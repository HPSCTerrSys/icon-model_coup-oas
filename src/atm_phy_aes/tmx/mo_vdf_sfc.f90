! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Classes and functions for the turbulent mixing package (tmx)

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_vdf_sfc

  USE mo_kind,              ONLY: wp, vp, i1, i4
  USE mo_exception,         ONLY: message, finish
  USE mo_fortran_tools,     ONLY: init, copy
  USE mtime,                ONLY: t_datetime => datetime
  USE mo_timer,             ONLY: timer_start, timer_stop, ltimer
  USE mo_master_config,     ONLY: isrestart  ! TODO: use config instead of USE
  USE mo_tmx_process_class, ONLY: t_tmx_process
  USE mo_tmx_field_class,   ONLY: t_tmx_field, t_domain, isfc_oce, isfc_ice, isfc_lnd
  USE mo_vdf_sfc_memory,    ONLY: t_vdf_sfc_config, t_vdf_sfc_inputs, t_vdf_sfc_diags, &
                                  build_vdf_sfc_config, build_vdf_sfc_inputs, build_vdf_sfc_diags
  USE mo_aes_vdf_config,    ONLY: aes_vdf_config
#ifndef __NO_JSBACH__
  USE mo_cuda_graphs,       ONLY: t_cuda_graphs, id_captured, create_graphs, &
                                  begin_capture, end_capture, replay, reset
  USE mo_jsb_interface,     ONLY: invalidate_cuda_graphs
  USE mo_jsb_time,          ONLY: is_time_ltrig_rad_m1
#endif
  USE, INTRINSIC :: iso_c_binding, ONLY: c_loc

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_vdf_sfc, t_vdf_sfc_inputs, t_vdf_sfc_config, t_vdf_sfc_diags, t_vdf_aggregator

  INTEGER, PARAMETER :: MAX_NO_STATES = 2 !< Maximum number of states in t_vdf_sfc

  TYPE, EXTENDS(t_tmx_process) :: t_vdf_sfc
    TYPE(t_vdf_sfc_config), POINTER :: config => NULL()
    TYPE(t_vdf_sfc_inputs), POINTER :: inputs => NULL()
    TYPE(t_vdf_sfc_diags),  POINTER :: diagnostics  => NULL()
    !
    ! Supported diffusion variables
    !
    INTEGER :: tsfc_idx = 1
    INTEGER :: qsat_idx = 2
#ifndef __NO_JSBACH__
    TYPE(t_cuda_graphs) :: graphs
#endif
  CONTAINS
    PROCEDURE :: Init => Init_vdf_sfc
    PROCEDURE :: Compute
    PROCEDURE :: Compute_diagnostics
    PROCEDURE :: Update_diagnostics
  END TYPE t_vdf_sfc

  INTERFACE t_vdf_sfc
    MODULE PROCEDURE t_vdf_sfc_construct
  END INTERFACE

  TYPE t_vdf_aggregator
    INTEGER :: aggregation_queue
  CONTAINS
    PROCEDURE :: BeginAggregate => begin_averaging
    PROCEDURE :: EndAggregate   => end_averaging
    PROCEDURE :: Aggregate      => average_tiles
  END TYPE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_vdf_sfc'

CONTAINS

  FUNCTION t_vdf_sfc_construct(name, dt, domain) RESULT(result)

    CHARACTER(len=*), INTENT(in) :: name
    REAL(wp),         INTENT(in) :: dt
    TYPE(t_domain),      POINTER :: domain
    TYPE(t_vdf_sfc),     POINTER :: result

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_vdf_sfc_construct'

    ! CALL message(routine, '')

    ALLOCATE(t_vdf_sfc::result)
    result%max_no_states = MAX_NO_STATES
    !$ACC ENTER DATA COPYIN(result)
    ! Call Init of abstract parent class
    CALL result%Init_process(dt=dt, name=name, domain=domain)
    __acc_attach(result%domain)

    ! Initialize memory structures
    ALLOCATE(result%config)
    ALLOCATE(result%inputs)
    ALLOCATE(result%diagnostics)

    ! Initialize the data structures
    CALL build_vdf_sfc_config(result%config, result%domain)
    CALL build_vdf_sfc_inputs(result%inputs, result%domain)
    CALL build_vdf_sfc_diags(result%diagnostics, result%domain)

  END FUNCTION t_vdf_sfc_construct

  SUBROUTINE Init_vdf_sfc(this)

    CLASS(t_vdf_sfc), INTENT(inout), TARGET :: this

    REAL(wp), POINTER :: co2(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init'

    !CALL message(routine, '')

    ! Initialize structure for config variables (second scan)
    CALL build_vdf_sfc_config(this%config, this%domain)
    !$ACC ENTER DATA COPYIN(this%config)

    ! Initialize structure for input variables (second scan)
    CALL build_vdf_sfc_inputs(this%inputs, this%domain)
    !$ACC ENTER DATA COPYIN(this%inputs)

    ! Initialize structure for diagnostic variables (second scan)
    CALL build_vdf_sfc_diags(this%diagnostics, this%domain)
    !$ACC ENTER DATA COPYIN(this%diagnostics)

    ! TODO: simple initialization of CO2
    co2 => this%inputs%co2%Get_ptr_r2d()
    co2(:,:) = 348.0e-06_wp * 1.51919227_wp

  END SUBROUTINE Init_vdf_sfc

  SUBROUTINE Compute(this, datetime)

    USE mo_tmx_surface_interface, ONLY: &
      & update_land, update_sea_ice, compute_lw_rad_net, compute_sw_rad_net, compute_albedo, &
      & compute_sfc_fluxes, compute_sfc_sat_spec_humidity, compute_energy_fluxes
    USE mo_physical_constants, ONLY: albedoW ! TODO

    CLASS(t_vdf_sfc), INTENT(inout), TARGET :: this
    TYPE(t_datetime), OPTIONAL, INTENT(in), POINTER :: datetime     !< date and time at beginning of time step

    ! Local variables for data access
    TYPE(t_vdf_sfc_config), POINTER :: conf
    TYPE(t_vdf_sfc_inputs), POINTER :: ins
    TYPE(t_vdf_sfc_diags),  POINTER :: diags

    TYPE(t_vdf_aggregator) :: aggregator

    INTEGER :: jg, jtile, isfc, jc, jcl, jb, graph_id
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & old_tsfc, tend_tsfc, new_tsfc, &
      & new_qsfc, &
      & fract_tile, tsfc_tile

    ! Local pointers to config variables
    REAL(wp), POINTER :: dtime, cvd, cvv, wind_gustiness
    INTEGER,  POINTER :: nice_thickness_classes
    LOGICAL,  POINTER :: l_co2

    ! Local pointers to input variables
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & psfc, rlds, rsds, &
      & rvds_dir, rnds_dir, rpds_dir, &
      & rvds_dif, rnds_dif, rpds_dif, &
      & emissivity, cosmu0, &
      & ta, qa, pa, ua, va, rho_atm, &
      & rsfl, ssfl, co2, co2flx_ant, &
      & ice_thickness, &
      & ocean_u, ocean_v, ice_u, ice_v
    REAL(vp), POINTER, DIMENSION(:,:) :: &
      & dz

    ! Local pointers to diagnostic variables
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & tsfc, tsfc_rad, &
      & evapotrans, lhfl, shfl, &
      & ustress, vstress, &
      & ufts, ufvs, lwfl_up, swfl_up, &
      & q_snocpymlt_lnd, &
      & albvisdir, albvisdif, &
      & albnirdir, albnirdif, albedo, &
      & co2flx_nat, co2flx, &
      & q_ice_top, q_ice_bot, snow_thickness

    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & albvisdir_tile, albvisdif_tile, &
      & albnirdir_tile, albnirdif_tile, &
      & kh_tile, km_tile, &
      & interp_fac_2m_tile, interp_fac_10m_tile, &
      & interp_fac_tsfc_tile, &
      & wind_rel_tile, &
      & lhfl_tile, shfl_tile, &
      & ustress_tile, vstress_tile, &
      & evapotrans_tile, &
      & lwfl_net_tile, swfl_net_tile, &
      & co2flx_nat_tile, &
      & albedo_tile, &
      & rho_tile, wind10m_tile

    REAL(wp) :: &
      & new_tsfc_rad(this%domain%nproma,this%domain%nblks_c,this%domain%ntiles), &
      & new_tsfc_eff(this%domain%nproma,this%domain%nblks_c,this%domain%ntiles), &
      & lwfl_net    (this%domain%nproma,this%domain%nblks_c), &
      & swfl_net    (this%domain%nproma,this%domain%nblks_c)

    ! Indices for valid points
    INTEGER, POINTER, DIMENSION(:,:) :: nvalid
    INTEGER, POINTER, DIMENSION(:,:,:) :: indices

    INTEGER :: acc_async_queues(this%domain%ntiles) ! OACC queues to process tiles in parallel

    CHARACTER(len=*), PARAMETER :: routine = modname//':Compute'

    IF (ltimer) CALL timer_start(this%timer_compute)

    jg = this%domain%patch%id

    ! Access memory structures
    conf => this%config
    ins => this%inputs
    diags => this%diagnostics

    ! Get pointers to config variables
    dtime => conf%dtime%Get_ptr_r0d()
    cvd => conf%cvd%Get_ptr_r0d()
    cvv => conf%cvv%Get_ptr_r0d()
    wind_gustiness => conf%wind_gustiness%Get_ptr_r0d()
    nice_thickness_classes => conf%nice_thickness_classes%Get_ptr_i0d()
    l_co2 => conf%l_co2%Get_ptr_l0d()

    ! Get pointers to input variables
    ta => ins%ta%Get_ptr_r2d()
    qa => ins%qa%Get_ptr_r2d()
    ua => ins%ua%Get_ptr_r2d()
    va => ins%va%Get_ptr_r2d()
    rho_atm => ins%rho_atm%Get_ptr_r2d()
    pa => ins%pa%Get_ptr_r2d()
    psfc => ins%psfc%Get_ptr_r2d()
    dz => ins%dz%Get_ptr_v2d()
    rlds => ins%rlds%Get_ptr_r2d()
    rsds => ins%rsds%Get_ptr_r2d()
    rvds_dir => ins%rvds_dir%Get_ptr_r2d()
    rnds_dir => ins%rnds_dir%Get_ptr_r2d()
    rpds_dir => ins%rpds_dir%Get_ptr_r2d()
    rvds_dif => ins%rvds_dif%Get_ptr_r2d()
    rnds_dif => ins%rnds_dif%Get_ptr_r2d()
    rpds_dif => ins%rpds_dif%Get_ptr_r2d()
    emissivity => ins%emissivity%Get_ptr_r2d()
    cosmu0 => ins%cosmu0%Get_ptr_r2d()
    rsfl => ins%rsfl%Get_ptr_r2d()
    ssfl => ins%ssfl%Get_ptr_r2d()
    co2 => ins%co2%Get_ptr_r2d()
    co2flx_ant => ins%co2flx_ant%Get_ptr_r2d()
    co2flx_nat => diags%co2flx_nat%Get_ptr_r2d()
    co2flx     => diags%co2flx%Get_ptr_r2d()
    ice_thickness => ins%ice_thickness%Get_ptr_r2d()
    ocean_u => ins%ocean_u%Get_ptr_r2d()
    ocean_v => ins%ocean_v%Get_ptr_r2d()
    ice_u => ins%ice_u%Get_ptr_r2d()
    ice_v => ins%ice_v%Get_ptr_r2d()

    ! Get pointers to tile-based inputs
    fract_tile => ins%fract_tile%Get_ptr_r3d()
    tsfc_tile => ins%tsfc_tile%Get_ptr_r3d()

    ! Get pointers to diagnostic variables
    tsfc => diags%tsfc%Get_ptr_r2d()
    tsfc_rad => diags%tsfc_rad%Get_ptr_r2d()
    evapotrans => diags%evapotrans%Get_ptr_r2d()
    lhfl => diags%lhfl%Get_ptr_r2d()
    shfl => diags%shfl%Get_ptr_r2d()
    ustress => diags%ustress%Get_ptr_r2d()
    vstress => diags%vstress%Get_ptr_r2d()
    ufts => diags%ufts%Get_ptr_r2d()
    ufvs => diags%ufvs%Get_ptr_r2d()
    lwfl_up => diags%lwfl_up%Get_ptr_r2d()
    swfl_up => diags%swfl_up%Get_ptr_r2d()
    albvisdir => diags%albvisdir%Get_ptr_r2d()
    albvisdif => diags%albvisdif%Get_ptr_r2d()
    albnirdir => diags%albnirdir%Get_ptr_r2d()
    albnirdif => diags%albnirdif%Get_ptr_r2d()
    albedo => diags%albedo%Get_ptr_r2d()
    q_ice_top => diags%q_ice_top%Get_ptr_r2d()
    q_ice_bot => diags%q_ice_bot%Get_ptr_r2d()
    snow_thickness => diags%snow_thickness%Get_ptr_r2d()
    q_snocpymlt_lnd => diags%q_snocpymlt_lnd%Get_ptr_r2d()

    ! Get pointers to tile-based diagnostics
    albvisdir_tile => diags%albvisdir_tile%Get_ptr_r3d()
    albvisdif_tile => diags%albvisdif_tile%Get_ptr_r3d()
    albnirdir_tile => diags%albnirdir_tile%Get_ptr_r3d()
    albnirdif_tile => diags%albnirdif_tile%Get_ptr_r3d()
    albedo_tile => diags%albedo_tile%Get_ptr_r3d()
    kh_tile => diags%kh_tile%Get_ptr_r3d()
    km_tile => diags%km_tile%Get_ptr_r3d()
    interp_fac_2m_tile => diags%interp_fac_2m_tile%Get_ptr_r3d()
    interp_fac_10m_tile => diags%interp_fac_10m_tile%Get_ptr_r3d()
    interp_fac_tsfc_tile => diags%interp_fac_tsfc_tile%Get_ptr_r3d()
    lhfl_tile => diags%lhfl_tile%Get_ptr_r3d()
    shfl_tile => diags%shfl_tile%Get_ptr_r3d()
    wind_rel_tile => diags%wind_rel_tile%Get_ptr_r3d()
    ustress_tile => diags%ustress_tile%Get_ptr_r3d()
    vstress_tile => diags%vstress_tile%Get_ptr_r3d()
    evapotrans_tile => diags%evapotrans_tile%Get_ptr_r3d()
    lwfl_net_tile => diags%lwfl_net_tile%Get_ptr_r3d()
    swfl_net_tile => diags%swfl_net_tile%Get_ptr_r3d()
    rho_tile => diags%rho_tile%Get_ptr_r3d()
    wind10m_tile => diags%wind10m_tile%Get_ptr_r3d()
    co2flx_nat_tile => diags%co2flx_nat_tile%Get_ptr_r3d()

    ! Get pointers to index arrays
    nvalid => diags%nvalid%Get_ptr_i2d()
    indices => diags%indices%Get_ptr_i3d()

    ! Get pointers to state variables
    old_tsfc  => this%states    (this%tsfc_idx)%p%Get_ptr_r3d()
    tend_tsfc => this%tendencies(this%tsfc_idx)%p%Get_ptr_r3d()
    new_tsfc  => this%new_states(this%tsfc_idx)%p%Get_ptr_r3d()
    new_qsfc  => this%new_states(this%qsat_idx)%p%Get_ptr_r3d()

    graph_id = -1
#ifndef __NO_JSBACH__
    IF (aes_vdf_config(jg)%lcuda_graph_vdf .AND. .NOT. this%is_initial_time) THEN
      IF (.NOT. this%graphs%initialized) THEN
        CALL create_graphs(this%graphs, 3, modname)
      END IF
      IF (invalidate_cuda_graphs) THEN
        CALL reset(this%graphs)
      ELSE
        graph_id = id_captured( this%graphs, ptr_keys=(/ C_LOC(ta(1,1)) /), &
          int_keys=(/ 1, merge(1, 0, is_time_ltrig_rad_m1(datetime, dtime, jg, .FALSE.)) /) )
        IF (graph_id > 0) THEN
          CALL replay(this%graphs, graph_id, 1)
          !$ACC WAIT(1)
          IF (ltimer) CALL timer_stop(this%timer_compute)
          RETURN
        ELSE
          CALL begin_capture( this%graphs, 1, ptr_keys=(/ C_LOC(ta(1,1)) /), &
            int_keys=(/ 1, merge(1, 0, is_time_ltrig_rad_m1(datetime, dtime, jg, .FALSE.)) /) )
        END IF
      END IF
    END IF
#endif

    !$ACC DATA CREATE(new_tsfc_rad, new_tsfc_eff, lwfl_net, swfl_net) &
    !$ACC   PRESENT(old_tsfc, tend_tsfc, new_tsfc, new_qsfc, tsfc_rad, lwfl_up, swfl_up, rlds, rsds) ASYNC(1)

    ! Process tiles in multiple OpenACC queues
    DO jtile=1,this%domain%ntiles
      isfc = this%domain%sfc_types(jtile)

      SELECT CASE(isfc)
      CASE(isfc_oce)
        acc_async_queues(jtile) = 98
      CASE(isfc_ice)
        acc_async_queues(jtile) = 99
      CASE(isfc_lnd)
        acc_async_queues(jtile) = 1 ! land tile should run on queue 1 to match jsbach
      END SELECT

      !$ACC WAIT(1) ASYNC(acc_async_queues(jtile))
    END DO

    DO jtile=1,this%domain%ntiles
      isfc = this%domain%sfc_types(jtile)

      SELECT CASE(isfc)
      CASE(isfc_oce)
        ! Ocean surface temperature is calculated outside of this, set tendency to zero
!$OMP PARALLEL
        CALL init(tend_tsfc(:,:,jtile),                         lacc=.TRUE., opt_acc_async_queue=acc_async_queues(jtile))
        CALL copy(old_tsfc(:,:,jtile), new_tsfc(:,:,jtile),     lacc=.TRUE., opt_acc_async_queue=acc_async_queues(jtile))
        CALL copy(old_tsfc(:,:,jtile), new_tsfc_rad(:,:,jtile), lacc=.TRUE., opt_acc_async_queue=acc_async_queues(jtile))
        CALL copy(old_tsfc(:,:,jtile), new_tsfc_eff(:,:,jtile), lacc=.TRUE., opt_acc_async_queue=acc_async_queues(jtile))
!$OMP END PARALLEL
        CALL compute_sfc_sat_spec_humidity(.FALSE., this%domain, this%domain%sfc_types(jtile), &
          & nvalid(:,jtile), indices(:,:,jtile), &
          & psfc(:,:), new_tsfc(:,:,jtile), new_qsfc(:,:,jtile), &
          & opt_acc_async_queue=acc_async_queues(jtile))
        ! TODO: This should be replaced by routine mo_surface_ocean:update_albedo_ocean from ECHAM6.2
!$OMP PARALLEL
        CALL init(albvisdir_tile(:,:,jtile), albedoW, lacc=.TRUE., opt_acc_async_queue=acc_async_queues(jtile))
        CALL init(albvisdif_tile(:,:,jtile), albedoW, lacc=.TRUE., opt_acc_async_queue=acc_async_queues(jtile))
        CALL init(albnirdir_tile(:,:,jtile), albedoW, lacc=.TRUE., opt_acc_async_queue=acc_async_queues(jtile))
        CALL init(albnirdif_tile(:,:,jtile), albedoW, lacc=.TRUE., opt_acc_async_queue=acc_async_queues(jtile))
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jc, jcl, jb) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = this%domain%i_startblk_c, this%domain%i_endblk_c
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR PRIVATE(jc) ASYNC(acc_async_queues(jtile))
          DO jcl = 1, nvalid(jb,jtile)
            jc = indices(jcl,jb,jtile)
            interp_fac_tsfc_tile(jc,jb,jtile) = interp_fac_2m_tile(jc,jb,jtile) * new_tsfc(jc,jb,jtile)
          END DO
          !$ACC END PARALLEL LOOP
        END DO
!$OMP END PARALLEL DO

      CASE(isfc_ice)
        IF (nice_thickness_classes /= 1) CALL finish(routine, 'Only one ice thickness class (kice) implemented!')
!$OMP PARALLEL
        CALL init(tend_tsfc(:,:,jtile), lacc=.TRUE., opt_acc_async_queue=acc_async_queues(jtile))
!$OMP END PARALLEL
        CALL update_sea_ice(this%domain, this%dt, &
          & nvalid(:,jtile), indices(:,:,jtile), &
          & old_tsfc(:,:,jtile), &
          & lwfl_net_tile(:,:,jtile), swfl_net_tile(:,:,jtile), &
          & lhfl_tile(:,:,jtile), shfl_tile(:,:,jtile), &
          & ssfl, ice_thickness, &
          & emissivity, &
          ! inout &
          & snow_thickness, &
          ! out &
          & new_tsfc(:,:,jtile), &
          & q_ice_top(:,:), q_ice_bot(:,:), &
          & albvisdir_tile(:,:,jtile), albvisdif_tile(:,:,jtile), &
          & albnirdir_tile(:,:,jtile), albnirdif_tile(:,:,jtile),  &
          & opt_acc_async_queue=acc_async_queues(jtile) &
          & )

!$OMP PARALLEL DO PRIVATE(jc, jcl, jb) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = this%domain%i_startblk_c, this%domain%i_endblk_c
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR PRIVATE(jc) ASYNC(acc_async_queues(jtile))
          DO jcl = 1, nvalid(jb,jtile)
            jc = indices(jcl,jb,jtile)
            new_tsfc_rad(jc,jb,jtile) = new_tsfc(jc,jb,jtile)
            new_tsfc_eff(jc,jb,jtile) = new_tsfc(jc,jb,jtile)
            tend_tsfc(jc,jb,jtile) = (new_tsfc(jc,jb,jtile) - old_tsfc(jc,jb,jtile)) / dtime
            interp_fac_tsfc_tile(jc,jb,jtile) = interp_fac_2m_tile(jc,jb,jtile) * new_tsfc(jc,jb,jtile)
          END DO
          !$ACC END PARALLEL LOOP
        END DO
!$OMP END PARALLEL DO

        CALL compute_sfc_sat_spec_humidity(.FALSE., this%domain, this%domain%sfc_types(jtile), &
          & nvalid(:,jtile), indices(:,:,jtile), &
          & psfc(:,:), new_tsfc(:,:,jtile), new_qsfc(:,:,jtile), &
          & opt_acc_async_queue=acc_async_queues(jtile))
      CASE(isfc_lnd)
!$OMP PARALLEL
        CALL init(tend_tsfc(:,:,jtile), lacc=.TRUE., opt_acc_async_queue=acc_async_queues(jtile))
!$OMP END PARALLEL
        CALL update_land(jg, this%domain, datetime, this%dt, cvd, &
          & dz, psfc, ta, qa, pa, &
          & rsfl, ssfl, &
          & rlds, &
          & rvds_dir, rnds_dir, rpds_dir, &
          & rvds_dif, rnds_dif, rpds_dif, &
          & cosmu0, &
          & wind_rel_tile(:,:,jtile), wind10m_tile(:,:,jtile), rho_tile(:,:,jtile), co2, &
          ! out
          & new_tsfc(:,:,jtile), new_tsfc_rad(:,:,jtile), new_tsfc_eff(:,:,jtile), q_snocpymlt_lnd, &
          & new_qsfc(:,:,jtile), &
          & albvisdir_tile(:,:,jtile), albvisdif_tile(:,:,jtile), &
          & albnirdir_tile(:,:,jtile), albnirdif_tile(:,:,jtile), &
          & kh_tile(:,:,jtile), km_tile(:,:,jtile), &
          & interp_fac_2m_tile(:,:,jtile), interp_fac_10m_tile(:,:,jtile), &
          & interp_fac_tsfc_tile(:,:,jtile), co2flx_nat_tile(:,:,jtile) &
          & )

!$OMP PARALLEL DO PRIVATE(jc, jcl, jb) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = this%domain%i_startblk_c, this%domain%i_endblk_c
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR PRIVATE(jc) ASYNC(acc_async_queues(jtile))
          DO jcl = 1, nvalid(jb,jtile)
            jc = indices(jcl,jb,jtile)
            tend_tsfc(jc,jb,jtile) = (new_tsfc(jc,jb,jtile) - old_tsfc(jc,jb,jtile)) / dtime
          END DO
          !$ACC END PARALLEL LOOP
        END DO
!$OMP END PARALLEL DO

      END SELECT

!$OMP PARALLEL DO PRIVATE(jc, jcl, jb) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = this%domain%i_startblk_c, this%domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR PRIVATE(jc) ASYNC(acc_async_queues(jtile))
        DO jcl = 1, nvalid(jb,jtile)
          jc = indices(jcl,jb,jtile)
          new_tsfc_rad(jc,jb,jtile) = new_tsfc_rad(jc,jb,jtile)**4._wp
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO

      ! Compute surface fluxes for heat, water vapor and momentum from new state
      ! Note: for land, latent and sensible heat fluxes are directly taken from land model
      ! because they are part of the surface energy balance equation that is solved there
      CALL compute_sfc_fluxes( &
        ! Input
        & this%domain, this%domain%sfc_types(jtile), nvalid(:,jtile), indices(:,:,jtile), &
        & wind_gustiness, &
        & cvd, &
        & ua, va, &
        & ta, qa, wind_rel_tile(:,:,jtile), &
        & ocean_u, ocean_v, ice_u, ice_v, &
        & rho_tile(:,:,jtile), &
        & new_qsfc(:,:,jtile), new_tsfc(:,:,jtile), &
        & kh_tile(:,:,jtile), km_tile(:,:,jtile), &
        ! Output
        & evapotrans_tile(:,:,jtile), lhfl_tile(:,:,jtile), shfl_tile(:,:,jtile), &
        & ustress_tile(:,:,jtile), vstress_tile(:,:,jtile), &
        & opt_acc_async_queue=acc_async_queues(jtile))

      CALL compute_lw_rad_net(this%domain, nvalid(:,jtile), indices(:,:,jtile), &
        & emissivity(:,:), rlds(:,:), new_tsfc_eff(:,:,jtile), lwfl_net_tile(:,:,jtile), &
        & opt_acc_async_queue=acc_async_queues(jtile) )

      CALL compute_sw_rad_net(this%domain, nvalid(:,jtile), indices(:,:,jtile), &
        & rvds_dir(:,:), rvds_dif(:,:), rnds_dir(:,:), rnds_dif(:,:), &
        & albvisdir_tile(:,:,jtile), albvisdif_tile(:,:,jtile), &
        & albnirdir_tile(:,:,jtile), albnirdif_tile(:,:,jtile), &
        & swfl_net_tile(:,:,jtile), &
        & opt_acc_async_queue=acc_async_queues(jtile))

        CALL compute_albedo(this%domain, nvalid(:,jtile), indices(:,:,jtile), &
        & rsds(:,:), rvds_dir(:,:), rvds_dif(:,:), rnds_dir(:,:), rnds_dif(:,:), &
        & albvisdir_tile(:,:,jtile), albvisdif_tile(:,:,jtile), &
        & albnirdir_tile(:,:,jtile), albnirdif_tile(:,:,jtile), &
        & albedo_tile(:,:,jtile), &
        & opt_acc_async_queue=acc_async_queues(jtile))

    END DO

    ! Join streams before aggregating
    DO jtile=1,this%domain%ntiles
      !$ACC WAIT(acc_async_queues(jtile)) ASYNC(1)
    END DO

    CALL aggregator%BeginAggregate()

    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, new_tsfc(:,:,:), tsfc, 'tsfc')
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, new_tsfc_rad(:,:,:), tsfc_rad, 'tsfc_rad4')

    ! Aggregate surface fluxes
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, evapotrans_tile, evapotrans,'evapotrans')
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, lhfl_tile,       lhfl, 'lhfl')
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, shfl_tile,       shfl, 'shfl')
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, ustress_tile,    ustress, 'ustress')
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, vstress_tile,    vstress, 'vstress')
    IF (l_co2) THEN
      CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, co2flx_nat_tile, co2flx_nat, 'CO2 nat')
    END IF
    !
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, albvisdir_tile,  albvisdir, 'albvisdir')
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, albvisdif_tile,  albvisdif, 'albvisdif')
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, albnirdir_tile,  albnirdir, 'albnirdir')
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, albnirdif_tile,  albnirdif, 'albnirdif')
    !
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, albedo_tile,     albedo, 'albedo')
    !
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, lwfl_net_tile,   lwfl_net, 'lwfl_net')
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, swfl_net_tile,   swfl_net, 'swfl_net')

    CALL aggregator%EndAggregate()

!$OMP PARALLEL DO PRIVATE(jc, jb) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = this%domain%i_startblk_c, this%domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COPYIN(l_co2)
      DO jc = this%domain%i_startidx_c(jb), this%domain%i_endidx_c(jb)
        tsfc_rad(jc,jb) = tsfc_rad(jc,jb)**0.25_wp
        lwfl_up(jc,jb) = rlds(jc,jb) - lwfl_net(jc,jb)
        swfl_up(jc,jb) = rsds(jc,jb) - swfl_net(jc,jb)
        IF (l_co2) THEN
          co2flx(jc,jb)  = co2flx_nat(jc,jb) + co2flx_ant(jc,jb)
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

    CALL compute_energy_fluxes( &
      & this%domain, cvv, cvd, &
      & shfl, evapotrans, ta, rho_atm, &
      & ufts, ufvs)

    !$ACC END DATA

#ifndef __NO_JSBACH__
    IF (graph_id == 0) THEN
      graph_id = end_capture(this%graphs)
      CALL replay(this%graphs, graph_id, 1)
    END IF
#endif

    !$ACC WAIT(1)

    IF (ltimer) CALL timer_stop(this%timer_compute)

  END SUBROUTINE Compute

  SUBROUTINE Compute_diagnostics(this, datetime)

    USE mo_vdf_diag_smag, ONLY: compute_wind_speed, compute_atm_potential_temperature, compute_sfc_density, &
      & compute_sfc_potential_temperature, compute_sfc_exchange_coefficients, compute_moist_richardson
    USE mo_tmx_surface_interface, ONLY: &
      & compute_lw_rad_net, compute_sw_rad_net, compute_sfc_fluxes, compute_sfc_sat_spec_humidity, compute_sfc_roughness, &
      & update_land, compute_10m_wind

    CLASS(t_vdf_sfc), INTENT(inout), TARGET :: this
    TYPE(t_datetime), OPTIONAL, INTENT(in), POINTER :: datetime

    ! Local variables for data access
    TYPE(t_vdf_sfc_config), POINTER :: conf
    TYPE(t_vdf_sfc_inputs), POINTER :: ins
    TYPE(t_vdf_sfc_diags),  POINTER :: diags

    TYPE(t_vdf_aggregator) :: aggregator

    INTEGER :: jg, jtile, graph_id

    REAL(wp) :: rough_min

    ! Local pointers to config variables
    REAL(wp), POINTER :: fsl, min_rough, wind_gustiness, &
      & rough_m_oce, rough_m_ice, min_sfc_wind, cvd
    INTEGER, POINTER :: nice_thickness_classes

    ! Local pointers to input variables
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & ta, tv, qa, ua, va, pa, psfc, zf, zh, &
      & rlds, rsds, rvds_dir, rnds_dir, rpds_dir, &
      & rvds_dif, rnds_dif, rpds_dif, &
      & emissivity, cosmu0, rsfl, ssfl, co2, &
      & ocean_u, ocean_v, ice_u, ice_v
    REAL(vp), POINTER, DIMENSION(:,:) :: &
      & dz

    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & fract_tile, tsfc_tile

    ! Local pointers to diagnostic variables
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & theta_atm, thetav_atm, rough_m, rough_h

    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & qsat_tile, rho_tile, theta_tile, thetav_tile, &
      & moist_rich_tile, rough_h_tile, rough_m_tile, &
      & km_tile, kh_tile, interp_fac_2m_tile, interp_fac_10m_tile, &
      & evapotrans_tile, lhfl_tile, shfl_tile, &
      & wind_rel_tile, ustress_tile, vstress_tile, &
      & u10m_tile, v10m_tile, wind10m_tile, &
      & albvisdir_tile, albvisdif_tile, &
      & albnirdir_tile, albnirdif_tile, &
      & lwfl_net_tile, swfl_net_tile

    ! Indices for valid points
    INTEGER, POINTER, DIMENSION(:,:) :: nvalid
    INTEGER, POINTER, DIMENSION(:,:,:) :: indices

    ! Get pointers to state variables
    REAL(wp), POINTER, DIMENSION(:,:,:) :: old_tsfc, old_qsat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Compute_diagnostics'

    ! CALL message(routine, '')

    IF (ltimer) CALL timer_start(this%timer_diagnostics)

    jg = this%domain%patch%id

    ! Access memory structures
    conf => this%config
    ins => this%inputs
    diags => this%diagnostics

    ! Get pointers to config variables
    fsl => conf%fsl%Get_ptr_r0d()
    min_rough => conf%min_rough%Get_ptr_r0d()
    wind_gustiness => conf%wind_gustiness%Get_ptr_r0d()
    rough_m_oce => conf%rough_m_oce%Get_ptr_r0d()
    rough_m_ice => conf%rough_m_ice%Get_ptr_r0d()
    min_sfc_wind => conf%min_sfc_wind%Get_ptr_r0d()
    cvd => conf%cvd%Get_ptr_r0d()
    nice_thickness_classes => conf%nice_thickness_classes%Get_ptr_i0d()

    ! Get pointers to input variables
    ta => ins%ta%Get_ptr_r2d()
    tv => ins%tv%Get_ptr_r2d()
    qa => ins%qa%Get_ptr_r2d()
    ua => ins%ua%Get_ptr_r2d()
    va => ins%va%Get_ptr_r2d()
    pa => ins%pa%Get_ptr_r2d()
    psfc => ins%psfc%Get_ptr_r2d()
    dz => ins%dz%Get_ptr_v2d()
    zf => ins%zf%Get_ptr_r2d()
    zh => ins%zh%Get_ptr_r2d()
    rlds => ins%rlds%Get_ptr_r2d()
    rsds => ins%rsds%Get_ptr_r2d()
    rvds_dir => ins%rvds_dir%Get_ptr_r2d()
    rnds_dir => ins%rnds_dir%Get_ptr_r2d()
    rpds_dir => ins%rpds_dir%Get_ptr_r2d()
    rvds_dif => ins%rvds_dif%Get_ptr_r2d()
    rnds_dif => ins%rnds_dif%Get_ptr_r2d()
    rpds_dif => ins%rpds_dif%Get_ptr_r2d()
    emissivity => ins%emissivity%Get_ptr_r2d()
    cosmu0 => ins%cosmu0%Get_ptr_r2d()
    rsfl => ins%rsfl%Get_ptr_r2d()
    ssfl => ins%ssfl%Get_ptr_r2d()
    co2 => ins%co2%Get_ptr_r2d()
    ocean_u => ins%ocean_u%Get_ptr_r2d()
    ocean_v => ins%ocean_v%Get_ptr_r2d()
    ice_u => ins%ice_u%Get_ptr_r2d()
    ice_v => ins%ice_v%Get_ptr_r2d()

    ! Get pointers to tile-based inputs
    fract_tile => ins%fract_tile%Get_ptr_r3d()
    tsfc_tile => ins%tsfc_tile%Get_ptr_r3d()

    ! Get pointers to diagnostic variables
    theta_atm => diags%theta_atm%Get_ptr_r2d()
    thetav_atm => diags%thetav_atm%Get_ptr_r2d()
    rough_m => diags%rough_m%Get_ptr_r2d()
    rough_h => diags%rough_h%Get_ptr_r2d()

    ! Get pointers to tile-based diagnostics
    qsat_tile => diags%qsat_tile%Get_ptr_r3d()
    rho_tile => diags%rho_tile%Get_ptr_r3d()
    theta_tile => diags%theta_tile%Get_ptr_r3d()
    thetav_tile => diags%thetav_tile%Get_ptr_r3d()
    moist_rich_tile => diags%moist_rich_tile%Get_ptr_r3d()
    rough_h_tile => diags%rough_h_tile%Get_ptr_r3d()
    rough_m_tile => diags%rough_m_tile%Get_ptr_r3d()
    km_tile => diags%km_tile%Get_ptr_r3d()
    kh_tile => diags%kh_tile%Get_ptr_r3d()
    interp_fac_2m_tile => diags%interp_fac_2m_tile%Get_ptr_r3d()
    interp_fac_10m_tile => diags%interp_fac_10m_tile%Get_ptr_r3d()
    evapotrans_tile => diags%evapotrans_tile%Get_ptr_r3d()
    lhfl_tile => diags%lhfl_tile%Get_ptr_r3d()
    shfl_tile => diags%shfl_tile%Get_ptr_r3d()
    wind_rel_tile => diags%wind_rel_tile%Get_ptr_r3d()
    ustress_tile => diags%ustress_tile%Get_ptr_r3d()
    vstress_tile => diags%vstress_tile%Get_ptr_r3d()
    u10m_tile => diags%u10m_tile%Get_ptr_r3d()
    v10m_tile => diags%v10m_tile%Get_ptr_r3d()
    wind10m_tile => diags%wind10m_tile%Get_ptr_r3d()
    albvisdir_tile => diags%albvisdir_tile%Get_ptr_r3d()
    albvisdif_tile => diags%albvisdif_tile%Get_ptr_r3d()
    albnirdir_tile => diags%albnirdir_tile%Get_ptr_r3d()
    albnirdif_tile => diags%albnirdif_tile%Get_ptr_r3d()
    lwfl_net_tile => diags%lwfl_net_tile%Get_ptr_r3d()
    swfl_net_tile => diags%swfl_net_tile%Get_ptr_r3d()

    ! Get pointers to index arrays
    nvalid => diags%nvalid%Get_ptr_i2d()
    indices => diags%indices%Get_ptr_i3d()

    old_tsfc => this%states(this%tsfc_idx)%p%Get_ptr_r3d()
    old_qsat => this%states(this%qsat_idx)%p%Get_ptr_r3d()

    graph_id = -1
#ifndef __NO_JSBACH__
    IF (aes_vdf_config(jg)%lcuda_graph_vdf .AND. .NOT. this%is_initial_time) THEN
      IF (.NOT. this%graphs%initialized) THEN
        CALL create_graphs(this%graphs, 3, modname)
      END IF
      IF (invalidate_cuda_graphs) THEN
        CALL reset(this%graphs)
      ELSE
        graph_id = id_captured( this%graphs, ptr_keys=(/ C_LOC(ta(1,1)) /), int_keys=(/ 2, 0 /) )
        IF (graph_id > 0) THEN
          CALL replay(this%graphs, graph_id, 1)
          !$ACC WAIT(1)
          IF (ltimer) CALL timer_stop(this%timer_diagnostics)
          RETURN
        ELSE
          CALL begin_capture( this%graphs, 1, ptr_keys=(/ C_LOC(ta(1,1)) /), int_keys=(/ 2, 0 /) )
        END IF
      END IF
    END IF
#endif

    CALL compute_valid_indices(this%domain, fract_tile, nvalid, indices)

    CALL compute_atm_potential_temperature(this%domain, ta, tv, pa, &
      & theta_atm, thetav_atm)

    ! Process tiles in multiple OpenACC queues
    ! Since jsbach here is only called once in the beginning,
    ! we don't have to run the land tile on queue 1 like in `Compute`
    DO jtile=1,this%domain%ntiles
      !$ACC WAIT(1) ASYNC(jtile)
    END DO

    DO jtile=1,this%domain%ntiles

      ! relative wind speed
      CALL compute_wind_speed(this%domain, nvalid(:,jtile), indices(:,:,jtile), &
        & min_sfc_wind, this%domain%sfc_types(jtile), &
        & ua, va, ocean_u, ocean_v, ice_u, ice_v, wind_rel_tile(:,:,jtile), &
        & opt_acc_async_queue=jtile)

      ! Surface saturated humidity
      CALL compute_sfc_sat_spec_humidity(this%is_initial_time .AND. .NOT. isrestart(), this%domain, this%domain%sfc_types(jtile), &
        & nvalid(:,jtile), indices(:,:,jtile), &
        & psfc(:,:), tsfc_tile(:,:,jtile), qsat_tile(:,:,jtile), &
        & opt_acc_async_queue=jtile)

      ! Density at surface
      CALL compute_sfc_density(this%domain, nvalid(:,jtile), indices(:,:,jtile), &
        ! & qsat_tile(:,:,jtile), psfc(:,:), ta(:,:), rho_tile(:,:,jtile))
        & qsat_tile(:,:,jtile), psfc(:,:), tsfc_tile(:,:,jtile), rho_tile(:,:,jtile), &
        & opt_acc_async_queue=jtile)

      IF (this%is_initial_time .AND. .NOT. isrestart()) THEN

        ! This only happens once, no need to use multiple streams here
        !$ACC WAIT(jtile) ASYNC(1)

        IF (this%domain%sfc_types(jtile) == isfc_lnd) THEN
          ! Compute inital 10m wind for update_land
          CALL compute_sfc_potential_temperature(this%domain, nvalid(:,jtile), indices(:,:,jtile), &
            & psfc, tsfc_tile(:,:,jtile), qsat_tile(:,:,jtile), &
            & theta_tile(:,:,jtile), thetav_tile(:,:,jtile))
          CALL compute_moist_richardson(this%domain, nvalid(:,jtile), indices(:,:,jtile), &
            & fsl, zf, thetav_atm, thetav_tile(:,:,jtile), wind_rel_tile(:,:,jtile), &
            & moist_rich_tile(:,:,jtile))
          CALL compute_sfc_roughness(.TRUE., this%domain, this%domain%sfc_types(jtile), &
            & nvalid(:,jtile), indices(:,:,jtile), min_rough, rough_m_oce, rough_m_ice, &
            & wind_rel_tile(:,:,jtile), km_tile(:,:,jtile), &
            & rough_h_tile(:,:,jtile), rough_m_tile(:,:,jtile))
          CALL compute_sfc_exchange_coefficients( &
            ! Input
            & this%domain, nvalid(:,jtile), indices(:,:,jtile), dz, &
            & qa, &
            & theta_atm, wind_rel_tile(:,:,jtile), rough_m_tile(:,:,jtile), &
            & theta_tile(:,:,jtile), qsat_tile(:,:,jtile), &
            ! Output
            & km_tile(:,:,jtile), kh_tile(:,:,jtile), &
            & interp_fac_2m_tile(:,:,jtile), interp_fac_10m_tile(:,:,jtile))
          CALL compute_10m_wind( &
            & this%domain, nvalid(:,jtile), indices(:,:,jtile), &
            & ua, va, &
            & interp_fac_10m_tile(:,:,jtile), &
            & u10m_tile(:,:,jtile), v10m_tile(:,:,jtile), wind10m_tile(:,:,jtile) &
            & )

          ! Call land in quasi-diagnostic mode with a time step of 1 second
          CALL update_land(jg, this%domain, datetime, 1._wp, cvd, &
            & dz, psfc, ta, qa, pa, &
            & rsfl, ssfl, &
            & rlds, &
            & rvds_dir, rnds_dir, rpds_dir, &
            & rvds_dif, rnds_dif, rpds_dif, &
            & cosmu0, &
            & wind_rel_tile(:,:,jtile), wind10m_tile(:,:,jtile), rho_tile(:,:,jtile), co2, &
            ! out
            & tsfc=tsfc_tile(:,:,jtile), &
            & qsat=qsat_tile(:,:,jtile) &
            & )

          !$ACC WAIT(1) ASYNC(jtile)
        END IF
      END IF

      ! Surface potential temperature
      CALL compute_sfc_potential_temperature(this%domain, nvalid(:,jtile), indices(:,:,jtile), &
        & psfc(:,:), tsfc_tile(:,:,jtile), qsat_tile(:,:,jtile), &
        & theta_tile(:,:,jtile), thetav_tile(:,:,jtile), &
        & opt_acc_async_queue=jtile)

      ! Moist Richardson number
      CALL compute_moist_richardson(this%domain, nvalid(:,jtile), indices(:,:,jtile), &
        & fsl, zf(:,:), thetav_atm(:,:), thetav_tile(:,:,jtile), wind_rel_tile(:,:,jtile), &
        & moist_rich_tile(:,:,jtile), &
        & opt_acc_async_queue=jtile)

      ! Surface roughness length
      IF (this%domain%sfc_types(jtile) == isfc_oce) THEN
        rough_min = min_rough
      ELSE
        rough_min = 0._wp
      END IF
      ! Uses old value of km_tile before computation of new exchange coefficients (only in case of ocean)
      CALL compute_sfc_roughness(this%is_initial_time .AND. .NOT. isrestart(), this%domain, this%domain%sfc_types(jtile), &
        & nvalid(:,jtile), indices(:,:,jtile), rough_min, rough_m_oce, rough_m_ice, &
        & wind_rel_tile(:,:,jtile), km_tile(:,:,jtile), &
        & rough_h_tile(:,:,jtile), rough_m_tile(:,:,jtile), &
        & opt_acc_async_queue=jtile)

      ! Surface exchange coefficients
      ! Note: coefficients for land will be computed in the land model itself
      IF (this%domain%sfc_types(jtile) /= isfc_lnd) THEN
        CALL compute_sfc_exchange_coefficients( &
            ! Input
          & this%domain, nvalid(:,jtile), indices(:,:,jtile), dz, &
          & qa, &
          & theta_atm, wind_rel_tile(:,:,jtile), rough_m_tile(:,:,jtile), &
          & theta_tile(:,:,jtile), qsat_tile(:,:,jtile), &
            ! Output
          & km_tile(:,:,jtile), kh_tile(:,:,jtile), &
          & interp_fac_2m_tile(:,:,jtile), interp_fac_10m_tile(:,:,jtile), &
          & opt_acc_async_queue=jtile)
      END IF

      IF (this%domain%sfc_types(jtile) == isfc_ice) THEN
        ! The ice model needs the latent and sensible heat fluxes as input (based on old state)
        CALL compute_sfc_fluxes( &
          ! Input
          & this%domain, isfc_ice, nvalid(:,jtile), indices(:,:,jtile), &
          & wind_gustiness, &
          & cvd, &
          & ua, va, &
          & ta, qa, wind_rel_tile(:,:,jtile), &
          & ocean_u, ocean_v, ice_u, ice_v, &
          & rho_tile(:,:,jtile), &
          & qsat_tile(:,:,jtile), tsfc_tile(:,:,jtile), &
          & kh_tile(:,:,jtile), km_tile(:,:,jtile), &
          ! Output
          & evapotrans_tile(:,:,jtile), lhfl_tile(:,:,jtile), shfl_tile(:,:,jtile), &
          & ustress_tile(:,:,jtile), vstress_tile(:,:,jtile), &
          & opt_acc_async_queue=jtile)

          ! The ice model needs the shortwave net surface flux as input (based on old state)
        CALL compute_sw_rad_net(this%domain, nvalid(:,jtile), indices(:,:,jtile), &
          & rvds_dir, rvds_dif, rnds_dir, rnds_dif, &
          & albvisdir_tile(:,:,jtile), albvisdif_tile(:,:,jtile), &
          & albnirdir_tile(:,:,jtile), albnirdif_tile(:,:,jtile), &
          ! Out
          & swfl_net_tile(:,:,jtile), &
          & opt_acc_async_queue=jtile)

          ! The ice model needs the longwave net surface flux as input (based on old state)
        CALL compute_lw_rad_net(this%domain, nvalid(:,jtile), indices(:,:,jtile), &
          & emissivity(:,:), rlds(:,:), tsfc_tile(:,:,jtile), lwfl_net_tile(:,:,jtile), &
          & opt_acc_async_queue=jtile )
      END IF

    END DO

    DO jtile=1,this%domain%ntiles
      !$ACC WAIT(jtile) ASYNC(1)
    END DO

    CALL aggregator%BeginAggregate()
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, rough_m_tile, rough_m, 'rough_m')
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, rough_h_tile, rough_h, 'rough_h')
    CALL aggregator%EndAggregate()

#ifndef __NO_JSBACH__
    IF (graph_id == 0) THEN
      graph_id = end_capture(this%graphs)
      CALL replay(this%graphs, graph_id, 1)
    END IF
#endif

    !$ACC WAIT(1)

    IF (ltimer) CALL timer_stop(this%timer_diagnostics)

  END SUBROUTINE Compute_diagnostics

  SUBROUTINE Update_diagnostics(this)

    CLASS(t_vdf_sfc), INTENT(inout), TARGET :: this

    ! Local variables for data access
    TYPE(t_vdf_sfc_config), POINTER :: conf
    TYPE(t_vdf_sfc_inputs), POINTER :: ins
    TYPE(t_vdf_sfc_diags),  POINTER :: diags

    ! Local pointers to tile-based inputs and diagnostics
    REAL(wp), POINTER, DIMENSION(:,:,:) :: fract_tile
    REAL(wp), POINTER, DIMENSION(:,:,:) :: km_tile, kh_tile

    ! Local pointers to grid-mean diagnostics
    REAL(wp), POINTER, DIMENSION(:,:) :: km, kh

    ! Indices for valid points
    INTEGER, POINTER, DIMENSION(:,:) :: nvalid
    INTEGER, POINTER, DIMENSION(:,:,:) :: indices

    TYPE(t_vdf_aggregator) :: aggregator

    CHARACTER(len=*), PARAMETER :: routine = modname//':Update_diagnostics'

    IF (ltimer) CALL timer_start(this%timer_diagnostics)

    ! Access memory structures
    conf => this%config
    ins => this%inputs
    diags => this%diagnostics

    ! Get pointers to tile-based inputs
    fract_tile => ins%fract_tile%Get_ptr_r3d()

    ! Get pointers to tile-based diagnostics
    km_tile => diags%km_tile%Get_ptr_r3d()
    kh_tile => diags%kh_tile%Get_ptr_r3d()

    ! Get pointers to grid-mean diagnostics
    km => diags%km%Get_ptr_r2d()
    kh => diags%kh%Get_ptr_r2d()

    ! Get pointers to index arrays
    nvalid => diags%nvalid%Get_ptr_i2d()
    indices => diags%indices%Get_ptr_i3d()

    ! TODO: Update roughness length for heat and momentum
    CALL aggregator%BeginAggregate()
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, km_tile, km, 'km')
    CALL aggregator%Aggregate(this%domain, fract_tile, nvalid, indices, kh_tile, kh, 'kh')
    CALL aggregator%EndAggregate()

    IF (this%is_initial_time) this%is_initial_time = .FALSE.

    IF (ltimer) CALL timer_stop(this%timer_diagnostics)

  END SUBROUTINE Update_diagnostics

  SUBROUTINE compute_valid_indices(domain, fract_tile, nvalid, indices)

    USE mo_index_list, ONLY: generate_index_list_batched

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp),       INTENT(in)          :: fract_tile(:,:,:)
    INTEGER,        INTENT(out)         :: nvalid(:,:), indices(:,:,:)

    INTEGER     :: ib, ic, ics, ice, jsfc, ntiles
    INTEGER(i1) :: pfrc_test(domain%nproma,domain%nblks_c, domain%ntiles)
    INTEGER     :: loidx    (domain%nproma,domain%ntiles), &
      &            is       (domain%ntiles)

    ntiles = domain%ntiles

    !$ACC DATA CREATE(pfrc_test, loidx, is) PRESENT(fract_tile, indices, nvalid) ASYNC(1)

!$OMP PARALLEL
    CALL init(nvalid, lacc=.TRUE.)
    CALL init(indices, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP DO PRIVATE(ib, ic, ics, ice, jsfc) ICON_OMP_RUNTIME_SCHEDULE
    DO ib = domain%i_startblk_c, domain%i_endblk_c
      ics = domain%i_startidx_c(ib)
      ice = domain%i_endidx_c  (ib)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jsfc = 1, ntiles
        DO ic = ics, ice
          pfrc_test(ic,ib,jsfc) = MERGE(1_i1, 0_i1, fract_tile(ic,ib,jsfc) > 0.0_wp)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END DO

!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(ib, ics, ice, jsfc, ic, loidx, is) ICON_OMP_RUNTIME_SCHEDULE
    DO ib = domain%i_startblk_c, domain%i_endblk_c
      ics = domain%i_startidx_c(ib)
      ice = domain%i_endidx_c  (ib)

      CALL generate_index_list_batched(pfrc_test(:,ib,:), loidx(:,:), ics, ice, &
        &   is(:), lacc=.TRUE., opt_acc_async_queue=1)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jsfc = 1, ntiles
        nvalid(ib,jsfc) = is(jsfc)
        !$ACC LOOP GANG VECTOR
        DO ic = 1, is(jsfc)
          indices(ic,ib,jsfc) = loidx(ic,jsfc)
        END DO
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END PARALLEL DO

    !$ACC END DATA

  END SUBROUTINE compute_valid_indices

  SUBROUTINE begin_averaging(this)
    CLASS(t_vdf_aggregator), INTENT(inout) :: this
    this%aggregation_queue = 1
  END SUBROUTINE

  SUBROUTINE end_averaging(this)
    CLASS(t_vdf_aggregator), INTENT(in) :: this
    INTEGER :: queue

    DO queue = 2, this%aggregation_queue
      !$ACC WAIT(queue) ASYNC(1)
    END DO
  END SUBROUTINE

  SUBROUTINE average_tiles(this, domain, fract_tile, nvalid, indices, var_in, var_out, msg)

    CLASS(t_vdf_aggregator), INTENT(inout) :: this
    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), INTENT(in)  :: &
      & fract_tile(:,:,:)
    REAL(wp), INTENT(in) :: &
      & var_in(:,:,:)
    INTEGER,  INTENT(in)  :: &
      & nvalid(:,:),         &
      & indices(:,:,:)
    REAL(wp), INTENT(out) :: &
      & var_out(:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: msg

    INTEGER :: jb, jbs, jbe, jls, js, ntiles, jsfc

    CHARACTER(len=*), PARAMETER :: routine = modname//':average_tiles'

    ! IF (PRESENT(msg)) CALL message(routine, 'working on '//msg)

    IF (domain%ntiles < 1) CALL finish(routine, 'This should not happen - ntiles < 1')

    ntiles = domain%ntiles
    jbs = domain%i_startblk_c
    jbe = domain%i_endblk_c

    this%aggregation_queue = this%aggregation_queue + 1
    !$ACC WAIT(1) ASYNC(this%aggregation_queue)

    IF (ntiles == 1) THEN

!$OMP PARALLEL
      CALL copy(var_in(:,:,1), var_out(:,:), lacc=.TRUE., opt_acc_async_queue=this%aggregation_queue)
!$OMP END PARALLEL

    ELSE

!$OMP PARALLEL
      CALL init(var_out, lacc=.TRUE., opt_acc_async_queue=this%aggregation_queue)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb,jls,js,jsfc) ICON_OMP_RUNTIME_SCHEDULE
      DO jb = jbs, jbe
        DO jsfc = 1, ntiles
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(this%aggregation_queue) PRIVATE(js)
          DO jls = 1, nvalid(jb,jsfc)
            js=indices(jls,jb,jsfc)
            var_out(js,jb) = var_out(js,jb) + fract_tile(js,jb,jsfc) * var_in(js,jb,jsfc)
          END DO
          !$ACC END PARALLEL LOOP
        END DO
      END DO
!$OMP END PARALLEL DO

    END IF

    ! IF (PRESENT(msg)) CALL message(routine, 'average complete for '//msg)

  END SUBROUTINE average_tiles

END MODULE mo_vdf_sfc
