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

! Subroutine interface_aes_tmx calls the turbulent mixing scheme
! and the surface schemes using tmx.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_interface_aes_tmx

  USE mo_kind                ,ONLY: wp, vp, dp
  USE mtime                  ,ONLY: datetime, OPERATOR(>)

  USE mo_exception           ,ONLY: finish, warning

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: ntracer
  USE mo_dynamics_config     ,ONLY: nnow, nnew, nnow_rcf
  USE mo_aes_phy_config      ,ONLY: aes_phy_config, aes_phy_tc, dt_zero
  USE mo_aes_phy_memory      ,ONLY: t_aes_phy_field, prm_field, prm_field_list, &
    &                               t_aes_phy_tend,  prm_tend, prm_tend_list, &
    &                               cdimissval

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_tmx

  USE mo_ccycle_config       ,ONLY: ccycle_config
  USE mo_physical_constants  ,ONLY: amco2, amd
  USE mo_bc_greenhouse_gases ,ONLY: ghg_co2mmr

  USE mo_run_config          ,ONLY: iqv, iqc, iqi, iqr, iqs, iqg, iqnc, iqni, iqt, ico2
  USE mo_vdf                 ,ONLY: t_vdf, new_vdf

  USE mo_aes_sfc_indices     ,ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_surface_diag        ,ONLY: nsurf_diag
  USE mo_aes_vdf_config      ,ONLY: aes_vdf_config
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_impl_constants_grf  ,ONLY: grf_bdywidth_c
  USE mo_impl_constants      ,ONLY: min_rlcell_int, min_rlcell,max_dom
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_nh_testcases_nml    ,ONLY: nh_test_name

  USE memman, ONLY: var_descriptor

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_aes_tmx, init_tmx, vdf_dom

  TYPE :: t_vdf_p
    TYPE(t_vdf), POINTER :: p
  END TYPE t_vdf_p

  TYPE(t_vdf_p) :: vdf_dom(max_dom)

  LOGICAL, SAVE :: l_init_or_restart = .TRUE.

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interface_aes_tmx'

CONTAINS

  SUBROUTINE interface_aes_tmx  (patch                ,&
       &                         is_in_sd_ed_interval ,&
       &                         is_active            ,&
       &                         datetime_old         ,&
       &                         dtime               )

    USE mo_physical_constants, ONLY: vmr_to_mmr_co2

    ! Arguments
    !
    TYPE(t_patch)   ,TARGET ,INTENT(inout) :: patch
    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    TYPE(datetime)          ,POINTER    :: datetime_old
    REAL(wp)                ,INTENT(in) :: dtime

    ! Shortcuts
    !
    LOGICAL :: lparamcpl, l2moment, ldtrad_gt0
    INTEGER :: fc_vdf
    TYPE(t_aes_phy_field)   ,POINTER    :: field
    TYPE(t_aes_phy_tend)    ,POINTER    :: tend

    TYPE(t_vdf), POINTER :: vdf

    ! Local variables
    !
    INTEGER  :: jg, jc, jk, jsfc, jt, jice, copy_nblks_c
    INTEGER  :: jb,jbs,jbe,jcs,jce,ncd,rls,rle
    INTEGER  :: nlev, nlevm1, nlevp1
    INTEGER  :: ntrac
    INTEGER  :: nice        ! for simplicity (ice classes)
    !
    ! Pointers to results (nproma,patch%nlev,patch%nblks_c)
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & tend_ta_vdf, tend_ua_vdf, tend_va_vdf
    ! Pointers to results (nproma,patch%nlev,patch%nblks_c,ntracer)
    REAL(wp), POINTER, DIMENSION(:,:,:,:) :: &
      & tend_tracer_vdf
    ! Pointers to results (nproma,patch%nlev+1,patch%nblks_c)
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & tend_wa_vdf
    ! Pointers to results (nproma,patch%nblks_c,nsfc_type)
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & tend_ts

    ! Local varaibles
    REAL(wp) :: q_rlw_impl(nproma,patch%nblks_c), &
      &         tend_ta_rlw_impl(nproma,patch%nblks_c), &
      &         zco2(nproma,patch%nblks_c)
    REAL(wp)          :: mmr_co2

    REAL(wp), POINTER :: ptr_r2d(:,:)

    !
    INTEGER, POINTER :: turb
    LOGICAL :: l_use_rad, l_co2

    CHARACTER(len=*), PARAMETER :: routine = modname//':interface_aes_tmx'

    IF (ltimer) CALL timer_start(timer_tmx)

    jg           = patch%id
    nlev         = patch%nlev

    ! associate pointers
    ! lparamcpl =  aes_phy_config(jg)%lparamcpl
    ! l2moment  =  aes_phy_config(jg)%l2moment
    ! ldtrad_gt0=  aes_phy_tc(jg)%dt_rad > dt_zero
    ! fc_vdf    =  aes_phy_config(jg)%fc_vdf
    field     => prm_field(jg)
    tend      => prm_tend (jg)
    vdf       => vdf_dom  (jg)%p

    nlevm1 = nlev-1
    nlevp1 = nlev+1
    ntrac  = ntracer-iqt+1  ! number of tracers excluding water vapour and hydrometeors

    nice   = prm_field(jg)%kice
    turb => aes_vdf_config(jg)%turb

    l_co2 = (iqt <= ico2 .AND. ico2 <= ntracer)

    IF ( is_in_sd_ed_interval ) THEN
      !
      IF ( is_active ) THEN

        l_use_rad = aes_phy_tc(jg)%dt_rad > dt_zero

        rls = grf_bdywidth_c + 1
        rle = min_rlcell_int

        jbs = patch%cells%start_block(rls)
        jbe = patch%cells%end_block  (rle)

        IF (.NOT. aes_vdf_config(jg)%use_tmx) THEN
          CALL finish(routine, 'ERROR: namelist parameter aes_vdf_config%use_tmx=.FALSE.!')
        END IF

        !$ACC DATA CREATE(q_rlw_impl, tend_ta_rlw_impl, zco2) &
        !$ACC   PRESENT(ccycle_config)

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce, mmr_co2) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs, jbe

          CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

          SELECT CASE (ccycle_config(jg)%iccycle)
            !
          CASE (0) ! no c-cycle
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              zco2(jc,jb) = 348.0e-06_wp * vmr_to_mmr_co2
            END DO
          CASE (1) ! c-cycle with interactive atm. co2 concentration
            !$ACC LOOP GANG VECTOR
            DO jc = jcs,jce
              zco2(jc,jb) = field% qtrc_phy(jc,nlev,jb,ico2)
            END DO
          CASE (2) ! c-cycle with prescribed atm. co2 concentration

            SELECT CASE (ccycle_config(jg)%ico2conc)
            CASE (2)
              mmr_co2 = ccycle_config(jg)%vmr_co2 * vmr_to_mmr_co2
              !
              !$ACC LOOP GANG VECTOR
              DO jc = jcs,jce
                zco2(jc,jb) = mmr_co2
              END DO
            CASE (4)
              !$ACC LOOP GANG VECTOR
              DO jc = jcs,jce
                zco2(jc,jb) = ghg_co2mmr
              END DO
            END SELECT

          END SELECT

          !$ACC END PARALLEL

        END DO
!$OMP END PARALLEL DO

        !
        ! Some variable pointers are swapped in mo_interface_iconam_aes.f90:interface_iconam_aes
        ! before and after the physics is called so that they are potentially pointing to targets
        ! that are different from the targets set during initialization of tmx.
        ! Here we need to make sure that the tmx variables point to the correct targets at each time step.
        !
        CALL vdf%atmo%states(vdf%atmo%tracer_idx)%p%Update(time_id=nnow_rcf(jg))
        CALL vdf%atmo%states(vdf%atmo%temp_idx)  %p%Update(time_id=nnow(jg))
        CALL vdf%atmo%states(vdf%atmo%wwind_idx) %p%Update(time_id=nnow(jg))

        CALL vdf%atmo%inputs%tracer_c%Update  (time_id=nnow_rcf(jg))
        CALL vdf%atmo%inputs%temp_c%Update    (time_id=nnow(jg))
        CALL vdf%atmo%inputs%w_wind_ic%Update (time_id=nnow(jg))
        CALL vdf%atmo%inputs%rho_c%Update     (time_id=nnew(jg))

        ! TODO: using vn from prog state instead of computing it in tmx leads to differences. Why?
        ! CALL vdf%atmo%inputs%vn_e%Update(name='vn', time_id=nnew(jg))

        CALL vdf%sfc%inputs%ta%Update(time_id=nnow(jg))
        CALL vdf%sfc%inputs%rho_atm%Update(time_id=nnew(jg))
        CALL vdf%sfc%inputs%qa%Update(time_id=nnow_rcf(jg))

        ptr_r2d => vdf%sfc%inputs%co2%Get_ptr_r2d()
!$OMP PARALLEL DO PRIVATE(jb, jc, jcs, jce) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs, jbe
          CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
          DO jc = jcs, jce
            ptr_r2d(jc,jb) = zco2(jc,jb)
          END DO
        END DO
!$OMP END PARALLEL DO

        !
        ! Move forward one time step
        !
        CALL vdf%Compute(datetime_old)

        ! Retrieve computed tendency for temperature
        tend_ta_vdf => vdf%atmo%Get_tendency_r3d(vdf%atmo%temp_idx)
        ! Retrieve computed tendencies for horizontal velocity
        tend_ua_vdf => vdf%atmo%Get_tendency_r3d(vdf%atmo%uwind_idx)
        tend_va_vdf => vdf%atmo%Get_tendency_r3d(vdf%atmo%vwind_idx)
        ! Retrieve computed tendency for vertical velocity
        tend_wa_vdf => vdf%atmo%Get_tendency_r3d(vdf%atmo%wwind_idx)
        ! Retrieve computed tendency for surface temperature on tiles (for output only)
        tend_ts => vdf%sfc%Get_tendency_r3d(vdf%sfc%tsfc_idx)
        ! Retrieve computed tendency for tracers
        tend_tracer_vdf => vdf%atmo%Get_tendency_r4d(vdf%atmo%tracer_idx)

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,jk) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs, jbe

          CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)

          !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(field%qtrc_phy) ASYNC(1)
          !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO jk = 1, nlev
            DO jc = jcs, jce

              ! Add tendency from turbulent transport to physics tendency
              ! (is used in iconam_aes interface)
              tend%ta_phy(jc,jk,jb) = tend%ta_phy(jc,jk,jb) + tend_ta_vdf(jc,jk,jb)
              ! Update physics state
              field%ta(jc,jk,jb) = field%ta(jc,jk,jb) + tend_ta_vdf(jc,jk,jb) * dtime

              ! Add tendency from turbulent transport to physics tendency
              ! (is used in iconam_aes interface)
              tend%qtrc_phy (jc,jk,jb,iqv) = tend%qtrc_phy (jc,jk,jb,iqv) + tend_tracer_vdf(jc,jk,jb,iqv)
              ! Update physics state
              field%qtrc_phy(jc,jk,jb,iqv) = field%qtrc_phy(jc,jk,jb,iqv) + tend_tracer_vdf(jc,jk,jb,iqv) * dtime

              ! Add tendency from turbulent transport to physics tendency
              ! (is used in iconam_aes interface)
              tend%qtrc_phy(jc,jk,jb,iqc) = tend%qtrc_phy(jc,jk,jb,iqc) + tend_tracer_vdf(jc,jk,jb,iqc)
              ! Update physics state
              field%qtrc_phy(jc,jk,jb,iqc) = field%qtrc_phy(jc,jk,jb,iqc) + tend_tracer_vdf(jc,jk,jb,iqc) * dtime

              ! Add tendency from turbulent transport to physics tendency
              ! (is used in iconam_aes interface)
              tend%qtrc_phy(jc,jk,jb,iqi) = tend%qtrc_phy(jc,jk,jb,iqi) + tend_tracer_vdf(jc,jk,jb,iqi)
              ! Update physics state
              field%qtrc_phy(jc,jk,jb,iqi) = field%qtrc_phy(jc,jk,jb,iqi) + tend_tracer_vdf(jc,jk,jb,iqi) * dtime
              !
              IF (l_co2) THEN
                IF (ccycle_config(jg)%iccycle /= 2) THEN
                  ! Add tendency from turbulent transport to physics tendency
                  ! (is used in iconam_aes interface)
                  tend%qtrc_phy(jc,jk,jb,ico2) = tend%qtrc_phy(jc,jk,jb,ico2) + tend_tracer_vdf(jc,jk,jb,ico2)
                  ! Update physics state
                  field%qtrc_phy(jc,jk,jb,ico2) = field%qtrc_phy(jc,jk,jb,ico2) + tend_tracer_vdf(jc,jk,jb,ico2) * dtime
                END IF
              END IF
              !
              tend%ua_phy(jc,jk,jb) = tend%ua_phy(jc,jk,jb) + tend_ua_vdf(jc,jk,jb)
              tend%va_phy(jc,jk,jb) = tend%va_phy(jc,jk,jb) + tend_va_vdf(jc,jk,jb)
              ! Update physics state
              field%ua(jc,jk,jb) = field%ua(jc,jk,jb) + tend_ua_vdf(jc,jk,jb) * dtime
              field%va(jc,jk,jb) = field%va(jc,jk,jb) + tend_va_vdf(jc,jk,jb) * dtime
              !
              tend%wa_phy(jc,jk,jb) = tend%wa_phy(jc,jk,jb) + tend_wa_vdf(jc,jk,jb)
              ! Update physics state
              field%wa(jc,jk,jb) = field%wa(jc,jk,jb) + tend_wa_vdf(jc,jk,jb) * dtime

            END DO
          END DO
          !$ACC END LOOP

          IF (ASSOCIATED(tend%ta_vdf)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = jcs, jce
                tend%ta_vdf(jc,jk,jb) = tend_ta_vdf(jc,jk,jb)
              END DO
            END DO
            !$ACC END LOOP
          END IF

          IF (ASSOCIATED(tend%ua_vdf)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = jcs, jce
                tend%ua_vdf(jc,jk,jb) = tend_ua_vdf(jc,jk,jb)
              END DO
            END DO
            !$ACC END LOOP
          END IF

          IF (ASSOCIATED(tend%va_vdf)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = jcs, jce
                tend%va_vdf(jc,jk,jb) = tend_va_vdf(jc,jk,jb)
              END DO
            END DO
            !$ACC END LOOP
          END IF

          IF (ASSOCIATED(tend%qtrc_vdf)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = jcs, jce
                tend%qtrc_vdf(jc,jk,jb,iqv) = tend_tracer_vdf(jc,jk,jb,iqv)
                tend%qtrc_vdf(jc,jk,jb,iqc) = tend_tracer_vdf(jc,jk,jb,iqc)
                tend%qtrc_vdf(jc,jk,jb,iqi) = tend_tracer_vdf(jc,jk,jb,iqi)
                IF (l_co2) tend%qtrc_vdf(jc,jk,jb,ico2) = tend_tracer_vdf(jc,jk,jb,ico2)
              END DO
            END DO
            !$ACC END LOOP
          END IF
          !$ACC END PARALLEL

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          IF (ASSOCIATED(tend%wa_vdf)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlevp1
              DO jc = jcs, jce
                tend%wa_vdf(jc,jk,jb) = tend_wa_vdf(jc,jk,jb)
              END DO
            END DO
            !$ACC END LOOP
          END IF

          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = jcs, jce

            q_rlw_impl(jc,jb) =                                               &
              &  ( (field%rld_rt(jc,nlev,jb)-field%rlu_rt(jc,nlev,jb))  & ! ( rln  from "radiation", at top of layer nlev
              &   -(field%rlds  (jc,jb)     -field%rlus  (jc,jb)     )) & !  -rlns from "radheating" and "update_surface")
              & -field%q_rlw_nlev(jc,jb)                                       ! -old heating in layer nlev from "radheating"

            ! Correction related to implicitness, due to the fact that surface model only used
            ! part of longwave radiation to compute new surface temperature
            !
            IF (ASSOCIATED(field%q_rlw_impl)) THEN
              field%q_rlw_impl(jc,jb) = q_rlw_impl(jc,jb)
            END IF

            ! convert heating
            tend_ta_rlw_impl(jc,jb) = q_rlw_impl(jc,jb) / field%cvair(jc,nlev,jb)
            !
            IF (ASSOCIATED(tend%ta_rlw_impl)) THEN
              tend%ta_rlw_impl(jc,jb) = tend_ta_rlw_impl(jc,jb)
            END IF

            ! for output: accumulate heating
            IF (ASSOCIATED(field% q_phy)) THEN
              field%q_phy(jc,nlev,jb) = field%q_phy(jc,nlev,jb) + q_rlw_impl(jc,jb)
            END IF
            IF (ASSOCIATED(field% q_phy_vi)) THEN
              field%q_phy_vi(jc,jb) = field%q_phy_vi(jc,jb) + q_rlw_impl(jc,jb)
            END IF
            !
            ! accumulate surface lw increment for later updating the model state
            ! but only if radiative processes are active.
            IF (l_use_rad) THEN
              ! use tendency to update the model state
              tend%ta_phy(jc,nlev,jb) = tend%ta_phy(jc,nlev,jb) + tend_ta_rlw_impl(jc,jb)
              !
              ! update physics state for input to the next physics process
              ! use tendency to update the physics state
              field%ta(jc,nlev,jb) = field%ta(jc,nlev,jb) + tend_ta_rlw_impl(jc,jb) * dtime
            END IF

          ENDDO
          !$ACC END LOOP

          !$ACC END PARALLEL

        END DO

!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,jsfc) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs, jbe

          CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO jsfc = 1, nsfc_type
            DO jc = jcs, jce
              IF (field%frac_tile(jc,jb,jsfc) > 0.0_wp) THEN
                field%ts_tile(jc,jb,jsfc) = field%ts_tile(jc,jb,jsfc) + tend_ts(jc,jb,jsfc) * dtime
              ELSE
                ! Set non-valid cells for each tile (fraction <= 0) to missing value for output purposes.
                ! They are also used in mo_aes_phy_diag:surface_fractions at the beginning of the next time step to initialize the
                ! sea ice surface temperature for newly formed sea ice.
                field%ts_tile(jc,jb,jsfc) = cdimissval
              END IF
            END DO
          END DO
          !$ACC END LOOP
          !$ACC END PARALLEL

        END DO
!$OMP END PARALLEL DO

        !$no EXIT DATA COPYOUT(tend_qtrc_vdf(iqv)%p)
        !$no EXIT DATA COPYOUT(tend_qtrc_vdf(iqc)%p)
        !$no EXIT DATA COPYOUT(tend_qtrc_vdf(iqi)%p)
        !$no EXIT DATA COPYOUT(tend_qtrc_vdf(1:3))

        !$ACC WAIT(1)

        !$ACC END DATA

        ! Retrieve diagnostics (mostly for output or restart only)
        ! IF (ASSOCIATED(field%q_vdf)) field%q_vdf(:,:,:) = vdf%atmo%Get_diagnostic_r3d('layer heating')

        ! field%cptgz(:,:,:) = vdf%atmo%Get_diagnostic_r3d('static energy')

        ! field%ts (:,:) = vdf%sfc%Get_diagnostic_r2d('sfc temperature')
        ! !
        ! field%lhflx (:,:) = vdf%sfc%Get_diagnostic_r2d('sfc latent heat flux')
        ! field%shflx (:,:) = vdf%sfc%Get_diagnostic_r2d('sfc sensible heat flux')
        ! field%u_stress (:,:) = vdf%sfc%Get_diagnostic_r2d('sfc zonal wind stress')
        ! field%v_stress (:,:) = vdf%sfc%Get_diagnostic_r2d('sfc mer. wind stress')
        ! !
        ! field%lhflx_tile (:,:,:) = vdf%sfc%Get_diagnostic_r3d('sfc latent heat flux, tile')
        ! field%shflx_tile (:,:,:) = vdf%sfc%Get_diagnostic_r3d('sfc sensible heat flux, tile')
        ! field%u_stress_tile (:,:,:) = vdf%sfc%Get_diagnostic_r3d('sfc zonal wind stress, tile')
        ! field%v_stress_tile (:,:,:) = vdf%sfc%Get_diagnostic_r3d('sfc mer. wind stress, tile')
        ! !
        ! field%cfm (:,:,:) = vdf%atmo%Get_diagnostic_r3d('exchange coefficient momentum')
        ! field%cfh (:,:,:) = vdf%atmo%Get_diagnostic_r3d('drag coefficient scalar')
        ! !
        ! field%z0m(:,:) = vdf%sfc%Get_diagnostic_r2d('roughness length momentum')
        ! field%z0h(:,:) = vdf%sfc%Get_diagnostic_r2d('roughness length heat')

      END IF
    END IF

    ! disassociate pointers
    NULLIFY(field)
    NULLIFY(tend)

    IF (l_init_or_restart) THEN
      l_init_or_restart = .FALSE.
    END IF

    IF (ltimer) CALL timer_stop(timer_tmx)

  END SUBROUTINE interface_aes_tmx

  !New vdf
  SUBROUTINE init_tmx(p_patch, dtime)

    USE mo_vdf,      ONLY: heat_type, momentum_type
    USE mo_tmx_field_class, ONLY: isfc_oce, isfc_ice, isfc_lnd

    USE mo_nonhydro_state,     ONLY: p_nh_state
    USE mo_nonhydro_types,     ONLY: t_nh_metrics, t_nh_diag
    USE mo_dynamics_config,    ONLY: nnow, nnow_rcf
    USE mo_physical_constants, ONLY: cpd, cvd, cvv

    USE mo_run_config,    ONLY: lmemman

    TYPE(t_patch), INTENT(inout), TARGET :: p_patch
    REAL(wp),      INTENT(in)    :: dtime

    TYPE(t_aes_phy_field),POINTER :: field
    TYPE(t_aes_phy_tend), POINTER :: tend
    TYPE(t_nh_metrics),   POINTER :: p_nh_metrics
    TYPE(t_nh_diag),      POINTER :: p_nh_diag

    TYPE(t_vdf), POINTER :: vdf

    TYPE(t_patch), POINTER :: patch
    INTEGER :: nlev, nlevp1, jg
    INTEGER :: rls, rle, jbs, jbe, jcs, jce, jb, jc
    INTEGER, ALLOCATABLE :: sfc_types(:)

    REAL(wp), POINTER :: dz_srf(:,:)

    TYPE(var_descriptor) :: var_desc

    LOGICAL           :: l_co2

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_tmx'

    IF (.NOT. lmemman) THEN
      CALL finish(routine, 'ERROR: Memory manager must be activated (lmemman=.TRUE.)')
    END IF

    patch => p_patch

    jg = patch%id

    nlev = patch%nlev
    nlevp1 = nlev + 1

    rls = grf_bdywidth_c + 1
    rle = min_rlcell_int

    jbs = patch%cells%start_block(rls)
    jbe = patch%cells%end_block  (rle)

    field     => prm_field(jg)
    tend      => prm_tend (jg)
    p_nh_metrics => p_nh_state(jg)%metrics
    p_nh_diag => p_nh_state(jg)%diag
    ! Question: use fields from AES field or from e.g. p_nh_state_lists(jg)%metrics p_nh_state_lists(jg)%diag? !!!!!!!!!!

    ! ALLOCATE(dz_srf(nproma,patch%nblks_c))

    l_co2 = (iqt <= ico2 .AND. ico2 <= ntracer)

    IF (ccycle_config(jg)%iccycle /= 0 .AND. .NOT. l_co2) THEN
      CALL finish(routine,'The C-cycle cannot be used without CO2 tracer (ico2<iqt or ntracer<ico2)')
    END IF

    ! The order of sfc_types must be consistent with how variables are defined in aes memory!
    ALLOCATE(sfc_types(0))
    IF (iwtr <= nsfc_type) sfc_types = [sfc_types, isfc_oce]
    IF (iice <= nsfc_type) sfc_types = [sfc_types, isfc_ice]
    IF (ilnd <= nsfc_type) THEN
#ifndef __NO_JSBACH__
      sfc_types = [sfc_types, isfc_lnd]
#else
      CALL finish(routine, 'JSBACH is not available, no land surface type defined')
#endif
    END IF

    vdf => new_vdf(patch, nproma, nlev=nlev, nsfc_tiles=nsfc_type, sfc_types=sfc_types, dt=dtime)
    __acc_attach(vdf)

    CALL vdf%atmo%Add_state(vdf%atmo%temp_idx,   diffusion_type=heat_type,     &
      & var_desc=var_descriptor('theta_v',jg,1,1,nnow(jg)), rank=3)
    CALL vdf%atmo%Add_state(vdf%atmo%uwind_idx,    diffusion_type=momentum_type, &
      & var_desc=var_descriptor('u',jg,1,1,-1), rank=3)
    CALL vdf%atmo%Add_state(vdf%atmo%vwind_idx,    diffusion_type=momentum_type, &
      & var_desc=var_descriptor('v',jg,1,1,-1), rank=3)
    CALL vdf%atmo%Add_state(vdf%atmo%wwind_idx,    diffusion_type=momentum_type, &
      & var_desc=var_descriptor('w',jg,1,1,nnew(jg)), rank=3)

    IF (l_co2) THEN
      CALL vdf%atmo%Add_state(vdf%atmo%tracer_idx, diffusion_type=heat_type, &
        & var_desc=var_descriptor('tracer',jg,1,1,nnow_rcf(jg)), rank=4, ref_pos=4, ref_idx=[iqv,iqc,iqi,ico2])
    ELSE
      CALL vdf%atmo%Add_state(vdf%atmo%tracer_idx, diffusion_type=heat_type, &
        & var_desc=var_descriptor('tracer',jg,1,1,nnow_rcf(jg)), rank=4, ref_pos=4, ref_idx=[iqv,iqc,iqi])
    END IF

    ! Bind variables to atmo config list

    CALL vdf%atmo%config%cpd%Assign_r0d(cpd)
    CALL vdf%atmo%config%cvd%Assign_r0d(cvd)
    CALL vdf%atmo%config%smag_constant%Assign_r0d(aes_vdf_config(jg)%smag_constant)
    CALL vdf%atmo%config%max_turb_scale%Assign_r0d(aes_vdf_config(jg)%max_turb_scale)
    CALL vdf%atmo%config%km_min%Assign_r0d(aes_vdf_config(jg)%km_min)
    CALL vdf%atmo%config%rturb_prandtl%Assign_r0d(aes_vdf_config(jg)%rturb_prandtl)
    CALL vdf%atmo%config%turb_prandtl%Assign_r0d(aes_vdf_config(jg)%turb_prandtl)
    CALL vdf%atmo%config%use_louis%Assign_l0d(aes_vdf_config(jg)%use_louis)
    CALL vdf%atmo%config%use_louis_land%Assign_l0d(aes_vdf_config(jg)%use_louis_land)
    CALL vdf%atmo%config%use_louis_ice%Assign_l0d(aes_vdf_config(jg)%use_louis_ice)
    CALL vdf%atmo%config%louis_constant_b%Assign_r0d(aes_vdf_config(jg)%louis_constant_b)
    CALL vdf%atmo%config%use_km_const%Assign_l0d(aes_vdf_config(jg)%use_km_const)
    CALL vdf%atmo%config%km_const%Assign_r0d(aes_vdf_config(jg)%km_const)
    CALL vdf%atmo%config%use_scale_turb_energy_flux%Assign_l0d(aes_vdf_config(jg)%use_scale_turb_energy_flux)
    CALL vdf%atmo%config%scale_turb_energy_flux%Assign_r0d(aes_vdf_config(jg)%scale_turb_energy_flux)
    CALL vdf%atmo%config%dtime%Assign_r0d(dtime)
    CALL vdf%atmo%config%solver_type%Assign_i0d(aes_vdf_config(jg)%solver_type)
    CALL vdf%atmo%config%energy_type%Assign_i0d(aes_vdf_config(jg)%energy_type)
    CALL vdf%atmo%config%dissipation_factor%Assign_r0d(aes_vdf_config(jg)%dissipation_factor)
    CALL vdf%atmo%config%l_co2%Assign_l0d(l_co2)

    ! Bind variables to atmo input list
    ! 3d
    CALL vdf%atmo%inputs%tracer_c%Update      ('tracer',             patch_id=jg, time_id=nnow_rcf(jg), ref_pos=4)
    CALL vdf%atmo%inputs%temp_c%Update        ('theta_v',            patch_id=jg, time_id=nnow(jg))
    CALL vdf%atmo%inputs%u_wind_c%Update      ('u',                  patch_id=jg)
    CALL vdf%atmo%inputs%v_wind_c%Update      ('v',                  patch_id=jg)
    CALL vdf%atmo%inputs%w_wind_ic%Update     ('w',                  patch_id=jg, time_id=nnow(jg))
    CALL vdf%atmo%inputs%temp_virt_c%Update   ('tempv',              patch_id=jg)
    CALL vdf%atmo%inputs%moist_mass_c%Update  ('airmass_new',        patch_id=jg)
    CALL vdf%atmo%inputs%cv_air_c%Update      ('cvair',              patch_id=jg)
    CALL vdf%atmo%inputs%z_c%Update           ('z_mc',               patch_id=jg)
    CALL vdf%atmo%inputs%z_ic%Update          ('z_ifc',              patch_id=jg)
    CALL vdf%atmo%inputs%pres_c%Update        ('pres',               patch_id=jg)
    CALL vdf%atmo%inputs%pres_ic%Update       ('pres_ifc',           patch_id=jg)
    CALL vdf%atmo%inputs%dz_c%Update          ('ddqz_z_full',        patch_id=jg)
    CALL vdf%atmo%inputs%dz_ic%Update         ('ddqz_z_half',        patch_id=jg)
    CALL vdf%atmo%inputs%inv_dz_c%Update      ('inv_ddqz_z_full',    patch_id=jg)
    CALL vdf%atmo%inputs%inv_dz_ic%Update     ('inv_ddqz_z_half',    patch_id=jg)
    CALL vdf%atmo%inputs%geo_height_c%Update  ('z_mc',               patch_id=jg)
    CALL vdf%atmo%inputs%geo_height_ic%Update ('z_ifc',              patch_id=jg)
    CALL vdf%atmo%inputs%geopot_agl_ic%Update ('geopot_agl_ifc',     patch_id=jg)
    CALL vdf%atmo%inputs%rho_c%Update         ('rho',                patch_id=jg, time_id=nnew(jg))
    ! CALL vdf%atmo%inputs%vn_e%Update          ('vn',                 patch_id=jg, time_id=nnew(jg))

    ! Bind diagnostic variables
    CALL vdf%atmo%diagnostics%ctgz                 %Update('cptgz',   patch_id=jg)
    CALL vdf%atmo%diagnostics%km                   %Update('cfm',     patch_id=jg)
    CALL vdf%atmo%diagnostics%kh                   %Update('cfh',     patch_id=jg)
    IF (ASSOCIATED(field%q_vdf)) THEN
      CALL vdf%atmo%diagnostics%heating            %Update('q_vdf',   patch_id=jg)
    END IF
    IF (ASSOCIATED(field%kedisp)) THEN
      CALL vdf%atmo%diagnostics%dissip_ke_vi       %Update('kedisp',  patch_id=jg)
    END IF
    IF (ASSOCIATED(field%cptgzvi)) THEN
      CALL vdf%atmo%diagnostics%ctgzvi             %Update('cptgzvi', patch_id=jg)
    END IF
    IF (ASSOCIATED(field%utmxvi)) THEN
      CALL vdf%atmo%diagnostics%int_energy_vi      %Update('utmxvi',  patch_id=jg)
    END IF
    IF (ASSOCIATED(tend%utmxvi)) THEN
      CALL vdf%atmo%diagnostics%int_energy_vi_tend %Update('tend_utmxvi', patch_id=jg)
    END IF

    !
    ! Surface
    !
    CALL vdf%sfc%Add_state(vdf%sfc%qsat_idx, dims=[nproma,patch%nblks_c,SIZE(sfc_types)], diffusion_type=heat_type)
    CALL vdf%sfc%Add_state(vdf%sfc%tsfc_idx, diffusion_type=heat_type, &
      &                     var_desc=var_descriptor('ts_tile', jg, 1, 1, -1), rank=3)

    ! Bind variables to sfc config list
    CALL vdf%sfc%config%cpd%Assign_r0d(cpd)
    CALL vdf%sfc%config%cvd%Assign_r0d(cvd)
    CALL vdf%sfc%config%cvv%Assign_r0d(cvv)
    CALL vdf%sfc%config%dtime%Assign_r0d(dtime)
    CALL vdf%sfc%config%min_sfc_wind%Assign_r0d(aes_vdf_config(jg)%min_sfc_wind)
    CALL vdf%sfc%config%wind_gustiness%Assign_r0d(aes_vdf_config(jg)%wind_g)
    CALL vdf%sfc%config%rough_m_oce%Assign_r0d(aes_vdf_config(jg)%z0m_oce)
    CALL vdf%sfc%config%rough_m_ice%Assign_r0d(aes_vdf_config(jg)%z0m_ice)
    CALL vdf%sfc%config%min_rough%Assign_r0d(aes_vdf_config(jg)%z0m_min)
    CALL vdf%sfc%config%fsl%Assign_r0d(aes_vdf_config(jg)%fsl)
    CALL vdf%sfc%config%nice_thickness_classes%Assign_i0d(field%kice)
    CALL vdf%sfc%config%l_co2%Assign_l0d(l_co2)

    ! Bind variables to sfc input list
    CALL vdf%sfc%inputs%qa%Update('tracer', time_id=nnow(jg), patch_id=jg, &
      &                           ref_pos=2, ref=nlev, ref2_pos=4, ref2=iqv)
    CALL vdf%sfc%inputs%ta%Update('theta_v', time_id=nnow(jg), patch_id=jg, ref_pos=2, ref=nlev)
    CALL vdf%sfc%inputs%tv%Update('tempv', patch_id=jg, ref_pos=2, ref=nlev)
    CALL vdf%sfc%inputs%ua%Update('u', patch_id=jg, ref_pos=2, ref=nlev)
    CALL vdf%sfc%inputs%va%Update('v', patch_id=jg, ref_pos=2, ref=nlev)
    CALL vdf%sfc%inputs%rho_atm%Update('rho', time_id=nnow(jg), patch_id=jg, ref_pos=2, ref=nlev)
    CALL vdf%sfc%inputs%pa%Update('pres', patch_id=jg, ref_pos=2, ref=nlev)
    CALL vdf%sfc%inputs%psfc%Update('pres_ifc', patch_id=jg, ref_pos=2, ref=nlevp1)
    ! CALL bind_variable(vdf%sfc%inputs%list%Search('surface pressure'), p_nh_diag%pres_sfc)
    CALL vdf%sfc%inputs%zf%Update('z_mc', patch_id=jg, ref_pos=2, ref=nlev)
    CALL vdf%sfc%inputs%zh%Update('z_ifc', patch_id=jg, ref_pos=2, ref=nlevp1)
    !
    CALL vdf%sfc%inputs%rsfl%Update('prlr', patch_id=jg)
    CALL vdf%sfc%inputs%ssfl%Update('prls', patch_id=jg)
    !
    CALL vdf%sfc%inputs%rlds%Update('rlds', patch_id=jg)
    CALL vdf%sfc%inputs%rsds%Update('rsds', patch_id=jg)
    CALL vdf%sfc%inputs%rvds_dir%Update('rvds_dir', patch_id=jg)
    CALL vdf%sfc%inputs%rnds_dir%Update('rnds_dir', patch_id=jg)
    CALL vdf%sfc%inputs%rpds_dir%Update('rpds_dir', patch_id=jg)
    CALL vdf%sfc%inputs%rvds_dif%Update('rvds_dif', patch_id=jg)
    CALL vdf%sfc%inputs%rnds_dif%Update('rnds_dif', patch_id=jg)
    CALL vdf%sfc%inputs%rpds_dif%Update('rpds_dif', patch_id=jg)
    !
    CALL vdf%sfc%inputs%tsfc_tile%Update('ts_tile', patch_id=jg)
    CALL vdf%sfc%inputs%fract_tile%Update('frac_tile', patch_id=jg)
    !
    ! TODO: lw surface emissivity should be tile-specific and, for land, should be returned from land model
    CALL vdf%sfc%inputs%emissivity%Update('emissivity', patch_id=jg)

    CALL vdf%sfc%inputs%ocean_u%Update('ocean_u', patch_id=jg)
    CALL vdf%sfc%inputs%ocean_v%Update('ocean_v', patch_id=jg)
    CALL vdf%sfc%inputs%ice_u%Update('ice_u', patch_id=jg)
    CALL vdf%sfc%inputs%ice_v%Update('ice_v', patch_id=jg)

    !
    CALL vdf%sfc%inputs%ice_thickness%Update('sit_icecl', patch_id=jg, ref_pos=2, ref=1)
    !
    CALL vdf%sfc%inputs%cosmu0%Update('cosmu0', patch_id=jg)
    !
    CALL vdf%sfc%inputs%co2flx_ant%Update('fco2ant', patch_id=jg)
    CALL vdf%sfc%diagnostics%co2flx_nat%Update('fco2nat', patch_id=jg)
    !
    CALL vdf%sfc%inputs%dz%Update('ddqz_z_half', patch_id=jg, ref_pos=2, ref=nlevp1)
    ! dz_srf(:,:) = 2._wp * (field%zh(:,nlev,:) - field%zh(:,nlevp1,:))
    ! dz_srf(:,:) = 2._wp * (p_nh_metrics%z_mc(:,nlev,:) - p_nh_metrics%z_ifc(:,nlevp1,:))
    ! CALL bind_variable(vdf%sfc%inputs%list%Search('reference height in surface layer times 2'), dz_srf)
    !
    CALL vdf%sfc%diagnostics%tsfc%Update('ts', patch_id=jg)
    CALL vdf%sfc%diagnostics%tsfc_rad%Update('ts_rad', patch_id=jg)
    CALL vdf%sfc%diagnostics%lwfl_net_tile%Update('rlns_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%swfl_net_tile%Update('rsns_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%lwfl_up%Update('rlus', patch_id=jg)
    CALL vdf%sfc%diagnostics%swfl_up%Update('rsus', patch_id=jg)
    !
    CALL vdf%sfc%diagnostics%evapotrans%Update('evspsbl', patch_id=jg)
    CALL vdf%sfc%diagnostics%lhfl%Update('hfls', patch_id=jg)
    CALL vdf%sfc%diagnostics%shfl%Update('hfss', patch_id=jg)
    CALL vdf%sfc%diagnostics%ustress%Update('tauu', patch_id=jg)
    CALL vdf%sfc%diagnostics%vstress%Update('tauv', patch_id=jg)
    CALL vdf%sfc%diagnostics%ufts%Update('ufts', patch_id=jg)
    CALL vdf%sfc%diagnostics%ufvs%Update('ufvs', patch_id=jg)
    CALL vdf%sfc%diagnostics%evapotrans_tile%Update('evspsbl_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%lhfl_tile%Update('hfls_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%shfl_tile%Update('hfss_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%ustress_tile%Update('tauu_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%vstress_tile%Update('tauv_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%co2flx_nat_tile%Update('co2_flux_tile', patch_id=jg)

    IF (ASSOCIATED(field%z0m)) &
      & CALL vdf%sfc%diagnostics%rough_m%Update('z0m', patch_id=jg)
    IF (ASSOCIATED(field%z0h)) &
      & CALL vdf%sfc%diagnostics%rough_h%Update('z0h', patch_id=jg)
    IF (ASSOCIATED(field%z0m_tile)) &
      & CALL vdf%sfc%diagnostics%rough_m_tile%Update('z0m_tile', patch_id=jg)
    IF (ASSOCIATED(field%z0h_tile)) &
      & CALL vdf%sfc%diagnostics%rough_h_tile%Update('z0h_tile', patch_id=jg)
    IF (ASSOCIATED(field%cfm_tile)) &
      & CALL vdf%sfc%diagnostics%km_tile%Update('cfm_tile', patch_id=jg)
    IF (ASSOCIATED(field%cfh_tile)) &
      & CALL vdf%sfc%diagnostics%kh_tile%Update('cfh_tile', patch_id=jg)
    !
    CALL vdf%sfc%diagnostics%albvisdir_tile%Update('albvisdir_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%albvisdif_tile%Update('albvisdif_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%albnirdir_tile%Update('albnirdir_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%albnirdif_tile%Update('albnirdif_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%albedo_tile%Update('albedo_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%albvisdir%Update('albvisdir', patch_id=jg)
    CALL vdf%sfc%diagnostics%albvisdif%Update('albvisdif', patch_id=jg)
    CALL vdf%sfc%diagnostics%albnirdir%Update('albnirdir', patch_id=jg)
    CALL vdf%sfc%diagnostics%albnirdif%Update('albnirdif', patch_id=jg)
    CALL vdf%sfc%diagnostics%albedo%Update('albedo', patch_id=jg)

    CALL vdf%sfc%diagnostics%q_ice_top%Update('qtop_icecl', patch_id=jg, ref_pos=2, ref=1)
    CALL vdf%sfc%diagnostics%q_ice_bot%Update('qbot_icecl', patch_id=jg, ref_pos=2, ref=1)
    CALL vdf%sfc%diagnostics%snow_thickness%Update('hs_icecl', patch_id=jg, ref_pos=2, ref=1)

    CALL vdf%sfc%diagnostics%t2m%Update('tas', patch_id=jg)
    CALL vdf%sfc%diagnostics%t2m_tile%Update('tas_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%hus2m%Update('qv2m', patch_id=jg)
    CALL vdf%sfc%diagnostics%hus2m_tile%Update('qv2m_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%dew2m%Update('dew2', patch_id=jg)
    CALL vdf%sfc%diagnostics%dew2m_tile%Update('dew2_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%u10m%Update('uas', patch_id=jg)
    CALL vdf%sfc%diagnostics%u10m_tile%Update('uas_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%v10m%Update('vas', patch_id=jg)
    CALL vdf%sfc%diagnostics%v10m_tile%Update('vas_tile', patch_id=jg)
    CALL vdf%sfc%diagnostics%wind10m%Update('sfcwind', patch_id=jg)
    CALL vdf%sfc%diagnostics%wind10m_tile%Update('sfcwind_tile', patch_id=jg)

    vdf_dom(jg)%p => vdf

    CALL vdf%Init()

  END SUBROUTINE init_tmx

END MODULE mo_interface_aes_tmx
