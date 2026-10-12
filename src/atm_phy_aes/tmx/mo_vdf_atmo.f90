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

MODULE mo_vdf_atmo

  USE mo_kind,              ONLY: wp, vp, sp
  USE mo_exception,         ONLY: message, finish
  USE mtime,                ONLY: t_datetime => datetime
  USE mo_timer,             ONLY: timer_start, timer_stop, ltimer
  USE mo_tmx_process_class, ONLY: t_tmx_process
  USE mo_tmx_field_class,   ONLY: t_tmx_field, t_domain
  USE mo_tmx_var,           ONLY: t_tmx_var
  USE mo_vdf_atmo_memory,   ONLY: t_vdf_atmo_config, t_vdf_atmo_inputs, t_vdf_atmo_diags, &
    &                             build_vdf_atmo_config, build_vdf_atmo_inputs, build_vdf_atmo_diags, &
    &                             t_vel_grad_tensor
  USE mo_tmx_smagorinsky,   ONLY: Smagorinsky_init, Smagorinsky_model
  USE mo_model_domain,      ONLY: t_patch
  USE mo_intp_data_strc,    ONLY: t_int_state, p_int_state
  USE mo_nonhydro_types,    ONLY: t_nh_metrics
  USE mo_nonhydro_state,    ONLY: p_nh_state
  USE mo_run_config,        ONLY: ntracer, iqv, iqc, iqi, iqr, iqs, iqg, iqnc, iqni, iqt, ico2
  USE mo_aes_sfc_indices,   ONLY: nsfc_type
  USE mo_physical_constants,ONLY: grav, rd, cpd, cpv, cvd, rd_o_cpd, &
    &                             p0ref, rgrav
  USE mo_impl_constants,    ONLY: min_rlcell, min_rledge_int, min_rlcell_int, min_rlvert_int
  USE mo_nh_vert_interp_les,ONLY: brunt_vaisala_freq, vert_intp_full2half_cell_3d
  USE mo_intp,              ONLY: cells2verts_scalar, cells2edges_scalar
  USE mo_sync,              ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array,     &
  &                             sync_patch_array_mult
  USE mo_loopindices,       ONLY: get_indices_e, get_indices_c
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_intp_rbf,          ONLY: rbf_vec_interpol_vertex, rbf_vec_interpol_edge
  USE mo_fortran_tools,     ONLY: init

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_vdf_atmo, t_vdf_atmo_config, &
    & prepare_diffusion_matrix ! , compute_temp_from_static_energy

  !Parameters for surface layer parameterizations: From Zeng_etal 1997 J. Clim
  REAL(wp), PARAMETER :: bsm = 5.0_wp  !Businger Stable Momentum
  REAL(wp), PARAMETER :: bum = 16._wp  !Businger Untable Momentum
  REAL(wp), PARAMETER :: bsh = 5.0_wp  !Businger Stable Heat
  REAL(wp), PARAMETER :: buh = 16._wp  !Businger Untable Heat

  INTEGER, PARAMETER :: MAX_NO_STATES = 5 !< Maximum number of states in t_vdf_atmo

  TYPE, EXTENDS(t_tmx_process) :: t_vdf_atmo
    TYPE(t_vdf_atmo_config), POINTER :: config => NULL()
    TYPE(t_vdf_atmo_inputs), POINTER :: inputs => NULL()
    TYPE(t_vdf_atmo_diags),  POINTER :: diagnostics => NULL()
    !
    ! Supported diffusion variables
    !
    INTEGER :: temp_idx   = 1
    INTEGER :: tracer_idx = 2
    INTEGER :: uwind_idx  = 3
    INTEGER :: vwind_idx  = 4
    INTEGER :: wwind_idx  = 5
  CONTAINS
    PROCEDURE :: Init => Init_vdf_atmo
    PROCEDURE :: Compute
    PROCEDURE :: Compute_diagnostics
    PROCEDURE :: temp_to_energy
    PROCEDURE :: energy_to_temp
    PROCEDURE :: compute_flux_x
    PROCEDURE :: Update_diagnostics
  END TYPE t_vdf_atmo

  INTERFACE prepare_diffusion_matrix
    MODULE PROCEDURE prepare_diffusion_matrix_wp
#ifdef __MIXED_PRECISION
    MODULE PROCEDURE prepare_diffusion_matrix_mp
#endif
  END INTERFACE prepare_diffusion_matrix

  INTERFACE t_vdf_atmo
    MODULE PROCEDURE t_vdf_atmo_construct
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_vdf_atmo'

CONTAINS

  !-----------------------------------------------------------------
  ! Function: t_vdf_atmo_construct
  ! Purpose: Constructs and initializes an instance of the t_vdf_atmo type.
  ! Inputs:
  !   - name: A character string representing the name of the process.
  !   - domain: A pointer to the t_domain type representing the model domain.
  ! Outputs:
  !   - result: A pointer to the newly constructed t_vdf_atmo object.
  ! Behavior:
  !   - Allocates memory for the t_vdf_atmo object.
  !   - Initializes the object by calling its parent class's Init_process method.
  !   - Sets up configuration and diagnostic structures.
  !-----------------------------------------------------------------
  FUNCTION t_vdf_atmo_construct(name, dt, domain) RESULT(result)

    CHARACTER(len=*), INTENT(in) :: name
    REAL(wp),         INTENT(in) :: dt
    TYPE(t_domain),      POINTER :: domain
    TYPE(t_vdf_atmo),    POINTER :: result

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_vdf_atmo_construct'

    CALL message(routine, '')

    ALLOCATE(t_vdf_atmo::result)
    result%max_no_states = MAX_NO_STATES
    !$ACC ENTER DATA COPYIN(result)
    ! Call Init of abstract parent class
    CALL result%Init_process(dt=dt, name=name, domain=domain)
    __acc_attach(result%domain)

    ! Initialize structure for config variables (first scan)
    ALLOCATE(t_vdf_atmo_config :: result%config)
    CALL build_vdf_atmo_config(result%config, result%domain)

    ! Initialize structure for input variables (first scan)
    ALLOCATE(t_vdf_atmo_inputs :: result%inputs)
    CALL build_vdf_atmo_inputs(result%inputs, result%domain)

    ! Initialize structure for diagnostic variables (first scan)
    ALLOCATE(t_vdf_atmo_diags :: result%diagnostics)
    CALL build_vdf_atmo_diags(result%diagnostics, result%domain)

  END FUNCTION t_vdf_atmo_construct
  !
  !============================================================================
  !
  SUBROUTINE Init_vdf_atmo(this)
    CLASS(t_vdf_atmo), INTENT(inout), TARGET :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init'

    CALL message(routine, '')

    ! Initialize structure for config variables (second scan)
    CALL build_vdf_atmo_config(this%config, this%domain)
    !$ACC ENTER DATA COPYIN(this%config)

    ! Initialize structure for input variables (second scan)
    CALL build_vdf_atmo_inputs(this%inputs, this%domain)
    !$ACC ENTER DATA COPYIN(this%inputs)

    ! Initialize structure for diagnostic variables (second scan)
    ! ACC copyin moved inside build_vdf_atmo_diags to avoid
    ! "partially present on device" OpenACC error
    CALL build_vdf_atmo_diags(this%diagnostics, this%domain)

  END SUBROUTINE Init_vdf_atmo
  !
  !============================================================================
  !
  SUBROUTINE Compute(this, datetime)

    CLASS(t_vdf_atmo), INTENT(inout), TARGET :: this
    TYPE(t_datetime), OPTIONAL, INTENT(in),   POINTER :: datetime     !< date and time at beginning of time step

    CHARACTER(len=*), PARAMETER :: routine = modname//':Compute'

    CALL message(routine, '')

  END SUBROUTINE Compute
  !
  !============================================================================
  !
  SUBROUTINE Compute_diagnostics(this, datetime)

    CLASS(t_vdf_atmo), INTENT(inout), TARGET :: this
    TYPE(t_datetime), OPTIONAL, INTENT(in), POINTER :: datetime

    TYPE(t_vdf_atmo_config),      POINTER :: config
    TYPE(t_vdf_atmo_inputs),      POINTER :: inputs
    TYPE(t_vdf_atmo_diags),       POINTER :: diags
    TYPE(t_int_state),            POINTER :: p_int         !< interpolation state
    TYPE(t_nh_metrics),           POINTER :: p_nh_metrics
    TYPE(t_domain),               POINTER :: domain
    TYPE(t_patch),                POINTER :: patch

    ! Pointers to configuration variables
    REAL(wp), POINTER :: &
      & cpd, cvd, rturb_prandtl, louis_constant_b, km_min, km_const
    LOGICAL, POINTER :: &
      & use_louis, use_louis_land, use_louis_ice, use_km_const

    ! Pointers to input variables
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & zf, zh, dz, inv_dzh, &
      & ptm1, pum1, pvm1, pwp1, &
      & ptvm1, rho, mair, cvair, papm1, paphm1, &
      & pqm1, pxlm1, pxim1, pxrm1, pxsm1, pxgm1 !, vn

    REAL(vp), POINTER, DIMENSION(:,:,:) :: &
      & dzh, inv_dzf

    ! Pointers to diagnostic variables
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & ghf, ctgz, div_c, theta_v, pprfac, km, kh, km_c, heating, &  ! 3D full level cell diagnostics
      & rho_ic, bruvais, mech_prod, km_ic, kh_ic, mix_len_sq, &  ! 3D half level cell diagnostics
      & vn, shear, div_stress, &  ! 3D full level edge diagnostics
      & vn_ie, vt_ie, w_ie, km_ie, &  ! 3D half level edge diagnostics
      & u_vert, v_vert, w_vert, km_iv ! 3D vertex diagnostics
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & louis_factor, &
      & fract_land, &
      & fract_ice

    INTEGER :: jg
    INTEGER :: rl_start, rl_end
    INTEGER :: istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Compute_diagnostics'

    ! Initialize Smagorinsky model
    IF (this%is_initial_time) THEN
      CALL Smagorinsky_init(this%domain, this%config, this%inputs, this%diagnostics)
    END IF

    ! Get pointer to structure for config variables
    config => this%config
    __acc_attach(config)

    ! Get pointer to structure for input variables
    inputs => this%inputs
    __acc_attach(inputs)

    ! Get pointer to structure for diagnostic variables
    diags => this%diagnostics
    __acc_attach(diags)

    ! Prepare variables and pointers
    domain => this%domain
    patch  => domain%patch

    jg = patch%id
    rl_start = 3
    rl_end = min_rlcell_int

    p_int         => p_int_state(jg)
    p_nh_metrics  => p_nh_state(jg)%metrics

    ! Get pointers to configuration variables
    cpd              => config%cpd%Get_ptr_r0d()
    cvd              => config%cvd%Get_ptr_r0d()
    rturb_prandtl    => config%rturb_prandtl%Get_ptr_r0d()
    use_louis        => config%use_louis%Get_ptr_l0d()
    use_louis_land   => config%use_louis_land%Get_ptr_l0d()
    use_louis_ice    => config%use_louis_ice%Get_ptr_l0d()
    louis_constant_b => config%louis_constant_b%Get_ptr_r0d()
    km_min           => config%km_min%Get_ptr_r0d()
    use_km_const     => config%use_km_const%Get_ptr_l0d()
    km_const         => config%km_const%Get_ptr_r0d()

    ! Get pointers to input variables
    zf            => inputs%geo_height_c%Get_ptr_r3d()
    zh            => inputs%geo_height_ic%Get_ptr_r3d()
    dz            => inputs%dz_c%Get_ptr_r3d()
    inv_dzf       => inputs%inv_dz_c%Get_ptr_v3d()
    dzh           => inputs%dz_ic%Get_ptr_v3d()
    inv_dzh       => inputs%inv_dz_ic%Get_ptr_r3d()
    ptm1          => inputs%temp_c%Get_ptr_r3d()
    ptvm1         => inputs%temp_virt_c%Get_ptr_r3d()
    pum1          => inputs%u_wind_c%Get_ptr_r3d()
    pvm1          => inputs%v_wind_c%Get_ptr_r3d()
    pwp1          => inputs%w_wind_ic%Get_ptr_r3d()
    rho           => inputs%rho_c%Get_ptr_r3d()
    mair          => inputs%moist_mass_c%Get_ptr_r3d()
    cvair         => inputs%cv_air_c%Get_ptr_r3d()
    pqm1          => inputs%tracer_c%Get_ptr_r3d(ref=iqv)
    pxlm1         => inputs%tracer_c%Get_ptr_r3d(ref=iqc)
    pxim1         => inputs%tracer_c%Get_ptr_r3d(ref=iqi)
    pxrm1         => inputs%tracer_c%Get_ptr_r3d(ref=iqr)
    pxsm1         => inputs%tracer_c%Get_ptr_r3d(ref=iqs)
    pxgm1         => inputs%tracer_c%Get_ptr_r3d(ref=iqg)
    papm1         => inputs%pres_c%Get_ptr_r3d()
    paphm1        => inputs%pres_ic%Get_ptr_r3d()
    ! vn            => inputs%vn_e%Get_ptr_r3d()
    fract_land    => inputs%fract_land%Get_ptr_r2d()
    fract_ice     => inputs%fract_ice%Get_ptr_r2d()

    ! Get pointers to diagnostic variables from the diagnostics structure
    ! 3D full level cell diagnostics
    ghf           => diags%ghf%Get_ptr_r3d()
    ctgz          => diags%ctgz%Get_ptr_r3d()
    div_c         => diags%div_c%Get_ptr_r3d()
    theta_v       => diags%theta_v%Get_ptr_r3d()
    pprfac        => diags%pprfac%Get_ptr_r3d()
    km            => diags%km%Get_ptr_r3d()
    kh            => diags%kh%Get_ptr_r3d()
    km_c          => diags%km_c%Get_ptr_r3d()
    heating       => diags%heating%Get_ptr_r3d()

    ! 3D half level cell diagnostics
    rho_ic        => diags%rho_ic%Get_ptr_r3d()
    bruvais       => diags%bruvais%Get_ptr_r3d()
    mech_prod     => diags%mech_prod%Get_ptr_r3d()
    km_ic         => diags%km_ic%Get_ptr_r3d()
    kh_ic         => diags%kh_ic%Get_ptr_r3d()
    mix_len_sq    => diags%mix_len_sq%Get_ptr_r3d()

    ! 3D full level edge diagnostics
    vn            => diags%vn%Get_ptr_r3d()
    shear         => diags%shear%Get_ptr_r3d()
    div_stress    => diags%div_stress%Get_ptr_r3d()

    ! 3D half level edge diagnostics
    vn_ie         => diags%vn_ie%Get_ptr_r3d()
    vt_ie         => diags%vt_ie%Get_ptr_r3d()
    w_ie          => diags%w_ie%Get_ptr_r3d()
    km_ie         => diags%km_ie%Get_ptr_r3d()

    ! 3D vertex diagnostics
    u_vert        => diags%u_vert%Get_ptr_r3d()
    v_vert        => diags%v_vert%Get_ptr_r3d()
    w_vert        => diags%w_vert%Get_ptr_r3d()
    km_iv         => diags%km_iv%Get_ptr_r3d()

    ! 2D diagnostics
    louis_factor  => diags%louis_factor%Get_ptr_r2d()

    IF (ltimer) CALL timer_start(this%timer_diagnostics)

    !----------------------------------------------------------------------------
    ! Get static energy
    !----------------------------------------------------------------------------
    CALL compute_geopotential_height_above_ground(domain, zf, zh, ghf)

    CALL compute_static_energy(domain, cpd, ptm1, ghf, ctgz)

    !----------------------------------------------------------------------------
    ! Get virtual potential temperature
    !----------------------------------------------------------------------------
    CALL get_virtual_potential_temperature(patch, ptvm1, papm1, theta_v, &
                                           rl_start, rl_end)

    !Get rho at interfaces
    CALL vert_intp_full2half_cell_3d(patch, p_nh_metrics, rho, rho_ic, &
                                     2, min_rlcell_int-2, lacc=.TRUE.)

    !----------------------------------------------------------------------------
    ! Get Brunt-Vaisala frequency
    !----------------------------------------------------------------------------
    CALL brunt_vaisala_freq(patch, p_nh_metrics, domain%nproma, theta_v, bruvais, &
                            opt_rlstart=3, lacc=.TRUE.)

    !----------------------------------------------------------------------------
    ! Compute velocities normal to edges
    !----------------------------------------------------------------------------
    CALL sync_patch_array(SYNC_C, patch, pum1, lacc=.TRUE.)
    CALL sync_patch_array(SYNC_C, patch, pvm1, lacc=.TRUE.)

    CALL compute_normal_velocity_edge(pum1, pvm1, patch, p_int, &
                                      grf_bdywidth_e+1, min_rledge_int, vn)

    CALL sync_patch_array(SYNC_E, patch, vn, lacc=.TRUE.)

    !----------------------------------------------------------------------------
    ! Interpolate velocities at required locations to compute velocity
    ! gradient tensor and turbulent exchange coefficients.
    ! ->  assumes that prognostic values are all synced, while diagnostic
    !     values might not
    !----------------------------------------------------------------------------
!$OMP PARALLEL
    CALL init(u_vert, lacc=.TRUE.)
    CALL init(v_vert, lacc=.TRUE.)
    CALL init(w_vert, lacc=.TRUE.)
!$OMP END PARALLEL

    CALL cells2verts_scalar(pwp1, patch, p_int%cells_aw_verts, w_vert, &
                            lacc=.TRUE., opt_rlend=min_rlvert_int, opt_acc_async=.TRUE.)

    CALL cells2edges_scalar(pwp1, patch, p_int%c_lin_e, w_ie, lacc=.TRUE., &
                            opt_rlend=min_rledge_int-2)

    ! RBF reconstruction of velocity at vertices: include halos
    CALL rbf_vec_interpol_vertex(vn, patch, p_int, u_vert, v_vert, &
                                 lacc=.TRUE., opt_rlend=min_rlvert_int)

    !sync them
    CALL sync_patch_array_mult(SYNC_V, patch, 3, lacc=.TRUE., f3din1=w_vert, f3din2=u_vert, f3din3=v_vert)

    !Get vn at interfaces and then get vt at interfaces
    !Boundary values are extrapolated like dynamics although
    !they are not required in current implementation
    CALL interpolate_normal_velocity_edge_interface(vn, patch, p_nh_metrics, 2, &
                                                    min_rledge_int-3, vn_ie)

    CALL rbf_vec_interpol_edge(vn_ie, patch, p_int, vt_ie, lacc=.TRUE., &
                               opt_rlstart=3, opt_rlend=min_rledge_int-2)

    !----------------------------------------------------------------------------
    ! Compute velocity gradient tensor
    !----------------------------------------------------------------------------
    CALL compute_velocity_gradient_tensor(u_vert, v_vert, w_vert,         &
                                          pwp1, vn_ie, vt_ie, w_ie,       &
                                          patch, p_nh_metrics,            &
                                          4, min_rledge_int-2,            &
                                          diags%vel_grad_e)

    !----------------------------------------------------------------------------
    ! Compute strain rate and interpolate to required positions
    !----------------------------------------------------------------------------
    CALL compute_shear(diags%vel_grad_e, patch, 4, min_rledge_int-2,        &
                       shear, div_stress)

    !Interpolate mech production term from mid level edge to interface level cell
    !except top and bottom boundaries
    CALL get_horizontal_divergence_strain_rate_cell(div_stress, patch, p_int, &
                                                    grf_bdywidth_c+1, min_rlcell_int-1, &
                                                    div_c)

    ! Interpolate mech. production term from mid level edge to interface level cell: mech_prod = 2 * |S|^2
    ! except top and bottom boundaries
    CALL interpolate_rate_of_strain_full2half_edge2cell(shear, patch, p_nh_metrics, p_int, &
                                                        3, min_rlcell_int-1, mech_prod)

    !----------------------------------------------------------------------------
    ! Compute turbulent exchange coefficient Km & Kh according to classical
    ! Smagorinsky model with stability correction (Lilly 1962) at interface
    ! cell centers
    !----------------------------------------------------------------------------
!$OMP PARALLEL
    CALL init(km_iv, lacc=.TRUE.)
    CALL init(km_c,  lacc=.TRUE.)
    CALL init(km_ie, lacc=.TRUE.)
    CALL init(kh_ic, lacc=.TRUE.)
    CALL init(km_ic, lacc=.TRUE.)
!$OMP END PARALLEL

    IF (.NOT. use_km_const) THEN
      CALL Smagorinsky_model(domain, mech_prod, bruvais, rho_ic,                   &
                             mix_len_sq, rturb_prandtl, use_louis, use_louis_land, use_louis_ice, &
                             louis_constant_b, louis_factor, fract_land, fract_ice, &
                             km_ic, kh_ic)
    ELSE
      CALL Assign_constant_eddy_viscosity(domain,  rho_ic, km_const, &
                                          rturb_prandtl,             &
                                          km_ic, kh_ic)
    END IF

    !----------------------------------------------------------------------------
    ! Interpolate turbulent exchange coefficients to required positions for
    ! flux calculation.
    ! -> halos also computed since they are used in diffusion later
    !----------------------------------------------------------------------------
    ! visc at cell center
    CALL interpolate_eddy_viscosity2cell(km_ic, km_min, patch, grf_bdywidth_c, &
                                         min_rlcell_int-1, km_c)

    ! visc at vertices
    CALL interpolate_eddy_viscosity2half_vertex(km_ic, km_min, patch, p_int, &
                                                km_iv)

    ! Now calculate visc at half levels at edge
    CALL interpolate_eddy_viscosity2half_edge(km_ic, km_min, patch, p_int, &
                                              km_ie)

    IF (ltimer) CALL timer_stop(this%timer_diagnostics)

  END SUBROUTINE Compute_diagnostics
  !============================================================================
  !
  !============================================================================
  !
  SUBROUTINE Update_diagnostics(this)

    CLASS(t_vdf_atmo), INTENT(inout), TARGET :: this

    TYPE(t_vdf_atmo_config),      POINTER :: config
    TYPE(t_vdf_atmo_inputs),      POINTER :: inputs
    TYPE(t_vdf_atmo_diags),       POINTER :: diags

    TYPE(t_tmx_field), POINTER :: field

    INTEGER :: jb, jk, jc

    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & state_ta, state_qv, state_qc, state_qi, &
      & tend_ta, tend_qv, tend_qc, tend_qi, &
      & new_state_ta, new_state_qv, new_state_qc, new_state_qi

    ! Pointers to configuration variables
    REAL(wp), POINTER :: &
      & cpd, dtime
    ! Pointers to input variables
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & ghf, dz, rho, pxrm1, pxsm1, pxgm1
    ! Pointers to diagnostic variables
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & ctgz, dissip_ke
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & ctgzvi, dissip_ke_vi, int_energy_vi, int_energy_vi_tend
    REAL(wp), DIMENSION(this%domain%nproma,this%domain%nblks_c) :: &
      & int_energy_vi_old

    CHARACTER(len=*), PARAMETER :: routine = modname//':Update_diagnostics'

    ! CALL message(routine, 'start')

    config => this%config
    inputs => this%inputs
    diags  => this%diagnostics

    state_ta     => this%states    (this%temp_idx)%p%Get_ptr_r3d() ! old state
    tend_ta      => this%tendencies(this%temp_idx)%p%Get_ptr_r3d() ! tendency
    new_state_ta => this%new_states(this%temp_idx)%p%Get_ptr_r3d() ! new state

    field => this%states(this%tracer_idx)%p
    state_qv     => this%states    (this%tracer_idx)%p%Get_ptr_r3d(ref=field%ref_idx(1)) ! water vapor old state
    tend_qv      => this%tendencies(this%tracer_idx)%p%Get_ptr_r3d(ref=field%ref_idx(1)) ! water vapor tendency
    new_state_qv => this%new_states(this%tracer_idx)%p%Get_ptr_r3d(ref=field%ref_idx(1)) ! water vapor new state
    state_qc     => this%states    (this%tracer_idx)%p%Get_ptr_r3d(ref=field%ref_idx(2)) ! cloud water old state
    tend_qc      => this%tendencies(this%tracer_idx)%p%Get_ptr_r3d(ref=field%ref_idx(2)) ! cloud water tendency
    new_state_qc => this%new_states(this%tracer_idx)%p%Get_ptr_r3d(ref=field%ref_idx(2)) ! cloud water new state
    state_qi     => this%states    (this%tracer_idx)%p%Get_ptr_r3d(ref=field%ref_idx(3)) ! cloud ice old state
    tend_qi      => this%tendencies(this%tracer_idx)%p%Get_ptr_r3d(ref=field%ref_idx(3)) ! cloud ice tendency
    new_state_qi => this%new_states(this%tracer_idx)%p%Get_ptr_r3d(ref=field%ref_idx(3)) ! cloud ice new state

    ! Get pointers to configuration variables
    cpd   => config%cpd%Get_ptr_r0d()
    dtime => config%dtime%Get_ptr_r0d()

    ! Get pointers to input variables
    dz    => inputs%dz_c%Get_ptr_r3d()
    rho   => inputs%rho_c%Get_ptr_r3d()
    pxrm1 => inputs%tracer_c%Get_ptr_r3d(ref=iqr)
    pxsm1 => inputs%tracer_c%Get_ptr_r3d(ref=iqs)
    pxgm1 => inputs%tracer_c%Get_ptr_r3d(ref=iqg)

    ! Get pointers to diagnostic variables
    ghf                => diags%ghf%Get_ptr_r3d()
    ctgz               => diags%ctgz%Get_ptr_r3d()
    ctgzvi             => diags%ctgzvi%Get_ptr_r2d()
    dissip_ke          => diags%dissip_ke%Get_ptr_r3d()
    dissip_ke_vi       => diags%dissip_ke_vi%Get_ptr_r2d()
    int_energy_vi      => diags%int_energy_vi%Get_ptr_r2d()
    int_energy_vi_tend => diags%int_energy_vi_tend%Get_ptr_r2d()

    IF (ltimer) CALL timer_start(this%timer_diagnostics)

    ASSOCIATE( &
      domain => this%domain &
      & )

    !$ACC DATA CREATE(int_energy_vi_old)

    CALL compute_static_energy( &
      & domain, &
      & cpd, new_state_ta(:,:,:), ghf(:,:,:), &
      & ctgz(:,:,:) &
      & )

    ! Vertical integrals

    ! TODO: include hydrometeors from microphysics?
    CALL compute_internal_energy_vi( &
      & domain, &
      & rho(:,:,:), dz(:,:,:), pxrm1(:,:,:), pxsm1(:,:,:), pxgm1(:,:,:), &
      & state_ta(:,:,:), state_qv(:,:,:), state_qc(:,:,:), state_qi(:,:,:), &
      & int_energy_vi_old(:,:) &
      & )
    CALL compute_internal_energy_vi( &
      & domain, &
      & rho(:,:,:), dz(:,:,:), pxrm1(:,:,:), pxsm1(:,:,:), pxgm1(:,:,:), &
      & new_state_ta(:,:,:), new_state_qv(:,:,:), new_state_qc(:,:,:), new_state_qi(:,:,:), &
      & int_energy_vi(:,:) &
      & )

!$OMP PARALLEL
      CALL init(ctgzvi, lacc=.TRUE.)
      CALL init(dissip_ke_vi, lacc=.TRUE.)
      CALL init(int_energy_vi_tend, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = domain%i_startidx_c(jb),domain%i_endidx_c(jb)

          ! static energy (reference air path)
          ctgzvi(jc,jb) = ctgzvi(jc,jb) + ctgz(jc,jk,jb) * rho(jc,jk,jb) * dz(jc,jk,jb)
          ! kinetic energy dissipation
          dissip_ke_vi(jc,jb) = dissip_ke_vi(jc,jb) + dissip_ke(jc,jk,jb)

          ! internal energy tendency
          int_energy_vi_tend(jc,jb) = (int_energy_vi(jc,jb) - int_energy_vi_old(jc,jb)) / dtime

        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)
    !$ACC END DATA

    END ASSOCIATE

    IF (this%is_initial_time) this%is_initial_time = .FALSE.

    IF (ltimer) CALL timer_stop(this%timer_diagnostics)

    ! CALL message(routine, 'end')

  END SUBROUTINE Update_diagnostics
  !
  !============================================================================
  !
  SUBROUTINE temp_to_energy(this, temperature, energy, use_new_moisture_state)

    CLASS(t_vdf_atmo), INTENT(in) :: this
    LOGICAL,  INTENT(in), OPTIONAL :: use_new_moisture_state
    REAL(wp), INTENT(in) :: temperature(:,:,:)
    REAL(wp), INTENT(out) :: energy(:,:,:)

    INTEGER :: jb, jk, jc
    INTEGER,  POINTER :: energy_type
    LOGICAL :: use_updated_moisture
    REAL(wp), POINTER :: cpd
    REAL(wp), POINTER :: geo_height(:,:,:)
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & qr, qs, qg, qv, qc, qi
    TYPE(t_tmx_field), POINTER :: field

    CHARACTER(len=*), PARAMETER :: routine = modname//':temp_to_energy'

    energy_type => this%config%energy_type%Get_ptr_i0d()

    ! No effect for this%energy_type=1
    use_updated_moisture = .TRUE.
    IF (PRESENT(use_new_moisture_state)) use_updated_moisture = use_new_moisture_state

    geo_height => this%diagnostics%ghf%Get_ptr_r3d()

    SELECT CASE(energy_type)
    CASE (1)
      cpd => this%config%cpd%Get_ptr_r0d()
      CALL compute_static_energy(this%domain, cpd, temperature, geo_height, energy)
    CASE (2)
      qr  => this%inputs%tracer_c%Get_ptr_r3d(ref=iqr)
      qs  => this%inputs%tracer_c%Get_ptr_r3d(ref=iqs)
      qg  => this%inputs%tracer_c%Get_ptr_r3d(ref=iqg)

      field => this%states(this%tracer_idx)%p
      IF (use_updated_moisture) THEN
        qv  => this%new_states(this%tracer_idx)%p%get_ptr_r3d(ref=field%ref_idx(1)) ! water vapor
        qc  => this%new_states(this%tracer_idx)%p%get_ptr_r3d(ref=field%ref_idx(2)) ! cloud water
        qi  => this%new_states(this%tracer_idx)%p%get_ptr_r3d(ref=field%ref_idx(3)) ! cloud ice
      ELSE
        qv  => this%states(this%tracer_idx)%p%get_ptr_r3d(ref=field%ref_idx(1)) ! water vapor
        qc  => this%states(this%tracer_idx)%p%get_ptr_r3d(ref=field%ref_idx(2)) ! cloud water
        qi  => this%states(this%tracer_idx)%p%get_ptr_r3d(ref=field%ref_idx(3)) ! cloud ice
      END IF

      CALL compute_internal_energy( &
        & this%domain, &
        & geo_height,  &
        & qr, qs, qg,  &
        & temperature, &
        & qv, qc, qi,  &
        & energy       &
        & )
    END SELECT

  END SUBROUTINE temp_to_energy
  !
  !============================================================================
  !
  SUBROUTINE energy_to_temp(this, energy, temperature, use_new_moisture_state)

    CLASS(t_vdf_atmo), INTENT(in) :: this
    LOGICAL,  INTENT(in), OPTIONAL :: use_new_moisture_state
    REAL(wp), INTENT(in)  :: energy(:,:,:)
    REAL(wp), INTENT(out) :: temperature(:,:,:)

    INTEGER :: jb, jk, jc
    INTEGER,  POINTER :: energy_type
    LOGICAL :: use_updated_moisture
    REAL(wp), POINTER :: cpd
    REAL(wp), POINTER :: geo_height(:,:,:)
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & qr, qs, qg, qv, qc, qi
    TYPE(t_tmx_field), POINTER :: field

    CHARACTER(len=*), PARAMETER :: routine = modname//':energy_to_temp'

    energy_type => this%config%energy_type%Get_ptr_i0d()

    ! No effect for this%energy_type=1
    use_updated_moisture = .TRUE.
    IF (PRESENT(use_new_moisture_state)) use_updated_moisture = use_new_moisture_state

    geo_height => this%diagnostics%ghf%Get_ptr_r3d()

    SELECT CASE(energy_type)
    CASE (1)
      cpd => this%config%cpd%Get_ptr_r0d()
      CALL compute_temp_from_static_energy(this%domain, cpd, energy, geo_height, temperature)
    CASE (2)
      qr  => this%inputs%tracer_c%Get_ptr_r3d(ref=iqr)
      qs  => this%inputs%tracer_c%Get_ptr_r3d(ref=iqs)
      qg  => this%inputs%tracer_c%Get_ptr_r3d(ref=iqg)
      field => this%states(this%tracer_idx)%p
      IF (use_updated_moisture) THEN
        qv  => this%new_states(this%tracer_idx)%p%get_ptr_r3d(ref=field%ref_idx(1)) ! water vapor
        qc  => this%new_states(this%tracer_idx)%p%get_ptr_r3d(ref=field%ref_idx(2)) ! cloud water
        qi  => this%new_states(this%tracer_idx)%p%get_ptr_r3d(ref=field%ref_idx(3)) ! cloud ice
      ELSE
        qv  => this%states(this%tracer_idx)%p%get_ptr_r3d(ref=field%ref_idx(1)) ! water vapor
        qc  => this%states(this%tracer_idx)%p%get_ptr_r3d(ref=field%ref_idx(2)) ! cloud water
        qi  => this%states(this%tracer_idx)%p%get_ptr_r3d(ref=field%ref_idx(3)) ! cloud ice
      END IF

      CALL compute_temperature_from_internal_energy( &
        & this%domain, &
        & geo_height,  &
        & qr, qs, qg,  &
        & energy,      &
        & qv, qc, qi,  &
        & temperature  &
        & )
      END SELECT

  END SUBROUTINE energy_to_temp
  !
  !============================================================================
  !
  SUBROUTINE compute_flux_x(this, shflx, ufts, ufvs, flux_x)

    CLASS(t_vdf_atmo), INTENT(in) :: this
    REAL(wp), INTENT(in) :: &
      & shflx(:,:), & !< sensible heat flux
      & ufts(:,:),  & !< energy flux at surface from thermal exchange
      & ufvs(:,:)     !< energy flux at surface from vapor exchange
    REAL(wp), INTENT(out) :: flux_x(:,:)

    INTEGER :: jb, jc
    INTEGER,  POINTER :: energy_type
    REAL(wp), POINTER :: cpd, cvd

    CHARACTER(len=*), PARAMETER :: routine = modname//':energy_flux_to_flux_x'

    energy_type      => this%config%energy_type%Get_ptr_i0d()

    ASSOCIATE(domain => this%domain)

!$OMP PARALLEL
    CALL init(flux_x, lacc=.TRUE.)
!$OMP END PARALLEL

    SELECT CASE(energy_type)
    CASE (1)
      cpd => this%config%cpd%Get_ptr_r0d()
      cvd => this%config%cvd%Get_ptr_r0d()

!$OMP PARALLEL DO PRIVATE(jb, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c, domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          flux_x(jc,jb) = shflx(jc,jb) * cpd / cvd
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO

    CASE (2)

!$OMP PARALLEL DO PRIVATE(jb, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c, domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          flux_x(jc,jb) = ufts(jc,jb) + ufvs(jc,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO

    END SELECT

    END ASSOCIATE

    !$ACC WAIT(1)

  END SUBROUTINE compute_flux_x
  !
  !============================================================================
  !
  ! Subroutines for diagnostics
  !
  !============================================================================
  SUBROUTINE compute_static_energy(domain, spec_heat, temperature, geo_height, static_energy)

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), INTENT(in) :: spec_heat
    REAL(wp), DIMENSION(:,:,:), INTENT(in) :: &
      & temperature, geo_height
    REAL(wp), DIMENSION(:,:,:), INTENT(out) :: &
      & static_energy

    INTEGER :: jb, jk, jc

!$OMP PARALLEL
    CALL init(static_energy, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          static_energy(jc,jk,jb) = spec_heat * temperature(jc,jk,jb) + grav * geo_height(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_static_energy
  !
  !============================================================================
  !
  SUBROUTINE compute_temp_from_static_energy(domain, spec_heat, static_energy, geo_height, temperature)

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), INTENT(in) :: spec_heat
    REAL(wp), DIMENSION(:,:,:), INTENT(in) :: &
      & static_energy, geo_height
    REAL(wp), DIMENSION(:,:,:), INTENT(out) :: &
      & temperature

    INTEGER :: jb, jk, jc

!$OMP PARALLEL
    CALL init(temperature, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          temperature(jc,jk,jb) = (static_energy(jc,jk,jb) - grav * geo_height(jc,jk,jb)) / spec_heat
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_temp_from_static_energy
  !
  !============================================================================
  !
  SUBROUTINE compute_geopotential_height_above_ground(domain, zf, zh, ghf)

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), DIMENSION(:,:,:), INTENT(in) :: &
      & zf, zh
    REAL(wp), DIMENSION(:,:,:), INTENT(out) :: &
      & ghf

    INTEGER :: jb, jk, jc, nlevp1

    nlevp1 = domain%nlev + 1

!$OMP PARALLEL
    CALL init(ghf, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          ghf(jc,jk,jb) = zf(jc,jk,jb) - zh(jc,nlevp1,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_geopotential_height_above_ground
  !============================================================================
  !
  ! Compute mass specific internal energy + geopotential from temperature
  ! and moisture state
  !
  SUBROUTINE compute_internal_energy( &
    & domain,      &
    & geo_height,  &
    & qr, qs, qg,  &
    & temperature, &
    & qv, qc, qi,  &
    & energy       &
    & )

    USE mo_aes_thermo, ONLY: internal_energy

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), DIMENSION(:,:,:), INTENT(in) :: &
      & geo_height,  &
      & qr,          &
      & qs,          &
      & qg,          &
      & temperature, &
      & qv,          &
      & qc,          &
      & qi
    REAL(wp), DIMENSION(:,:,:), INTENT(out) :: &
      & energy

    INTEGER :: jb, jk, jc
    REAL(wp) :: q_liquid, q_solid

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_internal_energy'

!$OMP PARALLEL
    CALL init(energy, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, q_liquid, q_solid) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(q_liquid, q_solid)
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          q_liquid = qc(jc,jk,jb) + qr(jc,jk,jb)
          q_solid  = qi(jc,jk,jb) + qs(jc,jk,jb) + qg(jc,jk,jb)
          energy(jc,jk,jb) = &
            & internal_energy(           &
            &   temperature(jc,jk,jb), & ! temperature
            &   qv         (jc,jk,jb), & ! qv
            &   q_liquid,              & ! liquid
            &   q_solid,               & ! solid
            &   1._wp,                 & ! density
            &   1._wp                  & ! delta z
            &) + grav * geo_height(jc,jk,jb) * cvd/cpd
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_internal_energy
  !============================================================================
  !
  ! Compute temperature from mass specific internal energy + geopotential
  ! and moisture state
  !
  SUBROUTINE compute_temperature_from_internal_energy( &
    & domain,      &
    & geo_height,  &
    & qr, qs, qg,  &
    & energy,      &
    & qv, qc, qi,  &
    & temperature  &
    & )

    USE mo_aes_thermo, ONLY: T_from_internal_energy

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), DIMENSION(:,:,:), INTENT(in) :: &
      & geo_height, &
      & qr,         &
      & qs,         &
      & qg,         &
      & energy,     &
      & qv,         &
      & qc,         &
      & qi
    REAL(wp), DIMENSION(:,:,:), INTENT(out) :: &
      & temperature

    INTEGER :: jb, jk, jc
    REAL(wp) :: q_liquid, q_solid, u

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_temperature_from_internal_energy'

!$OMP PARALLEL
    CALL init(temperature, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, q_liquid, q_solid, u) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(q_liquid, q_solid, u)
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          q_liquid = qc(jc,jk,jb) + qr(jc,jk,jb)
          q_solid  = qi(jc,jk,jb) + qs(jc,jk,jb) + qg(jc,jk,jb)
          u        = energy(jc,jk,jb) - grav * geo_height(jc,jk,jb) * cvd/cpd
          temperature(jc,jk,jb) = &
            & T_from_internal_energy(  &
            &   u,                     & ! internal energy
            &   qv         (jc,jk,jb), & ! qv
            &   q_liquid,              & ! liquid
            &   q_solid,               & ! solid
            &   1._wp,                 & ! density
            &   1._wp                  & ! delta z
            &)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_temperature_from_internal_energy
  !
  !============================================================================
  !
  SUBROUTINE compute_internal_energy_vi( &
    & domain,      &
    & rho,         &
    & dz,          &
    & qr, qs, qg,  &
    & temperature, &
    & qv, qc, qi,  &
    & uvi          &
    & )

    USE mo_aes_thermo, ONLY: internal_energy

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), DIMENSION(:,:,:), INTENT(in) :: &
      & rho, &
      & dz, &
      & qr, &
      & qs, &
      & qg, &
      & temperature, &
      & qv, &
      & qc, &
      & qi
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & uvi

    INTEGER :: jb, jk, jc
    REAL(wp) :: q_liquid, q_solid

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_internal_energy_vi'

!$OMP PARALLEL
    CALL init(uvi, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, q_liquid, q_solid) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(q_liquid, q_solid)
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          q_liquid = qc(jc,jk,jb) + qr(jc,jk,jb)
          q_solid  = qi(jc,jk,jb) + qs(jc,jk,jb) + qg(jc,jk,jb)
          uvi(jc,jb) = uvi(jc,jb) + &
            & internal_energy(           &
            &   temperature(jc,jk,jb), & ! temperature
            &   qv         (jc,jk,jb), & ! qv
            &   q_liquid,              & ! liquid
            &   q_solid,               & ! solid
            &   rho(jc,jk,jb),         & ! density
            &   dz (jc,jk,jb)          & ! delta z
            &)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_internal_energy_vi
  !
  !============================================================================
  !
  SUBROUTINE  get_virtual_potential_temperature(patch, ptvm1, papm1, theta_v, &
                                                rl_start, rl_end)

    REAL(wp), INTENT(in), POINTER      :: ptvm1(:,:,:), papm1(:,:,:)
    TYPE(t_patch), INTENT(in), POINTER :: patch
    REAL(wp), INTENT(in), POINTER      :: theta_v(:,:,:)
    INTEGER,  INTENT(in)               :: rl_start, rl_end

    INTEGER  :: jb, jk, jc, nlev
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx

    nlev = SIZE(theta_v,2)
    i_startblk = patch%cells%start_block(rl_start)
    i_endblk   = patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          theta_v(jc,jk,jb) = ptvm1(jc,jk,jb)*(p0ref/papm1(jc,jk,jb))**rd_o_cpd
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE get_virtual_potential_temperature
  !============================================================================
  !
  ! This subroutine computes normal velocity component at the edges, based on
  ! zonal and meridional wind components at cell centers?
  !
  SUBROUTINE compute_normal_velocity_edge(                  &
    pum1, pvm1, patch, p_int, rl_start, rl_end, vn )

    REAL(wp), INTENT(in), POINTER :: pum1(:,:,:), pvm1(:,:,:)
    TYPE(t_patch), INTENT(in), POINTER :: patch
    TYPE(t_int_state), TARGET, INTENT(IN)  ::  p_int

    INTEGER, INTENT(in) :: rl_start,rl_end
    REAL(wp), INTENT(in), POINTER :: vn(:,:,:)

    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: jb, je, jk, jbn, jcn, nlev
    REAL(wp) :: zvn1, zvn2

    nlev = SIZE(vn,2)

    i_startblk = patch%edges%start_block(rl_start)
    i_endblk   = patch%edges%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, je, i_startidx, i_endidx, jcn, jbn, zvn1, zvn2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk ! domain%i_startblk_e, domain%i_endblk_e
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2) PRIVATE(jcn, jbn, zvn1, zvn2)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx ! domain%i_startidx_e(jb), domain%i_endidx_e(jb)
          jcn  = patch%edges%cell_idx(je,jb,1)
          jbn  = patch%edges%cell_blk(je,jb,1)

          zvn1 =   pum1(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,1)%v1 &
            &    + pvm1(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,1)%v2
          !
          jcn  =   patch%edges%cell_idx(je,jb,2)
          jbn  =   patch%edges%cell_blk(je,jb,2)
          zvn2 =   pum1(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,2)%v1 &
            &    + pvm1(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,2)%v2
          !
          vn(je,jk,jb) = p_int%c_lin_e(je,1,jb)*zvn1 + p_int%c_lin_e(je,2,jb)*zvn2
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_normal_velocity_edge
  !
  !============================================================================
  !
  SUBROUTINE interpolate_normal_velocity_edge_interface(                  &
    vn, patch, p_nh_metrics, rl_start, rl_end, vn_ie )

    REAL(wp), INTENT(in), POINTER :: vn(:,:,:)
    TYPE(t_patch), INTENT(in), POINTER :: patch
    TYPE(t_nh_metrics),INTENT(in) :: p_nh_metrics

    INTEGER, INTENT(in) :: rl_start,rl_end
    REAL(wp), INTENT(in), POINTER :: vn_ie(:,:,:)

    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: jb, je, jk, nlev, nlevp1

    nlev = SIZE(vn,2)
    nlevp1 = nlev + 1

    i_startblk = patch%edges%start_block(rl_start)
    i_endblk   = patch%edges%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 2, nlev
        DO je = i_startidx, i_endidx
          vn_ie(je,jk,jb) = p_nh_metrics%wgtfac_e(je,jk,jb) * vn(je,jk,jb) +            &
                            ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * vn(je,jk-1,jb)
        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO je = i_startidx, i_endidx
        vn_ie(je,1,jb)      =   p_nh_metrics%wgtfacq1_e(je,1,jb) * vn(je,1,jb) &
                              + p_nh_metrics%wgtfacq1_e(je,2,jb) * vn(je,2,jb) &
                              + p_nh_metrics%wgtfacq1_e(je,3,jb) * vn(je,3,jb)

        vn_ie(je,nlevp1,jb) =   p_nh_metrics%wgtfacq_e(je,1,jb) * vn(je,nlev,jb)   &
                              + p_nh_metrics%wgtfacq_e(je,2,jb) * vn(je,nlev-1,jb) &
                              + p_nh_metrics%wgtfacq_e(je,3,jb) * vn(je,nlev-2,jb)
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE interpolate_normal_velocity_edge_interface
  !============================================================================
  !
  ! This subroutine computes the velocity gradient tensor at the cell edge.
  !
  ! Output:
  !   vel_grad_e[i,j]:   velocity gradient tensor, 3x3 tensor
  !                       -1. index: wind component at cell edge (1: normal, 2: tangential, 3: vertial)
  !                       -2. index: derivative direction (1: normal, 2: tangential, 3: vertial)
  !
  SUBROUTINE compute_velocity_gradient_tensor(                 &
    u_vert, v_vert, w_vert, pwp1, vn_ie, vt_ie, w_ie,          &
    patch, p_nh_metrics, rl_start, rl_end, vel_grad_e          &
    )

    REAL(wp), INTENT(in) :: u_vert(:,:,:), v_vert(:,:,:), w_vert(:,:,:), pwp1(:,:,:)
    REAL(wp), INTENT(in) :: vn_ie(:,:,:), vt_ie(:,:,:), w_ie(:,:,:)

    TYPE(t_patch), INTENT(in), POINTER :: patch
    TYPE(t_nh_metrics),INTENT(in) :: p_nh_metrics

    INTEGER, INTENT(in) :: rl_start,rl_end

    TYPE(t_vel_grad_tensor), INTENT(inout) :: vel_grad_e(3,3)

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, jb, jk, je, nlev

    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4
    REAL(wp) :: vt_vert1, vt_vert2, vt_vert3, vt_vert4
    REAL(wp) :: w_full_c1, w_full_c2, w_full_v1, w_full_v2

    nlev = SIZE(u_vert,2)

    i_startblk = patch%edges%start_block(rl_start)
    i_endblk   = patch%edges%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, je, i_startidx, i_endidx, vn_vert1, vn_vert2, vn_vert3, vn_vert4,  &
!$OMP                     vt_vert1, vt_vert2, vt_vert3, vt_vert4, w_full_c1, w_full_c2, w_full_v1, &
!$OMP                     w_full_v2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR TILE(32, 4) &
      !$ACC   PRIVATE(vn_vert1, vn_vert2, vn_vert3, vn_vert4, vt_vert1, vt_vert2, vt_vert3, vt_vert4) &
      !$ACC   PRIVATE(w_full_c1, w_full_c2, w_full_v1, w_full_v2)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          ! Get normal velocity component at vertices for target edge
          !TODO: can we modify the subroutine get_normal_velocity_vertex from tmx_numerics
          !      to use within OpenACC loop?
          vn_vert1 =    u_vert(patch%edges%vertex_idx(je,jb,1),jk,patch%edges%vertex_blk(je,jb,1)) &
                        * patch%edges%primal_normal_vert(je,jb,1)%v1                               &
                      + v_vert(patch%edges%vertex_idx(je,jb,1),jk,patch%edges%vertex_blk(je,jb,1)) &
                        * patch%edges%primal_normal_vert(je,jb,1)%v2

          vn_vert2 =    u_vert(patch%edges%vertex_idx(je,jb,2),jk,patch%edges%vertex_blk(je,jb,2)) &
                        * patch%edges%primal_normal_vert(je,jb,2)%v1                               &
                      + v_vert(patch%edges%vertex_idx(je,jb,2),jk,patch%edges%vertex_blk(je,jb,2)) &
                        * patch%edges%primal_normal_vert(je,jb,2)%v2

          vn_vert3 =    u_vert(patch%edges%vertex_idx(je,jb,3),jk,patch%edges%vertex_blk(je,jb,3)) &
                        * patch%edges%primal_normal_vert(je,jb,3)%v1                               &
                      + v_vert(patch%edges%vertex_idx(je,jb,3),jk,patch%edges%vertex_blk(je,jb,3)) &
                        * patch%edges%primal_normal_vert(je,jb,3)%v2

          vn_vert4 =    u_vert(patch%edges%vertex_idx(je,jb,4),jk,patch%edges%vertex_blk(je,jb,4)) &
                        * patch%edges%primal_normal_vert(je,jb,4)%v1                               &
                      + v_vert(patch%edges%vertex_idx(je,jb,4),jk,patch%edges%vertex_blk(je,jb,4)) &
                        * patch%edges%primal_normal_vert(je,jb,4)%v2


          ! Get tangential velocity component at vertices for target edge
          vt_vert1 =    u_vert(patch%edges%vertex_idx(je,jb,1),jk,patch%edges%vertex_blk(je,jb,1)) &
                        * patch%edges%dual_normal_vert(je,jb,1)%v1                                 &
                      + v_vert(patch%edges%vertex_idx(je,jb,1),jk,patch%edges%vertex_blk(je,jb,1)) &
                        * patch%edges%dual_normal_vert(je,jb,1)%v2

          vt_vert2 =    u_vert(patch%edges%vertex_idx(je,jb,2),jk,patch%edges%vertex_blk(je,jb,2)) &
                        * patch%edges%dual_normal_vert(je,jb,2)%v1                                 &
                      + v_vert(patch%edges%vertex_idx(je,jb,2),jk,patch%edges%vertex_blk(je,jb,2)) &
                        * patch%edges%dual_normal_vert(je,jb,2)%v2

          vt_vert3 =    u_vert(patch%edges%vertex_idx(je,jb,3),jk,patch%edges%vertex_blk(je,jb,3)) &
                        * patch%edges%dual_normal_vert(je,jb,3)%v1                                 &
                      + v_vert(patch%edges%vertex_idx(je,jb,3),jk,patch%edges%vertex_blk(je,jb,3)) &
                        * patch%edges%dual_normal_vert(je,jb,3)%v2

          vt_vert4 =    u_vert(patch%edges%vertex_idx(je,jb,4),jk,patch%edges%vertex_blk(je,jb,4)) &
                        * patch%edges%dual_normal_vert(je,jb,4)%v1                                 &
                      + v_vert(patch%edges%vertex_idx(je,jb,4),jk,patch%edges%vertex_blk(je,jb,4)) &
                        * patch%edges%dual_normal_vert(je,jb,4)%v2

          ! W at full levels
          w_full_c1  = 0.5_wp * (  pwp1(patch%edges%cell_idx(je,jb,1),jk,  patch%edges%cell_blk(je,jb,1))    &
                                 + pwp1(patch%edges%cell_idx(je,jb,1),jk+1,patch%edges%cell_blk(je,jb,1)) )

          w_full_c2  = 0.5_wp * (  pwp1(patch%edges%cell_idx(je,jb,2),jk,  patch%edges%cell_blk(je,jb,2))    &
                                 + pwp1(patch%edges%cell_idx(je,jb,2),jk+1,patch%edges%cell_blk(je,jb,2)) )

          ! W at full levels vertices from w at vertices at interface levels
          w_full_v1  = 0.5_wp * (  w_vert(patch%edges%vertex_idx(je,jb,1),jk,  patch%edges%vertex_blk(je,jb,1))   &
                                 + w_vert(patch%edges%vertex_idx(je,jb,1),jk+1,patch%edges%vertex_blk(je,jb,1)) )

          w_full_v2  = 0.5_wp * (  w_vert(patch%edges%vertex_idx(je,jb,2),jk,  patch%edges%vertex_blk(je,jb,2))   &
                                 + w_vert(patch%edges%vertex_idx(je,jb,2),jk+1,patch%edges%vertex_blk(je,jb,2)) )

          ! Compute velocity gradient tensor at edge of full levels
          ! e.g. vgrad_e(1,2) = du_1/dx_2
          ! with index notation for triangles:
          !     1: normal, 2: tangential, 3: vertical direction
          ! TODO: can we modify get_velocity_gradient_tensor_edge to use within OpenACC loop?
          vel_grad_e(1,1)%ptr(je,jk,jb) = ( vn_vert4 - vn_vert3 ) * patch%edges%inv_vert_vert_length(je,jb)

          vel_grad_e(1,2)%ptr(je,jk,jb) = ( vn_vert2 - vn_vert1 )                      &
                                          * patch%edges%tangent_orientation(je,jb)     &
                                          * patch%edges%inv_primal_edge_length(je,jb)

          vel_grad_e(1,3)%ptr(je,jk,jb) = ( vn_ie(je,jk,jb) - vn_ie(je,jk+1,jb) )      &
                                          * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)

          vel_grad_e(2,1)%ptr(je,jk,jb) = ( vt_vert4-vt_vert3 )                        &
                                          * patch%edges%inv_vert_vert_length(je,jb)

          vel_grad_e(2,2)%ptr(je,jk,jb) = ( vt_vert2-vt_vert1 )                        &
                                          * patch%edges%tangent_orientation(je,jb)     &
                                          * patch%edges%inv_primal_edge_length(je,jb)

          vel_grad_e(2,3)%ptr(je,jk,jb) = ( vt_ie(je,jk,jb) - vt_ie(je,jk+1,jb) )      &
                                          * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)

          vel_grad_e(3,1)%ptr(je,jk,jb) = ( w_full_c2 - w_full_c1 )                    &
                                          * patch%edges%inv_dual_edge_length(je,jb)

          vel_grad_e(3,2)%ptr(je,jk,jb) = ( w_full_v2 - w_full_v1 )                    &
                                          * patch%edges%tangent_orientation(je,jb)     &
                                          * patch%edges%inv_primal_edge_length(je,jb)

          vel_grad_e(3,3)%ptr(je,jk,jb) = ( w_ie(je,jk,jb) - w_ie(je,jk+1,jb) )        &
                                          * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)

        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_velocity_gradient_tensor
  !============================================================================
  !
  ! This subroutine computes the shear (rate of strane |S|) and trace of the
  ! strain-rate tensor (S_ij).
  !
  !   shear = 2 * |S|^2 =  D_ij * D_ij
  !
  !     with  S_ij = 0.5 ( du_i/dx_j + du_j/dx_i )
  !
  !
  ! Output:
  !   shear:          2 * |S|^2 (|S|=sqrt(2*S_ij*S_ij))
  !   div_of_stress:  trace(S_ij) = S_jj = 0.5 * D_jj = du_j/dx_j
  !
  ! shear_e = = 2 * |S|^2  (|S|: rate-of-strain)
  ! with D_ij = 2 * S_ij = du_i/dx_j + du_j/dx_i
  ! and  |S| = sqrt(2 * S_ij * S_ij)
  !
  SUBROUTINE compute_shear(                                   &
    vel_grad_e, patch, rl_start, rl_end,        &
    shear, div_of_stress                                      &
    )

    TYPE(t_vel_grad_tensor), INTENT(in) :: vel_grad_e(3,3)

    TYPE(t_patch), INTENT(in), POINTER :: patch
    INTEGER, INTENT(in) :: rl_start,rl_end

    REAL(wp), INTENT(inout) :: shear(:,:,:), div_of_stress(:,:,:)

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: jb, jk, je, nlev

    REAL(wp) :: D_12,D_13,D_23

    nlev = SIZE(div_of_stress,2)

    i_startblk = patch%edges%start_block(rl_start)
    i_endblk   = patch%edges%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, je, i_startidx, i_endidx, D_12, D_13, D_23) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR TILE(32, 4) &
      !$ACC   PRIVATE(D_12, D_13, D_23)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
         ! Compute local shear at edge:
         !   shear_e = 2 * |S|^2 (|S|=sqrt(2*S_ij*S_ij))
         ! Mechanical prod is half of this value divided by km

          ! Strain rates at edge center
          D_12       = vel_grad_e(1,2)%ptr(je,jk,jb) + vel_grad_e(2,1)%ptr(je,jk,jb)
          D_13       = vel_grad_e(1,3)%ptr(je,jk,jb) + vel_grad_e(3,1)%ptr(je,jk,jb)
          D_23       = vel_grad_e(2,3)%ptr(je,jk,jb) + vel_grad_e(3,2)%ptr(je,jk,jb)

          ! Mechanical prod is half of this value divided by km
           shear(je,jk,jb) = 4._wp * (  vel_grad_e(1,1)%ptr(je,jk,jb)**2._wp  &
                                      + vel_grad_e(2,2)%ptr(je,jk,jb)**2._wp  &
                                      + vel_grad_e(3,3)%ptr(je,jk,jb)**2._wp) &
                            + 2._wp * ( D_12**2._wp + D_13**2._wp + D_23**2._wp)

          ! Trace of strain-rate tensor S_ij:
          !   trace(S_ij) = S_jj = 0.5 * D_jj = du_j/dx_j
          div_of_stress(je,jk,jb) =   vel_grad_e(1,1)%ptr(je,jk,jb) &
                                    + vel_grad_e(2,2)%ptr(je,jk,jb) &
                                    + vel_grad_e(3,3)%ptr(je,jk,jb)
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_shear
  !============================================================================
  !
  ! This function comuptes the local velocity gradient tensor at the edges of
  ! the full levels, based on central finite differences.
  !
  ! Normal gradient of the normal & tangential velocity components are computed
  ! between vertices 3 and 4 (vertex numbering see figure 1 in Zaengl et al.
  ! 2015, Q. J. R. Meteorol. Soc.).
  ! Normal gradients of vertical velocity component is computed between
  ! adjacent cell centers.
  ! Tangential gradients are computed between vertices 1 and 2, and vertical
  ! gradients between edges at adjacent half levels.
  !
  ! Output:
  !   vgrad_e:  velocity gradient tensor (1. index: velocity component;
  !                                       2. index: derivative direction)
  !
  FUNCTION get_velocity_gradient_tensor_edge(               &
    patch, p_nh_metrics, je, jb, jk,                        &
    vn_vert1,vn_vert2,vn_vert3,vn_vert4,                    &
    vt_vert1,vt_vert2,vt_vert3,vt_vert4,                    &
    w_full_c1,w_full_c2,w_full_v1,w_full_v2,                &
    vn_ie, vt_ie, w_ie                                      &
    ) RESULT(vgrad_e)

    TYPE(t_patch), INTENT(in), POINTER :: patch
    TYPE(t_nh_metrics),INTENT(in) :: p_nh_metrics

    INTEGER, INTENT(in) :: je, jb, jk

    REAL(wp), INTENT(in) :: vn_vert1,vn_vert2,vn_vert3,vn_vert4
    REAL(wp), INTENT(in) :: vt_vert1,vt_vert2,vt_vert3,vt_vert4
    REAL(wp), INTENT(in) :: w_full_c1,w_full_c2,w_full_v1,w_full_v2

    REAL(wp), INTENT(in), POINTER :: vn_ie(:,:,:), vt_ie(:,:,:), w_ie(:,:,:)

    REAL(wp) :: vgrad_e(3,3)

    vgrad_e(1,1) = ( vn_vert4 - vn_vert3 )                            &
                      * patch%edges%inv_vert_vert_length(je,jb)
    vgrad_e(1,2) = ( vn_vert2 - vn_vert1 )                            &
                      * patch%edges%tangent_orientation(je,jb)        &
                      * patch%edges%inv_primal_edge_length(je,jb)
    vgrad_e(1,3) =  ( vn_ie(je,jk,jb) - vn_ie(je,jk+1,jb) )           &
                      * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)

    vgrad_e(2,1) = ( vt_vert4-vt_vert3 )                              &
                      * patch%edges%inv_vert_vert_length(je,jb)
    vgrad_e(2,2) = ( vt_vert2-vt_vert1 )                              &
                      * patch%edges%tangent_orientation(je,jb)        &
                      * patch%edges%inv_primal_edge_length(je,jb)
    vgrad_e(2,3) = ( vt_ie(je,jk,jb) - vt_ie(je,jk+1,jb) )            &
                      * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)

    vgrad_e(3,1) = ( w_full_c2 - w_full_c1 )                          &
                      * patch%edges%inv_dual_edge_length(je,jb)
    vgrad_e(3,2) = ( w_full_v2 - w_full_v1 )                          &
                      * patch%edges%tangent_orientation(je,jb)        &
                      * patch%edges%inv_primal_edge_length(je,jb)
    vgrad_e(3,3) = ( w_ie(je,jk,jb) - w_ie(je,jk+1,jb) )              &
                      * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)

  END FUNCTION get_velocity_gradient_tensor_edge
  !============================================================================
  !
  ! Computes horizontal divergence at cell mass obtained by Gauss theorem and
  ! fluxes defined at edges. (see equation (19) in Zaengl et al.
  ! 2015, Q. J. R. Meteorol. Soc.).
  !
  SUBROUTINE get_horizontal_divergence_strain_rate_cell(                &
    flux_e,patch,ptr_int,rl_start,rl_end,hdiv_c)

    REAL(wp), INTENT(in), POINTER :: flux_e(:,:,:)
    TYPE(t_patch), INTENT(in), POINTER :: patch
    TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

    INTEGER, INTENT(in) :: rl_start,rl_end
    REAL(wp), INTENT(in), POINTER :: hdiv_c(:,:,:)

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: jb, jk, jc, nlev

    nlev = SIZE(hdiv_c,2)

    i_startblk = patch%cells%start_block(rl_start)
    i_endblk   = patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          hdiv_c(jc,jk,jb) =                                                                       &
              (   flux_e(patch%cells%edge_idx(jc,jb,1),jk,patch%cells%edge_blk(jc,jb,1))   &
                  * ptr_int%e_bln_c_s(jc,1,jb)                                                      &
                + flux_e(patch%cells%edge_idx(jc,jb,2),jk,patch%cells%edge_blk(jc,jb,2))   &
                  * ptr_int%e_bln_c_s(jc,2,jb)                                                      &
                + flux_e(patch%cells%edge_idx(jc,jb,3),jk,patch%cells%edge_blk(jc,jb,3))   &
                  * ptr_int%e_bln_c_s(jc,3,jb)                                                      &
              )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE get_horizontal_divergence_strain_rate_cell
  !============================================================================
  !
  ! Interpolates the rate of strain from edges at full levels to cell center at
  ! half levels.
  !
  SUBROUTINE interpolate_rate_of_strain_full2half_edge2cell(                &
    shear_e,                                                                &
    patch,                                                                  &
    p_nh_metrics,                                                           &
    ptr_int,                                                                &
    rl_start,rl_end,                                                        &
    shear_ic)

    REAL(wp), INTENT(in), POINTER :: shear_e(:,:,:)
    TYPE(t_patch), INTENT(in), POINTER :: patch
    TYPE(t_nh_metrics),INTENT(in) :: p_nh_metrics
    TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

    INTEGER, INTENT(in) :: rl_start,rl_end
    REAL(wp), INTENT(in), POINTER :: shear_ic(:,:,:)

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, jb, jk, jc, nlev

    nlev = SIZE(shear_e,2)

    i_startblk = patch%cells%start_block(rl_start)
    i_endblk   = patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
          shear_ic(jc,jk,jb) =                                                                  &
                  p_nh_metrics%wgtfac_c(jc,jk,jb)                                               &
                  * (   shear_e(patch%cells%edge_idx(jc,jb,1),jk,patch%cells%edge_blk(jc,jb,1))   &
                        * ptr_int%e_bln_c_s(jc,1,jb)                                            &
                      + shear_e(patch%cells%edge_idx(jc,jb,2),jk,patch%cells%edge_blk(jc,jb,2))   &
                        * ptr_int%e_bln_c_s(jc,2,jb)                                            &
                      + shear_e(patch%cells%edge_idx(jc,jb,3),jk,patch%cells%edge_blk(jc,jb,3))   &
                        * ptr_int%e_bln_c_s(jc,3,jb)                                            &
                    )                                                                           &
                + ( 1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb) )                                   &
                  * (   shear_e(patch%cells%edge_idx(jc,jb,1),jk-1,patch%cells%edge_blk(jc,jb,1)) &
                        * ptr_int%e_bln_c_s(jc,1,jb)                                            &
                      + shear_e(patch%cells%edge_idx(jc,jb,2),jk-1,patch%cells%edge_blk(jc,jb,2)) &
                        * ptr_int%e_bln_c_s(jc,2,jb)                                            &
                      + shear_e(patch%cells%edge_idx(jc,jb,3),jk-1,patch%cells%edge_blk(jc,jb,3)) &
                        * ptr_int%e_bln_c_s(jc,3,jb)                                            &
                    )
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO
  END SUBROUTINE interpolate_rate_of_strain_full2half_edge2cell
  !============================================================================
  !
  ! Interpolates the eddy viscosity from cell center at half levels to cell
  ! center of full levels.
  ! The minimum return value of km_c is km_min.
  !
  SUBROUTINE interpolate_eddy_viscosity2cell(                               &
    km_ic,                                                                  &
    km_min,                                                                 &
    patch,                                                                  &
    rl_start,rl_end,                                                        &
    km_c )

    REAL(wp), INTENT(in), POINTER :: km_ic(:,:,:)
    REAL(wp), INTENT(in) :: km_min
    TYPE(t_patch), INTENT(in), POINTER :: patch

    INTEGER, INTENT(in) :: rl_start,rl_end
    REAL(wp), INTENT(in), POINTER :: km_c(:,:,:)

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, jb, jk, jc, nlev

    nlev = SIZE(km_c,2)

    i_startblk = patch%cells%start_block(rl_start)
    i_endblk   = patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          km_c(jc,jk,jb) = MAX( km_min,                                   &
                                ( km_ic(jc,jk,jb) + km_ic(jc,jk+1,jb) ) * &
                                0.5_wp )
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE interpolate_eddy_viscosity2cell
  !============================================================================
  !
  ! Interpolates the eddy viscosity from cell center at half levels to vertices
  ! at half levels.
  ! The minimum return value of km_c is km_min.
  !
  SUBROUTINE interpolate_eddy_viscosity2half_vertex(                        &
    km_ic,                                                                  &
    km_min,                                                                 &
    patch,                                                                  &
    ptr_int,                                                                &
    km_iv)

    REAL(wp), INTENT(in) :: km_ic(:,:,:)
    REAL(wp), INTENT(in) :: km_min
    TYPE(t_patch), INTENT(in), POINTER :: patch
    TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

    REAL(wp), INTENT(inout) :: km_iv(:,:,:)

    INTEGER :: jb, jk, jc, nblks, nlev, nproma

    nblks = SIZE(km_iv,3)
    nproma = SIZE(km_iv, 1)
    nlev = SIZE(km_iv, 2)

    CALL cells2verts_scalar(km_ic, patch, ptr_int%cells_aw_verts, km_iv, &
                            lacc=.TRUE., opt_rlstart=5, opt_rlend=min_rlvert_int-1,   &
                            opt_acc_async=.TRUE.)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
    DO jb = 1, nblks
      DO jk = 1, nlev
        DO jc = 1, nproma
          km_iv(jc,jk,jb) = MAX( km_min,  km_iv(jc,jk,jb) )
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE interpolate_eddy_viscosity2half_vertex
  !============================================================================
  !
  ! Interpolates the eddy viscosity from cell center at half levels to edges
  ! at half levels.
  ! The minimum return value of km_c is km_min.
  !
  SUBROUTINE interpolate_eddy_viscosity2half_edge(                        &
    km_ic,                                                                &
    km_min,                                                               &
    patch,                                                                &
    ptr_int,                                                              &
    km_ie)

    REAL(wp), INTENT(in) :: km_ic(:,:,:)
    REAL(wp), INTENT(in) :: km_min
    TYPE(t_patch), INTENT(in), POINTER :: patch
    TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

    REAL(wp), INTENT(inout) :: km_ie(:,:,:)

    INTEGER :: jb, jk, jc, nblks, nlev, nproma

    nblks = SIZE(km_ie,3)
    nproma = SIZE(km_ie, 1)
    nlev = SIZE(km_ie, 2)

    CALL cells2edges_scalar(km_ic, patch, ptr_int%c_lin_e, km_ie, lacc=.TRUE.,      &
                            opt_rlstart=grf_bdywidth_e, opt_rlend=min_rledge_int-1  )

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
    DO jb = 1, nblks
      DO jk = 1, nlev
        DO jc = 1, nproma
          km_ie(jc,jk,jb) = MAX( km_min,  km_ie(jc,jk,jb) )
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE interpolate_eddy_viscosity2half_edge
  !============================================================================
  !
  ! working precision version of prepare_diffusion_matrix.
  ! The coefficients of the system of equations for the
  ! implicit calculation of the tendencies are set up here.
  ! The coefficients for explicit calculations differ from
  ! those for implicit calculations only by one term. This
  ! terms is depending on the time increment and is omitted
  ! at this point. This allows to use the subroutine for
  ! calculating both, the explicit and the implicit
  ! coefficients. In case of implicit treatment it is added
  ! to the coefficient b later (see module mo_tmx_numerics;
  ! subroutine diffuse_vertical_implicit).
  !
  SUBROUTINE prepare_diffusion_matrix_wp( &
    & ics, ice,              & ! in
    & minlvl, maxlvl,        & ! in
    & lhalflvl,              & ! in
    & inv_mair,              & ! in
    & inv_dz,                & ! in
    & zk,                    & ! in
    & zprefac,               & ! in
    & a, b, c                & ! out
    & )

    ! Logical variable that takes into account whether the
    ! calculation is done on half or full levels
    LOGICAL, INTENT(in) :: lhalflvl

    ! Iteration boundaries for blocks, cells, and level
    INTEGER, INTENT(in) :: ics, ice, minlvl, maxlvl

    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & inv_dz       ! inverse distance between cell centers/interfaces [1 / m]

    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & inv_mair,  & ! inverse moist air mass [m2 / kg]
      & zk           ! turbulent diffusion coefficient multiplied by density [kg / (m * s)]

    ! factor containing the turbulent diffusion coefficient
    REAL(wp), OPTIONAL, INTENT(in) ::  zprefac

    ! Set up the system of equations of shape
    ! a*x_(k-1) + b*x_(k) + c*x_(k+1) = rhs,
    ! where x is the variable at time step t+1.
    REAL(wp), INTENT(out), DIMENSION(:,:) :: a, b, c

    ! Iterators for blocks, cells, and levels
    ! The correction factors lvlcorr_a and lvlcorr_c
    ! are used to address the difference between computations
    ! on half levels and on full levels.
    INTEGER  :: jc, jk, lvlcorr_a, lvlcorr_c, jk_corr_a, jk_corr_c

    ! Multiplier requiered for some coefficients
    REAL(wp) :: zmulti

    CHARACTER(len=*), PARAMETER :: routine = modname//':prepare_diffusion_matrix_wp'

    ! For half levels the coefficient "a" is calclated using
    ! infomation on the upper half level, i.e. jk-1, and the
    ! coefficient "c" using information on the current level jk.
    ! For full levels this is shifted, so that the coefficient
    ! "a" is calculated using information on the current level jk
    ! and coefficient "c" using information on the lower full level,
    ! i.e. a cell with index jk+1.
    IF(lhalflvl) THEN
      lvlcorr_a = -1
      lvlcorr_c =  0
    ELSE
      lvlcorr_a =  0
      lvlcorr_c =  1
    ENDIF

    IF(PRESENT(zprefac)) THEN
      zmulti = zprefac
    ELSE
      zmulti= 1._wp
    END IF

    ! Set up the tri-diagonal matrix
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR PRIVATE(jk_corr_a, jk_corr_c) COLLAPSE(2) ASYNC(1)
    DO jk=minlvl+1,maxlvl-1
      DO jc = ics, ice
        jk_corr_a = jk + lvlcorr_a
        jk_corr_c = jk + lvlcorr_c
        a(jc,jk) = - zmulti * zk(jc,jk_corr_a) * inv_dz(jc,jk_corr_a) * inv_mair(jc,jk)
        c(jc,jk) = - zmulti * zk(jc,jk_corr_c) * inv_dz(jc,jk_corr_c) * inv_mair(jc,jk)
        b(jc,jk) = - a(jc,jk) - c(jc,jk)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    jk_corr_a = minlvl + lvlcorr_a
    jk_corr_c = minlvl + lvlcorr_c
    ! Set up the upper boundary condition
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jc = ics, ice
      a(jc,minlvl) = 0._wp
      c(jc,minlvl) = - zmulti * zk(jc,jk_corr_c) * inv_dz(jc,jk_corr_c) * inv_mair(jc,minlvl)
      b(jc,minlvl) = - c(jc,minlvl)
    END DO
    !$ACC END PARALLEL LOOP

    jk_corr_a = maxlvl + lvlcorr_a
    jk_corr_c = maxlvl + lvlcorr_c
    ! Set up the lower boundary condition
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jc = ics, ice
      a(jc,maxlvl) = - zmulti * zk(jc,jk_corr_a) * inv_dz(jc,jk_corr_a) * inv_mair(jc,maxlvl)
      c(jc,maxlvl) = 0._wp
      b(jc,maxlvl) = - a(jc,maxlvl)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

  END SUBROUTINE prepare_diffusion_matrix_wp

#ifdef __MIXED_PRECISION
  !============================================================================
  !
  ! mixed precision version of prepare_diffusion_matrix.
  ! The coefficients of the system of equations for the
  ! implicit calculation of the tendencies are set up here.
  ! The coefficients for explicit calculations differ from
  ! those for implicit calculations only by one term. This
  ! terms is depending on the time increment and is omitted
  ! at this point. This allows to use the subroutine for
  ! calculating both, the explicit and the implicit
  ! coefficients. In case of implicit treatment it is added
  ! to the coefficient b later (see module mo_tmx_numerics;
  ! subroutine diffuse_vertical_implicit).
  !
  SUBROUTINE prepare_diffusion_matrix_mp( &
    & ics, ice,              & ! in
    & minlvl, maxlvl,        & ! in
    & lhalflvl,              & ! in
    & inv_mair,              & ! in
    & inv_dz,                & ! in
    & zk,                    & ! in
    & zprefac,               & ! in
    & a, b, c                & ! out
    & )

    ! Logical variable that takes into account whether the
    ! calculation is done on half or full levels
    LOGICAL, INTENT(in) :: lhalflvl

    ! Iteration boundaries for blocks, cells, and level
    INTEGER, INTENT(in) :: ics, ice, minlvl, maxlvl

    REAL(vp), INTENT(in), DIMENSION(:,:) :: &
      & inv_dz       ! inverse distance between cell centers/interfaces [1 / m]

    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & inv_mair,  & ! inverse moist air mass [m2 / kg]
      & zk           ! turbulent diffusion coefficient multiplied by density [kg / (m * s)]

    ! factor containing the turbulent diffusion coefficient
    REAL(wp), OPTIONAL, INTENT(in) ::  zprefac

    ! Set up the system of equations of shape
    ! a*x_(k-1) + b*x_(k) + c*x_(k+1) = rhs,
    ! where x is the variable at time step t+1.
    REAL(wp), INTENT(out), DIMENSION(:,:) :: a, b, c

    ! Iterators for blocks, cells, and levels
    ! The correction factors lvlcorr_a and lvlcorr_c
    ! are used to address the difference between computations
    ! on half levels and on full levels.
    INTEGER  :: jc, jk, lvlcorr_a, lvlcorr_c, jk_corr_a, jk_corr_c

    ! Multiplier requiered for some coefficients
    REAL(wp) :: zmulti

    CHARACTER(len=*), PARAMETER :: routine = modname//':prepare_diffusion_matrix'

    ! For half levels the coefficient "a" is calclated using
    ! infomation on the upper half level, i.e. jk-1, and the
    ! coefficient "c" using information on the current level jk.
    ! For full levels this is shifted, so that the coefficient
    ! "a" is calculated using information on the current level jk
    ! and coefficient "c" using information on the lower full level,
    ! i.e. a cell with index jk+1.
    IF(lhalflvl) THEN
      lvlcorr_a = -1
      lvlcorr_c =  0
    ELSE
      lvlcorr_a =  0
      lvlcorr_c =  1
    ENDIF

    IF(PRESENT(zprefac)) THEN
      zmulti = zprefac
    ELSE
      zmulti= 1._wp
    END IF

    ! Set up the tri-diagonal matrix
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR PRIVATE(jk_corr_a, jk_corr_c) COLLAPSE(2) ASYNC(1)
    DO jk=minlvl+1,maxlvl-1
      DO jc = ics, ice
        jk_corr_a = jk + lvlcorr_a
        jk_corr_c = jk + lvlcorr_c
        a(jc,jk) = - zmulti * zk(jc,jk_corr_a) * inv_dz(jc,jk_corr_a) * inv_mair(jc,jk)
        c(jc,jk) = - zmulti * zk(jc,jk_corr_c) * inv_dz(jc,jk_corr_c) * inv_mair(jc,jk)
        b(jc,jk) = - a(jc,jk) - c(jc,jk)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    jk_corr_a = minlvl + lvlcorr_a
    jk_corr_c = minlvl + lvlcorr_c
    ! Set up the upper boundary condition
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jc = ics, ice
      a(jc,minlvl) = 0._wp
      c(jc,minlvl) = - zmulti * zk(jc,jk_corr_c) * inv_dz(jc,jk_corr_c) * inv_mair(jc,minlvl)
      b(jc,minlvl) = - c(jc,minlvl)
    END DO
    !$ACC END PARALLEL LOOP

    jk_corr_a = maxlvl + lvlcorr_a
    jk_corr_c = maxlvl + lvlcorr_c
    ! Set up the lower boundary condition
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jc = ics, ice
      a(jc,maxlvl) = - zmulti * zk(jc,jk_corr_a) * inv_dz(jc,jk_corr_a) * inv_mair(jc,maxlvl)
      c(jc,maxlvl) = 0._wp
      b(jc,maxlvl) = - a(jc,maxlvl)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

  END SUBROUTINE prepare_diffusion_matrix_mp
#endif
  !
  ! This subroutine assigns a constant eddy viscosity and diffusivity (Km=const=Kh).
  ! For validation purposes of turbunce model.
  !
  SUBROUTINE Assign_constant_eddy_viscosity( &
    domain,                     &
    rho_ic,                     &
    km_const,                   &
    rturb_prandtl,              &
    km_ic,                      &
    kh_ic                       &
    )

    TYPE(t_domain), INTENT(in)    :: domain
    REAL(wp), INTENT(in), DIMENSION(:,:,:)  :: rho_ic
    REAL(wp), INTENT(in) :: km_const, rturb_prandtl

    REAL(wp), INTENT(inout), DIMENSION(:,:,:) :: km_ic, kh_ic

    INTEGER :: jb,jc,jk,nlev,nlevp1
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end

    nlev = domain%nlev
    nlevp1 = nlev+1

    rl_start   = 3
    rl_end     = min_rlcell_int
    i_startblk = domain%patch%cells%start_block(rl_start)
    i_endblk   = domain%patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(domain%patch, jb, i_startblk, i_endblk, &
                              i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2 , nlev
#else
      DO jk = 2 , nlev
        DO jc = i_startidx, i_endidx
#endif
          km_ic(jc,jk,jb) = rho_ic(jc,jk,jb) * km_const

          kh_ic(jc,jk,jb) = km_ic(jc,jk,jb) * rturb_prandtl

        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jc = i_startidx, i_endidx
        kh_ic(jc,1,jb)      = kh_ic(jc,2,jb)
        kh_ic(jc,nlevp1,jb) = kh_ic(jc,nlev,jb)
        km_ic(jc,1,jb)      = km_ic(jc,2,jb)
        km_ic(jc,nlevp1,jb) = km_ic(jc,nlev,jb)
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE Assign_constant_eddy_viscosity
  !
  !=================================================================
  !
END MODULE mo_vdf_atmo
