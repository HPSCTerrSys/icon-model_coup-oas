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
!! This module is an additional interface between nwp_gw_interface to the
!! gravity wave drag parametrization developed at the Goethe Uni:
!! Lagrangian WKB raytracer in phase space.
!!
!! Revision History
!! Initial version by Gergely Boeloeni, Goethe Uni Frankfurt (2016-03-02)
!!
MODULE mo_nwp_msgwam_interface
#ifdef __MSGWAM

  USE mo_kind,                  ONLY: wp, vp, sp
  USE mo_mpi,                   ONLY: my_process_is_stdio, p_wait, work_mpi_barrier
  USE mo_exception,             ONLY: message, finish, message_text, &
                                      debug_on, debug_off
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_diag, t_nh_prog, t_nh_metrics
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_tend
  USE mo_impl_constants,        ONLY: min_rlcell_int, min_rledge_int, success
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_physical_constants,    ONLY: grav, rd, cpd, cvd, rd_o_cpd, earth_angular_velocity
  USE mo_grid_config,           ONLY: grid_sphere_radius, is_plane_torus
  USE mo_sync,                  ONLY: sync_patch_array, sync_patch_array_mult, SYNC_E, &
                                      SYNC_C, SYNC_C1
  USE mo_parallel_config,       ONLY: nproma
  USE mo_run_config,            ONLY: msg_level, ldynamics
  USE mo_dynamics_config,       ONLY: lcoriolis
  USE mo_math_constants,        ONLY: pi, pi2, pi_2, rad2deg, deg2rad
  USE mo_math_gradients,        ONLY: grad_green_gauss_cell
  USE mo_util_vgrid_types,      ONLY: vgrid_buffer
  USE mtime,                    ONLY: datetime, timeDelta, newTimedelta, &
                                      deallocateTimedelta, getTimedeltaFromDatetime, &
                                      getTotalMillisecondsTimedelta
  USE mo_vertical_grid,         ONLY: nrdmax
  USE mo_intp,                  ONLY: cell_avg
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_intp_lonlat_baryctr,   ONLY: inside_triangle
  USE mo_delaunay_types,        ONLY: t_point
  USE mo_math_divrot,           ONLY: div, div_avg
  USE mo_nh_testcases_nml,      ONLY: nh_test_name
  USE mo_timer,                 ONLY: timers_level, timer_start, timer_stop, &
                                      timer_msgwam_fieldsgrads, timer_msgwam_remove_rays, &
                                      timer_msgwam_init_gw_orretal, timer_msgwam_saturation, &
                                      timer_msgwam_wave2grid, timer_msgwam_diagprof, &
                                      timer_msgwam_propagate_wave, timer_msgwam_smooth_hori, &
                                      timer_msgwam_tendency, timer_msgwam_sync_wave, &
                                      timer_msgwam_split_merge
  USE mo_msgwam_config,         ONLY: alpha_sat, dt_add, imethod_merge, imethod_split, &
                                      jklaunch_bg, jray_offset_bg, jray_offset_cv, l1ray, &
                                      lcalc_flux_4dir_bg, ldiagprof, lhsmooth, llimittend, &
                                      lmsgwam_noforce, lmsgwam_offline, lmsgwam_pmomflux, &
                                      lsaturation, lsteady, ltest_hprop, ltest_restart, &
                                      nhsmooth, nrays, nrays_add_bg, &
                                      nrays_add_cv, nrays_bg, nrays_cv
  USE mo_setup_msgwam_interface,ONLY: iout_msgwam, msgwam_write_restartfiles, nlaunch_max_bg, &
                                      nlaunch_max_cv, p_rwork, t_msgwam, p_ray, p_gridinfo4ray, &
                                      ndiag_msgwam, jb_diag, p_msgwam, p_mgmgrid, t_gridinfo4ray, &
                                      jc_diag, p_spl, p_lfluxbg
  USE mo_msgwam_util,           ONLY: datout, idx_rayedge, smooth_hori, test_blowup_uvt
  USE mo_gw_source_config,      ONLY: gws_conv_config
  USE mo_gw_source_bg,          ONLY: init_gw_orretal
  USE mo_gw_source_conv,        ONLY: init_gw_conv
  USE mo_msgwam_diagnostics,    ONLY: project_action, project_diagnostics
  USE mo_msgwam_splitmerge,     ONLY: split_rays, remove_rays
  USE mo_msgwam,                ONLY: saturation, wave2grid, propagate_wave
  USE mo_msgwam_stst,           ONLY: gwdrag_msgwam_stst
  USE mo_nwp_msgwam_utils,      ONLY: fieldsgrads, regrid_wave, sync_wave, tendency


  ! Modules for computing the Orr et al. 2010 launch spectrum
  USE data_gwd,    ONLY : nslope, gfluxlaun, &
    &                     ggaussa, ggaussb, ngauss, gcoeff, lozpr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nwp_msgwam_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
  SUBROUTINE nwp_msgwam_interface( dt_call,                   & !>input
        &                          mtime_datetime,            & !>input
        &                          p_sim_time,                & !>input
        &                          p_patch, p_metrics,        & !>input
        &                          p_int_state,               & !>input
        &                          p_diag,p_prog,             & !>input
        &                          prm_nwp_tend               ) !>inout
    TYPE(datetime),      POINTER,INTENT(IN)   :: mtime_datetime       !< date/time information
    TYPE(t_patch),        TARGET,INTENT(INOUT):: p_patch              ! grid/patch info.
    TYPE(t_nh_metrics)          ,INTENT(IN)   :: p_metrics
    TYPE(t_int_state),    TARGET,INTENT(IN)   :: p_int_state
    TYPE(t_nh_prog),      TARGET,INTENT(IN)   :: p_prog               ! the dyn prog vars
    TYPE(t_nh_diag),      TARGET,INTENT(IN)   :: p_diag               ! the dyn diag vars
    TYPE(t_nwp_phy_tend), TARGET,INTENT(INOUT):: prm_nwp_tend         ! atm tend vars
    REAL(wp),                    INTENT(IN)   :: dt_call              ! time step
    REAL(wp),                    INTENT(IN)   :: p_sim_time           ! elapsed simulation time

    INTEGER  :: jg
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk ! blocks
    INTEGER  :: i_startidx, i_endidx ! slices
    INTEGER  :: jb                   ! block indeces


  ! Later use this ifdef
  !#ifdef msgwam

    IF (nrays(p_patch%id) == 0)  RETURN

    ! Domain ID
    jg = p_patch%id

    ! Decide whether we run MS-GWaM or its steady state version
    IF (.NOT. lsteady) THEN

      IF ( ltest_restart(2) )  &   ! flag to test restart-file reading
        &  CALL msgwam_write_restartfiles( mtime_datetime, p_patch )

      CALL gwdrag_msgwam ( dt_call,                     & !>input
        &                  mtime_datetime,              & !>input
        &                  p_sim_time,                  & !>input
        &                  p_patch,                     & !>input
        &                  p_metrics,                   & !>input
        &                  p_int_state,                 & !>input
        &                  p_prog% rho,                 & !>input
        &                  p_diag% temp,                & !>input
        &                  p_prog% theta_v,             & !>input
        &                  p_diag% u,                   & !>input
        &                  p_diag% v,                   & !>input
        &                  p_msgwam(jg)                 ) !>inout

    ELSE

      CALL gwdrag_msgwam_stst ( mtime_datetime,&
        &                       p_sim_time ,                   & !>input
        &                       jg,                            & !>input
        &                       p_metrics,                     & !>input
        &                       p_prog% rho,                   & !>input
        &                       p_diag% u,                     & !>input
        &                       p_diag% v,                     & !>input
        &                       p_diag% temp,                  & !>input
        &                       p_msgwam(jg),                  & !>inout
        &                       p_patch%nlev,                  & !>input
        &                       p_patch%nlevp1,                & !>input
        &                       p_patch%cells%start_block,     & !>input
        &                       p_patch%cells%end_block,       & !>input
        &                       p_patch%nshift_total,          &
        &                       p_patch%cells%center(:,:)%lat, &
        &                       p_patch%cells%f_c(:,:)         )

    ENDIF

    ! Transfer the resultant GWD to the NWP physics tendency array prm_nwp_tend
    IF (lmsgwam_offline .OR. lmsgwam_noforce) THEN
      CALL message('nwp_msgwam_interface', 'MS-GWaM in offline / noforce mode!')
    ELSE
      ! Exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

      IF (lmsgwam_pmomflux) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx)
        DO jb = i_startblk, i_endblk
          CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
            &                 i_startidx, i_endidx, rl_start, rl_end )

          prm_nwp_tend%        ddt_u_gwd         (i_startidx:i_endidx,:,jb)  &
            &  = p_msgwam(jg)% ddt_u_gwd_pmom_mgm(i_startidx:i_endidx,:,jb)
          prm_nwp_tend%        ddt_v_gwd         (i_startidx:i_endidx,:,jb)  &
            &  = p_msgwam(jg)% ddt_v_gwd_pmom_mgm(i_startidx:i_endidx,:,jb)
          prm_nwp_tend%        ddt_temp_drag     (i_startidx:i_endidx,:,jb) = 0._wp
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ELSE
!$OMP PARALLEL
!$OMP DO PRIVATE (jb, i_startidx, i_endidx)
        DO jb = i_startblk, i_endblk
          CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
            &                 i_startidx, i_endidx, rl_start, rl_end )

          prm_nwp_tend%        ddt_u_gwd    (i_startidx:i_endidx,:,jb)  &
            &  = p_msgwam(jg)% ddt_u_gwd_mgm(i_startidx:i_endidx,:,jb)
          prm_nwp_tend%        ddt_v_gwd    (i_startidx:i_endidx,:,jb)  &
            &  = p_msgwam(jg)% ddt_v_gwd_mgm(i_startidx:i_endidx,:,jb)
          prm_nwp_tend%        ddt_temp_drag(i_startidx:i_endidx,:,jb)  &
            &  = p_msgwam(jg)% ddt_t_gwd_mgm(i_startidx:i_endidx,:,jb)
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ENDIF
    END IF

  !#endif

  END SUBROUTINE nwp_msgwam_interface
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE gwdrag_msgwam ( dt_call,                   & ! input
                        &   mtime_datetime,            & ! input
                        &   p_sim_time,                & ! input
                        &   p_patch,p_metrics,         & ! input
                        &   p_int_state,               & ! input
                        &   rho,                       & ! input
                        &   temp    ,                  & ! input
                        &   theta,                     & ! input
                        &   u       ,                  & ! input
                        &   v       ,                  & ! input
                        &   p_fld                      ) ! inout

    TYPE(datetime),      POINTER,INTENT(IN)    :: mtime_datetime       ! date/time information
    TYPE(t_patch),        TARGET,INTENT(INOUT) :: p_patch              ! grid/patch info.
    TYPE(t_nh_metrics)          ,INTENT(IN)    :: p_metrics
    TYPE(t_int_state),    TARGET,INTENT(IN)    :: p_int_state          ! interpolation state
    TYPE(t_msgwam),              INTENT(INOUT) :: p_fld                ! the atm phys vars
    REAL(wp),                    INTENT(IN)    :: dt_call              ! time step
    REAL(wp),                    INTENT(IN)    :: p_sim_time           ! elapsed simulation time on this grid level
    REAL(wp),            POINTER,INTENT(IN)    :: rho (:,:,:)
    REAL(wp),            POINTER,INTENT(IN)    :: temp(:,:,:)
    REAL(wp),            POINTER,INTENT(IN)    :: theta(:,:,:)
    REAL(wp),            POINTER,INTENT(IN)    :: u   (:,:,:)
    REAL(wp),            POINTER,INTENT(IN)    :: v   (:,:,:)
    ! Local array bounds:
    INTEGER                                    :: nlev, nlevp1         ! number of full and half levels
    INTEGER                                    :: rl_start, rl_end
    INTEGER                                    :: i_startblk, i_endblk ! blocks
    INTEGER                                    :: i_startidx, i_endidx ! slices
    INTEGER                                    :: jk,jc,jb             ! block indeces
    INTEGER                                    :: maxindexu(3), maxindexv(3)
    ! Background
    REAL(wp)                                   :: rho_half(nproma,p_patch%nlevp1, p_patch%nblks_c)      ! rho at half levels
    REAL(wp)                                   :: bvf2_full(nproma, p_patch%nlev, p_patch%nblks_c)      ! N**2 at full levels
    REAL(wp)                                   :: bvf2_half(nproma, p_patch%nlevp1, p_patch%nblks_c)    ! N**2 at half levels
    REAL(wp)                                   :: u_fld(nproma, p_patch%nlevp1, p_patch%nblks_c)        ! u at half levels
    REAL(wp)                                   :: v_fld(nproma, p_patch%nlevp1, p_patch%nblks_c)        ! v at half levels
    REAL(wp)                                   :: kvisc(nproma, p_patch%nlevp1, p_patch%nblks_c)        ! kinematic viscosity
    REAL(wp)                                   :: gammash2_full(nproma, p_patch%nlev, p_patch%nblks_c)  ! inverse pinc scale height
                                                                                                        ! pinc:
                                                                                                        ! pseudo-incompressible
    REAL(wp)                                   :: gammash2_half(nproma, p_patch%nlevp1, p_patch%nblks_c)! inverse pinc scale height
                                                                                                        ! pinc:
                                                                                                        ! pseudo-incompressible
    REAL(wp)                                   :: dn2dz(nproma, p_patch%nlevp1, p_patch%nblks_c)        ! vertical grad of bvf2 /2
    REAL(wp)                                   :: dg2dz(nproma, p_patch%nlevp1, p_patch%nblks_c)        ! vertical grad of gammash2 /2
    REAL(wp)                                   :: dudz (nproma, p_patch%nlevp1, p_patch%nblks_c)        ! vertical (radial) grad of u
    REAL(wp)                                   :: dvdz (nproma, p_patch%nlevp1, p_patch%nblks_c)        ! vertical (radial) grad of v
    REAL(wp)                                   :: n2_hgrad(2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! zonal/meridional grad of N**2
    REAL(wp)                                   :: g2_hgrad(2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! zonal/meridional grad of gammash2
    REAL(wp)                                   :: u_hgrad (2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! zonal/meridional grad of u
    REAL(wp)                                   :: v_hgrad (2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! zonal/meridional grad of v
    REAL(wp)                                   :: dn2dz_hgrad(2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! vert & zonal/meridional grad of N**2
    REAL(wp)                                   :: dg2dz_hgrad(2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! vert & zonal/meridional grad of gammash2
    REAL(wp)                                   :: dudz_hgrad (2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! vert & zonal/meridional grad of u
    REAL(wp)                                   :: dvdz_hgrad (2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! vert & zonal/meridional grad of v
    REAL(wp)                                   :: fc2(nproma, p_patch%nblks_c)                          ! f_c**2
    REAL(wp)                                   :: fdfdlat(nproma, p_patch%nblks_c)                      ! meridional gradient of f**2
    ! 1D Diagnostics
    REAL(wp)                                   :: kd(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: kd2(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: ld(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: ld2(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: md(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: md2(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: cgz_diag(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: B2(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: B2_save(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: A(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: A_save(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: mB2(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: mB2_save(nproma,p_patch%nlevp1)
    REAL(wp)                                   :: prec(nproma)                                  ! total precip
  ! REAL(wp),                      PARAMETER   :: lim_tend = 0.05_wp
    REAL(wp),                      PARAMETER   :: lim_tend_tmp = 1._wp
    !
    LOGICAL                                    :: ltimeadd_bg
    LOGICAL                                    :: ltimeadd_cv
    INTEGER                                    :: iter_sm
    INTEGER                                    :: idiag
    INTEGER                                    :: jg
    REAL(wp)                                   :: coef_sphere, coef_torus

#ifndef __msgwam1d
    CALL message('gwdrag_msgwam', 'Transient MS-GWaM called')
#else
    CALL message('gwdrag_msgwam', 'Transient MS-GWaM-1D called')
    ! off : lhsmooth, hard-coded tendency limiter
    ! need to turn on  lmsgwam_pmomflux
    !              off imethod_split/merge (as not implemented)
#endif

    IF (is_plane_torus) THEN
      coef_sphere = 0.     ;  coef_torus = 1._wp
    ELSE
      coef_sphere = 1._wp  ;  coef_torus = 0.
    END IF

    ! Number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Domain ID
    jg     = p_patch%id

    ! Exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c + 1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    ! Decide whether new waves are emmited (new ray volumes are added)
    ! ltimeadd_bg = .T. --> Desaubies type background source launched
    ! ltimeadd_cv = .T. --> convective source launched
    IF (nrays_add_bg(jg) /= 0) THEN
      IF (ltest_hprop .AND. l1ray) THEN ! horizontal propagation test with single ray
                                                          ! or simulation with idealized source
        ltimeadd_bg = (MOD(p_sim_time,dt_add)<dt_call*0.999_wp .AND. p_sim_time <= 0._wp)
      ELSEIF (ltest_hprop) THEN ! horizontal propagation test
        ltimeadd_bg = (MOD(p_sim_time,dt_add)<dt_call*0.999_wp .AND. p_sim_time <= 43200._wp) ! half day emission
      ELSE
        ltimeadd_bg = (MOD(p_sim_time,dt_add)<dt_call*0.999_wp)
      ENDIF
    ELSE
      ltimeadd_bg = .FALSE.
    ENDIF
    IF ( gws_conv_config%n_source(jg) > 0 ) THEN
      ltimeadd_cv = ( p_sim_time /= 0._wp )
    ELSE
      ltimeadd_cv = .FALSE.
    ENDIF

    ! Idealized gwp test case: no background source, no convective source
    IF (nh_test_name=='gwp') THEN
      ltimeadd_bg = .FALSE.
      ltimeadd_cv = .FALSE.
    ENDIF

  ! Print whether source is launched
    IF (msg_level >= 12) THEN
      WRITE(message_text,'(a,L8)') 'ltimeadd_bg:', ltimeadd_bg
      IF ( gws_conv_config%n_source(jg) > 0 )  &
        &  WRITE(message_text,'(a,L8)') 'ltimeadd_cv:', ltimeadd_cv
      CALL message('', TRIM(message_text))
    ENDIF

    !=========================== FIELDS/GRADS ===========================
    ! Calculation of all resolved variables and their vertical gradients.
    !====================================================================

    IF (timers_level > 4) CALL timer_start(timer_msgwam_fieldsgrads)

    ! Calculate derived fields from the resolved flow (e.g. Brunt-Vaeisaelae freq and density scale height)
    ! and their vertical gradients (needed for solving the ray equations)
    CALL fieldsgrads(p_patch       = p_patch,                         & ! grid information                    (in)
                    p_int_state   = p_int_state,                     & ! interpolation data of patch         (in)
                    rl_start      = rl_start,                        & !                                     (in)
                    rl_end        = rl_end,                          & !                                     (in)
                    nlev          = nlev,                            & ! no. of full levels                  (in)
                    dz            = p_metrics%ddqz_z_full(:, :, :),  & ! full layer thickness                (in)
                    wgtfac_c      = p_metrics%wgtfac_c(:, :, :),     & !                                     (in)
                    rho           = rho (:, :, :),                   & ! density at full levels              (in)
                    temp          = temp(:, :, :),                   & ! temperature at full levels          (in)
                    u             = u(:, :, :),                      & ! zonal wind (full levels)            (in)
                    v             = v(:, :, :),                      & ! meridional wind (full levels)       (in)
                    fc2           = fc2 (:, :),                      & ! f_c**2                              (out)
                    fdfdlat       = fdfdlat(:, :),                   & ! meridional gradient of f**2         (out)
                    rho_half      = rho_half(:, :, :),               & ! density at half levels              (out)
                    u_fld         = u_fld(:, :, :),                  & ! zonal wind                          (out)
                    v_fld         = v_fld(:, :, :),                  & ! meridional wind                     (out)
                    bvf2_full     = bvf2_full(:, :, :),              & ! Brunt-Vaeisaelae frequency**2          (out)
                    bvf2_half     = bvf2_half(:, :, :),              & ! Brunt-Vaeisaelae frequency**2          (out)
                    gammash2_full = gammash2_full(:, :, :),          & ! inverse pinc scaleheight**2         (out)
                    gammash2_half = gammash2_half(:, :, :),          & ! inverse pinc scaleheight**2         (out)
                    kvisc         = kvisc(:, :, :),                  & ! Kinematic viscosity at half levels  (out)
                    dn2dz         = dn2dz(:, :, :),                  & ! vertical grad of bvf2               (out)
                    dg2dz         = dg2dz(:, :, :),                  & ! vertical grad of gammash2           (out)
                    dudz          = dudz (:, :, :),                  & ! vertical (radial) grad of u         (out)
                    dvdz          = dvdz (:, :, :),                  & ! vertical (radial) grad of v         (out)
                    n2_hgrad      = n2_hgrad(:, :, :, :),            & ! zonal / meridional grad of BVF**2   (out)
                    g2_hgrad      = g2_hgrad(:, :, :, :),            & ! zonal / meridional grad of gamma**2 (out)
                    u_hgrad       = u_hgrad(:, :, :, :),             & ! zonal / meridional grad of u        (out)
                    v_hgrad       = v_hgrad(:, :, :, :),             & ! zonal / meridional grad of v        (out)
                    dn2dz_hgrad   = dn2dz_hgrad(:, :, :, :),         & ! vert & zonal / meridional grad of BVF**2   (out)
                    dg2dz_hgrad   = dg2dz_hgrad(:, :, :, :),         & ! vert & zonal / meridional grad of gamma**2 (out)
                    dudz_hgrad    = dudz_hgrad(:, :, :, :),          & ! vert & zonal / meridional grad of u        (out)
                    dvdz_hgrad    = dvdz_hgrad(:, :, :, :)           ) ! vert & zonal / meridional grad of v        (out)

    IF (timers_level > 4) CALL timer_stop(timer_msgwam_fieldsgrads)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, prec,                      &
!$OMP            kd, kd2, ld, ld2, md, md2, cgz_diag, A, A_save, B2, B2_save, mB2, mB2_save)


    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

      ! Considering N^2 changes between time steps, remove ray volumes where N^2 <= f^2.
      ! This will not happen as long as bvf2_min > f^2 everywhere.
  !   DO jc = i_startidx, i_endidx
  !     DO jk = 1, nlevp1
  !       IF (bvf2_half(jc,jk,jb) > fc2(jc))  CYCLE
  !       WHERE ( p_ray(jg)% iexist(jc,:,jb) == jk )
  !         p_ray(jg)% iexist(jc,:,jb) = 0
  !         p_ray(jg)% specid(jc,:,jb) = 0
  !       END WHERE
  !     ENDDO
  !   ENDDO

      !======================== MERGE/REMOVE RAYS ==========================
      ! Remove ray volumes if there are more than allowed (maxrays_bg(jg)
      ! or maxrays_cv(jg)) in order to keep computational costs under control.
      ! The removal is done separately for ray volumes coming from the
      ! background and the convective sources.
      ! TODO: the removal should be replaced by a merging procedure
      !=====================================================================

      IF (timers_level > 4) CALL timer_start(timer_msgwam_remove_rays)

      !======================== Diagnostics ==================================
      ! WA diagnostic
      !=======================================================================

      CALL idx_rayedge( nlev,i_startidx,i_endidx,1,nrays(jg),& ! (in)
                        p_metrics%z_mc     (:,:,jb),        & ! (in)
                        p_metrics%z_ifc    (:,:,jb),        & ! (in)
                        p_ray(jg)%iexist   (:,:,jb),        & ! (in)
                        p_ray(jg)%specid   (:,:,jb),        & ! (in)
                        p_ray(jg)%jk_active(:,:,jb),        & ! (in)
                        p_ray(jg)%z        (:,:,jb),        & ! (in)
                        p_ray(jg)%dz       (:,:,jb),        & ! (in)
                        p_ray(jg)%jk_full_rtop(:,:,jb),     & ! (inout)
                        p_ray(jg)%jk_full_rbot(:,:,jb),     & ! (inout)
                        p_ray(jg)%jk_half_rtop(:,:,jb),     & ! (inout)
                        p_ray(jg)%jk_half_rbot(:,:,jb)      ) ! (inout)

      CALL project_action(nlev        = nlev,                             & !
                          i_startidx  = i_startidx,                       & !
                          i_endidx    = i_endidx,                         & !
                          jray_start  = 1,                                & !
                          jray_end    = nrays(jg),                        & !
                          z           = p_metrics%z_mc(:,:,jb),           & !
                          zhalf       = p_metrics%z_ifc(:,:,jb),          & !
                          cellarea    = p_patch%cells%area(:,jb),         & !
                          zray        = p_ray(jg)%z(:,:,jb),              & !
                          dzray       = p_ray(jg)%dz(:,:,jb),             & !
                          dyray       = p_ray(jg)%dy(:,:,jb),             & !
                          dxray       = p_ray(jg)%dx(:,:,jb),             & !
                          dkray       = p_ray(jg)%dk(:,:,jb),             & !
                          dlray       = p_ray(jg)%dl(:,:,jb),             & !
                          dmray       = p_ray(jg)%dm(:,:,jb),             & !
                          dens        = p_ray(jg)%wadens(:,:,jb),         & !
                          iexist      = p_ray(jg)%iexist(:,:,jb),         & !
                          specid      = p_ray(jg)%specid(:,:,jb),         & !
                          jk_active   = p_ray(jg)%jk_active(:,:,jb),      & !
                          jkmin_half  = p_ray(jg)%jk_half_rtop(:,:,jb),   & !
                          jkmax_half  = p_ray(jg)%jk_half_rbot(:,:,jb),   & !
                          action      = p_fld%action_mgm_1(:,:,jb),       & !
                          diag        = 0,                                & !
                          jg          = jg                                ) !


      ! Remove ray volumes from convective sources when merge-scheme is disabled
      IF ( ltimeadd_cv .AND. ( imethod_merge <= 0 ) ) THEN
        CALL remove_rays(      nlev         = nlev,                            & ! no. of full levels       (in)
                              i_startidx   = i_startidx,                      & ! first index of the block (in)
                              i_endidx     = i_endidx,                        & ! last index of the block  (in)
                              jray_start   = jray_offset_cv(jg)+1,            &
                              jray_end     = jray_offset_cv(jg)+nrays_cv(jg), &
                              nrays_vacate = nrays_add_cv(jg)*nlaunch_max_cv, &
                              dxray        = p_ray(jg)%dx(:,:,jb),            & ! ray size in lon          (in)
                              dyray        = p_ray(jg)%dy(:,:,jb),            & ! ray size in lat          (in)
                              dzray        = p_ray(jg)%dz(:,:,jb),            & ! ray size in z            (in)
                              kray         = p_ray(jg)%k(:,:,jb),             & ! horiz (lon) wavenumber   (in)
                              dkray        = p_ray(jg)%dk(:,:,jb),            & ! horiz (lon) wavenumber   (in)
                              lray         = p_ray(jg)%l(:,:,jb),             & ! horiz (lat) wavenumber   (in)
                              dlray        = p_ray(jg)%dl(:,:,jb),            & ! horiz (lat) wavenumber   (in)
                              mray         = p_ray(jg)%m(:,:,jb),             & ! vertical wavenumber      (in)
                              dmray        = p_ray(jg)%dm(:,:,jb),            & ! vertical wavenumber      (in)
                              dens         = p_ray(jg)%wadens(:,:,jb),        & ! wave action density      (inout)
                              iexist       = p_ray(jg)%iexist(:,:,jb),        & ! existence of ray         (inout)
                              specid       = p_ray(jg)%specid(:,:,jb),        & ! spectral id of rays      (inout)
                              bvf2         = bvf2_half(:,:,jb),               & ! Brunt-Vaeisala freq**2    (in)
                              gammash2     = gammash2_half(:,:,jb),           & ! inverse pinc scale height(in)
                              fc2          = fc2(:,jb)                        ) ! Coriolis parameter**2    (in)
      END IF

      ! Remove ray volumes from the Desaubies type background source when merge-scheme is disabled
      IF ( ltimeadd_bg .AND. ( imethod_merge <= 0 ) ) THEN
        CALL remove_rays(      nlev         = nlev,                            & ! no. of full levels       (in)
                              i_startidx   = i_startidx,                      & ! first index of the block (in)
                              i_endidx     = i_endidx,                        & ! last index of the block  (in)
                              jray_start   = jray_offset_bg(jg)+1,            &
                              jray_end     = jray_offset_bg(jg)+nrays_bg(jg), &
                              nrays_vacate = nrays_add_bg(jg)*nlaunch_max_bg, &
                              dxray        = p_ray(jg)%dx(:,:,jb),            & ! ray size in lon          (in)
                              dyray        = p_ray(jg)%dy(:,:,jb),            & ! ray size in lat          (in)
                              dzray        = p_ray(jg)%dz(:,:,jb),            & ! ray size in z            (in)
                              kray         = p_ray(jg)%k(:,:,jb),             & ! horiz (lon) wavenumber   (in)
                              dkray        = p_ray(jg)%dk(:,:,jb),            & ! horiz (lon) wavenumber   (in)
                              lray         = p_ray(jg)%l(:,:,jb),             & ! horiz (lat) wavenumber   (in)
                              dlray        = p_ray(jg)%dl(:,:,jb),            & ! horiz (lat) wavenumber   (in)
                              mray         = p_ray(jg)%m(:,:,jb),             & ! vertical wavenumber      (in)
                              dmray        = p_ray(jg)%dm(:,:,jb),            & ! vertical wavenumber      (in)
                              dens         = p_ray(jg)%wadens(:,:,jb),        & ! wave action density      (inout)
                              iexist       = p_ray(jg)%iexist(:,:,jb),        & ! existence of ray         (inout)
                              specid       = p_ray(jg)%specid(:,:,jb),        & ! spectral id of rays      (inout)
                              bvf2         = bvf2_half(:,:,jb),               & ! Brunt-Vaeisala freq**2    (in)
                              gammash2     = gammash2_half(:,:,jb),           & ! inverse pinc scale height(in)
                              fc2          = fc2(:,jb)                        ) ! Coriolis parameter**2    (in)
      END IF

      IF (timers_level > 4) CALL timer_stop(timer_msgwam_remove_rays)


      !========================= INIT GW SOURCES ============================
      ! Launch GW sources (currently 3 kinds):
      ! 1) convective (init_gw_conv)
      ! 2) Desaubies type background (init_gw_orretal)
      ! 3) Idealized gwp (set up in src/testcases/mo_nh_isotherm_rest_atm_gwp.f90)
      ! TODO: mountain waves, jets/fronts
      !======================================================================

      ! Convective source
      IF ( ltimeadd_cv )                                                          &
        &  CALL init_gw_conv( nlev,jg,jb,i_startidx,i_endidx,jray_offset_cv(jg),  &
        &                     nlaunch_max_cv,                                     &
        &                     p_metrics%z_mc(:,:,jb),p_metrics%z_ifc(:,:,jb),     &
        &                     p_patch%cells%center(:,jb)%lon,                     &
        &                     p_patch%cells%center(:,jb)%lat,                     &
        &                     p_patch%cells%area(:,jb),                           &
        &                     p_fld% flag_cgw(:,jb) )

      ! Desaubies type background source
      IF (ltimeadd_bg) THEN

        IF (timers_level > 4) CALL timer_start(timer_msgwam_init_gw_orretal)

        ! Caculate total precipitation (input for the Orr et al., 2010 source)
  !     IF ( lozpr )  prec(:) =   prm_diag%rain_gsp_rate(:,jb) &  ! rain_gsp
  !       &                     + prm_diag%snow_gsp_rate(:,jb) &  ! snow_gsp
  !       &                     + prm_diag%rain_con_rate(:,jb) &  ! rain_con
  !       &                     + prm_diag%snow_con_rate(:,jb)    ! snow con
        IF ( lozpr )  prec(:) = 0._wp

        CALL init_gw_orretal(nlev      = nlev,                            & ! no. of full levels       (in)
                            i_startidx = i_startidx,                      & ! first index of the block (in)
                            i_endidx   = i_endidx,                        & ! last index of the block  (in)
                            jray_start = jray_offset_bg(jg)+1,            &
                            jray_end   = jray_offset_bg(jg)+nrays_bg(jg), &
                            mdatetime  = mtime_datetime,                  & ! last index of the block   (in)
                            zhalf      = p_metrics%z_ifc(:,:,jb),         & ! half level heights        (in)
                            z          = p_metrics%z_mc(:,:,jb),          & ! full level heights        (in)
                            specid     = p_ray(jg)%specid(:,:,jb),        & ! spectral id of rays       (inout)
                            lonray     = p_ray(jg)%lon(:,:,jb),           & ! ray position lon          (inout)
                            dxray      = p_ray(jg)%dx(:,:,jb),            & ! ray size in lon           (inout)
                            latray     = p_ray(jg)%lat(:,:,jb),           & ! ray position lat          (inout)
                            dyray      = p_ray(jg)%dy(:,:,jb),            & ! ray size in lat           (inout)
                            zray       = p_ray(jg)%z(:,:,jb),             & ! ray position z            (inout)
                            dzray      = p_ray(jg)%dz(:,:,jb),            & ! ray size in z             (inout)
                            coslatray  = p_ray(jg)%coslat(:,:,jb),        & ! cosine of ray lat         (inout)
                            kray       = p_ray(jg)%k(:,:,jb),             & ! ray horiz (lon) wavenumber(inout)
                            dkray      = p_ray(jg)%dk(:,:,jb),            & ! ray size in k             (inout)
                            lray       = p_ray(jg)%l(:,:,jb),             & ! ray horiz (lat) wavenumber(inout)
                            dlray      = p_ray(jg)%dl(:,:,jb),            & ! ray size in l             (inout)
                            mray       = p_ray(jg)%m(:,:,jb),             & ! ray vertical wavenumber   (inout)
                            dmray      = p_ray(jg)%dm(:,:,jb),            & ! ray size in m             (inout)
                            dens       = p_ray(jg)%wadens(:,:,jb),        & ! ray wave action density   (inout)
                            iexist     = p_ray(jg)%iexist(:,:,jb),        & ! existence of ray          (inout)
                            jk_active  = p_ray(jg)%jk_active(:,:,jb),     & ! jk index for launching    (inout)
                            jr_last    = p_ray(jg)%jr_last(:,:,jb),       & ! last-launched ray index   (inout)
                            jklaunch   = jklaunch_bg(jg),                 & ! launch level index        (in)
                            nlaunch_max= nlaunch_max_bg,                  & ! launch level index        (in)
                            cellarea   = p_patch%cells%area(:,jb),        & ! area of grid cel          (in)
                            lon        = p_patch%cells%center(:,jb)%lon,  & ! latitude of cell center   (in)
                            lat        = p_patch%cells%center(:,jb)%lat,  & ! latitude of cell center   (in)
                            coslat     = p_mgmgrid(jg)%coslat(:,jb),      & ! cosine of latitude        (in)
                            rhol       = rho_half(:,jklaunch_bg(jg),jb),  & ! rho at launch level       (in)
                            bvfl2      = bvf2_half(:,jklaunch_bg(jg),jb), & ! N**2 at launch level      (in)
                            lat_prof_bw= p_lfluxbg(jg)%lat_prof_bw(:,jb), & ! latitudinal factor        (in)
                            lat_prof_bs= p_lfluxbg(jg)%lat_prof_bs(:,jb), & ! latitudinal factor        (in)
                            lat_prof   = p_lfluxbg(jg)%lat_prof(:,jb),    & ! latitudinal factor        (inout)
                            prec       = prec(:),                         & ! total precipitation       (in)
                            jg         = jg                               )

        IF (timers_level > 4) CALL timer_stop(timer_msgwam_init_gw_orretal)

      ENDIF


      !======================== INDEX OF RAY EDGES ===========================
      ! Calculate the vertical level index (jk) closest to the bottom and top
      ! of the ray volume.
      ! TODO: something similar for horizontal edges if horizontal propagation
      !=======================================================================


      IF (timers_level > 4) CALL timer_start(timer_msgwam_saturation)    ! temporary

      CALL idx_rayedge( nlev,i_startidx,i_endidx,1,nrays(jg),& ! (in)
                        p_metrics%z_mc     (:,:,jb),        & ! (in)
                        p_metrics%z_ifc    (:,:,jb),        & ! (in)
                        p_ray(jg)%iexist   (:,:,jb),        & ! (in)
                        p_ray(jg)%specid   (:,:,jb),        & ! (in)
                        p_ray(jg)%jk_active(:,:,jb),        & ! (in)
                        p_ray(jg)%z        (:,:,jb),        & ! (in)
                        p_ray(jg)%dz       (:,:,jb),        & ! (in)
                        p_ray(jg)%jk_full_rtop(:,:,jb),     & ! (inout)
                        p_ray(jg)%jk_full_rbot(:,:,jb),     & ! (inout)
                        p_ray(jg)%jk_half_rtop(:,:,jb),     & ! (inout)
                        p_ray(jg)%jk_half_rbot(:,:,jb)      ) ! (inout)

      IF (timers_level > 4) CALL timer_stop(timer_msgwam_saturation)     !  temporary

      CALL project_action(nlev        = nlev,                             & !
                          i_startidx  = i_startidx,                       & !
                          i_endidx    = i_endidx,                         & !
                          jray_start  = 1,                                & !
                          jray_end    = nrays(jg),                        & !
                          z           = p_metrics%z_mc(:,:,jb),           & !
                          zhalf       = p_metrics%z_ifc(:,:,jb),          & !
                          cellarea    = p_patch%cells%area(:,jb),         & !
                          zray        = p_ray(jg)%z(:,:,jb),              & !
                          dzray       = p_ray(jg)%dz(:,:,jb),             & !
                          dyray       = p_ray(jg)%dy(:,:,jb),             & !
                          dxray       = p_ray(jg)%dx(:,:,jb),             & !
                          dkray       = p_ray(jg)%dk(:,:,jb),             & !
                          dlray       = p_ray(jg)%dl(:,:,jb),             & !
                          dmray       = p_ray(jg)%dm(:,:,jb),             & !
                          dens        = p_ray(jg)%wadens(:,:,jb),         & !
                          iexist      = p_ray(jg)%iexist(:,:,jb),         & !
                          specid      = p_ray(jg)%specid(:,:,jb),         & !
                          jk_active   = p_ray(jg)%jk_active(:,:,jb),      & !
                          jkmin_half  = p_ray(jg)%jk_half_rtop(:,:,jb),   & !
                          jkmax_half  = p_ray(jg)%jk_half_rbot(:,:,jb),   & !
                          action      = p_fld%action_mgm_2(:,:,jb),       & !
                          diag        = 0,                                & !
                          jg          = jg                                ) !


      !============================ SATURATION ==============================
      ! Wave breaking scheme based on static instability criterion: if the GW
      ! turns the potential temperture gradient to negative at a certain
      ! height, the wave action of ray volumes are reduced so that static
      ! stability sets in again. Saturation is (should be!) calculated for
      ! GWs from all sources (background, convective) together
      !======================================================================

      IF (lsaturation) THEN

        IF (timers_level > 4) CALL timer_start(timer_msgwam_saturation)

        CALL saturation(cellarea   = p_patch%cells%area(:,jb),    & !                          (in)
                        nlev       = nlev,                        & ! no. of full levels       (in)
                        i_startidx = i_startidx,                  & ! first index of the block (in)
                        i_endidx   = i_endidx,                    & ! last index of the block  (in)
                        jray_start = 1,                           & !                          (in)
                        jray_end   = nrays(jg),                   & !                          (in)
                        z          = p_metrics%z_mc(:,:,jb),      & ! full level heights       (in)
                        jk_active  = p_ray(jg)%jk_active(:,:,jb), & ! launch level index       (in)
                        zray       = p_ray(jg)%z(:,:,jb),         & ! position                 (in)
                        dzray      = p_ray(jg)%dz(:,:,jb),        & ! size in z-dir            (in)
                        dyray      = p_ray(jg)%dy(:,:,jb),        & ! meridional extent        (in)
                        dxray      = p_ray(jg)%dx(:,:,jb),        & ! zonal extent             (in)
                        kray       = p_ray(jg)%k(:,:,jb),         & ! horiz (lon) wavenumber   (in)
                        dkray      = p_ray(jg)%dk(:,:,jb),        & ! size in k-dir            (in)
                        lray       = p_ray(jg)%l(:,:,jb),         & ! horiz (lat) wavenumber   (in)
                        dlray      = p_ray(jg)%dl(:,:,jb),        & ! size in l-dir            (in)
                        mray       = p_ray(jg)%m(:,:,jb),         & ! vertical wavenumber      (in)
                        dmray      = p_ray(jg)%dm(:,:,jb),        & ! size in m-dir            (in)
                        dens       = p_ray(jg)%wadens(:,:,jb),    & ! wave action density      (inout)
                        iexist     = p_ray(jg)%iexist(:,:,jb),    & ! existence of ray         (inout)
                        specid     = p_ray(jg)%specid(:,:,jb),    & ! spectral id of rays      (inout)
                        jkmin      = p_ray(jg)%jk_full_rtop(:,:,jb), & ! jk closest to ray-v top  (in)
                        jkmax      = p_ray(jg)%jk_full_rbot(:,:,jb), & ! jk closest to ray-v bot  (in)
                        bvf2       = bvf2_half(:,:,jb),           & ! Brunt-Vaeisala freq**2    (in)
                        gammash2   = gammash2_half(:,:,jb),       & ! inverse pinc scale height(in)
                        fc2        = fc2(:,jb),                   & ! Coriolis parameter**2    (in)
                        rho        = rho_half(:,:,jb),            & ! rho at half levels       (in)
                        kd         = kd(:,:),                     & ! k diagnostic             (out)
                        kd2        = kd2(:,:),                    & ! k^2 diagnostic           (out)
                        ld         = ld(:,:),                     & ! l diagnostic             (out)
                        ld2        = ld2(:,:),                    & ! l^2 diagnostic           (out)
                        md         = md(:,:),                     & ! m^2 diagnostic           (out)
                        md2        = md2(:,:),                    & ! m^2 diagnostic           (out)
                        cgz_diag   = cgz_diag(:,:),               & ! cgz diagnostic           (out)
                        A          = A(:,:),                      & ! B^2 diagnostic           (out)
                        A_save     = A_save(:,:),                 & ! B^2 diagnostic           (out)
                        B2         = B2(:,:),                     & ! B^2 diagnostic           (out)
                        B2_save    = B2_save(:,:),                & ! B^2 diagnostic           (out)
                        mB2_new    = mB2(:,:),                    & ! m^2B^2 diagnostic        (out)
                        mB2        = mB2_save(:,:),               & ! m^2B^2 diagnostic        (out)
                        jg=jg)

        IF (timers_level > 4) CALL timer_stop(timer_msgwam_saturation)

      ENDIF


      !============================ WAVE2GRID ===============================
      ! Project quantities of the Lagrangian ray volumes to the Eulerian grid
      ! and calculate momentum, pseudo-momentum, and temperature fluxes.
      ! The projection is done separately for GWs from the background and
      ! from the convective sources
      !======================================================================

      IF (timers_level > 4) CALL timer_start(timer_msgwam_wave2grid)

      ! wave2grid
      IF (nrays_add_bg(jg) /= 0) THEN
        CALL wave2grid(nlev            = nlev,                            & ! no. of full levels       (in)
                      i_startidx      = i_startidx,                      & ! first index of the block (in)
                      i_endidx        = i_endidx,                        & ! last index of the block  (in)
                      jray_start      = 1,                               & !                          (in)
                      jray_end        = nrays(jg),                       & !                          (in)
                      z               = p_metrics%z_mc(:,:,jb),          & ! full level heights       (in)
                      zhalf           = p_metrics%z_ifc(:,:,jb),         & ! half level heights       (in)
                      fc              = p_patch%cells%f_c(:,jb),         & ! Coriolis parameter       (in)
                      cellarea        = p_patch%cells%area(:,jb),        & ! area of grid cell        (in)
                      zray            = p_ray(jg)%z(:,:,jb),             & ! position                 (in)
                      dzray           = p_ray(jg)%dz(:,:,jb),            & ! size in z-dir            (in)
                      dyray           = p_ray(jg)%dy(:,:,jb),            & ! meridional extent        (in)
                      dxray           = p_ray(jg)%dx(:,:,jb),             & ! zonal extent             (in)
                      kray            = p_ray(jg)%k(:,:,jb),             & ! horiz (lon) wavenumber   (in)
                      dkray           = p_ray(jg)%dk(:,:,jb),            & ! size in k-dir            (in)
                      lray            = p_ray(jg)%l(:,:,jb),             & ! horiz (lat) wavenumber   (in)
                      dlray           = p_ray(jg)%dl(:,:,jb),            & ! size in l-dir            (in)
                      mray            = p_ray(jg)%m(:,:,jb),             & ! vertical wavenumber      (in)
                      dmray           = p_ray(jg)%dm(:,:,jb),            & ! size in m-dir            (in)
                      dens            = p_ray(jg)%wadens(:,:,jb),        & ! wave action density      (in)
                      iexist          = p_ray(jg)%iexist(:,:,jb),        & ! existence of ray         (in)
                      specid          = p_ray(jg)%specid(:,:,jb),        & ! spectral id of rays      (in)
                      jk_active       = p_ray(jg)%jk_active(:,:,jb),     & ! launch level index       (in)
                      jkmin_full      = p_ray(jg)%jk_full_rtop(:,:,jb),  & ! jk closest to ray-v top  (in)
                      jkmax_full      = p_ray(jg)%jk_full_rbot(:,:,jb),  & ! jk closest to ray-v bot  (in)
                      jkmin_half      = p_ray(jg)%jk_half_rtop(:,:,jb),  & ! jk closest to ray-v top  (in)
                      jkmax_half      = p_ray(jg)%jk_half_rbot(:,:,jb),  & ! jk closest to ray-v bot  (in)
                      theta           = theta(:,:,jb),                   & ! potential temperature    (in)
                      bvf2_full       = bvf2_full(:,:,jb),               & ! Brunt-Vaeisala freq**2    (in)
                      bvf2_half       = bvf2_half(:,:,jb),               & ! Brunt-Vaeisala freq**2    (in)
                      gammash2_full   = gammash2_full(:,:,jb),           & ! inverse pinc scale height(in)
                      gammash2_half   = gammash2_half(:,:,jb),           & ! inverse pinc scale height(in)
                      fc2             = fc2(:,jb),                       & ! Coriolis parameter**2    (in)
                      uuflux          = p_fld% uufl_mgm    (:,:,jb),     & ! uu mometum flux          (out)
                      uvflux          = p_fld% uvfl_mgm    (:,:,jb),     & ! uv mometum flux          (out)
                      uwflux          = p_fld% uwfl_mgm    (:,:,jb),     & ! uw mometum flux          (out)
                      vvflux          = p_fld% vvfl_mgm    (:,:,jb),     & ! vv mometum flux          (out)
                      vwflux          = p_fld% vwfl_mgm    (:,:,jb),     & ! vw mometum flux          (out)
                      uupflux         = p_fld% uupfl_mgm   (:,:,jb),     & ! uu pseudo mometum flux   (out)
                      uvpflux         = p_fld% uvpfl_mgm   (:,:,jb),     & ! uv pseudo mometum flux   (out)
                      uwpflux         = p_fld% uwpfl_mgm   (:,:,jb),     & ! uw pseudo mometum flux   (out)
                      vvpflux         = p_fld% vvpfl_mgm   (:,:,jb),     & ! vv pseudo mometum flux   (out)
                      vwpflux         = p_fld% vwpfl_mgm   (:,:,jb),     & ! vw pseudo mometum flux   (out)
                      utflux          = p_fld% utfl_mgm    (:,:,jb),     & ! u theta flux             (out)
                      vtflux          = p_fld% vtfl_mgm    (:,:,jb)      ) ! v theta flux             (out)
      ELSE
        p_fld% uufl_mgm(:,:,jb) = 0._wp   ; p_fld% uvfl_mgm(:,:,jb) = 0._wp  ; p_fld% uwfl_mgm(:,:,jb) = 0._wp
        p_fld% vvfl_mgm(:,:,jb) = 0._wp   ; p_fld% vwfl_mgm(:,:,jb) = 0._wp
        p_fld% uupfl_mgm(:,:,jb) = 0._wp  ; p_fld% uvpfl_mgm(:,:,jb) = 0._wp ; p_fld% uwpfl_mgm(:,:,jb) = 0._wp
        p_fld% vvpfl_mgm(:,:,jb) = 0._wp  ; p_fld% vwpfl_mgm(:,:,jb) = 0._wp
        p_fld% utfl_mgm(:,:,jb) = 0._wp   ; p_fld% vtfl_mgm(:,:,jb) = 0._wp
      END IF

      IF (timers_level > 4) CALL timer_stop(timer_msgwam_wave2grid)

      !============================ DIAGNOSTICS ===============================
      ! Project quantities of the Lagrangian ray volumes to the Eulerian grid
      ! and calculate everything that the model does not necessarily require.
      ! The projection is done separately for GWs from the background and
      ! from the convective sources
      !======================================================================

      IF (nrays_add_bg(jg) /= 0 .OR. gws_conv_config%n_source(jg) > 0) THEN
        call project_diagnostics(nlev           = nlev,                             & ! no. of full levels       (in)
                                i_startidx     = i_startidx,                       & ! first index of the block (in)
                                i_endidx       = i_endidx,                         & ! last index of the block  (in)
                                jray_start     = 1,                                & ! first index for ray      (in)
                                jray_end       = nrays(jg),                        & ! last index for ray       (in)
                                z              = p_metrics%z_mc(:,:,jb),           & ! full level heights       (in)
                                zhalf          = p_metrics%z_ifc(:,:,jb),          & ! half level heights       (in)
                                cellarea       = p_patch%cells%area(:,jb),         & ! area of grid cell        (in)
                                jkmin_full     = p_ray(jg)%jk_full_rtop(:,:,jb),   & ! jk closest to ray-v top  (in)
                                jkmax_full     = p_ray(jg)%jk_full_rbot(:,:,jb),   & ! jk closest to ray-v bot  (in)
                                jkmin_half     = p_ray(jg)%jk_half_rtop(:,:,jb),   & ! jk closest to ray-v top  (in)
                                jkmax_half     = p_ray(jg)%jk_half_rbot(:,:,jb),   & ! jk closest to ray-v bot  (in)
                                iexist         = p_ray(jg)%iexist(:,:,jb),         & ! existence of ray         (in)
                                jk_active      = p_ray(jg)%jk_active(:,:,jb),      & ! launch level index       (in)
                                specid         = p_ray(jg)%specid(:,:,jb),         & ! spectral id of rays      (in)
                                zray           = p_ray(jg)%z(:,:,jb),              & ! position                 (in)
                                dzray          = p_ray(jg)%dz(:,:,jb),             & ! size in z-dir            (in)
                                dyray          = p_ray(jg)%dy(:,:,jb),             & ! meridional extent        (in)
                                dxray          = p_ray(jg)%dx(:,:,jb),             & ! zonal extent             (in)
                                kray           = p_ray(jg)%k(:,:,jb),              & ! horiz (lon) wavenumber   (in)
                                dkray          = p_ray(jg)%dk(:,:,jb),             & ! size in k-dir            (in)
                                lray           = p_ray(jg)%l(:,:,jb),              & ! horiz (lat) wavenumber   (in)
                                dlray          = p_ray(jg)%dl(:,:,jb),             & ! size in l-dir            (in)
                                mray           = p_ray(jg)%m(:,:,jb),              & ! vertical wavenumber      (in)
                                dmray          = p_ray(jg)%dm(:,:,jb),             & ! size in m-dir            (in)
                                dens           = p_ray(jg)%wadens(:,:,jb),         & ! wave action density      (in)
                                bvf2_full      = bvf2_full(:,:,jb),                & ! Brunt-Vaeisala freq**2    (in)
                                bvf2_half      = bvf2_half(:,:,jb),                & ! Brunt-Vaeisala freq**2    (in)
                                gammash2_full  = gammash2_full(:,:,jb),            & ! inverse pinc scale height(in)
                                gammash2_half  = gammash2_half(:,:,jb),            & ! inverse pinc scale height(in)
                                fc             = p_patch%cells%f_c(:,jb),          & ! Coriolis parameter       (in)
                                fc2            = fc2(:,jb),                        & ! Coriolis parameter**2    (in)
                                u              = u(:, :, jb),                      & ! zonal wind (full levels)      (in)
                                v              = v(:, :, jb),                      & ! meridional wind (full levels) (in)
                                theta          = theta(:,:,jb),                    & ! potential temperature    (in)
                                pmflux_e       = p_fld% pmfl_mgm_e   (:,:,jb),      & ! pseudo mometum flux E-ward      (out)
                                pmflux_w       = p_fld% pmfl_mgm_w   (:,:,jb),      & ! pseudo mometum flux W-ward      (out)
                                pmflux_s       = p_fld% pmfl_mgm_s   (:,:,jb),      & ! pseudo mometum flux S-ward      (out)
                                pmflux_n       = p_fld% pmfl_mgm_n   (:,:,jb),      & ! pseudo mometum flux N-ward      (out)
                                apmflux        = p_fld% apmfl_mgm   (:,:,jb),      & ! absolute p-momentum flux (out)
                                mflux_e        = p_fld% mfl_mgm_e   (:,:,jb),      & ! mometum flux E-ward      (out)
                                mflux_w        = p_fld% mfl_mgm_w   (:,:,jb),      & ! mometum flux W-ward      (out)
                                mflux_s        = p_fld% mfl_mgm_s   (:,:,jb),      & ! mometum flux S-ward      (out)
                                mflux_n        = p_fld% mfl_mgm_n   (:,:,jb),      & ! mometum flux N-ward      (out)
                                amflux         = p_fld% amfl_mgm    (:,:,jb),      & ! absolute momentum flux   (out)
                                waflux_u       = p_fld% wafl_mgm_u  (:,:,jb),      & ! wave action flux U-ward  (out)
                                waflux_d       = p_fld% wafl_mgm_d  (:,:,jb),      & ! wave action flux D-ward  (out)
                                waflux_e       = p_fld% wafl_mgm_e  (:,:,jb),      & ! wave action flux E-ward  (out)
                                waflux_w       = p_fld% wafl_mgm_w  (:,:,jb),      & ! wave action flux W-ward  (out)
                                waflux_s       = p_fld% wafl_mgm_s  (:,:,jb),      & ! wave action flux S-ward  (out)
                                waflux_n       = p_fld% wafl_mgm_n  (:,:,jb),      & ! wave action flux N-ward  (out)
                                ptflux_e       = p_fld% ptfl_mgm_e  (:,:,jb),      & ! theta flux E-ward        (out)
                                ptflux_w       = p_fld% ptfl_mgm_w  (:,:,jb),      & ! theta flux W-ward        (out)
                                ptflux_s       = p_fld% ptfl_mgm_s  (:,:,jb),      & ! theta flux S-ward        (out)
                                ptflux_n       = p_fld% ptfl_mgm_n  (:,:,jb),      & ! theta flux N-ward        (out)
                                aptflux        = p_fld% aptfl_mgm   (:,:,jb),      & ! absolute pot temp flux   (out)
                                energy         = p_fld% energy_mgm  (:,:,jb),      & ! GW energy                (out)
                                energy_p       = p_fld% energy_p_mgm(:,:,jb),      & ! GW potential energy      (out)
                                waction        = p_fld% action_mgm_3(:,:,jb),      & ! GW action                (out)
                                active_rays    = p_fld% active_rays_mgm(:,jb)      ) ! number of rays per cell  (out)

      END IF

      !============================== DATOUT ================================
      ! Output profile of subgrid-scale GW quantities for specific lat,lon
      ! coordinates. Mostly used for idealized tests but can be useful for
      ! visualizing how ray volumes travel in vertical, etc.
      !======================================================================

      IF ( ldiagprof .AND. jg == 1 ) THEN

        IF (timers_level > 4) CALL timer_start(timer_msgwam_diagprof)

        IF ( ndiag_msgwam > 0 ) THEN
        IF ( ANY(jb_diag(1:ndiag_msgwam) == jb) ) THEN

        CALL debug_on

        DO idiag = 1,ndiag_msgwam

          IF ( jb /= jb_diag(idiag) ) CYCLE

          jc = jc_diag(idiag)

          WRITE(message_text,'(a,i8)') 'No. of existing rays before propagate:',&
                                        COUNT(p_ray(jg)%iexist(jc,:,jb)/=0)
          CALL message('', TRIM(message_text))

          CALL datout(nlev     = nlev,                         &
                      jg       = jg,                           &
                      idiag    = idiag,                        &
                      zray     = p_ray(jg)%z(jc,:,jb),         &
                      mray     = p_ray(jg)%m(jc,:,jb),         &
                      dzray    = p_ray(jg)%dz(jc,:,jb),        &
                      dmray    = p_ray(jg)%dm(jc,:,jb),        &
                      dens     = p_ray(jg)%wadens(jc,:,jb),    &
                      specid   = p_ray(jg)%specid(jc,:,jb),    &
                      zz_half  = p_metrics%z_ifc(jc,:,jb),     &
                      zz       = p_metrics%z_mc(jc,:,jb),      &
                      rho      = rho_half(jc,:,jb),            &
                      fld_u    = u(jc,:,jb),                   &
                      fld_v    = v(jc,:,jb),                   &
                      uw_wr    = p_fld% uwfl_mgm(jc,:,jb),     &
                      vw_wr    = p_fld% vwfl_mgm(jc,:,jb),     &
                      bvf2     = bvf2_half(jc,:,jb),           &
                      kd       = kd(jc,:),                     &
                      kd2      = kd2(jc,:),                    &
                      ld       = ld(jc,:),                     &
                      ld2      = ld2(jc,:),                    &
                      md       = md(jc,:),                     &
                      md2      = md2(jc,:),                    &
                      cgz_diag = cgz_diag(jc,:),               &
                      A        = A(jc,:),                      &
                      A_save   = A_save(jc,:),                 &
                      B2       = B2(jc,:),                     &
                      B2_save  = B2_save(jc,:),                &
                      mB2      = mB2(jc,:),                    &
                      mB2_save = mB2_save(jc,:),               &
                      ener     = p_fld% energy_mgm(jc,:,jb)    )

          ! Diagnostics for the wave breaking parametrization (saturation)
          IF (lsaturation) THEN
            DO jk = nlevp1,1,-1
              IF (mB2_save(jc,jk) > (alpha_sat*bvf2_half(jc,jk,jb))**2) THEN
                CALL message('', 'lsaturation=.T. ==> Saturation Parametrized')
                WRITE(message_text,'(a,3E12.4)') 'alpha, N2, height:', &
                                    alpha_sat, bvf2_half(jc,jk,jb), p_metrics%z_ifc(jc,jk,jb)
                CALL message('', TRIM(message_text))
                WRITE(message_text,'(a,3E12.4)') 'N^4*alpha^2, m^2B^2 ori, m^2B^2 reduced:', &
                                    (bvf2_half(jc,jk,jb)*alpha_sat)**2, mB2_save(jc,jk), mB2(jc,jk)
                CALL message('', TRIM(message_text))
                IF (mB2(jc,jk) > mB2_save(jc,jk)) THEN
                  WRITE(message_text,'(a,i4)') 'Warning: saturation did not work well for layer:', jk
                  CALL message('', TRIM(message_text))
                ENDIF
              ENDIF
            ENDDO
          ENDIF

          IF (msg_level >= 12 ) THEN
            DO jk = nlev,1,-1
              WRITE(message_text,'(a,i6,8E12.4)') 'jk, utflux, vtflux, '// &
                                                      'uuflux, uvflux, uwflux, '// &
                                                      'vvflux, vwflux:', &
                  jk, &
                  p_fld% utfl_mgm(jc,jk,jb), p_fld% vtfl_mgm(jc,jk,jb), &
                  p_fld% uufl_mgm(jc,jk,jb), p_fld% uvfl_mgm(jc,jk,jb), p_fld% uwfl_mgm(jc,jk,jb), &
                  p_fld% vvfl_mgm(jc,jk,jb), p_fld% vwfl_mgm(jc,jk,jb)
              CALL message('', TRIM(message_text))
            ENDDO
          ENDIF

        ENDDO  ! idiag

        CALL debug_off

        ENDIF ! ANY(jb_diag(1:ndiag_msgwam) == jb)
        ENDIF ! ndiag_msgwam > 0

        IF (timers_level > 4) CALL timer_stop(timer_msgwam_diagprof)

      ENDIF ! ldiagprof .AND. jg == 1


      !======================== PROPAGATE GW FIELD ==========================
      ! Ray equations are solved in 6D (x-k space) in order to calculate the
      ! propagation and deformation of ray volumes in phase-space.
      !======================================================================

      IF (timers_level > 4) CALL timer_start(timer_msgwam_propagate_wave)

          CALL propagate_wave(    dt_prop    = dt_call,                         & ! dt for prop. routine call  (in)
                                  jb         = jb,                              & ! block index (debug)        (in)
                                  nlev       = nlev,                            & ! no. of full levels         (in)
                                  i_startidx = i_startidx,                      & ! first index of the block   (in)
                                  i_endidx   = i_endidx,                        & ! last index of the block    (in)
                                  jray_start = 1,                               &
                                  jray_end   = nrays(jg),                       &
                                  z          = p_metrics%z_mc(:, :, jb),        & ! full level heights         (in)
                                  dz         = p_metrics%ddqz_z_full(:, :, jb), & ! full level layer thick.    (in)
                                  zhalf      = p_metrics%z_ifc(:, :, jb),       & ! half level heights         (in)
                                  dzhalf     = p_metrics%ddqz_z_half(:, :, jb), & ! half level layer thick.    (in)
                                  lat        = p_patch%cells%center(:, jb)%lat, & ! latitude of cell center    (in)
                                  lon        = p_patch%cells%center(:, jb)%lon, & ! longitude of cell center   (in)
                                  fc2        = fc2(:, jb),                      & ! Coriolis parameter**2      (in)
                                  fdfdlat    = fdfdlat(:, jb),                  & ! gradiend of fc2            (in)
                                  kvisc      = kvisc(:, :, jb),                 & ! Kinematic viscosity        (in)
                                  bvf2       = bvf2_half(:, :, jb),             & ! Brunt-Vaeisala freq**2      (in)
                                  gammash2   = gammash2_half(:, :, jb),         & ! inverse pinc scale height  (in)
                                  u          = u_fld(:, :, jb),                 & ! zonal wind                 (in)
                                  v          = v_fld(:, :, jb),                 & ! meridional wind            (in)
                                  dn2dz      = dn2dz(:, :, jb),                 & ! vertical grad of bvf2      (in)
                                  dg2dz      = dg2dz(:, :, jb),                 & ! vertical grad of gammash2  (in)
                                  dudz       = dudz (:, :, jb),                 & ! vertical grad of u         (in)
                                  dvdz       = dvdz (:, :, jb),                 & ! vertical grad of v         (in)
                                  n2_hgrad   = n2_hgrad(:, :, :, jb),           & ! horizontal grads of bvf**2 (in)
                                  g2_hgrad   = g2_hgrad(:, :, :, jb),           & ! horizontal grads of ga**2  (in)
                                  u_hgrad    = u_hgrad(:, :, :, jb),            & ! horizontal grads of u      (in)
                                  v_hgrad    = v_hgrad(:, :, :, jb),            & ! horizontal grads of v      (in)
                                  dn2dz_hgrad= dn2dz_hgrad(:, :, :, jb),        & ! vert & hori grads of bvf2  (in)
                                  dg2dz_hgrad= dg2dz_hgrad(:, :, :, jb),        & ! vert & hori grads of ga**2 (in)
                                  dudz_hgrad = dudz_hgrad (:, :, :, jb),        & ! vert & hori grads of u     (in)
                                  dvdz_hgrad = dvdz_hgrad (:, :, :, jb),        & ! vert & hori grads of v     (in)
                                  lonray     = p_ray(jg)%lon(:, :, jb),         & ! ray position lon           (inout)
                                  dxray      = p_ray(jg)%dx(:, :, jb),          & ! ray size in lon            (inout)
                                  latray     = p_ray(jg)%lat(:, :, jb),         & ! ray position lat           (inout)
                                  dyray      = p_ray(jg)%dy(:, :, jb),          & ! ray size in lat            (inout)
                                  zray       = p_ray(jg)%z(:, :, jb),           & ! position                   (inout)
                                  dzray      = p_ray(jg)%dz(:, :, jb),          & ! size in z-dir              (inout)
                                  coslatray  = p_ray(jg)%coslat(:, :, jb),      & ! cosine of ray lat          (inout)
                                  kray       = p_ray(jg)%k(:, :, jb),           & ! horiz (lon) wavenumber     (in)
                                  dkray      = p_ray(jg)%dk(:, :, jb),          & ! size in k-dir              (in)
                                  lray       = p_ray(jg)%l(:, :, jb),           & ! horiz (lat) wavenumber     (in)
                                  dlray      = p_ray(jg)%dl(:, :, jb),          & ! size in l-dir              (in)
                                  mray       = p_ray(jg)%m(:, :, jb),           & ! vertical wavenumber        (inout)
                                  dmray      = p_ray(jg)%dm(:, :, jb),          & ! size in m-dir              (inout)
                                  dens       = p_ray(jg)%wadens(:, :, jb),      & ! wave action density        (inout)
                                  iexist     = p_ray(jg)%iexist(:, :, jb),      & ! existence of ray           (inout)
                                  specid     = p_ray(jg)%specid(:, :, jb),      & ! spectral id of rays        (inout)
                                  jk_active  = p_ray(jg)%jk_active(:, :, jb),   & ! launch level index         (in)
                                  jr_last    = p_ray(jg)%jr_last(:, :, jb),     & ! last-launched ray index    (inout)
                                  jg=jg,&
                                  coef_sphere=coef_sphere, coef_torus=coef_torus)

      IF (timers_level > 4) CALL timer_stop(timer_msgwam_propagate_wave)

      !======================== Diagnostics ==================================
      ! WA diagnostic
      !=======================================================================

      CALL idx_rayedge( nlev,i_startidx,i_endidx,1,nrays(jg),& ! (in)
                        p_metrics%z_mc     (:,:,jb),        & ! (in)
                        p_metrics%z_ifc    (:,:,jb),        & ! (in)
                        p_ray(jg)%iexist   (:,:,jb),        & ! (in)
                        p_ray(jg)%specid   (:,:,jb),        & ! (in)
                        p_ray(jg)%jk_active(:,:,jb),        & ! (in)
                        p_ray(jg)%z        (:,:,jb),        & ! (in)
                        p_ray(jg)%dz       (:,:,jb),        & ! (in)
                        p_ray(jg)%jk_full_rtop(:,:,jb),     & ! (inout)
                        p_ray(jg)%jk_full_rbot(:,:,jb),     & ! (inout)
                        p_ray(jg)%jk_half_rtop(:,:,jb),     & ! (inout)
                        p_ray(jg)%jk_half_rbot(:,:,jb)      ) ! (inout)

      CALL project_action(nlev        = nlev,                             & !
                          i_startidx  = i_startidx,                       & !
                          i_endidx    = i_endidx,                         & !
                          jray_start  = 1,                                & !
                          jray_end    = nrays(jg),                        & !
                          z           = p_metrics%z_mc(:,:,jb),           & !
                          zhalf       = p_metrics%z_ifc(:,:,jb),          & !
                          cellarea    = p_patch%cells%area(:,jb),         & !
                          zray        = p_ray(jg)%z(:,:,jb),              & !
                          dzray       = p_ray(jg)%dz(:,:,jb),             & !
                          dyray       = p_ray(jg)%dy(:,:,jb),             & !
                          dxray       = p_ray(jg)%dx(:,:,jb),             & !
                          dkray       = p_ray(jg)%dk(:,:,jb),             & !
                          dlray       = p_ray(jg)%dl(:,:,jb),             & !
                          dmray       = p_ray(jg)%dm(:,:,jb),             & !
                          dens        = p_ray(jg)%wadens(:,:,jb),         & !
                          iexist      = p_ray(jg)%iexist(:,:,jb),         & !
                          specid      = p_ray(jg)%specid(:,:,jb),         & !
                          jk_active   = p_ray(jg)%jk_active(:,:,jb),      & !
                          jkmin_half  = p_ray(jg)%jk_half_rtop(:,:,jb),   & !
                          jkmax_half  = p_ray(jg)%jk_half_rbot(:,:,jb),   & !
                          action      = p_fld%action_mgm_4(:,:,jb),       & !
                          diag        = 0,                                & !
                          jg          = jg                                ) !

    ENDDO ! jb

!$OMP END DO
!$OMP END PARALLEL


    !================ AVERAGE GW FIELD IN HORIZONTAL ======================
    ! ...
    ! TODO:
    !======================================================================

#ifndef __msgwam1d
    IF (lhsmooth) THEN

      IF (timers_level > 4) CALL timer_start(timer_msgwam_smooth_hori)

      ! Background source related fields
      DO iter_sm = 1, ABS(nhsmooth)
        CALL sync_patch_array_mult(SYNC_C, p_patch, 5, lacc=.FALSE., f3din1=p_fld% uwfl_mgm, f3din2=p_fld% vwfl_mgm,  &
          &                           f3din3=p_fld% uufl_mgm, f3din4=p_fld% uvfl_mgm, f3din5=p_fld% vvfl_mgm)
        CALL sync_patch_array_mult(SYNC_C, p_patch, 5, lacc=.FALSE., f3din1=p_fld% uwpfl_mgm, f3din2=p_fld% vwpfl_mgm,  &
          &                          f3din3=p_fld% uupfl_mgm, f3din4=p_fld% uvpfl_mgm, f3din5=p_fld% vvpfl_mgm)
        CALL sync_patch_array_mult(SYNC_C, p_patch, 2, lacc=.FALSE., f3din1=p_fld% utfl_mgm, f3din2=p_fld% vtfl_mgm)
        !
        CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%uwfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%vwfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlev,  p_fld%uufl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlev,  p_fld%uvfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlev,  p_fld%vvfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%uwpfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%vwpfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlev,  p_fld%uupfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlev,  p_fld%uvpfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlev,  p_fld%vvpfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlev,  p_fld%utfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlev,  p_fld%vtfl_mgm)
      ENDDO

      IF (nhsmooth > 0) THEN
        ! Fluxes that do not affect the tendencies are smoothed only once or not ever,
        ! depending on the sign of nhsmooth.

        CALL sync_patch_array_mult(SYNC_C, p_patch, 5, lacc=.FALSE., f3din1=p_fld% apmfl_mgm    , &
          &                                            f3din2=p_fld% amfl_mgm     , &
          &                                            f3din3=p_fld% aptfl_mgm    , &
          &                                            f3din4=p_fld% energy_mgm   , &
          &                                            f3din5=p_fld% energy_p_mgm )
        !
        CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%apmfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%amfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlev,  p_fld%aptfl_mgm)
        CALL smooth_hori_wrapper(p_patch,nlev,  p_fld%energy_mgm)
        CALL smooth_hori_wrapper(p_patch,nlev,  p_fld%energy_p_mgm)

        IF ( lcalc_flux_4dir_bg(jg) ) THEN
          CALL sync_patch_array_mult(SYNC_C, p_patch, 4, lacc=.FALSE., f3din1=p_fld% mfl_mgm_e  , &
            &                                            f3din2=p_fld% mfl_mgm_w  , &
            &                                            f3din3=p_fld% mfl_mgm_n  , &
            &                                            f3din4=p_fld% mfl_mgm_s  )
          CALL sync_patch_array_mult(SYNC_C, p_patch, 4, lacc=.FALSE., f3din1=p_fld% pmfl_mgm_e  , &
            &                                            f3din2=p_fld% pmfl_mgm_w  , &
            &                                            f3din3=p_fld% pmfl_mgm_n  , &
            &                                            f3din4=p_fld% pmfl_mgm_s  )
          CALL sync_patch_array_mult(SYNC_C, p_patch, 4, lacc=.FALSE., f3din1=p_fld% ptfl_mgm_e , &
            &                                            f3din2=p_fld% ptfl_mgm_w , &
            &                                            f3din3=p_fld% ptfl_mgm_n , &
            &                                            f3din4=p_fld% ptfl_mgm_s )
          !
          CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%mfl_mgm_e)
          CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%mfl_mgm_w)
          CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%mfl_mgm_n)
          CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%mfl_mgm_s)
          CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%pmfl_mgm_e)
          CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%pmfl_mgm_w)
          CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%pmfl_mgm_n)
          CALL smooth_hori_wrapper(p_patch,nlevp1,p_fld%pmfl_mgm_s)
          CALL smooth_hori_wrapper(p_patch,nlev,p_fld%ptfl_mgm_e)
          CALL smooth_hori_wrapper(p_patch,nlev,p_fld%ptfl_mgm_w)
          CALL smooth_hori_wrapper(p_patch,nlev,p_fld%ptfl_mgm_n)
          CALL smooth_hori_wrapper(p_patch,nlev,p_fld%ptfl_mgm_s)
        ENDIF

      ENDIF  ! nhsmooth > 0

      IF (timers_level > 4) CALL timer_stop(timer_msgwam_smooth_hori)

    ENDIF ! lhsmooth
#endif


    !============================= TENDENCY ===============================
    ! Calculate horizontal wind tendencies based on the vertical divergence
    ! of pseudo-momentum fluxes. The total tendency and the tendency due to
    ! convective GWs is calculated separately (the latter for diagnostics).
    ! TODO: calculate the 3D divergence of momentum fluxes in case of
    !       horizontal propagation
    !======================================================================

    IF (timers_level > 4) CALL timer_start(timer_msgwam_tendency)

    CALL tendency(p_patch,p_metrics,p_int_state,rho,temp,theta,p_fld)

    IF (timers_level > 4) CALL timer_stop(timer_msgwam_tendency)


    !====================== SYNCHRONIZE GW FIELD ==========================
    ! ...
    ! TODO:
    !======================================================================

#ifndef __msgwam1d
    IF (timers_level > 4) CALL timer_start(timer_msgwam_sync_wave)
    CALL sync_wave(p_ray(jg),p_patch)
    IF (timers_level > 4) CALL timer_stop(timer_msgwam_sync_wave)
#endif


    !========================== SPLIT GW FIELD ============================
    ! ...
    ! TODO:
    !======================================================================

    ! inform user about the disabled merging
    IF ( imethod_merge > 0 ) THEN
      WRITE(message_text,'(a)') 'MS-GWaM: Ray volume merging is currently not available. Disabling.'
      CALL message('', TRIM(message_text))

      imethod_merge = 0
    END IF

    IF ( imethod_split > 0 ) THEN
      IF (timers_level > 4) CALL timer_start(timer_msgwam_split_merge)

        !TODO: only call split merge per source if source is enabled (performance issues)

        CALL split_rays( p_ray  = p_ray(jg),                        &  ! the ray properties     (INOUT)
              &          p_spl = p_spl(jg),                         &  ! info for the splitting (INOUT)
              &          p_rwork = p_rwork(jg),                     &  ! work array             (INOUT)
              &          jray_start = jray_offset_cv(jg)+1,         &  ! first ray index        (IN)
              &          gammash2 = gammash2_half,                  &  ! pink scale height      (IN)
              &          fc2 = fc2,                                 &  ! Coriolis parameter     (IN)
              &          bvf2 = bvf2_half,                          &  ! nuoyancy frequency     (IN)
              &          jray_end = jray_offset_cv(jg)+nrays_cv(jg), &   ! last ray index         (IN)
              &          nzray=p_patch%nlevp1, &
              &          jg=jg, start_block=p_patch%cells%start_block,&
              & end_block=p_patch%cells%end_block, nlev=p_patch%nlev)

        CALL split_rays( p_ray  = p_ray(jg),                        &  ! the ray properties     (INOUT)
              &          p_spl = p_spl(jg),                         &  ! info for the splitting (INOUT)
              &          p_rwork = p_rwork(jg),                     &  ! work array             (INOUT)
              &          gammash2 = gammash2_half,                  &  ! pink scale height      (IN)
              &          fc2 = fc2,                                 &  ! Coriolis parameter     (IN)
              &          bvf2 = bvf2_half,                          &  ! nuoyancy frequency     (IN)
              &          jray_start = jray_offset_bg(jg)+1,         &  ! first ray index        (IN)
              &          jray_end = jray_offset_bg(jg)+nrays_bg(jg), &  ! last ray index         (IN)
              &           nzray=p_patch%nlevp1, &
              & jg=jg, start_block=p_patch%cells%start_block, &
              &end_block=p_patch%cells%end_block, nlev=p_patch%nlev)


      IF (timers_level > 4) CALL timer_stop(timer_msgwam_split_merge)
    END IF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

      !======================== Diagnostics ==================================
      ! WA diagnostic
      !=======================================================================

      CALL idx_rayedge( nlev,i_startidx,i_endidx,1,nrays(jg),& ! (in)
                        p_metrics%z_mc     (:,:,jb),        & ! (in)
                        p_metrics%z_ifc    (:,:,jb),        & ! (in)
                        p_ray(jg)%iexist   (:,:,jb),        & ! (in)
                        p_ray(jg)%specid   (:,:,jb),        & ! (in)
                        p_ray(jg)%jk_active(:,:,jb),        & ! (in)
                        p_ray(jg)%z        (:,:,jb),        & ! (in)
                        p_ray(jg)%dz       (:,:,jb),        & ! (in)
                        p_ray(jg)%jk_full_rtop(:,:,jb),     & ! (inout)
                        p_ray(jg)%jk_full_rbot(:,:,jb),     & ! (inout)
                        p_ray(jg)%jk_half_rtop(:,:,jb),     & ! (inout)
                        p_ray(jg)%jk_half_rbot(:,:,jb)      ) ! (inout)

      CALL project_action(nlev        = nlev,                             & !
                          i_startidx  = i_startidx,                       & !
                          i_endidx    = i_endidx,                         & !
                          jray_start  = 1,                                & !
                          jray_end    = nrays(jg),                        & !
                          z           = p_metrics%z_mc(:,:,jb),           & !
                          zhalf       = p_metrics%z_ifc(:,:,jb),          & !
                          cellarea    = p_patch%cells%area(:,jb),         & !
                          zray        = p_ray(jg)%z(:,:,jb),              & !
                          dzray       = p_ray(jg)%dz(:,:,jb),             & !
                          dxray       = p_ray(jg)%dy(:,:,jb),             & !
                          dyray       = p_ray(jg)%dx(:,:,jb),             & !
                          dkray       = p_ray(jg)%dk(:,:,jb),             & !
                          dlray       = p_ray(jg)%dl(:,:,jb),             & !
                          dmray       = p_ray(jg)%dm(:,:,jb),             & !
                          dens        = p_ray(jg)%wadens(:,:,jb),         & !
                          iexist      = p_ray(jg)%iexist(:,:,jb),         & !
                          specid      = p_ray(jg)%specid(:,:,jb),         & !
                          jk_active   = p_ray(jg)%jk_active(:,:,jb),      & !
                          jkmin_half  = p_ray(jg)%jk_half_rtop(:,:,jb),   & !
                          jkmax_half  = p_ray(jg)%jk_half_rbot(:,:,jb),   & !
                          action      = p_fld%action_mgm_5(:,:,jb),       & !
                          diag        = 0,                                & !
                          jg          = jg                                ) !

    ENDDO  ! jb loop for diagnostic
!$OMP END DO
!$OMP END PARALLEL

    !======================== RE-GRID GW FIELD ============================
    ! ...
    ! TODO:
    !======================================================================

#ifndef __msgwam1d
    IF (timers_level > 4) CALL timer_start(timer_msgwam_split_merge)

    CALL regrid_wave( p_patch=p_patch,                          & ! grid/patch info. (inout)
      &               p_gridinfo4ray=p_gridinfo4ray(jg),        & ! infos on grid cell and neighbors(in)
      &               z_ifc=p_metrics%z_ifc,                    & ! cell interface heights(in)
      &               z_mc=p_metrics%z_mc,                      & ! cell center heights (in)
      &               gammash2=gammash2_full,                   & ! pink scale height squared (in)
      &               fc2=fc2,                                  & ! Coriolis parameter squared (in)
      &               bvf2=bvf2_full,                           & ! buoyancy frequency squared (in)
      &               jray_start=jray_offset_bg(jg)+1,          & ! jray lower limit (in)
      &               jray_end=jray_offset_bg(jg)+nrays_bg(jg), & ! jray upper limit (in)
      &               p_ray=p_ray(jg)                           ) ! ray properties (inout)

    CALL regrid_wave( p_patch=p_patch,                          & ! grid/patch info. (inout)
      &               p_gridinfo4ray=p_gridinfo4ray(jg),        & ! infos on grid cell and neighbors(in)
      &               z_ifc=p_metrics%z_ifc,                    & ! cell interface heights(in)
      &               z_mc=p_metrics%z_mc,                      & ! cell center heights (in)
      &               gammash2=gammash2_full,                   & ! pink scale height squared (in)
      &               fc2=fc2,                                  & ! Coriolis parameter squared (in)
      &               bvf2=bvf2_full,                           & ! buoyancy frequency squared (in)
      &               jray_start=jray_offset_cv(jg)+1,          & ! jray lower limit (in)
      &               jray_end=jray_offset_cv(jg)+nrays_cv(jg), & ! jray upper limit (in)
      &               p_ray=p_ray(jg)                           ) ! ray properties (inout)

    IF (timers_level > 4) CALL timer_stop(timer_msgwam_split_merge)
#endif

    !======================== TENDENCY LIMITER ============================
    ! Tendency limiter to stabilize high-top runs. This limiter is taken
    ! from mo_nwp_gw_interface.f90. Normally we would not like to use it as
    ! MS-GWaM fluxes with a direct wave-meanflow interaction (+ wave
    ! breaking) should not be out of realistic range. We still keep the
    ! option...
    !======================================================================

    ! hard-coded tendency checker :  to be removed later
    CALL test_blowup_uvt(p_patch%id, p_fld%ddt_u_gwd_mgm, p_fld%ddt_v_gwd_mgm, p_fld%ddt_t_gwd_mgm,  &
      &                  lim_tend_tmp, lim_tend_tmp*0.1_wp, &
      &                  p_patch%nlev, p_patch%cells%start_block, p_patch%cells%end_block, &
      &                  p_patch%cells%center(:,:)%lat,p_patch%cells%center(:,:)%lon,'(checkpoint 1)')

    IF (llimittend) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb)
      DO jb = i_startblk, i_endblk

        p_fld%ddt_u_gwd_mgm(:,:,jb) = MAX(-lim_tend_tmp, MIN(lim_tend_tmp,  &
          &                           p_fld%ddt_u_gwd_mgm(:,:,jb)))
        p_fld%ddt_v_gwd_mgm(:,:,jb) = MAX(-lim_tend_tmp, MIN(lim_tend_tmp,  &
          &                           p_fld%ddt_v_gwd_mgm(:,:,jb)))
        p_fld%ddt_u_gwd_pmom_mgm(:,:,jb) = MAX(-lim_tend_tmp, MIN(lim_tend_tmp,  &
          &                                p_fld%ddt_u_gwd_pmom_mgm(:,:,jb)))
        p_fld%ddt_v_gwd_pmom_mgm(:,:,jb) = MAX(-lim_tend_tmp, MIN(lim_tend_tmp,  &
          &                                p_fld%ddt_v_gwd_pmom_mgm(:,:,jb)))

      ENDDO ! jb
!$OMP END DO
!$OMP END PARALLEL

    ENDIF ! llimittend


    !======================== TENDENCY DIAGNOSTIC ============================
    ! Print out maximum of absolute of tendencies
    ! TODO: one would need a syncronization and an mpi_allreduce operation
    !       here in order to get full domain diagnostics. This is here only a
    !       preliminary printout to have a first feeling how often the
    !       tendencies are limited via llimittend = .true.
    !=========================================================================
    IF (msg_level >= 12) THEN
      maxindexu = MAXLOC(ABS(p_fld%ddt_u_gwd_mgm))
      maxindexv = MAXLOC(ABS(p_fld%ddt_v_gwd_mgm))
      WRITE(message_text,'(a,E12.4,a,i4)') 'MAXABS of GW u-tendencies:', &
        &  MAXVAL(ABS(p_fld% ddt_u_gwd_mgm)), ' at level:', maxindexu(2)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,E12.4,a,i4)') 'MAXABS of GW v-tendencies:', &
        &  MAXVAL(ABS(p_fld% ddt_v_gwd_mgm)), ' at level:', maxindexv(2)
      CALL message('', TRIM(message_text))
    ENDIF

    ! Record number for subroutine datout (ldiagprof=.T.)
    iout_msgwam = iout_msgwam + 1

#ifndef __msgwam1d
    CALL message('gwdrag_msgwam', 'Transient MS-GWaM finished')
#else
    CALL message('gwdrag_msgwam', 'Transient MS-GWaM-1D finished')
#endif

  END SUBROUTINE gwdrag_msgwam
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE smooth_hori_wrapper( p_patch,        & !inout
    &                     nlev,           & !in
    &                     var             ) !inout

    ! In/out variables
    TYPE(t_patch)  , TARGET, INTENT(in)    ::  p_patch
    INTEGER,              INTENT(IN)    :: nlev
    REAL(wp),             INTENT(INOUT) :: var(:,:,:)

    CALL smooth_hori(p_patch%id,nlev,var, &
            & p_patch%nblks_c, p_patch%cells%start_block, p_patch%cells%end_block,&
            & p_patch%cells%vertex_idx, p_patch%cells%vertex_blk, p_patch%verts%cell_idx, &
            & p_patch%verts%cell_blk, p_patch%verts%num_edges)
  END SUBROUTINE smooth_hori_wrapper
#endif
END MODULE mo_nwp_msgwam_interface
