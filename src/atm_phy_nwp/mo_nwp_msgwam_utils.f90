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
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_nwp_msgwam_utils
#ifdef __MSGWAM

  USE mo_kind,                   ONLY: wp, vp, sp
  USE mo_mpi,                    ONLY: my_process_is_stdio, p_wait, work_mpi_barrier
  USE mo_exception,              ONLY: message, finish, message_text, &
                                      debug_on, debug_off
  USE mo_model_domain,           ONLY: t_patch
  USE mo_nonhydro_types,         ONLY: t_nh_diag, t_nh_prog, t_nh_metrics
  USE mo_intp_data_strc,         ONLY: t_int_state
  USE mo_nwp_phy_types,          ONLY: t_nwp_phy_tend
  USE mo_impl_constants,         ONLY: min_rlcell_int, min_rledge_int, success
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,            ONLY: get_indices_c, get_indices_e
  USE mo_physical_constants,     ONLY: grav, rd, cpd, cvd, rd_o_cpd, earth_angular_velocity
  USE mo_grid_config,            ONLY: grid_sphere_radius, is_plane_torus
  USE mo_sync,                   ONLY: sync_patch_array, sync_patch_array_mult, SYNC_E, &
                                      SYNC_C, SYNC_C1
  USE mo_parallel_config,        ONLY: nproma
  USE mo_run_config,             ONLY: msg_level, ldynamics
  USE mo_dynamics_config,        ONLY: lcoriolis
  USE mo_math_constants,         ONLY: pi, pi2, pi_2, rad2deg, deg2rad
  USE mo_math_gradients,         ONLY: grad_green_gauss_cell
  USE mo_util_vgrid_types,       ONLY: vgrid_buffer
  USE mtime,                     ONLY: datetime, timeDelta, newTimedelta, &
                                      deallocateTimedelta, getTimedeltaFromDatetime, &
                                      getTotalMillisecondsTimedelta
  USE mo_vertical_grid,          ONLY: nrdmax
  USE mo_intp,                   ONLY: cell_avg
  USE mo_intp_data_strc,         ONLY: t_int_state
  USE mo_intp_lonlat_baryctr,    ONLY: inside_triangle
  USE mo_delaunay_types,         ONLY: t_point
  USE mo_math_divrot,            ONLY: div, div_avg
  USE mo_nh_testcases_nml,       ONLY: nh_test_name
  USE mo_msgwam_config,          ONLY: bvf2_min, facgamma, imethod_merge, l1ray, &
                                      lcorrlongwaves, lmvisc, lsmoothb, lsmootht, &
                                      ltest_hprop, msgw_lower_bound_opt, nlev_lbnd, &
                                      nlevsmootht, nrays, nrays_coll, nsmooth
  USE mo_setup_msgwam_interface, ONLY: filename_ltest_hprop, iout_msgwam, msgwam_write_restartfiles, nlaunch_max_bg, &
                                      nlaunch_max_cv, p_rwork, t_msgwam, t_ray, p_ray, p_gridinfo4ray, &
                                      ndiag_msgwam, jb_diag, p_msgwam, p_mgmgrid, t_gridinfo4ray, &
                                      jc_diag, p_spl, p_lfluxbg,ray_coll
  USE mo_msgwam_util,            ONLY: dyn_visc_sutherland, get_jkhalf_closest_gen, heuristic_check_3d, &
                                      idx_rayedge, inside_triangle_torus, smooth_vert
  USE mo_gw_source_config,       ONLY: gws_conv_config
  USE mo_gw_source_bg,           ONLY: init_gw_orretal
  USE mo_gw_source_conv,         ONLY: init_gw_conv
  USE mo_msgwam_diagnostics,     ONLY: project_action, project_diagnostics
  USE mo_msgwam_splitmerge,      ONLY: split_rays, remove_rays, remove_rays_in_cell
  USE mo_msgwam,                 ONLY: saturation, wave2grid


  ! Modules for computing the Orr et al. 2010 launch spectrum
  USE data_gwd,    ONLY : nslope, gfluxlaun, &
    &                     ggaussa, ggaussb, ngauss, gcoeff, lozpr

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: fieldsgrads, regrid_wave, sync_wave, tendency

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
  SUBROUTINE fieldsgrads( p_patch, p_int_state, rl_start, rl_end, nlev, dz, wgtfac_c, &
                          rho, temp, u, v, fc2, fdfdlat, rho_half, u_fld, v_fld, &
                          bvf2_full, bvf2_half, gammash2_full, gammash2_half, kvisc, &
                          dn2dz, dg2dz, dudz, dvdz, &
                          n2_hgrad, g2_hgrad, u_hgrad, v_hgrad, &
                          dn2dz_hgrad, dg2dz_hgrad, dudz_hgrad, dvdz_hgrad )

    TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch                 ! grid/patch info.
    TYPE(t_int_state),     INTENT(IN)    :: p_int_state             ! interpolation struct
    INTEGER,               INTENT(IN)    :: rl_start, rl_end        ! iteration bounds
    INTEGER,               INTENT(IN)    :: nlev                    ! number of full levels
    REAL(wp),              INTENT(IN)    :: dz(:, :, :)             ! thickness of full levels
    REAL(vp),              INTENT(IN)    :: wgtfac_c(:, :, :)       !
    REAL(wp),              INTENT(IN)    :: rho(:, :, :)            ! density at full levels
    REAL(wp),              INTENT(IN)    :: temp(:, :, :)           ! temperature at full levels
    REAL(wp),              INTENT(IN)    :: u(:, :, :)              ! zonal wind (full levels)
    REAL(wp),              INTENT(IN)    :: v(:, :, :)              ! meridional wind (full levels)
    REAL(wp),              INTENT(  OUT) :: fc2(:, :)               ! fc^2
    REAL(wp),              INTENT(  OUT) :: fdfdlat(:, :)           ! d(fc^2)/dlat
    REAL(wp),              INTENT(  OUT) :: rho_half(:, :, :)       ! density at half levels
    REAL(wp),              INTENT(  OUT) :: u_fld(:, :, :)          ! zonal wind
    REAL(wp),              INTENT(  OUT) :: v_fld(:, :, :)          ! meridional wind
    REAL(wp),              INTENT(  OUT) :: bvf2_full(:, :, :)      ! Brunt-Vaeisaelae freq**2 at full levels
    REAL(wp),              INTENT(  OUT) :: bvf2_half(:, :, :)      ! Brunt-Vaeisaelae freq**2 at half levels
    REAL(wp),              INTENT(  OUT) :: gammash2_full(:, :, :)  ! inverse pinc scale height**2 (i.e. gamma**2) at full levels
    REAL(wp),              INTENT(  OUT) :: gammash2_half(:, :, :)  ! inverse pinc scale height**2 (i.e. gamma**2) at half levels
    REAL(wp),              INTENT(  OUT) :: kvisc(:, :, :)          ! Kinematic viscosity
    REAL(wp),              INTENT(  OUT) :: dn2dz(:, :, :)          ! vertical B-V freq gradients at half levels
    REAL(wp),              INTENT(  OUT) :: dg2dz(:, :, :)          ! vertical gradients of inverse pinc scale height at half levels
    REAL(wp),              INTENT(  OUT) :: dudz (:, :, :)          ! vertical (radial) wind u gradients at half levels
    REAL(wp),              INTENT(  OUT) :: dvdz (:, :, :)          ! vertical (radial) wind v gradients at half levels
    REAL(wp),              INTENT(  OUT) :: n2_hgrad(:, :, :, :)    ! zonal / meridional gradients of BVF**2
    REAL(wp),              INTENT(  OUT) :: g2_hgrad(:, :, :, :)    ! zonal / meridional gradients of gamma**2 at full levels
    REAL(wp),              INTENT(  OUT) :: u_hgrad(:, :, :, :)     ! zonal / meridional wind u gradients at full levels
    REAL(wp),              INTENT(  OUT) :: v_hgrad(:, :, :, :)     ! zonal / meridional wind v gradients at full levels
    REAL(wp),              INTENT(  OUT) :: dn2dz_hgrad(:, :, :, :) ! vert & zonal / meridional gradients of BVF**2
    REAL(wp),              INTENT(  OUT) :: dg2dz_hgrad(:, :, :, :) ! vert & zonal / meridional gradients of gamma**2 at full levels
    REAL(wp),              INTENT(  OUT) :: dudz_hgrad(:, :, :, :)  ! vert & zonal / meridional wind u gradients at full levels
    REAL(wp),              INTENT(  OUT) :: dvdz_hgrad(:, :, :, :)  ! vert & zonal / meridional wind v gradients at full levels

    REAL(wp)                             :: temp_half(nproma, p_patch%nlevp1)
    REAL(vp)                             :: hgrad(2, nproma, p_patch%nlevp1, p_patch%nblks_c)
    INTEGER                              :: jk, jc                  ! vertical and block indices
    INTEGER                              :: i_startidx, i_endidx    ! first and last block index
    INTEGER                              :: i_startblk, i_endblk    ! first and last block index
    INTEGER                              :: rl_end_h1               ! iteration bound with 1 halo line
    INTEGER                              :: i_endblk_h1             ! last block index including 1 halo line
    INTEGER                              :: jb                      ! block index
    INTEGER                              :: nlevp1
    REAL(wp)                             :: inv_dz
    REAL(wp)                             :: wgtfac_c_jkm1           ! 1. - wgtfac_c
    REAL(wp)                             :: Hrho                    ! density scale height
    REAL(wp),    PARAMETER               :: grav_o_cpd = grav / cpd !
    INTEGER  ::  jg
      ! Domain ID
    jg = p_patch%id

    IF (msg_level >= 12) CALL message('fieldsgrads', 'MS-GWaM: prepare resolved fields and grads')

    !----------------------------------------------------------------------
    ! Purpose:
    !         Calculate all necessary resolved fields appearing as input to MS-GWaM
    ! Method:
    !         -- smooth wind and temperature (OMP jb loop)
    !         -- calculate Brunt-Vaeisaelae frequency (OMP jb loop)
    !         -- calculate pinc scale height correction term gamma (OMP jb loop)
    !         -- calculate vertcial gradients (OMP jb loop)
    !         -- calculate the horizontal wind gradients (Green-Gauss method)
    !
    !----------------------------------------------------------------------

    ! get indices of first and last block
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

#ifdef __msgwam1d
    rl_end_h1   = rl_end
    i_endblk_h1 = i_endblk
#else
    rl_end_h1   = rl_end - 1   ! to include one inner line of halo
  ! rl_end_h1   = rl_end - 2   ! to include one pair (inner + outer) of halo line
    i_endblk_h1 = p_patch%cells%end_block(rl_end_h1)
#endif

    ! set number of vertical levels + 1
    nlevp1 = nlev+1

    IF (is_plane_torus)  fdfdlat(:, :) = 0._wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,wgtfac_c_jkm1,Hrho,temp_half)

    DO jb = i_startblk, i_endblk

      ! get cell indices
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

      ! Initialize Coriolis frequency
      fc2(:, jb) = p_patch%cells%f_c(:, jb)**2
      IF (.NOT. is_plane_torus)  fdfdlat(:, jb) = (4._wp*earth_angular_velocity)  &
        &            * p_patch%cells%f_c(:, jb) * p_mgmgrid(jg)%coslat(:, jb)

      ! Interpolate the resolved fields to half levels
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
          wgtfac_c_jkm1 = 1._wp - wgtfac_c(jc, jk, jb)
          rho_half (jc, jk, jb) = wgtfac_c(jc, jk, jb) * rho (jc, jk ,jb) &
                                      + wgtfac_c_jkm1 * rho (jc, jk - 1, jb)
          temp_half(jc, jk) = wgtfac_c(jc, jk, jb) * temp(jc, jk, jb) &
                                      + wgtfac_c_jkm1 * temp(jc, jk - 1, jb)
          u_fld(jc, jk, jb) = wgtfac_c(jc, jk, jb) * u(jc, jk, jb) &
                                  + wgtfac_c_jkm1 * u(jc, jk - 1, jb)
          v_fld(jc, jk, jb) = wgtfac_c(jc, jk, jb) * v(jc, jk, jb) &
                                  + wgtfac_c_jkm1 * v(jc, jk - 1, jb)
        ENDDO  ! jc
      ENDDO  ! jk

      DO jc = i_startidx, i_endidx
        rho_half (jc, 1, jb) = rho (jc, 1, jb)
        rho_half (jc, nlevp1, jb) = rho (jc, nlev, jb)
        temp_half(jc, 1) = temp(jc, 1, jb)
        temp_half(jc, nlevp1) = temp(jc, nlev, jb)
      ENDDO

      DO jc = i_startidx, i_endidx
        ! dA/dz = 0 to prevent extrapolation
        u_fld(jc, 1     , jb) = u_fld(jc, 2   , jb)
        u_fld(jc, nlevp1, jb) = u_fld(jc, nlev, jb)
        v_fld(jc, 1     , jb) = v_fld(jc, 2   , jb)
        v_fld(jc, nlevp1, jb) = v_fld(jc, nlev, jb)
      ENDDO

      ! Smoothing over 2*nsmooth+1 points
      IF (lsmoothb) THEN
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .TRUE., temp_half(:, :))
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .TRUE., u_fld(:, :, jb))
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .TRUE., v_fld(:, :, jb))
        DO jc = i_startidx, i_endidx
          ! dA/dz = 0 to prevent extrapolation
          u_fld(jc, 1     , jb) = u_fld(jc, 2   , jb)
          u_fld(jc, nlevp1, jb) = u_fld(jc, nlev, jb)
          v_fld(jc, 1     , jb) = v_fld(jc, 2   , jb)
          v_fld(jc, nlevp1, jb) = v_fld(jc, nlev, jb)
        ENDDO
      ENDIF

      ! Fields on full levels
      DO jk = 2, nlev - 1

          ! Diag printout
          IF (msg_level >= 15) THEN
            WRITE(message_text, '(a,3i6,E12.4, a,3i6,E12.4, a,3i6,E12.4)') &
          'jc, jk, jb, temp(jc,jk,jb)            :', i_startidx, jk, jb, temp(i_startidx,jk,jb), &
          ' ;           ', (i_endidx)/2, jk, jb, temp((i_endidx)/2,jk,jb), &
          ' ;           ', i_endidx, jk, jb, temp(i_endidx,jk,jb)
            CALL message('', TRIM(message_text))
          ENDIF

        DO jc = i_startidx, i_endidx

          ! Calculate Brunt-Vaeisaelae freq (full levels)
          ! negative pot temperature gradient is not allowed,
          ! i.e. in that case bvf2 is set to bvf2_min
  !       bvf2_full(jc, jk, jb) = MAX(bvf2_min, grav / temp(jc, jk, jb) * (grav_o_cpd + &
  !             (temp(jc, jk - 1, jb) - temp(jc, jk + 1, jb)) / (dzhalf(jc, jk, jb) + dzhalf(jc, jk + 1, jb))))
          bvf2_full(jc, jk, jb) = MAX(bvf2_min, grav / temp(jc, jk, jb) * (grav_o_cpd + &
                (temp_half(jc, jk) - temp_half(jc, jk + 1)) / dz(jc, jk, jb)))

          ! Calculate gamma (half levels): inverse pseudo incompressible scaleheight
          IF (lcorrlongwaves) THEN
            ! Density scaleheight assuming "locally" isothermal atmosphere
            Hrho = rd * temp(jc, jk, jb) / grav
            ! Assuming "locally" isothermal atmosphere
            gammash2_full(jc, jk, jb) = (facgamma / Hrho)**2
            ! facgamma is a namelist parameter with the following typical values:
            ! Pseudo-incompressible correction:         facgamma ~= 0.214
            ! Anelastic correction:                     facgamma  = 0.5
            ! Further decrease effects from long waves: facgamma  > 0.5
          ENDIF

        ENDDO ! jc
        ! Diag printout
        IF (lcorrlongwaves .AND. msg_level>= 15) THEN
              WRITE(message_text, '(a,3i6,2E12.4,a,3i6,2E12.4,a,3i6,2E12.4)') &
                'jc, jk, jb, bvf2_full, 1/gammash2_full:', &
                i_startidx, jk, jb, bvf2_full(i_startidx,jk,jb), 1._wp/gammash2_full(i_startidx,jk,jb), &
                ' ', (i_endidx)/2, jk, jb, bvf2_full((i_endidx)/2,jk,jb), 1._wp/gammash2_full((i_endidx)/2,jk,jb), &
                ' ', i_endidx, jk, jb, bvf2_full(i_endidx,jk,jb), 1._wp/gammash2_full(i_endidx,jk,jb)
              CALL message('', TRIM(message_text))
            ENDIF

      ENDDO ! jk

      ! fields on half levels
      DO jk = 2, nlev

        ! Diag printout
        IF (msg_level >= 15) THEN
          WRITE(message_text, '(a,3i6,E12.4,a,3i6,E12.4,a,3i6,E12.4)') 'jc, jk, jb, temp_half(jc, jk):', &
                i_startidx, jk, jb, temp_half(i_startidx, jk), i_endidx/2, jk, jb, temp_half(i_endidx/2, jk),&
                i_endidx,jk, jb, temp_half(i_endidx, jk)
          CALL message('', TRIM(message_text))
        ENDIF

        DO jc = i_startidx, i_endidx

          ! Calculate Brunt-Vaeisaelae freq (half levels)
          ! negative pot temperature gradient is not allowed,
          ! i.e. in that case bvf2 is set to bvf2_min
          bvf2_half(jc, jk, jb) = MAX(bvf2_min, grav / temp_half(jc, jk) * (grav_o_cpd + &
                (temp_half(jc, jk - 1) - temp_half(jc, jk + 1)) / (dz(jc, jk - 1, jb) + dz(jc, jk, jb))))

          ! Calculate gamma (half levels): inverse pseudo incompressible scaleheight
          IF (lcorrlongwaves) THEN
            ! Density scaleheight assuming "locally" isothermal atmosphere
            Hrho = rd * temp_half(jc, jk) / grav
            ! Assuming "locally" isothermal atmosphere
            gammash2_half(jc, jk, jb) = (facgamma / Hrho)**2
            ! facgamma is a namelist parameter with the following typical values:
            ! Pseudo-incompressible correction:         facgamma ~= 0.214
            ! Anelastic correction:                     facgamma  = 0.5
            ! Further decrease effects from long waves: facgamma  > 0.5
          ENDIF

          ! Calculate kinematic viscosity profile
          kvisc(jc, jk, jb) = dyn_visc_sutherland(temp_half(jc, jk)) / rho_half(jc, jk, jb)

        ENDDO ! jc

          ! Diag printout
        IF (msg_level >= 15) THEN
          WRITE(message_text, '(a,3i6,2E12.4,a,3i6,2E12.4,a,3i6,2E12.4)') 'jc, jk, jb, bvf2_half(jc, jk, jb), 1 / gammash2_half(jc, jk, jb):', &
                i_startidx, jk, jb, bvf2_half(i_startidx, jk, jb), 1._wp / gammash2_half(i_startidx, jk, jb), &
                i_endidx/2, jk, jb, bvf2_half(i_endidx/2, jk, jb), 1._wp / gammash2_half(i_endidx/2, jk, jb), &
                i_endidx, jk, jb, bvf2_half(i_endidx, jk, jb), 1._wp / gammash2_half(i_endidx, jk, jb)
          CALL message('', TRIM(message_text))
        ENDIF

      ENDDO ! jk

      IF (.NOT. lmvisc)  kvisc(:,:,jb) = 0.

      ! Values at vertical boundaries
      DO jc = i_startidx, i_endidx
        ! dA/dz = 0 to prevent extrapolation
        bvf2_half(jc, 1, jb)          = bvf2_half(jc, 2, jb)
        bvf2_half(jc, nlevp1, jb)     = bvf2_half(jc, nlev, jb)
        gammash2_half(jc, 1, jb)      = gammash2_half(jc, 2, jb)
        gammash2_half(jc, nlevp1, jb) = gammash2_half(jc, nlev, jb)

        bvf2_full(jc, 1, jb)        = bvf2_full(jc, 2, jb)
        bvf2_full(jc, nlev, jb)     = bvf2_full(jc, nlev-1, jb)
        gammash2_full(jc, 1, jb)    = gammash2_full(jc, 2, jb)
        gammash2_full(jc, nlev, jb) = gammash2_full(jc, nlev-1, jb)

        kvisc(jc, 1, jb)    = kvisc(jc, 2, jb)
        kvisc(jc, nlevp1, jb) = kvisc(jc, nlev, jb)
      ENDDO ! jc

    ENDDO ! jb

!$OMP END DO
!$OMP END PARALLEL

#ifndef __msgwam1d
    ! Synchronize input fields of horizontal gradient calculations
    CALL sync_patch_array_mult(SYNC_C1, p_patch, 4, lacc=.FALSE., f3din1=u_fld, f3din2=v_fld, f3din3=bvf2_half, f3din4=gammash2_half)
#endif

    ! The vertical gradients calculated below will be used in the equation for Dm/Dt,
    ! after interpolation to the ray-volume center.
    ! It has turned out that the use of vertically interpolated gradient fields
    ! helps to obtain numerical stable results, compared to the use of gradient fields
    ! directly calculated between the top and bottom of a ray volume.

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,inv_dz)
    DO jb = i_startblk, i_endblk_h1   ! including one inner line of halo, assuming that sync has
                                      ! already been made above with SYNC_C*, for the 4 variables
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk_h1, &
            & i_startidx, i_endidx, rl_start, rl_end_h1)

      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx

          inv_dz = 1._wp / (dz(jc, jk-1, jb) + dz(jc, jk, jb))  ! two full layers

          ! vertical gradients of N^2, Gamma^2
          dn2dz(jc, jk, jb) = (bvf2_half(jc, jk-1, jb) - bvf2_half(jc, jk+1, jb)) * inv_dz
          dg2dz(jc, jk, jb) = (gammash2_half(jc, jk-1, jb) - gammash2_half(jc, jk+1, jb)) * inv_dz

          ! vertical (radial) gradients of the horizontal wind u, v
          dudz(jc, jk, jb) = (u_fld(jc, jk-1, jb) - u_fld(jc, jk+1, jb)) * inv_dz
          dvdz(jc, jk, jb) = (v_fld(jc, jk-1, jb) - v_fld(jc, jk+1, jb)) * inv_dz

        ENDDO ! jc
      ENDDO ! jk

      ! Values at vertical boundaries
      dn2dz(:, 1, jb) = 0._wp  ;  dn2dz(:, nlevp1, jb) = 0._wp
      dg2dz(:, 1, jb) = 0._wp  ;  dg2dz(:, nlevp1, jb) = 0._wp
      dudz (:, 1, jb) = 0._wp  ;  dudz (:, nlevp1, jb) = 0._wp
      dvdz (:, 1, jb) = 0._wp  ;  dvdz (:, nlevp1, jb) = 0._wp

    ENDDO ! jb
!$OMP END DO
!$OMP END PARALLEL

#ifndef __msgwam1d
    ! get the gradients in lat and lon
    ! the green gauss algorithms are parallelized with OMP internally
    hgrad(:,:,:,:) = 0.
    CALL grad_green_gauss_cell(p_cc = u_fld,           & ! cell centered variable to interpolate
                              ptr_patch = p_patch,    & ! patch info
                              ptr_int = p_int_state,  & ! interpolation info
                              p_grad = hgrad,         & ! horizontal gradients
                              lacc=.FALSE.,           & ! OpenAcc flag
                              opt_slev = 2,           & ! lowest vertical level
                              opt_elev = nlev,        & ! highest vertical level
                              opt_rlstart = rl_start, &
                              opt_rlend   = rl_end    )
    u_hgrad = REAL(hgrad, KIND=wp)

    hgrad(:,:,:,:) = 0.
    CALL grad_green_gauss_cell(p_cc = v_fld,           & ! cell centered variable to interpolate
                              ptr_patch = p_patch,    & ! patch info
                              ptr_int = p_int_state,  & ! interpolation info
                              p_grad = hgrad,         & ! horizontal gradients
                              lacc=.FALSE.,           & ! OpenAcc flag
                              opt_slev = 2,           & ! lowest vertical level
                              opt_elev = nlev,        & ! highest vertical level
                              opt_rlstart = rl_start, &
                              opt_rlend   = rl_end    )
    v_hgrad = REAL(hgrad, KIND=wp)

    hgrad(:,:,:,:) = 0.
    CALL grad_green_gauss_cell(p_cc = bvf2_half,       & ! cell centered variable to interpolate
                              ptr_patch = p_patch,    & ! patch info
                              ptr_int = p_int_state,  & ! interpolation info
                              p_grad = hgrad,         & ! horizontal gradients
                              lacc=.FALSE.,           & ! OpenAcc flag
                              opt_slev = 2,           & ! lowest vertical level
                              opt_elev = nlev,        & ! highest vertical level
                              opt_rlstart = rl_start, &
                              opt_rlend   = rl_end    )
    n2_hgrad = REAL(hgrad, KIND=wp)

    hgrad(:,:,:,:) = 0.
    CALL grad_green_gauss_cell(p_cc = gammash2_half,   & ! cell centered variable to interpolate
                              ptr_patch = p_patch,    & ! patch info
                              ptr_int = p_int_state,  & ! interpolation info
                              p_grad = hgrad,         & ! horizontal gradients
                              lacc=.FALSE.,           & ! OpenAcc flag
                              opt_slev = 2,           & ! lowest vertical level
                              opt_elev = nlev,        & ! highest vertical level
                              opt_rlstart = rl_start, &
                              opt_rlend   = rl_end    )
    g2_hgrad = REAL(hgrad, KIND=wp)

    hgrad(:,:,:,:) = 0.
    CALL grad_green_gauss_cell(p_cc = dudz,            & ! cell centered variable to interpolate
                              ptr_patch = p_patch,    & ! patch info
                              ptr_int = p_int_state,  & ! interpolation info
                              p_grad = hgrad,         & ! horizontal gradients
                              lacc=.FALSE.,           & ! OpenAcc flag
                              opt_slev = 2,           & ! lowest vertical level
                              opt_elev = nlev,        & ! highest vertical level
                              opt_rlstart = rl_start, &
                              opt_rlend   = rl_end    )
    dudz_hgrad = REAL(hgrad, KIND=wp)

    hgrad(:,:,:,:) = 0.
    CALL grad_green_gauss_cell(p_cc = dvdz,            & ! cell centered variable to interpolate
                              ptr_patch = p_patch,    & ! patch info
                              ptr_int = p_int_state,  & ! interpolation info
                              p_grad = hgrad,         & ! horizontal gradients
                              lacc=.FALSE.,           & ! OpenAcc flag
                              opt_slev = 2,           & ! lowest vertical level
                              opt_elev = nlev,        & ! highest vertical level
                              opt_rlstart = rl_start, &
                              opt_rlend   = rl_end    )
    dvdz_hgrad = REAL(hgrad, KIND=wp)

    hgrad(:,:,:,:) = 0.
    CALL grad_green_gauss_cell(p_cc = dn2dz,           & ! cell centered variable to interpolate
                              ptr_patch = p_patch,    & ! patch info
                              ptr_int = p_int_state,  & ! interpolation info
                              p_grad = hgrad,         & ! horizontal gradients
                              lacc=.FALSE.,           & ! OpenAcc flag
                              opt_slev = 2,           & ! lowest vertical level
                              opt_elev = nlev,        & ! highest vertical level
                              opt_rlstart = rl_start, &
                              opt_rlend   = rl_end    )
    dn2dz_hgrad = REAL(hgrad, KIND=wp)

    hgrad(:,:,:,:) = 0.
    CALL grad_green_gauss_cell(p_cc = dg2dz,           & ! cell centered variable to interpolate
                              ptr_patch = p_patch,    & ! patch info
                              ptr_int = p_int_state,  & ! interpolation info
                              p_grad = hgrad,         & ! horizontal gradients
                              lacc=.FALSE.,           & ! OpenAcc flag
                              opt_slev = 2,           & ! lowest vertical level
                              opt_elev = nlev,        & ! highest vertical level
                              opt_rlstart = rl_start, &
                              opt_rlend   = rl_end    )
    dg2dz_hgrad = REAL(hgrad, KIND=wp)

    IF (.NOT. is_plane_torus) THEN   ! convert to dA/dlon/cos(lat) or dA/dlat
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx
            u_hgrad (:,jc,jk,jb) = u_hgrad (:,jc,jk,jb)*grid_sphere_radius
            v_hgrad (:,jc,jk,jb) = v_hgrad (:,jc,jk,jb)*grid_sphere_radius
            n2_hgrad(:,jc,jk,jb) = n2_hgrad(:,jc,jk,jb)*grid_sphere_radius
            g2_hgrad(:,jc,jk,jb) = g2_hgrad(:,jc,jk,jb)*grid_sphere_radius
            dudz_hgrad (:,jc,jk,jb) = dudz_hgrad (:,jc,jk,jb)*grid_sphere_radius
            dvdz_hgrad (:,jc,jk,jb) = dvdz_hgrad (:,jc,jk,jb)*grid_sphere_radius
            dn2dz_hgrad(:,jc,jk,jb) = dn2dz_hgrad(:,jc,jk,jb)*grid_sphere_radius
            dg2dz_hgrad(:,jc,jk,jb) = dg2dz_hgrad(:,jc,jk,jb)*grid_sphere_radius
          ENDDO
        ENDDO
      ENDDO
    END IF

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        ! dA/dz = 0 to prevent extrapolation
        u_hgrad (:,jc,1     ,jb) = u_hgrad (:,jc,2   ,jb)
        u_hgrad (:,jc,nlevp1,jb) = u_hgrad (:,jc,nlev,jb)
        v_hgrad (:,jc,1     ,jb) = v_hgrad (:,jc,2   ,jb)
        v_hgrad (:,jc,nlevp1,jb) = v_hgrad (:,jc,nlev,jb)
        n2_hgrad(:,jc,1     ,jb) = n2_hgrad(:,jc,2   ,jb)
        n2_hgrad(:,jc,nlevp1,jb) = n2_hgrad(:,jc,nlev,jb)
        g2_hgrad(:,jc,1     ,jb) = g2_hgrad(:,jc,2   ,jb)
        g2_hgrad(:,jc,nlevp1,jb) = g2_hgrad(:,jc,nlev,jb)
        dudz_hgrad (:,jc,1     ,jb) = 0.
        dudz_hgrad (:,jc,nlevp1,jb) = 0.
        dvdz_hgrad (:,jc,1     ,jb) = 0.
        dvdz_hgrad (:,jc,nlevp1,jb) = 0.
        dn2dz_hgrad(:,jc,1     ,jb) = 0.
        dn2dz_hgrad(:,jc,nlevp1,jb) = 0.
        dg2dz_hgrad(:,jc,1     ,jb) = 0.
        dg2dz_hgrad(:,jc,nlevp1,jb) = 0.
      ENDDO
    ENDDO
#endif

  END SUBROUTINE fieldsgrads
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE sync_wave(p_ray,p_patch)
    TYPE(t_ray),            INTENT(INOUT) :: p_ray                ! the ray properties
    TYPE(t_patch),   TARGET,INTENT(IN) :: p_patch              ! grid/patch info.

    IF (msg_level >= 12) CALL message('sync_wave', 'MS-GWaM: synchronize ray volumes')

    !----------------------------------------------------------------------
    ! Purpose:
    !         Synchronize derived type elements of ray volume properties
    !
    ! Method:
    !         -- based on sync_patch_array
    !            (src/parallel_infrastructure/mo_sync.f90)
    !
    !----------------------------------------------------------------------

    !CALL message('', TRIM(message_text))
    !      WRITE(message_text,'(a,i6)') 'shape p_ray%iexist:', SHAPE(p_ray%iexist)
    !print*, 'SHAPE p_ray%iexist:', SHAPE(p_ray%iexist)

    ! Synchronize integer properties
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%iexist, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%specid, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%jk_active, lacc = .FALSE.)
    ! For the 6 properties below no synchronization is needed
    ! because they matter either only in the source treatment or
    ! they are anyway recalculated in the new position by idx_rayedge
    ! p_ray_bg%jk_source
    ! p_ray_conv%jk_source
    ! p_ray%jr_last
    ! p_ray%jk_full_rtop
    ! p_ray%jk_full_rbot
    ! p_ray%jk_half_rtop
    ! p_ray%jk_half_rbot

    ! Synchronize real properties
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%lon, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%lat, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%z, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%dx, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%dy, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%dz, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%k, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%l, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%m, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%dk, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%dl, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%dm, lacc = .FALSE.)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%wadens, lacc = .FALSE.)

    ! Later, directly calculate the cosine here rather than perform MPI-comm. (if it is faster)
    CALL sync_patch_array(SYNC_C, p_patch, p_ray%coslat, lacc = .FALSE.)

  END SUBROUTINE sync_wave
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE regrid_wave( p_patch,        & !in(out)
    &                     p_gridinfo4ray, & !in
    &                     z_ifc,          & !in
    &                     z_mc,           & !in
    &                     gammash2,       & !in
    &                     fc2,            & !in
    &                     bvf2,           & !in
    &                     jray_start,     & !in
    &                     jray_end,       & !in
    &                     p_ray           ) !inout

    ! In/out variables
    TYPE(t_patch),        INTENT(INOUT) :: p_patch              ! grid/patch info.
    TYPE(t_gridinfo4ray), INTENT(IN)    :: p_gridinfo4ray       ! infos on grid cell and neighbors
    REAL(wp),             INTENT(IN)    :: z_ifc(:,:,:)         ! cell interface heights
    REAL(wp),             INTENT(IN)    :: z_mc (:,:,:)         ! cell center heights
    REAL(wp),             INTENT(IN)    :: gammash2(:,:,:)      ! pink scale height squared
    REAL(wp),             INTENT(IN)    :: fc2(:,:)             ! Coriolis parameter squared
    REAL(wp),             INTENT(IN)    :: bvf2(:,:,:)          ! buoyancy frequency squared
    INTEGER,              INTENT(IN   ) :: jray_start           ! jray lower limit
    INTEGER,              INTENT(IN   ) :: jray_end             ! jray upper limit
    TYPE(t_ray),          INTENT(INOUT) :: p_ray                ! ray properties

    ! Local variables
    TYPE(t_point)                       :: v1,v2,v3                           ! vertex longitudes/latitudes

    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: jb, jc, jstencil
    INTEGER  :: jray_loc, jray_nb
    INTEGER  :: jray_work, jray_launch, jray_insert, nrays_vacate
    INTEGER  :: specid
    INTEGER  :: icn, ibn
    INTEGER  :: jk, jk0, jk0m1
    INTEGER  :: num_rays_before  ! ray counters (per layer, before regridding)

    REAL(wp) :: ray_xyz(3)

    INTEGER  :: intriangle                        ! checking if a ray volume is in the current triangle
    INTEGER  :: nrays_per_layer(p_patch%nlevp1)   ! counter of ray volumes per layer in an individual column

    LOGICAL  :: l_warning(nproma)                  ! warn if all rays are to be removed in a column
    INTEGER  ::  jg
    jg = p_patch%id
    !----------------------------------------------------------------------
    ! Purpose:
    !         -- moves all propagated ray volumes to the next cell's ray
    !            data structures
    !         -- rays which are not launched yet, are not moved
    !         -- when propagation leads to overflowing arrays, as many rays
    !            as needed are removed
    !         -- alternatively, a merge routine could be implemented here
    !
    ! Method:
    !         -- move all but launching rays to a work array
    !         -- check if there are propagated rays in neighboring cells
    !         -- move all transported rays to the local cell
    !         -- check if the number of rays is now too large
    !         -- if yes: remove just as many rays as needed from work array
    !         -- finally fill the local ray volume array with the work
    !            array and take care of the slote filled with launching rays
    !
    ! TODO:
    !         -- speed up regridding for cells which need no regridding
    !            (currently all cells undergo the copying procedures)
    !
    !----------------------------------------------------------------------

    ! Prognostic domain
    rl_start   = grf_bdywidth_c + 1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    IF (.NOT. is_plane_torus) THEN

      ! Find ray volumes propagated to this cell from neighboring
      ! cells and initialize them in the local cell
      DO jb = i_startblk, i_endblk ! jb index of the local cell

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        l_warning(:) = .FALSE.

        DO jc = i_startidx, i_endidx ! jc index of the local cell

          ! Vertex 1 in Cartesian coordinates
          v1%x = p_gridinfo4ray%cellvertices_x(1,jc,jb)
          v1%y = p_gridinfo4ray%cellvertices_y(1,jc,jb)
          v1%z = p_gridinfo4ray%cellvertices_z(1,jc,jb)

          ! Vertex 2 in Cartesian coordinates
          v2%x = p_gridinfo4ray%cellvertices_x(2,jc,jb)
          v2%y = p_gridinfo4ray%cellvertices_y(2,jc,jb)
          v2%z = p_gridinfo4ray%cellvertices_z(2,jc,jb)

          ! Vertex 3 in Cartesian coordinates
          v3%x = p_gridinfo4ray%cellvertices_x(3,jc,jb)
          v3%y = p_gridinfo4ray%cellvertices_y(3,jc,jb)
          v3%z = p_gridinfo4ray%cellvertices_z(3,jc,jb)

          ! clean auxiliary arrays
          nrays_per_layer(:) = 0
          ! initialization is omitted due to performance problems
          ! care must be taken when handliong ray_coll to avoid
          ! previously added or uninitialized values
          ! ray_coll(:, :)%iexist = 0
          ! ray_coll(:, :)%wadens = 0

          ! set counter of rays in launch process to zero
          jray_launch = 0

          ! fill auxiliary array with local ray volumes
          DO jray_loc=jray_start, jray_end

            ! skip if there is no ray in current slot or spoectral element cannot be identified
            IF (p_ray%iexist(jc, jray_loc, jb) == 0) CYCLE

            ! skip if ray is not launched yet, or is in jr_last
            ! note: rays in jr_last(jc, jspec, jb) are not moved in the array but can be regridded below
            ! increment counter to estimate nrays_vacate (number of rays for deletion)
            specid = ABS(p_ray%specid(jc, jray_loc, jb))
            IF (specid /= 0) THEN
              IF (p_ray%specid(jc, jray_loc, jb) < 0 .OR. p_ray%jr_last(jc, specid, jb) == jray_loc) THEN
                jray_launch = jray_launch + 1
                CYCLE
              ENDIF
            ENDIF

            ! get vertical index
            jk = p_ray%iexist(jc, jray_loc, jb)

            ! increment the number of rays per layer
            nrays_per_layer(jk) = nrays_per_layer(jk) + 1

            ! copy ray to aucilliary array
            ray_coll(nrays_per_layer(jk), jk)%iexist    = p_ray%iexist   (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%specid    = p_ray%specid   (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%jk_active = p_ray%jk_active(jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%lon       = p_ray%lon      (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%lat       = p_ray%lat      (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%z         = p_ray%z        (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%dx        = p_ray%dx       (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%dy        = p_ray%dy       (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%dz        = p_ray%dz       (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%coslat    = p_ray%coslat   (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%k         = p_ray%k        (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%l         = p_ray%l        (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%m         = p_ray%m        (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%dk        = p_ray%dk       (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%dl        = p_ray%dl       (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%dm        = p_ray%dm       (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%wadens    = p_ray%wadens   (jc, jray_loc, jb)

            ! delete ray from p_ray
            p_ray%iexist(jc, jray_loc, jb) = 0
          END DO

          ! store number of rays in cell before regridding
          num_rays_before = INT(SUM(nrays_per_layer))

          ! Loop over neighboring cell indices
          DO jstencil = 1, p_gridinfo4ray%cellneighbors_nstencil(jc,jb)  ! 11 or 12

            icn = p_gridinfo4ray%cellneighbors_idx(jstencil,jc,jb) ! jc index of neighboring cell
            ibn = p_gridinfo4ray%cellneighbors_blk(jstencil,jc,jb) ! jb index of neighboring cell

            ! Loop over ray volumes in neighboring cell
            DO jray_nb = jray_start, jray_end

              ! Skip calculations if this ray volume in the neighboring cell
              ! does not exist or is below the launch level z_src
              IF (p_ray%iexist(icn,jray_nb,ibn) == 0 .OR. &
                  p_ray%specid(icn,jray_nb,ibn) < 0) CYCLE

              ! Ray volume center point in Cartesian coordinates
              ray_xyz(1) = p_ray%coslat(icn,jray_nb,ibn)*COS(p_ray%lon(icn,jray_nb,ibn))
              ray_xyz(2) = p_ray%coslat(icn,jray_nb,ibn)*SIN(p_ray%lon(icn,jray_nb,ibn))
              ray_xyz(3) = SIGN(SQRT(1. - p_ray%coslat(icn,jray_nb,ibn)**2), p_ray%lat(icn,jray_nb,ibn))

              intriangle = inside_triangle(ray_xyz, v1,v2,v3)

              IF (intriangle >= 0) THEN

                jk0 = p_ray%iexist(icn,jray_nb,ibn)
                jk0m1 = MAX(1, jk0 - 1)

                IF ( z_mc(icn,jk0,ibn) /= z_mc(jc,jk0,jb) .OR. z_mc(icn,jk0m1,ibn) /= z_mc(jc,jk0m1,jb) ) THEN
                  !
                  ! Some conditions used at the end of the integration (in propagate_wave) need
                  ! to be applied again here, as the vertical grids are changed.
                  !  - lower bound used in propagate_wave
                  IF (p_ray%z(icn,jray_nb,ibn) < z_ifc(jc, p_patch% nlevp1 - nlev_lbnd(jg), jb)) THEN
                    ! outofdomain
                    p_ray%iexist(icn,jray_nb,ibn) = 0      ! TODO: check whether it is problematic for OpenMP ?
                    CYCLE
                  END IF
                  !
                  ! update jk0
                  CALL get_jkhalf_closest_gen( jk0, p_ray% z(icn,jray_nb,ibn),  &
                    &                          p_patch% nlev, z_mc(:,:,jb), jc )
                END IF

                ! Initialize new ray volume in work array
                ! Synchronized properties: take them from the ray volume that moved
                ! here from neighboring cell
                nrays_per_layer(jk0) = nrays_per_layer(jk0) + 1
                ray_coll(nrays_per_layer(jk0), jk0)%iexist    = p_ray%iexist   (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%specid    = p_ray%specid   (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%jk_active = p_ray%jk_active(icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%lon       = p_ray%lon      (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%lat       = p_ray%lat      (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%z         = p_ray%z        (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%dx        = p_ray%dx       (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%dy        = p_ray%dy       (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%dz        = p_ray%dz       (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%coslat    = p_ray%coslat   (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%k         = p_ray%k        (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%l         = p_ray%l        (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%m         = p_ray%m        (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%dk        = p_ray%dk       (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%dl        = p_ray%dl       (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%dm        = p_ray%dm       (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%wadens    = p_ray%wadens   (icn, jray_nb, ibn)

                ! Remove from the neighboring cell
                p_ray%iexist(icn, jray_nb, ibn) = 0      ! TODO: check whether it is problematic for OpenMP ?

                ! if index of propagating ray is in jr_last, delete it from there
                ! - rays propagated away from the source are no longer referenced as last launched ray
                ! - this should strictly speaking not happen as specid of ghost rays should be negative
                !   and are excluded above
                specid = ABS(p_ray%specid(icn, jray_nb, ibn))
                IF (specid /=0) THEN    ! ony check for last launched rays if spectral ID can be identified
                  IF (p_ray%jr_last(icn, specid, ibn) == jray_nb) THEN
                    CALL message('', 'MS-GWaM WARNING: regridding ray associated to source (is in jr_last), dissociating.')
                    p_ray%jr_last(icn, specid, ibn) = 0
                  ENDIF
                ENDIF

                ! Horizontal propagation test
                IF (ltest_hprop .AND. l1ray) THEN
                  ! Open output file for horizontal propagation test
                  OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
                  WRITE(987,'(a,3i4,2x,6F14.5)') &
                        'Regrid ray taken from lon, lat:', jstencil, icn, ibn, &
                                            p_gridinfo4ray%cellvertices_lon(1,icn,ibn)*rad2deg,&
                                            p_gridinfo4ray%cellvertices_lat(1,icn,ibn)*rad2deg,&
                                            p_gridinfo4ray%cellvertices_lon(2,icn,ibn)*rad2deg,&
                                            p_gridinfo4ray%cellvertices_lat(2,icn,ibn)*rad2deg,&
                                            p_gridinfo4ray%cellvertices_lon(3,icn,ibn)*rad2deg,&
                                            p_gridinfo4ray%cellvertices_lat(3,icn,ibn)*rad2deg
                  WRITE(987,'(a,2i4,2x,6F14.5)') &
                        'Regrid ray taken   to lon, lat:', jc, jb, &
                                            p_gridinfo4ray%cellvertices_lon(1,jc,jb)*rad2deg,&
                                            p_gridinfo4ray%cellvertices_lat(1,jc,jb)*rad2deg,&
                                            p_gridinfo4ray%cellvertices_lon(2,jc,jb)*rad2deg,&
                                            p_gridinfo4ray%cellvertices_lat(2,jc,jb)*rad2deg,&
                                            p_gridinfo4ray%cellvertices_lon(3,jc,jb)*rad2deg,&
                                            p_gridinfo4ray%cellvertices_lat(3,jc,jb)*rad2deg
                  CLOSE(987)
                ENDIF ! ltest_hprop .AND. l1ray

              ENDIF ! intriangle

            ENDDO  ! jray_nb

          ENDDO  ! jstencil

          ! now `ray_coll` contains all ray volums associated to the current column
          ! the dimensions of `pray_col` are: (rays, #vertical_layers).
          ! `nrays_per_layer(#vertical_layers)` describes the number of rays per vert. layer
          ! this construction should allow for running a merging algorithm per layer on the current column
          !TODO: reintroduce a merging algorithm
          IF (imethod_merge > 0) THEN
            CALL finish('regrid_wave', 'MS-GWaM: Merging ray volumes is currently disabled and needs reimplementation. Aborting.')
          ELSE
            ! rather than merging, we now move all rays into the first (highest) vertical layer of the auxiliary array
            ! ray_coll has `nrays_coll = 39 * nrays` ray slots -> it should always fit all rays
            jray_work = nrays_per_layer(1)
            DO jk = 2, p_patch%nlevp1
              DO jray_loc = 1, nrays_per_layer(jk)

                jray_work = jray_work + 1

                ray_coll(jray_work, 1)%iexist    = ray_coll(jray_loc, jk)%iexist
                ray_coll(jray_work, 1)%specid    = ray_coll(jray_loc, jk)%specid
                ray_coll(jray_work, 1)%jk_active = ray_coll(jray_loc, jk)%jk_active
                ray_coll(jray_work, 1)%lon       = ray_coll(jray_loc, jk)%lon
                ray_coll(jray_work, 1)%lat       = ray_coll(jray_loc, jk)%lat
                ray_coll(jray_work, 1)%z         = ray_coll(jray_loc, jk)%z
                ray_coll(jray_work, 1)%dx        = ray_coll(jray_loc, jk)%dx
                ray_coll(jray_work, 1)%dy        = ray_coll(jray_loc, jk)%dy
                ray_coll(jray_work, 1)%dz        = ray_coll(jray_loc, jk)%dz
                ray_coll(jray_work, 1)%coslat    = ray_coll(jray_loc, jk)%coslat
                ray_coll(jray_work, 1)%k         = ray_coll(jray_loc, jk)%k
                ray_coll(jray_work, 1)%l         = ray_coll(jray_loc, jk)%l
                ray_coll(jray_work, 1)%m         = ray_coll(jray_loc, jk)%m
                ray_coll(jray_work, 1)%dk        = ray_coll(jray_loc, jk)%dk
                ray_coll(jray_work, 1)%dl        = ray_coll(jray_loc, jk)%dl
                ray_coll(jray_work, 1)%dm        = ray_coll(jray_loc, jk)%dm
                ray_coll(jray_work, 1)%wadens    = ray_coll(jray_loc, jk)%wadens

              ENDDO
            ENDDO

            ! deactivate all other rays after moving to array of first layer
            ray_coll(1:MAXVAL(nrays_per_layer), 2:p_patch%nlevp1)%iexist = 0

            ! remove ray volumes only if needed
            ! if nrays_vacate > 0 there are more rays in the aux array
            ! than should be in the ray array range -> remove the surplus rays
            nrays_vacate = jray_work + jray_launch - (jray_end - jray_start + 1)

            if (nrays_vacate > 0) THEN


              ! perform removal on auxiliary array
              CALL remove_rays_in_cell(jray_start=1,                            & ! start of ray loop
                                      jray_end=jray_work,                      & !   end of ray loop
                                      nlev=p_patch%nlev,                       & ! number of vertical levels
                                      nrays_vacate=nrays_vacate,               & ! number of ray slots to be empty after removing
                                      specid=ray_coll(1:jray_work, 1)%specid,  & ! ray spectral ID
                                      iexist=ray_coll(1:jray_work, 1)%iexist,  & ! ray vertical layer / existence
                                      dens=ray_coll(1:jray_work, 1)%wadens,    & ! ray phase space wave action density
                                      kray=ray_coll(1:jray_work, 1)%k,         & ! ray wave numbers
                                      lray=ray_coll(1:jray_work, 1)%l,         & ! ray wave numbers
                                      mray=ray_coll(1:jray_work, 1)%m,         & ! ray wave numbers
                                      dkray=ray_coll(1:jray_work, 1)%dk,       & ! ray size in wave number direction
                                      dlray=ray_coll(1:jray_work, 1)%dl,       & ! ray size in wave number direction
                                      dmray=ray_coll(1:jray_work, 1)%dm,       & ! ray size in wave number direction
                                      dxray=ray_coll(1:jray_work, 1)%dx,       & ! ray size in physical directions
                                      dyray=ray_coll(1:jray_work, 1)%dy,       & ! ray size in physical directions
                                      dzray=ray_coll(1:jray_work, 1)%dz,       & ! ray size in physical directions
                                      gammash2=gammash2(jc,:,jb),              & ! pink scale height squared
                                      fc2=fc2(jc,jb),                          & ! Coriolis parameter squared
                                      bvf2=bvf2(jc,:,jb),                      & ! buoyancy frequency squared
                                      l_warning=l_warning(jc)                  ) ! warning if attempting to remove everything
            ENDIF

            ! copy all alive rays back to `p_ray`
            jray_insert = -1  ! starts at minus one so the first insertion is done at jray_start
            DO jray_loc = 1, jray_work

              ! dont account for removed rays
              IF (ray_coll(jray_loc, 1)%iexist == 0) THEN
                CYCLE
              ENDIF

              ! increment ray counter
              jray_insert = jray_insert + 1

              ! keep incrementing if there is an existing ray in the location of the array
              ! these rays must be part of the launch process (specid<0 / in jr_last)
              DO WHILE (p_ray%iexist(jc, jray_start+jray_insert, jb) /= 0)
                jray_insert = jray_insert + 1
              ENDDO

              ! replace rays within p_ray region with `ray_coll` members
              p_ray%iexist   (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%iexist
              p_ray%specid   (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%specid
              p_ray%jk_active(jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%jk_active
              p_ray%lon      (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%lon
              p_ray%lat      (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%lat
              p_ray%z        (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%z
              p_ray%dx       (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%dx
              p_ray%dy       (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%dy
              p_ray%dz       (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%dz
              p_ray%coslat   (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%coslat
              p_ray%k        (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%k
              p_ray%l        (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%l
              p_ray%m        (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%m
              p_ray%dk       (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%dk
              p_ray%dl       (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%dl
              p_ray%dm       (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%dm
              p_ray%wadens   (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%wadens

            ENDDO  ! jray_loc

          ENDIF  ! imethod_merge <= 0  -> merging is not activated

        ENDDO  ! jc

        ! add warning in case the removal deleted all local rays somewhere
        IF ( ANY(l_warning(:)) )  CALL message('',  &
          &  'regrid_wave: trying to remove too many local ray volumes')

      ENDDO  ! jb

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jray_loc, i_startidx, i_endidx, jstencil, icn, ibn, &
!$OMP v1, v2, v3, ray_xyz, intriangle, specid ) ICON_OMP_GUIDED_SCHEDULE

      ! Find ray volumes that propagated out of this
      ! cell and remove them
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          ! Vertex 1 in Cartesian coordinates
          v1%x = p_gridinfo4ray%cellvertices_x(1,jc,jb)
          v1%y = p_gridinfo4ray%cellvertices_y(1,jc,jb)
          v1%z = p_gridinfo4ray%cellvertices_z(1,jc,jb)

          ! Vertex 2 in Cartesian coordinates
          v2%x = p_gridinfo4ray%cellvertices_x(2,jc,jb)
          v2%y = p_gridinfo4ray%cellvertices_y(2,jc,jb)
          v2%z = p_gridinfo4ray%cellvertices_z(2,jc,jb)

          ! Vertex 3 in Cartesian coordinates
          v3%x = p_gridinfo4ray%cellvertices_x(3,jc,jb)
          v3%y = p_gridinfo4ray%cellvertices_y(3,jc,jb)
          v3%z = p_gridinfo4ray%cellvertices_z(3,jc,jb)

          DO jray_loc = jray_start, jray_end

            IF (p_ray%iexist(jc,jray_loc,jb) == 0 .OR. &
                p_ray%specid(jc,jray_loc,jb) < 0) CYCLE

            ! Remainder: ray volume either being inside this cell or escaped to a halo cell

            ! Find ray volumes that propagated out of this cell and
            ! remove them
            ! Ray volume center point in Cartesian coordinates
            ray_xyz(1) = p_ray%coslat(jc,jray_loc,jb)*COS(p_ray%lon(jc,jray_loc,jb))
            ray_xyz(2) = p_ray%coslat(jc,jray_loc,jb)*SIN(p_ray%lon(jc,jray_loc,jb))
            ray_xyz(3) = SIGN(SQRT(1. - p_ray%coslat(jc,jray_loc,jb)**2), p_ray%lat(jc,jray_loc,jb))

            intriangle = inside_triangle(ray_xyz, v1,v2,v3)

            IF (intriangle < 0) THEN  ! not in triangle

              ! Remove ray volume
              p_ray%iexist(jc,jray_loc,jb) = 0

              ! remove from last launched ray volumes if applicable
              specid = ABS(p_ray%specid(jc, jray_loc, jb))
              IF (specid /= 0) THEN
                IF (p_ray%jr_last(jc, specid, jb) == jray_loc) THEN
                  p_ray%jr_last(jc, specid, jb) = 0
                ENDIF
              ENDIF

              ! Horizontal propagation test
              IF (ltest_hprop .AND. l1ray) THEN
                ! Open output file for horizontal propagation test
                OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
                WRITE(987,'(a,2i4,2x,6F14.5)') &
                      'Regrid ray removed at lon, lat:', jc, jb, &
                                          p_gridinfo4ray%cellvertices_lon(1,jc,jb)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lat(1,jc,jb)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lon(2,jc,jb)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lat(2,jc,jb)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lon(3,jc,jb)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lat(3,jc,jb)*rad2deg
                CLOSE(987)
              ENDIF ! ltest_hprop .AND l1ray

            ENDIF ! not intriangle

          ENDDO  ! jray_loc

        ENDDO  ! jc

      ENDDO  ! jb

!$OMP END DO
!$OMP END PARALLEL

    ELSE ! is_plane_torus

      ! Find ray volumes propagated to this cell from neighboring
      ! cells and initialize them in the local cell
      DO jb = i_startblk, i_endblk ! jb index of the local cell

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        l_warning(:) = .FALSE.

        DO jc = i_startidx, i_endidx ! jc index of the local cell

          ! Vertex 1 in Cartesian coordinates
          v1%x = p_gridinfo4ray%cellvertices_x(1,jc,jb)
          v1%y = p_gridinfo4ray%cellvertices_y(1,jc,jb)
          v1%z = p_gridinfo4ray%cellvertices_z(1,jc,jb)

          ! Vertex 2 in Cartesian coordinates
          v2%x = p_gridinfo4ray%cellvertices_x(2,jc,jb)
          v2%y = p_gridinfo4ray%cellvertices_y(2,jc,jb)
          v2%z = p_gridinfo4ray%cellvertices_z(2,jc,jb)

          ! Vertex 3 in Cartesian coordinates
          v3%x = p_gridinfo4ray%cellvertices_x(3,jc,jb)
          v3%y = p_gridinfo4ray%cellvertices_y(3,jc,jb)
          v3%z = p_gridinfo4ray%cellvertices_z(3,jc,jb)

          ! clean auxiliary arrays
          nrays_per_layer(:) = 0
          ! initialization is omitted due to performance problems
          ! care must be taken when handliong ray_coll to avoid
          ! previously added or uninitialized values
          ! ray_coll(:, :)%iexist = 0
          ! ray_coll(:, :)%wadens = 0

          ! set counter of rays in launch process to zero
          jray_launch = 0

          ! fill auxiliary array with local ray volumes
          DO jray_loc=jray_start, jray_end

            ! skip if there is no ray in current slot or spoectral element cannot be identified
            IF (p_ray%iexist(jc, jray_loc, jb) == 0) CYCLE

            ! skip if ray is not launched yet, or is in jr_last
            ! note: rays in jr_last(jc, jspec, jb) are not moved in the array but can be regridded below
            ! increment counter to estimate nrays_vacate (number of rays for deletion)
            specid = ABS(p_ray%specid(jc, jray_loc, jb))
            IF (specid /= 0) THEN
              IF (p_ray%specid(jc, jray_loc, jb) < 0 .OR. p_ray%jr_last(jc, specid, jb) == jray_loc) THEN
                jray_launch = jray_launch + 1
                CYCLE
              ENDIF
            ENDIF

            ! get vertical index
            jk = p_ray%iexist(jc, jray_loc, jb)

            ! increment the number of rays per layer
            nrays_per_layer(jk) = nrays_per_layer(jk) + 1

            ! copy ray to aucilliary array
            ray_coll(nrays_per_layer(jk), jk)%iexist    = p_ray%iexist   (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%specid    = p_ray%specid   (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%jk_active = p_ray%jk_active(jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%lon       = p_ray%lon      (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%lat       = p_ray%lat      (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%z         = p_ray%z        (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%dx        = p_ray%dx       (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%dy        = p_ray%dy       (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%dz        = p_ray%dz       (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%coslat    = p_ray%coslat   (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%k         = p_ray%k        (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%l         = p_ray%l        (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%m         = p_ray%m        (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%dk        = p_ray%dk       (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%dl        = p_ray%dl       (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%dm        = p_ray%dm       (jc, jray_loc, jb)
            ray_coll(nrays_per_layer(jk), jk)%wadens    = p_ray%wadens   (jc, jray_loc, jb)

            ! delete ray from p_ray
            p_ray%iexist(jc, jray_loc, jb) = 0
          END DO

          ! store number of rays in cell before regridding
          num_rays_before = INT(SUM(nrays_per_layer))

          ! Loop over neighboring cell indices
          DO jstencil = 1, p_gridinfo4ray%cellneighbors_nstencil(jc,jb)  ! 11 or 12

            icn = p_gridinfo4ray%cellneighbors_idx(jstencil,jc,jb) ! jc index of neighboring cell
            ibn = p_gridinfo4ray%cellneighbors_blk(jstencil,jc,jb) ! jb index of neighboring cell

            ! Loop over ray volumes in neighboring cell
            DO jray_nb = jray_start, jray_end

              ! Skip calculations if this ray volume in the neighboring cell
              ! does not exist or is below the launch level z_src
              IF (p_ray%iexist(icn,jray_nb,ibn) == 0 .OR. &
                  p_ray%specid(icn,jray_nb,ibn) < 0) CYCLE

              ! Ray volume center point in Cartesian coordinates
              ray_xyz(1) = p_ray%lon(icn,jray_nb,ibn)
              ray_xyz(2) = p_ray%lat(icn,jray_nb,ibn)
              ray_xyz(3) = p_ray%z  (icn,jray_nb,ibn) ! not used on the plane

              intriangle = inside_triangle_torus(ray_xyz, v1,v2,v3)

              IF (intriangle >= 0) THEN

                jk0 = p_ray%iexist(icn,jray_nb,ibn)
                jk0m1 = MAX(1, jk0 - 1)

                IF ( z_mc(icn,jk0,ibn) /= z_mc(jc,jk0,jb) .OR. z_mc(icn,jk0m1,ibn) /= z_mc(jc,jk0m1,jb) ) THEN
                  !
                  ! Some conditions used at the end of the integration (in propagate_wave) need
                  ! to be applied again here, as the vertical grids are changed.
                  !  - lower bound used in propagate_wave
                  IF (p_ray%z(icn,jray_nb,ibn) < z_ifc(jc, p_patch% nlevp1 - nlev_lbnd(jg), jb)) THEN
                    ! outofdomain
                    p_ray%iexist(icn,jray_nb,ibn) = 0      ! TODO: check whether it is problematic for OpenMP ?
                    CYCLE
                  END IF
                  !
                  ! update jk0
                  CALL get_jkhalf_closest_gen( jk0, p_ray% z(icn,jray_nb,ibn),  &
                    &                          p_patch% nlev, z_mc(:,:,jb), jc )
                END IF

                ! Initialize new ray volume in work array
                ! Synchronized properties: take them from the ray volume that moved
                ! here from neighboring cell
                nrays_per_layer(jk0) = nrays_per_layer(jk0) + 1
                ray_coll(nrays_per_layer(jk0), jk0)%iexist    = p_ray%iexist   (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%specid    = p_ray%specid   (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%jk_active = p_ray%jk_active(icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%lon       = p_ray%lon      (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%lat       = p_ray%lat      (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%z         = p_ray%z        (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%dx        = p_ray%dx       (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%dy        = p_ray%dy       (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%dz        = p_ray%dz       (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%coslat    = p_ray%coslat   (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%k         = p_ray%k        (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%l         = p_ray%l        (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%m         = p_ray%m        (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%dk        = p_ray%dk       (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%dl        = p_ray%dl       (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%dm        = p_ray%dm       (icn, jray_nb, ibn)
                ray_coll(nrays_per_layer(jk0), jk0)%wadens    = p_ray%wadens   (icn, jray_nb, ibn)

                ! Remove from the neighboring cell
                p_ray%iexist(icn, jray_nb, ibn) = 0      ! TODO: check whether it is problematic for OpenMP ?

                ! if index of propagating ray is in jr_last, delete it from there
                ! - rays propagated away from the source are no longer referenced as last launched ray
                ! - this should strictly speaking not happen as specid of ghost rays should be negative
                !   and are excluded above
                specid = ABS(p_ray%specid(icn, jray_nb, ibn))
                IF (specid /=0) THEN    ! ony check for last launched rays if spectral ID can be identified
                  IF (p_ray%jr_last(icn, specid, ibn) == jray_nb) THEN
                    CALL message('', 'MS-GWaM WARNING: regridding ray associated to source (is in jr_last), dissociating.')
                    p_ray%jr_last(icn, specid, ibn) = 0
                  ENDIF
                ENDIF

                ! Horizontal propagation test
                IF (ltest_hprop .AND. l1ray) THEN
                  ! Open output file for horizontal propagation test
                  OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
                  WRITE(987,'(a,3i4,2x,6F14.5)') &
                        'Regrid ray taken from x, y:', jstencil, icn, ibn, &
                                            p_gridinfo4ray%cellvertices_x(1,icn,ibn),&
                                            p_gridinfo4ray%cellvertices_y(1,icn,ibn),&
                                            p_gridinfo4ray%cellvertices_x(2,icn,ibn),&
                                            p_gridinfo4ray%cellvertices_y(2,icn,ibn),&
                                            p_gridinfo4ray%cellvertices_x(3,icn,ibn),&
                                            p_gridinfo4ray%cellvertices_y(3,icn,ibn)
                  WRITE(987,'(a,2i4,2x,6F14.5)') &
                        'Regrid ray taken   to x, y:', jc, jb, &
                                            p_gridinfo4ray%cellvertices_x(1,jc,jb),&
                                            p_gridinfo4ray%cellvertices_y(1,jc,jb),&
                                            p_gridinfo4ray%cellvertices_x(2,jc,jb),&
                                            p_gridinfo4ray%cellvertices_y(2,jc,jb),&
                                            p_gridinfo4ray%cellvertices_x(3,jc,jb),&
                                            p_gridinfo4ray%cellvertices_y(3,jc,jb)
                  CLOSE(987)
                ENDIF ! ltest_hprop .AND. l1ray

              ENDIF ! intriangle

            ENDDO  ! jray_nb

          ENDDO  ! jstencil

          ! now `ray_coll` contains all ray volums associated to the current column
          ! the dimensions of `pray_col` are: (rays, #vertical_layers).
          ! `nrays_per_layer(#vertical_layers)` describes the number of rays per vert. layer
          ! this construction should allow for running a merging algorithm per layer on the current column
          !TODO: reintroduce a merging algorithm
          IF (imethod_merge > 0) THEN
            CALL finish('regrid_wave', 'MS-GWaM: Merging ray volumes is currently disabled and needs reimplementation. Aborting.')
          ELSE
            ! rather than merging, we now move all rays into the first (highest) vertical layer of the auxiliary array
            ! ray_coll has `nrays_coll = 39 * nrays` ray slots -> it should always fit all rays
            jray_work = nrays_per_layer(1)
            DO jk = 2, p_patch%nlevp1
              DO jray_loc = 1, nrays_per_layer(jk)

                jray_work = jray_work + 1

                ray_coll(jray_work, 1)%iexist    = ray_coll(jray_loc, jk)%iexist
                ray_coll(jray_work, 1)%specid    = ray_coll(jray_loc, jk)%specid
                ray_coll(jray_work, 1)%jk_active = ray_coll(jray_loc, jk)%jk_active
                ray_coll(jray_work, 1)%lon       = ray_coll(jray_loc, jk)%lon
                ray_coll(jray_work, 1)%lat       = ray_coll(jray_loc, jk)%lat
                ray_coll(jray_work, 1)%z         = ray_coll(jray_loc, jk)%z
                ray_coll(jray_work, 1)%dx        = ray_coll(jray_loc, jk)%dx
                ray_coll(jray_work, 1)%dy        = ray_coll(jray_loc, jk)%dy
                ray_coll(jray_work, 1)%dz        = ray_coll(jray_loc, jk)%dz
                ray_coll(jray_work, 1)%coslat    = ray_coll(jray_loc, jk)%coslat
                ray_coll(jray_work, 1)%k         = ray_coll(jray_loc, jk)%k
                ray_coll(jray_work, 1)%l         = ray_coll(jray_loc, jk)%l
                ray_coll(jray_work, 1)%m         = ray_coll(jray_loc, jk)%m
                ray_coll(jray_work, 1)%dk        = ray_coll(jray_loc, jk)%dk
                ray_coll(jray_work, 1)%dl        = ray_coll(jray_loc, jk)%dl
                ray_coll(jray_work, 1)%dm        = ray_coll(jray_loc, jk)%dm
                ray_coll(jray_work, 1)%wadens    = ray_coll(jray_loc, jk)%wadens

              ENDDO
            ENDDO

            ! deactivate all other rays after moving to array of first layer
            ray_coll(1:MAXVAL(nrays_per_layer), 2:p_patch%nlevp1)%iexist = 0

            ! remove ray volumes only if needed
            ! if nrays_vacate > 0 there are more rays in the aux array
            ! than should be in the ray array range -> remove the surplus rays
            nrays_vacate = jray_work + jray_launch - (jray_end - jray_start + 1)

            if (nrays_vacate > 0) THEN


              ! perform removal on auxiliary array
              CALL remove_rays_in_cell(jray_start=1,                            & ! start of ray loop
                                      jray_end=jray_work,                      & !   end of ray loop
                                      nlev=p_patch%nlev,                       & ! number of vertical levels
                                      nrays_vacate=nrays_vacate,               & ! number of ray slots to be empty after removing
                                      specid=ray_coll(1:jray_work, 1)%specid,  & ! ray spectral ID
                                      iexist=ray_coll(1:jray_work, 1)%iexist,  & ! ray vertical layer / existence
                                      dens=ray_coll(1:jray_work, 1)%wadens,    & ! ray phase space wave action density
                                      kray=ray_coll(1:jray_work, 1)%k,         & ! ray wave numbers
                                      lray=ray_coll(1:jray_work, 1)%l,         & ! ray wave numbers
                                      mray=ray_coll(1:jray_work, 1)%m,         & ! ray wave numbers
                                      dkray=ray_coll(1:jray_work, 1)%dk,       & ! ray size in wave number direction
                                      dlray=ray_coll(1:jray_work, 1)%dl,       & ! ray size in wave number direction
                                      dmray=ray_coll(1:jray_work, 1)%dm,       & ! ray size in wave number direction
                                      dxray=ray_coll(1:jray_work, 1)%dx,       & ! ray size in physical directions
                                      dyray=ray_coll(1:jray_work, 1)%dy,       & ! ray size in physical directions
                                      dzray=ray_coll(1:jray_work, 1)%dz,       & ! ray size in physical directions
                                      gammash2=gammash2(jc,:,jb),              & ! pink scale height squared
                                      fc2=fc2(jc,jb),                          & ! Coriolis parameter squared
                                      bvf2=bvf2(jc,:,jb),                      & ! buoyancy frequency squared
                                      l_warning=l_warning(jc)                  ) ! warning if attempting to remove everything
            ENDIF

            ! copy all alive rays back to `p_ray`
            jray_insert = -1  ! starts at minus one so the first insertion is done at jray_start
            DO jray_loc = 1, jray_work

              ! dont account for removed rays
              IF (ray_coll(jray_loc, 1)%iexist == 0) THEN
                CYCLE
              ENDIF

              ! increment ray counter
              jray_insert = jray_insert + 1

              ! keep incrementing if there is an existing ray in the location of the array
              ! these rays must be part of the launch process (specid<0 / in jr_last)
              DO WHILE (p_ray%iexist(jc, jray_start+jray_insert, jb) /= 0)
                jray_insert = jray_insert + 1
              ENDDO

              ! replace rays within p_ray region with `ray_coll` members
              p_ray%iexist   (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%iexist
              p_ray%specid   (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%specid
              p_ray%jk_active(jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%jk_active
              p_ray%lon      (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%lon
              p_ray%lat      (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%lat
              p_ray%z        (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%z
              p_ray%dx       (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%dx
              p_ray%dy       (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%dy
              p_ray%dz       (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%dz
              p_ray%coslat   (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%coslat
              p_ray%k        (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%k
              p_ray%l        (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%l
              p_ray%m        (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%m
              p_ray%dk       (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%dk
              p_ray%dl       (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%dl
              p_ray%dm       (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%dm
              p_ray%wadens   (jc, jray_start+jray_insert, jb) = ray_coll(jray_loc, 1)%wadens

            ENDDO  ! jray_loc

          ENDIF  ! imethod_merge <= 0  -> merging is not activated

        ENDDO  ! jc

        ! add warning in case the removal deleted all local rays somewhere
        IF ( ANY(l_warning(:)) )  CALL message('',  &
          &  'regrid_wave: trying to remove too many local ray volumes')

      ENDDO  ! jb


!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jray_loc, i_startidx, i_endidx, jstencil, icn, ibn, &
!$OMP v1, v2, v3, ray_xyz, intriangle, specid ) ICON_OMP_GUIDED_SCHEDULE

      ! Find ray volumes that propagated out of this
      ! cell and remove them
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          ! Vertex 1 in Cartesian coordinates
          v1%x = p_gridinfo4ray%cellvertices_x(1,jc,jb)
          v1%y = p_gridinfo4ray%cellvertices_y(1,jc,jb)
          v1%z = p_gridinfo4ray%cellvertices_z(1,jc,jb)

          ! Vertex 2 in Cartesian coordinates
          v2%x = p_gridinfo4ray%cellvertices_x(2,jc,jb)
          v2%y = p_gridinfo4ray%cellvertices_y(2,jc,jb)
          v2%z = p_gridinfo4ray%cellvertices_z(2,jc,jb)

          ! Vertex 3 in Cartesian coordinates
          v3%x = p_gridinfo4ray%cellvertices_x(3,jc,jb)
          v3%y = p_gridinfo4ray%cellvertices_y(3,jc,jb)
          v3%z = p_gridinfo4ray%cellvertices_z(3,jc,jb)

          DO jray_loc = jray_start, jray_end

            IF (p_ray%iexist(jc,jray_loc,jb) == 0 .OR. &
                p_ray%specid(jc,jray_loc,jb) < 0) CYCLE

            ! Remainder: ray volume either being inside this cell or escaped to a halo cell

            ! Find ray volumes that propagated out of this cell and
            ! remove them
            ! Ray volume center point in Cartesian coordinates
            ray_xyz(1) = p_ray%lon(jc,jray_loc,jb)
            ray_xyz(2) = p_ray%lat(jc,jray_loc,jb)
            ray_xyz(3) = p_ray%  z(jc,jray_loc,jb) ! not used on the plane

            intriangle = inside_triangle_torus(ray_xyz, v1,v2,v3)

            IF (intriangle < 0) THEN  ! not in triangle

              ! Remove ray volume
              p_ray%iexist(jc,jray_loc,jb) = 0

              ! remove from last launched ray volumes if applicable
              specid = ABS(p_ray%specid(jc, jray_loc, jb))
              IF (specid /= 0) THEN
                IF (p_ray%jr_last(jc, specid, jb) == jray_loc) THEN
                  p_ray%jr_last(jc, specid, jb) = 0
                ENDIF
              ENDIF

              ! Horizontal propagation test
              IF (ltest_hprop .AND. l1ray) THEN
                ! Open output file for horizontal propagation test
                OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
                WRITE(987,'(a,2i4,2x,6F14.5)') &
                      'Regrid ray removed at x, y:', jc, jb, &
                                          p_gridinfo4ray%cellvertices_x(1,jc,jb),&
                                          p_gridinfo4ray%cellvertices_y(1,jc,jb),&
                                          p_gridinfo4ray%cellvertices_x(2,jc,jb),&
                                          p_gridinfo4ray%cellvertices_y(2,jc,jb),&
                                          p_gridinfo4ray%cellvertices_x(3,jc,jb),&
                                          p_gridinfo4ray%cellvertices_y(3,jc,jb)
                CLOSE(987)
              ENDIF ! ltest_hprop .AND l1ray

            ENDIF ! not intriangle

          ENDDO  ! jray_loc

        ENDDO  ! jc

      ENDDO  ! jb

!$OMP END DO
!$OMP END PARALLEL

    ENDIF ! is_plane_torus?

  END SUBROUTINE regrid_wave
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE tendency(p_patch, p_metrics, p_int_state, rho, temp, theta, p_fld)
    TYPE(t_patch),  TARGET,INTENT(IN)    :: p_patch
    TYPE(t_nh_metrics)    ,INTENT(IN)    :: p_metrics
    TYPE(t_int_state)     ,INTENT(IN)    :: p_int_state
    REAL(wp),      POINTER,INTENT(IN)    :: rho (:,:,:)   ! full lev
    REAL(wp),      POINTER,INTENT(IN)    :: temp (:,:,:)  ! full lev
    REAL(wp),      POINTER,INTENT(IN)    :: theta(:,:,:)  ! full lev
    TYPE(t_msgwam),        INTENT(INOUT) :: p_fld

    ! Local variables
    INTEGER               :: nlev
    INTEGER               :: rl_start, rl_end
    INTEGER               :: i_startblk, i_endblk
    INTEGER               :: i_startidx, i_endidx
    INTEGER               :: jk,jc,jce,jb,ist
    INTEGER               :: nblks_c,nblks_e
    INTEGER,      POINTER :: iidx(:,:,:), iblk(:,:,:)
    REAL(wp)              :: inv_rho, inv_dzrho
    REAL(wp)              :: fc_o_rho_th
    REAL(wp)              :: uufl_c1, uvfl_c1, vvfl_c1
    REAL(wp)              :: uufl_c2, uvfl_c2, vvfl_c2
    REAL(wp)              :: uupfl_c1, uvpfl_c1, vvpfl_c1
    REAL(wp)              :: uupfl_c2, uvpfl_c2, vvpfl_c2
    REAL(wp)              :: utfl_c1, vtfl_c1
    REAL(wp)              :: utfl_c2, vtfl_c2
    REAL(wp)              :: ddt_theta_gwd_mgm
    INTEGER               :: jk0
    REAL(wp)              :: inv_depth, flxgrad
    REAL(wp), ALLOCATABLE :: ufl_hor_e(:,:,:)  ! horizontal flux of u
                                              ! at cell edge midpoints
    REAL(wp), ALLOCATABLE :: vfl_hor_e(:,:,:)  ! horizontal flux of v
                                              ! at cell edge midpoints
    REAL(wp), ALLOCATABLE :: upfl_hor_e(:,:,:) ! horizontal pflux of u
                                              ! at cell edge midpoints
    REAL(wp), ALLOCATABLE :: vpfl_hor_e(:,:,:) ! horizontal pflux of v
                                              ! at cell edge midpoints
    REAL(wp), ALLOCATABLE :: ptfl_hor_e(:,:,:) ! horizontal flux of pot temp
                                              ! at cell edge midpoints
    REAL(wp), ALLOCATABLE :: ufldiv_hor(:,:,:) ! horizontal div of u fluxes
    REAL(wp), ALLOCATABLE :: vfldiv_hor(:,:,:) ! horizontal div of v fluxes
    REAL(wp), ALLOCATABLE :: upfldiv_hor(:,:,:) ! horizontal div of u fluxes
    REAL(wp), ALLOCATABLE :: vpfldiv_hor(:,:,:) ! horizontal div of v fluxes
    REAL(wp), ALLOCATABLE :: ptfldiv_hor(:,:,:)! horizontal div of pot temp fluxes

    REAL(wp), PARAMETER :: gam = cpd/cvd
    INTEGER :: jg
    jg = p_patch%id

    IF (msg_level >= 12) CALL message('tendency', 'MS-GWaM: wind tendency calculation')

    !----------------------------------------------------------------------
    ! Purpose:
    !         Calculate tendencies based on momentum flux convergence and
    !         the elastic terms
    !
    ! Method:
    !         tendencies = 1/rho*(vertical centered differences of pseudo-
    !                      momentum fluxes)
    !         staggering:
    !         -- pseudo-momentum fluxes defined on half levels
    !         -- tendencies defined on full levels
    !
    ! TODO:
    !
    !
    !----------------------------------------------------------------------

    nlev = p_patch%nlev

    IF ( msgw_lower_bound_opt == 1 ) THEN

      ! Exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c + 1
      rl_end   = min_rlcell_int

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

      jk0 = nlev+1 - nlev_lbnd(jg) - 1

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jk, i_startidx, i_endidx, inv_depth, flxgrad)

      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
          &                i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          inv_depth = 1._wp/(p_metrics%z_mc(jc,jk0,jb) - p_metrics%z_mc(jc,nlev,jb))
          !
#ifndef __msgwam1d
          flxgrad = p_fld%uufl_mgm(jc,jk0,jb)*inv_depth
          DO jk = jk0+1, nlev-1
            p_fld%uufl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
          ENDDO
          flxgrad = p_fld%uvfl_mgm(jc,jk0,jb)*inv_depth
          DO jk = jk0+1, nlev-1
            p_fld%uvfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
          ENDDO
          flxgrad = p_fld%vvfl_mgm(jc,jk0,jb)*inv_depth
          DO jk = jk0+1, nlev-1
            p_fld%vvfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
          ENDDO
          flxgrad = p_fld%utfl_mgm(jc,jk0,jb)*inv_depth
          DO jk = jk0+1, nlev-1
            p_fld%utfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
          ENDDO
          flxgrad = p_fld%vtfl_mgm(jc,jk0,jb)*inv_depth
          DO jk = jk0+1, nlev-1
            p_fld%vtfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
          ENDDO
          flxgrad = p_fld%uupfl_mgm(jc,jk0,jb) * inv_depth
          DO jk = jk0+1, nlev-1
            p_fld%uupfl_mgm(jc,jk,jb) = flxgrad * (p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
          ENDDO
          flxgrad = p_fld%uvpfl_mgm(jc,jk0,jb) * inv_depth
          DO jk = jk0+1, nlev-1
            p_fld%uvpfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
          ENDDO
          flxgrad = p_fld%vvpfl_mgm(jc,jk0,jb) * inv_depth
          DO jk = jk0+1, nlev-1
            p_fld%vvpfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
          ENDDO

          p_fld%uupfl_mgm(jc,nlev,jb) = 0._wp
          p_fld%uvpfl_mgm(jc,nlev,jb) = 0._wp
          p_fld%vvpfl_mgm(jc,nlev,jb) = 0._wp
          p_fld%uufl_mgm(jc,nlev,jb) = 0._wp
          p_fld%uvfl_mgm(jc,nlev,jb) = 0._wp
          p_fld%vvfl_mgm(jc,nlev,jb) = 0._wp
          p_fld%utfl_mgm(jc,nlev,jb) = 0._wp
          p_fld%vtfl_mgm(jc,nlev,jb) = 0._wp
          !
#endif
          flxgrad = 0.5_wp*(p_fld%uwfl_mgm(jc,jk0,jb) + p_fld%uwfl_mgm(jc,jk0+1,jb))*inv_depth
          DO jk = jk0+1, nlev
            p_fld%uwfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
          ENDDO
          flxgrad = 0.5_wp*(p_fld%vwfl_mgm(jc,jk0,jb) + p_fld%vwfl_mgm(jc,jk0+1,jb))*inv_depth
          DO jk = jk0+1, nlev
            p_fld%vwfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
          ENDDO
          flxgrad = 0.5_wp * (p_fld%uwpfl_mgm(jc,jk0,jb) + p_fld%uwpfl_mgm(jc,jk0+1,jb)) * inv_depth
          DO jk = jk0+1, nlev
            p_fld%uwpfl_mgm(jc,jk,jb) = flxgrad * (p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
          ENDDO
          flxgrad = 0.5_wp * (p_fld%vwpfl_mgm(jc,jk0,jb) + p_fld%vwpfl_mgm(jc,jk0+1,jb)) * inv_depth
          DO jk = jk0+1, nlev
            p_fld%vwpfl_mgm(jc,jk,jb) = flxgrad * (p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
          ENDDO
          p_fld%uwfl_mgm(jc,nlev+1,jb) = 0._wp
          p_fld%vwfl_mgm(jc,nlev+1,jb) = 0._wp
          p_fld%uwpfl_mgm(jc,nlev+1,jb) = 0._wp
          p_fld%vwpfl_mgm(jc,nlev+1,jb) = 0._wp
        ENDDO  ! jc
      ENDDO  ! jb

!$OMP END DO
!$OMP END PARALLEL

    END IF  ! msgw_lower_bound_opt == 1

#ifndef __msgwam1d
    !define pointers
    iidx  => p_patch%edges%cell_idx
    iblk  => p_patch%edges%cell_blk

    ! Local dimensions
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! Allocate auxiliary fields
    ALLOCATE( ufl_hor_e  (nproma,nlev,nblks_e), &
              vfl_hor_e  (nproma,nlev,nblks_e), &
              upfl_hor_e (nproma,nlev,nblks_e), &
              vpfl_hor_e (nproma,nlev,nblks_e), &
              ptfl_hor_e (nproma,nlev,nblks_e), &
              ufldiv_hor (nproma,nlev,nblks_c), &
              vfldiv_hor (nproma,nlev,nblks_c), &
              upfldiv_hor(nproma,nlev,nblks_c), &
              vpfldiv_hor(nproma,nlev,nblks_c), &
              ptfldiv_hor(nproma,nlev,nblks_c), &
              STAT=ist                        )
    IF (ist /= success) CALL finish ('tendency', 'Allocation of auxiliary fields failed!')

    ! Synchronize fields and fluxes
    CALL sync_patch_array_mult(SYNC_C1, p_patch, 5, lacc=.FALSE., f3din1=p_fld%uufl_mgm, f3din2=p_fld%uvfl_mgm,  &
      &                        f3din3=p_fld%vvfl_mgm, f3din4=p_fld%utfl_mgm, f3din5=p_fld%vtfl_mgm)
    CALL sync_patch_array_mult(SYNC_C1, p_patch, 3, lacc=.FALSE., f3din1=p_fld%uupfl_mgm, f3din2=p_fld%uvpfl_mgm, f3din3=p_fld%vvpfl_mgm)

    ! Initialize horizontal fluxes at cell edges
    ufl_hor_e = 0._wp ; vfl_hor_e = 0._wp ; ptfl_hor_e = 0._wp
    upfl_hor_e = 0._wp ; vpfl_hor_e = 0._wp

    ! Calculate normal horizontal fluxes at edge midpoints
    ! Exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_e + 1
    rl_end   = min_rledge_int

    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jce, i_startidx, i_endidx,                     &
!$OMP            uufl_c1, uvfl_c1, uufl_c2, uvfl_c2, vvfl_c1, vvfl_c2,  &
!$OMP            utfl_c1, vtfl_c1, utfl_c2, vtfl_c2)

    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
      DO jce = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO jce = i_startidx, i_endidx
#endif

          ! Horizontal flux of u at cell center points divided by rho
          uufl_c1 = p_fld%uufl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
          uvfl_c1 = p_fld%uvfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
          uufl_c2 = p_fld%uufl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
          uvfl_c2 = p_fld%uvfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
          ! Interpolation from cell center to cell edge midpoints
          ufl_hor_e(jce,jk,jb) = p_int_state%c_lin_e(jce,1,jb) &
                                * ( uufl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
                                +   uvfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
                              + p_int_state%c_lin_e(jce,2,jb) &
                                * ( uufl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
                                +   uvfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

          ! Horizontal flux of v at cell center points divided by rho
          vvfl_c1 = p_fld%vvfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
          vvfl_c2 = p_fld%vvfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
          ! Interpolation from cell center to cell edge midpoints
          vfl_hor_e(jce,jk,jb) = p_int_state%c_lin_e(jce,1,jb) &
                                * ( uvfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
                                +   vvfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
                              + p_int_state%c_lin_e(jce,2,jb) &
                                * ( uvfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
                                +   vvfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

          ! Horizontal pflux of u at cell center points divided by rho
          uupfl_c1 = p_fld%uupfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
          uvpfl_c1 = p_fld%uvpfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
          uupfl_c2 = p_fld%uupfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
          uvpfl_c2 = p_fld%uvpfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
          ! Interpolation from cell center to cell edge midpoints
          upfl_hor_e(jce,jk,jb) = p_int_state%c_lin_e(jce,1,jb) &
                                  * ( uupfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
                                  +   uvpfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
                                + p_int_state%c_lin_e(jce,2,jb) &
                                  * ( uupfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
                                  +   uvpfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

          ! Horizontal pflux of v at cell center points divided by rho
          vvpfl_c1 = p_fld%vvpfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
          vvpfl_c2 = p_fld%vvpfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
          ! Interpolation from cell center to cell edge midpoints
          vpfl_hor_e(jce,jk,jb) = p_int_state%c_lin_e(jce,1,jb) &
                                  * ( uvpfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
                                  +   vvpfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
                                + p_int_state%c_lin_e(jce,2,jb) &
                                  * ( uvpfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
                                  +   vvpfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

          ! Horizontal flux of pot. temperature at cell center points divided by rho
          utfl_c1 = p_fld%utfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
          vtfl_c1 = p_fld%vtfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
          utfl_c2 = p_fld%utfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
          vtfl_c2 = p_fld%vtfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
          ! Interpolation from cell center to cell edge midpoints
          ptfl_hor_e(jce,jk,jb) = p_int_state%c_lin_e(jce,1,jb) &
                                * ( utfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
                                +   vtfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
                                + p_int_state%c_lin_e(jce,2,jb) &
                                * ( utfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
                                +   vtfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

          ENDDO
        ENDDO
    ENDDO ! jb

!$OMP END DO
!$OMP END PARALLEL

    ! Synchronize horizontal fluxes at edge mid-points
    CALL sync_patch_array_mult(SYNC_E, p_patch, 3, lacc=.FALSE., f3din1=ufl_hor_e, f3din2=vfl_hor_e, f3din3=ptfl_hor_e)
    CALL sync_patch_array_mult(SYNC_E, p_patch, 2, lacc=.FALSE., f3din1=upfl_hor_e, f3din2=vpfl_hor_e)

    ! Calculate horizontal divergence of u flux
    CALL div_avg(ufl_hor_e,             & !in
                p_patch,               & !in
                p_int_state,           & !in
                p_int_state%c_bln_avg, & !in
                ufldiv_hor,            & !out
                lacc = .FALSE. )         !in

    ! Calculate horizontal divergence of v flux
    CALL div_avg(vfl_hor_e,             & !in
                p_patch,               & !in
                p_int_state,           & !in
                p_int_state%c_bln_avg, & !in
                vfldiv_hor,            & !out
                lacc = .FALSE. )         !in

    ! Calculate horizontal divergence of u pflux
    CALL div_avg(upfl_hor_e,            & !in
                p_patch,               & !in
                p_int_state,           & !in
                p_int_state%c_bln_avg, & !in
                upfldiv_hor,           & !out
                lacc = .FALSE. )         !in

    ! Calculate horizontal divergence of v pflux
    CALL div_avg(vpfl_hor_e,            & !in
                p_patch,               & !in
                p_int_state,           & !in
                p_int_state%c_bln_avg, & !in
                vpfldiv_hor,           & !out
                lacc = .FALSE. )         !in

    ! Calculate horizontal divergence of pot temperature flux
    CALL div_avg(ptfl_hor_e,            & !in
                p_patch,               & !in
                p_int_state,           & !in
                p_int_state%c_bln_avg, & !in
                ptfldiv_hor,           & !out
                lacc = .FALSE. )         !in
#endif

    ! Exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c + 1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx,  &
!$OMP            inv_rho, inv_dzrho, fc_o_rho_th, ddt_theta_gwd_mgm)

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          ! 1/rho, 1/(rho*dz)
          inv_rho   = 1._wp/rho(jc,jk,jb)
          inv_dzrho = 1._wp/(p_metrics%ddqz_z_full(jc,jk,jb)*rho(jc,jk,jb))

          ! f/(rho*theta)
          fc_o_rho_th = p_patch%cells%f_c(jc,jb)/(rho(jc,jk,jb)*theta(jc,jk,jb))

          ! Elastic terms (pot temperature fluxes are at full levels)
          p_fld%etx_mgm(jc,jk,jb) = - fc_o_rho_th * p_fld%vtfl_mgm(jc,jk,jb)
          p_fld%ety_mgm(jc,jk,jb) =   fc_o_rho_th * p_fld%utfl_mgm(jc,jk,jb)

          ! Calculate u tendency at cell centers on full levels
          p_fld%mfcxz_mgm(jc,jk,jb) = -(  p_fld%uwfl_mgm(jc,jk,jb)  &
                                        - p_fld%uwfl_mgm(jc,jk+1,jb) ) * inv_dzrho
#ifdef __msgwam1d
          p_fld%ddt_u_gwd_mgm(jc,jk,jb) = + p_fld%mfcxz_mgm(jc,jk,jb) & ! vert flux conv.
                + p_fld%etx_mgm(jc,jk,jb)                               ! elastic term
#else
          p_fld%ddt_u_gwd_mgm(jc,jk,jb) = + p_fld%mfcxz_mgm(jc,jk,jb) & ! vert flux conv.
                - ufldiv_hor(jc,jk,jb)*inv_rho                        & ! hori flux conv.
                + p_fld%etx_mgm(jc,jk,jb)                               ! elastic term
#endif

          ! Calculate v tendency at cell centers on full levels
          p_fld%mfcyz_mgm(jc,jk,jb) = -(  p_fld%vwfl_mgm(jc,jk,jb)  &
                                        - p_fld%vwfl_mgm(jc,jk+1,jb) ) * inv_dzrho

#ifdef __msgwam1d
          p_fld%ddt_v_gwd_mgm(jc,jk,jb) = p_fld%mfcyz_mgm(jc,jk,jb) & ! vert flux conv.
                + p_fld%ety_mgm(jc,jk,jb)                             ! elastic term
#else
          p_fld%ddt_v_gwd_mgm(jc,jk,jb) = p_fld%mfcyz_mgm(jc,jk,jb) & ! vert flux conv.
                - vfldiv_hor(jc,jk,jb)*inv_rho                      & ! hori flux conv.
                + p_fld%ety_mgm(jc,jk,jb)                             ! elastic term
#endif

          ! Calculate u tendency for pmom fluxes at cell centers on full levels
          p_fld%pmfcxz_mgm(jc,jk,jb) = -(  p_fld%uwpfl_mgm(jc,jk,jb)  &
                                        - p_fld%uwpfl_mgm(jc,jk+1,jb) ) * inv_dzrho

#ifdef __msgwam1d
          p_fld%ddt_u_gwd_pmom_mgm(jc,jk,jb) = p_fld%pmfcxz_mgm(jc,jk,jb)   ! vert flux conv.
#else
          p_fld%ddt_u_gwd_pmom_mgm(jc,jk,jb) = p_fld%pmfcxz_mgm(jc,jk,jb) & ! vert flux conv.
                - upfldiv_hor(jc,jk,jb)*inv_rho                             ! hori flux conv.
#endif

          ! Calculate v tendency at cell centers on full levels
          p_fld%pmfcyz_mgm(jc,jk,jb) = -(  p_fld%vwpfl_mgm(jc,jk,jb)  &
                                        - p_fld%vwpfl_mgm(jc,jk+1,jb) ) * inv_dzrho

#ifdef __msgwam1d
          p_fld%ddt_v_gwd_pmom_mgm(jc,jk,jb) = p_fld%pmfcyz_mgm(jc,jk,jb)   ! vert flux conv.
#else
          p_fld%ddt_v_gwd_pmom_mgm(jc,jk,jb) = p_fld%pmfcyz_mgm(jc,jk,jb) & ! vert flux conv.
                - vpfldiv_hor(jc,jk,jb)*inv_rho                             ! hori flux conv.
#endif

#ifndef __msgwam1d
          ! Calculate pot temp tendency at cell centers on full levels
          ddt_theta_gwd_mgm = - ptfldiv_hor(jc,jk,jb)*inv_rho   ! hori flux conv.
  !!!!! TEST !!!!! TO BE CLEANED LATER
          p_fld%ddt_pt_gwd_mgm(jc,jk,jb) = ddt_theta_gwd_mgm
  !!!!! TEST !!!!! TO BE CLEANED LATER
          ! Calculate temperature tendency from pot temp tendency
          p_fld%ddt_t_gwd_mgm(jc,jk,jb) = ddt_theta_gwd_mgm &
                      *temp(jc,jk,jb)/theta(jc,jk,jb)*gam  ! * pt_prog%exner(jc,jk,jb)*gam
#else
          p_fld%ddt_t_gwd_mgm(jc,jk,jb) = 0._wp
#endif

        ENDDO ! jc
      ENDDO ! jk

      ! Smoothing over 2*nsmooth+1 points
      IF (lsmootht) THEN
        ! Note: usually not necessary to smooth tendencies (defalult is .F.)
        IF (nlevsmootht >= nlev) nlevsmootht = nlev
        CALL smooth_vert(nlevsmootht,i_startidx,i_endidx,nsmooth,.FALSE.,  &
          &              p_fld%ddt_u_gwd_mgm(:,1:nlevsmootht+1,jb))
        CALL smooth_vert(nlevsmootht,i_startidx,i_endidx,nsmooth,.FALSE.,  &
          &              p_fld%ddt_v_gwd_mgm(:,1:nlevsmootht+1,jb))
        CALL smooth_vert(nlevsmootht,i_startidx,i_endidx,nsmooth,.FALSE.,  &
          &              p_fld%ddt_u_gwd_pmom_mgm(:,1:nlevsmootht+1,jb))
        CALL smooth_vert(nlevsmootht,i_startidx,i_endidx,nsmooth,.FALSE.,  &
          &              p_fld%ddt_v_gwd_pmom_mgm(:,1:nlevsmootht+1,jb))
        CALL smooth_vert(nlevsmootht,i_startidx,i_endidx,nsmooth,.FALSE.,  &
          &              p_fld%ddt_t_gwd_mgm(:,1:nlevsmootht+1,jb))
      ENDIF

    ENDDO ! jb

!$OMP END DO
!$OMP END PARALLEL

#ifndef __msgwam1d
    DEALLOCATE( ufl_hor_e, vfl_hor_e, ptfl_hor_e,    &
                ufldiv_hor, vfldiv_hor, ptfldiv_hor, &
                STAT=ist                             )
    IF (ist /= success) CALL finish ('tendency', 'Deallocation of auxiliary fields failed!')
    DEALLOCATE( upfl_hor_e, vpfl_hor_e,   &
                upfldiv_hor, vpfldiv_hor, &
                STAT=ist                  )
    IF (ist /= success) CALL finish ('tendency', 'Deallocation of auxiliary fields failed!')
#endif

    ! do heuristic checks of tendencies
    CALL heuristic_check_3d(p_patch%id, p_patch%nlev, p_patch%cells%start_block, p_patch%cells%end_block, &
      &  p_patch%cells%center(:,:)%lat, p_patch%cells%center(:,:)%lon, p_fld%ddt_u_gwd_mgm, 'ddt_u_gwd_mgm', 1)
    CALL heuristic_check_3d(p_patch%id, p_patch%nlev, p_patch%cells%start_block, p_patch%cells%end_block, &
      &  p_patch%cells%center(:,:)%lat, p_patch%cells%center(:,:)%lon, p_fld%ddt_v_gwd_mgm, 'ddt_v_gwd_mgm', 1)
    CALL heuristic_check_3d(p_patch%id, p_patch%nlev, p_patch%cells%start_block, p_patch%cells%end_block, &
      &  p_patch%cells%center(:,:)%lat, p_patch%cells%center(:,:)%lon, p_fld%ddt_t_gwd_mgm, 'ddt_t_gwd_mgm', 1)
    CALL heuristic_check_3d(p_patch%id, p_patch%nlev, p_patch%cells%start_block, p_patch%cells%end_block, &
      &  p_patch%cells%center(:,:)%lat, p_patch%cells%center(:,:)%lon, p_fld%ddt_u_gwd_pmom_mgm, 'ddt_u_gwd_pmom_mgm', 1)
    CALL heuristic_check_3d(p_patch%id, p_patch%nlev, p_patch%cells%start_block, p_patch%cells%end_block, &
      &  p_patch%cells%center(:,:)%lat, p_patch%cells%center(:,:)%lon, p_fld%ddt_v_gwd_pmom_mgm, 'ddt_v_gwd_pmom_mgm', 1)

  END SUBROUTINE tendency
#endif
END MODULE mo_nwp_msgwam_utils
