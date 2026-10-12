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

! Subroutine to initialize the wave test case

MODULE mo_wave_adv_exp

  USE mo_kind,                 ONLY: wp
  USE mo_model_domain,         ONLY: t_patch
  USE mo_wave_forcing_types,   ONLY: t_wave_forcing
  USE mo_wave_config,          ONLY: t_wave_config
  USE mo_physical_constants,   ONLY: earth_radius
  USE mo_math_constants,       ONLY: pi, pi2, deg2rad, dbl_eps
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, min_rlcell
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_wave_td_update,       ONLY: update_ice_free_mask

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_analytic_forcing

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_adv_exp'

  !--------------------------------------------------------------------

CONTAINS

  !
  ! Prescribes a time-constant, analytically defined wind field with a
  ! fixed direction, from which an analytic JONSWAP-type wave spectrum
  ! is constructed.
  ! The resulting spectral wave variance field may be used in
  ! transport-only tests for evaluating the horizontal transport
  ! operator.
  !
  ! This routine initializes the following forcing fields:
  ! u10m, v10m, sp10m, dir10m, sea_ice_c, ice_free_mask
  !
  SUBROUTINE init_analytic_forcing(p_patch, wave_config, p_forcing)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//'::init_analytic_forcing'

    TYPE(t_patch),        INTENT(IN)         :: p_patch
    TYPE(t_wave_config),  TARGET, INTENT(IN) :: wave_config
    TYPE(t_wave_forcing), INTENT(INOUT)      :: p_forcing

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: jc, jb
    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx

    REAL(wp) :: sin_tmp, cos_tmp, zlat, zlon, d1, r

    REAL(wp) :: theta_w  ! wind direction measured clockwise from true north in radians
    REAL(wp) :: lambda0, phi0  ! coordinates of center point
    REAL(wp) :: RR ! horizontal half width on a unit sphere

    wc => wave_config

    lambda0 = wc%peak_lon  * deg2rad
    phi0    = wc%peak_lat  * deg2rad
    theta_w = wc%dir_wsp10 * deg2rad
    RR      = wave_config%r_wsp10 / earth_radius

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,zlon,zlat,sin_tmp,cos_tmp,r,d1)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,      &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        zlon = p_patch%cells%center(jc,jb)%lon
        zlat = p_patch%cells%center(jc,jb)%lat

        sin_tmp = SIN(zlat) * SIN(phi0)
        cos_tmp = COS(zlat) * COS(phi0)
        r  = ACOS (sin_tmp + cos_tmp*COS(zlon-lambda0))    ! great circle distance on a unit sphere
        d1 = MIN( 1._wp, (r/RR) )

        ! cosine bell distribution, following Williamson (1992)
        p_forcing%sp10m(jc,jb) = 0.5_wp * wc%peak_wsp10 * (1._wp + COS(pi*d1))
        p_forcing%u10m(jc,jb)  = p_forcing%sp10m(jc,jb) * sin(theta_w)
        p_forcing%v10m(jc,jb)  = p_forcing%sp10m(jc,jb) * cos(theta_w)
        p_forcing%dir10m(jc,jb)= theta_w

        ! no seaice
        p_forcing%sea_ice_c = 0._wp
      END DO ! cell loop
    END DO
!$OMP END DO
!$OMP END PARALLEL

    ! compute mask of ice-free points
    CALL update_ice_free_mask(                       &
      &     p_patch       = p_patch,                 & ! IN
      &     sea_ice_c     = p_forcing%sea_ice_c,     & ! IN
      &     ice_free_mask = p_forcing%ice_free_mask_c) ! OUT

  END SUBROUTINE init_analytic_forcing

END MODULE mo_wave_adv_exp
