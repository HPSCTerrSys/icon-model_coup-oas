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

! Interface for the surface component of the turbulent mixing package (tmx)

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_tmx_surface_interface

  USE mo_kind, ONLY: wp, vp
  USE mo_exception, ONLY: finish
  USE mo_fortran_tools, ONLY: init, set_acc_async_queue
  USE mtime, ONLY: datetime
  USE mo_physical_constants, ONLY: grav, rgrav, tmelt, stbo, rhos, alf, cvv, clw, ci, vtmpc1, rd, nu
  USE mo_sea_ice_nml,        ONLY: Tf
  USE mo_aes_thermo, ONLY: &
    & lvc, lsc, &
    & sat_pres_water, sat_pres_ice, specific_humidity, dewpoint_temperature
  USE mo_turb_vdiff_params, ONLY: ckap
  USE mo_coupling_config,   ONLY: is_coupled_to_ocean
  USE mo_aes_phy_config,    ONLY: aes_phy_config  ! TODO: replace USE
  USE mo_tmx_field_class, ONLY: t_domain, isfc_oce, isfc_ice, isfc_lnd
  USE mo_turb_vdiff_params, ONLY: cchar, viscous_coeff
  USE mo_nh_testcases_nml,  ONLY: isrfc_type, shflx, lhflx

  ! If land is present, JSBACH is currently the only surface scheme supported by AES physcis package
#ifndef __NO_JSBACH__
  ! USE mo_jsb_interface, ONLY: jsbach_interface
  USE mo_jsb_interface,     ONLY: jsbach_interface, jsbach_get_var
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_land, update_sea_ice, &
    & compute_sfc_sat_spec_humidity, compute_sfc_fluxes, compute_sfc_roughness, &
    & compute_lw_rad_net, compute_sw_rad_net, compute_albedo, compute_energy_fluxes, &
    & compute_2m_temperature, compute_2m_humidity_and_dewpoint, compute_10m_wind

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_surface_interface'

CONTAINS

  ! Note: update_land is only called on OpenACC queue 1.
  !
  SUBROUTINE update_land(jg, domain, datetime_old, dtime, cvd, &
    & dz, pres_srf, ptemp, pq, pres_air, rsfl, ssfl, &
    & rlds, rvds_dir, rnds_dir, rpds_dir, rvds_dif, rnds_dif, rpds_dif, &
    & cosmu0, wind, wind10m, rho, co2, &
    ! out
    & tsfc, tsfc_rad, tsfc_eff, q_snocpymlt, qsat, &
    & alb_vis_dir, &
    & alb_vis_dif, &
    & alb_nir_dir, &
    & alb_nir_dif, &
    & kh, km, interp_fac_2m, interp_fac_10m, interp_fac_tsfc, &
    & co2flx)

    INTEGER, INTENT(in) :: &
      & jg
    TYPE(t_domain), INTENT(in), POINTER :: domain
    TYPE(datetime), INTENT(in), POINTER :: datetime_old ! date and time at beginning of this time step
    REAL(wp),       INTENT(in)          :: dtime
    REAL(wp),       INTENT(in)          :: cvd
    REAL(vp), INTENT(IN) :: &
      & dz        (:,:)     ! reference height in surface layer times 2
    REAL(wp), INTENT(IN) :: &
      & pres_srf  (:,:), &  ! surface pressure
      & ptemp     (:,:), &  ! temperature of lowest atmospheric level
      & pq        (:,:), &  ! humidity of lowest atmospheric level
      & pres_air  (:,:), &  ! pressure at lowest atmospheric level
      & rsfl      (:,:), &  ! surface rain flux, large-scale
      & ssfl      (:,:), &  ! surface snow flux, large-scale
      & rlds      (:,:), &
      & rvds_dir  (:,:), &
      & rnds_dir  (:,:), &
      & rpds_dir  (:,:), &
      & rvds_dif  (:,:), &
      & rnds_dif  (:,:), &
      & rpds_dif  (:,:), &
      & cosmu0    (:,:), &  ! cosine of zenith angle
      & wind      (:,:), &  ! wind speed at lowest level
      & wind10m   (:,:), &  ! wind speed at 10m
      & rho       (:,:), &  ! air density at surface
      & co2       (:,:)     ! CO2 at lowest atmospheric level

    REAL(wp), INTENT(out), OPTIONAL :: &
      & tsfc (:,:),           & ! new surface temperature
      & tsfc_rad (:,:),       & ! new surface radiative temperature
      & tsfc_eff (:,:),       & ! new surface effective temperature for rad heating
      & q_snocpymlt(:,:),     & ! heat used to melt snow on canopy
      & qsat(:,:),            & ! saturated specific humidity at surface
      & alb_vis_dir(:,:),     &
      & alb_vis_dif(:,:),     &
      & alb_nir_dir(:,:),     &
      & alb_nir_dif(:,:),     &
      & kh         (:,:),     & ! surface exchange coefficient (heat)
      & km         (:,:),     & ! surface exchange coefficient (momentum)
      & interp_fac_2m (:,:),  & ! neutral surface exchange coefficient (heat)
      & interp_fac_10m(:,:),  & ! neutral surface exchange coefficient (momentum)
      & interp_fac_tsfc(:,:), &
      & co2flx     (:,:)    ! CO2 flux into the atmosphere from natural sources

    INTEGER :: jb, jc, jcs, jce
    REAL(wp), DIMENSION(domain%nproma) :: &
      & dz_srf, rain_tmp, snow_tmp, &
      & rvds, rnds, rpds, fract_par_diffuse, &
      & t_acoef, t_bcoef, q_acoef, q_bcoef

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_land'

!$OMP PARALLEL
    CALL init(tsfc, lacc=.TRUE.)
    CALL init(qsat, lacc=.TRUE.)
    IF (PRESENT(km)) THEN
      CALL init(tsfc_rad, lacc=.TRUE.)
      CALL init(tsfc_eff, lacc=.TRUE.)
      CALL init(q_snocpymlt, lacc=.TRUE.)
      CALL init(alb_vis_dir, lacc=.TRUE.)
      CALL init(alb_vis_dif, lacc=.TRUE.)
      CALL init(alb_nir_dir, lacc=.TRUE.)
      CALL init(alb_nir_dif, lacc=.TRUE.)
      CALL init(kh, lacc=.TRUE.)
      CALL init(km, lacc=.TRUE.)
      CALL init(interp_fac_2m, lacc=.TRUE.)
      CALL init(interp_fac_10m, lacc=.TRUE.)
      CALL init(interp_fac_tsfc, lacc=.TRUE.)
      CALL init(co2flx, lacc=.TRUE.)
    END IF
!$OMP END PARALLEL

#ifndef __NO_JSBACH__

    !$ACC DATA CREATE(dz_srf, rain_tmp, snow_tmp, rvds, rnds, rpds, fract_par_diffuse, t_acoef, t_bcoef, q_acoef, q_bcoef) ASYNC(1)

!$OMP PARALLEL DO PRIVATE(jb, jcs, jce, jc, dz_srf, rain_tmp, snow_tmp, rvds, rnds, rpds, fract_par_diffuse, &
!$OMP                     t_acoef, t_bcoef, q_acoef, q_bcoef) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c

      jcs = domain%i_startidx_c(jb)
      jce = domain%i_endidx_c(jb)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jc = jcs, jce
        dz_srf(jc) = dz(jc,jb) * 0.5_wp
        rain_tmp(jc) = rsfl(jc,jb)
        snow_tmp(jc) = ssfl(jc,jb)
        rvds(jc) = rvds_dir(jc,jb) + rvds_dif(jc,jb)
        rnds(jc) = rnds_dir(jc,jb) + rnds_dif(jc,jb)
        rpds(jc) = rpds_dir(jc,jb) + rpds_dif(jc,jb)
        IF (rpds(jc) > 0._wp) THEN
          fract_par_diffuse(jc) = rpds_dif(jc,jb) / rpds(jc)
        ELSE
          fract_par_diffuse(jc) = 0._wp
        END IF
        t_acoef(jc) = 0._wp
        t_bcoef(jc) = cvd * ptemp(jc,jb) ! dry static energy
        q_acoef(jc) = 0._wp
        q_bcoef(jc) = pq(jc,jb)
      END DO
      !$ACC END PARALLEL LOOP

      IF (PRESENT(km)) THEN
        CALL jsbach_interface ( jg, jb, jcs, jce,                                         & ! in
          & datetime_old, dtime,                                                          & ! in
          & t_air             = ptemp(jcs:jce,jb),                                        & ! in
          & q_air             = pq(jcs:jce,jb),                                           & ! in
          & press_air         = pres_air(jcs:jce,jb),                                     & ! in
          & rain              = rain_tmp(jcs:jce),                                        & ! in
          & snow              = snow_tmp(jcs:jce),                                        & ! in
          & wind_air          = wind(jcs:jce,jb),                                         & ! in
          & wind_10m          = wind10m(jcs:jce,jb),                                      & ! in
          & lw_srf_down       = rlds(jcs:jce,jb),                                         & ! in
          & swvis_srf_down    = rvds(jcs:jce),                                            & ! in
          & swnir_srf_down    = rnds(jcs:jce),                                            & ! in
          & swpar_srf_down    = rpds(jcs:jce),                                            & ! in
          & fract_par_diffuse = fract_par_diffuse(jcs:jce),                               & ! in
          & dz_srf            = dz_srf(jcs:jce),                                          & ! in
          & press_srf         = pres_srf(jcs:jce,jb),                                     & ! in
          & rho_srf           = rho(jcs:jce,jb),                                          & ! in
          & t_acoef           = t_acoef(jcs:jce),                                         & ! in
          & t_bcoef           = t_bcoef(jcs:jce),                                         & ! in
          & q_acoef           = q_acoef(jcs:jce),                                         & ! in
          & q_bcoef           = q_bcoef(jcs:jce),                                         & ! in
          & cos_zenith_angle  = cosmu0(jcs:jce,jb),                                       & ! in
          & CO2_air           = co2(jcs:jce,jb),                                          & ! in
          & t_srf             = tsfc(jcs:jce,jb),                                         & ! out
          & t_srf_rad         = tsfc_rad(jcs:jce,jb),                                     & ! out
          & t_srf_eff         = tsfc_eff(jcs:jce,jb),                                     & ! out
          & q_snocpymlt       = q_snocpymlt(jcs:jce,jb),                                  & ! out
          & qsat_srf          = qsat(jcs:jce,jb),                                         & ! out
          & alb_vis_dir       = alb_vis_dir(jcs:jce,jb),                                  & ! out
          & alb_nir_dir       = alb_nir_dir(jcs:jce,jb),                                  & ! out
          & alb_vis_dif       = alb_vis_dif(jcs:jce,jb),                                  & ! out
          & alb_nir_dif       = alb_nir_dif(jcs:jce,jb),                                  & ! out
          & kh                = kh(jcs:jce,jb),                                           & ! out
          & km                = km(jcs:jce,jb),                                           & ! out
          & interp_fac_2m     = interp_fac_2m(jcs:jce,jb),                                & ! out
          & interp_fac_10m    = interp_fac_10m(jcs:jce,jb),                               & ! out
          & interp_fac_tsfc   = interp_fac_tsfc(jcs:jce,jb),                              & ! out
          & co2_flux          = co2flx(jcs:jce, jb)                                       & ! out
          ! & t_eff_srf         = ztsfc_lnd_eff(jcs:jce),                                   & ! out (T_s^eff) surface temp
          !                                                                                     ! (effective, for longwave rad)
          ! & s_srf             = zcpt_lnd(jcs:jce),                                        & ! out (s_s^star, for vdiff scheme)
          ! & fact_q_air        = pcair(jcs:jce),                                           & ! out
          ! & fact_qsat_srf     = pcsat(jcs:jce),                                           & ! out
          ! & evapotrans        = zevap_lnd(jcs:jce),                                       & ! out
          ! & latent_hflx       = zlhflx_lnd(jcs:jce),                                      & ! out
          ! & sensible_hflx     = zshflx_lnd(jcs:jce),                                      & ! out
          ! & grnd_hflx         = zgrnd_hflx(jcs:jce, idx_lnd),                             & ! out
          ! & grnd_hcap         = zgrnd_hcap(jcs:jce, idx_lnd),                             & ! out
          ! & rough_h_srf       = z0h_lnd(jcs:jce),                                         & ! out
          ! & rough_m_srf       = z0m_tile(jcs:jce, idx_lnd),                               & ! out
          ! & q_snocpymlt       = q_snocpymlt(jcs:jce),                                     & ! out
        )
      ELSE
        CALL jsbach_interface ( jg, jb, jcs, jce,                                         & ! in
          & datetime_old, dtime,                                                          & ! in
          & t_air             = ptemp(jcs:jce,jb),                                        & ! in
          & q_air             = pq(jcs:jce,jb),                                           & ! in
          & press_air         = pres_air(jcs:jce,jb),                                     & ! in
          & rain              = rain_tmp(jcs:jce),                                        & ! in
          & snow              = snow_tmp(jcs:jce),                                        & ! in
          & wind_air          = wind(jcs:jce,jb),                                         & ! in
          & wind_10m          = wind10m(jcs:jce,jb),                                      & ! in
          & lw_srf_down       = rlds(jcs:jce,jb),                                         & ! in
          & swvis_srf_down    = rvds(jcs:jce),                                            & ! in
          & swnir_srf_down    = rnds(jcs:jce),                                            & ! in
          & swpar_srf_down    = rpds(jcs:jce),                                            & ! in
          & fract_par_diffuse = fract_par_diffuse(jcs:jce),                               & ! in
          & dz_srf            = dz_srf(jcs:jce),                                          & ! in
          & press_srf         = pres_srf(jcs:jce,jb),                                     & ! in
          & rho_srf           = rho(jcs:jce,jb),                                          & ! in
          & t_acoef           = t_acoef(jcs:jce),                                         & ! in
          & t_bcoef           = t_bcoef(jcs:jce),                                         & ! in
          & q_acoef           = q_acoef(jcs:jce),                                         & ! in
          & q_bcoef           = q_bcoef(jcs:jce),                                         & ! in
          & cos_zenith_angle  = cosmu0(jcs:jce,jb),                                       & ! in
          & CO2_air           = co2(jcs:jce,jb),                                          & ! in
          & t_srf             = tsfc(jcs:jce,jb),                                         & ! out
          & qsat_srf          = qsat(jcs:jce,jb)                                          & ! out
        )
      END IF

    END DO
!$OMP END PARALLEL DO

    !$ACC END DATA

#else
    CALL finish(routine, "The JSBACH component is not activated")
#endif

  END SUBROUTINE update_land

  SUBROUTINE update_sea_ice(domain, dtime, &
    & nvalid, indices, &
    & old_tsfc, &
    & lwflx_net, swflx_net, lhflx, shflx, &
    & ssfl, ice_thickness, &
    & emissivity, &
    ! inout &
    & snow_thickness, &
    ! out &
    & new_tsfc, q_top, q_bot, &
    & albvisdir, albvisdif, albnirdir, albnirdif, &
    & opt_acc_async_queue )

  USE mo_ice_fast, ONLY: ice_fast

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), INTENT(in) :: dtime
    INTEGER,  INTENT(in)  :: &
      & nvalid(:), &
      & indices(:,:)
    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & old_tsfc,       &
      & lwflx_net,      &
      & swflx_net,      &
      & lhflx,          &
      & shflx,          &
      & ssfl,           &
      & ice_thickness,  &
      & emissivity

    REAL(wp), INTENT(inout), DIMENSION(:,:) :: &
      & snow_thickness

    REAL(wp), INTENT(out), DIMENSION(:,:) :: &
      & new_tsfc,  &  ! new surface temperature
      & q_top,     &
      & q_bot,     &
      & albvisdir, &
      & albvisdif, &
      & albnirdir, &
      & albnirdif

    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
    INTEGER :: acc_async_queue

    INTEGER :: jb, jc, jls, jcs, jce, kice
    REAL(wp), DIMENSION(domain%nproma) :: &
      & Tfw, nonsolar_flux, dnonsolar_flux_dt, &
      & T1, T2 !< Dummies, not used

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_sea_ice'

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

    !$ACC DATA CREATE(Tfw, nonsolar_flux, dnonsolar_flux_dt, T1, T2) ASYNC(acc_async_queue)

!$OMP PARALLEL
    CALL init(new_tsfc, Tf,      lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(q_top,             lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(q_bot,             lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(albvisdir,         lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(albvisdif,         lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(albnirdir,         lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(albnirdif,         lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL

    kice = 1

!$OMP PARALLEL DO PRIVATE(jb, jcs, jce, jls, jc, Tfw, T1, T2, nonsolar_flux, dnonsolar_flux_dt) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      jcs = domain%i_startidx_c(jb)
      jce = domain%i_endidx_c(jb)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue)
      DO jc = jcs, jce
        Tfw(jc) = Tf
        T1(jc) = 0._wp
        T2(jc) = 0._wp
        nonsolar_flux(jc) = 0._wp
        dnonsolar_flux_dt(jc) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(jc)
      DO jls = 1, nvalid(jb)
        jc = indices(jls,jb)
        nonsolar_flux(jc) = lwflx_net(jc,jb) + lhflx(jc,jb) + shflx(jc,jb)
        dnonsolar_flux_dt(jc) = -4._wp * emissivity(jc,jb) * stbo * old_tsfc(jc,jb)**3._wp
        new_tsfc(jc,jb) = old_tsfc(jc,jb) - tmelt
      END DO
      !$ACC END PARALLEL LOOP

      CALL ice_fast(jcs, jce, domain%nproma, kice, dtime, &
        &   new_tsfc(:,jb),        & ! inout
        &   T1(:),                 & ! inout, dummy
        &   T2(:),                 & ! inout, dummy
        &   ice_thickness(:,jb),   & ! in
        &   snow_thickness(:,jb),  & ! in
        &   q_top(:,jb),           & ! out
        &   q_bot(:,jb),           & ! out
        &   swflx_net(:,jb),       & ! in
        &   nonsolar_flux(:),      & ! in
        &   dnonsolar_flux_dt(:),  & ! in
        &   Tfw(:),                & ! in
        &   albvisdir(:,jb),       & ! out
        &   albvisdif(:,jb),       & ! out
        &   albnirdir(:,jb),       & ! out
        &   albnirdif(:,jb),       & ! out
        &   lacc=.TRUE.,           & ! in
        &   opt_acc_async_queue=acc_async_queue) ! in

      ! Update the thickness of snow on ice in atmosphere only simulation.
      ! In coupled experiments this is done by the ocean model in either
      ! ice_growth_zerolayer or ice_growth_winton.
      IF ( .NOT. is_coupled_to_ocean() ) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(jc)
        DO jls = 1, nvalid(jb)
          jc = indices(jls,jb)
          ! Snowfall on ice - no ice => no snow
          IF (ice_thickness(jc,jb) > 0._wp) THEN
            ! Snow only falls when it's below freezing
            IF (new_tsfc(jc,jb) < 0._wp) THEN
              snow_thickness(jc,jb) = snow_thickness(jc,jb) + ssfl(jc,jb) * dtime / rhos
            ENDIF
            ! Snow melt
            snow_thickness(jc,jb) = snow_thickness(jc,jb) - MIN( q_top(jc,jb) * dtime / (alf * rhos), snow_thickness(jc,jb) )
          ELSE
            snow_thickness(jc,jb) = 0._wp
          ENDIF
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue)
      DO jls = 1, nvalid(jb)
        jc = indices(jls,jb)
        new_tsfc(jc,jb) = new_tsfc(jc,jb) + tmelt
      END DO
      !$ACC END PARALLEL LOOP

    END DO

!$OMP END PARALLEL DO

    !$ACC END DATA

  END SUBROUTINE update_sea_ice

  SUBROUTINE compute_lw_rad_net( &
    & domain,                    &
    & nvalid, indices,           &
    & emissivity,                &
    & rlds,                      &
    & tsfc,                      &
    & lwfl_net,                  &
    & opt_acc_async_queue        &
    & )

    USE mo_physical_constants,ONLY: stbo

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & emissivity ,     &
      & rlds,            &
      & tsfc
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: lwfl_net
    !
    ! Optional ACC queue
    !
    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
    INTEGER :: acc_async_queue

    INTEGER  :: jb, jl, jls, js

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_lw_rad_net'

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

!$OMP PARALLEL
    CALL init(lwfl_net, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        lwfl_net(js,jb) = emissivity(js,jb) * (rlds(js,jb) - stbo * tsfc(js,jb)**4._wp)
      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_lw_rad_net
  !
  !=================================================================
  !
  SUBROUTINE compute_sw_rad_net( &
    & domain,                    &
    & nvalid, indices,           &
    & rvds_dir,                  &
    & rvds_dif,                  &
    & rnds_dir,                  &
    & rnds_dif,                  &
    & alb_vis_dir,               &
    & alb_vis_dif,               &
    & alb_nir_dir,               &
    & alb_nir_dif,               &
    & swfl_net,                  &
    & opt_acc_async_queue        &
    & )


    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & rvds_dir,            &
      & rvds_dif,            &
      & rnds_dir,            &
      & rnds_dif,            &
      & alb_vis_dir,         &
      & alb_vis_dif,         &
      & alb_nir_dir,         &
      & alb_nir_dif
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: swfl_net
    !
    ! Optional ACC queue
    !
    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
    INTEGER :: acc_async_queue

    INTEGER  :: jb, jl, jls, js

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sw_rad_net'

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

!$OMP PARALLEL
    CALL init(swfl_net, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        swfl_net(js,jb) = &
          & rvds_dif(js,jb) * (1._wp - alb_vis_dif(js,jb)) + &
          & rvds_dir(js,jb) * (1._wp - alb_vis_dir(js,jb)) + &
          & rnds_dif(js,jb) * (1._wp - alb_nir_dif(js,jb)) + &
          & rnds_dir(js,jb) * (1._wp - alb_nir_dir(js,jb))
      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_sw_rad_net
  !
  !=================================================================
  !
  SUBROUTINE compute_sfc_roughness( &
    & linit,                   &
    & domain,                  &
    & isfc,                    &
    & nvalid, indices,         &
    & rough_min, rough_oce, rough_ice, wind_rel, km,  &
    & rough_h, rough_m,        &
    & opt_acc_async_queue      &
    & )

    USE mo_turb_vdiff_params, ONLY: cchar

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    INTEGER,         INTENT(in) :: isfc
    !
    ! Input variables
    !
    LOGICAL :: linit
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), INTENT(in) :: rough_min, rough_oce, rough_ice
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & wind_rel,            &
      & km
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: rough_h, rough_m
    !
    ! Optional ACC queue
    !
    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
    INTEGER :: acc_async_queue

    INTEGER  :: jb, jl, jls, js

    REAL(wp), POINTER :: jsb_rough_m(:,:) => NULL(), jsb_rough_h(:,:) => NULL()
    REAL(wp) :: rough_tmp

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sfc_roughness'

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

!$OMP PARALLEL
    CALL init(rough_h, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(rough_m, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jl, jls, js, rough_tmp) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR PRIVATE(js, rough_tmp) ASYNC(acc_async_queue)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)

        IF (isfc == isfc_oce) THEN
          IF (linit) THEN
            rough_tmp = rough_oce
          ELSE
            rough_tmp = (wind_rel(js,jb)**2 * km(js,jb) * cchar * rgrav) &
                        + (viscous_coeff * MIN(0.01_wp, nu / (km(js,jb)**0.5_wp * wind_rel(js,jb)) ))
          END IF
          rough_tmp = MAX(rough_min, rough_tmp)
          rough_m(js,jb) = rough_tmp
          rough_h(js,jb) = rough_tmp
        ELSE IF (isfc == isfc_ice) THEN
          ! Nothing to do
          rough_m(js,jb) = rough_ice
          rough_h(js,jb) = rough_ice
        END IF

      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

    IF (isfc == isfc_lnd) THEN
#ifndef __NO_JSBACH__
      CALL jsbach_get_var('turb_rough_m', domain%patch%id, ptr2d=jsb_rough_m, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
      CALL jsbach_get_var('turb_rough_h', domain%patch%id, ptr2d=jsb_rough_h, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c,domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
        DO jls = 1, nvalid(jb)
          js = indices(jls,jb)
          rough_m(js,jb) = jsb_rough_m(js,jb)
          rough_h(js,jb) = jsb_rough_h(js,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO
      NULLIFY(jsb_rough_m, jsb_rough_h)
#else
      CALL finish(routine, "The JSBACH component is not activated")
#endif
    END IF

  END SUBROUTINE compute_sfc_roughness
  !
  !=================================================================
  !
  SUBROUTINE compute_sfc_sat_spec_humidity( &
    & linit, domain, isfc,     &
    & nvalid, indices,         &
    & ppsfc, ptsfc,            &
    & qsat,                    &
    & opt_acc_async_queue      &
    )

    ! Domain information
    TYPE(t_domain),        INTENT(in), POINTER :: domain
    INTEGER,         INTENT(in) :: isfc
    !
    ! Input variables
    !
    LOGICAL, INTENT(in) :: linit
    INTEGER, INTENT(in) :: &
      & nvalid(:),         &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      ptsfc ,     &
      ppsfc
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: qsat
    !
    ! Optional ACC queue
    !
    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
    INTEGER :: acc_async_queue

    INTEGER  :: jb, jls, js
    REAL(wp), POINTER :: jsb_qsat(:,:) => NULL()

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sfc_sat_spec_humidity'

#ifdef __NO_JSBACH__
    IF (isfc == isfc_lnd) CALL finish(routine, "The JSBACH component is not activated")
#endif

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c

      ! DO jl = 1, domain%nproma
      !   qsat(jl,jb) = 0._wp
      ! END DO

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        IF (isfc == isfc_oce) THEN
          qsat(js,jb) = specific_humidity(sat_pres_water(ptsfc(js,jb)),ppsfc(js,jb))
        ELSE IF (isfc == isfc_ice) THEN
          qsat(js,jb) = specific_humidity(sat_pres_ice(ptsfc(js,jb)),ppsfc(js,jb))
        ELSE IF (isfc == isfc_lnd .AND. linit) THEN
          qsat(js,jb) = specific_humidity(sat_pres_water(ptsfc(js,jb)),ppsfc(js,jb))
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

    IF (isfc == isfc_lnd .AND. .NOT. linit) THEN
#ifndef __NO_JSBACH__
      CALL jsbach_get_var('seb_qsat_star', domain%patch%id, ptr2d=jsb_qsat, &
        lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c,domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
        DO jls = 1, nvalid(jb)
          js = indices(jls,jb)
          qsat(js,jb) = jsb_qsat(js,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO
      NULLIFY(jsb_qsat)
#endif
    END IF

  END SUBROUTINE compute_sfc_sat_spec_humidity
  !
  !=================================================================
  !
  SUBROUTINE compute_sfc_fluxes( &
    & domain,                  &
    & isfc,                    &
    & nvalid, indices,         &
    & wind_g,                  &
    & cvd,                     &
    ! & ua, va, thetam1, qm1, wind, rho, qsat_sfc, theta_sfc, kh, km,  &
    & ua, va, ta, qm1, wind_rel, oce_u, oce_v, ice_u, ice_v, rho, qsat_sfc, t_sfc, kh, km,  &
    & evapotrans, latent_hflx, sensible_hflx, ustress, vstress,  &
    & opt_acc_async_queue      &
    & )

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    INTEGER,         INTENT(in) :: isfc
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), INTENT(in) :: cvd, wind_g
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & ua, &
      & va, &
      ! & thetam1, &
      & ta, &
      & qm1, &
      & kh ,     &
      & km,      &
      & rho,     &
      & qsat_sfc, &
      ! & theta_sfc, &
      & t_sfc, &
      & wind_rel, &
      & oce_u, &
      & oce_v, &
      & ice_u, &
      & ice_v
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: evapotrans, latent_hflx, sensible_hflx, ustress, vstress
    !
    ! Optional ACC queue
    !
    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
    INTEGER :: acc_async_queue

    INTEGER  :: jb, jl, jls, js

#ifndef __NO_JSBACH__
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & jsb_evapotrans_ptr => NULL(), jsb_latent_hflx_ptr => NULL(), jsb_sensible_hflx_ptr => NULL()
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sfc_fluxes'

    ! CALL message(routine, 'Start')

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

!$OMP PARALLEL
      CALL init(evapotrans,    lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
      CALL init(latent_hflx,   lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
      CALL init(sensible_hflx, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
      CALL init(ustress,       lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
      CALL init(vstress,       lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL

    IF (isrfc_type == 1) THEN
!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c,domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
        DO jls = 1, nvalid(jb)
          js = indices(jls,jb)
          latent_hflx(js,jb)   = -lhflx * rho(js,jb) * (lvc+(cvv-clw)*t_sfc(js,jb))
          evapotrans(js,jb)    = -lhflx * rho(js,jb)
          sensible_hflx(js,jb) = -shflx * cvd * rho(js,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO
      RETURN
    END IF

    IF (isfc == isfc_lnd) THEN
#ifndef __NO_JSBACH__
      CALL jsbach_get_var('hydro_evapotrans',  domain%patch%id, ptr2d=jsb_evapotrans_ptr, &
        & lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
      CALL jsbach_get_var('seb_latent_hflx',   domain%patch%id, ptr2d=jsb_latent_hflx_ptr, &
        & lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
      CALL jsbach_get_var('seb_sensible_hflx', domain%patch%id, ptr2d=jsb_sensible_hflx_ptr, &
        & lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
#else
      CALL finish(routine, "The JSBACH component is not activated")
#endif
    END IF

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        ! TODO: is the treatment of surface ocean current correct (cf. vdiff code)
        IF (isfc == isfc_oce) THEN
          ! FIXME: correct to use wind_g with relative wind speed?
          evapotrans(js,jb) = &
            & rho(js,jb) * kh(js,jb) * (SQRT(wind_g**2._wp + wind_rel(js,jb)**2._wp)) * (qm1(js,jb) - qsat_sfc(js,jb))
          latent_hflx(js,jb) = evapotrans(js,jb) * (lvc+(cvv-clw)*t_sfc(js,jb))
          sensible_hflx(js,jb) = cvd * rho(js,jb) * kh(js,jb) * &
            & (SQRT(wind_g**2._wp + wind_rel(js,jb)**2._wp))* (ta(js,jb) - t_sfc(js,jb))
          ustress(js,jb) = rho(js,jb) * km(js,jb) * wind_rel(js,jb) * (ua(js,jb) - oce_u(js,jb))
          vstress(js,jb) = rho(js,jb) * km(js,jb) * wind_rel(js,jb) * (va(js,jb) - oce_v(js,jb))
        ELSE IF (isfc == isfc_ice) THEN
          evapotrans(js,jb) = rho(js,jb) * wind_rel(js,jb) * kh(js,jb) * (qm1(js,jb) - qsat_sfc(js,jb))
          latent_hflx(js,jb) = evapotrans(js,jb) * (lsc+(cvv-ci)*t_sfc(js,jb))
          sensible_hflx(js,jb) = cvd * rho(js,jb) * wind_rel(js,jb) * kh(js,jb) * (ta(js,jb) - t_sfc(js,jb))
          ustress(js,jb) = rho(js,jb) * km(js,jb) * wind_rel(js,jb) * (ua(js,jb) - ice_u(js,jb))
          vstress(js,jb) = rho(js,jb) * km(js,jb) * wind_rel(js,jb) * (va(js,jb) - ice_v(js,jb))
        ELSE IF (isfc == isfc_lnd) THEN
#ifndef __NO_JSBACH__
          evapotrans(js,jb) = jsb_evapotrans_ptr(js,jb)
          latent_hflx(js,jb) = jsb_latent_hflx_ptr(js,jb)
          sensible_hflx(js,jb) = jsb_sensible_hflx_ptr(js,jb)
          ustress(js,jb) = rho(js,jb) * km(js,jb) * wind_rel(js,jb) * ua(js,jb)
          vstress(js,jb) = rho(js,jb) * km(js,jb) * wind_rel(js,jb) * va(js,jb)
#endif
        END IF
      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

#ifndef __NO_JSBACH__
    IF (isfc == isfc_lnd) THEN
      NULLIFY(jsb_evapotrans_ptr, jsb_latent_hflx_ptr, jsb_sensible_hflx_ptr)
    END IF
#endif

  END SUBROUTINE compute_sfc_fluxes
  !
  !=================================================================
  !
  SUBROUTINE compute_energy_fluxes( &
    & domain,     &
    & cvv, cvd,   &
    & shfl,       &
    & evapotrans, &
    & ta,         &
    & rho,        &
    & ufts,       &
    & ufvs)

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    REAL(wp), INTENT(in) :: cvv, cvd
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & shfl,       &
      & evapotrans, &
      & ta,         &
      & rho
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & ufts, &
      & ufvs

    INTEGER :: ic, ib

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_energy_fluxes'

!$OMP PARALLEL
    CALL init(ufts, lacc=.TRUE.)
    CALL init(ufvs, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(ib, ic) ICON_OMP_DEFAULT_SCHEDULE
    DO ib = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic = domain%i_startidx_c(ib), domain%i_endidx_c(ib)
        ufts(ic,ib) = shfl(ic,ib)
        ufvs(ic,ib) = ta(ic,ib) * evapotrans(ic,ib) * (cvv - cvd)
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_energy_fluxes
  !
  !=================================================================
  !
  SUBROUTINE compute_albedo(     &
    & domain,                    &
    & nvalid, indices,           &
    & rsds,                      &
    & rvds_dir,                  &
    & rvds_dif,                  &
    & rnds_dir,                  &
    & rnds_dif,                  &
    & alb_vis_dir,               &
    & alb_vis_dif,               &
    & alb_nir_dir,               &
    & alb_nir_dif,               &
    & albedo,                    &
    & opt_acc_async_queue        &
    & )


    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & rsds,                &
      & rvds_dir,            &
      & rvds_dif,            &
      & rnds_dir,            &
      & rnds_dif,            &
      & alb_vis_dir,         &
      & alb_vis_dif,         &
      & alb_nir_dir,         &
      & alb_nir_dif
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: albedo
    !
    ! Optional ACC queue
    !
    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
    INTEGER :: acc_async_queue

    INTEGER  :: jb, jls, js
    REAL(wp) :: zalbvis, zalbnir, rvds, rnds

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_albedo'

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

!$OMP PARALLEL
    CALL init(albedo, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js, zalbvis, zalbnir, rvds, rnds) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)

        zalbvis = 0._wp
        rvds = rvds_dir(js,jb) + rvds_dif(js,jb)
        IF(rvds > 0._wp) THEN
          zalbvis = &
            & (alb_vis_dir(js,jb) * rvds_dir(js,jb) + alb_vis_dif(js,jb) * rvds_dif(js,jb)) &
            & / rvds
        END IF
        zalbnir = 0._wp
        rnds = rnds_dir(js,jb) + rnds_dif(js,jb)
        IF(rnds > 0._wp) THEN
          zalbnir = &
            & (alb_nir_dir(js,jb) * rnds_dir(js,jb) + alb_nir_dif(js,jb) * rnds_dif(js,jb)) &
            & / rnds
        END IF
        IF(rsds(js,jb) > 0._wp) THEN
          albedo(js,jb) = &
            & (zalbvis * rvds + zalbnir * rnds) &
            & / rsds(js,jb)
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_albedo
  !
  !=================================================================
  !
  SUBROUTINE compute_10m_wind( &
    & domain,                  &
    & nvalid, indices,         &
    & ua, va, interp_fac,      &
    & u10m, v10m, wind10m,     &
    & opt_acc_async_queue)

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & interp_fac, ua, va
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & u10m, v10m, wind10m
    !
    ! Optional ACC queue
    !
    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
    INTEGER :: acc_async_queue

    INTEGER :: jb, jls, js

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_10m_wind'

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

!$OMP PARALLEL
    CALL init(u10m,    lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(v10m,    lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(wind10m, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        u10m(js,jb) = interp_fac(js,jb) * ua(js,jb)
        v10m(js,jb) = interp_fac(js,jb) * va(js,jb)
        wind10m(js,jb) = SQRT(u10m(js,jb)**2._wp + (v10m(js,jb)**2._wp))
      END DO
    !$ACC END PARALLEL LOOP
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_10m_wind
  !
  !=================================================================
  !
  SUBROUTINE compute_2m_temperature( &
    & domain, isfc,                  &
    & nvalid, indices, tatm, tsfc,   &
    & interp_fac, interp_tsfc,       &
    & t2m)

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:),        &
      & isfc
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & tatm, tsfc, interp_fac, interp_tsfc
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & t2m

    INTEGER :: jb, jls, js

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_2m_temperature'

!$OMP PARALLEL
    CALL init(t2m, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        t2m(js,jb) = tsfc(js,jb) + interp_fac(js,jb) * tatm(js,jb) - interp_tsfc(js,jb)
      END DO
    !$ACC END PARALLEL LOOP
    ENDDO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_2m_temperature
  !
  !=================================================================
  !
  SUBROUTINE compute_2m_humidity_and_dewpoint( &
    & domain, nvalid, indices,    &
    & patm, psfc,                 &
    & tatm, t2m,                  &
    & qv_atm, qc_atm, qi_atm,     &
    & hus2m, dew2m)

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & patm, psfc,             & ! pressure in lowest model level and at surface
      & tatm, t2m,              & ! temperature in lowest model level and at surface
      & qv_atm, qc_atm, qi_atm    ! water vapor, cloud water and cloud ice in lowest model level
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & hus2m,                  & ! specific humidity at 2 meter
      & dew2m                     ! dewpoint temperature at 2 meter

    INTEGER :: jb, jls, js
    REAL(wp) :: &
      & qsat_atm, qsat_2m, & ! saturated humidity in lowest model level and at 2 meter
      & qrel_atm,          & ! relative humidity in lowest model level
      & pres2m               ! pressure at 2 meter

    REAL(wp), PARAMETER :: zephum = 0.05_wp ! epsilon for rel. humidity
    !$ACC DECLARE COPYIN(zephum)

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_2m_humidity'

!$OMP PARALLEL
    CALL init(hus2m, lacc=.TRUE.)
!$OMP END PARALLEL

    ! Note: Below it is assumed that relative humidity is constant with height. This should
    ! be revisited, see comments in https://gitlab.dkrz.de/icon/icon-mpim/-/merge_requests/341
    !
!$OMP PARALLEL DO PRIVATE(jb, jls, js, qsat_atm, qrel_atm, pres2m, qsat_2m) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
      !$ACC   PRIVATE(js, qsat_atm, qrel_atm, pres2m, qsat_2m)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        qsat_atm = specific_humidity(sat_pres_water(tatm(js,jb)),patm(js,jb))
        qrel_atm = MAX(zephum, qv_atm(js,jb) / qsat_atm)
        pres2m = psfc(js,jb) * &
          &  (1._wp - 2._wp * grav / ( rd * t2m(js,jb) * (1._wp + vtmpc1 * qv_atm(js,jb) - qc_atm(js,jb) - qi_atm(js,jb))))
        qsat_2m = specific_humidity(sat_pres_water(t2m(js,jb)),pres2m)
        hus2m(js,jb) = qrel_atm * qsat_2m
        dew2m(js,jb) = dewpoint_temperature(t2m(js,jb), hus2m(js,jb), pres2m)
      END DO
      !$ACC END PARALLEL LOOP
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_2m_humidity_and_dewpoint
  !
  !=================================================================
  !

END MODULE mo_tmx_surface_interface
