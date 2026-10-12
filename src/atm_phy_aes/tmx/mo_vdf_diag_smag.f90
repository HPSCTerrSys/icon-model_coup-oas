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

MODULE mo_vdf_diag_smag

  USE mo_kind,              ONLY: wp, vp, i1
  USE mo_exception,         ONLY: message, finish
  USE mo_tmx_process_class, ONLY: t_tmx_process
  USE mo_tmx_field_class,   ONLY: t_tmx_field, t_domain, isfc_oce, isfc_ice, isfc_lnd
  USE mo_util_string,       ONLY: real2string
  USE mo_model_domain,      ONLY: t_patch
  USE mo_intp_data_strc,    ONLY: t_int_state, p_int_state
  USE mo_nonhydro_types,    ONLY: t_nh_metrics
  USE mo_nonhydro_state,    ONLY: p_nh_state
  USE mo_aes_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd
  ! USE mo_physical_constants,ONLY: grav, rd, vtmpc1, p0ref, rgrav, alv, als
  USE mo_physical_constants,ONLY: grav, rd, vtmpc1, rgrav
  USE mo_turb_vdiff_params, ONLY: ckap
  USE mo_math_constants,    ONLY: pi_2, ln2
  USE mo_impl_constants,    ONLY: min_rlcell, min_rledge_int, min_rlcell_int, min_rlvert_int
  USE mo_nh_vert_interp_les,ONLY: brunt_vaisala_freq, vert_intp_full2half_cell_3d
  USE mo_intp,              ONLY: cells2verts_scalar, cells2edges_scalar
  USE mo_sync,              ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array
  USE mo_loopindices,       ONLY: get_indices_e, get_indices_c
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_intp_rbf,          ONLY: rbf_vec_interpol_edge
  USE mo_fortran_tools,     ONLY: init, set_acc_async_queue

  USE mo_aes_thermo,        ONLY: potential_temperature, sat_pres_water, sat_pres_ice, specific_humidity
  ! USE mo_jsb_interface,     ONLY: jsbach_get_var

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sfc_exchange_coefficients, compute_sfc_exchange_coefficients, &
    &       compute_wind_speed, &
    &       compute_atm_potential_temperature, compute_moist_richardson, &
    ! &       compute_sfc_sat_spec_humidity, compute_sfc_fluxes, compute_sfc_roughness &
    &       compute_sfc_potential_temperature, compute_sfc_density

  !Parameters for surface layer parameterizations: From Zeng_etal 1997 J. Clim
  REAL(wp), PARAMETER :: bsm = 5.0_wp  !Businger Stable Momentum
  REAL(wp), PARAMETER :: bum = 16._wp  !Businger Untable Momentum
  REAL(wp), PARAMETER :: bsh = 5.0_wp  !Businger Stable Heat
  REAL(wp), PARAMETER :: buh = 16._wp  !Businger Untable Heat

  CHARACTER(len=*), PARAMETER :: modname = 'mo_vdf_diag_smag'

CONTAINS

  !
  !=================================================================
  !
  ! Subroutines for diagnostics
  !
  ! Note: will be re-factored
  !
  SUBROUTINE compute_sfc_density(                 &
    & domain,                                            &
    & nvalid, indices,                                   &
    & qsat_sfc, ppsfc, ptsfc,  &
    & rhos, &
    & opt_acc_async_queue )

    ! Domain information
    TYPE(t_domain),        INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & qsat_sfc, &
      ppsfc,     &
      ptsfc
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: rhos
    !
    ! Optional ACC queue
    !
    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
    INTEGER :: acc_async_queue

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sfc_density'

    INTEGER  :: jb, jls, js

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

!$OMP PARALLEL
    CALL init(rhos, 1._wp, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        rhos(js,jb) = ppsfc(js,jb) / (rd * ptsfc(js,jb) * (1._wp + vtmpc1 * qsat_sfc(js,jb)))
        ! TODO: in vdiff sfc_exchange_coeff 1._wp+vtmpc1*pqm1_b(jl)-pqxm1_b(jl) is used
        !       and temperature at lowest model level instead of ptsfc !?

      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_sfc_density  !
!
!=================================================================
!
  SUBROUTINE compute_wind_speed( &
    & domain,                  &
    & nvalid, indices,         &
    & min_sfc_wind, isfc,      &
    & pum1, pvm1,              &
    & oce_u, oce_v,            &
    & ice_u, ice_v,            &
    & wind_rel,                &
    & opt_acc_async_queue      &
    )

    ! Domain information
    TYPE(t_domain),        INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), INTENT(in) :: min_sfc_wind
    INTEGER, INTENT(IN) :: isfc
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      pum1 ,     &
      pvm1 ,     &
      oce_u,     &
      oce_v,     &
      ice_u,     &
      ice_v
    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: wind_rel

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_wind_speed'

    INTEGER :: jb, jls, js
    INTEGER :: acc_async_queue

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

!$OMP PARALLEL
    ! Initialize with min_sfc_wind to prevent division by zero later (from jsbach)
    CALL init(wind_rel, min_sfc_wind, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        IF (isfc == isfc_oce) THEN ! open ocean
          wind_rel(js,jb) = MAX( min_sfc_wind, SQRT( (pum1(js,jb) - oce_u(js,jb))**2._wp + (pvm1(js,jb) - oce_v(js,jb))**2._wp ) )
        ELSEIF (isfc == isfc_ice) THEN ! sea ice
          wind_rel(js,jb) = MAX( min_sfc_wind, SQRT( (pum1(js,jb) - ice_u(js,jb))**2._wp + (pvm1(js,jb) - ice_v(js,jb))**2._wp ) )
        ELSE ! land
          wind_rel(js,jb) = MAX( min_sfc_wind, SQRT(pum1(js,jb)**2._wp + pvm1(js,jb)**2._wp) )
        END IF
      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_wind_speed
  !
  !=================================================================
  !
  SUBROUTINE compute_atm_potential_temperature( &
    & domain,                  &
    & ptm1, ptvm1, papm1,      &
    & theta, thetav            &
    )

    ! Domain information
    TYPE(t_domain),        INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
    & ptm1 , &
    & ptvm1, &
    & papm1
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: theta, thetav

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_atm_potential_temperature'

    INTEGER  :: jb, jc

!$OMP PARALLEL
    CALL init(theta, lacc=.TRUE.)
    CALL init(thetav, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
        theta (jc,jb) = potential_temperature(ptm1(jc,jb),  papm1(jc,jb))
        thetav(jc,jb) = potential_temperature(ptvm1(jc,jb), papm1(jc,jb))
      END DO !jc
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_atm_potential_temperature
  !
  !=================================================================
  !
  SUBROUTINE compute_sfc_potential_temperature( &
    & domain,                  &
    & nvalid, indices,         &
    & ppsfc, ptsfc, qsat,      &
    & theta, thetav,           &
    & opt_acc_async_queue      &
    )

    ! Domain information
    TYPE(t_domain),        INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      ptsfc,     &
      qsat,      &
      ppsfc
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: theta, thetav
    !
    ! Optional ACC queue
    !
    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
    INTEGER :: acc_async_queue

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sfc_potential_temperature'

    INTEGER  :: jb, jls, js

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

!$OMP PARALLEL
    CALL init(theta,  lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(thetav, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        theta(js,jb)  = potential_temperature(ptsfc(js,jb), ppsfc(js,jb))
        thetav(js,jb) = theta(js,jb) * (1._wp + vtmpc1 * qsat(js,jb))
      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_sfc_potential_temperature
  !
  !=================================================================
  !
  SUBROUTINE compute_moist_richardson(     &
    & domain,                    &
    & nvalid, indices,           &
    & fsl, zf, &
    & thetav_atm, thetav_sfc, wind, &
    & richardson_number, &
    & opt_acc_async_queue )

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), INTENT(in) :: &
      & fsl
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & zf, thetav_atm, thetav_sfc, wind
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & richardson_number
      !
    ! Optional ACC queue
    !
    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
    INTEGER :: acc_async_queue

    REAL(wp) :: zthetav_mid, w1, w2
    INTEGER :: jb, jls, js

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_moist_richardson'

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

    w1 = fsl
    w2 = 1._wp - fsl

!$OMP PARALLEL
    CALL init(richardson_number, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js, zthetav_mid) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        zthetav_mid = w1 * thetav_atm(js,jb) + w2 * thetav_sfc(js,jb)
        richardson_number(js,jb) = zf(js,jb) * grav * (thetav_atm(js,jb) - thetav_sfc(js,jb)) / (zthetav_mid * wind(js,jb)**2)
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_moist_richardson
  !
  !=================================================================
  !
! #ifndef _OPENACC
!   ELEMENTAL &
! #endif
  PURE SUBROUTINE sfc_exchange_coefficients(              &
    & dz,                                                 &
    & pqm1, thetam1, mwind, rough_m, theta_sfc, qsat_sfc, &
    & km, kh, interp_fac_2m, interp_fac_10m               &
    & )

    REAL(wp), INTENT(in) :: &
      dz,        & ! height to be used as a reference height in surface layer
      thetam1,   &
      pqm1 ,     &
      mwind,     &
      rough_m,   &
      theta_sfc, &
      qsat_sfc
      !
    REAL(wp), INTENT(out) :: &
      kh,            &
      km,            &
      interp_fac_2m, &
      interp_fac_10m

    !$ACC ROUTINE SEQ

    REAL(wp) :: z_mc, RIB, tcn_mom, tcn_heat, &
      &         shfl_local, lhfl_local, bflx1, ustar, obukhov_length, inv_bus_mom
    REAL(wp) :: zrat_2m, zrat_10m, zbm, zbh, zbn, zcbn, zcbu, zcbs, zmerge
    REAL(wp) :: tch
    REAL(wp) :: tcm

    INTEGER  :: itr

    ! to prevent floating-point arithmetic inconsistencies later in
    ! the interpolation to u 10m and 2m T/T_d: has been 0.01 before
    ! (Cray FP instead of IEEE 754 FP format)
    REAL(wp),parameter :: zepsec = 0.028_wp
    REAL(wp),parameter :: zcons17 = 1._wp / ckap**2._wp

    z_mc = dz
    !First guess for tch and tcm using bulk approach
    RIB = grav * (thetam1-theta_sfc) * (z_mc-rough_m) / (theta_sfc * mwind**2._wp)
    tcn_mom = (ckap/LOG(z_mc/rough_m))**2._wp
    tcm     = tcn_mom * stability_function_mom(RIB,z_mc/rough_m,tcn_mom)

    tcn_heat        = ckap**2._wp/(LOG(z_mc/rough_m)*LOG(z_mc/rough_m))
    tch             = tcn_heat * stability_function_heat(RIB,z_mc/rough_m,tcn_heat)

    !now iterate
    DO itr = 1, 5
      shfl_local = tch*mwind * (theta_sfc - thetam1)
      lhfl_local = tch*mwind * (qsat_sfc - pqm1)
      bflx1= shfl_local + vtmpc1 * theta_sfc * lhfl_local
      ustar= SQRT(tcm)*mwind

      obukhov_length = -ustar**3._wp * theta_sfc * rgrav / (ckap * bflx1)

      inv_bus_mom = 1._wp / businger_mom(rough_m,z_mc,obukhov_length)
      tch         = inv_bus_mom / businger_heat(rough_m,z_mc,obukhov_length)
      tcm         = inv_bus_mom * inv_bus_mom
    END DO


    kh  = tch
    km  = tcm

    ! Interpolation heights
    zrat_2m  =  2._wp / z_mc
    zrat_10m = 10._wp / z_mc

    ! Transfer coefficients for momentum, heat, and neutral conditions
    zbm  = 1._wp / MAX(zepsec, SQRT(km) / ckap)
    zbh  = 1._wp / MAX(zepsec, kh * zbm / ckap**2._wp)
    zbn  = 1._wp / MAX(zepsec, SQRT(tcn_mom) / ckap)

    ! Compute 2m interpolation factor
    zcbn   = LOG(1._wp + (EXP (zbn) - 1._wp) * zrat_2m )
    zcbs   = -(zbn - zbh) * zrat_2m
    zcbu   = -LOG(1._wp + (EXP (zbn - zbh) - 1._wp) * zrat_2m)
    zmerge = MERGE(zcbs, zcbu, RIB > 0._wp)

    interp_fac_2m   = (zcbn + zmerge) / zbh

    ! Compute 10m interpolation factor

    zcbn   = LOG(1._wp + (EXP (zbn) - 1._wp) * zrat_10m )
    zcbs   = -(zbn - zbh) * zrat_10m
    zcbu   = -LOG(1._wp + (EXP (zbn - zbh) - 1._wp) * zrat_10m)
    zmerge = MERGE(zcbs, zcbu, RIB > 0._wp)

    interp_fac_10m   = (zcbn + zmerge) / zbh

  END SUBROUTINE sfc_exchange_coefficients
  !
  !=================================================================
  !
  ! Initial release by Junhong Lee; MPI-M (2020-02)
  ! Rewritten for the new framework of ICON-Sapphire by Junhong Lee, Mira Shevchenko; MPI-M (2022-10)
  !
  SUBROUTINE compute_sfc_exchange_coefficients(          &
    & domain,                                            &
    & nvalid, indices, dz,                               &
    & pqm1,                                              &
    & thetam1, mwind, rough_m, theta_sfc, qsat_sfc,      &
    & km, kh, interp_fac_2m, interp_fac_10m,             &
    & opt_acc_async_queue                                &
    )

    ! Domain information
    TYPE(t_domain),        INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(vp), DIMENSION(:,:), INTENT(in) :: &
      dz
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      thetam1,   &
      pqm1 ,     &
      mwind,     &
      rough_m,   &
      theta_sfc, &
      qsat_sfc
      !
      ! Output variables
      !
      REAL(wp), DIMENSION(:,:), INTENT(out) :: &
        kh,            &
        km,            &
        interp_fac_2m, &
        interp_fac_10m
      !
      ! Optional ACC queue
      !
      INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue
      INTEGER :: acc_async_queue

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sfc_exchange_coefficients'

    INTEGER  :: jb, jls, js
    REAL(wp) :: dz_temp

    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

!$OMP PARALLEL
    CALL init(km,         lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(kh,         lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(interp_fac_2m,  lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
    CALL init(interp_fac_10m, lacc=.TRUE., opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb,jls,js, dz_temp) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_async_queue) PRIVATE(js, dz_temp)
      DO jls = 1, nvalid(jb)
        js=indices(jls,jb)

        ! Note: dz=p_nh_metrics%ddqz_z_half(:,nlevp1,:) from atmosphere (src/atm_dyn_iconam/mo_vertical_grid.f90)
        ! is twice of what we want
        dz_temp = REAL(dz(js,jb), KIND=wp) * 0.5_wp
        CALL sfc_exchange_coefficients(                                                      &
          & dz_temp,                                                                         &
          & pqm1(js,jb),                                                                     &
          & thetam1(js,jb), mwind(js,jb), rough_m(js,jb), theta_sfc(js,jb), qsat_sfc(js,jb), &
          & km(js,jb), kh(js,jb), interp_fac_2m(js,jb), interp_fac_10m(js,jb)                &
          & )
      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_sfc_exchange_coefficients  !
  !
  !=================================================================
  !
  !! first guess stability_function_mom
  !! Taken from Holstag & Boville 1992 J.Clim, Eqs.2.9-2.11
  !!------------------------------------------------------------------------
#ifndef _OPENACC
  ELEMENTAL &
#endif
  PURE FUNCTION stability_function_mom(RIB, hz0, tc) RESULT(stab_fun)
     REAL(wp), INTENT(IN) :: RIB, hz0, tc

     REAL(wp) :: stab_fun, hz0_fac

     !$ACC ROUTINE SEQ

     IF(RIB.GE.0._wp)THEN
       stab_fun = 1._wp / ( 1._wp + 10._wp*RIB*(1._wp+8._wp*RIB) )
    ELSE
       hz0_fac = ( max(hz0, 1._wp)**(1._wp/3._wp) - 1._wp )**1.5_wp ! FLO - hz0 can be < 1 then, the **1.5 is invalid.
       !for water surface (z0/h)**(1/3)<<1 giving hz0_fac=SQRT(h/z0)
       !Generally it is explicitly written for water surface but i don't
       !see any reason to do that.
       stab_fun = 1._wp + 10._wp*ABS(RIB)/(1._wp + 75._wp*tc*hz0_fac*SQRT(ABS(RIB)))
     END IF

  END FUNCTION stability_function_mom
  !
  ! first guess stability_function_heat
  !------------------------------------------------------------------------
#ifndef _OPENACC
  ELEMENTAL &
#endif
  PURE FUNCTION stability_function_heat(RIB, hzh, tc) RESULT(stab_fun)
     REAL(wp), INTENT(IN) :: RIB, hzh, tc

     REAL(wp) :: stab_fun, hzh_fac

     !$ACC ROUTINE SEQ

     IF(RIB.GE.0._wp)THEN
       stab_fun = 1._wp / ( 1._wp + 10._wp*RIB*(1._wp+8._wp*RIB) )
     ELSE
       hzh_fac = ( max(hzh, 1._wp)**(1._wp/3._wp) - 1._wp )**1.5_wp
       !for water surface (zh/h)**(1/3)<<1 giving hzh_fac=SQRT(h/zh)
       !Generally it is explicitly written for water surface but i don't
       !see any reason to do that.
       stab_fun = 1._wp + 15._wp*ABS(RIB)/(1._wp + 75._wp*tc*hzh_fac*SQRT(ABS(RIB)))
     END IF
  END FUNCTION stability_function_heat

  !
  ! factor_mom
  !------------------------------------------------------------------------
  ! Businger Dyer similarity profile:
  ! Louis (1979) A Parametirc model of vertical eddy fluxes in the atmosphere
  ! and R. B. Stull's book
  !------------------------------------------------------------------------
#ifndef _OPENACC
  ELEMENTAL &
#endif
  PURE FUNCTION businger_mom(z0, z1, L) RESULT(factor)
    ! Use VALUE attribute in order to prevent 'Reference argument passing prevents parallelization' message by nvhpc
    REAL(wp), VALUE, INTENT(IN) :: z0, z1, L
    REAL(wp) :: factor, zeta, psi, lamda
    REAL(wp) :: zeta0, psi0, lamda0

    !$ACC ROUTINE SEQ

    IF(L > 0._wp)THEN !Stable
      zeta  = z1/L
      zeta0 = z0/L
      IF(zeta > 1._wp)THEN !Zeng etal 1997 J. Clim
        psi    = -bsm + bsm*zeta0 + (1._wp-bsm)*LOG(zeta) -  zeta + 1._wp
        factor = (LOG(z1/z0) - psi) / ckap
      ELSE
        psi  = -bsm*zeta
        psi0 = -bsm*zeta0
        factor = ( LOG(z1/z0) - psi + psi0 ) / ckap
      END IF
    ELSEIF(L < 0._wp)THEN !unstable
      zeta   = z1/L
      zeta0  = z0/L
      lamda  = SQRT(SQRT(1._wp - bum*zeta))
      lamda0 = SQRT(SQRT(1._wp - bum*zeta0))

      psi    = 2._wp * LOG(1._wp+lamda) + LOG(1._wp+lamda*lamda) - &
              2._wp * ATAN(lamda) + pi_2 - 3._wp*ln2

      psi0   = 2._wp * LOG(1._wp+lamda0) + LOG(1._wp+lamda0*lamda0) - &
              2._wp * ATAN(lamda0) + pi_2 - 3._wp*ln2

      factor = ( LOG(z1/z0) - psi + psi0 ) / ckap
    ELSE !neutral
      factor = LOG(z1/z0) / ckap
    END IF

  END FUNCTION businger_mom

  !
  ! factor_heat
  !------------------------------------------------------------------------
  ! Businger Dyer similarity profile:
  !------------------------------------------------------------------------
#ifndef _OPENACC
  ELEMENTAL &
#endif
  PURE FUNCTION businger_heat(z0, z1, L) RESULT(factor)
    ! Use VALUE attribute in order to prevent 'Reference argument passing prevents parallelization' message by nvhpc
    REAL(wp), VALUE, INTENT(IN) :: z0, z1, L
    REAL(wp) :: factor, zeta, lamda, psi
    REAL(wp) :: zeta0, lamda0, psi0

    !$ACC ROUTINE SEQ

    IF(L > 0._wp)THEN !Stable
      zeta   = z1/L
      zeta0  = z0/L
      IF(zeta > 1._wp)THEN !Zeng etal 1997 J. Clim
        psi    = -bsh + bsh*zeta0 + (1._wp-bsh)*LOG(zeta) -  zeta + 1._wp
        factor = (LOG(z1/z0) - psi) / ckap
      ELSE
        psi    = -bsh*zeta
        psi0   = -bsh*zeta0
        factor = (LOG(z1/z0) - psi + psi0) / ckap
      END IF
    ELSEIF(L < 0._wp)THEN !unstable
      zeta   = z1/L
      zeta0  = z0/L
      lamda  = SQRT(1._wp - buh*zeta)
      lamda0 = SQRT(1._wp - buh*zeta0)
      psi    = 2._wp * ( LOG(1._wp+lamda) - ln2 )
      psi0   = 2._wp * ( LOG(1._wp+lamda0) - ln2 )
      factor = (LOG(z1/z0) - psi + psi0) / ckap
    ELSE !Neutral
      factor = LOG(z1/z0) / ckap
    END IF

  END FUNCTION businger_heat


END MODULE mo_vdf_diag_smag
