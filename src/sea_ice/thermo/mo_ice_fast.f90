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

! This module contains "fast" processes in sea ice thermodynamics that are
! used by the atmosphere (uncoupled or coupled) or the ocean (uncoupled).

MODULE mo_ice_fast
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2007
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ltimer
  USE mo_exception,           ONLY: finish, message

  USE mo_timer,               ONLY: timer_start, timer_stop, timer_ice_fast

  USE mo_physical_constants,  ONLY: rhoi, ki, ci
  USE mo_sea_ice_nml,         ONLY: i_ice_therm, hci_layer, Tf

  USE mo_ice_winton,          ONLY: set_ice_temp_winton
  USE mo_ice_zerolayer,       ONLY: set_ice_temp_zerolayer, set_ice_temp_zerolayer_analytical
  USE mo_ice_parameterizations, ONLY: set_ice_albedo

  USE mo_fortran_tools,       ONLY: set_acc_host_or_device, set_acc_async_queue

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ice_fast

  CHARACTER(len=12)           :: str_module    = 'mo_ice_fast'  ! Output of module for 1 line debug

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! ! ice_fast: ice thermodynamics at atmospheric time-step.
  !!   Calculates ice/snow surface temp, air-ice fluxes and sets albedos.
  !!
  !! This function changes:
  !! Tsurf and T1, T2 (winton)  - temperature of snow/ice layer(s)
  !! Qtop, Qbot                 - heat flux available for surface/bottom melting
  !! alb{vis/nir}{dir/dif}      - albedos
  !!
  SUBROUTINE ice_fast(i_startidx_c, i_endidx_c, nbdim, kice, pdtime, &
            &   Tsurf,          & ! Surface temperature [degC]
            &   T1,             & ! Temperature of upper layer [degC]
            &   T2,             & ! Temperature of lower layer [degC]
            &   hi,             & ! Ice thickness
            &   hs,             & ! Snow thickness
            &   Qtop,           & ! Energy flux available for surface melting [W/m2]
            &   Qbot,           & ! Energy flux available for bottom melting [W/m2]
            &   SWnet,          & ! Net shortwave flux [W/m^2]
            &   nonsolar,       & ! Latent and sensible heat flux and longwave radiation [W/m^2]
            &   dnonsolardT,    & ! Derivative of non-solar fluxes w.r.t. temperature [W/m^2/K]
            &   Tfw,            & ! Freezing temperature of the ocean
            &   albvisdir,      & ! Albedo VIS, direct/parallel
            &   albvisdif,      & ! Albedo VIS, diffuse
            &   albnirdir,      & ! Albedo NIR, direct/parallel
            &   albnirdif,      & ! Albedo NIR, diffuse
            &   doy,            & ! Day of the year
            &   lacc,           &
            &   opt_acc_async_queue)

    INTEGER, INTENT(IN)    :: i_startidx_c, i_endidx_c, nbdim, kice
    REAL(wp),INTENT(IN)    :: pdtime
    REAL(wp),INTENT(INOUT) :: Tsurf      (nbdim,kice)
    REAL(wp),INTENT(INOUT) :: T1         (nbdim,kice)
    REAL(wp),INTENT(INOUT) :: T2         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hi         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hs         (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qtop       (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qbot       (nbdim,kice)
    REAL(wp),INTENT(IN)    :: SWnet      (nbdim,kice)
    REAL(wp),INTENT(IN)    :: nonsolar   (nbdim,kice)
    REAL(wp),INTENT(IN)    :: dnonsolardT(nbdim,kice)
    REAL(wp),INTENT(IN)    :: Tfw        (nbdim)
    REAL(wp),INTENT(OUT)   :: albvisdir  (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: albvisdif  (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: albnirdir  (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: albnirdif  (nbdim,kice)

    INTEGER, OPTIONAL,INTENT(IN)  :: doy
    LOGICAL, OPTIONAL,INTENT(IN)  :: lacc
    INTEGER, INTENT(IN), OPTIONAL :: opt_acc_async_queue

    INTEGER :: jk, ji
    LOGICAL :: lzacc
    INTEGER :: acc_async_queue

    !-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)
    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

    IF (ltimer) CALL timer_start(timer_ice_fast)

    SELECT CASE (i_ice_therm)

    CASE (1)
      CALL set_ice_temp_zerolayer(i_startidx_c, i_endidx_c, nbdim, kice, pdtime, &
                            &   Tsurf, hi, hs, Qtop, Qbot, SWnet, nonsolar, dnonsolardT, Tfw, &
                            &   lacc=lzacc, opt_acc_async_queue=acc_async_queue)

    CASE (2)
      CALL set_ice_temp_winton(i_startidx_c, i_endidx_c, nbdim, kice, pdtime, &
                    &   Tsurf, T1, T2, hi, hs, Qtop, Qbot, SWnet, nonsolar, dnonsolardT, Tfw, &
                    &   lacc=lzacc, opt_acc_async_queue=acc_async_queue)

    CASE (3)
      IF ( .NOT. PRESENT(doy) ) THEN
        CALL finish(TRIM('mo_ice_interface:ice_fast'),'i_ice_therm = 3 not allowed in this context')
      ENDIF
      CALL set_ice_temp_zerolayer_analytical(i_startidx_c, i_endidx_c, nbdim, kice, &
            &   Tsurf, hi, hs, Qtop, Qbot, Tfw, doy, &
            &   lacc=lzacc, opt_acc_async_queue=acc_async_queue)

    CASE (4)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO ji = 1, kice
        DO jk = 1, nbdim
          IF ( hi(jk,ji) > 0._wp ) THEN
            Tsurf(jk,ji) = min(0._wp, Tsurf(jk,ji) + (SWnet(jk,ji)+nonsolar(jk,ji) + ki/hi(jk,ji)*(Tf-Tsurf(jk,ji))) &
        &               / (ci*rhoi*hci_layer/pdtime-dnonsolardT(jk,ji)+ki/hi(jk,ji)))
          ELSE
            Tsurf(jk,ji) = Tf
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

    END SELECT

    ! New albedo based on the new surface temperature
    CALL set_ice_albedo(i_startidx_c, i_endidx_c, nbdim, kice, Tsurf, hi, hs, &
      & albvisdir, albvisdif, albnirdir, albnirdif, lacc=lzacc, opt_acc_async_queue=acc_async_queue)

    IF (ltimer) CALL timer_stop(timer_ice_fast)

   END SUBROUTINE ice_fast

END MODULE mo_ice_fast
