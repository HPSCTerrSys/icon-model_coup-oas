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

! Contains routines for computing the 1D JONSWAP spectrum,
! and potentially other helper routines for physics computations.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_phy_util

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_math_constants,       ONLY: pi2
  USE mo_physical_constants,   ONLY: grav

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: jonswap_nonblk
  PUBLIC :: fetch_law_nonblk

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_phy_util'

CONTAINS

  !>
  !! Calculation of the JONSWAP spectrum according to
  !! Hasselmann et al. 1973. Adaptation of WAM 4.5
  !! subroutine JONSWAP.
  !!
  !! Non-blocked version
  !!
  SUBROUTINE jonswap_nonblk(i_startidx, i_endidx, freqs, gamma, sa, sb, flmin, alphaj, fp, et)
    INTEGER,       INTENT(IN)    :: i_startidx, i_endidx
    REAL(wp),      INTENT(IN)    :: freqs(:)      !! FREQUENCIES.
    REAL(wp),      INTENT(IN)    :: gamma         !! OVERSHOOT FACTOR.
    REAL(wp),      INTENT(IN)    :: sa            !! LEFT PEAK WIDTH.
    REAL(wp),      INTENT(IN)    :: sb            !! RIGHT PEAK WIDTH.
    REAL(wp),      INTENT(IN)    :: flmin
    REAL(wp),      INTENT(IN)    :: alphaj(:)     !! OVERALL ENERGY LEVEL OF JONSWAP SPECTRA.
    REAL(wp),      INTENT(IN)    :: fp(:)         !! PEAK FREQUENCIES.
    REAL(wp),      INTENT(INOUT) :: et(:,:)       !! JONSWAP SPECTRA (is only OUT parameter)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':JONSWAP_nonblk'

    REAL(wp) :: ARG, sigma, G2ZPI4FRH5M

    INTEGER :: jc,jf

    DO jf = 1,SIZE(freqs)

      G2ZPI4FRH5M = grav**2 / pi2**4 * freqs(jf)**(-5)

      DO jc = i_startidx, i_endidx

        sigma = MERGE(sb,sa, freqs(jf)>fp(jc))
        ET(jc,jf) = 0._wp

        ARG = 1.25_wp*(FP(jc)/freqs(jf))**4
        IF (ARG.LT.50.0_wp) THEN
          ET(jc,jf) = ALPHAJ(jc) * G2ZPI4FRH5M * EXP(-ARG)
        END IF

        ARG = 0.5_wp*((freqs(jf)-FP(jc)) / (sigma*FP(jc)))**2
        IF (ARG.LT.99._wp) THEN
          ET(jc,jf) = ET(jc,jf)*exp(log(GAMMA)*EXP(-ARG))
        END IF

        ET(jc,jf) = MAX(ET(jc,jf),flmin)

      END DO  !jc
    END DO  !jf

  END SUBROUTINE jonswap_nonblk


  !>
  !! Calculation of JONSWAP parameters.
  !!
  !! Calculate the peak frequency from a fetch law
  !! and the JONSWAP alpha.
  !!
  !! Developted by S. Hasselmann (July 1990) and H. Guenther (December 1990).
  !! K.HASSELMAN,D.B.ROOS,P.MUELLER AND W.SWELL. A parametric wave prediction
  !! model. Journal of physical oceanography, Vol. 6, No. 2, March 1976.
  !!
  !! Adopted from WAM 4.5.
  !!
  !! non-blocked version
  !!
  SUBROUTINE fetch_law_nonblk(i_startidx, i_endidx, fetch, fpmax, sp10m, fp, alphaj)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':fetch_law_nonblk'
    !
    INTEGER,       INTENT(IN)    :: i_startidx, i_endidx
    REAL(wp),      INTENT(IN)    :: fetch
    REAL(wp),      INTENT(IN)    :: fpmax     ! maximum peak frequency (Hz)
    REAL(wp),      INTENT(IN)    :: sp10m(:)  ! wind speed at 10m (m/s)
    REAL(wp),      INTENT(INOUT) :: fp(:)     ! jonswap peak frequency (1/s) (is only OUT parameter)
    REAL(wp),      INTENT(INOUT) :: alphaj(:) ! jonswap alpha (-) (is only OUT parameter)

    INTEGER :: jc

    REAL(wp), PARAMETER :: A = 2.84_wp,  D = -(3._wp/10._wp) !! PEAK FREQUENCY FETCH LAW CONSTANTS
    REAL(wp), PARAMETER :: B = 0.033_wp, E = 2._wp/3._wp     !! ALPHA-PEAK FREQUENCY LAW CONSTANTS

    REAL(wp) :: ug

    ! ---------------------------------------------------------------------------- !
    !                                                                              !
    !     1. COMPUTE VALUES FROM FETCH LAWS.                                       !
    !        -------------------------------
    DO jc = i_startidx, i_endidx
      IF (sp10m(jc) > 0.1E-08_wp) THEN
        ug = grav / sp10m(jc)
        fp(jc) = MAX(0.13_wp, A*((grav*fetch)/(sp10m(jc)**2))**D)

        fp(jc) = MIN(fp(jc), fpmax/ug)

        alphaj(jc) = MAX(0.0081_wp, B * fp(jc)**E)
        fp(jc) = fp(jc) * ug
      ELSE
        alphaj(jc) = 0.0081_wp
        fp(jc) = fpmax
      END IF
    END DO

  END SUBROUTINE fetch_law_nonblk

END MODULE mo_wave_phy_util
