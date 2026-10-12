! Copyright (c) 2005-2019, SHTOOLS
! All rights reserved
!
! Contact: https://github.com/SHTOOLS/SHTOOLS
! Authors: https://github.com/SHTOOLS/SHTOOLS/blob/master/AUTHORS.md
! See LICENSES/BSD-3-Clause.txt for license information
! SPDX-License-Identifier: BSD-3-Clause
!
! This file has been modified for the use in ICON
! ---------------------------------------------------------------
!
! This module provides the orthonormal associated Legendre polynomials
! for spherical harmonics calculations in ICON. The subroutines are
! not optimized and it can be beneficial to use machine specific
! math libraries instead.
!
! The following minor changes have been made to comply with
! ICON coding standards:
! - The SAVE attribute has been removed, f1, f2 and sqr are
!   calculated all the time
! - Maybe the SAVE in combination with OMP threadprivate, which
!   was in the original SHTOOLS code, would also be an option
!   for ICON? But usually ICON does not use this approach.
! - dp has been renamed to wp and is imported from mo_kind
! - stop has been replaced by call finish
! - print statements have been removed
!
! DWD, 2024
!

MODULE mo_math_legendre

  USE mo_kind,               ONLY: wp
  USE mo_math_constants,     ONLY: pi
  USE mo_exception,          ONLY: finish
  USE iso_fortran_env,       ONLY: int32

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: PlmON, PlmIndex

  CHARACTER(len=*), PARAMETER :: routine = 'mo_math_legendre'
  REAL(wp),         PARAMETER :: sqr4pi = SQRT(4.0_wp*pi)

CONTAINS

  SUBROUTINE PlmON(p, lmax, z, csphase, cnorm, errstatus)
    !------------------------------------------------------------------------------
    !
    !   This function evalutates all of the normalized associated Legendre
    !   functions up to degree lmax. The functions are initially scaled by
    !   10^280 sin^m in order to minimize the effects of underflow at large m
    !   near the poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299).
    !   On a macOS system with a maximum allowable double precision value of
    !   2.225073858507203E-308 the scaled portion of the algorithm will not overflow
    !   for degrees less than or equal to 2800.
    !
    !   For each value of m, the rescaling factor is computed as
    !   rescalem=rescalem*sin(theta), with the intial value of rescalem being equal
    !   to 1/scalef (which is here equal to 10^280). This will gradually reduce this
    !   huge number to a tiny number, and will ultimately underflow. In order to
    !   prevent this underflow, when rescalem becomes less than 10^(-280), the
    !   subsequent rescaling factors of sin(theta) will be directly applied to Plm,
    !   and then this number will be multipled by the old value of rescalem.
    !
    !   Temporary variables in saved in an allocated array. In order to explicitly
    !   deallocate this memory, call this routine with a spherical harmonic degree
    !   of -1.
    !
    !   Calling Parameters
    !
    !       OUT
    !           p           A vector of all associated Legendgre polynomials
    !                       evaluated at z up to lmax. The length must by greater
    !                       or equal to (lmax+1)*(lmax+2)/2.
    !       OPTIONAL (IN)
    !           csphase     1: Do not include the phase factor of (-1)^m (default).
    !                       -1: Apply the phase factor of (-1)^m.
    !           cnorm       0: Use real normalization.
    !                       1: Use complex normalization.
    !       IN
    !           lmax        Maximum spherical harmonic degree to compute.
    !           z           cos(colatitude) or sin(latitude).
    !
    !       OPTIONAL (OUT)
    !           errstatus  If present, instead of executing a STOP when an error
    !                       is encountered, the variable errstatus will be
    !                       returned describing the error.
    !                       0 = No errors;
    !                       1 = Improper dimensions of input array;
    !                       2 = Improper bounds for input variable;
    !                       3 = Error allocating memory;
    !                       4 = File IO error.
    !   Notes:
    !   1.  The employed normalization is the "orthonormalized convention." The
    !       integral of (PlmON*cos(m theta))**2 or (PlmON*sin (m theta))**2 over
    !       all space is 1.
    !   2.  The integral of PlmON**2 over (-1,1) is (2 - delta(0,m))/2pi.
    !       If CNORM=1, then this is equal to 1/2pi.
    !   3.  The index of the array p corresponds to l*(l+1)/2 + m + 1. As such
    !       the array p should be dimensioned as (lmax+1)*(lmax+2)/2 in the
    !       calling routine.
    !   4.  The default is to exclude the Condon-Shortley phase of (-1)^m.
    !
    !   Copyright (c) 2005-2019, SHTOOLS
    !   All rights reserved.
    !
    !------------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER(int32), INTENT(in) :: lmax
    REAL(wp), INTENT(out) :: p(:)
    REAL(wp), INTENT(in) :: z
    INTEGER(int32), INTENT(in),  OPTIONAL :: csphase, cnorm
    INTEGER(int32), INTENT(out), OPTIONAL :: errstatus
    REAL(wp) :: pm2, pm1, pmm, plm, rescalem, u, scalef
    INTEGER(int32) :: k, kstart, m, l
    INTEGER(int32) :: phase

    ! In the original SHTOOLS code those work array where ALLOCATABLE
    ! with SAVE attribute. Then they can store precomputed results which
    ! makes the algorithm more efficient. They had an OMP threadprivate
    ! for OpenMP. Such a combination of ALLOCATE+SAVE+threadprivate is
    ! not used anywhere else in ICON, though. Therefore here this less
    ! efficient approach. For performance critical application the SHTOOLS
    ! code should maybe be replace by a high-performance math library.
    REAL(wp) :: sqr(2*lmax+1)
    REAL(wp) :: f1((lmax+1)*(lmax+2)/2)
    REAL(wp) :: f2((lmax+1)*(lmax+2)/2)

    IF (PRESENT(errstatus)) errstatus = 0

    IF (SIZE(p) < (lmax+1)*(lmax+2)/2) THEN
      IF (PRESENT(errstatus)) THEN
        errstatus = 1
        RETURN
      ELSE
        CALL finish(routine,'P must be dimensioned as (LMAX+1)*(LMAX+2)/2')
      END IF
    ELSE IF (lmax < 0) THEN
      IF (PRESENT(errstatus)) THEN
        errstatus = 2
        RETURN
      ELSE
        CALL finish(routine,"LMAX must be greater than or equal to 0.")
      END IF
    ELSE IF (ABS(z) > 1.0_wp) THEN
      IF (PRESENT(errstatus)) THEN
        errstatus = 2
        RETURN
      ELSE
        CALL finish(routine,"ABS(Z) must be less than or equal to 1.")
      END IF
    END IF
    IF (PRESENT(csphase)) THEN
      IF (csphase == -1) THEN
        phase = -1
      ELSE IF (csphase == 1) THEN
        phase = 1
      ELSE
        IF (PRESENT(errstatus)) THEN
          errstatus = 2
          RETURN
        ELSE
          CALL finish(routine,"CSPHASE must be 1 (exclude) or -1 (include).")
        END IF
      END IF
    ELSE
      phase = 1
    END IF

    scalef = 1.0e-280_wp

    !----------------------------------------------------------------------
    !
    !   Precompute square roots of integers that are used several times.
    !
    !----------------------------------------------------------------------
    DO l = 1, 2 * lmax+1
      sqr(l) = SQRT(DBLE(l))
    END DO

    !----------------------------------------------------------------------
    !
    !   Precompute multiplicative factors used in recursion relationships
    !       Plmbar(l,m) = x*f1(l,m)*Plmbar(l-1,m) - Plmbar(l-2,m)*f2(l,m)
    !       k = l*(l+1)/2 + m + 1
    !   Note that prefactors are not used for the case when m=l and m=l-1,
    !   as a different recursion is used for these two values.
    !
    !----------------------------------------------------------------------
    k = 3

    DO l = 2, lmax, 1
      k = k + 1
      f1(k) = sqr(2*l-1) * sqr(2*l+1) / DBLE(l)
      f2(k) = DBLE(l-1) * sqr(2*l+1) / sqr(2*l-3) / DBLE(l)
      DO m = 1, l-2
        k = k + 1
        f1(k) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
        f2(k) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
             / sqr(2*l-3) / sqr(l+m) / sqr(l-m)
      END DO
      k = k + 2
    END DO

    !--------------------------------------------------------------------------
    !
    !   Calculate P(l,0). These are not scaled.
    !
    !--------------------------------------------------------------------------

    u = SQRT((1.0_wp - z) * (1.0_wp + z)) ! sin(theta)
    pm2 = 1.0_wp / sqr4pi
    p(1) = pm2
    IF (lmax == 0) RETURN
    pm1 = sqr(3) * z / sqr4pi
    p(2) = pm1

    k = 2
    DO l = 2, lmax, 1
      k = k + l
      plm = f1(k) * z * pm1 - f2(k) * pm2
      p(k) = plm
      pm2 = pm1
      pm1 = plm
    END DO

    !--------------------------------------------------------------------------
    !
    !   Calculate P(m,m), P(m+1,m), and P(l,m)
    !
    !--------------------------------------------------------------------------
    IF (PRESENT(cnorm)) THEN
      IF (cnorm == 1) THEN
        pmm = scalef / sqr4pi
      ELSE
        pmm = sqr(2)*scalef / sqr4pi
      END IF
    ELSE
      pmm = sqr(2)*scalef / sqr4pi
    END IF
    rescalem = 1.0_wp / scalef
    kstart = 1

    DO m = 1, lmax - 1, 1
      rescalem = rescalem * u

      ! Calculate P(m,m)
      kstart = kstart + m + 1
      pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
      p(kstart) = pmm * rescalem
      pm2 = pmm

      ! Calculate P(m+1,m)
      k = kstart + m + 1
      pm1 = z * sqr(2 * m + 3) * pmm
      p(k) = pm1 * rescalem

      ! Calculate P(l,m)
      DO l = m + 2, lmax, 1
        k = k + l
        plm = z * f1(k) * pm1 - f2(k) * pm2
        p(k) = plm * rescalem
        pm2 = pm1
        pm1 = plm
      END DO
    END DO

    ! Calculate P(lmax,lmax)
    rescalem = rescalem * u
    kstart = kstart + m + 1
    pmm = phase * pmm * sqr(2*lmax+1) / sqr(2*lmax)
    p(kstart) = pmm * rescalem

  END SUBROUTINE PlmON

  FUNCTION PlmIndex(l, m)
    !-------------------------------------------------------------------------------
    !
    !   This function will return the index corresponding
    !   to a given l and m in the arrays of Legendre Polynomials
    !   generated by routines such as PlmBar and PlmSchmidt.
    !
    !   Calling Parameters
    !
    !       l   Spherical harmonic angular degree.
    !       m   Spherical harmonic angular order.
    !
    !   Copyright (c) 2005-2019, SHTOOLS
    !   All rights reserved.
    !
    !-------------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER(int32) :: PlmIndex
    INTEGER(int32), INTENT(in) :: l, m

    IF (l < 0) THEN
      CALL finish(routine,"L must be greater of equal to 0.")
    ELSE IF (m < 0 .OR.  m > l) THEN
      CALL finish("M must be greater than or equal to zero and less than or equal to L.")
    END IF
    PlmIndex = (l*(l+1))/2+m+1

  END FUNCTION PlmIndex

END MODULE mo_math_legendre
