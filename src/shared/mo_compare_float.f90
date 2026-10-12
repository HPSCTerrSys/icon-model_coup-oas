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

! This module contains routines to compare floating point values for some given
! tolerance.

MODULE mo_compare_float

  USE mo_kind, ONLY: dp, sp

  IMPLICIT NONE

  PUBLIC :: notEqual

  PRIVATE

  ! Scaling tolerance based on amplitude of the two terms. For all contiguous, real-values `a,b` we have an upper
  ! bound on the floating point tolerance when taking a-b. Assuming a,b not very small we can represent this as:
  !   floating_point_tolerance(a-b) <= machine_epsilon * max(|a|,|b|)
  INTERFACE notEqual
    MODULE PROCEDURE notEqual_1d_dp
    MODULE PROCEDURE notEqual_1d_sp
    MODULE PROCEDURE notEqual_0d_dp
    MODULE PROCEDURE notEqual_0d_sp
  END INTERFACE notEqual

CONTAINS

  ! Elementwise comparison with given tolerances
  LOGICAL FUNCTION notEqual_1d_dp(a, b, rel_tol, abs_tol) RESULT(notEqual)
    REAL(dp), INTENT(IN) :: a(:)
    REAL(dp), INTENT(IN) :: b(:)
    REAL(dp), OPTIONAL, INTENT(IN) :: rel_tol ! relative tolerance, default EPSILON(a)
    REAL(dp), OPTIONAL, INTENT(IN) :: abs_tol ! absolute tolerance, default TINY(a)
    ! local
    INTEGER :: i

    notEqual = .FALSE. ! start with equal, unless proven otherwise
    DO i = 1,SIZE(a)
      IF (notEqual_0d_dp(a(i), b(i), rel_tol, abs_tol)) THEN
        notEqual = .TRUE.
        RETURN
      END IF
    END DO
  END FUNCTION notEqual_1d_dp


  LOGICAL FUNCTION notEqual_1d_sp(a, b, rel_tol, abs_tol) RESULT(notEqual)
    REAL(sp), INTENT(IN) :: a(:)
    REAL(sp), INTENT(IN) :: b(:)
    REAL(sp), OPTIONAL, INTENT(IN) :: rel_tol ! relative tolerance, default EPSILON(a)
    REAL(sp), OPTIONAL, INTENT(IN) :: abs_tol ! absolute tolerance, default TINY(a)
    ! local
    INTEGER :: i

    notEqual = .FALSE. ! start with equal, unless proven otherwise
    DO i = 1,SIZE(a)
      IF (notEqual_0d_sp(a(i), b(i), rel_tol, abs_tol)) THEN
        notEqual = .TRUE.
        RETURN
      END IF
    END DO
  END FUNCTION notEqual_1d_sp


  LOGICAL FUNCTION notEqual_0d_dp(a, b, rel_tol, abs_tol) RESULT(notEqual)
    REAL(dp), INTENT(IN) :: a
    REAL(dp), INTENT(IN) :: b
    REAL(dp), OPTIONAL, INTENT(IN) :: rel_tol ! relative tolerance, default EPSILON(a)
    REAL(dp), OPTIONAL, INTENT(IN) :: abs_tol ! absolute tolerance, default TINY(a)
    ! local
    REAL(dp) :: amplitude, lrel_tol, labs_tol

    lrel_tol = EPSILON(a) ! Default
    IF (PRESENT(rel_tol)) lrel_tol = rel_tol

    labs_tol = TINY(a) ! Default
    IF (PRESENT(abs_tol)) labs_tol = abs_tol

    amplitude = MAX(ABS(a), ABS(b))
    notEqual = ABS(a - b) > MAX(labs_tol, lrel_tol * amplitude)
  END FUNCTION notEqual_0d_dp


  LOGICAL FUNCTION notEqual_0d_sp(a, b, rel_tol, abs_tol) RESULT(notEqual)
    REAL(sp), INTENT(IN) :: a
    REAL(sp), INTENT(IN) :: b
    REAL(sp), OPTIONAL, INTENT(IN) :: rel_tol ! relative tolerance, default EPSILON(a)
    REAL(sp), OPTIONAL, INTENT(IN) :: abs_tol ! absolute tolerance, default TINY(a)
    ! local
    REAL(sp) :: amplitude, lrel_tol, labs_tol

    lrel_tol = EPSILON(a) ! Default
    IF (PRESENT(rel_tol)) lrel_tol = rel_tol

    labs_tol = TINY(a) ! Default
    IF (PRESENT(abs_tol)) labs_tol = abs_tol

    amplitude = MAX(ABS(a), ABS(b))
    notEqual = ABS(a - b) > MAX(labs_tol, lrel_tol * amplitude)
  END FUNCTION notEqual_0d_sp

END MODULE mo_compare_float
