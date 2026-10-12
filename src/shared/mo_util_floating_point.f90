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

MODULE mo_util_floating_point
  USE :: mo_kind, ONLY: wp
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: is_fp_equivalent

CONTAINS

  ! Check if two floats are equivalent to a given tolerance (default machine epsilon)
  FUNCTION is_fp_equivalent(a,b, opt_tol)
    REAL(wp), INTENT(IN) :: a, b
    REAL(wp), OPTIONAL, INTENT(IN) :: opt_tol
    LOGICAL :: is_fp_equivalent
    REAL(wp) :: tol

    tol = EPSILON(1._wp)*MAX(ABS(a),ABS(b))
    IF (PRESENT(opt_tol)) tol = opt_tol

    IF ( ABS(a-b) >= tol ) THEN
      is_fp_equivalent = .FALSE.
    ELSE
      is_fp_equivalent = .TRUE.
    ENDIF
  END FUNCTION is_fp_equivalent

END MODULE mo_util_floating_point
