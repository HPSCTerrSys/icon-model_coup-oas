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

! Compile test to ensure all real and integer kinds provided in `mo_kind.f90`
! are available

PROGRAM test_kind

#ifdef __HAVE_QUAD_PRECISION
  USE mo_kind, ONLY: qp
#endif
  USE mo_kind, ONLY: wp, xwp, vp, dp, sp, &
      &                       i1, i2, i4, i8

  IMPLICIT NONE

  ! Real kinds
#ifdef __HAVE_QUAD_PRECISION
  REAL(qp) :: val_qp
#endif
  REAL(wp) :: val_wp
  REAL(xwp) :: val_xwp
  REAL(vp) :: val_vp
  REAL(dp) :: val_dp
  REAL(sp) :: val_sp

  ! Integer kinds
  INTEGER(i1) :: val_i1
  INTEGER(i2) :: val_i2
  INTEGER(i4) :: val_i4
  INTEGER(i8) :: val_i8

  ! Unused variables
#ifdef __HAVE_QUAD_PRECISION
  val_qp = 1._qp
#endif
  val_wp = 1._wp
  val_xwp = 1._xwp
  val_vp = 1._vp
  val_dp = 1._dp
  val_sp = 1._sp
  val_i1 = 1_i1
  val_i2 = 1_i2
  val_i4 = 1_i4
  val_i8 = 1_i8
END PROGRAM test_kind
