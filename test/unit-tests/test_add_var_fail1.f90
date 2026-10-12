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

PROGRAM test_add_var_fail1

  USE mo_kind, ONLY: dp, sp
  USE mo_impl_constants, ONLY: REAL_T, SINGLE_T
  USE mo_var, ONLY: t_var
  USE mo_var_list_register, ONLY: vlr_add
  USE mo_var_list, ONLY: add_var, t_var_list_ptr
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_test_common, ONLY: test_pass

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: listname = "foo"
  CHARACTER(*), PARAMETER :: dname = "bar"
  TYPE(t_var), POINTER :: vl_elem => NULL()
  TYPE(t_var_metadata) :: info
  TYPE(t_var_list_ptr) :: var_list_ptr
  REAL(dp) :: rval
  REAL(sp) :: sval

  ! Register var_list_ptr
  CALL vlr_add(var_list_ptr, TRIM(listname), model_type="test")

  ! FAIL 1: test_add_var_fail1 - data_type=SINGLE_T, and passing rvals
  CALL add_var(SINGLE_T, var_list_ptr, dname//"2", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, info%used_dimensions(1:info%ndims), vl_elem, &
      & initval_r=info%initval%rval, resetval_r=info%resetval%rval, &
      & missval_r=info%missval%rval)

  CALL test_pass()

END PROGRAM test_add_var_fail1
