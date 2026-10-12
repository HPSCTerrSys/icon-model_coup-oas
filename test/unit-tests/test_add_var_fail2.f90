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

PROGRAM test_add_var_fail2

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

  ! FAIL 2: test_add_var_fail2 - data_type=REAL_T, and passing svals
  CALL add_var(REAL_T, var_list_ptr, dname//"2", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, info%used_dimensions(1:info%ndims), vl_elem, &
      & initval_s=info%initval%sval, resetval_s=info%resetval%sval, &
      & missval_s=info%missval%sval)

  CALL test_pass()

END PROGRAM test_add_var_fail2
