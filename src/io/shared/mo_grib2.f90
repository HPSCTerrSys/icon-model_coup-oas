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

MODULE mo_grib2

  USE mo_kind,                  ONLY: dp

  IMPLICIT NONE

  PRIVATE

  ! max. number of additional GRIB2 integer keys per variable
  INTEGER, PARAMETER :: MAX_INT_KEYS = 15

  ! max. number of additional GRIB2 double keys per variable
  INTEGER, PARAMETER :: MAX_DBL_KEYS = 5

  ! max. array length for GRIB2 array keys
  INTEGER, PARAMETER :: MAX_ARR_LEN = 20


  TYPE t_grib2_global
    INTEGER :: centre
    INTEGER :: subcentre
    INTEGER :: generating_process
  END TYPE t_grib2_global

  TYPE t_grib2_int_key
    CHARACTER(len=50) :: key
    INTEGER           :: val
  END TYPE t_grib2_int_key

  TYPE t_grib2_dbl_key
    CHARACTER(len=50) :: key
    REAL(dp)          :: val
  END TYPE t_grib2_dbl_key

  TYPE t_grib2_intarr_key
    CHARACTER(len=50) :: key
    INTEGER           :: arr_len          ! actual array length
    INTEGER           :: vals(MAX_ARR_LEN)
  END TYPE t_grib2_intarr_key

  TYPE t_grib2_dblarr_key
    CHARACTER(len=50) :: key
    INTEGER           :: arr_len          ! actual array length
    REAL(dp)          :: vals(MAX_ARR_LEN)
  END TYPE t_grib2_dblarr_key

  TYPE t_grib2_key_list
    INTEGER                   :: nint_keys     ! no. of integer keys
    INTEGER                   :: ndbl_keys     ! no. of double keys
    INTEGER                   :: nintarr_keys  ! no. of integer array keys
    INTEGER                   :: ndblarr_keys  ! no. of double array keys
    TYPE (t_grib2_int_key)    :: int_key(MAX_INT_KEYS)
    TYPE (t_grib2_dbl_key)    :: dbl_key(MAX_DBL_KEYS)
    TYPE (t_grib2_intarr_key) :: intarr_key(MAX_INT_KEYS)
    TYPE (t_grib2_dblarr_key) :: dblarr_key(MAX_DBL_KEYS)
  END TYPE t_grib2_key_list

  TYPE t_grib2_var
    INTEGER :: discipline
    INTEGER :: category
    INTEGER :: number
    INTEGER :: bits
    INTEGER :: gridtype
    INTEGER :: subgridtype

    ! list of additional GRIB2 key/value pairs
    TYPE (t_grib2_key_list) :: additional_keys
  END TYPE t_grib2_var


  INTERFACE OPERATOR(+)
    MODULE PROCEDURE grib2_key_list_plus_int
    MODULE PROCEDURE grib2_key_list_plus_intarr
    MODULE PROCEDURE grib2_key_list_plus_dbl
    MODULE PROCEDURE grib2_key_list_plus_dblarr
  END INTERFACE OPERATOR(+)

  INTERFACE t_grib2_intarr_key
    MODULE PROCEDURE create__t_grib2_intarr_key
  END INTERFACE

  INTERFACE t_grib2_dblarr_key
    MODULE PROCEDURE create__t_grib2_dblarr_key
  END INTERFACE


  PUBLIC :: t_grib2_global
  PUBLIC :: t_grib2_var
  PUBLIC :: t_grib2_int_key, t_grib2_intarr_key
  PUBLIC :: t_grib2_dbl_key, t_grib2_dblarr_key
  PUBLIC :: t_grib2_key_list
  PUBLIC :: OPERATOR(+)

  ! constructor
  PUBLIC :: grib2_var

CONTAINS

  ! custom constructor for variable of type t_grib2_intarr_key
  !
  FUNCTION create__t_grib2_intarr_key(key, vals) RESULT(res)
    CHARACTER(len=*), INTENT(IN) :: key
    INTEGER         , INTENT(IN) :: vals(:)
    TYPE(t_grib2_intarr_key) :: res

    res%key                 = TRIM(key)
    res%arr_len             = SIZE(vals)
    res%vals(:)             = 0
    res%vals(1:res%arr_len) = vals
  END FUNCTION create__t_grib2_intarr_key

  ! custom constructor for variable of type t_grib2_dblarr_key
  !
  FUNCTION create__t_grib2_dblarr_key(key, vals) RESULT(res)
    CHARACTER(len=*), INTENT(IN) :: key
    REAL(dp)        , INTENT(IN) :: vals(:)
    TYPE(t_grib2_dblarr_key) :: res

    res%key                 = TRIM(key)
    res%arr_len             = SIZE(vals)
    res%vals(:)             = 0._dp
    res%vals(1:res%arr_len) = vals
  END FUNCTION create__t_grib2_dblarr_key

  ! constructor for GRIB2 derived data type
  !
  FUNCTION grib2_var(discipline, category, number, bits, gridtype, subgridtype)
    INTEGER, INTENT(IN) :: discipline
    INTEGER, INTENT(IN) :: category
    INTEGER, INTENT(IN) :: number
    INTEGER, INTENT(IN) :: bits
    INTEGER, INTENT(IN) :: gridtype
    INTEGER, INTENT(IN) :: subgridtype
    TYPE(t_grib2_var) :: grib2_var

    grib2_var%discipline  = discipline
    grib2_var%category    = category
    grib2_var%number      = number
    grib2_var%bits        = bits
    grib2_var%gridtype    = gridtype
    grib2_var%subgridtype = subgridtype

    grib2_var%additional_keys%nint_keys    = 0
    grib2_var%additional_keys%ndbl_keys    = 0
    grib2_var%additional_keys%nintarr_keys = 0
    grib2_var%additional_keys%ndblarr_keys = 0
  END FUNCTION grib2_var

  FUNCTION grib2_key_list_plus_int(a, b)
    TYPE(t_grib2_var) :: grib2_key_list_plus_int
    TYPE(t_grib2_var),     INTENT(IN) :: a
    TYPE(t_grib2_int_key), INTENT(IN) :: b
    ! local variables
    INTEGER :: i

    grib2_key_list_plus_int = a
    IF (a%additional_keys%nint_keys < MAX_INT_KEYS) THEN
      i = grib2_key_list_plus_int%additional_keys%nint_keys
      grib2_key_list_plus_int%additional_keys%nint_keys = i + 1
      grib2_key_list_plus_int%additional_keys%int_key(i+1) = b
    END IF
  END FUNCTION grib2_key_list_plus_int

  FUNCTION grib2_key_list_plus_intarr(a, b)
    TYPE(t_grib2_var) :: grib2_key_list_plus_intarr
    TYPE(t_grib2_var),        INTENT(IN) :: a
    TYPE(t_grib2_intarr_key), INTENT(IN) :: b
    ! local variables
    INTEGER :: i

    grib2_key_list_plus_intarr = a
    IF (a%additional_keys%nintarr_keys < MAX_INT_KEYS) THEN
      i = grib2_key_list_plus_intarr%additional_keys%nintarr_keys
      grib2_key_list_plus_intarr%additional_keys%nintarr_keys = i + 1
      grib2_key_list_plus_intarr%additional_keys%intarr_key(i+1) = b
    END IF
  END FUNCTION grib2_key_list_plus_intarr

  FUNCTION grib2_key_list_plus_dbl(a, b)
    TYPE(t_grib2_var) :: grib2_key_list_plus_dbl
    TYPE(t_grib2_var),     INTENT(IN) :: a
    TYPE(t_grib2_dbl_key), INTENT(IN) :: b
    ! local variables
    INTEGER :: i

    grib2_key_list_plus_dbl = a
    IF (a%additional_keys%ndbl_keys < MAX_DBL_KEYS) THEN
      i = grib2_key_list_plus_dbl%additional_keys%ndbl_keys
      grib2_key_list_plus_dbl%additional_keys%ndbl_keys = i + 1
      grib2_key_list_plus_dbl%additional_keys%dbl_key(i+1) = b
    END IF
  END FUNCTION grib2_key_list_plus_dbl

  FUNCTION grib2_key_list_plus_dblarr(a, b)
    TYPE(t_grib2_var) :: grib2_key_list_plus_dblarr
    TYPE(t_grib2_var),        INTENT(IN) :: a
    TYPE(t_grib2_dblarr_key), INTENT(IN) :: b
    ! local variables
    INTEGER :: i

    grib2_key_list_plus_dblarr = a
    IF (a%additional_keys%ndblarr_keys < MAX_DBL_KEYS) THEN
      i = grib2_key_list_plus_dblarr%additional_keys%ndblarr_keys
      grib2_key_list_plus_dblarr%additional_keys%ndblarr_keys = i + 1
      grib2_key_list_plus_dblarr%additional_keys%dblarr_key(i+1) = b
    END IF
  END FUNCTION grib2_key_list_plus_dblarr

END MODULE mo_grib2
