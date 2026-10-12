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

MODULE mo_ragnarok_f2c
  USE ISO_C_BINDING, ONLY: c_int, c_double, c_ptr, c_null_ptr, c_loc, c_funptr, &
      & c_null_funptr, c_bool
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_f2c_ftable
  PUBLIC :: t_dom_info
  PUBLIC :: t_f2c_patch_descr
  PUBLIC :: t_patch_info
  PUBLIC :: init_ragnarok_f2c
  PUBLIC :: t_comm_patch
  PUBLIC :: t_comm_pattern_cdescr
  PUBLIC :: ragnarok_sync_patch_array_r3_sync_c
  PUBLIC :: t_process_info

  TYPE, BIND(c) :: t_f2c_ftable
    TYPE(c_funptr) :: get_domain_info = c_null_funptr
    TYPE(c_funptr) :: get_patch_info = c_null_funptr
    TYPE(c_funptr) :: get_mo_model_domain_p_patch_descr = c_null_funptr
    TYPE(c_funptr) :: get_comm_patch = c_null_funptr
    TYPE(c_funptr) :: f2c_message = c_null_funptr
    TYPE(c_funptr) :: f2c_finish = c_null_funptr
    TYPE(c_funptr) :: f2c_exchange_data_r3d = c_null_funptr
    TYPE(c_funptr) :: get_process_info = c_null_funptr
    TYPE(c_funptr) :: f2c_new_timer = c_null_funptr
    TYPE(c_funptr) :: f2c_timer_start = c_null_funptr
    TYPE(c_funptr) :: f2c_timer_stop = c_null_funptr
    TYPE(c_funptr) :: f2c_timer_value = c_null_funptr
    TYPE(c_funptr) :: f2c_get_timer_config = c_null_funptr
  END TYPE t_f2c_ftable

  ! domain info:
  TYPE, BIND(c) :: t_dom_info
    INTEGER(c_int) :: nproma = 0
    INTEGER(c_int) :: id_min = -1
    INTEGER(c_int) :: id_max = -2
  END TYPE t_dom_info

  TYPE, BIND(c) :: t_f2c_patch_descr
    TYPE(c_ptr) :: cptr
  END TYPE t_f2c_patch_descr

  ! portable part of top level t_patch components
  TYPE, BIND(c) :: t_patch_info
    INTEGER(c_int) :: id = -1
    INTEGER(c_int) :: nlev = 0
    INTEGER(c_int) :: nblks_c = 0
    INTEGER(c_int) :: nblks_e = 0
    INTEGER(c_int) :: nblks_v = 0
  END TYPE t_patch_info

  TYPE, BIND(c) :: t_comm_pattern_cdescr
    TYPE(c_ptr) :: cptr
  END TYPE t_comm_pattern_cdescr

  TYPE, BIND(c) :: t_comm_patch
    TYPE(t_comm_pattern_cdescr) :: comm_pat_c
    TYPE(t_comm_pattern_cdescr) :: comm_pat_e
    TYPE(t_comm_pattern_cdescr) :: comm_pat_v
    TYPE(t_comm_pattern_cdescr) :: comm_pat_c1
  END TYPE t_comm_patch

  TYPE, BIND(c) :: t_process_info
    LOGICAL(c_bool) :: is_mpi_parallel
  END TYPE t_process_info

  INTERFACE
    SUBROUTINE init_ragnarok_f2c(funtab) BIND(c)
      IMPORT t_f2c_ftable
      TYPE(t_f2c_ftable), INTENT(IN) :: funtab
    END SUBROUTINE init_ragnarok_f2c


    SUBROUTINE ragnarok_sync_patch_array_r3_sync_c(pid, arr, arr_shape) BIND(c)
      IMPORT :: c_double, c_int
      INTEGER(c_int), VALUE :: pid
      REAL(c_double) :: arr(*)
      INTEGER(c_int), INTENT(IN) :: arr_shape(3)
    END SUBROUTINE ragnarok_sync_patch_array_r3_sync_c
  END INTERFACE

END MODULE mo_ragnarok_f2c
