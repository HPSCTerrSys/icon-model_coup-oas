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

PROGRAM test_insert_dimension
  USE mo_kind, ONLY: wp
  USE mo_test_common, ONLY: test_fail, test_pass
  USE mo_exception, ONLY: warning, message_text
  USE mo_fortran_tools, ONLY: insert_dimension
  USE ISO_C_BINDING, ONLY: c_double
  USE mo_util_timer, ONLY: util_gettimeofday
  IMPLICIT NONE
  INTEGER, PARAMETER :: num_test_iterations = 1000
  CHARACTER(*), PARAMETER :: modname = 'test_insert_dimension'
  INTEGER :: test_iteration
  CALL init_rng
  DO test_iteration = 1, num_test_iterations
    CALL test_wp
    CALL test_i
  END DO
CONTAINS
  SUBROUTINE test_wp
    CHARACTER(*), PARAMETER :: routine = modname//':test_wp'
    REAL(wp), TARGET, ALLOCATABLE :: a_wp(:, :)
    REAL(wp), POINTER :: p_wp(:, :, :)
    INTEGER :: a_shape_base(2), a_shape_use(2), a_offset(2), &
        a_stride(2), a_last_idx(2), dim_insert_pos
    INTEGER :: p_shape(3), slice_shape(2)
    REAL :: rand_val(2, 5)
    INTEGER :: ierr
    INTEGER, PARAMETER :: max_shape(2) = (/100, 1000/)
    CALL RANDOM_NUMBER(rand_val)
    a_shape_base = INT(max_shape * rand_val(:, 1) + 1)
    ALLOCATE(a_wp(a_shape_base(1), a_shape_base(2)), STAT=ierr)
    IF (ierr /= 0) CALL test_fail()

    a_shape_use = INT(a_shape_base * (0.5 * rand_val(:, 2) + 0.5) + 1)
    a_offset = INT((a_shape_base - a_shape_use) * rand_val(:, 3) + 1)
    dim_insert_pos = INT(rand_val(1, 4) * 3 + 1)
    IF (rand_val(2, 4) >= 0.5) THEN
      a_stride = INT(a_shape_use * rand_val(:, 5) + 1)
    ELSE
      a_stride = 1
    END IF
    a_last_idx = a_offset + a_shape_use - 1
#if 0
    WRITE(0, '(a,i0)') &
        "size(a, 1)=", SIZE(a_wp, 1), &
        "size(a, 2)=", SIZE(a_wp, 2), &
        "a_offset(1)=", a_offset(1), &
        "a_last_idx(1)=", a_last_idx(1), &
        "a_stride(1)=", a_stride(1), &
        "a_offset(2)=", a_offset(2), &
        "a_last_idx(2)=", a_last_idx(2), &
        "a_stride(2)=", a_stride(2)
#endif
    CALL RANDOM_NUMBER(a_wp)
    CALL insert_dimension(p_wp, a_wp(a_offset(1):a_last_idx(1):a_stride(1), &
        &                           a_offset(2):a_last_idx(2):a_stride(2)), &
        &                dim_insert_pos)
    p_shape = SHAPE(p_wp)
    slice_shape = SHAPE(a_wp(a_offset(1):a_last_idx(1):a_stride(1), &
        &                   a_offset(2):a_last_idx(2):a_stride(2)))
    SELECT CASE (dim_insert_pos)
    CASE (1)
      IF (ANY(SHAPE(p_wp(1, :, :)) /= slice_shape)) &
          CALL shape_problem(p_shape, slice_shape, dim_insert_pos)
      IF (ANY(p_wp(1, :, :) /= a_wp(a_offset(1):a_last_idx(1):a_stride(1), &
          &                        a_offset(2):a_last_idx(2):a_stride(2)))) &
          CALL incorrect_extraction(routine, dim_insert_pos)
    CASE (2)
      IF (ANY(SHAPE(p_wp(:, 1, :)) /= slice_shape)) &
          CALL shape_problem(p_shape, slice_shape, dim_insert_pos)
      IF (ANY(p_wp(:, 1, :) /= a_wp(a_offset(1):a_last_idx(1):a_stride(1), &
          &                        a_offset(2):a_last_idx(2):a_stride(2)))) &
          CALL incorrect_extraction(routine, dim_insert_pos)
    CASE (3)
      IF (ANY(SHAPE(p_wp(:, :, 1)) /= slice_shape)) &
          CALL shape_problem(p_shape, slice_shape, dim_insert_pos)
      IF (ANY(p_wp(:, :, 1) /= a_wp(a_offset(1):a_last_idx(1):a_stride(1), &
          &                        a_offset(2):a_last_idx(2):a_stride(2)))) THEN
        CALL incorrect_extraction(routine, dim_insert_pos)
      END IF
    END SELECT
    DEALLOCATE(a_wp)

    CALL test_pass()
  END SUBROUTINE test_wp


  SUBROUTINE test_i
    CHARACTER(*), PARAMETER :: routine = modname//':test_i'
    REAL, TARGET, ALLOCATABLE :: a_r(:, :)
    INTEGER, TARGET, ALLOCATABLE :: a_i(:, :)
    INTEGER, POINTER :: p_i(:, :, :)
    INTEGER :: a_shape_base(2), a_shape_use(2), a_offset(2), &
        a_stride(2), a_last_idx(2), dim_insert_pos
    INTEGER :: p_shape(3), slice_shape(2)
    REAL :: rand_val(2, 5)
    INTEGER :: ierr
    INTEGER, PARAMETER :: max_shape(2) = (/200, 1000/)
    CALL RANDOM_NUMBER(rand_val)
    a_shape_base = INT(max_shape * rand_val(:, 1) + 1)
    ALLOCATE(a_i(a_shape_base(1), a_shape_base(2)), &
        &   a_r(a_shape_base(1), a_shape_base(2)), &
        &   stat=ierr)
    IF (ierr /= 0) CALL test_fail()
    a_shape_use = INT(a_shape_base * (0.5 * rand_val(:, 2) + 0.5) + 1)
    a_offset = INT((a_shape_base - a_shape_use) * rand_val(:, 3) + 1)
    dim_insert_pos = INT(rand_val(1, 4) * 3 + 1)
    IF (rand_val(2, 4) >= 0.5) THEN
      a_stride = INT(a_shape_use * rand_val(:, 5) + 1)
    ELSE
      a_stride = 1
    END IF
    a_last_idx = a_offset + a_shape_use - 1
#if 0
    WRITE(0, '(a,i0)') &
        "size(a, 1)=", SIZE(a_i, 1), &
        "size(a, 2)=", SIZE(a_i, 2), &
        "a_offset(1)=", a_offset(1), &
        "a_last_idx(1)=", a_last_idx(1), &
        "a_stride(1)=", a_stride(1), &
        "a_offset(2)=", a_offset(2), &
        "a_last_idx(2)=", a_last_idx(2), &
        "a_stride(2)=", a_stride(2)
#endif
    CALL RANDOM_NUMBER(a_r)
    a_i = NINT(((a_r * 2.0) - 1.0) * REAL(HUGE(a_i)))
    CALL insert_dimension(p_i, a_i(a_offset(1):a_last_idx(1):a_stride(1), &
        &                         a_offset(2):a_last_idx(2):a_stride(2)), &
        &                dim_insert_pos)
    p_shape = SHAPE(p_i)
    slice_shape = SHAPE(a_i(a_offset(1):a_last_idx(1):a_stride(1), &
        &                  a_offset(2):a_last_idx(2):a_stride(2)))
    SELECT CASE (dim_insert_pos)
    CASE (1)
      IF (ANY(SHAPE(p_i(1, :, :)) /= slice_shape)) &
          CALL shape_problem(p_shape, slice_shape, dim_insert_pos)
      IF (ANY(p_i(1, :, :) /= a_i(a_offset(1):a_last_idx(1):a_stride(1), &
          &                      a_offset(2):a_last_idx(2):a_stride(2)))) &
          CALL incorrect_extraction(routine, dim_insert_pos)
    CASE (2)
      IF (ANY(SHAPE(p_i(:, 1, :)) /= slice_shape)) &
          CALL shape_problem(p_shape, slice_shape, dim_insert_pos)
      IF (ANY(p_i(:, 1, :) /= a_i(a_offset(1):a_last_idx(1):a_stride(1), &
          &                      a_offset(2):a_last_idx(2):a_stride(2)))) &
          CALL incorrect_extraction(routine, dim_insert_pos)
    CASE (3)
      IF (ANY(SHAPE(p_i(:, :, 1)) /= slice_shape)) &
          CALL shape_problem(p_shape, slice_shape, dim_insert_pos)
      IF (ANY(p_i(:, :, 1) /= a_i(a_offset(1):a_last_idx(1):a_stride(1), &
          &                      a_offset(2):a_last_idx(2):a_stride(2)))) &
          CALL incorrect_extraction(routine, dim_insert_pos)
    END SELECT
    DEALLOCATE(a_i, a_r)
  END SUBROUTINE test_i


  SUBROUTINE init_rng
    INTEGER :: rseed_size
    INTEGER, ALLOCATABLE :: rseed(:)
    INTEGER :: status
    CHARACTER(len=1) :: debug
    REAL(c_double) :: unix_time
    CALL RANDOM_SEED(size=rseed_size)
    ALLOCATE(rseed(rseed_size))
    rseed = 4711
    unix_time = util_gettimeofday()
    rseed(1) = IEOR(INT(unix_time), &
        &          INT((unix_time - FLOOR(unix_time)) * 1000000.0))
    CALL RANDOM_SEED(put=rseed)
    CALL GET_ENVIRONMENT_VARIABLE("DEBUG", debug, status=status)
    IF (0 == status .AND. "1" == debug) WRITE(0, '(a,i0)') 'rseed=', rseed(1)
  END SUBROUTINE init_rng


  SUBROUTINE shape_problem(p_shape, slice_shape, dim_insert_pos)
    INTEGER, INTENT(IN) :: p_shape(3), slice_shape(2), dim_insert_pos
    WRITE(0, '(a,i0)') 'dim_insert_pos=', dim_insert_pos
    WRITE(0, '(a,"(",i0,", ", i0, ", ", i0, ")")') 'p_shape=', p_shape
    WRITE(0, '(a,"(",i0,", ", i0, ")")') 'slice_shape=', slice_shape
    WRITE(0, '(a)') 'incorrect shape'
    FLUSH(0)
    CALL test_fail()
  END SUBROUTINE shape_problem


  SUBROUTINE incorrect_extraction(routine, call_id)
    CHARACTER(*), INTENT(IN) :: routine
    INTEGER, INTENT(IN) :: call_id
    WRITE(message_text, *) 'incorrect extraction call_id ', call_id
    CALL warning(routine, message_text)
    CALL test_fail()
  END SUBROUTINE incorrect_extraction
END PROGRAM test_insert_dimension
