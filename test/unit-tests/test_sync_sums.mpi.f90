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

! 4EXEC-JSON {
! 4EXEC-JSON "MPI_NTASKS": 2,
! 4EXEC-JSON "ENV": { "OMP_NUM_THREADS" : 2 }
! 4EXEC-JSON }

#if defined(_OPENMP)
#include "omp_definitions.inc"
#endif

PROGRAM test_sync_sums_parallel

  USE mo_kind, ONLY: wp
  USE mo_test_common, ONLY: test_fail, test_pass
  USE mo_mpi, ONLY: start_mpi, stop_mpi, get_my_mpi_work_comm_size, &
      &                             get_my_mpi_work_id
  USE mo_sync, ONLY: global_sum_array, global_sum, omp_global_sum_array, &
      &                             global_min, global_max
  USE mo_parallel_config, ONLY: l_fast_sum !< Controls which sum routine used in {,omp_}global_sum_array
  USE mo_exception, ONLY: warning

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: modname = 'test_sync_sums_parallel'
  !INTEGER, PARAMETER :: nproma = X, nblocks=Y, nlev=Z

  REAL(wp), TARGET, ALLOCATABLE :: wp_data_3d(:, :, :)
  REAL(wp) :: ref_wp_sum
  INTEGER, PARAMETER :: data_shape(3) = (/3, 2, 12/)
  INTEGER, PARAMETER :: data_size = PRODUCT(data_shape)
  INTEGER :: init_data(data_size)
  INTEGER :: i, comm_size, comm_rank

  ALLOCATE(wp_data_3d(data_shape(1), data_shape(2), data_shape(3)))
  DO i = 1, data_size
    init_data(i) = i
  END DO
  wp_data_3d = RESHAPE(REAL(init_data, wp), SHAPE=data_shape)

  CALL start_mpi(modname)
  comm_size = get_my_mpi_work_comm_size()
  comm_rank = get_my_mpi_work_id()

  ! global_sum_array: Simple example using array of init_data=[1,2,3,...]
  ref_wp_sum = REAL(comm_size * SUM(init_data), wp)
  CALL check_global_sum_array(wp_data_3d, ref_wp_sum, l_fast_sum_in=.TRUE.)
  CALL check_global_sum_array(wp_data_3d, ref_wp_sum, l_fast_sum_in=.FALSE.)

  ! global_sum: Simple example using array of init_data=[1,2,3,...]
  CALL check_global_sum(wp_data_3d, REAL(comm_size * init_data, wp), data_size)

#ifndef __SINGLE_PRECISION
  ! omp_global_sum_array: Simple example using array of init_data=[1,2,3,...]
  CALL check_omp_global_sum_array(RESHAPE(wp_data_3d, (/data_size/)), ref_wp_sum, l_fast_sum_in=.TRUE.)
  CALL check_omp_global_sum_array(RESHAPE(wp_data_3d, (/data_size/)), ref_wp_sum, l_fast_sum_in=.FALSE.)
#endif

  ! global_min
  CALL check_global_min(RESHAPE(wp_data_3d, (/data_size/)) + comm_rank, &
      &                   RESHAPE(wp_data_3d, (/data_size/)))

  ! global_max
  CALL check_global_max(RESHAPE(wp_data_3d, (/data_size/)) + comm_rank, &
      &                   RESHAPE(wp_data_3d, (/data_size/)) + comm_size - 1)

  DEALLOCATE(wp_data_3d)
  CALL stop_mpi ! Exits with code 0, similar to test_pass
  CALL test_pass

CONTAINS

  ! Total sum (local and global reduction)
  SUBROUTINE check_global_sum_array(wp_data_3d, ref_wp_sum, l_fast_sum_in)
    REAL(wp), INTENT(IN) :: wp_data_3d(:, :, :)
    REAL(wp), INTENT(IN) :: ref_wp_sum
    LOGICAL, INTENT(IN) :: l_fast_sum_in
    ! Local vars
    REAL(wp) :: wp_sum
    CHARACTER(LEN=256) :: message_text

    l_fast_sum = l_fast_sum_in

    wp_sum = global_sum_array(wp_data_3d)

    IF (wp_sum /= ref_wp_sum) THEN
      WRITE(message_text, '(a,L4)') "global_sum_array_wp does not match. l_fast_sum=", l_fast_sum
      CALL test_fail(modname, message_text)
    END IF

  END SUBROUTINE check_global_sum_array


  ! Total sum (global reduction, no local reduction)
  SUBROUTINE check_global_sum(wp_data_3d, ref_wp_sum, nsize)
    INTEGER, INTENT(IN) :: nsize
    REAL(wp), INTENT(IN) :: wp_data_3d(nsize)
    REAL(wp), INTENT(IN) :: ref_wp_sum(nsize)
    ! Local vars
    REAL(wp) :: wp_sum(nsize)
    CHARACTER(LEN=256) :: message_text

    wp_sum = global_sum(wp_data_3d)

    IF (ANY(wp_sum /= ref_wp_sum)) THEN
      WRITE(message_text, '(a,L4)') "global_sum_wp does not match. l_fast_sum=", l_fast_sum
      CALL test_fail(modname, message_text)
    END IF

  END SUBROUTINE check_global_sum

#ifndef __SINGLE_PRECISION
  ! Total sum (local and global reduction)
  SUBROUTINE check_omp_global_sum_array(wp_data_1d, ref_wp_sum, l_fast_sum_in)
    REAL(wp), INTENT(IN) :: wp_data_1d(:)
    REAL(wp), INTENT(IN) :: ref_wp_sum
    LOGICAL, INTENT(IN) :: l_fast_sum_in
    ! Local vars
    REAL(wp) :: wp_sum
    CHARACTER(LEN=256) :: message_text

    l_fast_sum = l_fast_sum_in

    !ICON_OMP_PARALLEL
    wp_sum = omp_global_sum_array(wp_data_1d)
    !ICON_OMP_END_PARALLEL

    IF (wp_sum /= ref_wp_sum) THEN
      WRITE(message_text, '(a,f10.3,a,f10.3)') "wp_sum=", wp_sum, ", ref_wp_sum=", ref_wp_sum
      CALL warning(modname, message_text)
      WRITE(message_text, '(a,L4)') "global_sum_array_wp does not match. l_fast_sum=", l_fast_sum
      CALL test_fail(modname, message_text)
    END IF

  END SUBROUTINE check_omp_global_sum_array
#endif

  ! global min
  SUBROUTINE check_global_min(wp_data_1d, ref_wp_min)
    REAL(wp), INTENT(IN) :: wp_data_1d(:)
    REAL(wp), INTENT(IN) :: ref_wp_min(:)
    ! Local vars
    REAL(wp), ALLOCATABLE :: wp_min(:)
    CHARACTER(LEN=256) :: message_text

    ALLOCATE(wp_min(SIZE(wp_data_1d)))
    wp_min = global_min(RESHAPE(wp_data_1d, (/data_size/)))

    IF (ANY(wp_min /= ref_wp_min)) THEN
      WRITE(message_text, '(a)') "global_min does not match"
      CALL test_fail(modname, message_text)
    END IF

    DEALLOCATE(wp_min)
  END SUBROUTINE check_global_min


  ! global max
  SUBROUTINE check_global_max(wp_data_1d, ref_wp_max)
    REAL(wp), INTENT(IN) :: wp_data_1d(:)
    REAL(wp), INTENT(IN) :: ref_wp_max(:)
    ! Local vars
    REAL(wp), ALLOCATABLE :: wp_max(:)
    CHARACTER(LEN=256) :: message_text

    ALLOCATE(wp_max(SIZE(wp_data_1d)))
    wp_max = global_max(RESHAPE(wp_data_1d, (/data_size/)))

    IF (ANY(wp_max /= ref_wp_max)) THEN
      WRITE(message_text, '(a)') "global_max does not match"
      CALL test_fail(modname, message_text)
    END IF

    DEALLOCATE(wp_max)
  END SUBROUTINE check_global_max

END PROGRAM test_sync_sums_parallel
