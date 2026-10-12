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

PROGRAM test_index_list

  USE mo_index_list, ONLY: generate_index_list, generate_index_list_batched
  USE mo_kind, ONLY: i4
  USE mo_test_common, ONLY: test_fail, test_pass
  USE mo_exception, ONLY: warning

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: modname = 'test_index_list'
  INTEGER, PARAMETER :: n = 20000
  INTEGER, PARAMETER :: nb = 7
  INTEGER(i4), TARGET :: conditions(n, nb)
  INTEGER, TARGET :: indices(n, nb), dev_indices(n, nb)
  INTEGER :: nvalid(nb), dev_nvalid(nb)
  REAL :: harvest(n, nb)

  INTEGER :: i, b, ic
  LOGICAL :: lacc

#ifdef _OPENACC
  lacc = .TRUE.
#else
  lacc = .FALSE.
#endif

  CALL RANDOM_NUMBER(harvest)
  conditions = INT(harvest * 2)

  !$ACC DATA COPYIN(conditions) CREATE(dev_indices, dev_nvalid)

  ! Test the non-batched version

  nvalid(1) = 0
  DO i = 1, n
    IF (conditions(i, 1) /= 0) THEN
      nvalid(1) = nvalid(1) + 1
      indices(nvalid(1), 1) = i
    END IF
  END DO

  CALL generate_index_list(conditions(:, 1), dev_indices(:, 1), 1, n, dev_nvalid(1), lacc=lacc, opt_acc_async_queue=1)
  !$ACC UPDATE HOST(dev_indices(:,1)) ASYNC(1)
  !$ACC WAIT(1)

  PRINT *, "CHECK NON-BATCHED: ", nvalid(1) == dev_nvalid(1), ALL(indices(:nvalid(1), 1) == dev_indices(:nvalid(1), 1))

  ! Test the non-batched async version

  CALL generate_index_list(conditions(:, 1), dev_indices(:, 1), 1, n, dev_nvalid(1), lacc=lacc, opt_acc_async_queue=1, &
      &                       opt_acc_copy_to_host=.FALSE.)
  !$ACC UPDATE HOST(dev_indices(:,1), dev_nvalid(1)) ASYNC(1)
  !$ACC WAIT(1)

  PRINT *, "CHECK NON-BATCHED ASYNC: ", nvalid(1) == dev_nvalid(1), &
      &                                    ALL(indices(:nvalid(1), 1) == dev_indices(:nvalid(1), 1))

  ! Test the non-batched version with a shift

  ic = 42
  nvalid(1) = 0
  DO i = ic, n
    IF (conditions(i, 1) /= 0) THEN
      nvalid(1) = nvalid(1) + 1
      indices(nvalid(1), 1) = i
    END IF
  END DO

  CALL generate_index_list(conditions(:, 1), dev_indices(:, 1), ic, n, dev_nvalid(1), lacc=lacc, opt_acc_async_queue=1)
  !$ACC UPDATE HOST(dev_indices(:,1)) ASYNC(1)
  !$ACC WAIT(1)

  PRINT *, "CHECK NON-BATCHED SHIFTED: ", nvalid(1) == dev_nvalid(1), &
      &                                      ALL(indices(:nvalid(1), 1) == dev_indices(:nvalid(1), 1))

  ! Test the batched version

  indices(:, :) = 0
  !$ACC KERNELS
  dev_indices(:, :) = 0
  !$ACC END KERNELS

  nvalid = 0
  DO b = 1, nb
    DO i = 1, n
      IF (conditions(i, b) /= 0) THEN
        nvalid(b) = nvalid(b) + 1
        indices(nvalid(b), b) = i
      END IF
    END DO
  END DO

  CALL generate_index_list_batched(conditions, dev_indices, 1, n, dev_nvalid, &
      &   lacc=.TRUE., opt_acc_async_queue=1)
  !$ACC UPDATE HOST(dev_indices, dev_nvalid) ASYNC(1)
  !$ACC WAIT(1)

  PRINT *, "CHECK BATCHED: ", ALL(nvalid == dev_nvalid), ALL(indices == dev_indices)

  ! Test the batched shifted version

  indices(:, :) = 0
  !$ACC KERNELS
  dev_indices(:, :) = 0
  !$ACC END KERNELS

  ic = 142
  nvalid = 0
  DO b = 1, nb
    DO i = ic, n
      IF (conditions(i, b) /= 0) THEN
        nvalid(b) = nvalid(b) + 1
        indices(nvalid(b), b) = i
      END IF
    END DO
  END DO

  CALL generate_index_list_batched(conditions, dev_indices, ic, n, dev_nvalid, lacc, opt_acc_async_queue=1)
  !$ACC UPDATE HOST(dev_indices, dev_nvalid) ASYNC(1)
  !$ACC WAIT(1)

  PRINT *, "CHECK BATCHED SHIFTED: ", ALL(nvalid == dev_nvalid), ALL(indices == dev_indices)

  !$ACC END DATA

  IF (ANY(nvalid /= dev_nvalid) .OR. ANY(indices /= dev_indices)) THEN
    CALL warning(modname, "")
    CALL test_fail()
  END IF
  CALL test_pass()

END PROGRAM test_index_list
