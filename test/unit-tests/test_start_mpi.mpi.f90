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

! 4EXEC-JSON { "MPI_NTASKS": 3 }

PROGRAM test_start_mpi_parallel

  USE mo_mpi, ONLY: start_mpi, stop_mpi, abort_mpi, &
      process_mpi_all_comm, split_global_mpi_communicator
  USE mpi
  USE mo_test_common, ONLY: test_fail, test_pass
  USE mo_io_units, ONLY: nerr

  IMPLICIT NONE

  INTEGER :: ierror
  INTEGER :: comm
  INTEGER :: comm_rank, comm_size
  CHARACTER(*), PARAMETER :: modname = 'test_start_mpi_parallel'

  CALL start_mpi('test_mpi_init')
  CALL split_global_mpi_communicator(1, 1)

  comm = process_mpi_all_comm
  CALL mpi_comm_size(comm, comm_size, ierror)
  IF (ierror /= mpi_success) THEN
    WRITE(nerr, '(a,a)') modname, ' mpi_comm_size failed.'
    WRITE(nerr, '(a,i4)') ' Error =  ', ierror
    CALL abort_mpi
  END IF
  CALL mpi_comm_rank(comm, comm_rank, ierror)
  IF (ierror /= mpi_success) THEN
    WRITE(nerr, '(a,a)') modname, ' mpi_comm_rank failed.'
    WRITE(nerr, '(a,i4)') ' Error =  ', ierror
    CALL abort_mpi
  END IF

  CALL stop_mpi

  CALL test_pass

END PROGRAM test_start_mpi_parallel
