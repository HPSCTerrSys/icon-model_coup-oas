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

! Utilities used for icon unit tests.
!   test_pass() -> Call if test passed
!   test_skip() -> Call to skip test
!   test_fail() -> Call if test failed

MODULE mo_test_common
  USE mo_util_system, ONLY: util_exit
  USE mo_io_units, ONLY: nerr
#ifndef NO_MPI
  USE mpi
#endif

  IMPLICIT NONE
  PUBLIC :: test_pass, test_skip, test_fail
CONTAINS

  ! Exits program and returns code 0.
  !   - fails if called in mpi region
  SUBROUTINE test_pass(name)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: name

    IF (is_mpi_live()) THEN
      CALL test_fail(name, "test_pass() called before mpi_finalize")
    END IF

    IF (PRESENT(name)) WRITE(nerr, '(a,a)') name, ' passed.'
    CALL util_exit(0)
  END SUBROUTINE test_pass


  SUBROUTINE test_skip()
    CALL util_exit(77)
  END SUBROUTINE test_skip


  SUBROUTINE test_fail(name, text)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: name
    CHARACTER(*), OPTIONAL, INTENT(IN) :: text
#ifndef NOMPI
    INTEGER :: p_error
#endif

    IF (PRESENT(name) .AND. PRESENT(text)) WRITE(nerr, '(a)') name//':'//text

#ifndef NOMPI
    IF (is_mpi_live()) THEN
      CALL MPI_ABORT(MPI_COMM_WORLD, 2, p_error)
    END IF
#endif
    CALL util_exit(2)
  END SUBROUTINE test_fail


  FUNCTION is_mpi_live()
    LOGICAL :: is_mpi_live
#ifndef NOMPI
    INTEGER :: p_error
    LOGICAL :: l_initialized, l_finalized

    CALL MPI_INITIALIZED(l_initialized, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE(nerr, '(a,i4)') ' MPI_INITIALIZED check failed. Error = ', p_error
      CALL util_exit(2)
    END IF

    CALL MPI_FINALIZED(l_finalized, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE(nerr, '(a,i4)') ' MPI_FINALIZED check failed. Error = ', p_error
      CALL util_exit(2)
    END IF

    is_mpi_live =(l_initialized .AND. .NOT.l_finalized)
#else
    is_mpi_live = .FALSE.
#endif

  END FUNCTION is_mpi_live

END MODULE mo_test_common
