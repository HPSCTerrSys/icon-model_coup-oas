! @brief Initialisation of oasis coupling
#ifdef COUP_OAS_ECLM

MODULE mo_oasis_coupling_frame

!  USE mo_kind                ,ONLY: wp
!  USE mo_mpi                 ,ONLY: global_mpi_communicator
  USE mo_oasis_define        ,ONLY: oas_comp_id,oas_comp_name,oas_error,kl_comm
  USE mod_oasis

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_oasis_coupling_frame' ! Output of module for debug

  PUBLIC :: construct_oasis_coupling, destruct_oasis_coupling, abort_oasis_icon

CONTAINS

  !>
  !! SUBROUTINE construct_oasis_coupling -- the initialisation for the coupling
  !! with OASIS3-MCT coupler
  SUBROUTINE construct_oasis_coupling

    ! init OASIS coupler
    CALL oasis_init_comp(oas_comp_id, oas_comp_name, oas_error)
    IF (oas_error /= 0) &
      CALL oasis_abort(oas_comp_id, str_module, 'Failure in oasis_init')
    CALL oasis_get_localcomm(kl_comm, oas_error)
    IF (oas_error /= 0) &
      CALL oasis_abort(oas_comp_id, str_module, 'Failure in oasis_get_localcomm')

!    CALL MPI_COMM_DUP (kl_comm, global_mpi_communicator, oas_error)
!    WRITE (nerr,'(a,a)') routine, ' DEBUG: MPI_COMM_DUP'
!    IF (oas_error /= 0) &
!      CALL oasis_abort(oas_comp_id, str_module, 'Failure in MPI_COMM_DUP')

!    yname = TRIM("icon_COUPLED")

  END SUBROUTINE construct_oasis_coupling

  !>
  !! SUBROUTINE destruct_oasis_coupling -- the termination for the coupling
  !! with OASIS3-MCT coupler
  SUBROUTINE destruct_oasis_coupling

    CALL oasis_terminate(oas_error)
    IF (oas_error /= 0) CALL oasis_abort(oas_comp_id, str_module, &
      'Failure in oasis_terminate')

  END SUBROUTINE destruct_oasis_coupling

  !>
  !! SUBROUTINE abort_oasis_icon -- calling oasis abort
  SUBROUTINE abort_oasis_icon(oas_message)

    CHARACTER(len=*), INTENT(in),optional :: oas_message   !< error message string

    CALL oasis_abort(oas_comp_id, oas_comp_name, oas_message)

  END SUBROUTINE abort_oasis_icon

END MODULE mo_oasis_coupling_frame

#endif
