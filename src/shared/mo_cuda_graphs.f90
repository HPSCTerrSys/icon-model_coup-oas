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

MODULE mo_cuda_graphs

  USE mo_impl_constants,        ONLY: max_dom, MAX_CHAR_LENGTH
  USE mo_exception,             ONLY: message, message_text, finish
  USE, INTRINSIC :: iso_c_binding

#ifdef ICON_USE_CUDA_GRAPH
  USE mo_run_config,            ONLY: msg_level
  USE mo_acc_device_management, ONLY: accGraph, accBeginCapture, accEndCapture, accGraphLaunch, accGraphDelete
#endif

  IMPLICIT NONE

  PUBLIC :: t_cuda_graphs, create_graphs, id_captured, &
            begin_capture, end_capture, replay, reset, reset_one, &
            reset_all_graphs

  TYPE t_cuda_graphs
#ifdef ICON_USE_CUDA_GRAPH
    TYPE(accGraph),      ALLOCATABLE :: cache(:)
    INTEGER(c_intptr_t), ALLOCATABLE :: keys(:,:)
    INTEGER :: total_captured, capturing_queue
#endif

    CHARACTER(len=MAX_CHAR_LENGTH) :: routine
    LOGICAL :: initialized = .FALSE.
  END TYPE t_cuda_graphs

  INTEGER :: num_registered_graphs = 0, num_registered_comm_graphs = 0
  TYPE p_cuda_graphs
    TYPE(t_cuda_graphs), POINTER :: p
  END TYPE p_cuda_graphs
  TYPE(p_cuda_graphs) :: registered_graphs(100)

  CONTAINS

  SUBROUTINE make_internal_keys(keys, routine, ptr_keys, int_keys)
    INTEGER(c_intptr_t),   INTENT(inout) :: keys(:)
    CHARACTER(len=*),      INTENT(in)    :: routine
    TYPE(c_ptr), OPTIONAL, INTENT(in)    :: ptr_keys(:)
    INTEGER,     OPTIONAL, INTENT(in)    :: int_keys(:)

    INTEGER :: tot_size, start

    tot_size = 0
    IF (PRESENT(ptr_keys)) tot_size = tot_size + SIZE(ptr_keys)
    IF (PRESENT(int_keys)) tot_size = tot_size + SIZE(int_keys)

    IF (SIZE(keys) /= tot_size) THEN
      WRITE(message_text,'(a,i2,a,i1)') 'wrong number of keys in CUDA graphs cache, expected ', &
        SIZE(keys), ' got ', tot_size
      CALL finish(routine, message_text)
    END IF

    start = 1
    IF (PRESENT(ptr_keys)) THEN
      keys(1:SIZE(ptr_keys)) = TRANSFER(ptr_keys, INT(1, c_intptr_t), SIZE(ptr_keys))
      start = SIZE(ptr_keys)+1
    END IF
    IF (PRESENT(int_keys)) THEN
      keys(start:) = INT(int_keys, c_intptr_t)
    END IF
  END SUBROUTINE make_internal_keys

  SUBROUTINE create_graphs(graphs, num_keys, routine)
    TYPE(t_cuda_graphs), TARGET, INTENT(out) :: graphs
    INTEGER,                     INTENT(in)  :: num_keys
    CHARACTER(len=*),            INTENT(in)  :: routine

    INTEGER :: max_graphs

    graphs%routine = routine

#ifdef ICON_USE_CUDA_GRAPH
    max_graphs = 5**num_keys * max_dom ! assume no more than 5 variants per dimension
    ALLOCATE( graphs%cache(max_graphs), graphs%keys(num_keys, max_graphs) )
    graphs%keys = 0
    graphs%total_captured = 0
    graphs%initialized = .TRUE.

    num_registered_graphs = num_registered_graphs + 1
    registered_graphs(num_registered_graphs)%p => graphs
#endif
  END SUBROUTINE create_graphs

  FUNCTION id_captured(graphs, ptr_keys, int_keys)
    TYPE(t_cuda_graphs),   INTENT(in) :: graphs
    TYPE(c_ptr), OPTIONAL, INTENT(in) :: ptr_keys(:)
    INTEGER,     OPTIONAL, INTENT(in) :: int_keys(:)
    INTEGER                           :: id_captured

#ifdef ICON_USE_CUDA_GRAPH
    INTEGER :: ig
    INTEGER(c_intptr_t) :: keys(SIZE(graphs%keys, 1))
#endif

#ifdef ICON_USE_CUDA_GRAPH
    id_captured = 0
    CALL make_internal_keys(keys, graphs%routine, ptr_keys, int_keys)
    DO ig=1,graphs%total_captured
      IF (ALL(keys == graphs%keys(:,ig))) THEN
        id_captured = ig
        EXIT
      END IF
    END DO
#else
    id_captured = -1
#endif
  END FUNCTION id_captured

  SUBROUTINE begin_capture(graphs, async_queue, ptr_keys, int_keys)
    TYPE(t_cuda_graphs),   INTENT(inout) :: graphs
    INTEGER,               INTENT(in)    :: async_queue
    TYPE(c_ptr), OPTIONAL, INTENT(in)    :: ptr_keys(:)
    INTEGER,     OPTIONAL, INTENT(in)    :: int_keys(:)

#ifdef ICON_USE_CUDA_GRAPH
    INTEGER(c_intptr_t) :: keys(SIZE(graphs%keys, 1))

    CALL make_internal_keys(keys, graphs%routine, ptr_keys, int_keys)
    graphs%total_captured = graphs%total_captured + 1

    IF (graphs%total_captured > SIZE(graphs%keys, 2)) THEN
      CALL finish(graphs%routine, 'captured too many CUDA graphs')
    END IF

    IF (msg_level >= 11) THEN
      WRITE(message_text,'(a,i2)') 'capturing CUDA graph id ', graphs%total_captured
      CALL message(graphs%routine, message_text)
    END IF

    graphs%keys(:,graphs%total_captured) = keys
    graphs%capturing_queue = async_queue
    CALL accBeginCapture(async_queue)
#endif
  END SUBROUTINE begin_capture

  FUNCTION end_capture(graphs)
    TYPE(t_cuda_graphs), INTENT(inout) :: graphs
    INTEGER                            :: end_capture

#ifdef ICON_USE_CUDA_GRAPH
    IF (msg_level >= 11) THEN
      WRITE(message_text,'(a,i2)') 'finished capturing CUDA graph id ', graphs%total_captured
      CALL message(graphs%routine, message_text)
    END IF

    CALL accEndCapture(graphs%capturing_queue, graphs%cache(graphs%total_captured))
    end_capture = graphs%total_captured
#else
    ! Always return dummy 0 here
    end_capture = 0
#endif
  END FUNCTION end_capture

  SUBROUTINE replay(graphs, id, async_queue)
    TYPE(t_cuda_graphs), INTENT(in) :: graphs
    INTEGER,             INTENT(in) :: id
    INTEGER,             INTENT(in) :: async_queue

#ifdef ICON_USE_CUDA_GRAPH
    IF (id < 1 .OR. id > graphs%total_captured) THEN
      WRITE(message_text,'(a,i2,a,i2)') 'passed invalid CUDA graph id ', id, &
        ', expected between 1 and ', graphs%total_captured
      CALL finish(graphs%routine, message_text)
    END IF

    IF (msg_level >= 12) THEN
      WRITE(message_text,'(a,i2)') 'executing CUDA graph id ', id
      CALL message(graphs%routine, message_text)
    END IF

    CALL accGraphLaunch(graphs%cache(id), async_queue)
#else
    CALL finish(graphs%routine, 'this routine should never be called when CUDA graphs are disabled')
#endif
  END SUBROUTINE replay

  SUBROUTINE reset(graphs)
    TYPE(t_cuda_graphs), INTENT(inout) :: graphs
    INTEGER :: ig

#ifdef ICON_USE_CUDA_GRAPH
    DO ig = 1, graphs%total_captured
      ! NVHPC WAR: removing graphs may trigger a segfault due to a prior memory error
      ! We will re-enable the cleanup when the bug is fixed
      ! CALL accGraphDelete(graphs%cache(ig))
    END DO
    graphs%total_captured = 0
    graphs%keys = 0
#endif
  END SUBROUTINE reset

  SUBROUTINE reset_one(graphs, id)
    TYPE(t_cuda_graphs), INTENT(inout) :: graphs
    INTEGER,             INTENT(in)    :: id

#ifdef ICON_USE_CUDA_GRAPH
    IF (id < 1 .OR. id > graphs%total_captured) THEN
      WRITE(message_text,'(a,i2,a,i2)') 'passed invalid CUDA graph id ', id, &
        ', expected between 1 and ', graphs%total_captured
      CALL finish(graphs%routine, message_text)
    END IF

    ! NVHPC WAR: removing graphs may trigger a segfault due to a prior memory error
    ! We will re-enable the cleanup when the bug is fixed
    ! CALL accGraphDelete(graphs%cache(id))
    IF (graphs%total_captured > 1) THEN
      graphs%keys(:,id) = graphs%keys(:,graphs%total_captured)
      graphs%cache(id) = graphs%cache(graphs%total_captured)
    END IF

    graphs%keys(:,graphs%total_captured) = 0
    graphs%total_captured = graphs%total_captured - 1
#endif
  END SUBROUTINE reset_one

  SUBROUTINE reset_all_graphs()
    INTEGER :: ig

#ifdef ICON_USE_CUDA_GRAPH
    CALL message("mo_cuda_graphs::reset_all_graphs", &
                 "resetting all the CUDA graphs (likely due to variables reallocation)")
    DO ig = 1, num_registered_graphs
      CALL reset(registered_graphs(ig)%p)
    END DO
#endif
  END SUBROUTINE reset_all_graphs

END MODULE mo_cuda_graphs
