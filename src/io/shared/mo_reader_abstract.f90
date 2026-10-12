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

MODULE mo_reader_abstract

  USE mo_kind,         ONLY: wp, i8
  USE mo_model_domain, ONLY: t_patch
  USE mtime,           ONLY: julianday, datetime, getDatetimeFromJulianDay, &
      & getJulianDayFromDatetime, OPERATOR(<), OPERATOR(<=), OPERATOR(>), OPERATOR(==)

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_abstract_reader
  PUBLIC :: t_abstract_indexed_reader

  PUBLIC :: t_reader_cyclic

  CHARACTER(len=*), PARAMETER :: modname = 'mo_abstract_reader'

  !> Abstract reader for time-series data.
  !!
  !! Models a bidirectional iterator over an unbounded sequence. Users must set
  !! a starting point by calling `goto` to get to a specified time stamp in the
  !! sequence. From there, `prev` and `next` can be used to move to adjacent
  !! elements. `is_valid` checks for validity of the iterator state (e.g. the
  !! actual bounds of the dataset may have been reached, or a file could not be
  !! opened). The `read` and `get_julian_day` methods return data and timestamps
  !! and may only be called on valid iterators.
  !!
  !! The `next` and `prev` operations can be used to back out of an invalid state
  !! unless that state resulted from a `goto` operation, in which case the iter-
  !! ator might remain in that invalid state until a `goto` operation succeeds.
  !!
  TYPE, ABSTRACT :: t_abstract_reader
  CONTAINS
    PROCEDURE(abstract_deinit),   DEFERRED :: deinit
    PROCEDURE(abstract_prev),     DEFERRED :: prev
    PROCEDURE(abstract_next),     DEFERRED :: next
    PROCEDURE(abstract_goto),     DEFERRED :: goto
    PROCEDURE(abstract_is_valid), DEFERRED :: is_valid

    PROCEDURE(abstract_read),           DEFERRED :: read
    PROCEDURE(abstract_get_julian_day), DEFERRED :: get_julian_day

    ! we need those, since we only assume the read data to be blocked
    ! we make no assumptions for the data to be placed on cells,
    ! edges, lonlat grids, ...
    PROCEDURE(abstract_get_nblks),  DEFERRED :: get_nblks
    PROCEDURE(abstract_get_npromz), DEFERRED :: get_npromz
  END TYPE t_abstract_reader

  ABSTRACT INTERFACE
    SUBROUTINE abstract_prev (this)
      IMPORT :: t_abstract_reader
      CLASS(t_abstract_reader), INTENT(INOUT) :: this
    END SUBROUTINE abstract_prev

    SUBROUTINE abstract_next (this)
      IMPORT :: t_abstract_reader
      CLASS(t_abstract_reader), INTENT(INOUT) :: this
    END SUBROUTINE abstract_next

    SUBROUTINE abstract_goto (this, target_datetime)
      IMPORT :: t_abstract_reader, datetime
      CLASS(t_abstract_reader), INTENT(INOUT) :: this
      TYPE(datetime), INTENT(IN) :: target_datetime
    END SUBROUTINE abstract_goto

    LOGICAL FUNCTION abstract_is_valid (this, msg)
      IMPORT :: t_abstract_reader
      CLASS(t_abstract_reader), INTENT(IN) :: this
      CHARACTER(len=*), INTENT(OUT), OPTIONAL :: msg
    END FUNCTION abstract_is_valid

    SUBROUTINE abstract_read (this, varname, dat)
      IMPORT :: t_abstract_reader, wp
      CLASS(t_abstract_reader), INTENT(INOUT) :: this
      CHARACTER(len=*), INTENT(IN) :: varname
      REAL(wp), ALLOCATABLE, INTENT(INOUT) :: dat(:,:,:,:)
    END SUBROUTINE abstract_read

    TYPE(julianday) FUNCTION abstract_get_julian_day (this)
      IMPORT :: t_abstract_reader, julianday
      CLASS(t_abstract_reader), INTENT(IN) :: this
    END FUNCTION

    SUBROUTINE abstract_deinit (this)
      IMPORT :: t_abstract_reader
      CLASS(t_abstract_reader), INTENT(inout) :: this
    END SUBROUTINE abstract_deinit

    FUNCTION abstract_get_nblks (this) RESULT(nblks)
      IMPORT :: t_abstract_reader
      CLASS(t_abstract_reader), INTENT(in   ) :: this
      INTEGER                                 :: nblks
    END FUNCTION

    FUNCTION abstract_get_npromz (this) RESULT(npromz)
      IMPORT :: t_abstract_reader
      CLASS(t_abstract_reader), INTENT(in   ) :: this
      INTEGER                                 :: npromz
    END FUNCTION
  END INTERFACE

  !> Reader interface that extends the underlying unbounded sequence with a notion of an absolute
  !! index. Adds a method that allows to seek to a specified index instead of a timestamp. Index
  !! values should be consecutive, starting at 1.
  TYPE, ABSTRACT, EXTENDS(t_abstract_reader) :: t_abstract_indexed_reader
  CONTAINS
    PROCEDURE(abstract_seek), DEFERRED :: seek
    PROCEDURE(abstract_get_index), DEFERRED :: get_index
  END TYPE

  ABSTRACT INTERFACE
    !> Seek to the specified index.
    SUBROUTINE abstract_seek (this, index)
      IMPORT :: t_abstract_indexed_reader, i8
      CLASS(t_abstract_indexed_reader), INTENT(INOUT) :: this
      INTEGER(i8), INTENT(IN) :: index
    END SUBROUTINE

    !> Get the current index.
    INTEGER(i8) FUNCTION abstract_get_index (this)
      IMPORT :: t_abstract_indexed_reader, i8
      CLASS(t_abstract_indexed_reader), INTENT(IN) :: this
    END FUNCTION
  END INTERFACE

  TYPE(julianday), PARAMETER, PRIVATE :: INVALID_JD = julianday(0, 100000000)

  !> Cyclic reader adaptor.
  !! Rolls up a reader into a cycle of specified length. The reader's current point is
  !! the starting point of the cycle, to which the cyclic reader will return after
  !! traversing the other `cycle_length - 1` elements.
  !!
  !! A `goto` operation moves to the closest point on the cycle before the given time,
  !!
  TYPE, EXTENDS(t_abstract_indexed_reader) :: t_reader_cyclic
    PRIVATE

    INTEGER :: cycle_length
    INTEGER :: index

    TYPE(julianday) :: start = INVALID_JD, end = INVALID_JD

    CLASS(t_abstract_reader), POINTER :: reader => NULL()
  CONTAINS

    PROCEDURE :: init => reader_cyclic_init

    PROCEDURE :: prev => reader_cyclic_prev
    PROCEDURE :: next => reader_cyclic_next
    PROCEDURE :: goto => reader_cyclic_goto
    PROCEDURE :: seek => reader_cyclic_seek
    PROCEDURE :: is_valid => reader_cyclic_is_valid

    PROCEDURE :: read => reader_cyclic_read
    PROCEDURE :: get_julian_day => reader_cyclic_get_julian_day
    PROCEDURE :: get_index => reader_cyclic_get_index

    PROCEDURE :: deinit => reader_cyclic_deinit

    PROCEDURE :: get_nblks => reader_cyclic_get_nblks
    PROCEDURE :: get_npromz => reader_cyclic_get_npromz

  END TYPE

CONTAINS

  SUBROUTINE reader_cyclic_init (this, reader, cycle_length)
    CLASS(t_reader_cyclic), INTENT(OUT) :: this
    CLASS(t_abstract_reader), TARGET, INTENT(INOUT) :: reader !< Underlying reader.
    INTEGER, INTENT(IN) :: cycle_length !< Length of the cycle.

    this%cycle_length = cycle_length
    this%index = 1
    this%start = reader%get_julian_day()
    this%reader => reader
  END SUBROUTINE reader_cyclic_init

  SUBROUTINE reader_cyclic_prev (this)
    CLASS(t_reader_cyclic), INTENT(INOUT) :: this

    INTEGER :: i
    TYPE(datetime) :: dt

    this%index = this%index - 1

    IF (this%index == 0) THEN
      IF (this%end == INVALID_JD) THEN
        DO i = 1, this%cycle_length - 1
          CALL this%reader%next()
        END DO
        this%end = this%reader%get_julian_day()
      ELSE
        CALL getDatetimeFromJulianDay(this%end, dt)
        CALL this%reader%goto(dt)
      END IF
    ELSE IF (this%index > 0) THEN
      CALL this%reader%prev()
    ELSE
      this%index = 0
    END IF
  END SUBROUTINE reader_cyclic_prev

  SUBROUTINE reader_cyclic_next (this)
    CLASS(t_reader_cyclic), INTENT(INOUT) :: this

    INTEGER :: i
    TYPE(datetime) :: dt

    this%index = this%index + 1

    IF (this%index > this%cycle_length) THEN
      CALL getDatetimeFromJulianDay(this%start, dt)
      CALL this%reader%goto(dt)
    ELSE
      CALL this%reader%next()
    END IF
  END SUBROUTINE reader_cyclic_next

  SUBROUTINE reader_cyclic_goto (this, target_datetime)
    CLASS(t_reader_cyclic), INTENT(INOUT) :: this
    TYPE(datetime), INTENT(IN) :: target_datetime

    TYPE(julianday) :: jd

    INTEGER :: i

    CALL getJulianDayFromDatetime (target_datetime, jd)

    IF (this%get_julian_day() > jd) THEN
      DO i = this%index - 1, 0, -1
        CALL this%reader%prev()
        this%index = i

        IF (this%get_julian_day() <= jd) RETURN
      END DO
    ELSE IF (this%get_julian_day() < jd) THEN
      DO i = this%index + 1, this%cycle_length
        CALL this%reader%next()
        this%index = i

        IF (this%get_julian_day() > jd) THEN
          CALL this%reader%prev()
          this%index = i - 1
          RETURN
        END IF
      END DO
    END IF

  END SUBROUTINE reader_cyclic_goto

  SUBROUTINE reader_cyclic_seek (this, index)
    CLASS(t_reader_cyclic), INTENT(INOUT) :: this
    INTEGER(i8), INTENT(IN) :: index

    INTEGER :: i

    IF (index > this%index) THEN
      DO i = this%index + 1, index
        CALL this%reader%next()
      END DO
    ELSE IF (index < this%index) THEN
      DO i = this%index - 1, index, -1
        CALL this%reader%prev()
      END DO
    END IF

    this%index = INT(index)

  END SUBROUTINE reader_cyclic_seek

  FUNCTION reader_cyclic_is_valid (this, msg) RESULT(valid)
    CLASS(t_reader_cyclic), INTENT(IN) :: this
    CHARACTER(len=*), INTENT(OUT), OPTIONAL :: msg
    LOGICAL :: valid

    IF (this%index <= 0) THEN
      valid = .FALSE.
      IF (PRESENT(msg)) msg = 'Index out of bounds'
      RETURN
    END IF

    valid = this%reader%is_valid(msg)
  END FUNCTION reader_cyclic_is_valid

  SUBROUTINE reader_cyclic_read (this, varname, dat)
    CLASS(t_reader_cyclic), INTENT(INOUT) :: this
    CHARACTER(len=*), INTENT(IN) :: varname
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: dat(:,:,:,:)

    CALL this%reader%read(varname, dat)
  END SUBROUTINE reader_cyclic_read

  FUNCTION reader_cyclic_get_julian_day (this) RESULT(jd)
    CLASS(t_reader_cyclic), INTENT(IN) :: this
    TYPE(julianday) :: jd

    jd = this%reader%get_julian_day()
  END FUNCTION reader_cyclic_get_julian_day

  FUNCTION reader_cyclic_get_index (this) RESULT(index)
    CLASS(t_reader_cyclic), INTENT(IN) :: this
    INTEGER(i8) :: index

    index = this%index

  END FUNCTION reader_cyclic_get_index

  SUBROUTINE reader_cyclic_deinit (this)
    CLASS(t_reader_cyclic), INTENT(INOUT) :: this

    CALL this%reader%deinit
  END SUBROUTINE reader_cyclic_deinit

  FUNCTION reader_cyclic_get_nblks (this) RESULT(nblks)
    CLASS(t_reader_cyclic), INTENT(IN) :: this
    INTEGER :: nblks

    nblks = this%reader%get_nblks()
  END FUNCTION reader_cyclic_get_nblks

  FUNCTION reader_cyclic_get_npromz (this) RESULT(npromz)
    CLASS(t_reader_cyclic), INTENT(IN) :: this
    INTEGER :: npromz

    npromz = this%reader%get_npromz()
  END FUNCTION reader_cyclic_get_npromz

END MODULE mo_reader_abstract
