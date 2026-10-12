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

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_interpolate_time

  USE mo_kind,              ONLY: wp, i8
  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
       &                               calculate_time_interpolation_weights
  USE mo_parallel_config,   ONLY: nproma
  USE mo_impl_constants,    ONLY: MAX_CHAR_LENGTH
  USE mtime,                ONLY: datetime, max_datetime_str_len,             &
       &                          julianday, getJulianDayFromDatetime,        &
       &                          datetimetostring, getDatetimeFromJulianDay, &
       &                          juliandelta, OPERATOR(<), OPERATOR(>),      &
       &                          OPERATOR(>=), OPERATOR(-), ASSIGNMENT(=),   &
       &                          no_of_ms_in_a_day

  USE mo_time_config,    ONLY: time_config
  USE mo_reader_abstract,   ONLY: t_abstract_reader, t_abstract_indexed_reader, t_reader_cyclic

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_time_intp
  PUBLIC :: t_time_intp_transient
  PUBLIC :: t_time_intp_monthlyclim
  PUBLIC :: intModeLinearMonthlyClim
  PUBLIC :: intModeLinear

  !> Abstract base class for common functionality of the interpolators.
  !!
  !! This class takes care of the interpolation. Behavior is modified by subclasses implementing
  !! functions deciding when to read new data, how to read the data, and what interpolation
  !! weight to use. Subclasses should provide the first two time steps in `dataa` and `datab` on
  !! initialization and push both to GPU via an `ENTER DATA CREATE`.
  TYPE, ABSTRACT :: t_time_intp
    REAL(wp), ALLOCATABLE :: dataa(:,:,:,:)
    REAL(wp), ALLOCATABLE :: datab(:,:,:,:)
  CONTAINS
    PROCEDURE :: intp => time_intp_intp

    PROCEDURE(need_new_data), PRIVATE, DEFERRED :: need_new_data
    PROCEDURE(read_new_data), PRIVATE, DEFERRED :: read_new_data
    PROCEDURE(get_weight), PRIVATE, DEFERRED :: get_weight
    PROCEDURE(get_npromz), PRIVATE, DEFERRED :: get_npromz

    PROCEDURE :: finalize => time_intp_final
  END TYPE

  ABSTRACT INTERFACE
    !> Signal that new data should be loaded.
    LOGICAL FUNCTION need_new_data(this, local_time)
      IMPORT :: t_time_intp, datetime
      CLASS(t_time_intp), INTENT(IN) :: this
      TYPE(datetime), INTENT(IN) :: local_time
    END FUNCTION

    !> Read the next time slice into `data`.
    SUBROUTINE read_new_data(this, local_time, data)
      IMPORT :: t_time_intp, datetime, wp
      CLASS(t_time_intp), INTENT(INOUT) :: this
      TYPE(datetime), INTENT(IN) :: local_time
      REAL(wp), ALLOCATABLE, INTENT(INOUT) :: data(:,:,:,:)
    END SUBROUTINE

    !> Compute the weight of the newer time step.
    REAL(wp) FUNCTION get_weight(this, local_time)
      IMPORT :: t_time_intp, datetime, wp
      CLASS(t_time_intp), INTENT(IN) :: this
      TYPE(datetime), INTENT(IN) :: local_time
    END FUNCTION

    !> Get length of the last block.
    INTEGER FUNCTION get_npromz(this)
      IMPORT :: t_time_intp
      CLASS(t_time_intp), INTENT(IN) :: this
    END FUNCTION
  END INTERFACE

  !> Interpolation for transient data.
  TYPE, EXTENDS(t_time_intp) :: t_time_intp_transient
    TYPE(julianday) :: time_old !< Timestamp of the old timestep.
    TYPE(julianday) :: time_new !< Timestamp of the new timestep.

    !> Name of variable to interpolate.
    CHARACTER(len=MAX_CHAR_LENGTH) :: var_name

    !> Interpolation mode. Options are constant, linear, and weird linear. Default is linear.
    INTEGER :: interpolation_mode

    CLASS(t_abstract_reader), POINTER :: reader !< Underlying reader.

  CONTAINS
    PROCEDURE :: init => time_intp_transient_init

    PROCEDURE :: need_new_data => time_intp_transient_need_new_data
    PROCEDURE :: read_new_data => time_intp_transient_read_new_data
    PROCEDURE :: get_weight => time_intp_transient_get_weight
    PROCEDURE :: get_npromz => time_intp_transient_get_npromz

    FINAL :: time_intp_transient_final_r0, time_intp_transient_final_r1
  END TYPE t_time_intp_transient

  !> Interpolation for monthly climatologies.
  TYPE, EXTENDS(t_time_intp) :: t_time_intp_monthlyclim
    !> Month that dataa refers to.
    INTEGER :: month1

    !> Name of variable to interpolate.
    CHARACTER(len=MAX_CHAR_LENGTH) :: var_name

    TYPE(t_reader_cyclic) :: reader !< Underlying reader.

  CONTAINS
    PROCEDURE :: init => time_intp_monthlyclim_init

    PROCEDURE :: need_new_data => time_intp_monthlyclim_need_new_data
    PROCEDURE :: read_new_data => time_intp_monthlyclim_read_new_data
    PROCEDURE :: get_weight => time_intp_monthlyclim_get_weight
    PROCEDURE :: get_npromz => time_intp_monthlyclim_get_npromz

    FINAL :: time_intp_monthlyclim_final_r0, time_intp_monthlyclim_final_r1
  END TYPE t_time_intp_monthlyclim

  INTEGER, PARAMETER :: intModeConstant          = 0
  INTEGER, PARAMETER :: intModeLinear            = 1
  INTEGER, PARAMETER :: intModeLinearMonthlyClim = 2
  INTEGER, PARAMETER :: intModeLinearWeird       = 11

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interpolate_time'

CONTAINS

  !> Interpolate field for `local_time`.
  SUBROUTINE time_intp_intp(this, local_time, interpolated, lacc)
    CLASS(t_time_intp), INTENT(INOUT), TARGET :: this
    TYPE(datetime), INTENT(IN) :: local_time !< Timestamp for interpolation.
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: interpolated(:,:,:,:) !< Interpolated data.
    LOGICAL, INTENT(IN) :: lacc !< OpenACC flag.

    REAL(wp) :: weight
    REAL(wp), POINTER, CONTIGUOUS :: dataa(:,:,:,:)
    REAL(wp), POINTER, CONTIGUOUS :: datab(:,:,:,:)

    INTEGER :: nblks, nlev, npromz, nlen
    INTEGER :: jw, jb, jk, jc

#ifndef _OPENACC
    ! Suppress unused warning.
    IF (lacc) THEN; END IF
#endif

    IF (this%need_new_data(local_time)) THEN
      BLOCK
        REAL(wp), ALLOCATABLE :: temp(:,:,:,:)

        ! Swap dataa and datab.
        CALL MOVE_ALLOC(this%dataa, temp)
        CALL MOVE_ALLOC(this%datab, this%dataa)
        CALL MOVE_ALLOC(temp, this%datab)
      END BLOCK

      CALL this%read_new_data(local_time, this%datab)
      !$ACC UPDATE DEVICE(this%datab) ASYNC(1)
    END IF

    weight = this%get_weight(local_time)

    IF (ALLOCATED(interpolated)) THEN
      IF (ANY(SHAPE(interpolated) /= SHAPE(this%dataa) )) THEN
        !$ACC EXIT DATA DELETE(interpolated) ASYNC(1)
        DEALLOCATE(interpolated)
      END IF
    END IF

    IF (.NOT. ALLOCATED(interpolated)) THEN
      ALLOCATE(interpolated(SIZE(this%dataa,1), SIZE(this%dataa,2), SIZE(this%dataa,3), SIZE(this%dataa,4)))
      !$ACC ENTER DATA CREATE(interpolated) ASYNC(1)
    END IF

    dataa => this%dataa
    datab => this%datab

    !$ACC DATA PRESENT(interpolated, dataa, datab) ASYNC(1) IF(lacc)

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lacc)
    interpolated(:,:,:,:) = 0.0_wp
    !$ACC END KERNELS

    npromz = this%get_npromz()
    nblks  = size(interpolated,3)
    nlev   = size(interpolated,2)

    ! we need this mess, since npromz == nproma is not garantueed
    DO jw = 1,size(interpolated,4)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1,nblks
        nlen = MERGE(nproma, npromz, jb /= nblks)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1,nlev
          DO jc = 1,nlen
            interpolated(jc,jk,jb,jw) = (1.0_wp-weight) * dataa(jc,jk,jb,jw) &
              &                                +weight  * datab(jc,jk,jb,jw)
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDDO

    !$ACC END DATA

  END SUBROUTINE time_intp_intp

  !> Finalize the base class. Deallocates data on GPU.
  SUBROUTINE time_intp_final (this)
    CLASS(t_time_intp), INTENT(INOUT) :: this

    IF (.NOT. ALLOCATED(this%dataa)) RETURN

    !$ACC EXIT DATA DELETE(this%dataa, this%datab) FINALIZE ASYNC(1)

    DEALLOCATE(this%dataa, this%datab)

  END SUBROUTINE time_intp_final

  !> Initialize transient interpolator. Checks data boundaries and loads the first two time steps.
  SUBROUTINE time_intp_transient_init(this, reader, local_time, var_name, int_mode)
    CLASS(t_time_intp_transient), TARGET, INTENT(OUT) :: this
    CLASS(t_abstract_reader), TARGET, INTENT(INOUT) :: reader !< Underlying reader.
    TYPE(datetime), INTENT(IN) :: local_time !< Initial time.
    CHARACTER(len=*), INTENT(IN) :: var_name !< Name of variable to interpolate.
    INTEGER, OPTIONAL, INTENT(IN) :: int_mode !< Interpolation mode (intModeConstant, intModeLinear)

    CHARACTER(len=*), PARAMETER :: routine = modname//"::time_intp_transient::init"

    this%reader   => reader
    this%var_name =  var_name

    IF (PRESENT(int_mode)) THEN
      this%interpolation_mode = int_mode
    ELSE
      this%interpolation_mode = intModeLinear
    ENDIF

    SELECT CASE (this%interpolation_mode)
    CASE (intModeConstant, intModeLinear)
    CASE (intModeLinearWeird)
      ! This should
      ! a) be renamed. (But what is this?)
      ! b) done similar to calculate_time_interpolation_weights in
      !    shared/mo_bcs_time_interpolation.f90 . Since I do not know how this
      !    works and whether it can be generalized to data with non-monthly
      !    intervals, I leave as this for now.
      CALL finish(routine, "You are weird")
    CASE (intModeLinearMonthlyClim)
      CALL finish(routine, 'Use a t_time_intp_monthlyclim for monthly climatologies')
    END SELECT

    BLOCK
      TYPE(julianday) :: jd
      CHARACTER(len=max_datetime_str_len) :: date_str

      CALL this%reader%goto(time_config%tc_stopdate)

      IF (.NOT. this%reader%is_valid(message_text)) THEN
        CALL datetimetostring(time_config%tc_stopdate, date_str)
        CALL finish(routine, "End of run ("//TRIM(date_str)//") seek failed: " // TRIM(message_text))
      END IF

      CALL getJulianDayFromDatetime(time_config%tc_stopdate, jd)

      IF (jd > this%reader%get_julian_day()) THEN
        ! We need the next time step to complete the interpolation for the last simulation step.
        CALL this%reader%next()

        IF (.NOT. this%reader%is_valid(message_text)) THEN
          CALL datetimetostring(time_config%tc_stopdate, date_str)
          CALL finish(routine, "End of run ("//TRIM(date_str)//") out of bounds: " // TRIM(message_text))
        END IF
      END IF

      CALL this%reader%goto(time_config%tc_startdate)

      IF (.NOT. this%reader%is_valid(message_text)) THEN
        CALL datetimetostring(time_config%tc_startdate, date_str)
        CALL finish(routine, "Start of run ("//TRIM(date_str)//") seek failed: " // TRIM(message_text))
      END IF

      CALL this%reader%goto(local_time)

      IF (.NOT. this%reader%is_valid(message_text)) THEN
        CALL datetimetostring(local_time, date_str)
        CALL finish(routine, "Initial ("//TRIM(date_str)//") seek failed: " // TRIM(message_text))
      END IF
    END BLOCK

    CALL reader%read(this%var_name, this%dataa)
    this%time_old = reader%get_julian_day()

    CALL reader%next()

    CALL reader%read(this%var_name, this%datab)
    this%time_new = reader%get_julian_day()

    !$ACC ENTER DATA CREATE(this%dataa, this%datab) ASYNC(1)
    !$ACC UPDATE DEVICE(this%dataa, this%datab) ASYNC(1)

    log_output: BLOCK
      TYPE(datetime) :: load_time_in_file, load_next_time_in_file
      CHARACTER(len=max_datetime_str_len) :: date_str1, date_str2

      CALL getDatetimeFromJulianDay(this%time_old, load_time_in_file)
      CALL getDatetimeFromJulianDay(this%time_new, load_next_time_in_file)

      CALL datetimetostring(load_time_in_file, date_str1)
      CALL datetimetostring(load_next_time_in_file, date_str2)
      WRITE(message_text,*) &
            &         "loading ", TRIM(this%var_name), " data for "//TRIM(date_str1)//" and "//TRIM(date_str2)
      CALL message(TRIM(routine),message_text)
    END BLOCK log_output

  END SUBROUTINE time_intp_transient_init

  SUBROUTINE time_intp_transient_final_r0 (this)
    TYPE(t_time_intp_transient), INTENT(INOUT) :: this

    ! Cannot put this on the abstract base class because final routines have to take a TYPE(class)
    ! parameter but abstract classes cannot be used as TYPE(). Catch 22.
    CALL this%finalize
  END SUBROUTINE

  SUBROUTINE time_intp_transient_final_r1 (this)
    TYPE(t_time_intp_transient), INTENT(INOUT) :: this(:)
    INTEGER :: i

    DO i = 1, SIZE(this)
      CALL this(i)%finalize
    END DO
  END SUBROUTINE

  !> Check if transient interpolator needs new data.
  !! New data is needed when the current timestamp has crossed the timestamp of the newer data slice.
  LOGICAL FUNCTION time_intp_transient_need_new_data (this, local_time)
    CLASS(t_time_intp_transient), INTENT(IN) :: this
    TYPE(datetime), INTENT(IN) :: local_time

    TYPE(julianday) :: current_jd

    CALL getJulianDayFromDatetime(local_time, current_jd)

    time_intp_transient_need_new_data = (current_jd > this%time_new)

  END FUNCTION time_intp_transient_need_new_data

  !> Read new data for the transient interpolation.
  !! Updates timestamps for weight calculation.
  SUBROUTINE time_intp_transient_read_new_data (this, local_time, data)
    CLASS(t_time_intp_transient), INTENT(INOUT) :: this
    TYPE(datetime), INTENT(IN) :: local_time
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: data(:,:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname // '::time_intp_transient::read_new_data'

    TYPE(datetime) :: dt
    CHARACTER(len=max_datetime_str_len) :: dt_str

    ! Suppress unused warning.
    IF (local_time%date%year /= 0) THEN; END IF

    CALL this%reader%next()

    this%time_old = this%time_new
    this%time_new = this%reader%get_julian_day()

    IF (this%time_old >= this%time_new) THEN
      BLOCK
        CHARACTER(len=max_datetime_str_len) :: dt_str_old, dt_str_new

        CALL getDateTimeFromJulianDay(this%time_old, dt)
        CALL datetimeToString(dt, dt_str_old)
        CALL getDateTimeFromJulianDay(this%time_new, dt)
        CALL datetimeToString(dt, dt_str_new)

        CALL finish(routine, 'timestamps are not monotonically increasing (' // &
            & TRIM(dt_str_old) // ' >= ' // TRIM(dt_str_new) // ')')
      END BLOCK
    END IF

    CALL getDateTimeFromJulianDay(this%time_new, dt)
    CALL datetimeToString(dt, dt_str)
    CALL message(routine, 'loaded ' // TRIM(this%var_name) // ' data for ' // TRIM(dt_str))

    CALL this%reader%read(this%var_name, data)

  END SUBROUTINE time_intp_transient_read_new_data

  !> Compute weight of the newer data.
  FUNCTION time_intp_transient_get_weight (this, local_time) RESULT(weight)
    CLASS(t_time_intp_transient), INTENT(IN) :: this
    TYPE(datetime), INTENT(IN) :: local_time
    REAL(wp) :: weight

    TYPE(julianday) :: current_jd

    CALL getJulianDayFromDatetime(local_time, current_jd)

    weight = 0._wp

    SELECT CASE (this%interpolation_mode)
    CASE (intModeConstant)
      weight = 0._wp ! Weight of the new data is zero.
    CASE (intModeLinear)
      BLOCK
        TYPE(juliandelta) :: delta_1, delta_2
        REAL(wp) :: ds1, ds2

        delta_1 = current_jd - this%time_old
        delta_2 = this%time_new - this%time_old
        ds1 = 1.0e-3_wp * (no_of_ms_in_a_day * delta_1%day + delta_1%ms)
        ds2 = 1.0e-3_wp * (no_of_ms_in_a_day * delta_2%day + delta_2%ms)
        weight = ds1/ds2
      END BLOCK
    END SELECT

  END FUNCTION time_intp_transient_get_weight

  FUNCTION time_intp_transient_get_npromz (this) RESULT(npromz)
    CLASS(t_time_intp_transient), INTENT(IN) :: this
    INTEGER :: npromz

    npromz = this%reader%get_npromz()
  END FUNCTION time_intp_transient_get_npromz


  !> Initialize the interpolator for monthly climatologies.
  SUBROUTINE time_intp_monthlyclim_init(this, timeseries_reader, local_time, var_name)
    CLASS(t_time_intp_monthlyclim), TARGET, INTENT(OUT) :: this
    CLASS(t_abstract_indexed_reader), TARGET, INTENT(INOUT) :: timeseries_reader !< Underlying reader.
    TYPE(datetime), INTENT(IN) :: local_time !< Initial time
    CHARACTER(len=*), INTENT(IN) :: var_name !< Name of variable to interpolate.

    CHARACTER(len=*), PARAMETER :: routine = modname//"::time_intp_monthlyclim_init"

    TYPE(t_time_interpolation_weights) :: tiw

    this%var_name = var_name

    CALL timeseries_reader%seek(1_i8)

    ! Set up for 12 monthly means.
    CALL this%reader%init(timeseries_reader, 12)

    tiw = calculate_time_interpolation_weights(local_time)

    ! Move to first month.
    CALL this%reader%seek(INT(tiw%month1, i8))

    WRITE (message_text, '(a,i2,a,i2)') 'Loading data for months ', tiw%month1, ' and ', tiw%month2
    CALL message(routine, message_text)

    CALL this%reader%read(this%var_name, this%dataa)
    this%month1 = tiw%month1

    CALL this%reader%next()

    CALL this%reader%read(this%var_name, this%datab)

    !$ACC ENTER DATA CREATE(this%dataa, this%datab) ASYNC(1)
    !$ACC UPDATE DEVICE(this%dataa, this%datab) ASYNC(1)

  END SUBROUTINE

  SUBROUTINE time_intp_monthlyclim_final_r0 (this)
    TYPE(t_time_intp_monthlyclim), INTENT(INOUT) :: this

    ! Cannot put this on the abstract base class because final routines have to take a TYPE(class)
    ! parameter but abstract classes cannot be used as TYPE(). Catch 22.
    CALL this%finalize
  END SUBROUTINE

  SUBROUTINE time_intp_monthlyclim_final_r1 (this)
    TYPE(t_time_intp_monthlyclim), INTENT(INOUT) :: this(:)
    INTEGER :: i

    DO i = 1, SIZE(this)
      CALL this(i)%finalize
    END DO
  END SUBROUTINE

  LOGICAL FUNCTION time_intp_monthlyclim_need_new_data (this, local_time)
    CLASS(t_time_intp_monthlyclim), INTENT(IN) :: this
    TYPE(datetime), INTENT(IN) :: local_time

    TYPE(t_time_interpolation_weights) :: tiw

    tiw = calculate_time_interpolation_weights(local_time)

    time_intp_monthlyclim_need_new_data = (tiw%month1 /= this%month1)

  END FUNCTION time_intp_monthlyclim_need_new_data

  SUBROUTINE time_intp_monthlyclim_read_new_data (this, local_time, data)
    CLASS(t_time_intp_monthlyclim), INTENT(INOUT) :: this
    TYPE(datetime), INTENT(IN) :: local_time
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: data(:,:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//'::time_intp_monthlyclim::read_new_data'

    TYPE(t_time_interpolation_weights) :: tiw

    tiw = calculate_time_interpolation_weights(local_time)

    this%month1 = tiw%month1

    WRITE (message_text, '(a,I0)') 'Loading data for month ', tiw%month2
    CALL message(routine, message_text)

    CALL this%reader%next()

    CALL this%reader%read(this%var_name, data)

  END SUBROUTINE time_intp_monthlyclim_read_new_data

  FUNCTION time_intp_monthlyclim_get_weight (this, local_time) RESULT(weight)
    CLASS(t_time_intp_monthlyclim), INTENT(IN) :: this
    TYPE(datetime), INTENT(IN) :: local_time
    REAL(wp) :: weight

    TYPE(t_time_interpolation_weights) :: tiw

    ! Suppress unused warning.
    IF (this%month1 == 0) THEN; END IF

    tiw = calculate_time_interpolation_weights(local_time)

    weight = tiw%weight2

  END FUNCTION time_intp_monthlyclim_get_weight

  FUNCTION time_intp_monthlyclim_get_npromz (this) RESULT(npromz)
    CLASS(t_time_intp_monthlyclim), INTENT(IN) :: this
    INTEGER :: npromz

    npromz = this%reader%get_npromz()
  END FUNCTION time_intp_monthlyclim_get_npromz

END MODULE mo_interpolate_time
