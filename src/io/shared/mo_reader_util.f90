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

MODULE mo_reader_util
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t
  USE mo_kind, ONLY: dp
  USE mo_exception, ONLY: finish, message_text

  USE mo_netcdf
  USE mtime, ONLY: julianday, juliandelta, datetime, timedelta, newDatetime, newTimedelta, &
      & getJulianDayFromDatetime, getDatetimeFromJulianDay, deallocateDatetime, deallocateTimedelta, &
      & timeDeltaToJulianDelta, OPERATOR(*), OPERATOR(+), OPERATOR(-), OPERATOR(>), OPERATOR(<=), OPERATOR(==), &
      & OPERATOR(<), no_of_ms_in_a_day

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_timestamps_from_netcdf
  PUBLIC :: shift_time
  PUBLIC :: divide_time

  !> Shift a timestamp by a number of units. This should be equal to `base + unit * len`,
  !! but mtime cannot deal with scaling where the scaled timedelta exceeds 24 hours.
  INTERFACE shift_time
    MODULE PROCEDURE shift_time_dp
    MODUlE PROCEDURE shift_time_i8
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_reader_util'

CONTAINS

  !> Read the timestamps from a NetCDF file, returning them as julian days.
  !! Conversion of time units is done automatically.
  FUNCTION read_timestamps_from_netcdf (filename, times, timevar, timedim) RESULT(err)
    CHARACTER(len=*), INTENT(IN) :: filename !< Name of the NetCDF file.
    TYPE(julianday), ALLOCATABLE, INTENT(OUT) :: times(:) !< Timestamps as julian days.
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: timevar !< Name of time variable (default: time).
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: timedim !< Name of the time dimension (default: time).
    INTEGER :: err !< NetCDF error code (0 if successful).

    INTEGER :: i
    INTEGER :: fileid
    INTEGER :: ntimes
    INTEGER :: tvid
    INTEGER :: tdid

    REAL(dp), ALLOCATABLE :: times_read(:)
    CHARACTER(len=NF90_MAX_NAME) :: timevar_s
    CHARACTER(len=NF90_MAX_NAME) :: timedim_s
    CHARACTER(len=NF90_MAX_NAME) :: cf_timeaxis_string
    TYPE(datetime), POINTER :: epoch
    TYPE(timedelta), POINTER :: time_unit

    IF (PRESENT(timevar)) THEN
      timevar_s = timevar
    ELSE
      timevar_s = 'time'
    END IF

    IF (PRESENT(timedim)) THEN
      timedim_s = timedim
    ELSE
      timedim_s = 'time'
    END IF

    err = nf90_open(TRIM(filename), nf90_nowrite, fileid); IF (err /= 0) RETURN
    err = nf90_inq_varid(fileid, timevar_s, tvid); IF (err /= 0) RETURN
    err = nf90_inq_dimid(fileid, timedim_s, tdid); IF (err /= 0) RETURN
    err = nf90_inquire_dimension(fileid, tdid, len=ntimes); IF (err /= 0) RETURN

    ALLOCATE(times_read(ntimes))

    err = nf90_get_var(fileid, tvid, times_read); IF (err /= 0) RETURN
    err = nf90_get_att(fileid, tvid, "units", cf_timeaxis_string); IF (err /= 0) RETURN
    err = nf90_close(fileid); IF (err /= 0) RETURN

    CALL get_cf_timeaxis_desc(TRIM(cf_timeaxis_string), epoch, time_unit)

    ALLOCATE(times(ntimes))

    DO i = 1, ntimes
      times(i) = shift_time(epoch, time_unit, times_read(i))
    END DO

    CALL deallocateTimedelta(time_unit)
    CALL deallocateDatetime(epoch)

  END FUNCTION read_timestamps_from_netcdf


  !> Shift a timestamp by a number of units. This should be equal to `base + unit * len`,
  !! but mtime cannot deal with scaling where the scaled timedelta exceeds 24 hours.
  !!
  !! The result is returned as a julianday for convenience of use in the SST/SIC reader.
  !!
  !! Fractional `len` is handled in two steps: first, the function adds the integer part
  !! of the total number of months `len * (12 * unit%year + unit%month)` to the base.
  !! Then, the remaining fractional part is handled by computing the length of the final
  !! month as a juliandelta and multiplying the result by the fractional part. The shorter-
  !! than-monthly components are always the same length independent of the starting point.
  !! They are added in a final step by calling the function with the same `len` but with
  !! a unit that has year and month zeroed out. So, a unit of `P1M2D` and `len` of 1.5
  !! results in `base + 1 month + 0.5 months + 3 days`.
  RECURSIVE FUNCTION shift_time_dp (base, unit, len) RESULT(jd)
    TYPE(datetime), INTENT(IN) :: base !< Base datetime.
    TYPE(timedelta), INTENT(IN) :: unit !< Unit duration of the shift.
    REAL(dp), INTENT(IN) :: len !< Number of units to shift.
    TYPE(julianday) :: jd

    IF (unit%year == 0 .AND. unit%month == 0) THEN
      BLOCK
        TYPE(juliandelta) :: junit, jdelta
        INTEGER(c_int64_t) :: ms_pre
        TYPE(julianday) :: base_jd

        ! All differences are of constant size. Convert to juliandelta, multiply, add
        CALL timeDeltaToJulianDelta(unit, base, junit)

        jdelta%sign = MERGE('+', '-', len >= 0._dp .EQV. junit%sign == '+')

        ms_pre = NINT(ABS(len) * (junit%ms + junit%day * no_of_ms_in_a_day), KIND=c_int64_t)
        jdelta%ms = MOD(ms_pre, INT(no_of_ms_in_a_day, KIND=c_int64_t))
        jdelta%day = ms_pre / no_of_ms_in_a_day

        CALL getJulianDayFromDatetime(base, base_jd)
        jd = base_jd + jdelta
     END BLOCK
    ELSE
      BLOCK
        TYPE(timedelta) :: big_step
        TYPE(timedelta) :: small_unit
        REAL(dp) :: step_in_months
        INTEGER :: whole_months
        REAL(dp) :: month_fraction
        TYPE(datetime) :: intermediate_dt

        ! Complicated case: Split into monthly and submonthly scales. Submonthly scales are independent of the starting point and
        ! are handled by the base case. So we first move by the total number of months, accounting for fractional months by linear
        ! interpolation, and then add the submonthly shift. Thus a unit of P1M2D and a length of 1.5 translates to `base + 1.5
        ! months + 3 days`.
        step_in_months = ABS(len) * (12 * unit%year + unit%month)
        whole_months = FLOOR(step_in_months)
        month_fraction = step_in_months - whole_months

        small_unit = unit
        small_unit%year = 0
        small_unit%month = 0

        big_step = timedelta( &
            & flag_std_form=1, &
            & sign=MERGE('+', '-', len >= 0._dp .EQV. unit%sign == '+'), &
            & year=whole_months / 12, &
            & month=MOD(whole_months, 12), &
            & day=0, hour=0, minute=0, second=0, ms=0)

        intermediate_dt = base + big_step

        IF (month_fraction > 1._dp / (no_of_ms_in_a_day * 31)) THEN
          BLOCK
            TYPE(timedelta) :: one_month
            TYPE(juliandelta) :: jone_month
            TYPE(juliandelta) :: partial_month
            TYPE(julianday) :: intermediate_jd
            INTEGER(c_int64_t) :: ms_pre

            one_month = timedelta( &
                & flag_std_form=1, &
                & sign=MERGE('+', '-', len >= 0._dp .EQV. unit%sign == '+'), &
                & year=0, month=1, &
                & day=0, hour=0, minute=0, second=0, ms=0)

            CALL timeDeltaToJulianDelta(one_month, intermediate_dt, jone_month)

            ms_pre = NINT(month_fraction * (jone_month%ms + jone_month%day * no_of_ms_in_a_day), KIND=c_int64_t)
            partial_month%sign = jone_month%sign
            partial_month%ms = MOD(ms_pre, INT(no_of_ms_in_a_day, KIND=c_int64_t))
            partial_month%day = ms_pre / no_of_ms_in_a_day

            CALL getJulianDayFromDatetime(intermediate_dt, intermediate_jd)
            intermediate_jd = intermediate_jd + partial_month
            CALL getDateTimeFromJulianDay(intermediate_jd, intermediate_dt)
          END BLOCK
        END IF

        jd = shift_time_dp(intermediate_dt, small_unit, len)

      END BLOCK
    END IF
  END FUNCTION


  !> Same as `shift_time_dp` but with an integer `len`.
  RECURSIVE FUNCTION shift_time_i8 (base, unit, len) RESULT(jd)
    TYPE(datetime), INTENT(IN) :: base !< Base datetime.
    TYPE(timedelta), INTENT(IN) :: unit !< Unit duration of the shift.
    INTEGER(c_int64_t), INTENT(IN) :: len !< Number of units to shift.
    TYPE(julianday) :: jd

    IF (unit%year == 0 .AND. unit%month == 0) THEN
      BLOCK
        TYPE(juliandelta) :: junit, jdelta
        INTEGER(c_int64_t) :: ms_pre
        TYPE(julianday) :: base_jd

        ! All differences are of constant size. Convert to juliandelta, multiply, add
        CALL timeDeltaToJulianDelta(unit, base, junit)

        jdelta%sign = MERGE('+', '-', len >= 0 .EQV. junit%sign == '+')

        ms_pre = INT(ABS(len) * (junit%ms + junit%day * no_of_ms_in_a_day), KIND=c_int64_t)
        jdelta%ms = MOD(ms_pre, INT(no_of_ms_in_a_day, KIND=c_int64_t))
        jdelta%day = ms_pre / no_of_ms_in_a_day

        CALL getJulianDayFromDatetime(base, base_jd)
        jd = base_jd + jdelta
     END BLOCK
    ELSE
      BLOCK
        TYPE(timedelta) :: big_step
        TYPE(timedelta) :: small_unit
        INTEGER(c_int64_t) :: step_in_months
        TYPE(datetime) :: intermediate_dt

        ! Complicated case: Split into monthly and submonthly scales. Submonthly scales are independent of the starting point and
        ! are handled by the base case. So we first move by the total number of months, and then add the submonthly shift.
        ! Thus a unit of P1M2D and a length of 2 translates to `base + 2 months + 4 days`.
        step_in_months = ABS(len) * (12 * unit%year + unit%month)

        small_unit = unit
        small_unit%year = 0
        small_unit%month = 0

        big_step = timedelta( &
            & flag_std_form=1, &
            & sign=MERGE('+', '-', len >= 0 .EQV. unit%sign == '+'), &
            & year=step_in_months / 12, &
            & month=INT(MOD(step_in_months, 12_c_int64_t)), &
            & day=0, hour=0, minute=0, second=0, ms=0)

        intermediate_dt = base + big_step

        jd = shift_time_i8(intermediate_dt, small_unit, len)

      END BLOCK
    END IF
  END FUNCTION


  !> Divide the delta between `time` and `base` by `unit` such that
  !! `time == base + quot * unit + rem`. Both the delta and `unit`
  !! have to be positive.
  SUBROUTINE divide_time(base, time, unit, quot, rem)
    TYPE(datetime), INTENT(IN) :: base !< Base datetime.
    TYPE(datetime), INTENT(IN) :: time !< Target datetime.
    TYPE(timedelta), INTENT(IN) :: unit !< Unit duration.
    INTEGER(c_int64_t), INTENT(OUT) :: quot !< Number of full unit durations between `time` and `base`.
    TYPE(timedelta), OPTIONAL, INTENT(OUT) :: rem !< Remaining duration such that `base + quot * unit + rem = time`.

    TYPE(julianday) :: jd_base, jd, jd_shift
    TYPE(datetime) :: nearest_multiple
    INTEGER :: highbit, bit

    CHARACTER(len=*), PARAMETER :: routine = modname // '::divide_time'

    CALL getJulianDayFromDateTime(time, jd)
    CALL getJulianDayFromDateTime(base, jd_base)

    IF (jd < jd_base .OR. unit%sign /= '+') CALL finish(routine, 'distance and unit must be positive')

    ! Find the highest set bit in the quotient.
    DO highbit = 0, INT(BIT_SIZE(quot)-1)
      jd_shift = shift_time_i8(base, unit, ISHFT(1_c_int64_t, highbit))
      IF (jd_shift > jd) EXIT
    END DO

    IF (highbit == BIT_SIZE(quot)) CALL finish(routine, 'Quotient too big. Divide by zero?')

    highbit = highbit - 1
    quot = ISHFT(1_c_int64_t, highbit)

    ! For each of the lower bits: if setting the bit overshoots, the bit is 0. Else, the bit is 1.
    ! This performs base-2 long division.
    DO bit = highbit - 1, 0, -1
      jd_shift = shift_time_i8(base, unit, quot + ISHFT(1_c_int64_t, bit))
      IF (jd_shift <= jd) quot = quot + ISHFT(1_c_int64_t, bit)
      IF (jd_shift == jd) EXIT
    END DO

    IF (PRESENT(rem)) THEN
      CALL getDateTimeFromJulianDay(shift_time_i8(base, unit, quot), nearest_multiple)
      rem = time - nearest_multiple
    END IF

  END SUBROUTINE divide_time


  SUBROUTINE get_cf_timeaxis_desc(cf_timeaxis_string, epoch, base_timeaxis_unit)
    CHARACTER(len=*), INTENT(IN) :: cf_timeaxis_string
    TYPE(datetime), POINTER, INTENT(OUT) :: epoch
    TYPE(timedelta), POINTER, INTENT(OUT) :: base_timeaxis_unit

    CHARACTER(len=*), PARAMETER :: routine = modname // '::get_cf_timeaxis_desc'

    ! The CF convention allows for a timezone to be included. We will
    ! ignore that one for all , but gets stored to word(5), if
    ! provided, to keep the algorithm simple.

    CHARACTER(len=16) :: word(5)
    INTEGER :: pos1, pos2, n
    INTEGER :: year, month, day
    INTEGER :: hour, minute, second, millis
    INTEGER :: hour_pre, minute_pre, second_pre
    REAL(dp) :: seconds

    pos1 = 1; pos2 = 0; n = 0;
    word(:) = ""

    DO
      pos2 = INDEX(cf_timeaxis_string(pos1:), " ")
      IF (pos2 == 0) THEN
        n = n + 1
        word(n) = cf_timeaxis_string(pos1:)
        EXIT
      ENDIF
      n = n + 1
      word(n) = cf_timeaxis_string(pos1:pos1+pos2-2)
      pos1 = pos2+pos1
    ENDDO

    ! correct the date part
    parse_date: BLOCK
      INTEGER :: idx1, idx2

      ! Exclude first character because year might be negative.
      idx1 = INDEX(word(3)(2:), '-') + 1
      idx2 = INDEX(word(3)(idx1+1:), '-') + idx1
      READ(word(3)(      :idx1-1),*) year
      READ(word(3)(idx1+1:idx2-1),*) month
      READ(word(3)(idx2+1:      ),*) day
    END BLOCK parse_date

    hour_pre = 0
    minute_pre = 0
    seconds = 0

    IF (word(4) /= "") THEN
      parse_time: BLOCK
        INTEGER :: idx1, idx2
        idx1 = INDEX(word(4), ':')
        idx2 = INDEX(word(4)(idx1+1:), ':')+idx1
        READ(word(4)(      :idx1-1),*) hour_pre
        READ(word(4)(idx1+1:idx2-1),*) minute_pre
        READ(word(4)(idx2+1:      ),*) seconds
      END BLOCK parse_time
    ENDIF

    millis = NINT(1000 * (seconds - FLOOR(seconds)))
    second_pre = FLOOR(seconds)
    minute_pre = minute_pre + floordiv(second_pre, 60, second)
    hour_pre = hour_pre + floordiv(minute_pre, 60, minute)
    day = day + floordiv(hour_pre, 24, hour)

    epoch => newDatetime(year, month, day, hour, minute, second, millis)
    IF (.NOT. ASSOCIATED(epoch)) THEN
      WRITE (message_text,'("Invalid epoch: ", 2(a, :, " "), &
          &", interpreted as ", i0, 2("-", i2.2), "T", 3(i2.2, :, ":"))') &
          & TRIM(word(3)), TRIM(word(4)), year, month, day, hour, minute, second
      CALL finish(routine, message_text)
    END IF

    SELECT CASE (TRIM(word(1)))
    CASE('years', 'year')
      base_timeaxis_unit => newTimedelta('+', 1, 0, 0, 0, 0, 0, 0)
    CASE('months', 'month')
      base_timeaxis_unit => newTimedelta('+', 0, 1, 0, 0, 0, 0, 0)
    CASE('days', 'day', 'd')
      base_timeaxis_unit => newTimedelta('+', 0, 0, 1, 0, 0, 0, 0)
    CASE('hours', 'hour', 'hr', 'h')
      base_timeaxis_unit => newTimedelta('+', 0, 0, 0, 1, 0, 0, 0)
    CASE('minutes', 'minute', 'min')
      base_timeaxis_unit => newTimedelta('+', 0, 0, 0, 0, 1, 0, 0)
    CASE('seconds', 'second', 'sec', 's')
      base_timeaxis_unit => newTimedelta('+', 0, 0, 0, 0, 0, 1, 0)
    CASE DEFAULT
      CALL finish(routine, 'Unknown time axis unit: ' // TRIM(word(1)))
      base_timeaxis_unit => NULL()
    END SELECT

  END SUBROUTINE get_cf_timeaxis_desc

  FUNCTION floordiv(a, b, rem) RESULT(quot)
    INTEGER, INTENT(IN) :: a, b
    INTEGER, INTENT(OUT) :: rem
    INTEGER :: quot

    rem = MODULO(a, b)
    quot = (a - rem) / b
  END FUNCTION

END MODULE mo_reader_util
