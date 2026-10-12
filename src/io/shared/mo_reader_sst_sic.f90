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

MODULE mo_reader_sst_sic

  USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t

  USE mo_kind,                    ONLY: wp
  USE mo_parallel_config,         ONLY: nproma
  USE mo_exception,               ONLY: finish
  USE mo_reader_abstract,         ONLY: t_abstract_reader
  USE mo_reader_util,             ONLY: read_timestamps_from_netcdf, shift_time, divide_time
  USE mo_io_units,                ONLY: FILENAME_MAX
  USE mo_model_domain,            ONLY: t_patch
  USE mo_netcdf_errhandler,       ONLY: nf
  USE mo_netcdf
  USE mtime,                      ONLY: julianday, juliandelta, getjuliandayfromdatetime, &
       &                                datetime, newdatetime, deallocatedatetime,        &
       &                                datetimeToString, max_datetime_str_len,           &
       &                                timedelta, newTimedelta, moduloTimedelta,         &
       &                                deallocateTimedelta, getDateTimeFromJulianDay,    &
       &                                OPERATOR(+), OPERATOR(-), OPERATOR(*),            &
       &                                OPERATOR(>), ASSIGNMENT(=)
  USE mo_mpi,                     ONLY: my_process_is_stdio, my_process_is_mpi_workroot, &
       &                                process_mpi_root_id, p_comm_work, p_bcast
  USE mo_read_netcdf_distributed, ONLY: distrib_nf_open, distrib_read, distrib_nf_close, &
       &                                idx_lvl_blk, t_distrib_read_data
  USE fortran_support,            ONLY: t_ptr_3d_wp
  USE mo_util_string,             ONLY: t_keyword_list, associate_keyword, with_keywords, &
       &                                int2string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_sst_sic_reader

  !> A single file from a sequence of files.
  !!
  !! Invariants:
  !!  - file is in a closed (default-constructed) state iff `times` is unallocated.
  !!  - fileid = -1 iff underlying NetCDF file not open.
  !!  - times is allocated and contains the timestamps present in file.
  TYPE :: t_sst_sic_file
    CHARACTER(len=FILENAME_MAX) :: filename = ''
    TYPE(datetime) :: timebase
    INTEGER :: fileid = -1
    TYPE(julianday), ALLOCATABLE :: times(:)
    TYPE(t_distrib_read_data), POINTER :: dist_io_data => NULL()
    INTEGER :: last_netcdf_error = 0
  CONTAINS
    PROCEDURE :: is_open => sst_sic_file_is_open
    PROCEDURE :: close => sst_sic_file_close
    PROCEDURE :: get_slice => sst_sic_file_get_slice
    FINAL :: sst_sic_file_finalize
  END TYPE

  INTERFACE t_sst_sic_file
    MODULE PROCEDURE sst_sic_file_open
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE sst_sic_file_assign
  END INTERFACE

  !> Reader for SST/SIC files.
  !!
  !! Invariants:
  !!  - Reader is valid iff `index > 0` and `file%is_open()`.
  !! Behavior:
  !!  - Moving past the end of data makes the reader invalid, but moving the opposite direction
  !!    makes it valid again.
  !!  - After a `goto`, the reader might be in two invalid states: if the target file exists,
  !!    the reader points to one before the first element or one beyond the last, and can be
  !!    brought into a valid state by moving into the opposite direction. If the file does not
  !!    exist, the reader can only be brought into a valid state by performing another `goto`.
  TYPE, EXTENDS(t_abstract_reader) :: t_sst_sic_reader

    TYPE(t_patch), POINTER :: p_patch => NULL()
    CHARACTER(len=FILENAME_MAX) :: filename_pattern
    TYPE(timedelta) :: file_interval

    TYPE(t_sst_sic_file) :: file

    INTEGER :: index = 0

  CONTAINS

    PROCEDURE :: init => sst_sic_reader_init
    PROCEDURE :: prev => sst_sic_reader_prev
    PROCEDURE :: next => sst_sic_reader_next
    PROCEDURE :: goto => sst_sic_reader_goto
    PROCEDURE :: is_valid => sst_sic_reader_is_valid

    PROCEDURE :: read => sst_sic_reader_read
    PROCEDURE :: get_julian_day => sst_sic_reader_get_julian_day

    PROCEDURE :: get_nblks       => sst_sic_reader_get_nblks
    PROCEDURE :: get_npromz      => sst_sic_reader_get_npromz

    PROCEDURE :: deinit          => sst_sic_reader_deinit

  END TYPE t_sst_sic_reader

  CHARACTER(len=*), PARAMETER :: modname = 'mo_reader_sst_sic'

CONTAINS

  FUNCTION sst_sic_file_open (filename, timebase, dist_io_data) RESULT(file)
    CHARACTER(len=*), INTENT(IN) :: filename
    TYPE(datetime), INTENT(IN) :: timebase
    TYPE(t_distrib_read_data), POINTER, INTENT(IN) :: dist_io_data
    TYPE(t_sst_sic_file) :: file

    INTEGER :: ntimes

    file%filename = filename
    file%timebase = timebase
    file%dist_io_data => dist_io_data

    IF (my_process_is_mpi_workroot()) THEN
      file%last_netcdf_error = read_timestamps_from_netcdf(TRIM(file%filename), file%times)
      ntimes = SIZE(file%times)
    END IF

    CALL p_bcast(file%last_netcdf_error, process_mpi_root_id, p_comm_work)

    IF (file%last_netcdf_error /= 0) RETURN

    CALL p_bcast(ntimes, process_mpi_root_id, p_comm_work)
    IF (.NOT. ALLOCATED(file%times)) THEN
      ALLOCATE(file%times(ntimes))
    END IF

    CALL p_bcast(file%times(:)%day, process_mpi_root_id, p_comm_work)
    CALL p_bcast(file%times(:)%ms, process_mpi_root_id, p_comm_work)
  END FUNCTION sst_sic_file_open

  LOGICAL FUNCTION sst_sic_file_is_open (this)
    CLASS(t_sst_sic_file), INTENT(IN) :: this

    sst_sic_file_is_open = ALLOCATED(this%times)
  END FUNCTION sst_sic_file_is_open

  SUBROUTINE sst_sic_file_close (this)
    CLASS(t_sst_sic_file), INTENT(INOUT) :: this

    IF (this%fileid /= -1) THEN
      CALL distrib_nf_close(this%fileid)
      this%fileid = -1
    END IF

    IF (ALLOCATED(this%times)) THEN
      DEALLOCATE(this%times)
    END IF
  END SUBROUTINE sst_sic_file_close

  SUBROUTINE sst_sic_file_finalize (this)
    TYPE(t_sst_sic_file), INTENT(INOUT) :: this

    CALL this%close()
  END SUBROUTINE sst_sic_file_finalize

  SUBROUTINE sst_sic_file_assign (to, from)
    TYPE(t_sst_sic_file), INTENT(OUT) :: to
    TYPE(t_sst_sic_file), INTENT(IN) :: from

    to%dist_io_data => from%dist_io_data
    to%fileid = -1
    to%filename = from%filename
    to%last_netcdf_error = from%last_netcdf_error
    to%timebase = from%timebase
    IF (ALLOCATED(from%times)) THEN
      to%times = from%times
    END IF
  END SUBROUTINE sst_sic_file_assign


  SUBROUTINE sst_sic_file_get_slice (this, index, varname, dat)
    CLASS(t_sst_sic_file), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: index
    CHARACTER(len=*), INTENT(IN) :: varname
    REAL(wp), INTENT(INOUT), TARGET :: dat(:,:,:)
    TYPE(t_ptr_3d_wp) :: ptr(1)

    IF (this%fileid == -1) THEN
      this%fileid = distrib_nf_open(TRIM(this%filename))
    END IF

    ptr(1)%p => dat
    CALL distrib_read(this%fileid, varname, ptr, [this%dist_io_data], &
        & edim=[1], dimo=idx_lvl_blk, start_ext_dim=[index], end_ext_dim=[index])

    WHERE (dat < -1e10_wp)
      dat = -1._wp
    END WHERE

  END SUBROUTINE


  !> Initialize a SST/SIC reader for the given patch.
  !!
  !! The filename pattern supports the replacements `<year>`, `<month>`, `<day>`, `<hh>`,
  !! `<mm>`, `<ss>`, referring to the current simulation time. The file interval gives
  !! the interval after which the next file gets opened, anchored at the start of the
  !! current simulation year. E.g., with an interval of 'P30D', a new file will be opened
  !! YYYY-01-31 00:00:00, YYYY-03-02 00:00:00, etc.
  SUBROUTINE sst_sic_reader_init(this, p_patch, filename_pattern, file_interval)
    CLASS(t_sst_sic_reader),    INTENT(inout) :: this
    TYPE(t_patch),      TARGET, INTENT(in   ) :: p_patch !< Domain patch.
    CHARACTER(len=*),           INTENT(in   ) :: filename_pattern !< Filename pattern.
    CHARACTER(len=*), OPTIONAL, INTENT(in   ) :: file_interval !< File interval as ISO 8601 duration.

    TYPE(timedelta), POINTER :: td

    this%filename_pattern = TRIM(filename_pattern)
    this%p_patch => p_patch

    IF (PRESENT(file_interval)) THEN
      td => newTimedelta(file_interval)
    ELSE
      td => newTimedelta('P1M')
    END IF

    this%file_interval = td
    CALL deallocateTimedelta(td)

  END SUBROUTINE sst_sic_reader_init

  !> Move reader to specific timestamp.
  !! The move will open the file that is supposed to contain the timestamp, e.g. for an interval
  !! of 'P30D' and a timestamp of 2022-03-01 12:00:00, the file for 2022-01-31 00:00:00 will be
  !! opened because it should contain timestamps up to 2022-03-02 00:00:00, exclusive.
  !! The reader points to the last step in the file that is earlier than the target timestamp.
  !!
  !! If the file cannot be opened, the reader is left in an invalid state from which it can only
  !! recover by performing another `goto` operation.
  SUBROUTINE sst_sic_reader_goto (this, target_datetime)
    CLASS(t_sst_sic_reader), INTENT(INOUT) :: this
    TYPE(datetime), INTENT(IN) :: target_datetime !< Target timestamp.

    CHARACTER(len=FILENAME_MAX) :: filename
    CHARACTER(len=max_datetime_str_len) :: target_datetime_str
    TYPE(datetime) :: base
    TYPE(julianday) :: jd

    INTEGER(c_int64_t) :: quot
    INTEGER :: i, ntimes

    CALL getJulianDayFromDatetime(target_datetime, jd)

    this%index = 0

    ! Base for filename generation is the current year.
    base = target_datetime
    base%date%day = 1
    base%date%month = 1
    base%time%hour = 0
    base%time%minute = 0
    base%time%second = 0
    base%time%ms = 0

    ! Round down to the nearest multiple of the file interval.
    CALL divide_time(base, target_datetime, this%file_interval, quot)
    CALL getDatetimeFromJulianDay(shift_time(base, this%file_interval, quot), base)

    filename = generate_filename(this%filename_pattern, base)

    ! No need to reopen the file if the new one is the same.
    IF (.NOT. this%file%is_open() .OR. this%file%filename /= filename) THEN
      this%file = t_sst_sic_file(TRIM(filename), base, this%p_patch%cells%dist_io_data)
    END IF

    IF (.NOT. this%file%is_open()) RETURN

    ntimes = SIZE(this%file%times)

    DO i = 1, ntimes
      IF (this%file%times(i) > jd) EXIT
    END DO

    this%index = i - 1

  END SUBROUTINE sst_sic_reader_goto

  !> Get the current file step's timestamp as julian day.
  FUNCTION sst_sic_reader_get_julian_day (this) RESULT(jd)
    CLASS(t_sst_sic_reader), INTENT(IN) :: this
    TYPE(julianday) :: jd

    jd = this%file%times(this%index)
  END FUNCTION sst_sic_reader_get_julian_day

  !> Get the data associated with the current file step.
  SUBROUTINE sst_sic_reader_read (this, varname, dat)
    CLASS(t_sst_sic_reader), INTENT(INOUT) :: this
    CHARACTER(len=*), INTENT(IN) :: varname
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: dat(:,:,:,:)

    IF (ALLOCATED(dat)) THEN
      IF (ANY(SHAPE(dat) /= [nproma, 1, this%p_patch%nblks_c, 1])) THEN
        DEALLOCATE(dat)
      END IF
    END IF

    IF (.NOT. ALLOCATED(dat)) THEN
      ALLOCATE(dat(nproma, 1, this%p_patch%nblks_c, 1))
      dat(:,:,:,:) = -1._wp
    END IF

    CALL this%file%get_slice(this%index, varname, dat(:,:,:,1))

  END SUBROUTINE sst_sic_reader_read

  !> Move the reader to the next file step.
  !! If the file is exhausted, opens the next file in the sequence.
  !!
  !! If the next file does not exist, the reader is in an invalid state, from which it can recover
  !! by either a `prev` or a `goto` operation.
  SUBROUTINE sst_sic_reader_next (this)
    CLASS(t_sst_sic_reader), INTENT(INOUT) :: this

    TYPE(datetime) :: next_timebase

    IF (this%file%is_open()) THEN
      IF (this%index < SIZE(this%file%times)) THEN
        this%index = this%index + 1
        RETURN
      END IF
    END IF

    next_timebase = this%file%timebase + this%file_interval

    this%file = t_sst_sic_file( &
        & generate_filename(this%filename_pattern, next_timebase), &
        & next_timebase, &
        & this%p_patch%cells%dist_io_data &
      )

    IF (this%file%is_open()) THEN
      this%index = 1
    ELSE
      this%index = 0
    END IF
  END SUBROUTINE sst_sic_reader_next

  !> Move the reader to the previous file step.
  !! If the file is exhausted, opens the previous file in the sequence.
  !!
  !! If the previous file does not exist, the reader is in an invalid state from which it can
  !! recover by either a `next` or a `goto` operation.
  SUBROUTINE sst_sic_reader_prev (this)
    CLASS(t_sst_sic_reader), INTENT(INOUT) :: this

    TYPE(datetime) :: next_timebase
    TYPE(timedelta) :: interval

    IF (this%file%is_open()) THEN
      IF (this%index > 1) THEN
        this%index = this%index - 1
        RETURN
      END IF
    END IF

    ! mtime timedelta does not support `-`.
    interval = this%file_interval
    interval%sign = '-'

    next_timebase = this%file%timebase + interval

    this%file = t_sst_sic_file( &
        & generate_filename(this%filename_pattern, next_timebase), &
        & next_timebase, &
        & this%p_patch%cells%dist_io_data &
      )

    IF (this%file%is_open()) THEN
      this%index = SIZE(this%file%times)
    ELSE
      this%index = 0
    END IF
  END SUBROUTINE sst_sic_reader_prev

  !> Check if the reader is valid, returning the last encountered error if the reader is currently
  !! invalid.
  FUNCTION sst_sic_reader_is_valid (this, msg) RESULT(valid)
    CLASS(t_sst_sic_reader), INTENT(IN) :: this
    CHARACTER(len=*), INTENT(OUT), OPTIONAL :: msg
    LOGICAL :: valid

    valid = this%file%is_open() .AND. this%index > 0
    IF (valid) RETURN

    IF (PRESENT(msg)) THEN
      msg = TRIM(this%file%filename) // ': ' // TRIM(nf90_strerror(this%file%last_netcdf_error))
    END IF

  END FUNCTION sst_sic_reader_is_valid

  FUNCTION sst_sic_reader_get_nblks (this) RESULT(nblks)
    CLASS(t_sst_sic_reader), INTENT(in   ) :: this
    INTEGER                                :: nblks
    nblks = this%p_patch%nblks_c
  END FUNCTION sst_sic_reader_get_nblks

  FUNCTION sst_sic_reader_get_npromz (this) RESULT(npromz)
    CLASS(t_sst_sic_reader), INTENT(in   ) :: this
    INTEGER                                :: npromz
    npromz = this%p_patch%npromz_c
  END FUNCTION sst_sic_reader_get_npromz

  SUBROUTINE sst_sic_reader_deinit (this)
    CLASS(t_sst_sic_reader), INTENT(INOUT) :: this

    CALL this%file%close()

  END SUBROUTINE sst_sic_reader_deinit

  FUNCTION generate_filename (pattern, dt) RESULT(filename)
    CHARACTER(len=*), INTENT(IN) :: pattern
    TYPE(datetime), INTENT(IN) :: dt
    CHARACTER(len=FILENAME_MAX) :: filename

    TYPE(t_keyword_list), POINTER :: kws

    kws => NULL()

    CALL associate_keyword('<year>', TRIM(int2string(INT(dt%date%year), '(i4.4)')), kws)
    CALL associate_keyword('<month>', TRIM(int2string(dt%date%month, '(i2.2)')), kws)
    CALL associate_keyword('<day>', TRIM(int2string(dt%date%day, '(i2.2)')), kws)
    CALL associate_keyword('<hh>', TRIM(int2string(dt%time%hour, '(i2.2)')), kws)
    CALL associate_keyword('<mm>', TRIM(int2string(dt%time%minute, '(i2.2)')), kws)
    CALL associate_keyword('<ss>', TRIM(int2string(dt%time%second, '(i2.2)')), kws)

    filename = with_keywords(kws, TRIM(pattern))

  END FUNCTION generate_filename

END MODULE mo_reader_sst_sic
