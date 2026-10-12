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
!
! Contains the implementation of the boundary conditions for the limited area mode
! for the ocean (ICON-O-LAM).
!
! Some of the more important keywords and concepts:
! * latbc - short for lateral boundary conditions
! * latbc zone - horizontally outermost region of LAM, which is where the boundary
!       conditions are applied; this region has two distinct sub-regions (zones):
! - boundary interpolation zone - outer part of the latbc zone where parent grid
!       data (usually global) is prescribed; this zone does not get modified by
!       solvers (is non-prognostic)
! - nudging zone - inner part of the latbc zone (between boundary interp. zone and
!       inner part of the LAM domain); this region is prognostic, but gets nudged
!       towards parent data to reduce discontinuties that may arise

MODULE mo_ocean_limarea
  !-------------------------------------------------------------------------
  ! USE section
  !-------------------------------------------------------------------------
  USE mo_kind,              ONLY: wp, i8
  USE mo_exception,         ONLY: finish
  USE mo_model_domain,      ONLY: t_patch, t_patch_3d, p_patch
  USE mo_sync,              ONLY: sync_c, sync_e, sync_patch_array
  USE mo_read_interface,    ONLY: read_2D_1time, read_3D_1time, on_cells, on_edges, t_stream_id, &
    & openInputFile, closeFile
  USE mo_grid_subset,       ONLY: t_subset_range, get_index_range
  USE mo_parallel_config,   ONLY: nproma
  USE mo_ocean_nml,         ONLY: n_zlev, vert_cor_type, no_tracer, ocean_latbc_bnd_intp_width, &
    &                             ocean_latbc_nudg_width, ocean_max_refin_c_ctrl, &
    &                             ocean_latbc_itype, ocean_latbc_nudg_func, ocean_latbc_exp_const_A, &
    &                             ocean_latbc_exp_const_b, ocean_latbc_dtime, &
    &                             ocean_latbc_start_datetime_str, ocean_latbc_filepattern, &
    &                             ocean_latbc_path, ocean_latbc_from_vn

  USE mo_io_units,          ONLY: filename_max, nerr
  USE mo_impl_constants,    ONLY: SUCCESS
  USE mtime,                ONLY: timedelta, datetime, getPTStringFromMS, &
    &                             newTimedelta, max_timedelta_str_len, newDatetime, &
    &                             deallocateDatetime, deallocateTimedelta, &
    &                             max_datetime_str_len, getTotalSecondsTimeDelta, &
    &                             OPERATOR(+), OPERATOR(-), OPERATOR(>), OPERATOR(<)
  USE mo_util_mtime,        ONLY: mtime_utils, FMT_DDDHH, FMT_DDHHMMSS, FMT_HHH
  USE mo_impl_constants,    ONLY: MAX_CHAR_LENGTH
  USE mo_util_string,       ONLY: t_keyword_list, int2string, &
    &                             associate_keyword, with_keywords
  USE mo_grid_config,       ONLY: nroot, n_dom
  USE mo_time_config,       ONLY: time_config
  USE mo_mpi,               ONLY: my_process_is_stdio
  USE mo_ocean_types,       ONLY: t_hydro_ocean_state
  USE mo_dynamics_config,   ONLY: nold, nnew
  USE mo_math_types,        ONLY: t_cartesian_coordinates
  USE mo_scalar_product,    ONLY: map_cell2edges_3D
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_math_utilities,    ONLY: gvec2cvec
  USE mo_math_constants,    ONLY: rad2deg
  USE mo_run_config,        ONLY: debug_check_level

  IMPLICIT NONE
  PRIVATE

  !-------------------------------------------------------------------------
  ! PUBLIC declarations
  !-------------------------------------------------------------------------
  ! Attributes
  PUBLIC :: ocean_latbc_data

  ! Methods
  PUBLIC :: init_ocean_latbc
  PUBLIC :: destruct_ocean_latbc
  PUBLIC :: preload_ocean_latbc
  PUBLIC :: apply_ocean_ssh_latbc
  PUBLIC :: apply_ocean_velocity_latbc
  PUBLIC :: apply_ocean_tracer_latbc

  !-------------------------------------------------------------------------
  ! TYPE definitions
  !-------------------------------------------------------------------------
  TYPE t_ocean_latbc_data
    ! the three levels that follow are used for two time levels loaded from
    ! the lateral boundary files plus a third for the interpolated level
    REAL(wp), ALLOCATABLE :: to(:,:,:,:) ! (3, nproma, level, blocks)
    REAL(wp), ALLOCATABLE :: so(:,:,:,:) ! (3, nproma, level, blocks)
    REAL(wp), ALLOCATABLE :: ssh(:,:,:) ! (3, nproma, blocks)
    REAL(wp), ALLOCATABLE :: stretch_c(:,:,:) ! (3, nproma, blocks)
    REAL(wp), ALLOCATABLE :: vn(:,:,:,:) ! (3, nproma, level, blocks)
    REAL(wp), ALLOCATABLE :: u(:,:,:,:) ! (3, nproma, level, blocks)
    REAL(wp), ALLOCATABLE :: v(:,:,:,:) ! (3, nproma, level, blocks)

    REAL(wp), ALLOCATABLE :: nudge_cell(:,:) ! (nproma, blocks)
    REAL(wp), ALLOCATABLE :: nudge_edge(:,:) ! (nproma, blocks)
  CONTAINS
    PROCEDURE :: alloc => t_ocean_latbc_data_allocate
    PROCEDURE :: finalize => t_ocean_latbc_data_finalize
  END TYPE t_ocean_latbc_data

  !-------------------------------------------------------------------------
  ! MODULE PARAMETERS
  !-------------------------------------------------------------------------
  CHARACTER(LEN=16), PARAMETER :: module_name = 'mo_ocean_limarea'

  ! CONSTANTS (for better readability):
  INTEGER, PARAMETER :: OCEAN_LATBC_TYPE_CONST       = 0 ! constant lateral boundary conditions
  INTEGER, PARAMETER :: OCEAN_LATBC_TYPE_VAR         = 1 ! time-dependent lateral boundary conditions

  !-------------------------------------------------------------------------
  ! MODULE VARIABLES
  !-------------------------------------------------------------------------
  TYPE(t_ocean_latbc_data), TARGET    :: ocean_latbc_data
  INTEGER                             :: ocean_latbc_timelevel ! timelevel of younger latbc. the other is calculated on the fly
  TYPE(datetime), POINTER             :: ocean_latbc_start_datetime ! starting date for the latbc
  TYPE(timedelta), POINTER            :: ocean_latbc_dtime_mtime ! dt between latbc levels
  CHARACTER(LEN=filename_max)         :: ocean_latbc_filename ! filename read by the module methods
  TYPE(datetime)                      :: ocean_latbc_prev_datetime ! mtime date of older latbc data
  TYPE(datetime)                      :: ocean_latbc_next_datetime ! mtime date of newer latbc data
  CHARACTER(LEN=max_datetime_str_len) :: ocean_latbc_next_datetime_str ! string version
  INTEGER                             :: ocean_latbc_desync ! 0 - both latbc OK, 1 - old latbc out of date, 2+ - both latbcs out of date

CONTAINS
  !-------------------------------------------------------------------------
  SUBROUTINE init_ocean_latbc(patch_3d, operators_coeff)
    TYPE(t_patch_3d), POINTER, INTENT(in) :: patch_3d
    TYPE(t_operator_coeff), INTENT(in), TARGET :: operators_coeff

    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: numCellBlocks, numEdgeBlocks
    REAL(wp) :: ocean_latbc_dtime_in_ms
    CHARACTER(LEN=max_timedelta_str_len) :: ocean_latbc_dtime_str
    INTEGER :: errno
    INTEGER :: start_index_c, end_index_c
    INTEGER :: start_index_e, end_index_e
    INTEGER :: jc, je, blockNo
    CHARACTER(*), PARAMETER :: method_name = module_name//"::init_ocean_latbc"
    TYPE(timedelta), POINTER :: ocean_time_step
    TYPE(datetime), POINTER :: ocean_start_date

    991 FORMAT(a,a)

    patch_2d => patch_3d%p_patch_2d(1)
    ocean_start_date => time_config%tc_exp_startdate
    ocean_time_step => time_config%tc_dt_model

    ! allocate the arrays for latbc data
    numCellBlocks = patch_2d%alloc_cell_blocks
    numEdgeBlocks = patch_2d%nblks_e
    CALL ocean_latbc_data%alloc(numCellBlocks, numEdgeBlocks)

    CALL precalculate_latbc_nudging_coefficients()

    ! convert latbc start date from string to mtime object (datetime)
    IF (TRIM(ocean_latbc_start_datetime_str) /= '') THEN
      ocean_latbc_start_datetime => newDatetime(TRIM(ocean_latbc_start_datetime_str), errno=errno)
      IF (errno /= SUCCESS) CALL finish(method_name, &
        & "LAM: Error allocating ocean datetime object (ocean_latbc_start_datetime_str is probably malformed)")
    ELSE
      ocean_latbc_start_datetime => newDatetime(ocean_start_date, errno=errno) ! default setting
      IF (errno /= SUCCESS) CALL finish(method_name, "LAM: Error allocating ocean datetime object")
    END IF

    IF (ocean_latbc_itype == OCEAN_LATBC_TYPE_CONST) THEN ! Preload latbcs into 3rd level
      ocean_latbc_filename = TRIM(ocean_latbc_path) &
        & // TRIM(generate_filename_ocean(ocean_latbc_start_datetime))
      IF (debug_check_level > 4 .AND. my_process_is_stdio()) &
        & WRITE(nerr,991) " reading next latbc file: ", TRIM(ocean_latbc_filename)
      CALL read_latbc_from_file(patch_3d, 3, operators_coeff)
      RETURN ! Nothing that follows is needed for constant latbcs
    END IF

    ! convert ocean_latbc_dtime (integer) into mtime object (timedelta)
    ocean_latbc_dtime_in_ms = 1000._wp * ocean_latbc_dtime
    CALL getPTStringFromMS(NINT(ocean_latbc_dtime_in_ms,i8), ocean_latbc_dtime_str)
    ocean_latbc_dtime_mtime => newTimedelta(ocean_latbc_dtime_str, errno)
    IF (errno /= SUCCESS) CALL finish(method_name, "LAM: Error in initialization of ocean_latbc_dtime_mtime")

    ! set the datetime of initial latbc files (no support for IAU yet)
    ocean_latbc_prev_datetime = ocean_latbc_start_datetime
    ocean_latbc_next_datetime = ocean_latbc_start_datetime + ocean_latbc_dtime_mtime
    IF (ocean_latbc_prev_datetime > ocean_start_date + ocean_time_step .OR. &
      & ocean_latbc_next_datetime < ocean_start_date) &
      & CALL finish(method_name, "LAM: The latbc starting date is not valid")

    ocean_latbc_desync = 2 ! mark full desync in order to load both latbcs on next preload phase
    ocean_latbc_timelevel = 1 ! mark next latbc level

    ! mapping for sparse grid would also be done here if it was there

  CONTAINS

    !-------------------------------------------------------------------------
    SUBROUTINE precalculate_latbc_nudging_coefficients()
      INTEGER, PARAMETER   :: discrete_parabolic_sigmoid = 1 ! it looks similar to sigmoid
      INTEGER, PARAMETER   :: exponential_decay = 2

      INTEGER :: idx ! iteration index for cell/edge nudging rows
      INTEGER :: N   ! total number of cell/edge nudging rows
      INTEGER :: jc, je, blockNo, etype
      REAL(wp) :: total_sum
      REAL(wp), DIMENSION(2*ocean_latbc_nudg_width,2) :: nudging_func
      INTEGER,  DIMENSION(2*ocean_latbc_nudg_width) :: deltas, partial_sum
      TYPE(t_subset_range), POINTER :: all_cells, all_edges

      nudging_func = 0.0_wp
      SELECT CASE (ocean_latbc_nudg_func)
      CASE (discrete_parabolic_sigmoid)
        ! this is a discrete approximation of: convex parabolic rise followed by C1 continuous
        ! concave parabolic rise forming an S-shaped smooth step function. the main building block
        ! is the discrete derivative (deltas  which goes eg. 1 2 3 4 3 2 1 or 1 2 3 3 2 1)
        DO etype = 1, 2 ! element type, 1=cells, 2=edges
          N = etype * ocean_latbc_nudg_width ! lucky choice of indices, as edges require 2x the cells
          total_sum = ((N+2)/2)*((N+3)/2)
          deltas = 1
          DO idx = 2, (N+1)/2
            deltas(idx) = deltas(idx-1) + 1
          END DO
          deltas((N+1)/2+1) = N/2+1
          DO idx = N/2+2, N
            deltas(idx) = deltas(idx-1) - 1
          END DO
          partial_sum = 1 ! prevent other indices from having NaNs when not using all indices
          nudging_func(1,etype) = 1 - partial_sum(1) / total_sum
          DO idx = 2, N
            partial_sum(idx) = partial_sum(idx-1) + deltas(idx)
            nudging_func(idx,etype) = 1 - partial_sum(idx) / total_sum
          END DO
        END DO
      CASE (exponential_decay)
        ! this is an exponential decay function y=A*exp(-b*x)
        ! cells
        DO idx = 1, ocean_latbc_nudg_width
          nudging_func(idx,1) = ocean_latbc_exp_const_A*exp(-ocean_latbc_exp_const_b*idx)
        END DO
        ! edges
        DO idx = 1, 2*ocean_latbc_nudg_width
          nudging_func(idx,2) = ocean_latbc_exp_const_A*exp(-ocean_latbc_exp_const_b/2*idx)
        END DO
      CASE DEFAULT ! should never happen due to the namelist crosscheck
        CALL finish(method_name, "LAM: invalid value for ocean_latbc_nudg_func")
      END SELECT

      ! for cells
      N = ocean_latbc_nudg_width
      all_cells => patch_2d%cells%all
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_index_c, end_index_c)
        DO jc = start_index_c, end_index_c ! do it later using ordering
          idx = patch_2d%cells%refin_ctrl(jc,blockNo) - ocean_latbc_bnd_intp_width
          IF (idx > 0 .AND. idx <= N) THEN
            ocean_latbc_data%nudge_cell(jc,blockNo) = nudging_func(idx,1)
          END IF
        END DO ! jc
      END DO ! blockNo

      ! for edges
      N = 2*ocean_latbc_nudg_width
      all_edges => patch_2d%edges%all
      DO blockNo = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, blockNo, start_index_e, end_index_e)
        DO je = start_index_e, end_index_e ! do it later using oredering
          idx = patch_2d%edges%refin_ctrl(je,blockNo) - 2*ocean_latbc_bnd_intp_width
          IF (idx > 0 .AND. idx <= N) THEN
            ocean_latbc_data%nudge_edge(je,blockNo) = nudging_func(idx,2)
          END IF
        END DO ! je
      END DO ! blockNo
    END SUBROUTINE precalculate_latbc_nudging_coefficients

  END SUBROUTINE init_ocean_latbc

  !-------------------------------------------------------------------------
  SUBROUTINE destruct_ocean_latbc()
    CALL ocean_latbc_data%finalize
    CALL deallocateDatetime(ocean_latbc_start_datetime)
    IF (ocean_latbc_itype /= OCEAN_LATBC_TYPE_CONST) CALL deallocateTimedelta(ocean_latbc_dtime_mtime)
  END SUBROUTINE destruct_ocean_latbc

  !-------------------------------------------------------------------------
  ! This function takes the time of the lateral boundary to be loaded from
  ! the caller, together with a user-provided template (through namelist),
  ! and generates the filename of the lat. boundary conditions to be read.
  ! Of course, this only works if the data has been properly preprocessed,
  ! which is here assumed to be the case (each time step in a separate file
  ! that includes that time in the filename in some way, shape or form).
  ! For example, user would enter (in ocean_limarea_nml):
  ! ocean_latbc_filepattern = 'LAM_east_pacific_r2b8_<y><m><d>T<h><min><sec>Z.nc'
  ! and this function would produce something like:
  ! result_str = 'LAM_east_pacific_r2b8_20200831T010500Z.nc'
  ! More than just year, month, day, hour, minute and second can be used.
  ! For full reference, check the function body.
  FUNCTION generate_filename_ocean(latbc_mtime, opt_mtime_begin) RESULT(result_str)
    CHARACTER(MAX_CHAR_LENGTH)                       :: result_str
    TYPE(datetime),   INTENT(IN)                     :: latbc_mtime
    ! Optional: Start date, which a time span in the filename is related to.
    TYPE(datetime),   INTENT(IN),  POINTER, OPTIONAL :: opt_mtime_begin
    ! Local variables
    CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: method_name = module_name//'::generate_filename_ocean'
    TYPE (t_keyword_list), POINTER        :: keywords => NULL()
    CHARACTER(MAX_CHAR_LENGTH)            :: str
    INTEGER :: jlev

    jlev = p_patch(1)%level ! not tested yet - does it work in the ocean the same as in the atmo

    WRITE(str,'(i4)')   latbc_mtime%date%year
    CALL associate_keyword("<y>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%date%month
    CALL associate_keyword("<m>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%date%day
    CALL associate_keyword("<d>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%hour
    CALL associate_keyword("<h>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%minute
    CALL associate_keyword("<min>",       TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%second !FLOOR(latbc_mtime%time%second)
    CALL associate_keyword("<sec>",       TRIM(str),                        keywords)

    CALL associate_keyword("<nroot>",     TRIM(int2string(nroot,'(i1)')),   keywords)
    CALL associate_keyword("<nroot0>",    TRIM(int2string(nroot,'(i2.2)')), keywords)
    CALL associate_keyword("<jlev>",      TRIM(int2string(jlev, '(i2.2)')), keywords)
    CALL associate_keyword("<dom>",       TRIM(int2string(1,'(i2.2)')),     keywords)

    IF (PRESENT(opt_mtime_begin)) THEN
      CALL associate_keyword("<ddhhmmss>", &
        &                    TRIM(mtime_utils%ddhhmmss(opt_mtime_begin, latbc_mtime, FMT_DDHHMMSS)), &
        &                    keywords)
      CALL associate_keyword("<dddhh>",    &
        &                    TRIM(mtime_utils%ddhhmmss(opt_mtime_begin, latbc_mtime, FMT_DDDHH)),    &
        &                    keywords)
      CALL associate_keyword("<hhh>",    &
        &                    TRIM(mtime_utils%ddhhmmss(opt_mtime_begin, latbc_mtime, FMT_HHH)),    &
        &                    keywords)
    END IF

    ! replace keywords in latbc_filename
    result_str = TRIM(with_keywords(keywords, TRIM(ocean_latbc_filepattern)))
  END FUNCTION generate_filename_ocean

  !-------------------------------------------------------------------------
  SUBROUTINE t_ocean_latbc_data_allocate(this, numCellBlocks, numEdgeBlocks)
    CLASS(t_ocean_latbc_data) :: this
    INTEGER, INTENT(in) :: numCellBlocks
    INTEGER, INTENT(in) :: numEdgeBlocks
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//'::allocate_ocean_latbc'
    INTEGER :: ierrstat=0

    IF (no_tracer >= 1) THEN
      ALLOCATE(this%to(3,nproma,n_zlev,numCellBlocks))
      this%to = 0.0_wp
      IF (no_tracer >= 2) THEN
        ALLOCATE(this%so(3,nproma,n_zlev,numCellBlocks))
        this%so = 0.0_wp
      END IF
    END IF
    ALLOCATE(this%ssh(3,nproma,numCellBlocks))
    ALLOCATE(this%stretch_c(3,nproma,numCellBlocks))
    this%ssh = 0.0_wp
    this%stretch_c = 0.0_wp
    IF (.NOT. ocean_latbc_from_vn) THEN
      ALLOCATE(this%u(3,nproma,n_zlev,numCellBlocks))
      ALLOCATE(this%v(3,nproma,n_zlev,numCellBlocks))
      this%u = 0.0_wp
      this%v = 0.0_wp
    END IF
    ALLOCATE(this%vn(3,nproma,n_zlev,numEdgeBlocks))
    ALLOCATE(this%nudge_cell(nproma,numCellBlocks))
    ALLOCATE(this%nudge_edge(nproma,numEdgeBlocks), STAT=ierrstat)
    this%vn = 0.0_wp
    this%nudge_cell = 0.0_wp
    this%nudge_edge = 0.0_wp

    IF (ierrstat /= SUCCESS) CALL finish(method_name, "LAM: an ALLOCATE failed!")
  END SUBROUTINE t_ocean_latbc_data_allocate

  !-------------------------------------------------------------------------
  SUBROUTINE t_ocean_latbc_data_finalize(this)
    CLASS(t_ocean_latbc_data) :: this
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//'::t_ocean_latbc_data_finalize'
    INTEGER :: ierrstat=0

    IF (no_tracer >= 1) THEN
      IF (ALLOCATED(this%to)) DEALLOCATE(this%to, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(method_name, "LAM: DEALLOCATE to failed!")
      IF (no_tracer >= 2) THEN
        IF (ALLOCATED(this%so)) DEALLOCATE(this%so, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(method_name, "LAM: DEALLOCATE so failed!")
      END IF
    END IF
    IF (ALLOCATED(this%ssh)) DEALLOCATE(this%ssh, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "LAM: DEALLOCATE ssh failed!")
    IF (ALLOCATED(this%stretch_c)) DEALLOCATE(this%stretch_c, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "LAM: DEALLOCATE stretch_c failed!")
    IF (.NOT. ocean_latbc_from_vn) THEN
      IF (ALLOCATED(this%u)) DEALLOCATE(this%u, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(method_name, "LAM: DEALLOCATE u failed!")
      IF (ALLOCATED(this%v)) DEALLOCATE(this%v, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(method_name, "LAM: DEALLOCATE v failed!")
    END IF
    IF (ALLOCATED(this%vn)) DEALLOCATE(this%vn, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "LAM: DEALLOCATE vn failed!")

    IF (ALLOCATED(this%nudge_cell)) DEALLOCATE(this%nudge_cell, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "LAM: DEALLOCATE nudge_cell failed!")
    IF (ALLOCATED(this%nudge_edge)) DEALLOCATE(this%nudge_edge, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "LAM: DEALLOCATE nudge_edge failed!")
  END SUBROUTINE t_ocean_latbc_data_finalize

  !-------------------------------------------------------------------------
  SUBROUTINE preload_ocean_latbc(patch_3d, current_time, operators_coeff, jstep)
    ! current_time is the time at the end of the timestep, which is n+1

    ! we never need to fix the initial state, since we believe it to be correct (allowing inconsistencies
    ! between initial state and initial boundary conditions is not something we need to consider).
    ! the main reason we even need to differentiate between starting points of experiment and latbc,
    ! is if we use IAU, for which we need to have states from the past (I assume)

    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3d
    TYPE(datetime), INTENT(in) :: current_time ! make it pointer?
    TYPE(t_operator_coeff), INTENT(in), TARGET :: operators_coeff
    INTEGER, INTENT(in) :: jstep

    TYPE(t_patch), POINTER :: patch_2d
    TYPE(timedelta) :: current_ocean_latbc_dtime_mtime
    REAL(wp) :: current_ocean_latbc_dtime
    REAL(wp) :: latbc_f ! interpolation factor
    INTEGER :: jg, l1, l2

    991 FORMAT(a,a)
    992 FORMAT(a,1pg26.18)

    IF (ocean_latbc_itype == OCEAN_LATBC_TYPE_CONST) RETURN ! Nothing to preload

    jg = n_dom ! no support for nested grids!
    patch_2d => patch_3d%p_patch_2d(1)

    ! make sure to find the right pair of prev and next latbc before the interpolation
    DO WHILE (ocean_latbc_next_datetime < current_time)
      ocean_latbc_prev_datetime = ocean_latbc_next_datetime
      ocean_latbc_next_datetime = ocean_latbc_next_datetime + ocean_latbc_dtime_mtime
      ocean_latbc_desync = ocean_latbc_desync + 1
    END DO

    ! If ocean_latbc_desync is greater than 1, newer latbc needs to be replaced
    IF (ocean_latbc_desync > 1) THEN
      ocean_latbc_filename = TRIM(ocean_latbc_path) &
        & // TRIM(generate_filename_ocean(ocean_latbc_prev_datetime))
      IF (debug_check_level > 4 .AND. my_process_is_stdio()) &
        & WRITE(nerr,991) " reading next latbc file: ", TRIM(ocean_latbc_filename)
      CALL read_latbc_from_file(patch_3d, ocean_latbc_timelevel, operators_coeff)
    END IF

    ! If ocean_latbc_desync is greater than 0, old latbc needs to be replaced
    IF (ocean_latbc_desync > 0) THEN
      ocean_latbc_filename = TRIM(ocean_latbc_path) &
        & // TRIM(generate_filename_ocean(ocean_latbc_next_datetime))
      ocean_latbc_timelevel = 3 - ocean_latbc_timelevel ! jump to old level, it's the new one now
      IF (debug_check_level > 4 .AND. my_process_is_stdio()) &
        & WRITE(nerr,991) " reading next latbc file: ", TRIM(ocean_latbc_filename)
      CALL read_latbc_from_file(patch_3d, ocean_latbc_timelevel, operators_coeff) ! overwrite old
    END IF
    ocean_latbc_desync = 0 ! both latbcs are now valid

    ! calculate the interpolation coefficient
    current_ocean_latbc_dtime_mtime = current_time - ocean_latbc_prev_datetime
    current_ocean_latbc_dtime = REAL(getTotalSecondsTimeDelta(current_ocean_latbc_dtime_mtime, &
      &                                                       ocean_latbc_prev_datetime))
    latbc_f = current_ocean_latbc_dtime / ocean_latbc_dtime

    IF (debug_check_level > 6 .AND. my_process_is_stdio()) WRITE(nerr,992) " latbc_f = ", latbc_f

    ! interpolate
    l1 = 3 - ocean_latbc_timelevel
    l2 = ocean_latbc_timelevel
    IF (no_tracer >= 1) THEN
      ocean_latbc_data%to(3,:,:,:) = &
      & (1-latbc_f)*ocean_latbc_data%to(l1,:,:,:) + &
      &    latbc_f *ocean_latbc_data%to(l2,:,:,:) ! 1 for temperature
      IF (no_tracer >= 2) THEN
        ocean_latbc_data%so(3,:,:,:) = &
        & (1-latbc_f)*ocean_latbc_data%so(l1,:,:,:) + &
        &    latbc_f *ocean_latbc_data%so(l2,:,:,:) ! 2 for salinity
      END IF
    END IF
    ocean_latbc_data%ssh(3,:,:) = &
    & (1-latbc_f)*ocean_latbc_data%ssh(l1,:,:) + &
    &    latbc_f *ocean_latbc_data%ssh(l2,:,:)
    ocean_latbc_data%vn(3,:,:,:) = &
    & (1-latbc_f)*ocean_latbc_data%vn(l1,:,:,:) + &
    &    latbc_f *ocean_latbc_data%vn(l2,:,:,:)

  ! At this point, apply routines are ready to be applied as needed

  END SUBROUTINE preload_ocean_latbc

  !-------------------------------------------------------------------------
  SUBROUTINE apply_ocean_ssh_latbc(patch_3d, ocean_state, mt)

    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    INTEGER :: mt ! model timelevel, either nold(jg) or nnew(jg)

    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: blockNo, jc, jg
    INTEGER :: start_index_c, end_index_c
    REAL(wp) :: latbc_nudge
    CHARACTER(*), PARAMETER :: method_name = module_name//"::apply_ocean_ssh_latbc"
    REAL(wp), POINTER, DIMENSION(:,:) :: ssh => NULL()

    jg = n_dom ! no support for nested grids!
    patch_2d => patch_3d%p_patch_2d(1)

    IF (vert_cor_type == 0) THEN ! z coord system
      CALL finish(method_name, "LAM: z coordinate system not yet supported")
      !ssh => ocean_state%p_prog(mt)%h
    ELSE ! IF (vert_cor_type == 1) THEN ! z* coord system
      ssh => ocean_state%p_prog(mt)%eta_c
    END IF

    ! assign to proper cells
    all_cells => patch_2d%cells%all
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index_c, end_index_c)
      DO jc = start_index_c, end_index_c ! do it later using ordering
        IF (patch_2d%cells%refin_ctrl(jc,blockNo) > 0 .AND. &
          & patch_2d%cells%refin_ctrl(jc,blockNo) <= ocean_latbc_bnd_intp_width) THEN
          ssh(jc,blockNo) = ocean_latbc_data%ssh(3,jc,blockNo)
        ELSE IF (patch_2d%cells%refin_ctrl(jc,blockNo) > ocean_latbc_bnd_intp_width .AND. &
          & patch_2d%cells%refin_ctrl(jc,blockNo) <= ocean_latbc_bnd_intp_width + ocean_latbc_nudg_width) THEN ! nudge cells
          latbc_nudge = ocean_latbc_data%nudge_cell(jc,blockNo)
          ssh(jc,blockNo) = &
          & (1-latbc_nudge)* ssh(jc,blockNo) + &
          &    latbc_nudge * ocean_latbc_data%ssh(3,jc,blockNo) ! 3 for interpolated layer
        END IF
      END DO ! jc
    END DO ! blockNo

  END SUBROUTINE apply_ocean_ssh_latbc

  !-------------------------------------------------------------------------
  SUBROUTINE apply_ocean_velocity_latbc(patch_3d, ocean_state, mt)

    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    INTEGER :: mt ! model timelevel, either nold(jg) or nnew(jg)

    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges
    INTEGER :: blockNo, je, jk, jg
    INTEGER :: start_index_e, end_index_e
    INTEGER :: start_level, end_level
    REAL(wp) :: latbc_nudge
    CHARACTER(*), PARAMETER :: method_name = module_name//"::apply_ocean_velocity_latbc"

    jg = n_dom ! no support for nested grids!
    patch_2d => patch_3d%p_patch_2d(1)
    start_level = 1

    ! assign to proper edges
    all_edges => patch_2d%edges%all
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_index_e, end_index_e)
      DO je = start_index_e, end_index_e ! do it later using oredering
        end_level = patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
        IF (patch_2d%edges%refin_ctrl(je,blockNo) > 0 .AND. &
          & patch_2d%edges%refin_ctrl(je,blockNo) <= 2*ocean_latbc_bnd_intp_width) THEN
          DO jk = start_level, end_level
            ocean_state%p_prog(mt)%vn(je,jk,blockNo) = &
            & ocean_latbc_data%vn(3,je,jk,blockNo) ! 3 for interpolated layer
          END DO
        ELSE IF (patch_2d%edges%refin_ctrl(je,blockNo) > 2*ocean_latbc_bnd_intp_width .AND. &
          & patch_2d%edges%refin_ctrl(je,blockNo) <= 2*(ocean_latbc_bnd_intp_width + ocean_latbc_nudg_width)) THEN ! nudge edges
          latbc_nudge = ocean_latbc_data%nudge_edge(je,blockNo)
          DO jk = start_level, end_level
            ocean_state%p_prog(mt)%vn(je,jk,blockNo) = &
            & (1-latbc_nudge)* ocean_state%p_prog(mt)%vn(je,jk,blockNo) + &
            &    latbc_nudge * ocean_latbc_data%vn(3,je,jk,blockNo) ! 3 for interpolated layer
          END DO
        END IF
      END DO ! je
    END DO ! blockNo

  END SUBROUTINE apply_ocean_velocity_latbc

  !-------------------------------------------------------------------------
  SUBROUTINE apply_ocean_tracer_latbc(patch_3d, ocean_state, mt)

    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    INTEGER, INTENT(in) :: mt ! model timelevel, either nold(jg) or nnew(jg)

    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: blockNo, jc, jk, jg
    INTEGER :: start_index_c, end_index_c
    INTEGER :: start_level, end_level
    REAL(wp) :: latbc_nudge
    CHARACTER(*), PARAMETER :: method_name = module_name//"::apply_ocean_tracer_latbc"

    jg = n_dom ! no support for nested grids!
    patch_2d => patch_3d%p_patch_2d(1)
    start_level = 1

    ! assign to proper cells
    all_cells => patch_2d%cells%all
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index_c, end_index_c)
      DO jc = start_index_c, end_index_c ! do it later using ordering
        end_level = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
        IF (patch_2d%cells%refin_ctrl(jc,blockNo) > 0 .AND. &
          & patch_2d%cells%refin_ctrl(jc,blockNo) <= ocean_latbc_bnd_intp_width) THEN
          IF (no_tracer >= 1) THEN
            DO jk = start_level, end_level
              ocean_state%p_prog(mt)%tracer(jc,jk,blockNo,1) = & ! 1 for temperature
              & ocean_latbc_data%to(3,jc,jk,blockNo) ! 3 for interpolated layer
            END DO
            IF (no_tracer >= 2) THEN
              DO jk = start_level, end_level
                ocean_state%p_prog(mt)%tracer(jc,jk,blockNo,2) = & ! 2 for salinity
                & ocean_latbc_data%so(3,jc,jk,blockNo) ! 3 for interpolated layer
              END DO
            END IF
          END IF
        ELSE IF (patch_2d%cells%refin_ctrl(jc,blockNo) > ocean_latbc_bnd_intp_width .AND. &
          & patch_2d%cells%refin_ctrl(jc,blockNo) <= ocean_latbc_bnd_intp_width + ocean_latbc_nudg_width) THEN ! nudge cells
          latbc_nudge = ocean_latbc_data%nudge_cell(jc,blockNo)
          IF (no_tracer >= 1) THEN
            DO jk = start_level, end_level
              ocean_state%p_prog(mt)%tracer(jc,jk,blockNo,1) = & ! 1 for temperature
              & (1-latbc_nudge)* ocean_state%p_prog(mt)%tracer(jc,jk,blockNo,1) + &
              &    latbc_nudge * ocean_latbc_data%to(3,jc,jk,blockNo) ! 3 for interpolated layer
            END DO
            IF (no_tracer >= 2) THEN
              DO jk = start_level, end_level
                ocean_state%p_prog(mt)%tracer(jc,jk,blockNo,2) = & ! 2 for salinity
                & (1-latbc_nudge)* ocean_state%p_prog(mt)%tracer(jc,jk,blockNo,2) + &
                &    latbc_nudge * ocean_latbc_data%so(3,jc,jk,blockNo) ! 3 for interpolated layer
              END DO
            END IF
          END IF
        END IF
      END DO ! jc
    END DO ! blockNo

  END SUBROUTINE apply_ocean_tracer_latbc

  !-------------------------------------------------------------------------
  SUBROUTINE read_latbc_from_file(patch_3d, timelevel, operators_coeff)
    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3d
    INTEGER, INTENT(in) :: timelevel
    TYPE(t_operator_coeff), INTENT(in), TARGET :: operators_coeff

    LOGICAL :: has_missValue(6) = (/ .false., .false., .false., .false., .false., .false. /)
    REAL(wp) :: missValue(6) = (/ -99999999.0_wp, -99999999.0_wp, -99999999.0_wp, -99999999.0_wp, -99999999.0_wp, -99999999.0_wp /)
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: varNames(6) = (/ "to ", "so ", "zos", "u  ", "v  ", "vn " /)

    TYPE(t_patch),POINTER :: patch_2d
    TYPE(t_stream_id) :: stream_id
    INTEGER :: blockNo, idx, level
    INTEGER :: start_cell_index, end_cell_index, start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_cartesian_coordinates), ALLOCATABLE ::cellVelocity_cc(:,:,:)
    CHARACTER(*), PARAMETER :: method_name = module_name//':read_latbc_from_file '
    REAL(wp) :: lat, lon
    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    CALL openInputFile(stream_id, TRIM(ocean_latbc_filename), patch_2d)

    IF (no_tracer >= 1) THEN
      CALL read_3D_1time(stream_id=stream_id, location=on_cells, &
        & variable_name=TRIM(varNames(1)), fill_array=ocean_latbc_data%to(timelevel,:,:,:), &
        & has_missValue=has_missValue(1),  missValue=missValue(1))
      IF (no_tracer >= 2) THEN
        CALL read_3D_1time(stream_id=stream_id, location=on_cells, &
          & variable_name=TRIM(varNames(2)), fill_array=ocean_latbc_data%so(timelevel,:,:,:), &
          & has_missValue=has_missValue(2),  missValue=missValue(2))
      END IF
    END IF
    CALL read_2D_1time(stream_id=stream_id, location=on_cells, &
      & variable_name=TRIM(varNames(3)), fill_array=ocean_latbc_data%ssh(timelevel,:,:), &
      & has_missValue=has_missValue(3),  missValue=missValue(3))
    IF (ocean_latbc_from_vn) THEN
      CALL read_3D_1time(stream_id=stream_id, location=on_edges, &
        & variable_name=TRIM(varNames(6)), fill_array=ocean_latbc_data%vn(timelevel,:,:,:), &
        & has_missValue=has_missValue(6),  missValue=missValue(6))
    ELSE
      CALL read_3D_1time(stream_id=stream_id, location=on_cells, &
        & variable_name=TRIM(varNames(4)), fill_array=ocean_latbc_data%u(timelevel,:,:,:), &
        & has_missValue=has_missValue(4),  missValue=missValue(4))
      CALL read_3D_1time(stream_id=stream_id, location=on_cells, &
        & variable_name=TRIM(varNames(5)), fill_array=ocean_latbc_data%v(timelevel,:,:,:), &
        & has_missValue=has_missValue(5),  missValue=missValue(5))
    ENDIF

    CALL closeFile(stream_id)

    IF (ocean_latbc_from_vn) THEN
      CONTINUE
      !CALL sync_patch_array(sync_e, patch_2d, ocean_latbc_data%vn(timelevel,:,:,:), lacc=.FALSE.)
    ELSE
      ALLOCATE(cellVelocity_cc(nproma,n_zlev,patch_2d%alloc_cell_blocks))

      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
        DO idx = start_cell_index, end_cell_index
          lon = patch_2d%cells%center(idx,blockNo)%lon
          lat = patch_2d%cells%center(idx,blockNo)%lat

          DO level = 1, n_zlev
            CALL gvec2cvec (ocean_latbc_data%u(timelevel, idx, level, blockNo), &
                    & ocean_latbc_data%v(timelevel, idx, level, blockNo), &
                    & lon, lat,                         &
                    & cellVelocity_cc(idx, level, blockNo)%x(1), &
                    & cellVelocity_cc(idx, level, blockNo)%x(2), &
                    & cellVelocity_cc(idx, level, blockNo)%x(3), &
                    & patch_2d%geometry_info)
          END DO
        END DO
      END DO

      CALL map_cell2edges_3D(patch_3d, cellVelocity_cc, ocean_latbc_data%vn(timelevel,:,:,:), operators_coeff)
      !CALL sync_patch_array(sync_e, patch_2d, ocean_latbc_data%vn(timelevel,:,:,:), lacc=.FALSE.)
      DEALLOCATE(cellVelocity_cc)
    ENDIF

    !CALL sync_patch_array(sync_c, patch_2d, ocean_latbc_data%to(timelevel,:,:,:), lacc=.FALSE.)
    !CALL sync_patch_array(sync_c, patch_2d, ocean_latbc_data%so(timelevel,:,:,:), lacc=.FALSE.)
    !CALL sync_patch_array(sync_c, patch_2d, ocean_latbc_data%ssh(timelevel,:,:), lacc=.FALSE.)

  END SUBROUTINE read_latbc_from_file
  !-------------------------------------------------------------------------
END MODULE mo_ocean_limarea
