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

! Reading of forcing data in standalone mode
!
! Routines for reading forcing data at user-defined time intervals in
! standalone mode.

MODULE mo_wave_forcing

  USE mo_kind,                     ONLY: wp
  USE mo_exception,                ONLY: finish, warning, message, message_text
  USE mo_impl_constants,           ONLY: SUCCESS
  USE mo_io_units,                 ONLY: filename_max
  USE mo_model_domain,             ONLY: t_patch
  USE mo_master_config,            ONLY: getModelBaseDir
  USE mo_grid_config,              ONLY: n_dom, nroot
  USE mo_run_config,               ONLY: msg_level
  USE mo_wave_config,              ONLY: t_wave_config, generate_filename
  USE mo_reader_sst_sic,           ONLY: t_sst_sic_reader
  USE mo_interpolate_time,         ONLY: t_time_intp_transient
  USE mo_fortran_tools,            ONLY: copy, DO_DEALLOCATE
  USE mtime,                       ONLY: datetime, datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_wave_td_update,           ONLY: update_ice_free_mask, &
    &                                    update_speed_and_direction
  USE mo_mpi,                      ONLY: my_process_is_mpi_workroot, p_io, p_bcast, &
    &                                    p_comm_work

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_forcing'

  ! wave-specific forcing data type
  !
  TYPE :: t_read_wave_forcing

    TYPE(t_patch), POINTER :: p_patch => NULL()  ! pointer to actual patch
    ! time interpolated fields
    !
    REAL(wp), ALLOCATABLE :: u10m_raw (:,:,:,:)  ! zonal wind in 10m (m/s)
    REAL(wp), ALLOCATABLE :: v10m_raw (:,:,:,:)  ! meridional wind in 10m (m/s)
    REAL(wp), ALLOCATABLE :: sic_raw  (:,:,:,:)  ! sea ice concentration at centers (1)
    REAL(wp), ALLOCATABLE :: slh_raw  (:,:,:,:)  ! sea level height at centers ( m )
    REAL(wp), ALLOCATABLE :: uosc_raw (:,:,:,:)  ! zonal ocean surface current at centers ( m/s )
    REAL(wp), ALLOCATABLE :: vosc_raw (:,:,:,:)  ! meridional ocean surface current at centers ( m/s )

    ! reader
    TYPE(t_sst_sic_reader) :: u10_reader, v10_reader
    TYPE(t_sst_sic_reader) :: sic_reader
    TYPE(t_sst_sic_reader) :: slh_reader
    TYPE(t_sst_sic_reader) :: uosc_reader, vosc_reader

    ! time interpolator
    TYPE(t_time_intp_transient) :: u10_intp, v10_intp
    TYPE(t_time_intp_transient) :: sic_intp
    TYPE(t_time_intp_transient) :: slh_intp
    TYPE(t_time_intp_transient) :: uosc_intp, vosc_intp

    ! forcing flags
    LOGICAL :: l_wind_exist, l_ice_exist, l_slh_exist, l_osc_exist

    ! initialization flag
    LOGICAL :: isInit = .FALSE.
  CONTAINS

   PROCEDURE :: init               => read_wave_forcing__init
   PROCEDURE :: deinit             => read_wave_forcing__deinit
   PROCEDURE :: get_new_rawdata    => read_wave_forcing__get_new_rawdata
   PROCEDURE :: update_forcing     => read_wave_forcing__update_forcing
  END TYPE t_read_wave_forcing


  !
  ! note that the following TARGET attribute is essential! Otherwise the pointer to the
  ! specific reader inside the time interpolator object (this%reader in time_intp_intp)
  ! will lose its association status.
  TYPE(t_read_wave_forcing), ALLOCATABLE, TARGET :: reader_wave_forcing(:)

  PUBLIC  :: t_read_wave_forcing
  PUBLIC  :: reader_wave_forcing

  PUBLIC  :: construct_reader_wave_forcing
  PUBLIC  :: destruct_reader_wave_forcing

CONTAINS

  !>
  !! Initialize reader of wave forcing data
  !!
  !! Initialize the reader object which reads wave forcing data from file
  !! and performs a linear interpolation in time to a specified target date.
  !!
  SUBROUTINE read_wave_forcing__init (self, p_patch, destination_time, wave_forc_wind_file, &
    &        wave_forc_ice_file, wave_forc_slh_file, wave_forc_osc_file)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_wave_forcing__init'

    CLASS(t_read_wave_forcing),  TARGET     :: self
    TYPE(t_patch),  TARGET,      INTENT(IN) :: p_patch
    TYPE(datetime), POINTER,     INTENT(IN) :: destination_time
    CHARACTER(LEN=filename_max), INTENT(IN) :: wave_forc_wind_file
    CHARACTER(LEN=filename_max), INTENT(IN) :: wave_forc_ice_file
    CHARACTER(LEN=filename_max), INTENT(IN) :: wave_forc_slh_file
    CHARACTER(LEN=filename_max), INTENT(IN) :: wave_forc_osc_file

  !-----------------------------------------

    ! set pointer to current patch
    self%p_patch => p_patch


    IF(my_process_is_mpi_workroot()) THEN
      INQUIRE (FILE=wave_forc_wind_file, EXIST=self%l_wind_exist)
      IF (self%l_wind_exist) THEN
        CALL message(routine,'10m wind from: '//wave_forc_wind_file)
      ELSE
        WRITE(message_text,'(a,a,a)') 'Instant wind data file ', TRIM(wave_forc_wind_file), &
          &                           ' is not found.'
        CALL finish(routine, message_text)
      ENDIF

      INQUIRE (FILE=wave_forc_ice_file, EXIST=self%l_ice_exist)
      IF (self%l_ice_exist) THEN
        CALL message(routine,'ice concentration from: '//wave_forc_ice_file)
      ELSE
        WRITE(message_text,'(a,a,a)') 'Instant sea-ice data file ', TRIM(wave_forc_ice_file), &
          &                           ' is not found. Run without sea ice'
        CALL warning(routine, message_text)
      ENDIF

      INQUIRE (FILE=wave_forc_slh_file, EXIST=self%l_slh_exist)
      IF (self%l_slh_exist) THEN
        CALL message(routine,'sea level height from: '//wave_forc_slh_file)
      ELSE
        WRITE(message_text,'(a,a,a)') 'Instant sea level data file ', TRIM(wave_forc_slh_file), &
          &                           ' is not found. Run without sea level height'
        CALL warning(routine, message_text)
      ENDIF

      INQUIRE (FILE=wave_forc_osc_file, EXIST=self%l_osc_exist)
      IF (self%l_osc_exist) THEN
        CALL message(routine,'ocean surface currents from: '//wave_forc_osc_file)
      ELSE
        WRITE(message_text,'(a,a,a)') 'Instant ocean surface current data file ', TRIM(wave_forc_osc_file), &
          &                           ' is not found. Run without ocean currents'
        CALL warning(routine, message_text)
      ENDIF
    ENDIF
    !
    ! broadcast l_wind_exist, l_ice_exist, l_slh_exist, l_osc_exist from I-PE to WORK PEs
    CALL p_bcast(self%l_wind_exist, p_io, p_comm_work)
    CALL p_bcast(self%l_ice_exist,  p_io, p_comm_work)
    CALL p_bcast(self%l_slh_exist,  p_io, p_comm_work)
    CALL p_bcast(self%l_osc_exist,  p_io, p_comm_work)



    ! Initialize reader and time interpolator
    IF (self%l_wind_exist) THEN
      CALL self%u10_reader%init(p_patch, wave_forc_wind_file)
      CALL self%v10_reader%init(p_patch, wave_forc_wind_file)
      CALL self%u10_intp  %init(self%u10_reader, destination_time, "u_10m")
      CALL self%v10_intp  %init(self%v10_reader, destination_time, "v_10m")
    END IF

    IF (self%l_ice_exist) THEN
      CALL self%sic_reader%init(p_patch, wave_forc_ice_file)
      CALL self%sic_intp  %init(self%sic_reader, destination_time, "fr_seaice")
    END IF

    IF (self%l_slh_exist) THEN
      CALL self%slh_reader%init(p_patch, wave_forc_slh_file)
      CALL self%slh_intp  %init(self%slh_reader, destination_time, "SLH")
    END IF

    IF (self%l_osc_exist) THEN
      CALL self%uosc_reader%init(p_patch, wave_forc_osc_file)
      CALL self%vosc_reader%init(p_patch, wave_forc_osc_file)
      CALL self%uosc_intp  %init(self%uosc_reader, destination_time, "UOSC")
      CALL self%vosc_intp  %init(self%vosc_reader, destination_time, "VOSC")
    END IF

    ! update initialization flag
    self%isInit = .TRUE.
  END SUBROUTINE read_wave_forcing__init


  !>
  !! Read new timelevel of wave forcing data
  !!
  !! A new timelevel of wave forcing data is read from file and
  !! a linear interpolation to the specified destination time is performed.
  !! The result is stored inside the read object itself.
  !!
  SUBROUTINE read_wave_forcing__get_new_rawdata (self, destination_time)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_wave_forcing__get_new_rawdata'

    CLASS(t_read_wave_forcing),   TARGET    :: self
    TYPE(datetime), POINTER,     INTENT(IN) :: destination_time  ! validity time of the new data

    IF (self%l_wind_exist) THEN
      ! get new u10m
      CALL self%u10_intp%intp(destination_time, self%u10m_raw, lacc=.FALSE.)
      ! get new v10m
      CALL self%v10_intp%intp(destination_time, self%v10m_raw, lacc=.FALSE.)
    END IF

    ! get new sic
    IF (self%l_ice_exist) THEN
      CALL self%sic_intp%intp(destination_time, self%sic_raw, lacc=.FALSE.)
    ENDIF

    ! get new slh
    IF (self%l_slh_exist) THEN
      CALL self%slh_intp%intp(destination_time, self%slh_raw, lacc=.FALSE.)
    ENDIF

    IF (self%l_osc_exist) THEN
      ! get new uosc
      CALL self%uosc_intp%intp(destination_time, self%uosc_raw, lacc=.FALSE.)
      ! get new vosc
      CALL self%vosc_intp%intp(destination_time, self%vosc_raw, lacc=.FALSE.)
    END IF

  END SUBROUTINE read_wave_forcing__get_new_rawdata


  !>
  !! The forcing state vector is updated
  !!
  !! The forcing state vector is updated by copying the raw data fields to the
  !! corresponding state vars. Additional diagnostic fields are updated on the basis
  !! of the recent raw data fields, which results in a fully updated forcing state vector.
  !!
  SUBROUTINE read_wave_forcing__update_forcing (self, destination_time, u10m, v10m, &
    &                                           sp10m, dir10m, sic, slh, uosc, vosc, &
    &                                           sp_osc, dir_osc, ice_free_mask_c)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_wave_forcing__update_forcing'

    CLASS(t_read_wave_forcing), TARGET     :: self
    TYPE(datetime), POINTER, INTENT(IN)    :: destination_time ! validity time of the new data
    REAL(wp),                INTENT(INOUT) :: u10m(:,:), v10m(:,:)     ! zonal and meridional wind
                                                                       ! components at 10m asl
    REAL(wp),                INTENT(INOUT) :: sp10m(:,:)               ! wind speed at 10m asl
    REAL(wp),                INTENT(INOUT) :: dir10m(:,:)              ! wind direction at 10m asl [rad]
    REAL(wp),                INTENT(INOUT) :: sic(:,:)                 ! sea ice fraction
    REAL(wp),                INTENT(INOUT) :: slh(:,:)                 ! sea level height
    REAL(wp),                INTENT(INOUT) :: uosc(:,:), vosc(:,:)     ! ocean surface currents
    REAL(wp),                INTENT(INOUT) :: sp_osc(:,:)              ! ocean surface current velocity
    REAL(wp),                INTENT(INOUT) :: dir_osc(:,:)             ! ocean surface current direction [rad]
    INTEGER,                 INTENT(INOUT) :: ice_free_mask_c(:,:)     ! ice mask

    ! local
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: destination_time_string

    ! Sanity check
    IF (.NOT. self%isInit) THEN
      CALL finish(routine, "Error: Forcing state reader has not been initialized!")
    ENDIF

    IF (msg_level > 12) THEN
      CALL datetimeToString(destination_time, destination_time_string)
      WRITE(message_text,'(a,a)') 'Update forcing data for ', TRIM(destination_time_string)
      CALL message(routine, message_text)
    ENDIF

    ! get new forcing data (read from file)
    CALL self%get_new_rawdata(destination_time)

    ! check if the size of the raw data field matches with the
    ! corresponding field in the forcing state
    IF (self%l_wind_exist) THEN
      CALL check_matching_size(arr1=self%u10m_raw, arr2=u10m)
      CALL check_matching_size(arr1=self%v10m_raw, arr2=v10m)
    END IF

    IF (self%l_ice_exist) THEN
      CALL check_matching_size(arr1=self%sic_raw,  arr2=sic)
    ENDIF

    IF (self%l_slh_exist) THEN
      CALL check_matching_size(arr1=self%slh_raw,  arr2=slh)
    ENDIF

    IF (self%l_osc_exist) THEN
      CALL check_matching_size(arr1=self%uosc_raw, arr2=uosc)
      CALL check_matching_size(arr1=self%vosc_raw, arr2=vosc)
    END IF

    ! copy ray data fields to the forcing state vector
    !
!$OMP PARALLEL
    IF (self%l_wind_exist) THEN
      CALL copy(src=self%u10m_raw(:,1,:,1), dest=u10m, lacc=.FALSE.)
      CALL copy(src=self%v10m_raw(:,1,:,1), dest=v10m, lacc=.FALSE.)
    END IF

    IF (self%l_ice_exist) THEN
      CALL copy(src=self%sic_raw  (:,1,:,1), dest=sic, lacc=.FALSE.)
    ENDIF

    IF (self%l_slh_exist) THEN
      CALL copy(src=self%slh_raw  (:,1,:,1), dest=slh, lacc=.FALSE.)
    ENDIF

    IF (self%l_osc_exist) THEN
      CALL copy(src=self%uosc_raw (:,1,:,1), dest=uosc, lacc=.FALSE.)
      CALL copy(src=self%vosc_raw (:,1,:,1), dest=vosc, lacc=.FALSE.)
    END IF
!$OMP END PARALLEL

    ! update additional diagnostic fields
    !

    ! update wind speed and direction
    CALL update_speed_and_direction(p_patch = self%p_patch, &  ! IN
      &                               u     = u10m,         &  ! IN
      &                               v     = v10m,         &  ! IN
      &                              sp     = sp10m,        &  ! OUT
      &                              dir    = dir10m)          ! OUT

    ! update ocean current velocity and direction
    CALL update_speed_and_direction(p_patch = self%p_patch, &  ! IN
      &                               u     = uosc,         &  ! IN
      &                               v     = vosc,         &  ! IN
      &                              sp     = sp_osc,       &  ! OUT
      &                              dir    = dir_osc)         ! OUT

    ! update ice-free mask
    CALL update_ice_free_mask(p_patch       = self%p_patch,   & ! IN
      &                       sea_ice_c     = sic,            & ! IN
      &                       ice_free_mask = ice_free_mask_c)  ! OUT

  CONTAINS

    SUBROUTINE check_matching_size (arr1, arr2)
      REAL(wp), INTENT(IN) :: arr1(:,:,:,:)
      REAL(wp), INTENT(IN) :: arr2(:,:)

      CHARACTER(len=*), PARAMETER :: routine = modname//':check_matching_size'

      IF (SIZE(arr1,1) /= SIZE(arr2,1) .OR. SIZE(arr1,3) /= SIZE(arr2,2)) THEN
        WRITE(message_text,'(a)') 'SIZE mismatch between source and destination array'
        CALL finish(routine, message_text)
      ENDIF
    END SUBROUTINE
  END SUBROUTINE read_wave_forcing__update_forcing


  !>
  !! Destruct reader of wave forcing data
  !!
  !! Destruct the reader object which reads in wave forcing data from file
  !!
  SUBROUTINE read_wave_forcing__deinit (self)

  !  CHARACTER(len=*), PARAMETER :: routine = modname//':read_wave_forcing__deinit'

    CLASS(t_read_wave_forcing)              :: self

    self%p_patch => NULL()

    ! Destruct reader
    IF (self%l_wind_exist) THEN
      CALL self%u10_reader%deinit()
      CALL self%v10_reader%deinit()
    END IF
    IF (self%l_ice_exist)  CALL self%sic_reader%deinit()
    IF (self%l_slh_exist)  CALL self%slh_reader%deinit()
    IF (self%l_osc_exist) THEN
      CALL self%uosc_reader%deinit()
      CALL self%vosc_reader%deinit()
    END IF

    CALL DO_DEALLOCATE(self%u10m_raw)
    CALL DO_DEALLOCATE(self%v10m_raw)
    CALL DO_DEALLOCATE(self%sic_raw)
    CALL DO_DEALLOCATE(self%slh_raw)
    CALL DO_DEALLOCATE(self%uosc_raw)

    self%isInit=.FALSE.
  END SUBROUTINE read_wave_forcing__deinit


  !>
  !! Wrapper for forcing reader construction
  !!
  SUBROUTINE construct_reader_wave_forcing (p_patch, wave_config, tc_start_date)
    TYPE(t_patch),           INTENT(IN) :: p_patch(:)
    TYPE(t_wave_config),     INTENT(IN) :: wave_config(:)
    TYPE(datetime), POINTER, INTENT(IN) :: tc_start_date

    ! local
    INTEGER :: jg, jlev
    INTEGER :: ierrstat
    CHARACTER(LEN=filename_max) :: wave_forc_wind_fn(n_dom) ! forc_file_prefix+'_wind' for U and V 10 meter wind (m/s)
    CHARACTER(LEN=filename_max) :: wave_forc_ice_fn(n_dom)  ! forc_file_prefix+'_ice'  for sea ice concentration (fraction of 1)
    CHARACTER(LEN=filename_max) :: wave_forc_slh_fn(n_dom)  ! forc_file_prefix+'_slh'  for sea level height (m)
    CHARACTER(LEN=filename_max) :: wave_forc_osc_fn(n_dom)  ! forc_file_prefix+'_osc'  for U and V ocean surface currents (m/s)

    CHARACTER(len=*), PARAMETER :: routine = modname//':construct_reader_wave_forcing'

    IF (msg_level > 6) THEN
      CALL message(routine,'Construct wave forcing reader for standalone run')
    ENDIF

    ALLOCATE(reader_wave_forcing(n_dom), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, 'Allocation failed for reader_wave_forcing')

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level

      wave_forc_wind_fn(jg) = generate_filename(TRIM(wave_config(jg)%forc_file_prefix)//"_wind.nc",&
        &                 getModelBaseDir(), nroot, jlev, jg)
      wave_forc_ice_fn(jg)  = generate_filename(TRIM(wave_config(jg)%forc_file_prefix)//"_ice.nc", &
        &                 getModelBaseDir(), nroot, jlev, jg)
      wave_forc_slh_fn(jg)  = generate_filename(TRIM(wave_config(jg)%forc_file_prefix)//"_slh.nc", &
        &                 getModelBaseDir(), nroot, jlev, jg)
      wave_forc_osc_fn(jg)  = generate_filename(TRIM(wave_config(jg)%forc_file_prefix)//"_osc.nc", &
        &                 getModelBaseDir(), nroot, jlev, jg)

      ! initialize reader of external forcing data
      CALL reader_wave_forcing(jg)%init(p_patch             = p_patch(jg),             & !in
        &                               destination_time    = tc_start_date,           & !in
        &                               wave_forc_wind_file = wave_forc_wind_fn(jg),   & !in
        &                               wave_forc_ice_file  = wave_forc_ice_fn(jg),    & !in
        &                               wave_forc_slh_file  = wave_forc_slh_fn(jg),    & !in
        &                               wave_forc_osc_file  = wave_forc_osc_fn(jg) )     !in
    ENDDO

  END SUBROUTINE


  !>
  !! Wrapper for forcing reader destruction
  !!
  SUBROUTINE destruct_reader_wave_forcing ()

    INTEGER :: jg
    INTEGER :: ierrstat
    CHARACTER(len=*), PARAMETER :: routine = modname//':destruct_reader_wave_forcing'

    IF (msg_level > 6) THEN
      CALL message(routine,'Destruct wave forcing reader')
    ENDIF

    IF (ALLOCATED(reader_wave_forcing)) THEN
      DO jg=1,n_dom
        CALL reader_wave_forcing(jg)%deinit()
      ENDDO
      DEALLOCATE(reader_wave_forcing, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, 'Deallocation failed for reader_wave_forcing')
    ENDIF

  END SUBROUTINE

END MODULE mo_wave_forcing
