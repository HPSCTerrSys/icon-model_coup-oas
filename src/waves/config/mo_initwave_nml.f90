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

! Namelist for wave model initialization
!
! The following routine is called by read_wave_namelist and controls
! specifics of the wave model initialization.

MODULE mo_initwave_nml

  USE mo_kind,                      ONLY: wp
  USE mo_exception,                 ONLY: finish, message, print_value
  USE mo_io_units,                  ONLY: nnml, nnml_output, filename_max
  USE mo_master_control,            ONLY: use_restart_namelists
  USE mo_namelist,                  ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_restart_nml_and_att,       ONLY: open_tmpfile, store_and_close_namelist,     &
    &                                     open_and_restore_namelist, close_tmpfile
  USE mo_mpi,                       ONLY: my_process_is_stdio
  USE mo_nml_annotate,              ONLY: temp_defaults, temp_settings
  USE mo_impl_constants,            ONLY: max_dom
  USE mo_wave_constants,            ONLY: MODE_COLD, MODE_ANA
  USE mo_time_config,               ONLY: set_tc_timeshift
  USE mo_initwave_config,           ONLY: initwave_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_initwave_namelist

CONTAINS

  !>
  !! Read namelist for wave model initialization
  !!
  !! This subroutine
  !! - reads the Namelist for wave model initialization
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  SUBROUTINE read_initwave_namelist (filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: istat, funit
    INTEGER :: iunit
    INTEGER :: jg

    REAL(wp):: dt_shift      !< time interval by which the actual start date (tc_start_date)
                             !< is shifted backwards in time. [s]
    INTEGER :: init_mode     !< MODE_ANA : read wave energy spectrum from analysis file
                             !< MODE_COLD: initialize by analytic wind-speed based parameterization
                             !             (such as JONSWAP)
    LOGICAL :: lskip_inv_post_op !< skip inverse post op (i.e. scaling of the input variable)
                                 !  when reading the state file
    CHARACTER(LEN=filename_max) :: initial_wave_spectrum_filename

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_initwave_nml: read_initwave_nml'

    NAMELIST /initwave_nml/ dt_shift, init_mode, &
      &      lskip_inv_post_op, initial_wave_spectrum_filename

    !-----------------------
    ! 1. default settings
    !-----------------------

    dt_shift  = 0._wp     ! no shift backwards in time.
                          ! => tc_current_date = tc_start_date at model start
    init_mode = MODE_COLD ! coldstart from JONSWAP

    lskip_inv_post_op = .FALSE. ! do not skip inverse post op when reading the state file

    initial_wave_spectrum_filename = ""   ! <path>/wave_DOM<dom>_<y><m><d>T<h><min><sec>Z.nc or <path>/wave_DOM<dom>_ML_<num>.nc

    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, initwave_nml)   ! write defaults to temporary text file
    END IF

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('initwave_nml')
      READ(funit,NML=initwave_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------

    CALL open_nml(TRIM(filename))
    CALL position_nml ('initwave_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, initwave_nml)       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, initwave_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------


    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------
    DO jg = 1,max_dom
      initwave_config(jg)%init_mode = init_mode
      initwave_config(jg)%initial_wave_spectrum_filename = initial_wave_spectrum_filename
      initwave_config(jg)%lskip_inv_post_op = lskip_inv_post_op
    ENDDO

    ! transfer dt_shift to time_config state
    CALL set_tc_timeshift(dt_shift)

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=initwave_nml)
      CALL store_and_close_namelist(funit, 'initwave_nml')
    ENDIF


    !-----------------------------------------------------
    ! 7. write the contents of the namelist to an ASCII file
    !-----------------------------------------------------
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=initwave_nml)

  END SUBROUTINE read_initwave_namelist

END MODULE mo_initwave_nml
