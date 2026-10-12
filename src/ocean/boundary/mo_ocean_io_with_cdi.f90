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


! This module provides an input instruction list for the ocean
! that is used to determine, which variables may be read from
! which input file. It works analogous to initicon.
!

MODULE mo_oce_io_with_cdi

  USE mo_exception,          ONLY: message, finish
  USE mo_run_config,         ONLY: check_uuid_gracefully
  USE mo_ocean_types,        ONLY: t_hydro_ocean_state
  USE mo_ocean_physics_types,ONLY: t_ho_params
  USE mo_impl_constants,     ONLY: SUCCESS, MODE_DWDANA_OCE, MODE_IAU_OCE, max_dom
  USE mo_dictionary,         ONLY: t_dictionary
  USE mo_model_domain,       ONLY: t_patch,t_patch_3d
  USE mo_ocean_nml,          ONLY: lread_ana_oce, init_mode_oce, lconsistency_checks_oce, n_zlev
  USE mo_input_instructions, ONLY: readInstructionListOce_make, t_readInstructionListPtr
  USE mo_initicon_config,    ONLY: fgFilename, anaFilename, lread_ana, ana_varnames_map_file, &
  &                                is_iau_active, iau_wgt_adv
  USE mo_input_request_list, ONLY: t_InputRequestList, InputRequestList_create
  USE mo_initicon_utils,     ONLY: initVarnamesDict
  USE mo_util_string,        ONLY: int2string
  USE mo_sea_ice_types,      ONLY: t_sea_ice
  USE mo_netcdf
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_grid_config,        ONLY: n_dom
  USE mo_parallel_config,    ONLY: nproma
  USE mo_io_units,           ONLY: filename_max
  USE mo_util_uuid_types,    ONLY: t_uuid
  USE mo_fortran_tools,      ONLY: init
  USE mo_cdi,                ONLY: cdiDefAdditionalKey
  USE mo_sea_ice_nml,        ONLY: kice
  USE mo_ocean_initicono,    ONLY: t_initicono_state, t_pi_oce, t_pi_seaice, t_pi_oce_in, t_pi_seaice_in, &
  &                                fetch_dwdfg_oce, fetch_dwdfg_seaice, fetch_dwdana_oce,                 &
  &                                fetch_dwdana_seaice, t_initicono_read
  USE mo_dynamics_config,    ONLY: nold
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range
  USE mo_ocean_nudging_types,ONLY: t_ocean_nudge
  USE mo_ocean_nudging,      ONLY: ocean_nudge
  USE mtime,                 ONLY: datetime

  !treat input_instructions first

  IMPLICIT NONE
  PRIVATE

  !variables

  PUBLIC :: ana_varnames_dict_oce

  !functions

  PUBLIC :: init_oce, apply_ocean_iau

  CHARACTER(LEN=14), PARAMETER :: modname = 'oceIOWithCDI'

  ! The possible RETURN values of readInstructionList_sourceOfVar().
  ! kInputSourceBoth  : First guess and analysis increment read from file
  ! kInputSourceAnaI  : Analysis interpolated from parent grid
  ! kInputSourceFgAnaI: First guess read from file, analysis increment interpolated from parent grid
  ENUM, BIND(C)
    ENUMERATOR :: kInputSourceUnset = 1, kInputSourceNone, kInputSourceFg, kInputSourceAna, &
    &           kInputSourceBoth, kInputSourceCold, kInputSourceAnaI, kInputSourceFgAnaI
  END ENUM

  ! The possible values for statusFg AND statusAna:
  ENUM, BIND(C)
    ENUMERATOR :: kStateNoFetch = 1, kStateFailedFetch, kStateRead, kStateFailedOptFetch
  END ENUM

  TYPE(t_initicono_state),   ALLOCATABLE, TARGET :: initicono(:)

  ! dictionary which maps internal variable names onto
  ! GRIB2 shortnames or NetCDF var names.
  TYPE (t_dictionary) :: ana_varnames_dict_oce


  TYPE :: t_fetchParams
    TYPE(t_readInstructionListPtr), ALLOCATABLE :: inputInstructions(:)
    CLASS(t_InputRequestList), POINTER :: requestList
    CHARACTER(LEN = :), ALLOCATABLE :: routine
    LOGICAL :: isFg
  END TYPE t_fetchParams

CONTAINS

  FUNCTION gridUuids(p_patch) RESULT(resultVar)
    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    TYPE(t_uuid), DIMENSION(n_dom) :: resultVar

    resultVar(:) = p_patch(:)%grid_uuid
  END FUNCTION gridUuids

  !-------------
  !>
  !! SUBROUTINE init_oce
  !! ICON-O initialization routine: Reads in ICON-O analysis
  !!
  SUBROUTINE init_oce (patch_3d, p_sea_ice, ocean_state, params_oce, read_initicono)

    TYPE(t_patch_3d), INTENT(INOUT), TARGET          :: patch_3d
    TYPE(t_hydro_ocean_state), INTENT(INOUT), TARGET :: ocean_state(:)
    TYPE(t_sea_ice), TARGET, INTENT(INOUT)           :: p_sea_ice
    TYPE(t_ho_params), INTENT(INOUT)                 :: params_oce
    TYPE(t_initicono_read), INTENT(IN)               :: read_initicono


    CHARACTER(LEN = *), PARAMETER :: routine = modname//':init_oce'
    INTEGER :: ist, jg
    TYPE(t_readInstructionListPtr) :: inputInstructions(n_dom)

    ! Allocate initicono data type
    ALLOCATE (initicono(n_dom), stat=ist)
    IF (ist /= SUCCESS)  CALL finish(routine,'allocation for initicon-oce failed')

    DO jg=1,n_dom
      CALL construct_initicono(initicono(jg), patch_3d%p_patch_2D(jg))
    END DO

    ! Read IN the dictionary for the variable names (IF we need it)
    CALL initVarnamesDict(ana_varnames_dict_oce, .FALSE.)

    ! -----------------------------------------------
    ! make the CDI aware of some custom GRIB keys
    ! -----------------------------------------------

    CALL cdiDefAdditionalKey("localInformationNumber")
    CALL cdiDefAdditionalKey("localNumberOfExperiment")
    CALL cdiDefAdditionalKey("typeOfFirstFixedSurface")
    CALL cdiDefAdditionalKey("typeOfGeneratingProcess")
    CALL cdiDefAdditionalKey("backgroundProcess")
    CALL cdiDefAdditionalKey("totalNumberOfTileAttributePairs")
    CALL cdiDefAdditionalKey("tileIndex")
    CALL cdiDefAdditionalKey("tileAttribute")

    ! -----------------------------------------------
    ! generate analysis/FG input instructions
    ! -----------------------------------------------
    DO jg=1,n_dom
      inputInstructions(jg)%ptr => readInstructionListOce_make(patch_3d%p_patch_2D(jg), init_mode_oce, ana_varnames_dict_oce)
    END DO

    ! -----------------------------------------------
    ! READ AND process the input DATA
    ! -----------------------------------------------
    SELECT CASE(init_mode_oce)
      CASE(MODE_DWDANA_OCE)
        CALL message(modname,'MODE_DWD: perform initialization with DWD analysis for oce')
      CASE (MODE_IAU_OCE)
        CALL message(modname,'MODE_IAU_OCE: perform initialization with incremental analysis update for oce')
      CASE DEFAULT
        CALL finish(modname, "Invalid operation mode!")
    END SELECT

    ! read and initialize ICON prognostic fields
    !
    CALL read_dwdfg_oce(patch_3d, patch_3d%p_patch_2D(:), inputInstructions, ocean_state, params_oce, p_sea_ice, read_initicono)
    IF(lread_ana_oce) CALL read_dwdana_oce(patch_3d, patch_3d%p_patch_2D(:), inputInstructions, ocean_state, p_sea_ice, read_initicono)

    CALL deallocate_initicono(initicono)

    DEALLOCATE (initicono, stat=ist)
    IF (ist /= success) CALL finish(routine,'deallocation for initicon-o failed')

    DO jg=1,n_dom
      IF(my_process_is_stdio()) CALL inputInstructions(jg)%ptr%printSummary(jg)
      CALL inputInstructions(jg)%ptr%destruct()
      DEALLOCATE(inputInstructions(jg)%ptr, stat=ist)
      IF(ist /= success) CALL finish(routine,'deallocation of an input instruction list failed')
    ENDDO

  END SUBROUTINE init_oce

  !-------------
  !>
  !! SUBROUTINE deallocate_initicono
  !! Deallocates the components of the initicono data type
  !!
  SUBROUTINE deallocate_initicono (initicono)

    TYPE(t_initicono_state), INTENT(INOUT) :: initicono(:)

    INTEGER :: jg
  !------------------------------------------------------------------

    ! call destructor
    DO jg=1,n_dom
      CALL initicono(jg)%finalize()
    ENDDO

    ! destroy variable name dictionaries:
    CALL ana_varnames_dict_oce%finalize()

  END SUBROUTINE deallocate_initicono

  !-------------
  !>  SUBROUTINE construct_initicono
  !! Ensures that all fields have a defined VALUE.
  !!   * resets all linitialized flags
  !!   * copies topography AND coordinate surfaces
  !!   * allocates the fields we USE
  !!       * zeros OUT these fields to ensure deteministic checksums
  !!   * nullificates all other pointers
  !!
  !! This initalizes all ALLOCATED memory to avoid nondeterministic
  !! checksums when ONLY a part of a field IS READ from file due to
  !! nonfull blocks.


  SUBROUTINE construct_initicono(initicono, p_patch)

    TYPE(t_initicono_state), INTENT(INOUT) :: initicono
    TYPE(t_patch)                          :: p_patch

    ! Local variables: loop control and dimensions
    INTEGER :: nlev, nblks_c, nblks_e

    nlev = p_patch%nlev
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    !WS 2017-04-12:  ORDERED was added here to work around a CCE 8.5.5 bug
    !$OMP ORDERED
    CALL construct_oce_in(initicono%oce_in)
    CALL construct_seaice_in(initicono%seaice_in)
    CALL construct_oce(initicono%oce)
    CALL construct_oce_inc(initicono%oce_inc)
    CALL construct_seaice(initicono%seaice)
    CALL construct_seaice_inc(initicono%seaice_inc)
    !$OMP END ORDERED

    !----------------------------------------------------------------------

    CONTAINS

    SUBROUTINE construct_oce_in(oce_in)
      TYPE(t_pi_oce_in),        INTENT(INOUT) :: oce_in

      ALLOCATE(oce_in%to (nproma,nlev,nblks_c),  &
      &       oce_in%so  (nproma,nlev,nblks_c),  &
      &       oce_in%u   (nproma,nlev,nblks_c),  &
      &       oce_in%v   (nproma,nlev,nblks_c),  &
      &       oce_in%vn  (nproma,nlev,nblks_e),  &
      &       oce_in%zos      (nproma,nblks_c),  &
      &       oce_in%stretch_c(nproma,nblks_c))
      oce_in%nlev         = 72
      oce_in%linitialized = .TRUE.

    END SUBROUTINE construct_oce_in

    SUBROUTINE construct_seaice_in(seaice_in)
      TYPE(t_pi_seaice_in), INTENT(INOUT) :: seaice_in

      seaice_in%lev          = 0
      seaice_in%linitialized = .FALSE.
    END SUBROUTINE construct_seaice_in

    ! Allocate ocean output data
    SUBROUTINE construct_oce(oce)
      TYPE(t_pi_oce), INTENT(INOUT) :: oce

      ALLOCATE(oce%to  (nproma,nlev,nblks_c), &
      &        oce%v   (nproma,nlev,nblks_c), &
      &        oce%u   (nproma,nlev,nblks_c), &
      &        oce%vn  (nproma,nlev,nblks_e), &
      &        oce%so  (nproma,nlev,nblks_c), &
      &        oce%zos       (nproma,nblks_c), &
      &        oce%stretch_c (nproma,nblks_c))

     !$OMP PARALLEL
      CALL init(oce%to(:,:,:), lacc=.FALSE.)
      CALL init(oce%u(:,:,:), lacc=.FALSE.)
      CALL init(oce%v(:,:,:), lacc=.FALSE.)
      CALL init(oce%vn(:,:,:), lacc=.FALSE.)
      CALL init(oce%so(:,:,:), lacc=.FALSE.)
      CALL init(oce%zos(:,:), lacc=.FALSE.)
      CALL init(oce%stretch_c(:,:), lacc=.FALSE.)
      !$OMP END PARALLEL

      oce%nlev         = nlev
      oce%linitialized = .TRUE.
    END SUBROUTINE construct_oce

    ! ocean assimilation increments
    SUBROUTINE construct_oce_inc(oce_inc)
      TYPE(t_pi_oce), INTENT(INOUT) :: oce_inc

      IF ( init_mode_oce == MODE_IAU_OCE ) THEN
        ALLOCATE(oce_inc%to (nproma,nlev,nblks_c), &
        &        oce_inc%so (nproma,nlev,nblks_c), &
        &        oce_inc%u   (nproma,nlev,nblks_c), &
        &        oce_inc%v   (nproma,nlev,nblks_c), &
        &        oce_inc%vn  (nproma,nlev,nblks_e), &
        &        oce_inc%zos       (nproma,nblks_c), &
        &        oce_inc%stretch_c (nproma,nblks_c) )
        !$OMP PARALLEL
        CALL init(oce_inc%to(:,:,:), lacc=.FALSE.)
        CALL init(oce_inc%so(:,:,:), lacc=.FALSE.)
        CALL init(oce_inc%u(:,:,:), lacc=.FALSE.)
        CALL init(oce_inc%v(:,:,:), lacc=.FALSE.)
        CALL init(oce_inc%vn(:,:,:), lacc=.FALSE.)
        CALL init(oce_inc%zos(:,:), lacc=.FALSE.)
        CALL init(oce_inc%stretch_c(:,:), lacc=.FALSE.)
        !$OMP END PARALLEL

        oce_inc%nlev         = nlev
        oce_inc%linitialized = .TRUE.
      ELSE
        oce_inc%nlev         = 0
        oce_inc%linitialized = .FALSE.
      ENDIF
    END SUBROUTINE construct_oce_inc

    ! Allocate surface output data
    SUBROUTINE construct_seaice(seaice)
      TYPE(t_pi_seaice), INTENT(INOUT) :: seaice

      ALLOCATE(seaice%hi   (nproma, kice, nblks_c ), &
      &        seaice%hs   (nproma, kice, nblks_c ), &
      &        seaice%conc (nproma, kice, nblks_c ))

      !$OMP PARALLEL
      CALL init(seaice%hi(:,:,:), lacc=.FALSE.)
      CALL init(seaice%hs(:,:,:), lacc=.FALSE.)
      CALL init(seaice%conc(:,:,:), lacc=.FALSE.)
      !$OMP END PARALLEL

      seaice%linitialized = .TRUE.

    END SUBROUTINE construct_seaice

    ! surface assimilation increments
    SUBROUTINE construct_seaice_inc(seaice_inc)
      TYPE(t_pi_seaice), INTENT(INOUT) :: seaice_inc

      IF ( init_mode_oce == MODE_IAU_OCE ) THEN
        ALLOCATE(seaice_inc%hi (nproma,kice, nblks_c ), &
        &        seaice_inc%hs (nproma,kice, nblks_c ), &
        &        seaice_inc%conc (nproma,kice, nblks_c ))
        !$OMP PARALLEL
        CALL init(seaice_inc%hi(:,:,:), lacc=.FALSE.)
        CALL init(seaice_inc%hs(:,:,:), lacc=.FALSE.)
        CALL init(seaice_inc%conc(:,:,:), lacc=.FALSE.)
        !$OMP END PARALLEL

        seaice_inc%linitialized = .TRUE.
      ELSE
        seaice_inc%linitialized = .FALSE.
      ENDIF
    END SUBROUTINE construct_seaice_inc

  END SUBROUTINE construct_initicono

  ! Read the data from the first-guess file.
  SUBROUTINE read_dwdfg_oce(patch_3d, p_patch, inputInstructions, ocean_state, params_oce, p_sea_ice, read_initicono)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_patch), INTENT(INOUT) :: p_patch(:)
    TYPE(t_readInstructionListPtr) :: inputInstructions(n_dom)
    TYPE(t_hydro_ocean_state), INTENT(INOUT) :: ocean_state(:)
    TYPE(t_ho_params), INTENT(INOUT) :: params_oce
    TYPE(t_sea_ice), INTENT(INOUT) :: p_sea_ice
    TYPE(t_initicono_read) :: read_initicono

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":read_dwdfg_oce"
    INTEGER :: jg, jg1
    CLASS(t_InputRequestList), POINTER :: requestList
    CHARACTER(LEN=filename_max) :: fgFilename_str(max_dom)

    !The input file paths & types are NOT initialized IN all modes, so we need to avoid creating InputRequestLists IN these cases.
    SELECT CASE(init_mode_oce)
      CASE(MODE_DWDANA_OCE, MODE_IAU_OCE)
      CASE DEFAULT
        CALL finish(routine, "assertion failed: unknown ocean init_mode")
    END SELECT

    DO jg = 1, n_dom

      ! Create a request list for all the relevant variable names.
      requestList => InputRequestList_create()
      CALL inputInstructions(jg)%ptr%fileRequests(requestList, lIsFg = .TRUE.)

      fgFilename_str(jg) = " "
      fgFilename_str(jg) = fgFilename(p_patch(jg))

      IF (my_process_is_stdio()) THEN
        ! consistency check: check for duplicate file names which may
        ! occur, for example, if the keyword pattern (namelist
        ! parameter) has been defined ambiguously by the user.
        DO jg1 = 1,(jg-1)
          IF (fgFilename_str(jg1) == fgFilename_str(jg)) THEN
            CALL finish(routine, "Error! Namelist parameter fgFilename has been defined ambiguously "//&
            & "for domains "//TRIM(int2string(jg1, '(i0)'))//" and "//TRIM(int2string(jg, '(i0)'))//"!")
          END IF
        END DO
      END IF
    END DO

    ! Scan the input files AND distribute the relevant variables across the processes.
    DO jg = 1, n_dom
      IF(my_process_is_stdio()) THEN
        CALL message(routine, 'read oce_FG fields from '//TRIM(fgFilename_str(jg)))
      ENDIF  ! p_io
      IF (ana_varnames_map_file /= ' ') THEN

        CALL requestList%readFile(p_patch(jg), TRIM(fgFilename_str(jg)), .TRUE., &
        &                       opt_dict = ana_varnames_dict_oce)
      ELSE
        CALL requestList%readFile(p_patch(jg), TRIM(fgFilename_str(jg)), .TRUE.)
      END IF
      IF(my_process_is_stdio()) THEN
        CALL requestList%printInventory()
        IF(lconsistency_checks_oce) THEN
          CALL requestList%checkRuntypeAndUuids([CHARACTER(LEN=1)::], gridUuids(p_patch), lIsFg=.TRUE., &
          &  lHardCheckUuids=.NOT.check_uuid_gracefully)
        END IF
      END IF
    END DO

    ! Fetch the input DATA from the request list.
    CALL fetch_dwdfg_oce(patch_3d, requestList, ocean_state, params_oce, inputInstructions, read_initicono)
    CALL fetch_dwdfg_seaice(requestList, p_sea_ice, inputInstructions, read_initicono)

    ! Cleanup.
    CALL requestList%destruct()
    DEALLOCATE(requestList)
  END SUBROUTINE read_dwdfg_oce

  ! Read data from analysis files.
  SUBROUTINE read_dwdana_oce(patch_3d, p_patch, inputInstructions, ocean_state, p_sea_ice, read_initicono)
    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_patch), INTENT(INOUT) :: p_patch(:)
    TYPE(t_readInstructionListPtr) :: inputInstructions(n_dom)
    TYPE(t_hydro_ocean_state), INTENT(INOUT), TARGET :: ocean_state(:)
    TYPE(t_sea_ice), TARGET, INTENT(INOUT) :: p_sea_ice
    TYPE(t_initicono_read) :: read_initicono

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":read_dwdana_oce"
#if !defined __GFORTRAN__ || __GNUC__ >= 6

    CHARACTER(LEN = :), ALLOCATABLE :: incrementsList(:)
#else
    CHARACTER(LEN = 3) :: incrementsList_IAU_OCE(5)
    CHARACTER(LEN = 1) :: incrementsList_DEFAULT(1)
#endif
    CLASS(t_InputRequestList), POINTER :: requestList
    CHARACTER(LEN=filename_max) :: anaFilename_str(max_dom)
    INTEGER :: jg, jg1

    !The input file paths & types are not initialized in all modes, so we need to avoid creating InputRequestLists in these cases.
    SELECT CASE(init_mode_oce)
      CASE(MODE_DWDANA_OCE, MODE_IAU_OCE)
      CASE DEFAULT
        CALL finish(routine, "assertion failed: unknown ocean init_mode")
    END SELECT

    ! Create a request list for all the relevant variable names.
    requestList => InputRequestList_create()
    DO jg = 1, n_dom
      CALL inputInstructions(jg)%ptr%fileRequests(requestList, lIsFg = .FALSE.)
    END DO

    ! Scan the input files AND distribute the relevant variables across the processes.
    DO jg = 1, n_dom
      anaFilename_str(jg) = ""
      anaFilename_str(jg) = anaFilename(p_patch(jg))

      IF (my_process_is_stdio()) THEN
        ! consistency check: check for duplicate file names which may
        ! occur, for example, if the keyword pattern (namelist
        ! parameter) has been defined ambiguously by the user.
        DO jg1 = 1,(jg-1)
          IF (anaFilename_str(jg1) == anaFilename_str(jg)) THEN
            CALL finish(routine, "Error! Namelist parameter anaFilename has been defined ambiguously "//&
            & "for domains "//TRIM(int2string(jg1, '(i0)'))//" and "//TRIM(int2string(jg, '(i0)'))//"!")
          END IF
        END DO
      END IF
    END DO

    DO jg = 1, n_dom
      IF(lread_ana_oce) THEN
        lread_ana = .TRUE.
        IF(my_process_is_stdio()) THEN
          CALL message(routine, 'read oce_ANA fields from '//TRIM(anaFilename_str(jg)))
        ENDIF  ! p_io
        IF (ana_varnames_map_file /= ' ') THEN
          CALL requestList%readFile(p_patch(jg), TRIM(anaFilename_str(jg)), .FALSE., &
          &                       opt_dict = ana_varnames_dict_oce)
        ELSE
          CALL requestList%readFile(p_patch(jg), TRIM(anaFilename_str(jg)), .FALSE.)
        END IF
      END IF
    END DO
    IF(my_process_is_stdio()) THEN
      CALL requestList%printInventory()
      IF(lconsistency_checks_oce) THEN
        ! Workaround for GNU compiler (<6.0), which still does not fully support deferred length character arrays
        ! Make use of deferred length character arrays if the GNU compiler is not used, or if
        ! its version number is at least equal to 6.0.
#if !defined __GFORTRAN__ || __GNUC__ >= 6
        SELECT CASE(init_mode_oce)
          CASE(MODE_IAU_OCE)
            incrementsList = [CHARACTER(LEN=3) :: 'u', 'v', 'normal_velocity', 'to', 'so', 'zos']
          CASE DEFAULT
            incrementsList = [CHARACTER(LEN=1) :: ]
        END SELECT
        CALL requestList%checkRuntypeAndUuids(incrementsList, gridUuids(p_patch), lIsFg = .FALSE., &
          &    lHardCheckUuids = .NOT.check_uuid_gracefully)
#else
        SELECT CASE(init_mode_oce)
          CASE(MODE_IAU_OCE)
            incrementsList_IAU_OCE = (/'u  ', 'v  ', 'normal_velocity', 'to', 'so ', 'zos' /)
            CALL requestList%checkRuntypeAndUuids(incrementsList_IAU_OCE, gridUuids(p_patch), lIsFg = .FALSE., &
              &    lHardCheckUuids = .NOT.check_uuid_gracefully)
            write(0,*) "incrementsList_IAU_OCE: ", incrementsList_IAU_OCE
          CASE DEFAULT
            incrementsList_DEFAULT = (/' '/)
            CALL requestList%checkRuntypeAndUuids(incrementsList_DEFAULT, gridUuids(p_patch), lIsFg = .FALSE., &
              &    lHardCheckUuids = .NOT.check_uuid_gracefully)
        END SELECT
#endif
      END IF
    END IF
    ! Fetch the input DATA from the request list.
    SELECT CASE(init_mode_oce)
      CASE(MODE_DWDANA_OCE, MODE_IAU_OCE)
        IF(lread_ana_oce) THEN
          CALL fetch_dwdana_oce(patch_3d, requestList, ocean_state, initicono, inputInstructions, read_initicono)
          CALL fetch_dwdana_seaice(requestList, p_sea_ice, initicono, inputInstructions, read_initicono)
        ENDIF
        IF(init_mode_oce .EQ. MODE_IAU_OCE ) THEN
          ocean_nudge%temp_incr = initicono(1)%oce_inc%to
          ocean_nudge%sal_incr = initicono(1)%oce_inc%so
          ocean_nudge%u_incr = initicono(1)%oce_inc%u
          ocean_nudge%v_incr = initicono(1)%oce_inc%v
          ocean_nudge%vn_incr = initicono(1)%oce_inc%vn
          ocean_nudge%zos_incr = initicono(1)%oce_inc%zos
          ocean_nudge%stretch_incr = initicono(1)%oce_inc%stretch_c
          ocean_nudge%hi_incr = initicono(1)%seaice_inc%hi
          ocean_nudge%hs_incr = initicono(1)%seaice_inc%hs
          ocean_nudge%conc_incr = initicono(1)%seaice_inc%conc
        ENDIF
    END SELECT

    ! Cleanup.
    CALL requestList%destruct()
    DEALLOCATE(requestList)
  END SUBROUTINE read_dwdana_oce


  SUBROUTINE apply_ocean_iau(patch_3d, ocean_state, ocean_nudge, p_sea_ice, read_initicono, current_time)
    TYPE(t_patch_3d), INTENT(INOUT) :: patch_3d
    TYPE(t_hydro_ocean_state), INTENT(INOUT), TARGET :: ocean_state(:)
    TYPE(t_sea_ice), TARGET, INTENT(INOUT) :: p_sea_ice
    TYPE(t_initicono_read) :: read_initicono
    TYPE(t_subset_range), POINTER :: all_cells, all_edges
    TYPE(t_ocean_nudge), TARGET,INTENT(inout)   :: ocean_nudge
    INTEGER :: b, idx, level, start_index, end_index, k
    TYPE(datetime), INTENT(IN), POINTER   :: current_time
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":apply_ocean_iau"


    IF(.NOT. is_iau_active) return

    start_index = 1
    end_index = 1

    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    !Treat cell centers
    DO b = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, b, start_index, end_index)
      DO idx = start_index, end_index
        DO level = 1, n_zlev
! Treat 3-d fields
          ocean_state(1)%p_prog(nold(1))%tracer(idx,level,b,1) = &
          & ocean_state(1)%p_prog(nold(1))%tracer(idx,level,b,1) + iau_wgt_adv*ocean_nudge%temp_incr(idx,level,b)
          ocean_state(1)%p_prog(nold(1))%tracer(idx,level,b,2) = &
          & ocean_state(1)%p_prog(nold(1))%tracer(idx,level,b,2) + iau_wgt_adv*ocean_nudge%sal_incr(idx,level,b)
          ocean_state(1)%p_diag%u(idx,level,b) = &
          & ocean_state(1)%p_diag%u(idx,level,b) + iau_wgt_adv*ocean_nudge%u_incr(idx,level,b)
          ocean_state(1)%p_diag%v(idx,level,b) = &
          & ocean_state(1)%p_diag%v(idx,level,b) + iau_wgt_adv*ocean_nudge%v_incr(idx,level,b)
        ENDDO

        DO k= 1, p_sea_ice%kice
          p_sea_ice%hi(idx,k,b) = p_sea_ice%hi(idx,k,b) + iau_wgt_adv*ocean_nudge%hi_incr(idx,k,b)
          p_sea_ice%hs(idx,k,b) = p_sea_ice%hs(idx,k,b) + iau_wgt_adv*ocean_nudge%hs_incr(idx,k,b)
          p_sea_ice%conc(idx,k,b) = p_sea_ice%conc(idx,k,b) + iau_wgt_adv*ocean_nudge%conc_incr(idx,k,b)
        ENDDO
        !Treat 2-d fields
        ocean_state(1)%p_prog(nold(1))%h(idx,b) = &
        & ocean_state(1)%p_prog(nold(1))%h(idx,b) + iau_wgt_adv*ocean_nudge%zos_incr(idx,b)
        ocean_state(1)%p_prog(nold(1))%stretch_c(idx,b) = &
        & ocean_state(1)%p_prog(nold(1))%stretch_c(idx,b) + iau_wgt_adv*ocean_nudge%stretch_incr(idx,b)
      ENDDO
    ENDDO

    all_edges => patch_3d%p_patch_2D(1)%edges%all
    !Treat edges
    DO b = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, b, start_index, end_index)
      DO idx = start_index, end_index
        DO level = 1, n_zlev
          ocean_state(1)%p_prog(nold(1))%vn(idx,level,b) = &
          & ocean_state(1)%p_prog(nold(1))%vn(idx,level,b) + iau_wgt_adv*ocean_nudge%vn_incr(idx,level,b)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE

END MODULE mo_oce_io_with_cdi
