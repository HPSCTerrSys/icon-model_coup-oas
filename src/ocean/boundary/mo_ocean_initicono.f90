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

! This module contains the I/O routines for initicon

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ocean_initicono

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: finish
  USE mo_impl_constants,          ONLY: SUCCESS, MODE_IAU_OCE, max_ntracer
  USE mo_impl_constants,          ONLY: max_ntracer
  USE mo_ocean_types,             ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,           ONLY: t_sea_ice
  USE mo_ocean_nml,               ONLY: init_mode_oce, vert_cor_type, fillValue, &
  &                                     n_zlev, fg_apply_fillValues
  USE mo_initicon_types,          ONLY: t_pi_tracer
  USE mtime,                      ONLY: datetime
  USE mo_run_config,              ONLY: ntracer
  USE mo_input_instructions,      ONLY: t_readInstructionListPtr, &
  &                               kInputSourceFg, kInputSourceAna, kInputSourceBoth
  USE mo_input_request_list,      ONLY: t_InputRequestList
  USE mo_initicon_types,          ONLY: t_init_state_const
  USE mo_dynamics_config,         ONLY: nold
  USE mo_fortran_tools,           ONLY: DO_DEALLOCATE, DO_PTR_DEALLOCATE
  USE mo_initicon_io,             ONLY: fetch3d, fetch3d_with_status, fetchSurface, t_fetchParams
  USE mo_ocean_physics_types,     ONLY: t_ho_params
  USE mo_grid_subset,             ONLY: t_subset_range, get_index_range
  USE mo_model_domain,            ONLY: t_patch_3d

  IMPLICIT NONE
  PRIVATE

  !
  !variables

  PUBLIC :: t_inito_state
  PUBLIC :: t_initicono_state  !> state vector for initicon-o
  PUBLIC :: t_pi_oce_in
  PUBLIC :: t_pi_seaice_in
  PUBLIC :: t_pi_oce
  PUBLIC :: t_pi_seaice
  PUBLIC :: t_pi_seaice_inc
  PUBLIC :: t_initicono_read
  PUBLIC :: read_initicono

  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_ocean_initicono'


  ! functions
  PUBLIC :: t_inito_state_finalize
  PUBLIC :: fetch_dwdfg_oce, fetch_dwdfg_seaice
  PUBLIC :: fetch_dwdana_oce, fetch_dwdana_seaice


  ! ocean input variables
  TYPE :: t_pi_oce_in

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    ! vertical dimension of 3D input fields
    INTEGER :: nlev

    REAL(wp), POINTER, DIMENSION(:,:,:) :: to    => NULL(), &
    &                                    so      => NULL(), &
    &                                    u       => NULL(), &
    &                                    v       => NULL(), &
    &                                    vn       => NULL()
    REAL(wp), POINTER, DIMENSION(:,:) :: stretch_c => NULL(), &
    &                                    zos     => NULL()

    TYPE (t_pi_tracer), DIMENSION(max_ntracer) :: tracer

    CONTAINS
      PROCEDURE :: finalize => t_pi_oce_in_finalize   !< destructor
  END TYPE t_pi_oce_in

  ! surface input variables
  TYPE :: t_pi_seaice_in

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    ! number of ice levels (should be one)
    INTEGER :: lev

    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: hs, hi, conc

    CONTAINS
      PROCEDURE :: finalize => t_pi_seaice_in_finalize !< destructor
  END TYPE t_pi_seaice_in


  TYPE :: t_pi_oce

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    ! vertical dimension of 3D input fields
    INTEGER :: nlev

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: u, v, vn, to, so
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: stretch_c, zos

    TYPE (t_pi_tracer), DIMENSION(max_ntracer) :: tracer

    CONTAINS
      PROCEDURE :: finalize => t_pi_oce_finalize   !< destructor
  END TYPE t_pi_oce

  !
  TYPE :: t_pi_seaice

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: hi, hs, conc

    CONTAINS
      PROCEDURE :: finalize => t_pi_seaice_finalize   !< destructor
  END TYPE t_pi_seaice


  ! seaice field increments
  TYPE :: t_pi_seaice_inc
    !
    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:)   :: hi
    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:)   :: hs
    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:)   :: conc

    CONTAINS
      PROCEDURE :: finalize => t_pi_seaice_inc_finalize   !< destructor
  END TYPE t_pi_seaice_inc

  ! state vector type: base class
  !
  TYPE :: t_inito_state

    TYPE (t_pi_oce_in)     :: oce_in
    TYPE (t_pi_seaice_in)  :: seaice_in
    TYPE (t_pi_oce)        :: oce
    TYPE (t_pi_seaice)     :: seaice

    TYPE (t_init_state_const), POINTER :: const => NULL()

    TYPE (datetime)        :: vDateTime  ! validity DateTime
    CONTAINS
      PROCEDURE, PUBLIC :: finalize => t_inito_state_finalize
  END TYPE t_inito_state

  ! Adapt the types of mo_initicon_types for the ocean case. Name it initicono.
  ! complete state vector type
  !
  TYPE, EXTENDS(t_inito_state) :: t_initicono_state

    TYPE (t_pi_oce)      :: oce_inc
    TYPE (t_pi_seaice)   :: seaice_inc

    CONTAINS
      PROCEDURE, PUBLIC :: finalize => t_initicono_state_finalize
  END TYPE t_initicono_state


  !This type contains all variables that can be read
  !with initicon-o. For now this needs to be extended
  !manually.
  TYPE :: t_initicono_read

    LOGICAL :: u = .FALSE., v = .FALSE., vn = .FALSE.
    LOGICAL :: to= .FALSE., so = .FALSE.
    LOGICAL :: zos = .FALSE., stretch_c = .FALSE.
    LOGICAL :: hi = .FALSE., hs = .FALSE., conc = .FALSE.

  END TYPE t_initicono_read

  TYPE(t_initicono_read) :: read_initicono

  CHARACTER(LEN=4) :: to_var


  CONTAINS

  SUBROUTINE t_pi_oce_in_finalize(oce_in)
    CLASS(t_pi_oce_in), INTENT(INOUT) :: oce_in

    INTEGER :: idx

    oce_in%linitialized = .FALSE.
    CALL DO_PTR_DEALLOCATE(oce_in%to)
    CALL DO_PTR_DEALLOCATE(oce_in%so)
    CALL DO_PTR_DEALLOCATE(oce_in%u)
    CALL DO_PTR_DEALLOCATE(oce_in%v)
    CALL DO_PTR_DEALLOCATE(oce_in%vn)
    CALL DO_PTR_DEALLOCATE(oce_in%zos)
    CALL DO_PTR_DEALLOCATE(oce_in%stretch_c)

    DO idx=1, ntracer
      CALL oce_in%tracer(idx)%finalize()
    END DO
  END SUBROUTINE t_pi_oce_in_finalize

  SUBROUTINE t_pi_seaice_in_finalize(seaice_in)
    CLASS(t_pi_seaice_in), INTENT(INOUT) :: seaice_in

    seaice_in%linitialized = .FALSE.
    CALL DO_DEALLOCATE(seaice_in%hi)
    CALL DO_DEALLOCATE(seaice_in%hs)
    CALL DO_DEALLOCATE(seaice_in%conc)
  END SUBROUTINE t_pi_seaice_in_finalize


  SUBROUTINE t_pi_oce_finalize(oce)
    CLASS(t_pi_oce), INTENT(INOUT) :: oce

    INTEGER :: idx

    oce%linitialized = .FALSE.
    CALL DO_DEALLOCATE(oce%to)
    CALL DO_DEALLOCATE(oce%so)
    CALL DO_DEALLOCATE(oce%u)
    CALL DO_DEALLOCATE(oce%v)
    CALL DO_DEALLOCATE(oce%vn)
    CALL DO_DEALLOCATE(oce%zos)
    CALL DO_DEALLOCATE(oce%stretch_c)

    DO idx=1, ntracer
      CALL oce%tracer(idx)%finalize()
    END DO
  END SUBROUTINE t_pi_oce_finalize


  SUBROUTINE t_pi_seaice_finalize(seaice)
    CLASS(t_pi_seaice), INTENT(INOUT) :: seaice

    seaice%linitialized = .FALSE.
    CALL DO_DEALLOCATE(seaice%hs)
    CALL DO_DEALLOCATE(seaice%hi)
    CALL DO_DEALLOCATE(seaice%conc)

  END SUBROUTINE t_pi_seaice_finalize


  SUBROUTINE t_pi_seaice_inc_finalize(seaice_inc)
    CLASS(t_pi_seaice_inc), INTENT(INOUT) :: seaice_inc

    seaice_inc%linitialized = .FALSE.
    CALL DO_DEALLOCATE(seaice_inc%hi)
    CALL DO_DEALLOCATE(seaice_inc%hs)
    CALL DO_DEALLOCATE(seaice_inc%conc)
  END SUBROUTINE t_pi_seaice_inc_finalize


  SUBROUTINE t_inito_state_finalize(inito_data)
    CLASS(t_inito_state), INTENT(INOUT) :: inito_data

    CALL inito_data%oce_in%finalize()
    CALL inito_data%seaice_in%finalize()
    CALL inito_data%oce%finalize()
    CALL inito_data%seaice%finalize()

  END SUBROUTINE t_inito_state_finalize


  SUBROUTINE t_initicono_state_finalize(inito_data)
    CLASS(t_initicono_state), INTENT(INOUT) :: inito_data

    ! call base class destructor
    CALL t_inito_state_finalize(inito_data)

    CALL inito_data%oce_inc%finalize()
    CALL inito_data%seaice_inc%finalize()
  END SUBROUTINE t_initicono_state_finalize

  !>
  !! Fetch the DWD first guess from the request list (ocean only)
  !! First guess (FG) is read for to, so, u, v (or vn), zos
  !! whereas DA output is read for to, so, u, v (or vn), zos
  SUBROUTINE fetch_dwdfg_oce(patch_3d, requestList, ocean_state, params_oce, inputInstructions, read_initicono)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_hydro_ocean_state), INTENT(INOUT), TARGET :: ocean_state(:)
    TYPE(t_ho_params), INTENT(INOUT)              :: params_oce
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)
    TYPE(t_initicono_read) :: read_initicono

    CHARACTER(*), PARAMETER :: routine = modname//':fetch_dwdfg_oce'
    TYPE(t_fetchParams)      :: params
    REAL(wp), POINTER :: my_ptr3d(:,:,:)
    LOGICAL :: lfound_u, lfound_v, lfound_vn, lfound_to, lfound_so, lfound_tke

    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)))
    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .TRUE.

      !request the first guess fields (ocean only)
      IF(vert_cor_type .eq. 0) THEN
        CALL fetchSurface(params, 'zos', 1, ocean_state(1)%p_prog(nold(1))%h)
        IF(fg_apply_fillValues) CALL add_missval_2d(patch_3d, ocean_state(1)%p_prog(nold(1))%h)
      ELSEIF(vert_cor_type .eq. 1) THEN
        CALL fetchSurface(params, 'zos', 1, ocean_state(1)%p_prog(nold(1))%eta_c)
        IF(fg_apply_fillValues) CALL add_missval_2d(patch_3d, ocean_state(1)%p_prog(nold(1))%eta_c)
      ENDIF
      CALL fetchSurface(params, 'stretch_c', 1, ocean_state(1)%p_prog(nold(1))%stretch_c)
      IF(fg_apply_fillValues) CALL add_missval_2d(patch_3d, ocean_state(1)%p_prog(nold(1))%stretch_c)

      !The following variables are read in even though they are in
      !the diagnostic group.

      CALL fetch3d_with_status(routine, 'dwdfg file', params, 'u', 1, ocean_state(1)%p_diag%u, lfound_u)
      IF(lfound_u .AND. fg_apply_fillValues) CALL add_missval_3d(patch_3d, ocean_state(1)%p_diag%u)
      CALL fetch3d_with_status(routine, 'dwdfg file', params, 'v', 1, ocean_state(1)%p_diag%v, lfound_v)
      IF(lfound_u .AND. fg_apply_fillValues) CALL add_missval_3d(patch_3d, ocean_state(1)%p_diag%v)
      CALL fetch3d_with_status(routine, 'dwdfg file', params, 'normal_velocity', 1, ocean_state(1)%p_prog(nold(1))%vn, &
                              & lfound_vn)

      CALL fetch3d_with_status(routine, 'dwdfg file', params, 'tke', 1, params_oce%vmix_params%tke, lfound_tke)
      IF(lfound_tke .AND. fg_apply_fillValues) CALL add_missval_3d(patch_3d, params_oce%vmix_params%tke)
      my_ptr3d => ocean_state(1)%p_diag%SWPT
      CALL fetch3d_with_status(routine, 'dwdfg file', params, 'to', 1, my_ptr3d, lfound_to)

      IF(lfound_to) THEN
        to_var = 'to'
        IF(fg_apply_fillValues) CALL add_missval_3d(patch_3d, ocean_state(1)%p_diag%SWPT)
      ELSE
        CALL fetch3d(params, 'SWPT', 1, my_ptr3d)
        to_var = 'SWPT'
        IF(fg_apply_fillValues) CALL add_missval_3d(patch_3d, ocean_state(1)%p_diag%SWPT)
      ENDIF

      !GRIB files contain temperature in Kelvin but ICON-O wants degrees Celsius, conversion is done with a post-op
      IF(MAXVAL(ocean_state(1)%p_diag%SWPT) >= 100._wp) THEN
        WRITE(0,*) MAXVAL(ocean_state(1)%p_diag%SWPT)
        CALL finish(routine, "Wrong temperature unit. Please check the naming of your temperature field in your input file.")
      ENDIF
      ocean_state(1)%p_prog(nold(1))%tracer(:,:,:,1) = ocean_state(1)%p_diag%SWPT

      my_ptr3d => ocean_state(1)%p_prog(nold(1))%tracer(:,:,:,2)
      CALL fetch3d_with_status(routine, 'dwdfg file', params, 'so', 1, my_ptr3d, lfound_so)
      IF(lfound_so .AND. fg_apply_fillValues) CALL add_missval_3d(patch_3d, my_ptr3d)

    END SUBROUTINE fetch_dwdfg_oce

  !>
  !! Fetch DA-analysis DATA from the request list (ocean only)
  !!
  !! Depending on the initialization mode, either full fields or increments
  !! are read (atmosphere only). The following full fields are read, if available:
  !! u, v, vn, to, so, zos, stretch_c
  !!
  SUBROUTINE fetch_dwdana_oce(patch_3d, requestList, ocean_state, initicono, inputInstructions, read_initicono)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_hydro_ocean_state), INTENT(INOUT), TARGET :: ocean_state(:)
    TYPE(t_initicono_state), INTENT(INOUT), TARGET :: initicono(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)
    TYPE(t_initicono_read) :: read_initicono

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':fetch_dwdana_oce'
    TYPE(t_pi_oce), POINTER :: my_ptr
    REAL(wp), POINTER :: my_ptr3d(:,:,:)
    REAL(wp), POINTER :: my_ptr2d(:,:)
    TYPE(t_fetchParams) :: params
    LOGICAL :: lHaveFg, lfound_to, lfound_so, lfound_u, lfound_v, lfound_vn

    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)))
    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .FALSE.

    ! Depending on the initialization mode chosen (incremental vs. non-incremental)
    ! input fields are stored in different locations.
    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      my_ptr => initicono(1)%oce_inc
    ELSE
      my_ptr => initicono(1)%oce
    ENDIF

    ! start reading DA output (ocean only)
    ! The dynamical variables temp, salinity, u and v can be directly taken from the analysis and are thus written
    ! to the prognostic state
    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar(to_var) == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, to_var, 1, my_ptr%to, lfound_to)
      ! check whether we are using DATA from both FG and ANA input, so that it's correctly listed in the input source table
      IF(lHaveFg .AND. inputInstructions(1)%ptr%sourceOfVar(to_var) == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource(to_var, kInputSourceBoth)
      END IF
      IF(fg_apply_fillValues) CALL add_missval_3d(patch_3d, my_ptr%to)
    ELSE
      my_ptr3d => ocean_state(1)%p_prog(nold(1))%tracer(:,:,:,1)
      CALL fetch3d_with_status(routine, 'dwdana file', params, to_var, 1, my_ptr3d, lfound_to)
      IF(fg_apply_fillValues) CALL add_missval_3d(patch_3d, my_ptr3d)
    ENDIF

    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('so') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'so', 1, my_ptr%so, lfound_so)
      ! check whether we are using DATA from both FG and ANA input, so that it's correctly listed in the input source table
      IF(lHaveFg .AND. inputInstructions(1)%ptr%sourceOfVar('so') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('so', kInputSourceBoth)
      END IF
      IF(lfound_so .AND. fg_apply_fillValues) CALL add_missval_3d(patch_3d, my_ptr%so)
    ELSE
      my_ptr3d => ocean_state(1)%p_prog(nold(1))%tracer(:,:,:,2)
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'so', 1, my_ptr3d, lfound_so)
      IF(fg_apply_fillValues) CALL add_missval_3d(patch_3d, my_ptr3d)
    ENDIF

    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('u') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'u', 1, my_ptr%u, lfound_u)
      !check whether we are using DATA from both FG and ANA input, so that it's correctly listed in the input source table
      IF(lHaveFg .AND. inputInstructions(1)%ptr%sourceOfVar('u') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('u', kInputSourceBoth)
      END IF
      IF(lfound_u .AND. fg_apply_fillValues) CALL add_missval_3d(patch_3d, my_ptr%u)
    ELSE
      my_ptr3d => ocean_state(1)%p_diag%u
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'u', 1, my_ptr3d, lfound_u)
      IF(fg_apply_fillValues) CALL add_missval_3d(patch_3d, my_ptr3d)
    ENDIF

    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('v') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'v', 1, my_ptr%v, lfound_v)
      !check whether we are using DATA from both FG and ANA input, so that it's correctly listed in the input source table
      IF(lHaveFg .AND. inputInstructions(1)%ptr%sourceOfVar('v') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('v', kInputSourceBoth)
      END IF
      IF(lfound_v .AND. fg_apply_fillValues) CALL add_missval_3d(patch_3d, my_ptr%v)
    ELSE
      my_ptr3d => ocean_state(1)%p_diag%v
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'v', 1, my_ptr3d, lfound_v)
      IF(fg_apply_fillValues) CALL add_missval_3d(patch_3d, my_ptr3d)
    ENDIF

    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('normal_velocity') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'normal_velocity', 1, my_ptr%vn, lfound_vn)
      !check whether we are using DATA from both FG and ANA input, so that it's correctly listed in the input source table
      IF(lHaveFg .AND. inputInstructions(1)%ptr%sourceOfVar('normal_velocity') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('normal_velocity', kInputSourceBoth)
      END IF
    ELSE
      my_ptr3d => ocean_state(1)%p_prog(nold(1))%vn
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'normal_velocity', 1, my_ptr3d, lfound_vn)
    ENDIF

    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('zos') == kInputSourceFg
      CALL fetchSurface(params, 'zos', 1, my_ptr%zos)
      ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
      IF(lHaveFg.AND.inputInstructions(1)%ptr%sourceOfVar('zos') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('zos', kInputSourceBoth)
      END IF
    ELSE
      IF(vert_cor_type .eq. 0) THEN
        my_ptr2d => ocean_state(1)%p_prog(nold(1))%h
      ELSEIF(vert_cor_type .eq. 1) THEN
        my_ptr2d => ocean_state(1)%p_prog(nold(1))%eta_c
      ENDIF
      CALL fetchSurface(params, 'zos', 1, my_ptr2d)
    ENDIF

    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('stretch_c') == kInputSourceFg
      CALL fetchSurface(params, 'stretch_c', 1, my_ptr%stretch_c)
      ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
      IF(lHaveFg.AND.inputInstructions(1)%ptr%sourceOfVar('stretch_c') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('stretch_c', kInputSourceBoth)
      END IF
    ELSE
      my_ptr2d => ocean_state(1)%p_prog(nold(1))%stretch_c
      CALL fetchSurface(params, 'stretch_c', 1, my_ptr2d)
    ENDIF

  END SUBROUTINE fetch_dwdana_oce

  !>
  !! Fetch DWD first guess DATA from request list (ocean only)
  SUBROUTINE fetch_dwdfg_seaice(requestList, p_sea_ice, inputInstructions, read_initicono)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_sea_ice),           INTENT(INOUT) :: p_sea_ice
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)
    TYPE(t_initicono_read) :: read_initicono
    TYPE(t_fetchParams) :: params
    LOGICAL :: lfound_hi, lfound_hs, lfound_conc
    INTEGER :: error

    CHARACTER(len=*), PARAMETER :: routine = modname//':fetch_dwdfg_seaice'

    ! Workaround for the intel compiler botching implicit allocation.
    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .TRUE.

    ! sea-ice related fields
    CALL fetch3d_with_status(routine, 'dwdfg file', params, 'hi', 1, p_sea_ice%hi(:,:,:), lfound_hi)
    CALL fetch3d_with_status(routine, 'dwdfg file', params, 'hs', 1, p_sea_ice%hs(:,:,:), lfound_hs)
    CALL fetch3d_with_status(routine, 'dwdfg file', params, 'conc', 1, p_sea_ice%conc(:,:,:), lfound_conc)

  END SUBROUTINE fetch_dwdfg_seaice

  !>
  !! Fetch DWD analysis DATA from the request list (seaice only)
  !!
  !! Analysis is read for:
  !! hi, hs, conc
  !!
  SUBROUTINE fetch_dwdana_seaice(requestList, p_sea_ice, initicono, inputInstructions, read_initicono)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_sea_ice), TARGET, INTENT(INOUT) :: p_sea_ice
    TYPE(t_initicono_state), INTENT(INOUT), TARGET :: initicono(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)
    TYPE(t_pi_seaice), POINTER :: my_ptr
    TYPE(t_initicono_read) :: read_initicono

    TYPE(t_fetchParams) :: params
    LOGICAL :: lHaveFg, lfound_hi, lfound_hs, lfound_conc

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':fetch_dwdana_seaice'

    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)))
    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .FALSE.

    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      my_ptr => initicono(1)%seaice_inc
    ENDIF

    ! hi
    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('hi') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'hi', 1, my_ptr%hi, lfound_hi)
      ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
      IF(lHaveFg .AND. inputInstructions(1)%ptr%sourceOfVar('hi') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('hi', kInputSourceBoth)
      END IF
    ELSE
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'hi', 1, p_sea_ice%hi(:,:,:), lfound_hi)
    ENDIF

    ! hs
    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('hs') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'hs', 1, my_ptr%hs, lfound_hs)
      ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
      IF(lHaveFg .AND. inputInstructions(1)%ptr%sourceOfVar('hs') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('hs', kInputSourceBoth)
      END IF
    ELSE
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'hs', 1, p_sea_ice%hs(:,:,:), lfound_hs)
    ENDIF

    ! conc
    IF ( init_mode_oce == MODE_IAU_OCE ) THEN
      lHaveFg = inputInstructions(1)%ptr%sourceOfVar('conc') == kInputSourceFg
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'conc', 1, my_ptr%conc, lfound_conc)
      ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
      IF(lHaveFg .AND. inputInstructions(1)%ptr%sourceOfVar('conc') == kInputSourceAna) THEN
        CALL inputInstructions(1)%ptr%setSource('conc', kInputSourceBoth)
      END IF
    ELSE
      CALL fetch3d_with_status(routine, 'dwdana file', params, 'conc', 1, p_sea_ice%conc(:,:,:), lfound_conc)
    ENDIF

  END SUBROUTINE fetch_dwdana_seaice

  SUBROUTINE add_missval_3d(patch_3d, var)

   TYPE(t_patch_3d ) :: patch_3d
   REAL(wp), INTENT(INOUT) :: var(:,:,:)
   INTEGER :: blocks, idx, start_cell_index, end_cell_index, level
   TYPE(t_subset_range), POINTER :: all_cells

   CHARACTER(LEN = *), PARAMETER :: routine = modname//':add_missval_3d'

    all_cells => patch_3d%p_patch_2d(1)%cells%ALL

    DO blocks = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blocks, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        DO level = patch_3d%p_patch_1d(1)%dolic_c(idx,blocks) + 1, n_zlev
            var(idx,level,blocks) = fillValue
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE add_missval_3d

  SUBROUTINE add_missval_2d(patch_3d, var)

   TYPE(t_patch_3d ) :: patch_3d
   REAL(wp), INTENT(INOUT) :: var(:,:)
   INTEGER :: blocks, idx, start_cell_index, end_cell_index
   TYPE(t_subset_range), POINTER :: all_cells

   CHARACTER(LEN = *), PARAMETER :: routine = modname//':add_missval2d'

   all_cells => patch_3d%p_patch_2d(1)%cells%ALL

   DO blocks = all_cells%start_block, all_cells%end_block
     CALL get_index_range(all_cells, blocks, start_cell_index, end_cell_index)
     DO idx = start_cell_index, end_cell_index
       IF (patch_3d%p_patch_1d(1)%dolic_c(idx,blocks) .EQ. 0) THEN
         var(idx,blocks) = fillValue
       ENDIF
     ENDDO
   ENDDO

  END SUBROUTINE add_missval_2d


END MODULE mo_ocean_initicono
