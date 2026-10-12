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

MODULE mo_ocean_era5_provider_coupling

  USE mo_kind,            ONLY: wp
  USE mo_parallel_config, ONLY: nproma
  USE mo_exception,       ONLY: finish
  USE mo_impl_constants,  ONLY: max_char_length
  USE mo_coupling_utils,  ONLY: cpl_def_cell_field_mask, cpl_def_field, &
    &                           cpl_get_field, cpl_update_field, &
    &                           cpl_get_field_action
  USE mo_model_domain,    ONLY: t_patch, t_patch_3d
  USE mo_mpi,             ONLY: p_comm_work, p_lor
  USE mo_sync,            ONLY: sync_c, sync_patch_array, sync_patch_array_mult
  USE mo_fortran_tools,   ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_ocean_era5_provider_coupling' ! Output of module for debug

  PUBLIC :: construct_ocean_era5_provider_coupling_post_sync, &
            couple_ocean_to_era5_provider

!  INTEGER, PARAMETER :: no_of_fields = 13
  INTEGER, PARAMETER :: no_of_fields = 12

  TYPE :: field_data
    INTEGER :: id
    CHARACTER(len=max_char_length) :: name
    REAL(wp), POINTER :: data(:,:)
  END TYPE field_data

  TYPE :: field_data_ptr
    TYPE(field_data), POINTER :: p
  END TYPE field_data_ptr

!  TYPE(field_data), TARGET :: field_10m_wind_speed
  TYPE(field_data), TARGET :: field_ustress
  TYPE(field_data), TARGET :: field_vstress
  TYPE(field_data), TARGET :: field_u10
  TYPE(field_data), TARGET :: field_v10
  TYPE(field_data), TARGET :: field_ldown
  TYPE(field_data), TARGET :: field_swdown
  TYPE(field_data), TARGET :: field_precip
  TYPE(field_data), TARGET :: field_runoff
  TYPE(field_data), TARGET :: field_slp
  TYPE(field_data), TARGET :: field_t2m
  TYPE(field_data), TARGET :: field_tcc
  TYPE(field_data), TARGET :: field_tdew
  TYPE(field_data_ptr) :: fields(no_of_fields)

CONTAINS

  SUBROUTINE construct_ocean_era5_provider_coupling_post_sync( &
    patch_3d, comp_id, grid_id, cell_point_id, timestepstring)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: grid_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    TYPE(t_patch_3d ), TARGET, INTENT(in) :: patch_3d
    TYPE(t_patch), POINTER :: patch_horz

    INTEGER :: patch_no
    INTEGER :: cell_mask_id

    LOGICAL :: use_mask
    INTEGER :: block_index, cell_index, i

    LOGICAL, ALLOCATABLE  :: is_valid(:)

    TYPE field_id_ptr
      INTEGER, POINTER :: p
    END TYPE field_id_ptr


    CHARACTER(LEN=max_char_length) :: field_name(no_of_fields)
    TYPE(field_id_ptr)             :: field_ids(no_of_fields)

    CHARACTER(LEN=*), PARAMETER   :: &
      routine = str_module // ':construct_ocean_era5_provider_coupling_post_sync'

    patch_no = 1
    patch_horz => patch_3d%p_patch_2d(patch_no)


    ! The integer land-sea mask:
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !
    ! This integer mask for the ocean is available in patch_3D%surface_cell_sea_land_mask(:,:)
    ! The logical mask for the coupler is set to .FALSE. for land points to exclude them from mapping by yac.
    ! These points are not touched by yac.

    use_mask = &
      p_lor( &
        ANY(patch_3d%surface_cell_sea_land_mask( &
              1:nproma,1:patch_horz%nblks_c) /= 0.0), p_comm_work)

    IF ( use_mask ) THEN

      ALLOCATE(is_valid(nproma*patch_horz%nblks_c))

!ICON_OMP_PARALLEL_DO PRIVATE(block_index, cell_index) ICON_OMP_DEFAULT_SCHEDULE
      DO block_index = 1, patch_horz%nblks_c
        DO cell_index = 1, nproma
          IF ( patch_3d%surface_cell_sea_land_mask(cell_index, block_index) < 0 ) THEN
            ! ocean and ocean-coast is valid (-2, -1)
            is_valid((block_index-1)*nproma+cell_index) = .TRUE.
          ELSE
            ! land is undef (1, 2)
            is_valid((block_index-1)*nproma+cell_index) = .FALSE.
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

      CALL cpl_def_cell_field_mask(routine, grid_id, is_valid, cell_mask_id)

      DEALLOCATE(is_valid)

    ELSE
      cell_mask_id = -1

    ENDIF

    ! Define the fields
!    field_10m_wind_speed = def_field("10m_wind_speed")
    field_ustress = def_field("ustress")
    field_vstress = def_field("vstress")
    field_ldown = def_field("ldown")
    field_swdown = def_field("swdown")
    field_precip = def_field("precip")
    field_runoff = def_field("river_runoff")
    field_slp = def_field("sea_level_pressure")
    field_t2m = def_field("t2m")
    field_tcc = def_field("tcc")
    field_tdew = def_field("tdew")
    field_u10 = def_field("u10")
    field_v10 = def_field("v10")

    ! Store field pointers in array for easier access
    ! (for improved performance the order of these fields should match the
    !  order in which they are sent by the provider)
    fields(1)%p=> field_u10
    fields(2)%p=> field_v10
    fields(3)%p=> field_ustress
    fields(4)%p=> field_vstress
    fields(5)%p=> field_ldown
    fields(6)%p=> field_swdown
    fields(7)%p=> field_precip
    fields(8)%p=> field_slp
    fields(9)%p=> field_t2m
    fields(10)%p=> field_tcc
    fields(11)%p=> field_tdew
!    fields(12)%p=> field_10m_wind_speed
    fields(12)%p=> field_runoff

  CONTAINS

    FUNCTION def_field(field_name) RESULT(field)

      CHARACTER(LEN=*), INTENT(IN) :: field_name
      TYPE(field_data) :: field

      INTEGER, PARAMETER :: collection_size = 1

      IF (use_mask) THEN
        CALL cpl_def_field( &
          comp_id, cell_point_id, cell_mask_id, timestepstring, &
          TRIM(field_name), collection_size, field%id)
      ELSE
        CALL cpl_def_field( &
          comp_id, cell_point_id, timestepstring, &
          TRIM(field_name), collection_size, field%id)
      END IF

      field%name = field_name
      field%data => NULL()

    END FUNCTION def_field

  END SUBROUTINE construct_ocean_era5_provider_coupling_post_sync

  !>
  !! Receives fields from the era5 provider in the ocean model
  !!
  SUBROUTINE couple_ocean_to_era5_provider(p_patch,  era5_t2m, era5_tdew, &
    era5_wind10, era5_ustress, era5_vstress, era5_ldown, era5_swdown, &
    era5_precip, era5_runoff, era5_slp, era5_tcc, era5_u10, era5_v10, lacc)

    TYPE(t_patch_3d), INTENT(in) :: p_patch
    TYPE(t_patch), POINTER :: patch_horz

    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    REAL(wp), TARGET, INTENT(inout) :: era5_t2m(:,:)
    REAL(wp), TARGET, INTENT(inout) :: era5_tdew(:,:)
    REAL(wp), TARGET, INTENT(inout) :: era5_wind10(:,:)
    REAL(wp), TARGET, INTENT(inout) :: era5_ustress(:,:)
    REAL(wp), TARGET, INTENT(inout) :: era5_vstress(:,:)
    REAL(wp), TARGET, INTENT(inout) :: era5_ldown(:,:)
    REAL(wp), TARGET, INTENT(inout) :: era5_swdown(:,:)
    REAL(wp), TARGET, INTENT(inout) :: era5_precip(:,:)
    REAL(wp), TARGET, INTENT(inout) :: era5_runoff(:,:)
    REAL(wp), TARGET, INTENT(inout) :: era5_slp(:,:)
    REAL(wp), TARGET, INTENT(inout) :: era5_tcc(:,:)
    REAL(wp), TARGET, INTENT(inout) :: era5_u10(:,:)
    REAL(wp), TARGET, INTENT(inout) :: era5_v10(:,:)

    LOGICAL :: received_data_array(no_of_fields)
    LOGICAL :: received_data
    LOGICAL :: lzacc

    ! Local declarations for coupling:
    INTEGER :: cell_index     ! nproma loop count
    INTEGER :: block_index    ! block loop count
    INTEGER :: field_idx      ! field loop count

    REAL(wp), PARAMETER :: dummy = 0.0_wp
    REAL(wp), ALLOCATABLE :: buffer(:,:,:)

    CHARACTER(LEN=*), PARAMETER   :: &
      routine = str_module // ':couple_ocean_to_era5_provider'

    CALL set_acc_host_or_device(lzacc, lacc)

    patch_horz => p_patch%p_patch_2d(1)

    ! check the action for the next get calls
    DO field_idx = 1, no_of_fields
      CALL cpl_get_field_action( &
          routine, fields(field_idx)%p%id, received_data_array(field_idx))
    END DO

    ! check whether this is a coupling step (all fields received data)
    received_data = ALL(received_data_array(:))

    ! check whether only some of the fields received data -> error
    IF (ANY(received_data_array(:)) .AND. .NOT. received_data) THEN
      CALL finish( &
        routine, &
        "Inconsistent data reception from ERA5 provider: only some fields received data.")
    END IF

    ! if this is a coupling step
    IF (received_data) THEN

      ! update host copies of the era5 arrays before using them

      !$ACC UPDATE SELF(era5_t2m, era5_tdew, era5_wind10, era5_ustress) &
      !$ACC   SELF(era5_vstress, era5_ldown, era5_swdown, era5_precip) &
      !$ACC   SELF(era5_runoff, era5_slp, era5_tcc, era5_u10, era5_v10) IF(lzacc)

      ! set pointers to data arrays
      field_t2m%data => era5_t2m
      field_tdew%data => era5_tdew
      field_ustress%data => era5_ustress
      field_vstress%data => era5_vstress
      field_u10%data => era5_u10
      field_v10%data => era5_v10
      field_ldown%data => era5_ldown
      field_swdown%data => era5_swdown
      field_precip%data => era5_precip
      field_runoff%data => era5_runoff
      field_slp%data => era5_slp
      field_tcc%data => era5_tcc

      ! receive data from era5 provider
      DO field_idx = 1, no_of_fields
        ! only the first get call needs first_get=.TRUE.
        ! (switches internal timer from cpl_get to cpl_1st_get)
        CALL cpl_get_field( &
            routine, fields(field_idx)%p%id, fields(field_idx)%p%name, &
            patch_horz%n_patch_cells, fields(field_idx)%p%data, &
            first_get=(field_idx == 1))
      END DO

      ALLOCATE(buffer(nproma,no_of_fields,patch_horz%alloc_cell_blocks))

      ! store all received fields in buffer for halo exchange
      ! (this allows the halo exchange of all fields in one call)
      DO field_idx = 1, no_of_fields
!ICON_OMP_PARALLEL_DO PRIVATE(block_index, cell_index) ICON_OMP_DEFAULT_SCHEDULE
        DO block_index = 1, patch_horz%nblks_c
          DO cell_index = 1, nproma
            buffer(cell_index,field_idx,block_index) = &
              fields(field_idx)%p%data(cell_index,block_index)
          ENDDO
        ENDDO
!ICON_OMP_END_PARALLEL_DO
      END DO

      ! the last block may be only partially filled with valid data
      ! initialise remaining block with dummy values
      buffer(patch_horz%npromz_c+1:nproma, :, patch_horz%nblks_c) = dummy

      ! do a halo exchange of all received fields
      CALL sync_patch_array(sync_c, patch_horz, buffer, lacc=.FALSE.)

      ! copy back halo-exchanged data into era5 arrays
      DO field_idx = 1, no_of_fields
!ICON_OMP_PARALLEL_DO PRIVATE(block_index, cell_index) ICON_OMP_DEFAULT_SCHEDULE
        DO block_index = 1, patch_horz%nblks_c
          DO cell_index = 1, nproma
            fields(field_idx)%p%data(cell_index,block_index) = &
              buffer(cell_index,field_idx,block_index)
          ENDDO
        ENDDO
!ICON_OMP_END_PARALLEL_DO
      END DO

      ! calculate 10m wind speed
      era5_wind10(:,:) = SQRT(era5_u10(:,:)**2 + era5_v10(:,:)**2)

      !$ACC UPDATE DEVICE(era5_t2m, era5_tdew, era5_wind10, era5_ustress) &
      !$ACC   DEVICE(era5_vstress, era5_ldown, era5_swdown, era5_precip) &
      !$ACC   DEVICE(era5_runoff, era5_slp, era5_tcc, era5_u10, era5_v10) IF(lzacc)

    ELSE ! not received_data

      ! update internal clock of all fields instead of call get
      DO field_idx = 1, no_of_fields
        CALL cpl_update_field(routine, fields(field_idx)%p%id)
      END DO

    END IF ! received_data

  END SUBROUTINE couple_ocean_to_era5_provider

END MODULE mo_ocean_era5_provider_coupling
