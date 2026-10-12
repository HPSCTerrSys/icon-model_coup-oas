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


#include "omp_definitions.inc"
!----------------------------

MODULE mo_ocean_wave_coupling

  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_mpi,                 ONLY: p_comm_work, p_lor
  USE mo_physical_constants,  ONLY: tmelt, rhoh2o
  USE mo_run_config,          ONLY: ltimer, msg_level
  USE mo_dynamics_config,     ONLY: nnew
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling
  USE mo_sync,                ONLY: sync_c, sync_patch_array, global_sum_array
  USE mo_exception,           ONLY: message, message_text
  USE mo_dbg_nml,             ONLY: idbg_mxmn, idbg_val
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE fortran_support,        ONLY: assert_acc_host_only

  !-------------------------------------------------------------
  ! For the coupling
  !
  USE mo_coupling_utils,      ONLY: cpl_def_cell_field_mask, cpl_def_field, &
    &                               cpl_put_field, cpl_get_field
  USE mo_parallel_config,     ONLY: nproma
  USE mo_coupling_config,     ONLY: is_coupled_to_waves

  !-------------------------------------------------------------

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_ocean_wave_coupling'  ! Output of module for 1 line debug

  PUBLIC :: construct_ocean_wave_coupling, couple_ocean_to_waves

  INTEGER, TARGET :: field_id_stokes_u
  INTEGER, TARGET :: field_id_stokes_v
  INTEGER, TARGET :: field_id_tau_w
  INTEGER, TARGET :: field_id_h_s
  INTEGER, TARGET :: field_id_tm02
  INTEGER, TARGET :: field_id_kp
  INTEGER, TARGET :: field_id_tauoc_x
  INTEGER, TARGET :: field_id_tauoc_y
  INTEGER, TARGET :: field_id_phioc
  INTEGER, TARGET :: field_id_oce_u
  INTEGER, TARGET :: field_id_oce_v
  INTEGER, TARGET :: field_id_ssh
  INTEGER, TARGET :: field_id_ssd

  INTEGER, SAVE :: nbr_inner_cells

CONTAINS

  !--------------------------------------------------------------------------
  !>
  !! Registers fields required for the coupling between ocean and
  !! wave
  !!
  !! This subroutine is called from construct_ocean_coupling.
  !!
  SUBROUTINE construct_ocean_wave_coupling( &
    patch_3d, comp_id, grid_id, cell_point_id, timestepstring, &
    ocean_grid_name)

    TYPE(t_patch_3d ), TARGET, INTENT(in) :: patch_3d
    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: grid_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring
    CHARACTER(LEN=*), INTENT(IN) :: ocean_grid_name

    INTEGER                :: patch_no
    TYPE(t_patch), POINTER :: patch_horz

    INTEGER :: cell_mask_id

!    INTEGER :: mask_checksum
    INTEGER :: blockNo, cell_index, i
    LOGICAL :: use_mask

    LOGICAL, ALLOCATABLE  :: is_valid(:)

    TYPE field_id_ptr
      INTEGER, POINTER :: p
    END TYPE field_id_ptr

    INTEGER, PARAMETER :: no_of_fields = 13
    CHARACTER(LEN=max_char_length) :: field_name(no_of_fields)
    INTEGER                        :: collection_size(no_of_fields)
    TYPE(field_id_ptr)             :: field_ids(no_of_fields)

    CHARACTER(LEN=*), PARAMETER   :: &
      routine = str_module // ':construct_ocean_wave_coupling'

    patch_no = 1
    patch_horz => patch_3d%p_patch_2d(patch_no)

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Receive fields from waves to ocean
    !  "zonal_stokes_drift"
    !  "meridional_stokes_drift"
    !  "wave_stress"
    !  "significant_wave_height"
    !  "m2_wave_period"
    !  "peak_wavenumber"
    !  "zonal_wave_to_ocean_stress"
    !  "meridional_wave_to_ocean_stress"
    !  "wave_to_ocean_energy_flux"
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from ocean to waves
    !  "zonal_sea_surface_current"
    !  "meridional_sea_surface_current"
    !  "sea_surface_height"
    !  "sea_surface_water_density"
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    field_name(1) = "zonal_stokes_drift"
    collection_size(1) = 1
    field_ids(1)%p => field_id_stokes_u
    field_name(2) = "meridional_stokes_drift"
    collection_size(2) = 1
    field_ids(2)%p => field_id_stokes_v
    field_name(3) = "wave_stress"
    collection_size(3) = 1
    field_ids(3)%p => field_id_tau_w
    field_name(4) = "significant_wave_height"
    collection_size(4) = 1
    field_ids(4)%p => field_id_h_s
    field_name(5) = "m2_wave_period"
    collection_size(5) = 1
    field_ids(5)%p => field_id_tm02
    field_name(6) = "peak_wavenumber"
    collection_size(6) = 1
    field_ids(6)%p => field_id_kp
    field_name(7) = "zonal_wave_to_ocean_stress"
    collection_size(7) = 1
    field_ids(7)%p => field_id_tauoc_x
    field_name(8) = "meridional_wave_to_ocean_stress"
    collection_size(8) = 1
    field_ids(8)%p => field_id_tauoc_y
    field_name(9) = "wave_to_ocean_energy_flux"
    collection_size(9) = 1
    field_ids(9)%p => field_id_phioc
    field_name(10) = "zonal_sea_surface_current"
    collection_size(10) = 1
    field_ids(10)%p => field_id_oce_u
    field_name(11) = "meridional_sea_surface_current"
    collection_size(11) = 1
    field_ids(11)%p => field_id_oce_v
    field_name(12) = "sea_surface_height"
    collection_size(12) = 1
    field_ids(12)%p => field_id_ssh
    field_name(13) =   "sea_surface_water_density"
    collection_size(13) = 1
    field_ids(13)%p => field_id_ssd

    !
    ! mask generation : ... not yet defined ...
    !
    ! We could use the patch_horz%cells%decomp_info%owner_local information
    ! e.g. to mask out halo points. We do we get the info about what is local and what
    ! is remote.
    !
    ! The integer land-sea mask:
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !
    ! This integer mask for the ocean is available in patch_3D%surface_cell_sea_land_mask(:,:)
    ! The logical mask for the coupler is set to .FALSE. for land points to exclude them from mapping by yac.
    ! These points are not touched by yac.

    use_mask = .FALSE.

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,cell_index, use_mask) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_c
      DO cell_index = 1, nproma
        use_mask = use_mask .AND. (patch_3d%surface_cell_sea_land_mask(cell_index, blockNo) /= 0)
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

    use_mask = p_lor(use_mask, comm=p_comm_work)

    IF ( use_mask ) THEN

      ALLOCATE(is_valid(nproma*patch_horz%nblks_c))

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        DO cell_index = 1, nproma
          IF ( patch_3d%surface_cell_sea_land_mask(cell_index, blockNo) < 0 ) THEN
            ! ocean and ocean-coast is valid (-2, -1)
            is_valid((blockNo-1)*nproma+cell_index) = .TRUE.
          ELSE
            ! land is undef (1, 2)
            is_valid((blockNo-1)*nproma+cell_index) = .FALSE.
          END IF
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

      CALL cpl_def_cell_field_mask(routine, grid_id, is_valid, cell_mask_id)

      DO i = 1, no_of_fields

        CALL cpl_def_field( &
          comp_id, cell_point_id, cell_mask_id, timestepstring, &
          field_name(i), collection_size(i), field_ids(i)%p)

      END DO

    ELSE ! no mask

      DO i = 1, no_of_fields

        CALL cpl_def_field( &
          comp_id, cell_point_id, timestepstring, &
          TRIM(field_name(i)), collection_size(i), field_ids(i)%p)

      END DO

    END IF

  END SUBROUTINE construct_ocean_wave_coupling

  !>
  !! Exchange fields between ocean and waves model
  !!
  SUBROUTINE couple_ocean_to_waves(patch_3d, ice_conc, cur_u, cur_v, ssh, ssd, u2d_stokes, v2d_stokes, tau_w, swh, Tm2, kp, &
                                  & tauoc_x, tauoc_y, phioc, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(IN)        :: patch_3d
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: ice_conc(:,:)   ! ice concentration (first ice class) (1)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: cur_u(:,:)      ! zonal sea surface current [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: cur_v(:,:)      ! meridional sea surface current [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: ssh(:,:)        ! sea surface height [m]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: ssd(:,:)        ! sea surface water density [kg/m3]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: u2d_stokes(:,:) ! zonal surface Stokes velocity component from surface waves [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: v2d_stokes(:,:) ! meridional surface Stokes velocity component from surface waves [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: tau_w(:,:)      ! wave stress [m2/s2]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: swh(:,:)        ! significant wave height [m]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: Tm2(:,:)        ! m2 wave period [s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: kp(:,:)         ! peak wavenumber [1/m]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: tauoc_x(:,:)    ! zonal wave-to-ocean stress [m2/s2]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: tauoc_y(:,:)    ! meridional wave-to-ocean stress [m2/s2]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: phioc(:,:)      ! wave-to-ocean energy flux [kg/s3]

    LOGICAL, OPTIONAL, INTENT(IN) :: lacc ! If true, use openacc

    ! Local declarations for coupling:
    INTEGER :: nbr_hor_cells  ! = inner and halo points
    INTEGER :: cell_index     ! nproma loop count
    INTEGER :: nn             ! block offset
    INTEGER :: blockNo        ! block loop count
    INTEGER :: nlen           ! nproma/npromz
    INTEGER :: nblks_c        ! number of blocks
    INTEGER :: no_arr         ! no of arrays in bundle for put/get calls
    TYPE(t_patch), POINTER:: patch_horz

    REAL(wp), PARAMETER :: dummy = 0.0_wp

    REAL(wp) :: diag_tmp
    LOGICAL  :: received_data

    CHARACTER(LEN=*), PARAMETER   :: routine = str_module // ':couple_ocean_to_waves'

    IF (.NOT. is_coupled_to_waves() ) RETURN

    CALL assert_acc_host_only(routine, lacc)

    patch_horz   => patch_3D%p_patch_2D(1)

    nbr_hor_cells = patch_horz%n_patch_cells
    nblks_c = patch_horz%nblks_c

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Receive fields from waves to ocean
    !  "zonal_stokes_drift"
    !  "meridional_stokes_drift"
    !  "wave_stress"
    !  "significant_wave_height"
    !  "m2_wave_period"
    !  "peak_wavenumber"
    !  "zonal_wave_to_ocean_stress"
    !  "meridional_wave_to_ocean_stress"
    !  "wave_to_ocean_energy_flux"
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from ocean to waves
    !  "zonal_sea_surface_current"
    !  "meridional_sea_surface_current"
    !  "sea_surface_height"
    !  "sea_surface_water_density"
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !

    !  Send fields from ocean to waves
    !   "zonal_sea_surface_current"
    CALL cpl_put_field(routine, field_id_oce_u, 'zonal_sea_surface_current', nbr_hor_cells, cur_u)
    !   "meridional_sea_surface_current"
    CALL cpl_put_field(routine, field_id_oce_v, 'meridional_sea_surface_current', nbr_hor_cells, cur_v)
    !   "sea_surface_height"
    CALL cpl_put_field(routine, field_id_ssh, 'sea_surface_height', nbr_hor_cells, ssh)
    !   "sea_surface_water_density"
    CALL cpl_put_field(routine, field_id_ssd, 'sea_surface_water_density', nbr_hor_cells, ssd)


    !   Receive fields from waves
    !   "zonal_stokes_drift"
    CALL cpl_get_field(routine, field_id_stokes_u, 'zonal_stokes_drift', &
      nbr_hor_cells, u2d_stokes, first_get=.TRUE., received_data=received_data)
    IF (received_data) CALL sync_patch_array(sync_c, patch_horz, u2d_stokes, lacc)

    !   "meridional_stokes_drift"
    CALL cpl_get_field(routine, field_id_stokes_v, 'meridional_stokes_drift', &
      nbr_hor_cells, v2d_stokes, first_get=.TRUE., received_data=received_data)
    IF (received_data) CALL sync_patch_array(sync_c, patch_horz, v2d_stokes, lacc)

    !   "wave_stress"
    CALL cpl_get_field(routine, field_id_tau_w, 'wave_stress', &
      nbr_hor_cells, tau_w, first_get=.TRUE., received_data=received_data)
    IF (received_data) CALL sync_patch_array(sync_c, patch_horz, tau_w, lacc)

    !   "significant_wave_height"
    CALL cpl_get_field(routine, field_id_h_s, 'significant_wave_height', &
      nbr_hor_cells, swh, first_get=.TRUE., received_data=received_data)
    IF (received_data) CALL sync_patch_array(sync_c, patch_horz, swh, lacc)

    !   "m2_wave_period"
    CALL cpl_get_field(routine, field_id_tm02, 'm2_wave_period', &
      nbr_hor_cells, Tm2, first_get=.TRUE., received_data=received_data)
    IF (received_data) CALL sync_patch_array(sync_c, patch_horz, Tm2, lacc)

    !   "peak_wavenumber"
    CALL cpl_get_field(routine, field_id_kp, 'peak_wavenumber', &
      nbr_hor_cells, kp, first_get=.TRUE., received_data=received_data)
    IF (received_data) CALL sync_patch_array(sync_c, patch_horz, kp, lacc)

    !   "zonal_wave_to_ocean_stress"
    CALL cpl_get_field(routine, field_id_tauoc_x, 'zonal_wave_to_ocean_stress', &
      nbr_hor_cells, tauoc_x, first_get=.TRUE., received_data=received_data)
    IF (received_data) CALL sync_patch_array(sync_c, patch_horz, tauoc_x, lacc)

    !   "meridional_wave_to_ocean_stress"
    CALL cpl_get_field(routine, field_id_tauoc_y, 'meridional_wave_to_ocean_stress', &
      nbr_hor_cells, tauoc_y, first_get=.TRUE., received_data=received_data)
    IF (received_data) CALL sync_patch_array(sync_c, patch_horz, tauoc_y, lacc)

    !   "wave_to_ocean_energy_flux"
    CALL cpl_get_field(routine, field_id_phioc, 'wave_to_ocean_energy_flux', &
      nbr_hor_cells, phioc, first_get=.TRUE., received_data=received_data)
    IF (received_data) CALL sync_patch_array(sync_c, patch_horz, phioc, lacc)


  END SUBROUTINE couple_ocean_to_waves
  !--------------------------------------------------------------------------

END MODULE mo_ocean_wave_coupling
