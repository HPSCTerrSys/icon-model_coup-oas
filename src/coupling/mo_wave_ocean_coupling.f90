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

! Interface between ocean surface waves and atmosphere, through a coupler

MODULE mo_wave_ocean_coupling

  USE mo_kind,           ONLY: wp
  USE mo_model_domain,   ONLY: t_patch
  USE mo_coupling_utils, ONLY: cpl_def_field, cpl_put_field, cpl_get_field
  USE mo_sync,           ONLY: SYNC_C, sync_patch_array
  USE fortran_support,   ONLY: assert_acc_host_only

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_wave_ocean_coupling_post_sync, couple_wave_to_ocean

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_ocean_coupling'

  INTEGER :: field_id_stokes_u
  INTEGER :: field_id_stokes_v
  INTEGER :: field_id_tau_w
  INTEGER :: field_id_h_s
  INTEGER :: field_id_tm02
  INTEGER :: field_id_kp
  INTEGER :: field_id_tauoc_x
  INTEGER :: field_id_tauoc_y
  INTEGER :: field_id_phioc
  INTEGER :: field_id_cur_u
  INTEGER :: field_id_cur_v
  INTEGER :: field_id_ssh
  INTEGER :: field_id_ssd

CONTAINS

  !>
  !! Registers fields required for the coupling between wave and ocean
  !!
  !! This subroutine is called from construct_wave_coupling.
  !!
  SUBROUTINE construct_wave_ocean_coupling_post_sync( &
    comp_id, cell_point_id, timestepstring, wave_comp_name, wave_grid_name)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring
    CHARACTER(LEN=*), INTENT(IN) :: wave_comp_name
    CHARACTER(LEN=*), INTENT(IN) :: wave_grid_name

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "zonal_stokes_drift", 1, field_id_stokes_u)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "meridional_stokes_drift", 1, field_id_stokes_v)

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "wave_stress", 1, field_id_tau_w)

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "significant_wave_height",1, field_id_h_s)

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "m2_wave_period",1,field_id_tm02)

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "peak_wavenumber", 1, field_id_kp)

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "zonal_wave_to_ocean_stress", 1, field_id_tauoc_x)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "meridional_wave_to_ocean_stress", 1, field_id_tauoc_y)

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "wave_to_ocean_energy_flux", 1, field_id_phioc)

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "zonal_sea_surface_current", 1, field_id_cur_u)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "meridional_sea_surface_current", 1, field_id_cur_v)

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "sea_surface_height", 1, field_id_ssh)

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "sea_surface_water_density", 1, field_id_ssd)

  END SUBROUTINE construct_wave_ocean_coupling_post_sync

  !>
  !! Exchange fields between the wave model and the ocean model
  !!
  !! Send fields to ocean:
  !!  "zonal_stokes_drift"
  !!  "meridional_stokes_drift"
  !!  "wave_stress"
  !!  "significant_wave_height"
  !!  "m2_wave_period"
  !!  "peak_wavenumber"
  !!  "zonal_wave_to_ocean_stress"
  !!  "meridional_wave_to_ocean_stress"
  !!  "wave_to_ocean_energy_flux"
  !!
  !! Receive fields from ocean:
  !!  "zonal_sea_surface_current"
  !!  "meridional_sea_surface_current"
  !!  "sea_surface_height"
  !!  "sea_surface_water_density"
  !!
  !! This subroutine is called from perform_wave_stepping.
  !!
  SUBROUTINE couple_wave_to_ocean(p_patch, stokes_u, stokes_v, tau_w, h_s, tm02, kp, tauoc_x, tauoc_y, phioc, cur_u, cur_v, ssh, ssd, lacc)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':couple_wave_to_ocean'

    TYPE(t_patch),                INTENT(IN)    :: p_patch
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: stokes_u(:,:) ! zonal stokes drift [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: stokes_v(:,:) ! meridional stokes drift [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: tau_w(:,:)    ! wave stress [m2/s2]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: h_s(:,:)      ! significant wave height [m]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: tm02(:,:)     ! m2 wave period [s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: kp(:,:)       ! total peak wavenumber [1/m]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: tauoc_x(:,:)  ! zonal wave-to-ocean stress [m2/s2]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: tauoc_y(:,:)  ! meridional wave-to-ocean stress [m2/s2]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: phioc(:,:)    ! wave-to-ocean energy flux [kg/s3]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: cur_u(:,:)    ! zonal sea surface current [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: cur_v(:,:)    ! meridional sea surface current [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: ssh(:,:)      ! sea surface height [m]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: ssd(:,:)      ! sea surface water density [kg/m3]
    LOGICAL,  OPTIONAL,           INTENT(IN)    :: lacc          ! if true, use openacc

    LOGICAL :: received_data

    CALL assert_acc_host_only(routine, lacc)

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from waves to ocean
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

    ! --------------------------------------------------
    !  Send Stokes drift components to the ocean
    !  'zonal_stokes_drift'
    !  'meridional_stokes_drift'
    ! --------------------------------------------------
    !
    CALL cpl_put_field( &
      routine, field_id_stokes_u, 'zonal_stokes_drift', p_patch%n_patch_cells, stokes_u)

    CALL cpl_put_field( &
      routine, field_id_stokes_v, 'meridional_stokes_drift', p_patch%n_patch_cells, stokes_v)

    ! --------------------------------------------------
    !  Send wave stress to the ocean
    !  'wave_stress'
    ! --------------------------------------------------
    !
    CALL cpl_put_field( &
      routine, field_id_tau_w, 'wave_stress', p_patch%n_patch_cells, tau_w)

    ! --------------------------------------------------
    !  Send significant wave height to the ocean
    !  'significant_wave_height'
    ! --------------------------------------------------
    !
    CALL cpl_put_field( &
      routine, field_id_h_s, 'significant_wave_height', p_patch%n_patch_cells, h_s)

    ! --------------------------------------------------
    !  Send Tm02 period to the ocean
    !  'm2_wave_period'
    ! --------------------------------------------------
    !
    CALL cpl_put_field( &
      routine, field_id_tm02, 'm2_wave_period', p_patch%n_patch_cells, tm02)

    ! --------------------------------------------------
    !  Send peak wavenumber to the ocean
    !  'peak_wavenumber'
    ! --------------------------------------------------
    !
    CALL cpl_put_field( &
      routine, field_id_kp, 'peak_wavenumber', p_patch%n_patch_cells, kp)

    ! --------------------------------------------------
    !  Send wave-to-ocean stress components to the ocean
    !  'zonal_wave_to_ocean_stress'
    !  'meridional_wave_to_ocean_stress'
    ! --------------------------------------------------
    !
    CALL cpl_put_field( &
      routine, field_id_tauoc_x, 'zonal_wave_to_ocean_stress', p_patch%n_patch_cells, tauoc_x)

    CALL cpl_put_field( &
      routine, field_id_tauoc_y, 'meridional_wave_to_ocean_stress', p_patch%n_patch_cells, tauoc_y)

    ! --------------------------------------------------
    !  Send wave-to-ocean energy flux to the ocean
    !  'wave_to_ocean_energy_flux'
    ! --------------------------------------------------
    !
    CALL cpl_put_field( &
      routine, field_id_phioc, 'wave_to_ocean_energy_flux', p_patch%n_patch_cells, phioc)


    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !   Receive fields from ocean
    !  "zonal_sea_surface_current"
    !  "meridional_sea_surface_current"
    !  "sea_surface_height"
    !  "sea_surface_water_density"
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    ! -----------------------------------------
    !  Receive surface currents
    !  'zonal_sea_surface_current'
    !  'meridional_sea_surface_current'
    ! -----------------------------------------
    !
    CALL cpl_get_field( &
      routine, field_id_cur_u, 'zonal_sea_surface_current', p_patch%n_patch_cells, cur_u, &
      first_get=.TRUE., received_data=received_data)
    IF (received_data) &
      CALL sync_patch_array( &
        SYNC_C, p_patch, cur_u, opt_varname='zonal_sea_surface_current', lacc=lacc)

    CALL cpl_get_field( &
      routine, field_id_cur_v, 'meridional_sea_surface_current', p_patch%n_patch_cells, cur_v, &
      first_get=.TRUE., received_data=received_data)
    IF (received_data) &
      CALL sync_patch_array( &
        SYNC_C, p_patch, cur_v, opt_varname='meridional_sea_surface_current', lacc=lacc)

    ! ----------------------------------------------
    !  Receive sea surface height
    !  'sea_surface_height'
    ! ----------------------------------------------
    !
    CALL cpl_get_field( &
      routine, field_id_ssh, 'sea_surface_height', p_patch%n_patch_cells, ssh, &
      received_data=received_data)
    IF (received_data) &
      CALL sync_patch_array( &
        SYNC_C, p_patch, ssh, opt_varname='sea_surface_height', lacc=lacc)

    ! ----------------------------------------------
    !  Receive sea surface density
    !  'sea_surface_water_density'
    ! ----------------------------------------------
    !
    CALL cpl_get_field( &
      routine, field_id_ssd, 'sea_surface_water_density', p_patch%n_patch_cells, ssd, &
      received_data=received_data)
    IF (received_data) &
      CALL sync_patch_array( &
        SYNC_C, p_patch, ssd, opt_varname='sea_surface_water_density', lacc=lacc)

  END SUBROUTINE couple_wave_to_ocean

END MODULE mo_wave_ocean_coupling
