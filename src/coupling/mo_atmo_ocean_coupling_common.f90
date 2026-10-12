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
! @brief Common field definition for the atmosphere-ocean coupling


!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_atmo_ocean_coupling_common

  USE mo_coupling_utils,  ONLY: cpl_def_field, cpl_get_field_datetime, &
                                cpl_get_field_is_source, cpl_get_field_name
  USE mo_grid_config,     ONLY: n_dom
  USE mtime,              ONLY: datetime, OPERATOR(==), OPERATOR(/=)
  USE mo_exception,       ONLY: finish
  USE mo_ccycle_config,   ONLY: CCYCLE_MODE_INTERACTIVE, ccycle_config

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_atmo_ocean_coupling_common' ! Output of module for debug

  PUBLIC :: construct_atmo_ocean_coupling_common
  PUBLIC :: construct_atmo_ocean_coupling_common_finalize
  PUBLIC :: destruct_atmo_ocean_coupling_common
  PUBLIC :: t_out_field_id, t_in_field_id
  PUBLIC :: out_field_ids, in_field_ids

  TYPE t_out_field_id
    INTEGER :: umfl
    INTEGER :: vmfl
    INTEGER :: freshflx
    INTEGER :: heatflx
    INTEGER :: seaice_atm
    INTEGER :: sp10m
    INTEGER :: co2_vmr
    INTEGER :: pres_msl
  END TYPE t_out_field_id

  TYPE t_in_field_id
    INTEGER :: sst
    INTEGER :: seaice_oce
    INTEGER :: surface_velocity
    INTEGER :: co2_flx
  END TYPE t_in_field_id

  TYPE(t_out_field_id) :: out_field_ids
  TYPE(t_in_field_id), ALLOCATABLE :: in_field_ids(:)

CONTAINS

  !>
  !! Registers fields required for the coupling between atmo and
  !! ocean
  !!
  !! This subroutine is called from construct_atmo_ocean_coupling.
  !!
  SUBROUTINE construct_atmo_ocean_coupling_common( &
    comp_id, cell_point_id, cell_mask_id, timestepstring, &
    use_ocean_velocity)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id(0:n_dom)
    INTEGER, INTENT(IN) :: cell_mask_id(0:n_dom)
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring
    LOGICAL, INTENT(IN) :: use_ocean_velocity

    INTEGER :: jg

    ! define outgoing fields
    ! (in case of nested coupling; cell_point_id(0) contains the combined grid)
    jg = MERGE(0, 1, n_dom > 1)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "surface_downward_eastward_stress", 2, out_field_ids%umfl)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "surface_downward_northward_stress", 2, out_field_ids%vmfl)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "surface_fresh_water_flux", 3, out_field_ids%freshflx)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "total_heat_flux", 4, out_field_ids%heatflx)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "atmosphere_sea_ice_bundle", 2, out_field_ids%seaice_atm)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "10m_wind_speed", 1, out_field_ids%sp10m)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "co2_mixing_ratio", 1, out_field_ids%co2_vmr)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "sea_level_pressure", 1, out_field_ids%pres_msl)

    ALLOCATE(in_field_ids(n_dom))

    ! define incoming fields
    ! (data is received on all domains)
    DO jg = 1, n_dom
      CALL cpl_def_field( &
        comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
        "sea_surface_temperature", 1, in_field_ids(jg)%sst)

      IF (use_ocean_velocity) THEN
        ! ocean and sea ice velocity bundle
        CALL cpl_def_field( &
          comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
          "surface_velocity_bundle", 4, in_field_ids(jg)%surface_velocity)

      ELSE
        in_field_ids(jg)%surface_velocity = -1
      END IF

      CALL cpl_def_field( &
        comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
        "ocean_sea_ice_bundle", 3, in_field_ids(jg)%seaice_oce)

      CALL cpl_def_field( &
        comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
        "co2_flux", 1, in_field_ids(jg)%co2_flx)
    END DO

  END SUBROUTINE construct_atmo_ocean_coupling_common

  !>
  !! This subroutine ensures consistency in the coupling definition and is
  !! called after the coupling definition phase
  SUBROUTINE construct_atmo_ocean_coupling_common_finalize()

    CHARACTER(len=*), PARAMETER :: &
    &  routine = str_module//':construct_atmo_ocean_coupling_common_finalize'

    CALL check_initial_datetime(out_field_ids%umfl)
    CALL check_initial_datetime(out_field_ids%vmfl)
    CALL check_initial_datetime(out_field_ids%freshflx)
    CALL check_initial_datetime(out_field_ids%heatflx)
    CALL check_initial_datetime(out_field_ids%seaice_atm)
    CALL check_initial_datetime(out_field_ids%sp10m)
    CALL check_initial_datetime(out_field_ids%pres_msl)

    !check for the global domain of the carbon cycle is interactive
    IF(ccycle_config(1)%iccycle == CCYCLE_MODE_INTERACTIVE) &
      CALL check_initial_datetime(out_field_ids%co2_vmr)

  CONTAINS

    SUBROUTINE check_initial_datetime(field_id)

      INTEGER, INTENT(IN) :: field_id

      TYPE(datetime), SAVE :: ref_field_datetime
      LOGICAL, SAVE :: ref_field_datetime_is_set = .FALSE.

      TYPE(datetime) ::  curr_field_datetime

      IF( cpl_get_field_is_source(routine, field_id) ) THEN
        curr_field_datetime = cpl_get_field_datetime(routine, field_id)
        IF (ref_field_datetime_is_set) THEN
          IF (ref_field_datetime /= curr_field_datetime) THEN
            CALL finish( &
              routine, &
              "inconsistent definition of field datetime in atm-oce-coupling &
              &for " // cpl_get_field_name(routine, field_id))
          END IF
        ELSE
          ref_field_datetime = curr_field_datetime
          ref_field_datetime_is_set = .TRUE.
        END IF
      END IF

    END SUBROUTINE check_initial_datetime

  END SUBROUTINE construct_atmo_ocean_coupling_common_finalize

  SUBROUTINE destruct_atmo_ocean_coupling_common()

    DEALLOCATE(in_field_ids)

  END SUBROUTINE destruct_atmo_ocean_coupling_common

END MODULE mo_atmo_ocean_coupling_common
