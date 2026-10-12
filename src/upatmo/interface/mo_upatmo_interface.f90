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

! Expose public types and procedures for upatmo physics.
!
! This module contains no implementations. The public types and procedures of
! other modules in src/upatmo/interface are used directly by ICON rather than
! through this interface due to dependency cycles.

MODULE mo_upatmo_interface

  USE mo_nwp_upatmo_interface, ONLY: nwp_upatmo_interface, nwp_upatmo_update

  USE mo_upatmo_flowevent_utils, ONLY: t_upatmoRestartAttributes, &
      upatmoRestartAttributesAssign, &
      upatmoRestartAttributesPack, &
      upatmoRestartAttributesSet, &
      upatmoRestartAttributesGet, &
      upatmoRestartAttributesPrepare, &
      upatmoRestartAttributesDeallocate

  USE mo_upatmo_nml, ONLY: read_upatmo_namelist

  USE mo_upatmo_phy_setup, ONLY: init_upatmo_phy_nwp, finalize_upatmo_phy_nwp

  USE mo_upatmo_state, ONLY: construct_upatmo_state, destruct_upatmo_state, &
      prm_upatmo

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
      nwp_upatmo_interface, nwp_upatmo_update, &
      t_upatmoRestartAttributes, &
      upatmoRestartAttributesAssign, &
      upatmoRestartAttributesPack, &
      upatmoRestartAttributesSet, &
      upatmoRestartAttributesGet, &
      upatmoRestartAttributesPrepare, &
      upatmoRestartAttributesDeallocate, &
      read_upatmo_namelist, &
      init_upatmo_phy_nwp, finalize_upatmo_phy_nwp, &
      construct_upatmo_state, destruct_upatmo_state, &
      prm_upatmo

END MODULE mo_upatmo_interface
