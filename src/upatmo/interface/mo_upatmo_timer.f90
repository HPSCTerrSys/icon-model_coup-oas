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

MODULE mo_upatmo_timer

#ifdef __SCT__
  USE sct, ONLY: new_timer => sct_new_timer
#else
  USE mo_real_timer, ONLY: new_timer
#endif

  USE mo_run_config, ONLY: ltimer

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_upatmo_timer        !< procedure of this module

  ! upper atmosphere
  PUBLIC :: timer_expol
  PUBLIC :: timer_upatmo, timer_upatmo_constr, timer_upatmo_destr, timer_upatmo_phy, &
      &       timer_upatmo_phy_init, timer_upatmo_phy_tend, timer_upatmo_phy_diag,     &
      &       timer_upatmo_phy_imf, timer_upatmo_phy_rad, timer_upatmo_phy_acc

  !-------------------
  ! Module variables
  !-------------------

  ! upper atmosphere
  INTEGER :: timer_expol
  INTEGER :: timer_upatmo, timer_upatmo_constr, timer_upatmo_destr, timer_upatmo_phy, &
      &        timer_upatmo_phy_init, timer_upatmo_phy_tend, timer_upatmo_phy_diag,     &
      &        timer_upatmo_phy_imf, timer_upatmo_phy_rad, timer_upatmo_phy_acc

CONTAINS

  SUBROUTINE init_upatmo_timer

    IF (.NOT.ltimer)  RETURN

    ! upper atmosphere
    timer_expol           = new_timer("upatmo_expol")
    timer_upatmo          = new_timer("upper_atmosphere")
    timer_upatmo_constr   = new_timer("upatmo_construction")
    timer_upatmo_destr    = new_timer("upatmo_destruction")
    timer_upatmo_phy      = new_timer("upatmo_physics")
    timer_upatmo_phy_init = new_timer("upatmo_phy_initialization")
    timer_upatmo_phy_tend = new_timer("upatmo_phy_update_tendencies")
    timer_upatmo_phy_diag = new_timer("upatmo_phy_update_diag_vars")
    timer_upatmo_phy_imf  = new_timer("upatmo_phy_group_imf")
    timer_upatmo_phy_rad  = new_timer("upatmo_phy_group_rad")
    timer_upatmo_phy_acc  = new_timer("upatmo_phy_accmlt_tendencies")

  END SUBROUTINE init_upatmo_timer

END MODULE mo_upatmo_timer
