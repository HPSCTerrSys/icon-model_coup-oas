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

MODULE mo_ragnarok_timer
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_ragnarok_timer

  INTERFACE
    SUBROUTINE init_ragnarok_timer() BIND(c)
    END SUBROUTINE init_ragnarok_timer
  END INTERFACE

END MODULE mo_ragnarok_timer
