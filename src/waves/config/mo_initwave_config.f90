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

! Initwave config

MODULE mo_initwave_config

  USE mo_impl_constants,       ONLY: max_dom
  USE mo_io_units,             ONLY: filename_max

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_initwave_config'

  ! types
  PUBLIC :: t_initwave_config

  ! objects
  PUBLIC :: initwave_config


  TYPE t_initwave_config
    INTEGER :: init_mode     !< MODE_ANA : read wave energy spectrum from analysis file
                             !< MODE_COLD: initialize by analytic wind-speed based parameterization
                             !             (such as JONSWAP)
    LOGICAL :: lskip_inv_post_op !< skip inverse post op (i.e. scaling of the input variable)
                                 !  when reading the state file
    CHARACTER(LEN=filename_max)::initial_wave_spectrum_filename

  END TYPE t_initwave_config

  TYPE(t_initwave_config), TARGET:: initwave_config(max_dom)

CONTAINS

END MODULE mo_initwave_config
