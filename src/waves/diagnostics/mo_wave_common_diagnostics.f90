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

! Shared diagnostics, which may be needed by any of the three states
! diag_dyn, diag_out, diag_cpl.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_common_diagnostics

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_fortran_tools,       ONLY: init

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_common_diagnostics'

  PUBLIC :: significant_wave_height

CONTAINS

  !>
  !! Calculation of total significant wave height
  !! based on WAM 4.5 formulation
  !!
  SUBROUTINE significant_wave_height(p_patch, emean, hs)

    TYPE(t_patch),     INTENT(IN)    :: p_patch
    REAL(wp),          INTENT(IN)    :: emean(:,:)  !< total energy [m^2]
    REAL(wp),          INTENT(INOUT) :: hs(:,:)     !< significant wave height [m]

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':significant_wave_height'

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb


    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx
        hs(jc,jb) = 4.0_wp * SQRT(emean(jc,jb))
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE significant_wave_height

END MODULE mo_wave_common_diagnostics
