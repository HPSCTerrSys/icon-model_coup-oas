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
MODULE mo_ragnarok_microphysics
#ifdef __SINGLE_PRECISION
  USE :: ISO_C_BINDING, ONLY: c_int, wp => c_float
#else
  USE :: ISO_C_BINDING, ONLY: c_int, wp => c_double
#endif

  IMPLICIT NONE
  PRIVATE

  INTERFACE

    SUBROUTINE run(nvec, ke, ivstart, ivend, kstart, dt, cia, dz, t, rho, p, qv,   &
        &                         qc, qi, qr, qs, qg, qnc, prr_gsp, pri_gsp, prs_gsp, prg_gsp, &
        &                         pre_gsp, pflx) BIND(c)

      IMPORT c_int, wp
      ! arguments
      INTEGER(c_int), VALUE, INTENT(IN) :: nvec, ke, ivstart, ivend, kstart
      REAL(wp), VALUE, INTENT(IN) :: dt
      REAL(wp), VALUE, INTENT(IN) :: cia
      REAL(wp), DIMENSION(*), INTENT(IN) :: dz, rho, p
      REAL(wp), DIMENSION(*), INTENT(INOUT) :: t
      REAL(wp), DIMENSION(*), INTENT(INOUT) :: qv, qc, qi, qr, qs, qg
      REAL(wp), DIMENSION(*), INTENT(IN) :: qnc
      REAL(wp), DIMENSION(*), INTENT(OUT) :: pflx
      REAL(wp), DIMENSION(*), INTENT(INOUT) :: prr_gsp, pri_gsp, prs_gsp, prg_gsp, pre_gsp

    END SUBROUTINE run

  END INTERFACE

  PUBLIC :: graupel_run

CONTAINS

  SUBROUTINE graupel_run(nvec, ke, ivstart, ivend, kstart, dt, cia, qnc, dz, rho, p, t, qv,   &
      &                         qc, qi, qr, qs, qg, pflx, prr_gsp, pri_gsp, prs_gsp, prg_gsp, &
      &                         pre_gsp)

    INTEGER, INTENT(IN) :: nvec, ke, ivstart, ivend, kstart
    REAL(wp), INTENT(IN) :: dt
    REAL(wp), INTENT(IN) :: cia
    REAL(wp), DIMENSION(:, :), INTENT(IN) :: dz, rho, p
    REAL(wp), DIMENSION(:, :), INTENT(INOUT) :: t
    REAL(wp), DIMENSION(:, :), INTENT(INOUT) :: qv, qc, qi, qr, qs, qg
    REAL(wp), DIMENSION(:), INTENT(IN) :: qnc
    REAL(wp), DIMENSION(:, :), INTENT(OUT) :: pflx
    REAL(wp), DIMENSION(:), INTENT(INOUT) :: prr_gsp, pri_gsp, prs_gsp, prg_gsp, pre_gsp

    CALL run(nvec, ke, ivstart - 1, ivend, kstart - 1, dt, cia, dz, t, rho, p, qv,  &
        &                         qc, qi, qr, qs, qg, qnc, prr_gsp, pri_gsp, prs_gsp, prg_gsp, &
        &                         pre_gsp, pflx)

  END SUBROUTINE graupel_run

END MODULE mo_ragnarok_microphysics
