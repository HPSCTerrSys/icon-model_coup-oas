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


#if (defined(_OPENMP) && defined(OCE_SOLVE_OMP))
#include "omp_definitions.inc"
#endif
! contains extension to solver backend type: CG

MODULE mo_ocean_solve_cg

  USE mo_kind, ONLY: wp, sp
  USE mo_exception, ONLY: finish
  USE mo_ocean_solve_backend, ONLY: t_ocean_solve_backend
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_ocean_solve_cg
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_solve_cg'

  TYPE, EXTENDS(t_ocean_solve_backend) :: t_ocean_solve_cg
    PRIVATE
! arrays only used by CG
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: z_wp, d_wp, r_wp, rsq_wp
    REAL(KIND=sp), ALLOCATABLE, DIMENSION(:,:) :: z_sp, d_sp, r_sp, rsq_sp
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: z_wp, d_wp, r_wp, rsq_wp
!DIR$ ATTRIBUTES ALIGN : 64 :: z_sp, d_sp, r_sp, rsq_sp
#endif
! interfaces
  CONTAINS
    PROCEDURE :: doit_wp => ocean_solve_cg_cal_wp ! override deferred
    PROCEDURE :: doit_sp => ocean_solve_cg_cal_sp ! override deferred
    PROCEDURE, PRIVATE :: recover_arrays_wp => ocean_solve_cg_recover_arrays_wp
    PROCEDURE, PRIVATE :: recover_arrays_sp => ocean_solve_cg_recover_arrays_sp
    GENERIC, PRIVATE :: recover_arrays => recover_arrays_wp, recover_arrays_sp
  END TYPE t_ocean_solve_cg

CONTAINS

! get solver arrays (alloc them, if not done so, yet) - wp-variant
  SUBROUTINE ocean_solve_cg_recover_arrays_wp(this, x, b, z, d, r, r2, lacc)
    CLASS(t_ocean_solve_cg), INTENT(INOUT), TARGET :: this
    REAL(KIND=wp), INTENT(INOUT), POINTER, DIMENSION(:,:) :: &
      & x, b, z, d, r, r2
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: nblk_e
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    nblk_e = MERGE(this%trans%nblk, 1, this%trans%nblk > 0)

    IF (.NOT.ALLOCATED(this%z_wp)) THEN
      ALLOCATE(this%z_wp(this%trans%nidx, nblk_e), &
        & this%d_wp(this%trans%nidx, this%trans%nblk_a), &
        & this%r_wp(this%trans%nidx, nblk_e), &
        & this%rsq_wp(this%trans%nidx, nblk_e))
      this%d_wp(:, this%trans%nblk+1:this%trans%nblk_a) = 0._wp
      !$ACC ENTER DATA COPYIN(this%z_wp, this%d_wp, this%r_wp, this%rsq_wp) ASYNC(1) IF(lzacc)
      !$ACC WAIT(1)
    END IF
    x => this%x_wp
    b => this%b_wp
    z => this%z_wp
    d => this%d_wp
    r => this%r_wp
    r2 => this%rsq_wp
  END SUBROUTINE ocean_solve_cg_recover_arrays_wp

! actual CG solve (vanilla) - wp-variant
SUBROUTINE ocean_solve_cg_cal_wp(this, lacc)
    CLASS(t_ocean_solve_cg), INTENT(INOUT) :: this
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    REAL(KIND=wp) :: alpha, beta, dz_glob, tol, tol2, rn, rn_last
    INTEGER :: nidx, nidx_e, nblk, nblk_e, iblk, idx, k, m, k_final
#if defined(_OPENACC) || defined(__NO_CONT_SOLV_OCE__)
    REAL(KIND=wp), POINTER, DIMENSION(:,:) :: x, b, z, d, r, r2
#else
    REAL(KIND=wp), POINTER, DIMENSION(:,:), CONTIGUOUS :: x, b, z, d, r, r2
#endif
    LOGICAL :: done, lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

! retrieve extends of vector to solve
    nblk = this%trans%nblk
    nblk_e = MERGE(this%trans%nblk, 1, this%trans%nblk > 0)
    nidx = this%trans%nidx
    nidx_e = this%trans%nidx_e
    m = this%par%m
    k_final = -1

! retrieve arrays
    CALL this%recover_arrays(x, b, z, d, r, r2, lacc=lzacc)

    !$ACC DATA PRESENT(x, b, z, d, r, r2) IF(lzacc)

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO idx = nidx_e+1, nidx
      b(idx, nblk_e) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

! compute initial residual and auxiliary vectors
    CALL this%trans%sync(x, lacc=lzacc)

    CALL this%lhs%apply(x, z, lacc=lzacc)

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO idx = nidx_e+1, nidx
      z(idx, nblk_e) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    ! The nblk_e upper limit is such that the loop runs at least once, even if the domain is empty.
    ! This ensures that r, d, and r2 are initialized.
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(lzacc)
    DO iblk = 1, nblk_e
      DO idx = 1, nidx
        r(idx, iblk) = b(idx, iblk) - z(idx, iblk)
        d(idx, iblk) = r(idx, iblk)
        r2(idx, iblk) = r(idx, iblk) * r(idx, iblk)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
!ICON_OMP END PARALLEL DO

    CALL this%trans%global_sum(r2, rn, lacc=lzacc)

    ! tolerance
    tol = this%abs_tol_wp
    IF (.NOT.this%par%use_atol) THEN
      tol = this%par%tol * SQRT(rn)
      this%abs_tol_wp = tol
    END IF
    tol2 = tol * tol
    done = .false.

! enter CG iteration
    DO k = 1, m
! check if done
      IF (done) CYCLE
!      IF (this%trans%is_leader_pe) PRINT*,"it, res",k,SQRT(rn)
! check if already reached desired tolerance
      IF (rn .LE. tol2) THEN
        done = .true.
        k_final = k
        CYCLE
      END IF
! correct search direction (in direction of gradient) / update search vector
      IF (k .GT. 1) THEN
        beta = rn / rn_last

        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO idx = nidx_e+1, nidx
          d(idx, nblk_e) = 0._wp
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC WAIT(1)

!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(lzacc)
        DO iblk = 1, nblk
          DO idx = 1, nidx
            d(idx, iblk) = r(idx, iblk) + beta * d(idx, iblk)
          END DO
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC WAIT(1)
!ICON_OMP END PARALLEL DO
      END IF

      CALL this%trans%sync(d, lacc=lzacc)

      CALL this%lhs%apply(d, z, lacc=lzacc)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO idx = nidx_e+1, nidx
        d(idx, nblk_e) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
! compute extrapolated location of minimum in direction of d
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(lzacc)
      DO iblk = 1, nblk
        DO idx = 1, nidx
          r2(idx, iblk) = d(idx, iblk) * z(idx, iblk)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
!ICON_OMP END PARALLEL DO
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO idx = nidx_e+1, nidx
        r2(idx, nblk_e) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      CALL this%trans%global_sum(r2, dz_glob, lacc=lzacc)
      alpha = rn / dz_glob
! update guess and residuum
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(lzacc)
      DO iblk = 1, nblk
        DO idx = 1, nidx
          x(idx, iblk) = x(idx, iblk) + alpha * d(idx, iblk)
          r(idx, iblk) = r(idx, iblk) - alpha * z(idx, iblk)
          r2(idx, iblk) = r(idx, iblk) * r(idx, iblk)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
!ICON_OMP END PARALLEL DO
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO idx = nidx_e+1, nidx
        r2(idx, nblk_e) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
! save old and compute new residual norm
      rn_last = rn

      CALL this%trans%global_sum(r2, rn, lacc=lzacc)
    END DO
    this%niter_cal(1) = k_final
    this%res_wp(1) = SQRT(rn)

    CALL this%trans%sync(x, lacc=lzacc)

    !$ACC END DATA
  END SUBROUTINE ocean_solve_cg_cal_wp

! get solver arrays (alloc them, if not done so, yet) - sp-variant
  SUBROUTINE ocean_solve_cg_recover_arrays_sp(this, x, b, z, d, r, r2)
    CLASS(t_ocean_solve_cg), INTENT(INOUT), TARGET :: this
    REAL(KIND=sp), INTENT(INOUT), POINTER, DIMENSION(:,:) :: &
      & x, b, z, d, r, r2

    IF (.NOT.ALLOCATED(this%z_sp)) THEN
      ALLOCATE(this%z_sp(this%trans%nidx, this%trans%nblk), &
        & this%d_sp(this%trans%nidx, this%trans%nblk_a), &
        & this%r_sp(this%trans%nidx, this%trans%nblk), &
        & this%rsq_sp(this%trans%nidx, this%trans%nblk))
      this%d_sp(:, this%trans%nblk+1:this%trans%nblk_a) = 0._sp
    END IF
    x => this%x_sp
    b => this%b_sp
    z => this%z_sp
    d => this%d_sp
    r => this%r_sp
    r2 => this%rsq_sp
  END SUBROUTINE ocean_solve_cg_recover_arrays_sp

! we should not get here...
  SUBROUTINE ocean_solve_cg_cal_sp(this)
    CLASS(t_ocean_solve_cg), INTENT(INOUT) :: this
    REAL(KIND=sp) :: alpha, beta, dz_glob, tol, tol2, rn, rn_last
    INTEGER :: nidx_e, nblk, nblk_e, iblk, k, m, k_final
#if defined(_OPENACC) || defined(__NO_CONT_SOLV_OCE__)
    REAL(KIND=sp), POINTER, DIMENSION(:,:) :: x, b, z, d, r, r2
#else
    REAL(KIND=sp), POINTER, DIMENSION(:,:), CONTIGUOUS :: x, b, z, d, r, r2
#endif
    LOGICAL :: done

#ifdef _OPENACC
    CALL finish("mo_ocean_solve_cg::ocean_solve_cg_cal_sp", "not ported to GPU")
#endif

! retrieve extends of vector to solve
    nblk = this%trans%nblk
    nblk_e = MERGE(this%trans%nblk, 1, this%trans%nblk > 0)
    nidx_e = this%trans%nidx_e
    m = this%par%m
    k_final = -1
! retrieve arrays
    CALL this%recover_arrays(x, b, z, d, r, r2)
    b(nidx_e+1:, nblk_e) = 0._sp
! compute initial residual and auxiliary vectors
    CALL this%trans%sync(x)
    CALL this%lhs%apply(x, z)
    z(nidx_e+1:, nblk_e) = 0._sp
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO iblk = 1, nblk
      r(:, iblk) = b(:, iblk) - z(:, iblk)
      d(:, iblk) = r(:, iblk)
      r2(:, iblk) = r(:, iblk) * r(:, iblk)
    END DO
!ICON_OMP END PARALLEL DO
    CALL this%trans%global_sum(r2, rn)
    ! tolerance
    tol = REAL(this%abs_tol_sp, sp)
    IF (.NOT.this%par%use_atol) THEN
      tol = this%par%tol * SQRT(rn)
      this%abs_tol_sp = REAL(tol, wp)
    END IF
    tol2 = tol * tol
    done = .false.
! enter CG iteration
    DO k = 1, m
! check if done
      IF (done) CYCLE
!      IF (this%trans%is_leader_pe) PRINT*,"it, res",k,SQRT(rn)
! check if already reached desired tolerance
      IF (rn .LE. tol2) THEN
        done = .true.
        k_final = k
        CYCLE
      END IF
! correct search direction (in direction of gradient) / update search vector
      IF (k .GT. 1) THEN
        beta = rn / rn_last
        d(nidx_e+1:, nblk_e) = 0._sp
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
        DO iblk = 1, nblk
          d(:, iblk) = r(:, iblk) + beta * d(:, iblk)
        END DO
!ICON_OMP END PARALLEL DO
      END IF
      CALL this%trans%sync(d)
      CALL this%lhs%apply(d, z)
      d(nidx_e+1:, nblk_e) = 0._sp
! compute extrapolated location of minimum in direction of d
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        r2(:, iblk) = d(:, iblk) * z(:, iblk)
      END DO
!ICON_OMP END PARALLEL DO
      r2(nidx_e+1:, nblk_e) = 0._sp
      CALL this%trans%global_sum(r2, dz_glob)
      alpha = rn / dz_glob
! update guess and residuum
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        x(:, iblk) = x(:, iblk) + alpha * d(:, iblk)
        r(:, iblk) = r(:, iblk) - alpha * z(:, iblk)
        r2(:, iblk) = r(:, iblk) * r(:, iblk)
      END DO
!ICON_OMP END PARALLEL DO
      r2(nidx_e+1:, nblk_e) = 0._sp
! save old and compute new residual norm
      rn_last = rn
      CALL this%trans%global_sum(r2, rn)
    END DO
    this%niter_cal(2) = k_final
    this%res_wp(2) = REAL(SQRT(rn), wp)
    CALL this%trans%sync(x)
  END SUBROUTINE ocean_solve_cg_cal_sp

END MODULE mo_ocean_solve_cg
