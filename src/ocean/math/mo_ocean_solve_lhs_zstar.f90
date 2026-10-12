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

! contains extended lhs-matrix-generator type
! provides surface height lhs for free ocean surface solve

#if (defined(_OPENMP) && defined(OCE_SOLVE_OMP))
#include "omp_definitions.inc"
#endif

MODULE mo_lhs_zstar

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: finish
  USE mo_ocean_solve_lhs_type, ONLY: t_lhs_agen
  USE mo_ocean_types, ONLY: t_solverCoeff_singlePrecision, t_operator_coeff
  USE mo_grid_subset, ONLY: t_subset_range, get_index_range
  USE mo_ocean_nml, ONLY: select_lhs, select_lhs_matrix, ab_gam, ab_beta, &
    & l_edge_based, l_lhs_direct, ocean_latbc_bnd_intp_width, &
    & is_ocean_limited_area, select_lhs_operators
  USE mo_communication, ONLY: exchange_data
  USE mo_parallel_config, ONLY: nproma
  USE mo_impl_constants, ONLY: sea_boundary
  USE mo_run_config, ONLY: dtime, debug_check_level
  USE mo_physical_constants, ONLY: grav
  USE mo_ocean_math_operators, ONLY: &
    & grad_fd_norm_oce_2d_3d, div_oce_2D_onTriangles_onBlock, &
    & div_oce_2D_general_onBlock
  USE mo_scalar_product, ONLY: map_edges2edges_viacell_2D_zstar
  USE mo_model_domain, ONLY: t_patch_3d, t_patch
  USE mo_mpi,  ONLY: my_process_is_mpi_parallel
  USE mo_sync, ONLY: sync_e, sync_patch_array
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

#ifdef _OPENACC
  USE openacc, ONLY: acc_is_present
#endif

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_surface_height_lhs_zstar'

  PUBLIC :: t_lhs_zstar

  TYPE, EXTENDS(t_lhs_agen) :: t_lhs_zstar
    PRIVATE
    TYPE(t_patch_3d), POINTER :: patch_3d => NULL()
    TYPE(t_patch), POINTER :: patch_2d => NULL()
    REAL(wp), POINTER :: thickness_e_wp(:,:) => NULL()
    REAL(wp), POINTER, PUBLIC :: stretch_e(:,:) => NULL()
    TYPE(t_operator_coeff), POINTER :: op_coeffs_wp => NULL()
    TYPE(t_solverCoeff_singlePrecision), POINTER :: op_coeffs_sp => NULL()
    REAL(wp), ALLOCATABLE, DIMENSION(:,:), PRIVATE :: z_grad_h_wp, z_e_wp
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: z_grad_h_wp, z_e_wp
#endif
  CONTAINS
    PROCEDURE :: lhs_wp => lhs_zstar
    PROCEDURE :: construct => lhs_zstar_construct
    PROCEDURE, PRIVATE :: internal_matrix_wp => lhs_zstar_matrix_wp
    PROCEDURE, PRIVATE :: internal_wp => lhs_zstar_wp
    PROCEDURE :: lhs_matrix_shortcut => lhs_zstar_matrix_shortcut
  END TYPE t_lhs_zstar

  INTEGER, PARAMETER :: topLevel = 1

  CONTAINS

  !init generator object
  SUBROUTINE lhs_zstar_construct(this, patch_3d, thick_e, &
      & op_coeffs_wp, op_coeffs_sp, lacc)
    CLASS(t_lhs_zstar), INTENT(INOUT) :: this
    TYPE(t_patch_3d), POINTER, INTENT(IN) :: patch_3d
    REAL(wp), POINTER, INTENT(IN) :: thick_e(:,:)
    TYPE(t_operator_coeff), TARGET, INTENT(IN) :: op_coeffs_wp
    TYPE(t_solverCoeff_singlePrecision), TARGET, INTENT(IN) :: op_coeffs_sp
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    this%patch_3d => patch_3d
    this%patch_2d => patch_3d%p_patch_2d(1)
    this%thickness_e_wp => thick_e
    this%op_coeffs_wp => op_coeffs_wp
    this%op_coeffs_sp => op_coeffs_sp
    this%is_const = .false.
    this%use_shortcut = .false. !(select_lhs .GT. select_lhs_matrix .AND. select_lhs .LE. select_lhs_matrix + 1)
    IF (this%patch_2d%cells%max_connectivity .NE. 3 .AND. .NOT.l_lhs_direct) &
      & CALL finish("t_surface_height_lhs::lhs_surface_height_construct", &
      &  "internal matrix implementation only works with triangular grids!")

    !$ACC ENTER DATA COPYIN(this, this%patch_3d, this%patch_2d, this%thickness_e_wp) &
    !$ACC   COPYIN(this%op_coeffs_wp, this%op_coeffs_sp) ASYNC(1) IF(lzacc)
    !$ACC WAIT(1)

    this%is_init = .true.
  END SUBROUTINE lhs_zstar_construct

! interface routine for the left hand side computation
  SUBROUTINE lhs_zstar(this, x, ax, lacc)
    CLASS(t_lhs_zstar), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN) :: x(:,:)
    REAL(wp), INTENT(INOUT) :: ax(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    LOGICAL :: lzacc
    TYPE(t_subset_range), POINTER :: edges_in_domain

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA PRESENT(this%patch_3d%p_patch_2d(1)%edges%cell_idx, this%patch_3d%p_patch_2d(1)%edges%cell_blk) &
    !$ACC   PRESENT(this%patch_3d%p_patch_2d(1)%cells%edge_idx, this%patch_3d%p_patch_2d(1)%cells%edge_blk) &
    !$ACC   PRESENT(this%patch_3d%p_patch_1d(1)%dolic_e, this%patch_2d%cells%in_domain) &
    !$ACC   PRESENT(this%op_coeffs_wp%grad_coeff, this%op_coeffs_wp%div_coeff, x) &
    !$ACC   PRESENT(ax) IF(lzacc)

    IF (this%use_shortcut) &
      & CALL finish("t_surface_height_lhs::lhs_surface_height_wp", &
        & "should not be here because of shortcut!")
    IF ( (select_lhs .NE. select_lhs_operators) ) & !.OR. (l_lhs_direct .EQV. .false.) ) &
      & CALL finish("zstar currently only supports operator based solver", &
        & "Set select_lhs=1!!") ! and l_lhs_direct=.true.!!")
! alloc internal arrays, if needed to
! if using operator implementation of lhs
    IF (select_lhs .EQ. select_lhs_operators) THEN
      IF (ALLOCATED(this%z_e_wp)) THEN
        IF (nproma .NE. SIZE(this%z_e_wp, 1) .OR. &
            & this%patch_2d%nblks_e .NE. SIZE(this%z_e_wp, 2)) THEN
          !$ACC EXIT DATA DELETE(this%z_e_wp, this%z_grad_h_wp) ASYNC(1) IF(lzacc)
          !$ACC WAIT(1)
          DEALLOCATE(this%z_e_wp, this%z_grad_h_wp)
        END IF
      END IF
      IF (.NOT.ALLOCATED(this%z_e_wp)) THEN
        ALLOCATE(this%z_e_wp(nproma, this%patch_2d%nblks_e), &
            & this%z_grad_h_wp(nproma, this%patch_2d%nblks_e))
#ifdef NAGFOR
          this%z_e_wp(:,:) = 0._wp
          this%z_e_wp(:,:)   = 0._wp
          this%z_grad_h_wp(:,:) = 0._wp
          this%z_grad_h_wp(:,:)   = 0._wp
#endif
          !$ACC ENTER DATA COPYIN(this%z_e_wp, this%z_grad_h_wp) ASYNC(1) IF(lzacc)
        !$ACC WAIT(1)
      END IF
      CALL this%internal_wp(x, ax, lacc=lzacc)
    ELSE if (select_lhs .EQ. select_lhs_matrix) THEN
      CALL this%internal_matrix_wp(x, ax, lacc=lzacc)
    END IF
    !$ACC WAIT(1)
    !$ACC END DATA
  END SUBROUTINE lhs_zstar

  ! internal backend routine to compute surface height lhs -- "matrix" implementation
  SUBROUTINE lhs_zstar_matrix_wp(this, x, lhs, lacc)
    CLASS(t_lhs_zstar), INTENT(INOUT) :: this
#if defined(_OPENACC) || defined(__NO_CONT_SOLV_OCE__)
    REAL(wp), INTENT(IN) :: x(:,:)
    REAL(wp), INTENT(INOUT) :: lhs(:,:)
#else
    REAL(wp), INTENT(IN), CONTIGUOUS :: x(:,:)
    REAL(wp), INTENT(INOUT), CONTIGUOUS :: lhs(:,:)
#endif
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    INTEGER :: start_index, end_index, jc, blkNo
    TYPE(t_subset_range), POINTER :: cells_in_domain
#if defined(_OPENACC) || defined(__NO_CONT_SOLV_OCE__)
    REAL(wp), POINTER, DIMENSION(:,:,:) :: lhs_coeffs
    INTEGER, POINTER, DIMENSION(:,:,:) :: idx, blk
#else
    REAL(wp), POINTER, DIMENSION(:,:,:), CONTIGUOUS :: lhs_coeffs
    INTEGER, POINTER, DIMENSION(:,:,:), CONTIGUOUS :: idx, blk
#endif
    LOGICAL :: lzacc

CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
        IF (lzacc) CALL finish("mo_surface_height_lhs_zstar::lhs_surface_height_ab_mim_matrix_wp", &
          & "OpenACC version currently not implemented")
#endif

    cells_in_domain => this%patch_2D%cells%in_domain
    lhs_coeffs => this%op_coeffs_wp%lhs_all
    idx => this%op_coeffs_wp%lhs_CellToCell_index
    blk => this%op_coeffs_wp%lhs_CellToCell_block
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO blkNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blkNo, start_index, end_index)
      lhs(:,blkNo) = 0.0_wp
      DO jc = start_index, end_index
        IF(.NOT.(this%patch_3d%lsm_c(jc,1,blkNo) > sea_boundary)) THEN
          lhs(jc,blkNo) = x(jc,blkNo) * lhs_coeffs(0, jc, blkNo) + &
               & x(idx(1, jc, blkNo), blk(1, jc, blkNo)) * lhs_coeffs(1, jc, blkNo) + &
               & x(idx(2, jc, blkNo), blk(2, jc, blkNo)) * lhs_coeffs(2, jc, blkNo) + &
               & x(idx(3, jc, blkNo), blk(3, jc, blkNo)) * lhs_coeffs(3, jc, blkNo) + &
               & x(idx(4, jc, blkNo), blk(4, jc, blkNo)) * lhs_coeffs(4, jc, blkNo) + &
               & x(idx(5, jc, blkNo), blk(5, jc, blkNo)) * lhs_coeffs(5, jc, blkNo) + &
               & x(idx(6, jc, blkNo), blk(6, jc, blkNo)) * lhs_coeffs(6, jc, blkNo) + &
               & x(idx(7, jc, blkNo), blk(7, jc, blkNo)) * lhs_coeffs(7, jc, blkNo) + &
               & x(idx(8, jc, blkNo), blk(8, jc, blkNo)) * lhs_coeffs(8, jc, blkNo) + &
               & x(idx(9, jc, blkNo), blk(9, jc, blkNo)) * lhs_coeffs(9, jc, blkNo)
        END IF
      END DO
    END DO ! blkNo
!ICON_OMP_END_PARALLEL_DO
    IF (debug_check_level > 20) THEN
      DO blkNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blkNo, start_index, end_index)
        DO jc = start_index, end_index
          IF(this%patch_3d%lsm_c(jc,1,blkNo) > sea_boundary) THEN
            IF (lhs(jc,blkNo) /= 0.0_wp) &
              & CALL finish("lhs_surface_height_ab_mim", "lhs(jc,blkNo) /= 0 on land")
          ENDIF
        END DO
      END DO
    ENDIF
  END SUBROUTINE lhs_zstar_matrix_wp

! internal backend routine to compute surface height lhs -- "operator" implementation
  SUBROUTINE lhs_zstar_wp(this, x, lhs, lacc)
    CLASS(t_lhs_zstar), INTENT(INOUT) :: this
#if defined(_OPENACC) || defined(__NO_CONT_SOLV_OCE__)
    REAL(wp), INTENT(IN) :: x(:,:)
    REAL(wp), INTENT(INOUT) :: lhs(:,:)
#else
    REAL(wp), INTENT(IN), CONTIGUOUS :: x(:,:)
    REAL(wp), INTENT(INOUT), CONTIGUOUS :: lhs(:,:)
#endif
    LOGICAL, INTENT(in), OPTIONAL :: lacc

    REAL(wp) :: gdt2_inv, gam_times_beta
    INTEGER :: start_index, end_index, jc, blkNo
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    CHARACTER(len=*), PARAMETER :: routine = modname//':lhs_zstar_wp'
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    if (.NOT.ASSOCIATED(this%stretch_e)) CALL finish(routine,"stretch_e not pointing somewhere...")
    cells_in_domain => this%patch_2D%cells%in_domain
    edges_in_domain => this%patch_2D%edges%in_domain

    gdt2_inv = REAL(1.0_wp / (grav*(dtime)**2),wp)
    gam_times_beta = REAL(ab_gam * ab_beta, wp)

#ifdef NAGFOR
    lhs(:,:) = 0._wp
#endif

    !$ACC DATA PRESENT(cells_in_domain%end_block) &
    !$ACC   COPY(x, lhs) IF(lzacc)

    CALL grad_fd_norm_oce_2d_3d( x, this%patch_2D, this%op_coeffs_wp%grad_coeff(:,1,:), &
      & this%z_grad_h_wp(:,:), subset_range=this%patch_2D%edges%gradIsCalculable, lacc=lzacc)

    IF (this%patch_2d%cells%max_connectivity /= 3 .AND. my_process_is_mpi_parallel()) &
      & CALL exchange_data(p_pat=this%patch_2D%comm_pat_e, lacc=lzacc, recv=this%z_grad_h_wp)

    !! Multiply with mass matrix and and sum up over column
    !! The coefficients are arrived at by having a constant for all levels >= 2
    !! and updating the first level in update_thickness_dependent_operator_coeff
    CALL map_edges2edges_viacell_2D_zstar( this%patch_3d, &
      & this%z_grad_h_wp(:,:), this%op_coeffs_wp, this%stretch_e(:, :), this%z_e_wp(:,:), &
      & lacc=lzacc)

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO blkNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blkNo, start_index, end_index)
      IF (this%patch_2d%cells%max_connectivity .EQ. 3) THEN
        CALL div_oce_2D_onTriangles_onBlock(this%z_e_wp, this%patch_2D, &
          & this%op_coeffs_wp%div_coeff, lhs(:,blkNo), level=topLevel, &
          & blockNo=blkNo, start_index=start_index, end_index=end_index, lacc=lzacc)
      ELSE
        CALL div_oce_2D_general_onBlock(this%z_e_wp, this%patch_2D, &
          & this%op_coeffs_wp%div_coeff, lhs(:,blkNo), level=topLevel, &
          & blockNo=blkNo, start_index=start_index, end_index=end_index, lacc=lzacc)
      END IF
      !Step 4) Finalize LHS calculations
      IF (.NOT. is_ocean_limited_area) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = start_index, end_index
          lhs(jc,blkNo) = x(jc,blkNo) * gdt2_inv - gam_times_beta * lhs(jc,blkNo)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE ! LAM modification - cells in the boundary interpolation zone have custom lhs values
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = start_index, end_index
          IF (this%patch_2d%cells%refin_ctrl(jc,blkNo) == 0 .OR. &
            & this%patch_2d%cells%refin_ctrl(jc,blkNo) > ocean_latbc_bnd_intp_width) THEN
            lhs(jc,blkNo) = x(jc,blkNo) * gdt2_inv - gam_times_beta * lhs(jc,blkNo)
          ELSE
            lhs(jc,blkNo) = x(jc,blkNo) * gdt2_inv
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      ENDIF
      !$ACC WAIT(1)
    END DO ! blkNo
!ICON_OMP_END_PARALLEL_DO
    IF (debug_check_level > 20) THEN
      !$ACC UPDATE SELF(this%patch_3d%lsm_c, lhs) IF(lzacc)
      DO blkNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blkNo, start_index, end_index)
        DO jc = start_index, end_index
          IF(this%patch_3d%lsm_c(jc,1,blkNo) > sea_boundary) THEN
            IF (lhs(jc,blkNo) /= 0.0_wp) &
              & CALL finish(routine, "lhs(jc,blkNo) /= 0 on land")
          END IF
        END DO
      END DO
    ENDIF
    !$ACC END DATA
  END SUBROUTINE lhs_zstar_wp

  SUBROUTINE lhs_zstar_matrix_shortcut(this, idx, blk, coeff, lacc)
    CLASS(t_lhs_zstar), INTENT(INOUT) :: this
    INTEGER, INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: idx, blk
    REAL(KIND=wp), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: coeff
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    INTEGER :: nidx, nblk, nnz, iidx, iblk, inz
#if defined(_OPENACC) || defined(__NO_CONT_SOLV_OCE__)
    INTEGER, POINTER, DIMENSION(:,:,:) :: opc_idx, opc_blk
    REAL(wp), POINTER, DIMENSION(:,:,:) :: opc_coeff
#else
    INTEGER, POINTER, DIMENSION(:,:,:), CONTIGUOUS :: opc_idx, opc_blk
    REAL(wp), POINTER, DIMENSION(:,:,:), CONTIGUOUS :: opc_coeff
#endif
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish("t_surface_height_lhs_zstar::lhs_surface_height_ab_mim_matrix_shortcut", &
      & "OpenACC version currently not implemented")
#endif

    IF (.NOT.this%use_shortcut) &
      & CALL finish( &
        & "t_surface_height_lhs::lhs_surface_height_ab_mim_matrix_shortcut", &
        & "wrong turn!")
    nidx = SIZE(coeff, 1)
    nblk = SIZE(coeff, 2)
    nnz = 10
    opc_coeff => this%op_coeffs_wp%lhs_all
    IF (.NOT.ALLOCATED(idx)) THEN
      DEALLOCATE(coeff)
      ALLOCATE(idx(nidx, nblk, nnz), blk(nidx, nblk, nnz), &
        & coeff(nidx, nblk, nnz))
      opc_idx => this%op_coeffs_wp%lhs_CellToCell_index
      opc_blk => this%op_coeffs_wp%lhs_CellToCell_block
!ICON_OMP PARALLEL PRIVATE(iblk, iidx, inz)
!ICON_OMP DO SCHEDULE(GUIDED)
      DO iblk = 1, nblk
        idx(:, iblk, 1) = (/(iidx, iidx = 1, nidx)/)
        blk(:, iblk, 1) = iblk
        coeff(:, iblk, 1) = opc_coeff(0, :, iblk)
      END DO
!ICON_OMP END DO
!ICON_OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
      DO inz = 2, nnz
        DO iblk = 1, nblk
          DO iidx = 1, nidx
            IF (ABS(opc_coeff(inz - 1, iidx, iblk)) .EQ. 0._wp) THEN
              idx(iidx, iblk, inz) = idx(iidx, iblk, 1)
              blk(iidx, iblk, inz) = blk(iidx, iblk, 1)
              coeff(iidx, iblk, inz) = 0._wp
            ELSE
              idx(iidx, iblk, inz) = opc_idx(inz - 1, iidx, iblk)
              blk(iidx, iblk, inz) = opc_blk(inz - 1, iidx, iblk)
              coeff(iidx, iblk, inz) = opc_coeff(inz - 1, iidx, iblk)
            END IF
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
!ICON_OMP END PARALLEL
    ELSE
!ICON_OMP PARALLEL DO SCHEDULE(GUIDED) COLLAPSE(2)
      DO inz = 1, nnz
        DO iblk = 1, nblk
          coeff(:, iblk, inz) = opc_coeff(inz - 1, :, iblk)
        END DO
      END DO
!ICON_OMP END PARALLEL DO
    END IF
  END SUBROUTINE lhs_zstar_matrix_shortcut

END MODULE mo_lhs_zstar
