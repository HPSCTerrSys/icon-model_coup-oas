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

! Classes and functions for the turbulent mixing package (tmx)

MODULE mo_tmx_field_class

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: finish
  USE mo_tmx_var, ONLY: t_tmx_var
  USE memman, ONLY: var_descriptor, add_var_dp, get_var_data
  ! Todo: refactor so that t_patch is not needed
  USE mo_model_domain      ,ONLY: t_patch

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_tmx_field, t_tmx_field_p, &! bind_tmx_field, &
    &       t_domain, isfc_oce, isfc_ice, isfc_lnd

  TYPE, EXTENDS(t_tmx_var) :: t_tmx_field
    INTEGER :: id = -1
    INTEGER :: diffusion_type = -1 ! Field type
    INTEGER, ALLOCATABLE :: ref_idx(:)
  END TYPE t_tmx_field

  INTERFACE t_tmx_field
    MODULE PROCEDURE t_tmx_field_constructor
  END INTERFACE

  TYPE :: t_tmx_field_p
    TYPE(t_tmx_field), POINTER :: p
  END TYPE t_tmx_field_p

  TYPE t_domain
    INTEGER ::              &
      & id,                 &
      & nlev = 0,           &
      & ntiles = 0,         &
      & nproma = 1,         &
      & npromz = 1,         &
      & nblks_c = 1,        & ! Number of blocks (cells)
      & i_startblk_c = 1,   & ! Start block on cells
      & i_endblk_c = 1,     & ! End block on cells
      & nblks_e = 1,        & ! Number of blocks (edges)
      & i_startblk_e = 1,   & ! Start block on edges
      & i_endblk_e = 1,     & ! End block on edges
      & nblks_v = 1           ! Number of blocks (vertice)
    INTEGER, ALLOCATABLE :: &
      & sfc_types(:),       & ! Surface type for each tile if ntiles > 0
      & i_startidx_c(:),    & ! Start indices on cells (for each block)
      & i_endidx_c(:),      & ! End indices on cells (for each block)
      & i_startidx_e(:),    & ! Start indices on edges (for each block)
      & i_endidx_e(:)         ! End indices on edges (for each block)
    REAL(wp), ALLOCATABLE :: &
      & lon(:,:), lat(:,:), area(:,:)
    TYPE(t_patch), POINTER :: patch
  END TYPE

  INTERFACE t_domain
    PROCEDURE t_domain_constructor
  END INTERFACE

  ! Surface types
  ! Todo: currently, 1-3 need to be consistent with surface types in aes!
  ENUM, BIND(C)
    ENUMERATOR :: isfc_oce=1, isfc_ice, isfc_lnd
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_field_class'

CONTAINS

  FUNCTION t_tmx_field_constructor(name, type_id, dims, diffusion_type, var_desc, ref_pos, ref_idx) RESULT(field)

    CHARACTER(len=*), INTENT(in) :: name
    CHARACTER(len=*), INTENT(in) :: type_id
    INTEGER,          INTENT(in) :: dims(:)
    INTEGER,          INTENT(in) :: diffusion_type
    TYPE(var_descriptor), INTENT(in) :: var_desc
    INTEGER, OPTIONAL,    INTENT(in) :: ref_pos
    INTEGER, OPTIONAL,    INTENT(in) :: ref_idx(:)
    TYPE(t_tmx_field), POINTER   :: field

    TYPE(t_tmx_var), POINTER :: var

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_tmx_field_constructor'

    var => t_tmx_var(name, type_id, var_desc, dims=dims, ref_pos=ref_pos)
    ALLOCATE(field)
    field%t_tmx_var = var
    CALL MOVE_ALLOC(var%dims, field%dims)
    DEALLOCATE(var)

    field%diffusion_type = diffusion_type

    IF (PRESENT(ref_pos)) THEN
      IF (ref_pos > field%rank) CALL finish(routine, name//' - ref_pos > rank')
      field%ref_pos = ref_pos
    END IF
    IF (PRESENT(ref_idx)) THEN
      IF (.NOT. PRESENT(ref_pos)) field%ref_pos = field%rank
      IF (ANY(ref_idx > field%dims(field%ref_pos))) CALL finish(routine, name//' - ref_idx > dims(ref_pos)')
      ALLOCATE(field%ref_idx(SIZE(ref_idx)))
      field%ref_idx(:) = ref_idx(:)
    END IF
    IF (field%ref_pos > 0 .AND. field%ref_pos < field%rank) CALL finish(routine, name//' - ref_pos < rank currently not supported')

  END FUNCTION t_tmx_field_constructor

  FUNCTION t_domain_constructor(patch, nproma, nlev, ntiles, sfc_types) RESULT(domain)

    USE mo_loopindices,        ONLY: get_indices_e, get_indices_c
    USE mo_impl_constants,     ONLY: min_rlcell_int, min_rledge_int
    USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
    USE mo_math_constants,     ONLY: rad2deg

    TYPE(t_domain), POINTER :: domain

    TYPE(t_patch), POINTER :: patch
    INTEGER, INTENT(in) :: &
      & nproma
    INTEGER, INTENT(in), OPTIONAL :: &
      & nlev, &
      & ntiles, &
      & sfc_types(:)

    INTEGER :: rl_start, rl_end, jb

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_domain_constructor'

    ALLOCATE(domain)

    domain%id = patch%id
    domain%nproma = nproma

    domain%patch => patch
    __acc_attach(domain%patch)

    domain%npromz = patch%npromz_c
    IF (PRESENT(nlev)) THEN
      domain%nlev = nlev
    END IF

    ALLOCATE(domain%lon(nproma,patch%nblks_c))
    ALLOCATE(domain%lat(nproma,patch%nblks_c))
    domain%lon(:,:) = rad2deg * patch%cells%center(:,:)%lon
    domain%lat(:,:) = rad2deg * patch%cells%center(:,:)%lat
    ALLOCATE(domain%area(nproma,patch%nblks_c))
    domain%area(:,:) = patch%cells%area(:,:)

    IF (PRESENT(ntiles)) THEN
      domain%ntiles = ntiles
    END IF
    IF (domain%ntiles > 0) THEN
      IF (.NOT. PRESENT(sfc_types)) CALL finish(routine, 'sfc_types required')
      ALLOCATE(domain%sfc_types, source=sfc_types)
    END IF

    domain%nblks_c = patch%nblks_c
    domain%nblks_e = patch%nblks_e
    domain%nblks_v = patch%nblks_v

    ALLOCATE(domain%i_startidx_c(patch%nblks_c))
    ALLOCATE(domain%i_endidx_c  (patch%nblks_c))
    ALLOCATE(domain%i_startidx_e(patch%nblks_e))
    ALLOCATE(domain%i_endidx_e  (patch%nblks_e))

    rl_start = grf_bdywidth_c + 1
    rl_end   = min_rlcell_int
    domain%i_startblk_c = patch%cells%start_block(rl_start)
    domain%i_endblk_c   = patch%cells%end_block(rl_end)
    DO jb=domain%i_startblk_c,domain%i_endblk_c
      CALL get_indices_c(patch, jb, domain%i_startblk_c, domain%i_endblk_c, &
        &                           domain%i_startidx_c(jb), domain%i_endidx_c(jb), rl_start, rl_end)
    END DO

    rl_start = grf_bdywidth_e + 1
    rl_end   = min_rledge_int
    domain%i_startblk_e = patch%edges%start_block(rl_start)
    domain%i_endblk_e   = patch%edges%end_block(rl_end)
    DO jb=domain%i_startblk_e,domain%i_endblk_e
      CALL get_indices_e(patch, jb, domain%i_startblk_e, domain%i_endblk_e, &
        &                           domain%i_startidx_e(jb), domain%i_endidx_e(jb), rl_start, rl_end)
    END DO

    !$ACC ENTER DATA COPYIN(domain)
    !$ACC ENTER DATA COPYIN(domain%lon, domain%lat, domain%area)
    !$ACC ENTER DATA COPYIN(domain%i_startidx_c, domain%i_endidx_c, domain%i_startidx_e, domain%i_endidx_e)
    !$ACC ENTER DATA COPYIN(domain%sfc_types)

  END FUNCTION t_domain_constructor

END MODULE mo_tmx_field_class
