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

! Variable management for the TMX package
!
! This module provides the variable and variable list data structures
! for the turbulent mixing (TMX) package. It enables type-safe memory
! management for scalar, vector, and multi-dimensional field variables.
!
! The module contains one main data structure:
! - t_tmx_var: Variable container that extends var_descriptor
!
! Variables are stored in memory manager (memman) and can be accessed
! via pointers of various ranks and types (real, integer, logical).
MODULE mo_tmx_var

  USE iso_c_binding
  USE memman, ONLY: &
    & var_descriptor,                           &
    & add_var_real, add_var_integer,            &
    & allocate_var_real, allocate_var_integer,  &
    & get_var_data,                             &
    & mm_host_device_uid,                       &
    & mm_invalid_device_uid,                    &
    & mm_get_gpu_device_uid !, &
    ! & copy_to_device
  USE mo_kind, ONLY: wp, vp, i4, i1
  USE mo_exception, ONLY: finish, message
  USE mo_timer, ONLY: new_timer, timer_start, timer_stop
  USE mo_util_string, ONLY: int2string

#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_tmx_var, t_tmx_var_p

  !----------------------------------------------------------------------
  ! Type definitions
  !----------------------------------------------------------------------

  ! Polymorphic pointer container for different kinds and ranks
  !
  ! This container holds pointers to data of different types and ranks,
  ! allowing centralized pointer management across different data types.
  TYPE :: t_ptr_container
    ! Real wp pointers
    REAL(wp), POINTER :: r0d => NULL()
    REAL(wp), POINTER :: r1d(:) => NULL()
    REAL(wp), POINTER :: r2d(:,:) => NULL()
    REAL(wp), POINTER :: r3d(:,:,:) => NULL()
    REAL(wp), POINTER :: r4d(:,:,:,:) => NULL()
    REAL(wp), POINTER :: r5d(:,:,:,:,:) => NULL()

    ! Real vp pointers
    REAL(vp), POINTER :: v0d => NULL()
    REAL(vp), POINTER :: v1d(:) => NULL()
    REAL(vp), POINTER :: v2d(:,:) => NULL()
    REAL(vp), POINTER :: v3d(:,:,:) => NULL()
    REAL(vp), POINTER :: v4d(:,:,:,:) => NULL()
    REAL(vp), POINTER :: v5d(:,:,:,:,:) => NULL()

    ! Integer pointers
    INTEGER, POINTER :: i0d => NULL()
    INTEGER, POINTER :: i1d(:) => NULL()
    INTEGER, POINTER :: i2d(:,:) => NULL()
    INTEGER, POINTER :: i3d(:,:,:) => NULL()
    INTEGER, POINTER :: i4d(:,:,:,:) => NULL()
    INTEGER, POINTER :: i5d(:,:,:,:,:) => NULL()

    ! Logical pointers
    INTEGER(i1), POINTER :: l0d => NULL()
    INTEGER(i1), POINTER :: l5d(:,:,:,:,:) => NULL()
  END TYPE t_ptr_container

  ! TMX variable container class that extends var_descriptor
  !
  ! This class provides a wrapper around the memory manager's var_descriptor,
  ! adding TMX-specific functionality for variable initialization, access,
  ! and pointer management. It supports scalar, vector, and multi-dimensional
  ! field variables of different types (real, integer, logical).
  TYPE, EXTENDS(var_descriptor) :: t_tmx_var
    ! From base type:
    ! CHARACTER(kind=c_char, len=:), ALLOCATABLE :: name
    ! INTEGER :: patch_id
    ! INTEGER :: hgrid_id
    ! INTEGER :: vgrid_id
    ! INTEGER :: time_id
    CHARACTER(len=:), ALLOCATABLE :: type_id     !< Variable type ('double', 'mixed', 'integer', or 'logical')
    INTEGER                       :: rank        !< Dimensionality of the variable (0-5)
    INTEGER, ALLOCATABLE          :: dims(:)     !< Dimensions of the variable
    INTEGER                       :: ref_pos = -1 !< Position of first reference dimension
    INTEGER                       :: ref = -1     !< Index of first reference dimension
    INTEGER                       :: ref2_pos = -1 !< Position of second reference dimension
    INTEGER                       :: ref2 = -1     !< Index of second reference dimension
    LOGICAL                       :: is_attached = .FALSE. !< Flag indicating if var is attached to memory
    LOGICAL                       :: is_inizialized = .FALSE. !< Flag indicating if var is initialized
  CONTAINS
    PROCEDURE :: Init                  => init_tmx_var
    PROCEDURE :: Update                => t_tmx_var_update_var
    PROCEDURE :: Attach                => t_tmx_var_attach_var
    PROCEDURE, PRIVATE :: Add_var      => t_tmx_var_add_to_mmgr
    PROCEDURE, PRIVATE :: Allocate     => t_tmx_var_allocate

    ! Assignment procedures
    PROCEDURE :: Assign_r0d            => t_tmx_var_assign_r0d
    PROCEDURE :: Assign_i0d            => t_tmx_var_assign_i0d
    PROCEDURE :: Assign_l0d            => t_tmx_var_assign_l0d

    ! Internal helper procedures for pointer management
    PROCEDURE, PRIVATE :: Get_base_ptr => t_tmx_var_get_base_ptr
    PROCEDURE, PRIVATE :: Slice_ptr    => t_tmx_var_slice_ptr

    ! Pointer-getting procedures for real values (wp kind)
    PROCEDURE :: Get_ptr_r0d           => t_tmx_var_get_ptr_r0d
    PROCEDURE :: Get_ptr_r1d           => t_tmx_var_get_ptr_r1d
    PROCEDURE :: Get_ptr_r2d           => t_tmx_var_get_ptr_r2d
    PROCEDURE :: Get_ptr_r3d           => t_tmx_var_get_ptr_r3d
    PROCEDURE :: Get_ptr_r4d           => t_tmx_var_get_ptr_r4d

    ! Pointer-getting procedures for real values (vp kind)
    PROCEDURE :: Get_ptr_v0d           => t_tmx_var_get_ptr_v0d
    PROCEDURE :: Get_ptr_v1d           => t_tmx_var_get_ptr_v1d
    PROCEDURE :: Get_ptr_v2d           => t_tmx_var_get_ptr_v2d
    PROCEDURE :: Get_ptr_v3d           => t_tmx_var_get_ptr_v3d
    PROCEDURE :: Get_ptr_v4d           => t_tmx_var_get_ptr_v4d

    ! Pointer-getting procedures for integer values
    PROCEDURE :: Get_ptr_i0d           => t_tmx_var_get_ptr_i0d
    PROCEDURE :: Get_ptr_i1d           => t_tmx_var_get_ptr_i1d
    PROCEDURE :: Get_ptr_i2d           => t_tmx_var_get_ptr_i2d
    PROCEDURE :: Get_ptr_i3d           => t_tmx_var_get_ptr_i3d
    PROCEDURE :: Get_ptr_i4d           => t_tmx_var_get_ptr_i4d

    ! Pointer-getting procedures for logical values
    PROCEDURE :: Get_ptr_l0d           => t_tmx_var_get_ptr_l0d
  END TYPE t_tmx_var

  INTERFACE t_tmx_var
    MODULE PROCEDURE t_tmx_var_constructor
  END INTERFACE

  ! Pointer wrapper for t_tmx_var
  !
  ! This type provides a container for pointers to t_tmx_var objects.
  ! It is used in the variable list to store references to variables.
  TYPE t_tmx_var_p
    TYPE(t_tmx_var), POINTER :: p !< Pointer to a t_tmx_var object
  END TYPE t_tmx_var_p

  LOGICAL, TARGET :: true_value = .TRUE.
  LOGICAL, TARGET :: false_value = .FALSE.

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_var'

  INTEGER, PARAMETER :: host_device = mm_host_device_uid
  INTEGER :: gpu_device = mm_invalid_device_uid
  LOGICAL :: mmgr_require_init = .TRUE.

  INTEGER :: timer_tmx_mmgr

CONTAINS

  ! preliminary helper:
  SUBROUTINE error_stop(line)
    INTEGER, INTENT(in) :: line
    PRINT*,'error at line', line
    CALL finish('mo_tmx_var::error_stop', 'Mmgr returned non-zero error code.')
  END SUBROUTINE error_stop

  ! preliminary helper:
  SUBROUTINE checkz(ierr, line)
    INTEGER, INTENT(in) :: ierr, line
    IF (ierr /= 0) CALL error_stop(line)
  END SUBROUTINE checkz

  ! preliminary helper:
  SUBROUTINE init_mmgr
    CHARACTER(len=*), PARAMETER :: routine = 'mo_tmx_var::init_mmgr'
    IF (.NOT. mmgr_require_init) RETURN
    gpu_device = mm_get_gpu_device_uid(0)

    timer_tmx_mmgr = new_timer('tmx_mmgr')

    mmgr_require_init = .FALSE.
  END SUBROUTINE init_mmgr

  !----------------------------------------------------------------------
  ! t_tmx_var implementation
  !----------------------------------------------------------------------

  ! Constructor for t_tmx_var type
  !
  ! Creates and initializes a new t_tmx_var object with the given name,
  ! type, and dimensions. If a var_descriptor is provided, the new variable
  ! will be attached to it.
  FUNCTION t_tmx_var_constructor(name, type_id, var_desc, rank, dims, patch_id, ref_pos) RESULT(var)

    CHARACTER(len=*),           INTENT(in)           :: name     !< Name of the variable
    CHARACTER(len=*),           INTENT(in)           :: type_id  !< Type of the variable ('double', 'mixed', 'integer', or 'logical')
    TYPE(var_descriptor),       INTENT(in), OPTIONAL :: var_desc !< Optional var_descriptor to attach to
    INTEGER,                    INTENT(in), OPTIONAL :: rank     !< Dimensionality of the variable (0-5)
    INTEGER,                    INTENT(in), OPTIONAL :: dims(:)  !< Dimensions of the variable
    INTEGER,                    INTENT(in), OPTIONAL :: patch_id !< ID of the patch the variable belongs to
    INTEGER,                    INTENT(in), OPTIONAL :: ref_pos  !< Position of first reference dimension
    TYPE(t_tmx_var), POINTER                         :: var      !< Pointer to the created t_tmx_var

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_tmx_var_constructor'

    ALLOCATE(var)
    CALL var%Init(name, type_id, rank, dims, patch_id, ref_pos)

    IF (PRESENT(var_desc)) THEN
      CALL var%Attach(var_desc, ref_pos=ref_pos)
      CALL var%Init(name, type_id)
    END IF

  END FUNCTION t_tmx_var_constructor

  ! Initialize a t_tmx_var object with the given name, type, and dimensions
  !
  ! This subroutine initializes the t_tmx_var object and allocates memory for
  ! it if necessary. It can be called either to set up the initial attributes
  ! or to create the actual memory allocation.
  !
  ! this    The t_tmx_var object to initialize
  ! name    Name of the variable
  ! type_id Type of the variable ('double', 'mixed', 'integer', or 'logical')
  ! rank    Dimensionality of the variable (0-5)
  ! dims    Dimensions of the variable
  ! patch_id ID of the patch the variable belongs to
  ! ref_pos Position of first reference dimension
  SUBROUTINE init_tmx_var(this, name, type_id, rank, dims, patch_id, ref_pos)

    CLASS(t_tmx_var),           INTENT(inout)        :: this    !< The t_tmx_var object to initialize
    CHARACTER(len=*),           INTENT(in)           :: name    !< Name of the variable
    CHARACTER(len=*),           INTENT(in)           :: type_id !< Type of the variable ('double', 'mixed', 'integer', or 'logical')
    INTEGER,                    INTENT(in), OPTIONAL :: rank    !< Dimensionality of the variable (0-5)
    INTEGER,                    INTENT(in), OPTIONAL :: dims(:) !< Dimensions of the variable
    INTEGER,                    INTENT(in), OPTIONAL :: patch_id !< ID of the patch the variable belongs to
    INTEGER,                    INTENT(in), OPTIONAL :: ref_pos !< Position of first reference dimension

    INTEGER :: istat, id
    REAL(wp), POINTER :: ptr_r5d(:,:,:,:,:)
    REAL(vp), POINTER :: ptr_v5d(:,:,:,:,:)
    INTEGER,  POINTER :: ptr_i5d(:,:,:,:,:)
    INTEGER(i1),  POINTER :: ptr_l5d(:,:,:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_tmx_var'

    IF (.NOT. this%is_inizialized) THEN
      ! First phase: Initialize attributes

      IF (mmgr_require_init) CALL init_mmgr

      id = 1
      IF (PRESENT(patch_id)) id = patch_id

      this%type_id = type_id

      CALL this%Update(name, patch_id=id, hgrid_id=1, vgrid_id=1, time_id=-1, ref_pos=ref_pos)

      IF (PRESENT(rank) .AND. PRESENT(dims)) THEN
        IF (SIZE(dims) == rank) THEN
          this%rank = rank
          ALLOCATE(this%dims(this%rank))
          IF (rank > 0) this%dims(:) = dims(:)
        ELSE
          CALL finish(routine, 'Inconsistent rank and dims arguments')
        END IF
      ELSE IF (PRESENT(dims)) THEN
        this%rank = SIZE(dims)
        ALLOCATE(this%dims(this%rank))
        IF (this%rank > 0) this%dims(:) = dims(:)
      ELSE IF (PRESENT(rank)) THEN
        this%rank = rank
        ALLOCATE(this%dims(rank))
      ELSE
        CALL finish(routine, 'rank and/or dims argument required')
      END IF

      this%is_inizialized = .TRUE.

    ELSE
      ! Second phase: Allocate memory if not already attached

      SELECT CASE (this%type_id)
      CASE ('double')
        istat = get_var_data(ptr_r5d, this%var_descriptor)
      CASE ('mixed')
        istat = get_var_data(ptr_v5d, this%var_descriptor)
      CASE ('integer')
        istat = get_var_data(ptr_i5d, this%var_descriptor)
      CASE ('logical')
        istat = get_var_data(ptr_l5d, this%var_descriptor)
      END SELECT

      IF (.NOT. this%is_attached .OR. istat /= 0) THEN

        CALL this%Add_var(this%dims)
        CALL this%Allocate()

        this%is_attached = .TRUE.

      END IF

      SELECT CASE (this%type_id)
      CASE ('double')
        istat = get_var_data(ptr_r5d, this%var_descriptor)
      CASE ('mixed')
        istat = get_var_data(ptr_v5d, this%var_descriptor)
      CASE ('integer')
        istat = get_var_data(ptr_i5d, this%var_descriptor)
      CASE ('logical')
        istat = get_var_data(ptr_l5d, this%var_descriptor)
      END SELECT
      IF (istat /= 0) CALL finish(routine, 'Construction of tmx_var '//TRIM(this%name)//' failed.')

    END IF

  END SUBROUTINE init_tmx_var

  ! Add variable to memory manager
  !
  ! Creates a variable in the memory manager with the appropriate type and dimensions
  SUBROUTINE t_tmx_var_add_to_mmgr(this, dims)

    CLASS(t_tmx_var), INTENT(in) :: this !< The t_tmx_var object
    INTEGER,          INTENT(in) :: dims(:) !< Dimensions of the variable

    TYPE(var_descriptor) :: var_desc
    INTEGER              :: mm_dims(5)
    INTEGER              :: ndims, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_tmx_var_add_to_mmgr'

    ndims = SIZE(dims)
    mm_dims(1:ndims) = dims(:)
    mm_dims(ndims+1:5) = 1

    SELECT CASE (this%type_id)
    CASE ('double')
      istat = add_var_real(this%var_descriptor, wp, mm_dims)
    CASE ('mixed')
      istat = add_var_real(this%var_descriptor, vp, mm_dims)
    CASE ('integer')
      istat = add_var_integer(this%var_descriptor, i4, mm_dims)
    CASE ('logical')
      istat = add_var_integer(this%var_descriptor, i1, mm_dims)
    END SELECT
    IF (istat /= 0) CALL finish(routine, 'Adding '//this%name//' failed')

  END SUBROUTINE t_tmx_var_add_to_mmgr

  ! Allocate memory for the variable
  !
  ! Allocates memory for the variable in the memory manager
  SUBROUTINE t_tmx_var_allocate(this)
#ifdef _OPENACC
    USE mo_openacc, ONLY: acc_map_data
#endif
    CLASS(t_tmx_var), INTENT(in) :: this !< The t_tmx_var object

    INTEGER :: istat
    REAL(wp), PARAMETER :: one_wp = 1.0
    INTEGER(i4), PARAMETER :: one_i4 = 1
    INTEGER(i1), PARAMETER :: one_i1 = 1
    INTEGER, PARAMETER :: wp_size = storage_size(one_wp)/8
    INTEGER, PARAMETER :: i4_size = bit_size(one_i4)/8
    INTEGER, PARAMETER :: i1_size = bit_size(one_i1)/8
    REAL(wp), POINTER :: ptr_r5d_h(:,:,:,:,:)
    REAL(wp), POINTER :: ptr_r5d_d(:,:,:,:,:)
    REAL(vp), POINTER :: ptr_v5d_h(:,:,:,:,:)
    REAL(vp), POINTER :: ptr_v5d_d(:,:,:,:,:)
    INTEGER, POINTER :: ptr_i5d_h(:,:,:,:,:)
    INTEGER, POINTER :: ptr_i5d_d(:,:,:,:,:)
    INTEGER(i1), POINTER :: ptr_l5d_h(:,:,:,:,:)
    INTEGER(i1), POINTER :: ptr_l5d_d(:,:,:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_tmx_var_allocate'

    SELECT CASE (this%type_id)
    CASE ('double')
      istat = allocate_var_real(this%var_descriptor, wp)
      CALL checkz(get_var_data(ptr_r5d_h, this%var_descriptor), __LINE__)
#ifdef _OPENACC
      CALL checkz(allocate_var_real(this%var_descriptor, wp, gpu_device), __LINE__)
      CALL checkz(get_var_data(ptr_r5d_d, this%var_descriptor, gpu_device), __LINE__)
      CALL acc_map_data(c_LOC(ptr_r5d_h), c_LOC(ptr_r5d_d), wp_size*SIZE(ptr_r5d_h))
#endif
    CASE ('mixed')
      istat = allocate_var_real(this%var_descriptor, vp)
      CALL checkz(get_var_data(ptr_v5d_h, this%var_descriptor), __LINE__)
#ifdef _OPENACC
      CALL checkz(allocate_var_real(this%var_descriptor, vp, gpu_device), __LINE__)
      CALL checkz(get_var_data(ptr_v5d_d, this%var_descriptor, gpu_device), __LINE__)
      CALL acc_map_data(c_LOC(ptr_v5d_h), c_LOC(ptr_v5d_d), wp_size*SIZE(ptr_v5d_h))
#endif
    CASE ('integer')
      istat = allocate_var_integer(this%var_descriptor, i4)
      CALL checkz(get_var_data(ptr_i5d_h, this%var_descriptor), __LINE__)
#ifdef _OPENACC
      CALL checkz(allocate_var_integer(this%var_descriptor, i4, gpu_device), __LINE__)
      CALL checkz(get_var_data(ptr_i5d_d, this%var_descriptor, gpu_device), __LINE__)
      CALL acc_map_data(c_LOC(ptr_i5d_h), c_LOC(ptr_i5d_d), i4_size*SIZE(ptr_i5d_h))
#endif
    CASE ('logical')
      istat = allocate_var_integer(this%var_descriptor, i1)
      CALL checkz(get_var_data(ptr_l5d_h, this%var_descriptor), __LINE__)
#ifdef _OPENACC
      CALL checkz(allocate_var_integer(this%var_descriptor, i1, gpu_device), __LINE__)
      CALL checkz(get_var_data(ptr_l5d_d, this%var_descriptor, gpu_device), __LINE__)
      CALL acc_map_data(c_LOC(ptr_l5d_h), c_LOC(ptr_l5d_d), i1_size*SIZE(ptr_l5d_h))
#endif
    END SELECT
    IF (istat /= 0) CALL finish(routine, 'Allocation of '//this%name//' failed')

  END SUBROUTINE t_tmx_var_allocate

  ! Update variable attributes
  !
  ! Updates the attributes of a t_tmx_var object
  SUBROUTINE t_tmx_var_update_var(this, name, patch_id, hgrid_id, vgrid_id, time_id, &
    &                             ref_pos, ref, ref2_pos, ref2)

    CLASS(t_tmx_var), INTENT(inout)        :: this      !< The t_tmx_var object to update
    CHARACTER(len=*), INTENT(in), OPTIONAL :: name      !< New name for the variable
    INTEGER,          INTENT(in), OPTIONAL :: patch_id  !< New patch ID
    INTEGER,          INTENT(in), OPTIONAL :: hgrid_id  !< New horizontal grid ID
    INTEGER,          INTENT(in), OPTIONAL :: vgrid_id  !< New vertical grid ID
    INTEGER,          INTENT(in), OPTIONAL :: time_id   !< New time ID
    INTEGER,          INTENT(in), OPTIONAL :: ref_pos   !< New position for first reference dimension
    INTEGER,          INTENT(in), OPTIONAL :: ref       !< New index for first reference dimension
    INTEGER,          INTENT(in), OPTIONAL :: ref2_pos  !< New position for second reference dimension
    INTEGER,          INTENT(in), OPTIONAL :: ref2      !< New index for second reference dimension

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_tmx_var_update_var'

    IF (PRESENT(name))     this%name     = name
    IF (PRESENT(patch_id)) this%patch_id = patch_id
    IF (PRESENT(hgrid_id)) this%hgrid_id = hgrid_id
    IF (PRESENT(vgrid_id)) this%vgrid_id = vgrid_id
    IF (PRESENT(time_id))  this%time_id  = time_id

    IF (PRESENT(ref_pos)) this%ref_pos = ref_pos
    IF (PRESENT(ref))     this%ref     = ref
    IF (PRESENT(ref2_pos)) THEN
      IF (.NOT. PRESENT(ref_pos)) THEN
        CALL finish(routine, '*ref2_pos* argument requires *ref_pos* argument')
      ELSE
        IF (ref2_pos == ref_pos) CALL finish(routine, 'Invalid *ref2_pos* argument')
        IF (ref2_pos < ref_pos) THEN
          this%ref_pos = ref2_pos
          IF (PRESENT(ref2)) this%ref = ref2
          this%ref2_pos = ref_pos
          IF (PRESENT(ref)) this%ref2 = ref
        ELSE
          this%ref2_pos = ref2_pos
          IF (PRESENT(ref2)) this%ref2 = ref2
        END IF
      END IF
    END IF

    this%is_attached = .TRUE.

  END SUBROUTINE t_tmx_var_update_var

  ! Attach a variable to a var_descriptor
  !
  ! Links the t_tmx_var to an existing var_descriptor
  SUBROUTINE t_tmx_var_attach_var(this, var_desc, ref_pos, ref, ref2_pos, ref2)

    CLASS(t_tmx_var),     INTENT(inout) :: this     !< The t_tmx_var object
    TYPE(var_descriptor), INTENT(in)    :: var_desc !< The var_descriptor to attach to
    INTEGER, OPTIONAL,    INTENT(in)    :: ref_pos  !< Position of first reference dimension
    INTEGER, OPTIONAL,    INTENT(in)    :: ref      !< Index for first reference dimension
    INTEGER, OPTIONAL,    INTENT(in)    :: ref2_pos !< Position of second reference dimension
    INTEGER, OPTIONAL,    INTENT(in)    :: ref2     !< Index for second reference dimension

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_tmx_var_attach_var'

    CALL this%Update(var_desc%name, var_desc%patch_id, var_desc%hgrid_id, &
      & var_desc%vgrid_id, var_desc%time_id, ref_pos=ref_pos, ref=ref, &
      & ref2_pos=ref2_pos, ref2=ref2)

    this%is_attached = .TRUE.

  END SUBROUTINE t_tmx_var_attach_var

  !----------------------------------------------------------------------
  ! Assignment methods for t_tmx_var
  !----------------------------------------------------------------------

  ! Assign real scalar value to variable
  SUBROUTINE t_tmx_var_assign_r0d(this, value)

    CLASS(t_tmx_var), INTENT(inout) :: this  !< The t_tmx_var object
    REAL(wp), INTENT(in) :: value           !< The real value to assign

    REAL(wp), POINTER :: ptr_r5d(:,:,:,:,:)
    INTEGER :: istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_tmx_var_assign_r0d'

    istat = -1
    ptr_r5d => NULL()

    istat = get_var_data(ptr_r5d, this%var_descriptor)
    IF (istat /= 0) THEN
      CALL this%Add_var(this%dims)
      CALL this%Allocate()
      this%is_attached = .TRUE.
    END IF
    istat = get_var_data(ptr_r5d, this%var_descriptor)
    ptr_r5d(1,1,1,1,1) = value
    !$ACC UPDATE DEVICE(ptr_r5d)

  END SUBROUTINE t_tmx_var_assign_r0d

  ! Assign integer scalar value to variable
  SUBROUTINE t_tmx_var_assign_i0d(this, value)

    CLASS(t_tmx_var), INTENT(inout) :: this  !< The t_tmx_var object
    INTEGER, INTENT(in) :: value             !< The integer value to assign

    INTEGER, POINTER :: ptr_i5d(:,:,:,:,:)
    INTEGER :: istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_tmx_var_assign_i0d'

    istat = -1
    ptr_i5d => NULL()

    istat = get_var_data(ptr_i5d, this%var_descriptor)
    IF (istat /= 0) THEN
      CALL this%Add_var(this%dims)
      CALL this%Allocate()
      this%is_attached = .TRUE.
    END IF
    istat = get_var_data(ptr_i5d, this%var_descriptor)
    ptr_i5d(1,1,1,1,1) = value
    !$ACC UPDATE DEVICE(ptr_i5d)

  END SUBROUTINE t_tmx_var_assign_i0d

  ! Assign logical scalar value to variable
  SUBROUTINE t_tmx_var_assign_l0d(this, value)

    CLASS(t_tmx_var), INTENT(inout) :: this  !< The t_tmx_var object
    LOGICAL, INTENT(in) :: value             !< The logical value to assign

    INTEGER(i1), POINTER :: ptr_l0d, ptr_l5d(:,:,:,:,:)
    INTEGER :: istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_tmx_var_assign_l0d'

    istat = -1
    ptr_l5d => NULL()

    istat = get_var_data(ptr_l5d, this%var_descriptor)
    IF (istat /= 0) THEN
      CALL this%Add_var(this%dims)
      CALL this%Allocate()
      this%is_attached = .TRUE.
    END IF
    istat = get_var_data(ptr_l5d, this%var_descriptor)
    IF (value) THEN
      ptr_l5d(1,1,1,1,1) = 1_i1
    ELSE
      ptr_l5d(1,1,1,1,1) = 0_i1
    END IF

    !$ACC UPDATE DEVICE(ptr_l5d)

  END SUBROUTINE t_tmx_var_assign_l0d

  !----------------------------------------------------------------------
  ! Internal helper procedures for centralized pointer management
  !----------------------------------------------------------------------

  ! Get base pointer from memory manager
  !
  ! Retrieves the 5D base pointer from the memory manager for the variable.
  ! This centralizes all interaction with the memory manager's get_var_data.
  SUBROUTINE t_tmx_var_get_base_ptr(this, ptr_container)

    CLASS(t_tmx_var), INTENT(in)       :: this          !< The t_tmx_var object
    TYPE(t_ptr_container), INTENT(out) :: ptr_container !< Container to receive pointers

    INTEGER :: istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_tmx_var_get_base_ptr'

    CALL timer_start(timer_tmx_mmgr)

    SELECT CASE (this%type_id)
    CASE ('double')
      istat = get_var_data(ptr_container%r5d, this%var_descriptor)
    CASE ('mixed')
      istat = get_var_data(ptr_container%v5d, this%var_descriptor)
    CASE ('integer')
      istat = get_var_data(ptr_container%i5d, this%var_descriptor)
    CASE ('logical')
      istat = get_var_data(ptr_container%l5d, this%var_descriptor)
    CASE DEFAULT
      istat = -1
    END SELECT

    IF (istat /= 0) CALL finish(routine, "Couldn't get data for variable "//this%name)

    CALL timer_stop(timer_tmx_mmgr)

  END SUBROUTINE t_tmx_var_get_base_ptr

  ! Slice pointer to target rank
  !
  ! Performs array slicing based on ref_pos and ref indices to reduce the
  ! dimensionality of the pointer from 5D to the target rank. This centralizes
  ! all the slicing logic that was previously duplicated across functions.
  SUBROUTINE t_tmx_var_slice_ptr(this, ptr_container, target_rank, ref, ref2)

    CLASS(t_tmx_var), INTENT(in)          :: this          !< The t_tmx_var object
    TYPE(t_ptr_container), INTENT(inout)  :: ptr_container !< Container with pointers to slice
    INTEGER, INTENT(in)                   :: target_rank   !< Target rank for the result
    INTEGER, OPTIONAL, INTENT(in)         :: ref           !< Optional reference index
    INTEGER, OPTIONAL, INTENT(in)         :: ref2          !< Optional second reference index

    INTEGER :: ref_idx, ref2_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_tmx_var_slice_ptr'

    ! Determine reference indices with safeguards
    ref_idx = -1
    ref2_idx = -1

    IF (PRESENT(ref)) THEN
      IF (this%ref /= -1) THEN
        ! If this%ref is already set, treat the ref parameter as ref2_idx
        ref2_idx = ref
        IF (PRESENT(ref2)) THEN
          CALL finish(routine, 'Cannot provide *ref2* when *ref* parameter is used as second reference for '//this%name)
        END IF
      ELSE
        ! Normal case: ref parameter sets ref_idx
        ref_idx = ref
      END IF
    END IF

    IF (PRESENT(ref2)) THEN
      IF (.NOT. PRESENT(ref)) THEN
        CALL finish(routine, '*ref2* argument requires *ref* argument')
      END IF
      IF (this%ref2 /= -1) THEN
        CALL finish(routine, 'Cannot override *ref2* for '//this%name//': already set to '// &
          &                  int2string(this%ref2))
      END IF
      ref2_idx = ref2
    END IF

    ! Validate slicing request
    IF (ref_idx /= -1 .AND. ref2_idx /= -1) THEN
      IF (this%rank /= target_rank + 2) THEN
        CALL finish(routine, 'Invalid slicing: rank mismatch with two refs for '//this%name//': '// &
          &                  int2string(this%rank)//' '//int2string(target_rank)//'+2')
      END IF
    ELSE IF (ref_idx /= -1) THEN
      IF (this%rank /= target_rank + 1) THEN
        CALL finish(routine, 'Invalid slicing: rank mismatch with one ref for '//this%name//': '// &
          &                  int2string(this%rank)//' '//int2string(target_rank)//'+1')
      END IF
    ELSE
      IF (this%rank /= target_rank) THEN
        CALL finish(routine, 'Invalid slicing: rank mismatch for '//this%name//': '// &
          &                  int2string(this%rank)//' '//int2string(target_rank))
      END IF
    END IF

    ! Use stored references if not provided (may still be -1)
    IF (ref_idx == -1) ref_idx = this%ref
    IF (ref2_idx == -1) ref2_idx = this%ref2

    ! Perform the actual slicing based on type
    SELECT CASE (this%type_id)
    CASE ('double')
      CALL slice_real_wp_ptr(ptr_container, target_rank, this%ref_pos, ref_idx, &
        & this%ref2_pos, ref2_idx)
    CASE ('mixed')
      CALL slice_real_vp_ptr(ptr_container, target_rank, this%ref_pos, ref_idx, &
        & this%ref2_pos, ref2_idx)
    CASE ('integer')
      CALL slice_integer_ptr(ptr_container, target_rank, this%ref_pos, ref_idx, &
        & this%ref2_pos, ref2_idx)
    CASE ('logical')
      CALL slice_logical_ptr(ptr_container, target_rank, this%ref_pos, ref_idx, &
        & this%ref2_pos, ref2_idx)
    END SELECT

  END SUBROUTINE t_tmx_var_slice_ptr

  ! Helper to slice real(wp) pointers
  SUBROUTINE slice_real_wp_ptr(ptrs, target_rank, ref_pos, ref_idx, ref2_pos, ref2_idx)

    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: target_rank, ref_pos, ref_idx, ref2_pos, ref2_idx

    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_real_wp_ptr'

    IF (ref_idx /= -1 .AND. ref2_idx /= -1) THEN
      ! Two-level slicing
      SELECT CASE (target_rank)
      CASE (1)
        CALL slice_wp_5d_to_1d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
      CASE (2)
        CALL slice_wp_5d_to_2d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
      CASE (3)
        CALL slice_wp_5d_to_3d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
      END SELECT
    ELSE IF (ref_idx /= -1) THEN
      ! Single-level slicing
      SELECT CASE (target_rank)
      CASE (0)
        ptrs%r0d => ptrs%r5d(1,1,1,1,1)
      CASE (1)
        CALL slice_wp_5d_to_1d_single_ref(ptrs, ref_pos, ref_idx)
      CASE (2)
        CALL slice_wp_5d_to_2d_single_ref(ptrs, ref_pos, ref_idx)
      CASE (3)
        CALL slice_wp_5d_to_3d_single_ref(ptrs, ref_pos, ref_idx)
      CASE (4)
        CALL slice_wp_5d_to_4d_single_ref(ptrs, ref_pos, ref_idx)
      END SELECT
    ELSE
      ! No slicing, just dimension reduction
      SELECT CASE (target_rank)
      CASE (0)
        ptrs%r0d => ptrs%r5d(1,1,1,1,1)
      CASE (1)
        ptrs%r1d => ptrs%r5d(:,1,1,1,1)
      CASE (2)
        ptrs%r2d => ptrs%r5d(:,:,1,1,1)
      CASE (3)
        ptrs%r3d => ptrs%r5d(:,:,:,1,1)
      CASE (4)
        ptrs%r4d => ptrs%r5d(:,:,:,:,1)
      END SELECT
    END IF

  END SUBROUTINE slice_real_wp_ptr

  ! Helper to slice real(vp) pointers
  SUBROUTINE slice_real_vp_ptr(ptrs, target_rank, ref_pos, ref_idx, ref2_pos, ref2_idx)

    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: target_rank, ref_pos, ref_idx, ref2_pos, ref2_idx

    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_real_vp_ptr'

    IF (ref_idx /= -1 .AND. ref2_idx /= -1) THEN
      ! Two-level slicing
      SELECT CASE (target_rank)
      CASE (1)
        CALL slice_vp_5d_to_1d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
      CASE (2)
        CALL slice_vp_5d_to_2d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
      CASE (3)
        CALL slice_vp_5d_to_3d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
      END SELECT
    ELSE IF (ref_idx /= -1) THEN
      ! Single-level slicing
      SELECT CASE (target_rank)
      CASE (0)
        ptrs%v0d => ptrs%v5d(1,1,1,1,1)
      CASE (1)
        CALL slice_vp_5d_to_1d_single_ref(ptrs, ref_pos, ref_idx)
      CASE (2)
        CALL slice_vp_5d_to_2d_single_ref(ptrs, ref_pos, ref_idx)
      CASE (3)
        CALL slice_vp_5d_to_3d_single_ref(ptrs, ref_pos, ref_idx)
      CASE (4)
        CALL slice_vp_5d_to_4d_single_ref(ptrs, ref_pos, ref_idx)
      END SELECT
    ELSE
      ! No slicing, just dimension reduction
      SELECT CASE (target_rank)
      CASE (0)
        ptrs%v0d => ptrs%v5d(1,1,1,1,1)
      CASE (1)
        ptrs%v1d => ptrs%v5d(:,1,1,1,1)
      CASE (2)
        ptrs%v2d => ptrs%v5d(:,:,1,1,1)
      CASE (3)
        ptrs%v3d => ptrs%v5d(:,:,:,1,1)
      CASE (4)
        ptrs%v4d => ptrs%v5d(:,:,:,:,1)
      END SELECT
    END IF

  END SUBROUTINE slice_real_vp_ptr

  ! Helper to slice integer pointers
  SUBROUTINE slice_integer_ptr(ptrs, target_rank, ref_pos, ref_idx, ref2_pos, ref2_idx)

    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: target_rank, ref_pos, ref_idx, ref2_pos, ref2_idx

    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_integer_ptr'

    IF (ref_idx /= -1 .AND. ref2_idx /= -1) THEN
      ! Two-level slicing
      SELECT CASE (target_rank)
      CASE (1)
        CALL slice_int_5d_to_1d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
      CASE (2)
        CALL slice_int_5d_to_2d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
      CASE (3)
        CALL slice_int_5d_to_3d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
      END SELECT
    ELSE IF (ref_idx /= -1) THEN
      ! Single-level slicing
      SELECT CASE (target_rank)
      CASE (0)
        ptrs%i0d => ptrs%i5d(1,1,1,1,1)
      CASE (1)
        CALL slice_int_5d_to_1d_single_ref(ptrs, ref_pos, ref_idx)
      CASE (2)
        CALL slice_int_5d_to_2d_single_ref(ptrs, ref_pos, ref_idx)
      CASE (3)
        CALL slice_int_5d_to_3d_single_ref(ptrs, ref_pos, ref_idx)
      CASE (4)
        CALL slice_int_5d_to_4d_single_ref(ptrs, ref_pos, ref_idx)
      END SELECT
    ELSE
      ! No slicing, just dimension reduction
      SELECT CASE (target_rank)
      CASE (0)
        ptrs%i0d => ptrs%i5d(1,1,1,1,1)
      CASE (1)
        ptrs%i1d => ptrs%i5d(:,1,1,1,1)
      CASE (2)
        ptrs%i2d => ptrs%i5d(:,:,1,1,1)
      CASE (3)
        ptrs%i3d => ptrs%i5d(:,:,:,1,1)
      CASE (4)
        ptrs%i4d => ptrs%i5d(:,:,:,:,1)
      END SELECT
    END IF

  END SUBROUTINE slice_integer_ptr

  ! Helper to slice logical pointers
  SUBROUTINE slice_logical_ptr(ptrs, target_rank, ref_pos, ref_idx, ref2_pos, ref2_idx)

    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: target_rank, ref_pos, ref_idx, ref2_pos, ref2_idx

    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_logical_ptr'

    IF (ref_idx /= -1 .OR. ref2_idx /= -1) THEN
      CALL finish(routine, 'Slicing not supported for logical pointers')
    ELSE
      ! No slicing, just dimension reduction
      SELECT CASE (target_rank)
      CASE (0)
        ptrs%l0d => ptrs%l5d(1,1,1,1,1)
      ! CASE (1)
      !   ptrs%l1d => ptrs%l5d(:,1,1,1,1)
      ! CASE (2)
      !   ptrs%l2d => ptrs%l5d(:,:,1,1,1)
      ! CASE (3)
      !   ptrs%l3d => ptrs%l5d(:,:,:,1,1)
      ! CASE (4)
      !   ptrs%l4d => ptrs%l5d(:,:,:,:,1)
      END SELECT
    END IF

  END SUBROUTINE slice_logical_ptr

  ! Detailed slicing subroutines for wp kind - single ref
  SUBROUTINE slice_wp_5d_to_1d_single_ref(ptrs, ref_pos, ref_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_wp_5d_to_1d_single_ref'

    IF (ref_pos < 1 .OR. ref_pos > 2) CALL finish(routine, 'Invalid *ref_pos*')
    IF (SIZE(ptrs%r5d,3) /= 1 .OR. SIZE(ptrs%r5d,4) /= 1 .OR. SIZE(ptrs%r5d,5) /= 1) THEN
      CALL finish(routine, 'Only 2-dim variable allowed with *ref*')
    END IF

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%r1d => ptrs%r5d(ref_idx,:,1,1,1)
    CASE(2)
      ptrs%r1d => ptrs%r5d(:,ref_idx,1,1,1)
    END SELECT
  END SUBROUTINE slice_wp_5d_to_1d_single_ref

  SUBROUTINE slice_wp_5d_to_2d_single_ref(ptrs, ref_pos, ref_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_wp_5d_to_2d_single_ref'

    IF (ref_pos < 1 .OR. ref_pos > 3) CALL finish(routine, 'Invalid *ref_pos*')
    IF (SIZE(ptrs%r5d,4) /= 1 .OR. SIZE(ptrs%r5d,5) /= 1) THEN
      CALL finish(routine, 'Only 3-dim variable allowed with *ref*')
    END IF

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%r2d => ptrs%r5d(ref_idx,:,:,1,1)
    CASE(2)
      ptrs%r2d => ptrs%r5d(:,ref_idx,:,1,1)
    CASE(3)
      ptrs%r2d => ptrs%r5d(:,:,ref_idx,1,1)
    END SELECT
  END SUBROUTINE slice_wp_5d_to_2d_single_ref

  SUBROUTINE slice_wp_5d_to_3d_single_ref(ptrs, ref_pos, ref_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_wp_5d_to_3d_single_ref'

    IF (ref_pos < 1 .OR. ref_pos > 4) CALL finish(routine, 'Invalid *ref_pos*')
    IF (SIZE(ptrs%r5d,5) /= 1) CALL finish(routine, 'Only 4-dim variable allowed with *ref*')

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%r3d => ptrs%r5d(ref_idx,:,:,:,1)
    CASE(2)
      ptrs%r3d => ptrs%r5d(:,ref_idx,:,:,1)
    CASE(3)
      ptrs%r3d => ptrs%r5d(:,:,ref_idx,:,1)
    CASE(4)
      ptrs%r3d => ptrs%r5d(:,:,:,ref_idx,1)
    END SELECT
  END SUBROUTINE slice_wp_5d_to_3d_single_ref

  SUBROUTINE slice_wp_5d_to_4d_single_ref(ptrs, ref_pos, ref_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_wp_5d_to_4d_single_ref'

    IF (SIZE(SHAPE(ptrs%r5d)) /= 5) CALL finish(routine, 'Only 5-dim variable allowed with *ref*')
    IF (ref_pos < 1 .OR. ref_pos > 5) CALL finish(routine, 'Invalid *ref_pos*')

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%r4d => ptrs%r5d(ref_idx,:,:,:,:)
    CASE(2)
      ptrs%r4d => ptrs%r5d(:,ref_idx,:,:,:)
    CASE(3)
      ptrs%r4d => ptrs%r5d(:,:,ref_idx,:,:)
    CASE(4)
      ptrs%r4d => ptrs%r5d(:,:,:,ref_idx,:)
    CASE(5)
      ptrs%r4d => ptrs%r5d(:,:,:,:,ref_idx)
    END SELECT
  END SUBROUTINE slice_wp_5d_to_4d_single_ref

  ! Detailed slicing subroutines for wp kind - double ref
  SUBROUTINE slice_wp_5d_to_1d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx, ref2_pos, ref2_idx
    REAL(wp), POINTER :: ptr_r2d(:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_wp_5d_to_1d_double_ref'

    IF (ref_pos < 1 .OR. ref_pos > 2) CALL finish(routine, 'Invalid *ref_pos*')
    IF (ref2_pos < 2 .OR. ref2_pos > 3) CALL finish(routine, 'Invalid *ref2_pos*')
    IF (SIZE(ptrs%r5d,5) /= 1 .OR. SIZE(ptrs%r5d,4) /= 1) THEN
      CALL finish(routine, 'Only 3-dim variable allowed with *ref* and *ref2*')
    END IF

    SELECT CASE (ref2_pos)
    CASE(2)
      ptr_r2d => ptrs%r5d(:,ref2_idx,:,1,1)
    CASE(3)
      ptr_r2d => ptrs%r5d(:,:,ref2_idx,1,1)
    END SELECT

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%r1d => ptr_r2d(ref_idx,:)
    CASE(2)
      ptrs%r1d => ptr_r2d(:,ref_idx)
    END SELECT
  END SUBROUTINE slice_wp_5d_to_1d_double_ref

  SUBROUTINE slice_wp_5d_to_2d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx, ref2_pos, ref2_idx
    REAL(wp), POINTER :: ptr_r3d(:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_wp_5d_to_2d_double_ref'

    IF (ref_pos < 1 .OR. ref_pos > 3) CALL finish(routine, 'Invalid *ref_pos*')
    IF (ref2_pos < 2 .OR. ref2_pos > 4) CALL finish(routine, 'Invalid *ref2_pos*')
    IF (SIZE(ptrs%r5d,5) /= 1) THEN
      CALL finish(routine, 'Only 4-dim variable allowed with *ref* and *ref2*')
    END IF

    SELECT CASE (ref2_pos)
    CASE(2)
      ptr_r3d => ptrs%r5d(:,ref2_idx,:,:,1)
    CASE(3)
      ptr_r3d => ptrs%r5d(:,:,ref2_idx,:,1)
    CASE(4)
      ptr_r3d => ptrs%r5d(:,:,:,ref2_idx,1)
    END SELECT

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%r2d => ptr_r3d(ref_idx,:,:)
    CASE(2)
      ptrs%r2d => ptr_r3d(:,ref_idx,:)
    CASE(3)
      ptrs%r2d => ptr_r3d(:,:,ref_idx)
    END SELECT
  END SUBROUTINE slice_wp_5d_to_2d_double_ref

  SUBROUTINE slice_wp_5d_to_3d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx, ref2_pos, ref2_idx
    REAL(wp), POINTER :: ptr_r4d(:,:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_wp_5d_to_3d_double_ref'

    IF (ref_pos < 1 .OR. ref_pos > 4) CALL finish(routine, 'Invalid *ref_pos*')
    IF (ref2_pos /= 5) CALL finish(routine, 'Invalid *ref2_pos*')
    IF (SIZE(SHAPE(ptrs%r5d)) /= 5) THEN
      CALL finish(routine, 'Only 5-dim variable allowed with *ref* and *ref2*')
    END IF

    SELECT CASE (ref2_pos)
    CASE(5)
      ptr_r4d => ptrs%r5d(:,:,:,:,ref2_idx)
    END SELECT

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%r3d => ptr_r4d(ref_idx,:,:,:)
    CASE(2)
      ptrs%r3d => ptr_r4d(:,ref_idx,:,:)
    CASE(3)
      ptrs%r3d => ptr_r4d(:,:,ref_idx,:)
    CASE(4)
      ptrs%r3d => ptr_r4d(:,:,:,ref_idx)
    END SELECT
  END SUBROUTINE slice_wp_5d_to_3d_double_ref

  ! Detailed slicing subroutines for vp kind - single ref
  SUBROUTINE slice_vp_5d_to_1d_single_ref(ptrs, ref_pos, ref_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_vp_5d_to_1d_single_ref'

    IF (ref_pos < 1 .OR. ref_pos > 2) CALL finish(routine, 'Invalid *ref_pos*')
    IF (SIZE(ptrs%v5d,3) /= 1 .OR. SIZE(ptrs%v5d,4) /= 1 .OR. SIZE(ptrs%v5d,5) /= 1) THEN
      CALL finish(routine, 'Only 2-dim variable allowed with *ref*')
    END IF

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%v1d => ptrs%v5d(ref_idx,:,1,1,1)
    CASE(2)
      ptrs%v1d => ptrs%v5d(:,ref_idx,1,1,1)
    END SELECT
  END SUBROUTINE slice_vp_5d_to_1d_single_ref

  SUBROUTINE slice_vp_5d_to_2d_single_ref(ptrs, ref_pos, ref_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_vp_5d_to_2d_single_ref'

    IF (ref_pos < 1 .OR. ref_pos > 3) CALL finish(routine, 'Invalid *ref_pos*')
    IF (SIZE(ptrs%v5d,4) /= 1 .OR. SIZE(ptrs%v5d,5) /= 1) THEN
      CALL finish(routine, 'Only 3-dim variable allowed with *ref*')
    END IF

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%v2d => ptrs%v5d(ref_idx,:,:,1,1)
    CASE(2)
      ptrs%v2d => ptrs%v5d(:,ref_idx,:,1,1)
    CASE(3)
      ptrs%v2d => ptrs%v5d(:,:,ref_idx,1,1)
    END SELECT
  END SUBROUTINE slice_vp_5d_to_2d_single_ref

  SUBROUTINE slice_vp_5d_to_3d_single_ref(ptrs, ref_pos, ref_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_vp_5d_to_3d_single_ref'

    IF (ref_pos < 1 .OR. ref_pos > 4) CALL finish(routine, 'Invalid *ref_pos*')
    IF (SIZE(ptrs%v5d,5) /= 1) CALL finish(routine, 'Only 4-dim variable allowed with *ref*')

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%v3d => ptrs%v5d(ref_idx,:,:,:,1)
    CASE(2)
      ptrs%v3d => ptrs%v5d(:,ref_idx,:,:,1)
    CASE(3)
      ptrs%v3d => ptrs%v5d(:,:,ref_idx,:,1)
    CASE(4)
      ptrs%v3d => ptrs%v5d(:,:,:,ref_idx,1)
    END SELECT
  END SUBROUTINE slice_vp_5d_to_3d_single_ref

  SUBROUTINE slice_vp_5d_to_4d_single_ref(ptrs, ref_pos, ref_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_vp_5d_to_4d_single_ref'

    IF (SIZE(SHAPE(ptrs%v5d)) /= 5) CALL finish(routine, 'Only 5-dim variable allowed with *ref*')
    IF (ref_pos < 1 .OR. ref_pos > 5) CALL finish(routine, 'Invalid *ref_pos*')

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%v4d => ptrs%v5d(ref_idx,:,:,:,:)
    CASE(2)
      ptrs%v4d => ptrs%v5d(:,ref_idx,:,:,:)
    CASE(3)
      ptrs%v4d => ptrs%v5d(:,:,ref_idx,:,:)
    CASE(4)
      ptrs%v4d => ptrs%v5d(:,:,:,ref_idx,:)
    CASE(5)
      ptrs%v4d => ptrs%v5d(:,:,:,:,ref_idx)
    END SELECT
  END SUBROUTINE slice_vp_5d_to_4d_single_ref

  ! Detailed slicing subroutines for vp kind - double ref
  SUBROUTINE slice_vp_5d_to_1d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx, ref2_pos, ref2_idx
    REAL(vp), POINTER :: ptr_v2d(:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_vp_5d_to_1d_double_ref'

    IF (ref_pos < 1 .OR. ref_pos > 2) CALL finish(routine, 'Invalid *ref_pos*')
    IF (ref2_pos < 2 .OR. ref2_pos > 3) CALL finish(routine, 'Invalid *ref2_pos*')
    IF (SIZE(ptrs%v5d,5) /= 1 .OR. SIZE(ptrs%v5d,4) /= 1) THEN
      CALL finish(routine, 'Only 3-dim variable allowed with *ref* and *ref2*')
    END IF

    SELECT CASE (ref2_pos)
    CASE(2)
      ptr_v2d => ptrs%v5d(:,ref2_idx,:,1,1)
    CASE(3)
      ptr_v2d => ptrs%v5d(:,:,ref2_idx,1,1)
    END SELECT

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%v1d => ptr_v2d(ref_idx,:)
    CASE(2)
      ptrs%v1d => ptr_v2d(:,ref_idx)
    END SELECT
  END SUBROUTINE slice_vp_5d_to_1d_double_ref

  SUBROUTINE slice_vp_5d_to_2d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx, ref2_pos, ref2_idx
    REAL(vp), POINTER :: ptr_v3d(:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_vp_5d_to_2d_double_ref'

    IF (ref_pos < 1 .OR. ref_pos > 3) CALL finish(routine, 'Invalid *ref_pos*')
    IF (ref2_pos < 2 .OR. ref2_pos > 4) CALL finish(routine, 'Invalid *ref2_pos*')
    IF (SIZE(ptrs%v5d,5) /= 1) THEN
      CALL finish(routine, 'Only 4-dim variable allowed with *ref* and *ref2*')
    END IF

    SELECT CASE (ref2_pos)
    CASE(2)
      ptr_v3d => ptrs%v5d(:,ref2_idx,:,:,1)
    CASE(3)
      ptr_v3d => ptrs%v5d(:,:,ref2_idx,:,1)
    CASE(4)
      ptr_v3d => ptrs%v5d(:,:,:,ref2_idx,1)
    END SELECT

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%v2d => ptr_v3d(ref_idx,:,:)
    CASE(2)
      ptrs%v2d => ptr_v3d(:,ref_idx,:)
    CASE(3)
      ptrs%v2d => ptr_v3d(:,:,ref_idx)
    END SELECT
  END SUBROUTINE slice_vp_5d_to_2d_double_ref

  SUBROUTINE slice_vp_5d_to_3d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx, ref2_pos, ref2_idx
    REAL(vp), POINTER :: ptr_v4d(:,:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_vp_5d_to_3d_double_ref'

    IF (ref_pos < 1 .OR. ref_pos > 4) CALL finish(routine, 'Invalid *ref_pos*')
    IF (ref2_pos /= 5) CALL finish(routine, 'Invalid *ref2_pos*')
    IF (SIZE(SHAPE(ptrs%v5d)) /= 5) THEN
      CALL finish(routine, 'Only 5-dim variable allowed with *ref* and *ref2*')
    END IF

    SELECT CASE (ref2_pos)
    CASE(5)
      ptr_v4d => ptrs%v5d(:,:,:,:,ref2_idx)
    END SELECT

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%v3d => ptr_v4d(ref_idx,:,:,:)
    CASE(2)
      ptrs%v3d => ptr_v4d(:,ref_idx,:,:)
    CASE(3)
      ptrs%v3d => ptr_v4d(:,:,ref_idx,:)
    CASE(4)
      ptrs%v3d => ptr_v4d(:,:,:,ref_idx)
    END SELECT
  END SUBROUTINE slice_vp_5d_to_3d_double_ref

  ! Detailed slicing subroutines for integer kind - single ref
  SUBROUTINE slice_int_5d_to_1d_single_ref(ptrs, ref_pos, ref_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_int_5d_to_1d_single_ref'

    IF (ref_pos < 1 .OR. ref_pos > 2) CALL finish(routine, 'Invalid *ref_pos*')
    IF (SIZE(ptrs%i5d,3) /= 1 .OR. SIZE(ptrs%i5d,4) /= 1 .OR. SIZE(ptrs%i5d,5) /= 1) THEN
      CALL finish(routine, 'Only 2-dim variable allowed with *ref*')
    END IF

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%i1d => ptrs%i5d(ref_idx,:,1,1,1)
    CASE(2)
      ptrs%i1d => ptrs%i5d(:,ref_idx,1,1,1)
    END SELECT
  END SUBROUTINE slice_int_5d_to_1d_single_ref

  SUBROUTINE slice_int_5d_to_2d_single_ref(ptrs, ref_pos, ref_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_int_5d_to_2d_single_ref'

    IF (ref_pos < 1 .OR. ref_pos > 3) CALL finish(routine, 'Invalid *ref_pos*')
    IF (SIZE(ptrs%i5d,4) /= 1 .OR. SIZE(ptrs%i5d,5) /= 1) THEN
      CALL finish(routine, 'Only 3-dim variable allowed with *ref*')
    END IF

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%i2d => ptrs%i5d(ref_idx,:,:,1,1)
    CASE(2)
      ptrs%i2d => ptrs%i5d(:,ref_idx,:,1,1)
    CASE(3)
      ptrs%i2d => ptrs%i5d(:,:,ref_idx,1,1)
    END SELECT
  END SUBROUTINE slice_int_5d_to_2d_single_ref

  SUBROUTINE slice_int_5d_to_3d_single_ref(ptrs, ref_pos, ref_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_int_5d_to_3d_single_ref'

    IF (ref_pos < 1 .OR. ref_pos > 4) CALL finish(routine, 'Invalid *ref_pos*')
    IF (SIZE(ptrs%i5d,5) /= 1) CALL finish(routine, 'Only 4-dim variable allowed with *ref*')

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%i3d => ptrs%i5d(ref_idx,:,:,:,1)
    CASE(2)
      ptrs%i3d => ptrs%i5d(:,ref_idx,:,:,1)
    CASE(3)
      ptrs%i3d => ptrs%i5d(:,:,ref_idx,:,1)
    CASE(4)
      ptrs%i3d => ptrs%i5d(:,:,:,ref_idx,1)
    END SELECT
  END SUBROUTINE slice_int_5d_to_3d_single_ref

  SUBROUTINE slice_int_5d_to_4d_single_ref(ptrs, ref_pos, ref_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_int_5d_to_4d_single_ref'

    IF (SIZE(SHAPE(ptrs%i5d)) /= 5) CALL finish(routine, 'Only 5-dim variable allowed with *ref*')
    IF (ref_pos < 1 .OR. ref_pos > 5) CALL finish(routine, 'Invalid *ref_pos*')

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%i4d => ptrs%i5d(ref_idx,:,:,:,:)
    CASE(2)
      ptrs%i4d => ptrs%i5d(:,ref_idx,:,:,:)
    CASE(3)
      ptrs%i4d => ptrs%i5d(:,:,ref_idx,:,:)
    CASE(4)
      ptrs%i4d => ptrs%i5d(:,:,:,ref_idx,:)
    CASE(5)
      ptrs%i4d => ptrs%i5d(:,:,:,:,ref_idx)
    END SELECT
  END SUBROUTINE slice_int_5d_to_4d_single_ref

  ! Detailed slicing subroutines for integer kind - double ref
  SUBROUTINE slice_int_5d_to_1d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx, ref2_pos, ref2_idx
    INTEGER, POINTER :: ptr_i2d(:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_int_5d_to_1d_double_ref'

    IF (ref_pos < 1 .OR. ref_pos > 2) CALL finish(routine, 'Invalid *ref_pos*')
    IF (ref2_pos < 2 .OR. ref2_pos > 3) CALL finish(routine, 'Invalid *ref2_pos*')
    IF (SIZE(ptrs%i5d,4) /= 1 .OR. SIZE(ptrs%i5d,5) /= 1) THEN
      CALL finish(routine, 'Only 3-dim variable allowed with *ref* and *ref2*')
    END IF

    SELECT CASE (ref2_pos)
    CASE(2)
      ptr_i2d => ptrs%i5d(:,ref2_idx,:,1,1)
    CASE(3)
      ptr_i2d => ptrs%i5d(:,:,ref2_idx,1,1)
    END SELECT

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%i1d => ptr_i2d(ref_idx,:)
    CASE(2)
      ptrs%i1d => ptr_i2d(:,ref_idx)
    END SELECT
  END SUBROUTINE slice_int_5d_to_1d_double_ref

  SUBROUTINE slice_int_5d_to_2d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx, ref2_pos, ref2_idx
    INTEGER, POINTER :: ptr_i3d(:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_int_5d_to_2d_double_ref'

    IF (ref_pos < 1 .OR. ref_pos > 3) CALL finish(routine, 'Invalid *ref_pos*')
    IF (ref2_pos < 2 .OR. ref2_pos > 4) CALL finish(routine, 'Invalid *ref2_pos*')
    IF (SIZE(ptrs%i5d,5) /= 1) THEN
      CALL finish(routine, 'Only 4-dim variable allowed with *ref* and *ref2*')
    END IF

    SELECT CASE (ref2_pos)
    CASE(2)
      ptr_i3d => ptrs%i5d(:,ref2_idx,:,:,1)
    CASE(3)
      ptr_i3d => ptrs%i5d(:,:,ref2_idx,:,1)
    CASE(4)
      ptr_i3d => ptrs%i5d(:,:,:,ref2_idx,1)
    END SELECT

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%i2d => ptr_i3d(ref_idx,:,:)
    CASE(2)
      ptrs%i2d => ptr_i3d(:,ref_idx,:)
    CASE(3)
      ptrs%i2d => ptr_i3d(:,:,ref_idx)
    END SELECT
  END SUBROUTINE slice_int_5d_to_2d_double_ref

  SUBROUTINE slice_int_5d_to_3d_double_ref(ptrs, ref_pos, ref_idx, ref2_pos, ref2_idx)
    TYPE(t_ptr_container), INTENT(inout) :: ptrs
    INTEGER, INTENT(in) :: ref_pos, ref_idx, ref2_pos, ref2_idx
    INTEGER, POINTER :: ptr_i4d(:,:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':slice_int_5d_to_3d_double_ref'

    IF (ref_pos < 1 .OR. ref_pos > 4) CALL finish(routine, 'Invalid *ref_pos*')
    IF (ref2_pos /= 5) CALL finish(routine, 'Invalid *ref2_pos*')
    IF (SIZE(SHAPE(ptrs%i5d)) /= 5) THEN
      CALL finish(routine, 'Only 5-dim variable allowed with *ref* and *ref2*')
    END IF

    SELECT CASE (ref2_pos)
    CASE(5)
      ptr_i4d => ptrs%i5d(:,:,:,:,ref2_idx)
    END SELECT

    SELECT CASE (ref_pos)
    CASE(1)
      ptrs%i3d => ptr_i4d(ref_idx,:,:,:)
    CASE(2)
      ptrs%i3d => ptr_i4d(:,ref_idx,:,:)
    CASE(3)
      ptrs%i3d => ptr_i4d(:,:,ref_idx,:)
    CASE(4)
      ptrs%i3d => ptr_i4d(:,:,:,ref_idx)
    END SELECT
  END SUBROUTINE slice_int_5d_to_3d_double_ref

  !----------------------------------------------------------------------
  ! Get pointer methods for real values
  !----------------------------------------------------------------------

  ! Get pointer to scalar real variable
  FUNCTION t_tmx_var_get_ptr_r0d(this) RESULT(ptr_r0d)

    CLASS(t_tmx_var), INTENT(in) :: this  !< The t_tmx_var object
    REAL(wp), POINTER            :: ptr_r0d  !< Pointer to the scalar real value

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 0)
    ptr_r0d => ptrs%r0d

  END FUNCTION t_tmx_var_get_ptr_r0d

  ! Get pointer to scalar real variable (vp kind)
  FUNCTION t_tmx_var_get_ptr_v0d(this) RESULT(ptr_v0d)

    CLASS(t_tmx_var), INTENT(in) :: this  !< The t_tmx_var object
    REAL(vp), POINTER            :: ptr_v0d  !< Pointer to the scalar real value

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 0)
    ptr_v0d => ptrs%v0d

  END FUNCTION t_tmx_var_get_ptr_v0d

  ! Get pointer to 1D real array
  FUNCTION t_tmx_var_get_ptr_r1d(this, ref, ref2) RESULT(ptr_r1d)

    CLASS(t_tmx_var), INTENT(in) :: this  !< The t_tmx_var object
    INTEGER, OPTIONAL, INTENT(in) :: ref   !< Optional reference index
    INTEGER, OPTIONAL, INTENT(in) :: ref2  !< Optional second reference index
    REAL(wp), POINTER            :: ptr_r1d(:)  !< Pointer to the 1D real array

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 1, ref, ref2)
    ptr_r1d => ptrs%r1d

  END FUNCTION t_tmx_var_get_ptr_r1d

  ! Get pointer to 1D real array (vp kind)
  FUNCTION t_tmx_var_get_ptr_v1d(this, ref, ref2) RESULT(ptr_v1d)

    CLASS(t_tmx_var), INTENT(in) :: this  !< The t_tmx_var object
    INTEGER, OPTIONAL, INTENT(in) :: ref   !< Optional reference index
    INTEGER, OPTIONAL, INTENT(in) :: ref2  !< Optional second reference index
    REAL(vp), POINTER            :: ptr_v1d(:)  !< Pointer to the 1D real array

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 1, ref, ref2)
    ptr_v1d => ptrs%v1d

  END FUNCTION t_tmx_var_get_ptr_v1d

  ! Helper function to get pointer to 2D real array (wp kind)
  FUNCTION t_tmx_var_get_ptr_r2d(this, ref, ref2) RESULT(ptr_r2d)

    CLASS(t_tmx_var), INTENT(in) :: this
    INTEGER, OPTIONAL, INTENT(in) :: ref, ref2
    REAL(wp), POINTER            :: ptr_r2d(:,:)

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 2, ref, ref2)
    ptr_r2d => ptrs%r2d

  END FUNCTION t_tmx_var_get_ptr_r2d

  ! Helper function to get pointer to 2D real array (vp kind)
  FUNCTION t_tmx_var_get_ptr_v2d(this, ref, ref2) RESULT(ptr_v2d)

    CLASS(t_tmx_var), INTENT(in) :: this
    INTEGER, OPTIONAL, INTENT(in) :: ref, ref2
    REAL(vp), POINTER            :: ptr_v2d(:,:)

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 2, ref, ref2)
    ptr_v2d => ptrs%v2d

  END FUNCTION t_tmx_var_get_ptr_v2d

  ! Get pointer to 3D real array
  FUNCTION t_tmx_var_get_ptr_r3d(this, ref, ref2) RESULT(ptr_r3d)

    CLASS(t_tmx_var),  INTENT(in) :: this  !< The t_tmx_var object
    INTEGER, OPTIONAL, INTENT(in) :: ref   !< Optional reference index
    INTEGER, OPTIONAL, INTENT(in) :: ref2  !< Optional second reference index
    REAL(wp), POINTER             :: ptr_r3d(:,:,:)  !< Pointer to the 3D real array

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 3, ref, ref2)
    ptr_r3d => ptrs%r3d

  END FUNCTION t_tmx_var_get_ptr_r3d

  ! Get pointer to 3D real array (vp kind)
  FUNCTION t_tmx_var_get_ptr_v3d(this, ref, ref2) RESULT(ptr_v3d)

    CLASS(t_tmx_var),  INTENT(in) :: this  !< The t_tmx_var object
    INTEGER, OPTIONAL, INTENT(in) :: ref   !< Optional reference index
    INTEGER, OPTIONAL, INTENT(in) :: ref2  !< Optional second reference index
    REAL(vp), POINTER             :: ptr_v3d(:,:,:)  !< Pointer to the 3D real array

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 3, ref, ref2)
    ptr_v3d => ptrs%v3d

  END FUNCTION t_tmx_var_get_ptr_v3d

  ! Get pointer to 4D real array
  FUNCTION t_tmx_var_get_ptr_r4d(this, ref) RESULT(ptr_r4d)

    CLASS(t_tmx_var),  INTENT(in) :: this  !< The t_tmx_var object
    INTEGER, OPTIONAL, INTENT(in) :: ref   !< Optional reference index
    REAL(wp), POINTER             :: ptr_r4d(:,:,:,:)  !< Pointer to the 4D real array

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 4, ref)
    ptr_r4d => ptrs%r4d

  END FUNCTION t_tmx_var_get_ptr_r4d

  ! Get pointer to 4D real array (vp kind)
  FUNCTION t_tmx_var_get_ptr_v4d(this, ref) RESULT(ptr_v4d)

    CLASS(t_tmx_var),  INTENT(in) :: this  !< The t_tmx_var object
    INTEGER, OPTIONAL, INTENT(in) :: ref   !< Optional reference index
    REAL(vp), POINTER             :: ptr_v4d(:,:,:,:)  !< Pointer to the 4D real array

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 4, ref)
    ptr_v4d => ptrs%v4d

  END FUNCTION t_tmx_var_get_ptr_v4d

  !----------------------------------------------------------------------
  ! Get pointer methods for integer values
  !----------------------------------------------------------------------

  ! Get pointer to scalar integer variable
  FUNCTION t_tmx_var_get_ptr_i0d(this) RESULT(ptr_i0d)

    CLASS(t_tmx_var), INTENT(in) :: this  !< The t_tmx_var object
    INTEGER, POINTER             :: ptr_i0d  !< Pointer to the scalar integer value

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 0)
    ptr_i0d => ptrs%i0d

  END FUNCTION t_tmx_var_get_ptr_i0d

  ! Get pointer to 1D integer array
  FUNCTION t_tmx_var_get_ptr_i1d(this, ref, ref2) RESULT(ptr_i1d)

    CLASS(t_tmx_var), INTENT(in) :: this  !< The t_tmx_var object
    INTEGER, OPTIONAL, INTENT(in) :: ref   !< Optional reference index
    INTEGER, OPTIONAL, INTENT(in) :: ref2  !< Optional second reference index
    INTEGER, POINTER             :: ptr_i1d(:)  !< Pointer to the 1D integer array

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 1, ref, ref2)
    ptr_i1d => ptrs%i1d

  END FUNCTION t_tmx_var_get_ptr_i1d

  ! Get pointer to 2D integer array
  FUNCTION t_tmx_var_get_ptr_i2d(this, ref, ref2) RESULT(ptr_i2d)

    CLASS(t_tmx_var), INTENT(in) :: this  !< The t_tmx_var object
    INTEGER, OPTIONAL, INTENT(in) :: ref   !< Optional reference index
    INTEGER, OPTIONAL, INTENT(in) :: ref2  !< Optional second reference index
    INTEGER, POINTER             :: ptr_i2d(:,:)  !< Pointer to the 2D integer array

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 2, ref, ref2)
    ptr_i2d => ptrs%i2d

  END FUNCTION t_tmx_var_get_ptr_i2d

  ! Get pointer to 3D integer array
  FUNCTION t_tmx_var_get_ptr_i3d(this, ref, ref2) RESULT(ptr_i3d)

    CLASS(t_tmx_var),  INTENT(in) :: this  !< The t_tmx_var object
    INTEGER, OPTIONAL, INTENT(in) :: ref   !< Optional reference index
    INTEGER, OPTIONAL, INTENT(in) :: ref2  !< Optional second reference index
    INTEGER, POINTER              :: ptr_i3d(:,:,:)  !< Pointer to the 3D integer array

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 3, ref, ref2)
    ptr_i3d => ptrs%i3d

  END FUNCTION t_tmx_var_get_ptr_i3d

  ! Get pointer to 4D integer array
  FUNCTION t_tmx_var_get_ptr_i4d(this, ref) RESULT(ptr_i4d)

    CLASS(t_tmx_var),  INTENT(in) :: this  !< The t_tmx_var object
    INTEGER, OPTIONAL, INTENT(in) :: ref   !< Optional reference index
    INTEGER, POINTER              :: ptr_i4d(:,:,:,:)  !< Pointer to the 4D integer array

    TYPE(t_ptr_container) :: ptrs

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 4, ref)
    ptr_i4d => ptrs%i4d

  END FUNCTION t_tmx_var_get_ptr_i4d

  !----------------------------------------------------------------------
  ! Get pointer methods for logical values
  !----------------------------------------------------------------------

  ! Get pointer to scalar logical variable (from integer representation)
  FUNCTION t_tmx_var_get_ptr_l0d(this) RESULT(ptr_l0d)

    CLASS(t_tmx_var), INTENT(in) :: this  !< The t_tmx_var object
    LOGICAL, POINTER             :: ptr_l0d  !< Pointer to the scalar logical value

    TYPE(t_ptr_container) :: ptrs
    INTEGER(i1), POINTER :: ptr_i0d

    CALL this%Get_base_ptr(ptrs)
    CALL this%Slice_ptr(ptrs, 0)
    ptr_i0d => ptrs%l0d
    IF (ASSOCIATED(ptr_i0d) .AND. ptr_i0d == 1_i1) THEN
      ptr_l0d => true_value
    ELSE
      ptr_l0d => false_value
    END IF

  END FUNCTION t_tmx_var_get_ptr_l0d

END MODULE mo_tmx_var
