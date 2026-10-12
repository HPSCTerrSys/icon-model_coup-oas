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

MODULE mo_tmx_process_class

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: message, finish
  USE mo_fortran_tools, ONLY: init_contiguous_dp, init_contiguous_sp
  USE mo_util_string, ONLY: int2string
  USE mtime,        ONLY: t_datetime => datetime
  USE mo_timer, ONLY: new_timer
  USE mo_surrogate_class, ONLY: t_surrogate
  USE mo_tmx_field_class, ONLY: t_tmx_field, t_tmx_field_p, t_domain !bind_tmx_field
  USE mo_tmx_time_integration_class, ONLY: t_time_scheme
  USE mo_tmx_var, ONLY: t_tmx_var, t_tmx_var_p
  USE memman, ONLY: var_descriptor, add_var_real, allocate_var_dp, get_var_data

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_tmx_process, t_tmx_process_p

  TYPE, ABSTRACT, EXTENDS(t_surrogate) :: t_tmx_process
    ! PRIVATE
    CHARACTER(LEN=:),      ALLOCATABLE :: name         !< Process name
    TYPE(t_domain),        POINTER     :: domain       !< Spatial domain
    REAL(wp)                           :: dt           !< Time step
    LOGICAL                            :: is_initial_time
    TYPE(t_tmx_field_p),   ALLOCATABLE :: states(:)       !< State variables
    TYPE(t_tmx_var_p),     ALLOCATABLE :: tendencies(:)   !< Tendency variables
    TYPE(t_tmx_var_p),     ALLOCATABLE :: new_states(:)   !< New state variables
    INTEGER                            :: max_no_states = 0 !< Maximum number of states
    CLASS(t_time_scheme),  ALLOCATABLE :: time_scheme  !< Time integration scheme
    TYPE(t_tmx_process_p), POINTER     :: processes(:) => NULL() !< Subprocesses
    INTEGER                            :: timer_compute, timer_diagnostics
  CONTAINS
    PROCEDURE                          :: Init_process => Init_tmx_process
    PROCEDURE(init_iface),         DEFERRED :: Init
    PROCEDURE                          :: Add_process
    ! PROCEDURE                          :: Add_state_r2d
    PROCEDURE                          :: Add_state_multi
    PROCEDURE                          :: Add_state_shape_real
    GENERIC                            :: Add_state => Add_state_multi, Add_state_shape_real
    PROCEDURE                          :: Set_time_scheme
    PROCEDURE                          :: Step_forward
    PROCEDURE(compute_iface),      DEFERRED :: Compute !< Compute
    PROCEDURE(compute_iface),      DEFERRED :: Compute_diagnostics !< Compute diagnostics
    PROCEDURE(compute_diag_iface), DEFERRED :: Update_diagnostics !< Update diagnostics
    PROCEDURE                          :: Get_tendency_r2d
    PROCEDURE                          :: Get_tendency_r3d
    PROCEDURE                          :: Get_tendency_r4d
  END TYPE t_tmx_process

  TYPE t_tmx_process_p
    CLASS(t_tmx_process), POINTER :: p
  END TYPE t_tmx_process_p

  ABSTRACT INTERFACE
    SUBROUTINE init_iface(this)
      IMPORT :: t_tmx_process
      CLASS(t_tmx_process), INTENT(inout), TARGET :: this
    END SUBROUTINE
    SUBROUTINE compute_iface(this, datetime)
      IMPORT :: t_tmx_process, t_datetime
      CLASS(t_tmx_process),     INTENT(inout), TARGET :: this
      TYPE(t_datetime), OPTIONAL, INTENT(in),   POINTER :: datetime     !< date and time at beginning of time step
    END SUBROUTINE compute_iface
    SUBROUTINE compute_diag_iface(this)
      IMPORT :: t_tmx_process
      CLASS(t_tmx_process),     INTENT(inout), TARGET :: this
    END SUBROUTINE compute_diag_iface
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_process_class'

CONTAINS

  SUBROUTINE Init_tmx_process(this, dt, name, domain)

    CLASS(t_tmx_process), INTENT(inout)        :: this
    REAL(wp),             INTENT(in)           :: dt
    CHARACTER(len=*),     INTENT(in), OPTIONAL :: name
    TYPE(t_domain),       POINTER,    OPTIONAL :: domain

    this%dt = dt
    this%is_initial_time = .TRUE.

    IF (PRESENT(name)) THEN
      this%name = name
    ELSE
      this%name = 'Unnamed'
    END IF

    IF (PRESENT(domain)) THEN
      this%domain => domain
      __acc_attach(this%domain)
    END IF

    ALLOCATE(this%states(this%max_no_states))
    !$ACC ENTER DATA COPYIN(this%states)
    ALLOCATE(this%new_states(this%max_no_states))
    !$ACC ENTER DATA COPYIN(this%new_states)
    ALLOCATE(this%tendencies(this%max_no_states))
    !$ACC ENTER DATA COPYIN(this%tendencies)

    this%timer_compute     = new_timer('tmx_'//name//'_compute')
    this%timer_diagnostics = new_timer('tmx_'//name//'_diag')

  END SUBROUTINE Init_tmx_process

  SUBROUTINE Add_process(this, process)

    CLASS(t_tmx_process), INTENT(inout) :: this
    CLASS(t_tmx_process), TARGET, INTENT(in) :: process

    TYPE(t_tmx_process_p), POINTER :: processes_(:)
    INTEGER :: n_processes, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_process'

    IF (.NOT. ASSOCIATED(this%processes)) THEN
      ALLOCATE(this%processes(1))
      this%processes(1)%p => process
    ELSE
      n_processes = SIZE(this%processes)
      ALLOCATE(processes_(n_processes+1))
      DO i=1,n_processes
        processes_(i)%p => this%processes(i)%p
      END DO
      processes_(n_processes+1)%p => process
      DEALLOCATE(this%processes)
      this%processes => processes_
    END IF

    CALL message(routine, 'Added sub-process '//process%name//' to '//this%name)

  END SUBROUTINE Add_process

  SUBROUTINE Add_state_shape_real(this, idx, dims, diffusion_type)

    CLASS(t_tmx_process), INTENT(inout), TARGET :: this
    INTEGER,              INTENT(in)            :: idx
    INTEGER,              INTENT(in)            :: dims(:)
    INTEGER,              INTENT(in)            :: diffusion_type

    INTEGER :: ndims, istat

    TYPE(var_descriptor) :: var_desc
    TYPE(t_tmx_field), POINTER :: field
    TYPE(t_tmx_var), POINTER :: var
    REAL(wp), POINTER :: ptr_r3d(:,:,:), ptr_r5d(:,:,:,:,:)
    CHARACTER(LEN=:), ALLOCATABLE :: idx_str

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_state_shape_real'

    ndims = SIZE(dims)

    IF (ndims < 2 .OR. ndims > 3) CALL finish(routine, 'Only 2d or 3d real states supported at this time')

    IF (ndims == 2) THEN
      IF (dims(1) /= this%domain%nproma .OR. dims(2) /= this%domain%nblks_c) THEN
        CALL finish(routine, 'Dimension mismatch')
      END IF
    ELSE IF (ndims == 3) THEN
      IF (this%domain%ntiles > 1) THEN
        IF (dims(1) /= this%domain%nproma .OR. dims(2) /= this%domain%nblks_c .OR. dims(3) /= this%domain%ntiles) THEN
          CALL finish(routine, 'Dimension mismatch')
        END IF
      ELSE IF (this%domain%ntiles == 1) THEN
        IF (dims(1) /= this%domain%nproma .OR. dims(2) /= this%domain%nblks_c .OR. dims(3) /= this%domain%ntiles) THEN
          CALL finish(routine, 'Dimension mismatch')
        END IF
      ELSE IF (this%domain%nlev > 1) THEN
        IF (dims(1) /= this%domain%nproma .OR. dims(2) /= this%domain%nlev .OR. dims(3) /= this%domain%nblks_c) THEN
          CALL finish(routine, 'Dimension mismatch')
        END IF
      ELSE
        CALL finish(routine, 'Dimension mismatch')
      END IF
    ELSE IF (ndims == 4) THEN
      IF (dims(1) /= this%domain%nproma .OR. dims(2) /= this%domain%nlev .OR. &
        & dims(3) /= this%domain%nblks_c .OR. dims(4) /= this%domain%ntiles) THEN
        CALL finish(routine, 'Dimension mismatch')
      END IF
    END IF

    idx_str = TRIM(ADJUSTL(int2string(idx)))

    var_desc = var_descriptor('TMX '//this%name//' state '//idx_str, 1, 1, 1, 1)
    field => t_tmx_field(var_desc%name, "double", dims, diffusion_type, var_desc)
    this%states(idx)%p => field

    ! Add variable to new_state varlist with the same attributes as state
    var_desc = var_descriptor('TMX '//this%name//' new state '//idx_str, &
      &field%patch_id, field%hgrid_id, field%vgrid_id, 5)

    var => t_tmx_var(var_desc%name, "double", var_desc, dims=dims)
    this%new_states(idx)%p => var
    istat = get_var_data(ptr_r5d, var_desc)
    IF (istat /= 0) ERROR STOP
!$OMP PARALLEL
#ifdef __SINGLE_PRECISION
    CALL init_contiguous_sp(ptr_r5d, PRODUCT(SHAPE(ptr_r5d)), 0._wp, lacc=.TRUE.)
#else
    CALL init_contiguous_dp(ptr_r5d, PRODUCT(SHAPE(ptr_r5d)), 0._wp, lacc=.TRUE.)
#endif
!$OMP END PARALLEL
    NULLIFY(var, ptr_r5d)

    ! Add variable to tendencies varlist with the same attributes as state
    var_desc = var_descriptor('TMX '//this%name//' tendency '//idx_str, &
      & field%patch_id, field%hgrid_id, field%vgrid_id, 6)
    var => t_tmx_var(var_desc%name, "double", var_desc, dims=dims)
    this%tendencies(idx)%p => var
    istat = get_var_data(ptr_r5d, var_desc)
    IF (istat /= 0) ERROR STOP
!$OMP PARALLEL
#ifdef __SINGLE_PRECISION
    CALL init_contiguous_sp(ptr_r5d, PRODUCT(SHAPE(ptr_r5d)), 0._wp, lacc=.TRUE.)
#else
    CALL init_contiguous_dp(ptr_r5d, PRODUCT(SHAPE(ptr_r5d)), 0._wp, lacc=.TRUE.)
#endif
!$OMP END PARALLEL
    NULLIFY(var, ptr_r5d)

    CALL message(routine, 'New state: TMX '//idx_str//' for process '//this%name)

  END SUBROUTINE Add_state_shape_real

  SUBROUTINE Add_state_multi(this, idx, diffusion_type, var_desc, rank, ref_pos, ref_idx)

    CLASS(t_tmx_process), INTENT(inout), TARGET :: this
    INTEGER,              INTENT(in)            :: idx
    INTEGER,              INTENT(in)            :: diffusion_type
    TYPE(var_descriptor), INTENT(in)            :: var_desc
    INTEGER,              INTENT(in)            :: rank
    INTEGER, OPTIONAL,    INTENT(in)            :: ref_pos
    INTEGER, OPTIONAL,    INTENT(in)            :: ref_idx(:)

    REAL(wp), POINTER :: ptr_r5d(:,:,:,:,:)
    TYPE(t_tmx_field), POINTER :: field
    TYPE(t_tmx_var), POINTER :: var
    TYPE(var_descriptor) :: tmp_var_desc
    INTEGER :: istat
    INTEGER :: dims(rank)
    INTEGER :: tmp_shape(5)
    CHARACTER(LEN=:), ALLOCATABLE :: idx_str

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_state_multi'

    idx_str = TRIM(ADJUSTL(int2string(idx)))

    IF (rank < 2 .OR. rank > 4) CALL finish(routine, &
      & 'State '//idx_str//' - only 2d, 3d or 4d variable supported')

    istat = get_var_data(ptr_r5d, var_desc)
    IF (istat /= 0) THEN
      CALL finish(routine, var_desc%name//' not found.')
    END IF

    tmp_shape(1:5) = SHAPE(ptr_r5d)
    dims(1:rank) = tmp_shape(1:rank)

    field => t_tmx_field('TMX '//this%name//' state '//idx_str, &
      & "double", dims, diffusion_type, var_desc, ref_pos, ref_idx)
    this%states(idx)%p => field

    SELECT CASE (field%rank)
    CASE (2)
      IF (SIZE(ptr_r5d,1) /= this%domain%nproma .OR. SIZE(ptr_r5d,2) /= this%domain%nblks_c) THEN
        CALL finish(routine, 'Dimension mismatch for '//var_desc%name//' in '//this%name)
      END IF
    CASE (3)
      IF (this%domain%ntiles > 1) THEN
        IF (SIZE(ptr_r5d,1) /= this%domain%nproma .OR. SIZE(ptr_r5d,2) /= this%domain%nblks_c) THEN
          CALL finish(routine, 'Dimension mismatch for '//var_desc%name//' in '//this%name)
        END IF
      ELSE IF (this%domain%nlev > 1) THEN
        IF (SIZE(ptr_r5d,1) /= this%domain%nproma .OR. SIZE(ptr_r5d,3) /= this%domain%nblks_c) THEN
          CALL finish(routine, 'Dimension mismatch for '//var_desc%name//' in '//this%name)
        END IF
      END IF
    CASE (4)
      IF (this%domain%ntiles > 1) THEN
        IF (SIZE(ptr_r5d,1) /= this%domain%nproma .OR. SIZE(ptr_r5d,2) /= this%domain%nblks_c) THEN
          CALL finish(routine, 'Dimension mismatch for '//var_desc%name//' in '//this%name)
        END IF
      ELSE IF (this%domain%nlev > 1) THEN
        IF (SIZE(ptr_r5d,1) /= this%domain%nproma .OR. SIZE(ptr_r5d,3) /= this%domain%nblks_c) THEN
          CALL finish(routine, 'Dimension mismatch for '//var_desc%name//' in '//this%name)
        END IF
      END IF
    END SELECT

    ! Add variable to new_states varlist with the same attributes as state
    tmp_var_desc = var_descriptor('TMX '//this%name//' new state '//idx_str, &
      & field%patch_id, field%hgrid_id, field%vgrid_id, 5)

    var => t_tmx_var(tmp_var_desc%name, "double", tmp_var_desc, dims=dims, ref_pos=ref_pos)
    this%new_states(idx)%p => var
    istat = get_var_data(ptr_r5d, tmp_var_desc)
!$OMP PARALLEL
#ifdef __SINGLE_PRECISION
    CALL init_contiguous_sp(ptr_r5d, PRODUCT(SHAPE(ptr_r5d)), 0._wp, lacc=.TRUE.)
#else
    CALL init_contiguous_dp(ptr_r5d, PRODUCT(SHAPE(ptr_r5d)), 0._wp, lacc=.TRUE.)
#endif
!$OMP END PARALLEL

    ! Add variable to tendencies varlist with the same attributes as state
    tmp_var_desc = var_descriptor('TMX '//this%name//' tendency '//idx_str, &
      & field%patch_id, field%hgrid_id, field%vgrid_id, 6)

    var => t_tmx_var(tmp_var_desc%name, "double", tmp_var_desc, dims=dims, ref_pos=ref_pos)
    this%tendencies(idx)%p => var
    istat = get_var_data(ptr_r5d, tmp_var_desc)
!$OMP PARALLEL
#ifdef __SINGLE_PRECISION
    CALL init_contiguous_sp(ptr_r5d, PRODUCT(SHAPE(ptr_r5d)), 0._wp, lacc=.TRUE.)
#else
    CALL init_contiguous_dp(ptr_r5d, PRODUCT(SHAPE(ptr_r5d)), 0._wp, lacc=.TRUE.)
#endif
!$OMP END PARALLEL

    CALL message(routine, 'New state: TMX '//idx_str//' ('//field%var_descriptor%name//') for process '//this%name)

  END SUBROUTINE Add_state_multi

  FUNCTION Get_tendency_r2d(this, idx) RESULT(result)

    CLASS(t_tmx_process), INTENT(in) :: this
    INTEGER,              INTENT(in) :: idx
    REAL(wp), POINTER                :: result(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_tendency_r2d'

    result => NULL()
    result => this%tendencies(idx)%p%Get_ptr_r2d()
    IF (.NOT. ASSOCIATED(result)) CALL finish(routine, 'Could not fetch tendency for state '// &
      & TRIM(ADJUSTL(int2string(idx)))//' of process '//this%name)

  END FUNCTION Get_tendency_r2d

  FUNCTION Get_tendency_r3d(this, idx) RESULT(result)

    CLASS(t_tmx_process), INTENT(in) :: this
    INTEGER,              INTENT(in) :: idx
    REAL(wp), POINTER                :: result(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_tendency_r3d'

    result => NULL()
    result => this%tendencies(idx)%p%Get_ptr_r3d()
    IF (.NOT. ASSOCIATED(result)) CALL finish(routine, 'Could not fetch tendency for '// &
      & TRIM(ADJUSTL(int2string(idx)))//' of process '//this%name)

  END FUNCTION Get_tendency_r3d

  FUNCTION Get_tendency_r4d(this, idx) RESULT(result)

    CLASS(t_tmx_process), INTENT(in) :: this
    INTEGER,              INTENT(in) :: idx
    REAL(wp), POINTER                :: result(:,:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_tendency_r4d'

    result => NULL()
    result => this%tendencies(idx)%p%Get_ptr_r4d()
    IF (.NOT. ASSOCIATED(result)) CALL finish(routine, 'Could not fetch tendency for '// &
      & TRIM(ADJUSTL(int2string(idx)))//' of process '//this%name)

  END FUNCTION Get_tendency_r4d

  SUBROUTINE Set_time_scheme(this, time_scheme)

    CLASS(t_tmx_process), INTENT(inout) :: this
    CLASS(t_time_scheme), INTENT(in)    :: time_scheme

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_time_scheme'

    IF (ALLOCATED(this%time_scheme)) DEALLOCATE(this%time_scheme)

    ALLOCATE(this%time_scheme, source=time_scheme)

  END SUBROUTINE Set_time_scheme

  SUBROUTINE Step_forward(this)

    CLASS(t_tmx_process) :: this
    ! REAL(wp), INTENT(in) :: dt

    CHARACTER(len=*), PARAMETER :: routine = modname//':Step_forward'

    IF (ALLOCATED(this%time_scheme)) THEN
      CALL this%time_scheme%Step_forward(this, this%dt)
    ELSE
      CALL finish(routine, 'No time scheme defined')
    END IF

    END SUBROUTINE Step_forward

  END MODULE mo_tmx_process_class
