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

! Write and read RBF coefficients

MODULE mo_rbf_coefficients_io

  USE mo_kind, ONLY: wp, sp
  USE mo_intp_data_strc, ONLY: t_int_state
  USE mo_exception, ONLY: message, message_text, warning, finish
  USE mo_mpi, ONLY: my_process_is_mpi_workroot, process_mpi_all_workroot_id, &
      &                               work_mpi_barrier, p_bcast, p_comm_work, &
      &                               p_comm_work_test
  USE mo_interpol_config, ONLY: rbf_vec_dim_c, rbf_c2grad_dim, rbf_vec_dim_v, &
      &                               rbf_vec_dim_e, rbf_vec_scale_c, rbf_vec_scale_e, &
      &                               rbf_vec_scale_v, rbf_coeffs_filename
  USE mo_parallel_config, ONLY: nproma, p_test_run
  USE mo_model_domain, ONLY: t_patch
  USE mo_netcdf_errhandler, ONLY: nf
  USE mo_netcdf, ONLY: NF90_CLOBBER, NF90_GLOBAL, NF90_NETCDF4, &
      & nf90_create, nf90_put_att, nf90_enddef, nf90_close, nf90_get_att, nf90_def_dim, &
      & nf90_def_var, nf90_put_var, nf90_noerr
#ifdef __SINGLE_PRECISION
  USE mo_netcdf, ONLY: NF90_FLOAT
#else
  USE mo_netcdf, ONLY: NF90_DOUBLE
#endif
  USE mo_communication, ONLY: exchange_data, t_comm_gather_pattern
  USE mo_read_interface, ONLY: openInputFile, closeFile, on_cells, on_edges, &
      &                               on_vertices, t_stream_id, read_3D
  USE mo_sync, ONLY: sync_patch_array, SYNC_C, SYNC_E, SYNC_V
  USE mo_read_netcdf_types, ONLY: t_alloc_3d
  USE mo_util_uuid_types, ONLY: t_uuid, UUID_STRING_LENGTH
  USE mo_util_uuid, ONLY: uuid_parse, uuid_unparse, OPERATOR(==)
  USE mo_io_units, ONLY: filename_max
  USE mo_util_floating_point, ONLY: is_fp_equivalent

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rbf_coefficients_write, rbf_coefficients_read
  PUBLIC :: rbf_coefficients_round

  ! pack src array into allocated dst array {nproma,nlev,nblks}
  INTERFACE allocate_and_pack_into_nlev
    MODULE PROCEDURE allocate_and_pack_into_nlev_3d_wp
    MODULE PROCEDURE allocate_and_pack_into_nlev_4d_wp
  END INTERFACE allocate_and_pack_into_nlev


  ! unpack src array {nproma,nlev,nblks} into allocated dst array
  INTERFACE unpack_from_nlev
    MODULE PROCEDURE unpack_from_nlev_3d_wp
    MODULE PROCEDURE unpack_from_nlev_4d_wp
  END INTERFACE unpack_from_nlev

  ! Module variables
  INTEGER, PARAMETER :: MAX_LEN_NAMES = 30
  INTEGER, PARAMETER :: MAX_NDIMS = 4
  CHARACTER(*), PARAMETER :: modname = "mo_rbf_coefficients_io"
#ifdef __SINGLE_PRECISION
  INTEGER, PARAMETER :: NF90_REAL_WP = NF90_FLOAT
#else
  INTEGER, PARAMETER :: NF90_REAL_WP = NF90_DOUBLE
#endif

  ! Used by rbf_coefficients_write() to store info metadata about each output variable
  TYPE t_rbf_netcdf_var
    CHARACTER(len=MAX_LEN_NAMES) :: varname
    INTEGER :: varid !< netcdf var id
    INTEGER :: vartype !< netcdf var {NF90_INT, NF90_REAL_WP}
    INTEGER :: ndims !< number of deblocked dimensions
    REAL(wp), POINTER :: ptr_wp(:, :, :) => NULL()
    INTEGER :: dim_indices(MAX_NDIMS) = -1 !< Indices for `dims` belonging to this variable in all_dims(:)
    !< used to populate dims and dimids from all_dims and all_dimids
    INTEGER :: dims(MAX_NDIMS) = -1 !< Dimensions for this 4D array
    INTEGER :: dimids(MAX_NDIMS) = -1 !< netcdf dimids for this 4D array
    INTEGER :: pat_type ! SYNC_C, SYNC_E, SYNC_V used to identify which t_comm_gather_pattern
  END TYPE t_rbf_netcdf_var

CONTAINS

  !
  ! Round rbf coefficients to single-precision
  !
  SUBROUTINE rbf_coefficients_round(ptr_int_state)

    TYPE(t_int_state), INTENT(INOUT) :: ptr_int_state

    CALL warning('RBF:', 'Rounding RBF coefficients!')

    ptr_int_state%rbf_vec_coeff_c = REAL(ptr_int_state%rbf_vec_coeff_c, KIND=sp)
    ptr_int_state%rbf_vec_coeff_e = REAL(ptr_int_state%rbf_vec_coeff_e, KIND=sp)
    ptr_int_state%rbf_vec_coeff_v = REAL(ptr_int_state%rbf_vec_coeff_v, KIND=sp)
    ptr_int_state%rbf_c2grad_coeff = REAL(ptr_int_state%rbf_vec_coeff_c, KIND=sp)

  END SUBROUTINE rbf_coefficients_round


  !
  ! Write rbf coefficients to single file
  !
  SUBROUTINE rbf_coefficients_write(ptr_int_state, ptr_patch, jg)

    TYPE(t_int_state), INTENT(INOUT) :: ptr_int_state
    TYPE(t_patch), INTENT(IN) :: ptr_patch
    INTEGER, INTENT(IN) :: jg
    ! Local vars
    INTEGER :: ncid ! Only used on root
    INTEGER :: i
    CHARACTER(len=filename_max) :: filename
    LOGICAL :: is_root ! root proc for gather and io
    CHARACTER(*), PARAMETER :: routine = modname//":rbf_coefficients_write"

    ! Describing each dimension used by rbf coeff
    INTEGER, PARAMETER :: ndimids = 13
    INTEGER :: dimids_def(ndimids)
    INTEGER :: dims_def(ndimids)
    CHARACTER(len=MAX_LEN_NAMES) :: dimnames_def(ndimids)

    ! Describing each rbf array
    INTEGER, PARAMETER :: nvars = 4 ! 4:coeffs-only, 15:coeffs+indices
    TYPE(t_rbf_netcdf_var) :: rbf_vars(nvars)
    TYPE(t_alloc_3d), TARGET :: buf_wp(nvars)
    CHARACTER(LEN=UUID_STRING_LENGTH) :: uuid_grid_string

    CALL message(routine, 'Writing RBF coefficients')

    ! gather-scatter patterns use work0 as root
    is_root = my_process_is_mpi_workroot()

    ! i = 1, ndimids
    i = 1
    dims_def(i) = 2
    dimnames_def(i) = 'two'
    i = 2
    dims_def(i) = rbf_vec_dim_c
    dimnames_def(i) = 'rbf_vec_dim_c'
    i = 3
    dims_def(i) = rbf_vec_dim_e
    dimnames_def(i) = 'rbf_vec_dim_e'
    i = 4
    dims_def(i) = rbf_vec_dim_v
    dimnames_def(i) = 'rbf_vec_dim_v'
    i = 5
    dims_def(i) = rbf_c2grad_dim
    dimnames_def(i) = 'rbf_c2grad_dim'
    i = 6
    dims_def(i) = ptr_patch%n_patch_cells_g
    dimnames_def(i) = 'ncells'
    i = 7
    dims_def(i) = ptr_patch%n_patch_edges_g
    dimnames_def(i) = 'nedges'
    i = 8
    dims_def(i) = ptr_patch%n_patch_verts_g
    dimnames_def(i) = 'nverts'
    i = 9
    dims_def(i) = 2 * rbf_vec_dim_c
    dimnames_def(i) = 'two_x_rbf_vec_dim_c'
    i = 10
    dims_def(i) = 2 * rbf_vec_dim_e
    dimnames_def(i) = 'two_x_rbf_vec_dim_e'
    i = 11
    dims_def(i) = 2 * rbf_vec_dim_v
    dimnames_def(i) = 'two_x_rbf_vec_dim_v'
    i = 12
    dims_def(i) = 2 * rbf_c2grad_dim
    dimnames_def(i) = 'two_x_rbf_c2grad_dim'
    i = 13
    dims_def(i) = 1
    dimnames_def(i) = 'one'

    ! Floating-point rbf coefficients only
    CALL allocate_and_pack_into_nlev(ptr_int_state%rbf_vec_coeff_c, buf_wp(1)%a, ptr_patch%nblks_c)
    CALL init_rbf_netcdf_var(rbf_vars(1), 'rbf_vec_coeff_c', dims_def, ndimids, (/6, 9, -1, -1/), 2, buf_wp(1)%a, SYNC_C)
    CALL allocate_and_pack_into_nlev(ptr_int_state%rbf_c2grad_coeff, buf_wp(2)%a, ptr_patch%nblks_c)
    CALL init_rbf_netcdf_var(rbf_vars(2), 'rbf_c2grad_coeff', dims_def, ndimids, (/6, 12, -1, -1/), 2, buf_wp(2)%a, SYNC_C)
    CALL allocate_and_pack_into_nlev(ptr_int_state%rbf_vec_coeff_v, buf_wp(3)%a, ptr_patch%nblks_v)
    CALL init_rbf_netcdf_var(rbf_vars(3), 'rbf_vec_coeff_v', dims_def, ndimids, (/8, 11, -1, -1/), 2, buf_wp(3)%a, SYNC_V)
    CALL allocate_and_pack_into_nlev(ptr_int_state%rbf_vec_coeff_e, buf_wp(4)%a, ptr_patch%nblks_e)
    CALL init_rbf_netcdf_var(rbf_vars(4), 'rbf_vec_coeff_e', dims_def, ndimids, (/7, 3, -1, -1/), 2, buf_wp(4)%a, SYNC_E)

    ! Create all {dimid,varids} and put in rbf_netcdf_var
    IF (is_root) THEN

      CALL get_filename(jg, filename)
      WRITE(message_text, '(A,A,A)') "filename '", TRIM(filename), "'"
      CALL message(routine, message_text)

      ! Create file
      CALL nf(nf90_create(filename, IOR(NF90_CLOBBER, NF90_NETCDF4), ncid), filename)

      ! Create dimids and put in rbf_netcdf_var
      CALL create_nc_dimids(ncid, dims_def, dimnames_def, dimids_def, ndimids)
      CALL put_dimids_rbf_netcdf_var(rbf_vars, nvars, dimids_def, ndimids)

      ! Create varids and put in rbf_netcdf_var
      CALL create_and_put_nc_varids(ncid, rbf_vars, nvars)

      ! Attributes
      CALL uuid_unparse(ptr_patch%grid_uuid, uuid_grid_string)
      CALL nf(nf90_put_att(ncid, NF90_GLOBAL, 'uuidOfHGrid', uuid_grid_string), "")
      CALL nf(nf90_put_att(ncid, NF90_GLOBAL, 'rbf_vec_scale_c', rbf_vec_scale_c(jg)), "")
      CALL nf(nf90_put_att(ncid, NF90_GLOBAL, 'rbf_vec_scale_e', rbf_vec_scale_e(jg)), "")
      CALL nf(nf90_put_att(ncid, NF90_GLOBAL, 'rbf_vec_scale_v', rbf_vec_scale_v(jg)), "")
      ! Fallback attributes if uuidOfHGrid is NULL
      CALL nf(nf90_put_att(ncid, NF90_GLOBAL, 'ncells', ptr_patch%n_patch_cells_g), "")
      CALL nf(nf90_put_att(ncid, NF90_GLOBAL, 'nedges', ptr_patch%n_patch_edges_g), "")
      CALL nf(nf90_put_att(ncid, NF90_GLOBAL, 'nverts', ptr_patch%n_patch_verts_g), "")

      CALL nf(nf90_enddef(ncid), "")
    END IF

    DO i = 1, nvars
      CALL gather_and_write_rbf_netcdf_var(ncid, ptr_patch, rbf_vars(i), is_root)

      IF (ALLOCATED(buf_wp(i)%a)) DEALLOCATE(buf_wp(i)%a)
    END DO

    ! Close file
    IF (is_root) THEN
      CALL nf(nf90_close(ncid), "Closing file: "//filename)
    END IF

    CALL work_mpi_barrier
    CALL message(routine, 'Writing RBF coefficients complete')

  END SUBROUTINE rbf_coefficients_write


  !
  ! Read rbf coefficients from single file
  !
  SUBROUTINE rbf_coefficients_read(ptr_int_state, ptr_patch, jg, rbf_read_status)

    TYPE(t_int_state), INTENT(INOUT) :: ptr_int_state
    TYPE(t_patch), INTENT(IN) :: ptr_patch
    INTEGER, INTENT(IN) :: jg
    INTEGER, INTENT(OUT) :: rbf_read_status
    ! Local vars
    CHARACTER(*), PARAMETER :: routine = modname//":rbf_coefficients_read"
    TYPE(t_stream_id) :: stream_id !< file stream_id on workroot proc
    CHARACTER(len=filename_max) :: filename
    REAL(wp), ALLOCATABLE :: buf_wp(:, :, :)
    REAL(wp) :: attrib_wp
    INTEGER :: attrib_int
    CHARACTER(len=UUID_STRING_LENGTH) :: attrib_str
    TYPE(t_uuid) :: attrib_grid_uuid ! unparsed from patch%grid_uuid
    INTEGER :: p_comm

    CALL message(routine, 'Reading RBF coefficients')

    CALL get_filename(jg, filename)
    WRITE(message_text, '(A,A)') 'Reading file ', TRIM(filename)
    CALL message('mo_rbf_coefficients', message_text)

    CALL openInputFile(stream_id, filename, ptr_patch)

    ! Verify attributes
    rbf_read_status = 0
    IF (my_process_is_mpi_workroot()) THEN
      CALL nf(nf90_get_att(stream_id%file_id, NF90_GLOBAL, 'rbf_vec_scale_c', attrib_wp), "")
      IF (.NOT.is_fp_equivalent(attrib_wp, rbf_vec_scale_c(jg))) THEN
        rbf_read_status = -1
        CALL warning(routine, "rbf_vec_scale_c does not match input file")
      END IF

      CALL nf(nf90_get_att(stream_id%file_id, NF90_GLOBAL, 'rbf_vec_scale_e', attrib_wp), "")
      IF (.NOT.is_fp_equivalent(attrib_wp, rbf_vec_scale_e(jg))) THEN
        rbf_read_status = -1
        CALL warning(routine, "rbf_vec_scale_e does not match input file")
      END IF

      CALL nf(nf90_get_att(stream_id%file_id, NF90_GLOBAL, 'rbf_vec_scale_v', attrib_wp), "")
      IF (.NOT.is_fp_equivalent(attrib_wp, rbf_vec_scale_v(jg))) THEN
        rbf_read_status = -1
        CALL warning(routine, "rbf_vec_scale_v does not match input file")
      END IF

      IF (nf90_get_att(stream_id%file_id, NF90_GLOBAL, 'uuidOfHGrid', attrib_str) == nf90_noerr) THEN
        ! uuid provided, check if it matches ptr_patch%grid_uuid
        CALL uuid_parse(attrib_str, attrib_grid_uuid)
        IF (.NOT.(ptr_patch%grid_uuid == attrib_grid_uuid)) THEN
          rbf_read_status = -1
          CALL warning(routine, "uuidOfHGrid does not match input file")
        END IF
      ELSE
        ! uuid not provided, check gridsize as fallback
        CALL warning(routine, "uuidOfHGrid not set as an attribute!")

        CALL nf(nf90_get_att(stream_id%file_id, NF90_GLOBAL, 'ncells', attrib_int), "")
        IF (attrib_int /= ptr_patch%n_patch_cells_g) THEN
          rbf_read_status = -1
          CALL warning(routine, "ncells does not match input file")
        END IF

        CALL nf(nf90_get_att(stream_id%file_id, NF90_GLOBAL, 'nedges', attrib_int), "")
        IF (attrib_int /= ptr_patch%n_patch_edges_g) THEN
          rbf_read_status = -1
          CALL warning(routine, "nedges does not match input file")
        END IF

        CALL nf(nf90_get_att(stream_id%file_id, NF90_GLOBAL, 'nverts', attrib_int), "")
        IF (attrib_int /= ptr_patch%n_patch_verts_g) THEN
          rbf_read_status = -1
          CALL warning(routine, "nverts does not match input file")
        END IF
      END IF
    END IF

    ! Broadcast rbf_read_status
    IF (p_test_run) THEN
      p_comm = p_comm_work_test
    ELSE
      p_comm = p_comm_work
    END IF
    CALL p_bcast(rbf_read_status, process_mpi_all_workroot_id, comm=p_comm)
    IF (rbf_read_status /= 0) RETURN

    ! Floating-point rbf coefficients only
    CALL read_3D(stream_id, on_cells, 'rbf_vec_coeff_c', alloc_array=buf_wp)
    CALL sync_patch_array(SYNC_C, ptr_patch, buf_wp, .FALSE.)
    CALL unpack_from_nlev(buf_wp, ptr_int_state%rbf_vec_coeff_c, ptr_patch%nblks_c)
    DEALLOCATE(buf_wp)

    CALL read_3D(stream_id, on_cells, 'rbf_c2grad_coeff', alloc_array=buf_wp)
    CALL sync_patch_array(SYNC_C, ptr_patch, buf_wp, .FALSE.)
    CALL unpack_from_nlev(buf_wp, ptr_int_state%rbf_c2grad_coeff, ptr_patch%nblks_c)
    DEALLOCATE(buf_wp)

    CALL read_3D(stream_id, on_edges, 'rbf_vec_coeff_e', alloc_array=buf_wp)
    CALL sync_patch_array(SYNC_E, ptr_patch, buf_wp, .FALSE.)
    CALL unpack_from_nlev(buf_wp, ptr_int_state%rbf_vec_coeff_e, ptr_patch%nblks_e)
    DEALLOCATE(buf_wp)

    CALL read_3D(stream_id, on_vertices, 'rbf_vec_coeff_v', alloc_array=buf_wp)
    CALL sync_patch_array(SYNC_V, ptr_patch, buf_wp, .FALSE.)
    CALL unpack_from_nlev(buf_wp, ptr_int_state%rbf_vec_coeff_v, ptr_patch%nblks_v)
    DEALLOCATE(buf_wp)

    CALL closeFile(stream_id)

    CALL message(routine, 'Reading RBF coefficients complete')

  END SUBROUTINE rbf_coefficients_read


  !
  ! Unpack routines for reading coefficients as 3d fields, and refitting to rbf arrays (used by rbf_coefficients_read())
  !
  ! Unpack from {nproma,nlev,nblks}->{dim1,nproma,nblks}
  SUBROUTINE unpack_from_nlev_3d_wp(src, dst, nblks)
    REAL(wp), INTENT(IN) :: src(:, :, :)
    REAL(wp), INTENT(OUT) :: dst(:, :, :)
    INTEGER, INTENT(IN) :: nblks
    ! Local vars
    CHARACTER(*), PARAMETER :: routine = modname//":unpack_from_nlev_3d_wp"
    INTEGER :: i, j, k

    IF (SIZE(dst, 2) /= nproma) CALL finish(routine, "dim(ndims-1)/=nproma")
    IF (SIZE(dst, 3) /= nblks) CALL finish(routine, "dim(ndims)/=nblks")

    DO i = 1, SIZE(dst, 1)
      DO j = 1, SIZE(dst, 2)
        DO k = 1, SIZE(dst, 3)
          dst(i, j, k) = src(j, i, k)
        END DO
      END DO
    END DO
  END SUBROUTINE unpack_from_nlev_3d_wp


  ! Unpack from {nproma,nlev,nblks}->{dim1,dim2,nproma,nblks}
  SUBROUTINE unpack_from_nlev_4d_wp(src, dst, nblks)
    REAL(wp), INTENT(IN) :: src(:, :, :)
    REAL(wp), INTENT(OUT) :: dst(:, :, :, :)
    INTEGER, INTENT(IN) :: nblks
    ! Local vars
    CHARACTER(*), PARAMETER :: routine = modname//":unpack_from_nlev_4d_wp"
    INTEGER :: i, j

    IF (SIZE(dst, 3) /= nproma) CALL finish(routine, "dim(ndims-1)/=nproma")
    IF (SIZE(dst, 4) /= nblks) CALL finish(routine, "dim(ndims)/=nblks")
    DO j = 1, SIZE(dst, 2)
      DO i = 1, SIZE(dst, 1)
        dst(i, j, :, :) = src(:, (j - 1) * SIZE(dst, 1) + i, :)
      END DO
    END DO
  END SUBROUTINE unpack_from_nlev_4d_wp


  !
  ! Pack routines for rbf arrays into 3d fields (dynamically allocated) - (used by rbf_coefficients_write())
  !
  ! Packing {dim1,nproma,nblks}->{nproma,dim1,nblks}
  SUBROUTINE allocate_and_pack_into_nlev_3d_wp(src, alloc_array, nblks)
    REAL(wp), INTENT(IN) :: src(:, :, :)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: alloc_array(:, :, :)
    INTEGER, INTENT(IN) :: nblks
    ! Local vars
    CHARACTER(*), PARAMETER :: routine = modname//":allocate_and_pack_into_nlev_3d_wp"
    INTEGER :: i, j, k

    IF (SIZE(src, 2) /= nproma) CALL finish(routine, "dim(ndims-1)/=nproma")
    IF (SIZE(src, 3) /= nblks) CALL finish(routine, "dim(ndims)/=nblks")
    ALLOCATE(alloc_array(nproma, SIZE(src, 1), nblks))

    DO i = 1, SIZE(src, 1)
      DO j = 1, SIZE(src, 2)
        DO k = 1, SIZE(src, 3)
          alloc_array(j, i, k) = src(i, j, k)
        END DO
      END DO
    END DO
  END SUBROUTINE allocate_and_pack_into_nlev_3d_wp


  ! Packing {dim1,dim2,nproma,nblks}->{nproma,dim1*dim2,nblks}
  SUBROUTINE allocate_and_pack_into_nlev_4d_wp(src, alloc_array, nblks)
    REAL(wp), INTENT(IN) :: src(:, :, :, :)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: alloc_array(:, :, :)
    INTEGER, INTENT(IN) :: nblks
    ! Local vars
    CHARACTER(*), PARAMETER :: routine = modname//":allocate_and_pack_into_nlev_4d_wp"
    INTEGER :: i, j, k, l

    IF (SIZE(src, 3) /= nproma) CALL finish(routine, "dim(ndims-1)/=nproma")
    IF (SIZE(src, 4) /= nblks) CALL finish(routine, "dim(ndims)/=nblks")
    ALLOCATE(alloc_array(nproma, SIZE(src, 1) * SIZE(src, 2), nblks))

    DO j = 1, SIZE(src, 2)
      DO i = 1, SIZE(src, 1)
        DO k = 1, SIZE(src, 3)
          DO l = 1, SIZE(src, 4)
            alloc_array(k, (j - 1) * SIZE(src, 1) + i, l) = src(i, j, k, l)
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE allocate_and_pack_into_nlev_4d_wp


  !
  ! Util routines for rbf_coefficients_write()
  !
  SUBROUTINE create_and_put_nc_varids(ncid, rbf_netcdf_vars, nvars)
    INTEGER, INTENT(IN) :: ncid
    INTEGER, INTENT(IN) :: nvars
    TYPE(t_rbf_netcdf_var), TARGET, INTENT(INOUT) :: rbf_netcdf_vars(nvars)
    ! Local vars
    INTEGER :: ivar
    TYPE(t_rbf_netcdf_var), POINTER :: var

    ! Define variable
    DO ivar = 1, nvars
      var => rbf_netcdf_vars(ivar)
      CALL nf(nf90_def_var(ncid, TRIM(var%varname), var%vartype, var%dimids(1:var%ndims), var%varid), "")
    END DO

  END SUBROUTINE create_and_put_nc_varids


  SUBROUTINE create_nc_dimids(ncid, dims, dimnames, dimids, ndimids)
    INTEGER, INTENT(IN) :: ncid
    INTEGER, INTENT(IN) :: ndimids !< global number of defined dimensions
    INTEGER, INTENT(IN) :: dims(ndimids)
    CHARACTER(len=MAX_LEN_NAMES), INTENT(IN) :: dimnames(ndimids)
    INTEGER, INTENT(OUT) :: dimids(ndimids)
    ! Local vars
    INTEGER :: i

    ! Define dimensions
    DO i = 1, ndimids
      CALL nf(nf90_def_dim(ncid, dimnames(i), dims(i), dimids(i)), "")
    END DO

  END SUBROUTINE create_nc_dimids


  SUBROUTINE put_dimids_rbf_netcdf_var(rbf_netcdf_vars, nvars, dimids_g, ndimids)
    INTEGER, INTENT(IN) :: nvars
    INTEGER, INTENT(IN) :: ndimids !< global number of defined dimensions
    TYPE(t_rbf_netcdf_var), TARGET, INTENT(INOUT) :: rbf_netcdf_vars(nvars)
    INTEGER, INTENT(IN) :: dimids_g(ndimids)
    ! Local vars
    INTEGER :: ivar, idim
    TYPE(t_rbf_netcdf_var), POINTER :: var

    ! Put dimids for var from global dimids_g list
    DO ivar = 1, nvars
      var => rbf_netcdf_vars(ivar)

      DO idim = 1, var%ndims
        var%dimids(idim) = dimids_g(var%dim_indices(idim))
      END DO
    END DO

  END SUBROUTINE put_dimids_rbf_netcdf_var


  SUBROUTINE init_rbf_netcdf_var(rbf_netcdf_var, varname, dims_def, ndims_def, dim_indices, ndims, &
      &                            ptr_wp, pat_type)

    TYPE(t_rbf_netcdf_var), INTENT(OUT) :: rbf_netcdf_var

    CHARACTER(len=*), INTENT(IN) :: varname
    INTEGER, INTENT(IN) :: ndims_def
    REAL(wp), TARGET, INTENT(IN) :: ptr_wp(:, :, :)
    INTEGER, INTENT(IN) :: dims_def(ndims_def)
    INTEGER, INTENT(IN) :: dim_indices(MAX_NDIMS)
    INTEGER, INTENT(IN) :: ndims
    INTEGER, INTENT(IN) :: pat_type
    ! Local var
    INTEGER :: idim
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":init_rbf_netcdf_var"

    ! Integer or real variable
    rbf_netcdf_var%vartype = NF90_REAL_WP
    rbf_netcdf_var%ptr_wp => ptr_wp

    rbf_netcdf_var%varname = TRIM(varname)
    rbf_netcdf_var%ndims = ndims
    rbf_netcdf_var%dim_indices(1:ndims) = dim_indices(1:ndims)
    DO idim = 1, ndims
      rbf_netcdf_var%dims(idim) = dims_def(dim_indices(idim))
    END DO
    rbf_netcdf_var%pat_type = pat_type
  END SUBROUTINE init_rbf_netcdf_var


  SUBROUTINE gather_and_write_rbf_netcdf_var(ncid, ptr_patch, var, is_root)
    INTEGER, INTENT(IN) :: ncid
    TYPE(t_patch), TARGET, INTENT(IN) :: ptr_patch
    TYPE(t_rbf_netcdf_var), INTENT(IN) :: var
    LOGICAL, INTENT(IN) :: is_root
    ! Local vars
    CHARACTER(*), PARAMETER :: routine = modname//":gather_and_write_rbf_netcdf_var"
    REAL(wp), ALLOCATABLE :: out_buf_wp(:, :) ! dims [ncells/edges/verts, nlev]
    INTEGER, PARAMETER :: fill_value = -999
    TYPE(t_comm_gather_pattern), POINTER :: gather_pattern

    IF (var%ndims /= 2) CALL finish(routine, "only implimented for (nproma,nlev,nblks)")

    ! Allocate buffer on io proc
    IF (is_root) THEN
      ALLOCATE(out_buf_wp(var%dims(1), var%dims(2)))
    ELSE
      ALLOCATE(out_buf_wp(0, 0))
    END IF

    ! Choose correct gather_pattern
    SELECT CASE (var%pat_type)
    CASE (SYNC_C)
      gather_pattern => ptr_patch%comm_pat_gather_c
    CASE (SYNC_E)
      gather_pattern => ptr_patch%comm_pat_gather_e
    CASE (SYNC_V)
      gather_pattern => ptr_patch%comm_pat_gather_v
    CASE DEFAULT
      CALL finish(routine, 'Illegal type parameter')
    END SELECT

    ! Cycle each field and do output
    ! 2D deblock gather:
    !INTEGER, INTENT(IN   ) :: in_array(:,:,:)  !! dimension (nproma, nlev, nblk)
    !INTEGER, INTENT(INOUT) :: out_array(:,:)   !! dimension (global length, nlev); only required on root
    CALL exchange_data(in_array=var%ptr_wp(:, :, :), out_array=out_buf_wp(:, :), &
        &                gather_pattern=gather_pattern, fill_value=REAL(fill_value, KIND=wp))

    ! Output
    IF (is_root) THEN
      CALL nf(nf90_put_var(ncid, var%varid, out_buf_wp), "")
    END IF
    DEALLOCATE(out_buf_wp)

  END SUBROUTINE gather_and_write_rbf_netcdf_var


  SUBROUTINE get_filename(jg, filename)
    INTEGER, INTENT(IN) :: jg
    CHARACTER(LEN=filename_max), INTENT(OUT) :: filename
    WRITE(filename, '(A)') rbf_coeffs_filename(jg)
  END SUBROUTINE get_filename

END MODULE mo_rbf_coefficients_io
