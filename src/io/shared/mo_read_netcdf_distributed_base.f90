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

#include "omp_definitions.inc"
! -----------------------------------------------------------------------------
! MODULE: mo_read_netcdf_distributed_base
!
! Purpose
! -------
! Provide a reusable, distributed NetCDF read facility. The module implements:
!  - partitioning of global 1-D data across a sub-group of MPI ranks that do IO
!    (IO ranks),
!  - local read to a compact buffer on IO ranks
!
! Design overview
! ---------------
! * IO ranks: A subset of p_comm_work processes are selected to perform NetCDF
!   reads. Selection is controlled by the parallel_config variables
!   (io_process_stride and io_process_rotate).
! * basic_data: Shared information describing the global array layout for a
!   particular global length (e,g, global number of cells, vertices, or edges).
!   Multiple distrib_read handles can share one
!   basic_data entry (reference counted). This avoids repeated partitioning.
! * Read pattern:
!     1. The designated IO ranks open the NetCDF file (collective in pNetCDF
!        mode, otherwise fallback to serial open on IO ranks).
!     2. On IO ranks: read the chunk assigned to that IO rank into a compact
!        local buffer (one-dimensional slice read whenever possible).
!
! Important runtime knobs and semantics
! ------------------------------------
! * IO-selection:
!     - io_process_stride: when > 0, directly specifies spacing between IO ranks.
!       Effective number of IO processes = ceil(p_n_work/io_stride) (with rotate).
!     - io_process_rotate: offset applied to the stride to pick the IO ranks.
!     - If io_process_stride <= 0, module heuristics choose a stride by
!       approximating a square-root decomposition of p_n_work and rounding the
!       stride up to the next power of two.
! * Collective behavior:
!     - When compiled with Parallel NetCDF (HAVE_PARALLEL_NETCDF and MPI), the
!       module calls nf90_open_par and sets variable access to NF90_COLLECTIVE.
!       This requires that the MPI communicator passed to this module (p_comm_work)
!       matches the communicator used by the application at startup for IO groups.
!     - If nf90_open_par returns an error, the code falls back to non-parallel
!       nf90_open (on IO ranks) and emits a warning message.
! * Memory footprint:
!     - IO ranks allocate compact read buffers sized to the IO chunk (ish).
!       Non-IO ranks allocate zero-sized io_chunk extents. The global `basic_data`
!       array caches index-mapping structures to avoid repeated allocation.
!
! Performance & tuning guidance
! -----------------------------
! 1. Choose io_process_stride so each IO rank reads a large contiguous block but
!    not so many that IO ranks are overloaded. Typical rule-of-thumb:
!      - Aim for each IO rank to read 1-4 GB concurrently (depends on memory).
!      - When in doubt, increase io_process_stride (fewer IO ranks).
! 2. If underlying filesystem benefits from fewer large readers (Lustre, GPFS),
!    prefer fewer IO ranks doing larger reads (bigger stride).
!
! -----------------------------------------------------------------------------
MODULE mo_read_netcdf_distributed_base
  USE mo_kind, ONLY: dp, sp
  USE mo_exception, ONLY: finish, message, warning
  USE mo_mpi, ONLY: p_n_work, p_pe_work, p_bcast, p_alltoall, p_comm_work, p_max
  USE ppm_extents, ONLY: extent, extent_intersect
  USE mo_decomposition_tools, ONLY: &
    & t_glb2loc_index_lookup, init_glb2loc_index_lookup, &
    & set_inner_glb_index, deallocate_glb2loc_index_lookup, &
    & uniform_partition, partidx_of_elem_uniform_deco
  USE fortran_support, ONLY: t_ptr_2d_dp, t_ptr_2d_int, t_ptr_2d_sp, &
    & t_ptr_3d_dp, t_ptr_3d_int, t_ptr_3d_sp, t_ptr_4d_dp, t_ptr_4d_int, &
    & t_ptr_4d_sp
  USE mo_parallel_config, ONLY: io_process_stride, io_process_rotate
  USE mo_netcdf_errhandler, ONLY: nf
  USE mo_netcdf
#if defined (HAVE_PARALLEL_NETCDF) && !defined (NOMPI)
  USE mpi, ONLY: MPI_INFO_NULL, MPI_UNDEFINED, MPI_Comm_split, MPI_COMM_NULL
#endif

  IMPLICIT NONE
  PRIVATE

  CHARACTER(*), PARAMETER :: modname = 'mo_read_netcdf_distributed_base'

  PUBLIC :: t_basic_distrib_read_data
  PUBLIC :: setup_distrib_read_base
  PUBLIC :: delete_distrib_read_base
  PUBLIC :: distrib_read_get_decomp
  PUBLIC :: distrib_read_get_owner
  PUBLIC :: distrib_read_get_alltoallv_args
  PUBLIC :: distrib_nf_open_base
  PUBLIC :: distrib_nf_close_base
  PUBLIC :: distrib_nf_inq_varexists_base
  PUBLIC :: distrib_nf_inq_attexists_base
  PUBLIC :: distrib_nf_inq_dimlen_base
  PUBLIC :: distrib_inq_var_dims_base
  PUBLIC :: distrib_nf_get_att_base
  PUBLIC :: distrib_read_base
  PUBLIC :: distrib_read_root_base

  TYPE t_basic_distrib_read_data
    INTEGER :: n_g = -1 ! global number of points (-1 if unused)
    TYPE(extent) :: io_chunk ! io decomposition for reading
    INTEGER :: n_ref = 0 ! number of times this data is referenced
  END TYPE t_basic_distrib_read_data

  TYPE(t_basic_distrib_read_data), TARGET, ALLOCATABLE :: basic_data(:)

  ! These module-level variables mirror the original module; they are safe to
  ! keep here because the io decomposition only depends on p_n_work and the
  ! io_process_* runtime knobs.
  INTEGER :: parRootRank = -1
#if defined (HAVE_PARALLEL_NETCDF) && !defined (NOMPI)
  INTEGER :: io_comm = MPI_COMM_NULL
#endif
  LOGICAL :: this_PE_does_IO = .FALSE.
  LOGICAL :: this_PE_does_meta_IO = .FALSE.
  INTEGER :: io_stride = -1, n_io_proc = -1, rotate = -1

  ! -------------------------------------------------------------------
  ! Reads a 1D/2D INT/FLOAT/DOUBLE array from file into the distributed
  ! buffer
  !
  ! Inputs:
  !  - ncid: NetCDF file id on IO ranks, -1 on non-IO ranks
  !          (returned by distrib_nf_open_base)
  !  - vname: variable name
  !  - buffer: 1D/2D INT/FLOAT/DOUBLE buffer whose first dimension is
  !            allocated according to io_chunk%size
  !  - basic_data_index: valid index into basic_data structure
  !                      (returned by setup_distrib_read_base)
  INTERFACE distrib_read_base
    MODULE PROCEDURE distrib_read_base_1D_int
    MODULE PROCEDURE distrib_read_base_2D_int
    MODULE PROCEDURE distrib_read_base_3D_int
    MODULE PROCEDURE distrib_read_base_1D_sp
    MODULE PROCEDURE distrib_read_base_2D_sp
    MODULE PROCEDURE distrib_read_base_3D_sp
    MODULE PROCEDURE distrib_read_base_1D_dp
    MODULE PROCEDURE distrib_read_base_2D_dp
    MODULE PROCEDURE distrib_read_base_3D_dp
  END INTERFACE distrib_read_base

  ! -------------------------------------------------------------------
  ! Reads a 1D INT array from file on the root rank and broadcasts
  ! it to all other ranks
  !
  ! Inputs:
  !  - ncid: NetCDF file id on IO ranks, -1 on non-IO ranks
  !          (returned by distrib_nf_open_base)
  !  - vname: variable name
  !  - buffer: 1D INT buffer (size should match with count argument)
  !  - start, count: defines slice to be read
  !
  ! Notes:
  !  - Usage only recommended for small arrays, whose contents is required by
  !    all ranks in p_comm_work
  INTERFACE distrib_read_root_base
    MODULE PROCEDURE distrib_read_root_base_1D_int
  END INTERFACE distrib_read_root_base

  ! -------------------------------------------------------------------
  ! Reads attribute data from file
  !
  ! Inputs:
  !  - ncid: NetCDF file id on IO ranks, -1 on non-IO ranks
  !          (returned by distrib_nf_open_base)
  !  - vname: variable name ("global" references global attributes)
  !  - name: attribute name
  !  - value: attribute value
  INTERFACE distrib_nf_get_att_base
    MODULE PROCEDURE distrib_nf_get_att_int_base
    MODULE PROCEDURE distrib_nf_get_att_string_base
  END INTERFACE distrib_nf_get_att_base

CONTAINS

  ! -----------------------------------------------
  ! Prepares data structures for distributed reading of
  ! a global 1-D array of length n_g (e.g. global number of cells).
  !
  ! Inputs:
  !  - n_g: length of the global array
  !
  ! Outputs / inout:
  !  - basic_data_index: handle to be used when reading data of size n_g
  !
  ! Side-effects:
  !  - Allocates or reuses an entry in internal basic_data cache (reference counted).
  !  - Determines IO ranks from io_process_stride / io_process_rotate and creates
  !    an MPI communicator for IO ranks if pNetCDF is enabled.
  SUBROUTINE setup_distrib_read_base(n_g, basic_data_index)

    INTEGER, INTENT(IN) :: n_g
    INTEGER, INTENT(OUT) :: basic_data_index

    CHARACTER(*), PARAMETER :: &
      routine = modname//'::setup_distrib_read_base'

    ! initialize module-level io variables lazily
    IF (n_io_proc == -1) CALL init_module_vars(routine)

    ! find matching or allocate basic_data entry
    basic_data_index = get_data_idx()

  CONTAINS

    INTEGER FUNCTION get_data_idx() RESULT(idx)
      INTEGER :: i, basic_data_size
      TYPE(t_basic_distrib_read_data), ALLOCATABLE :: temp_basic_data(:)
      INTEGER, PARAMETER :: basic_data_alloc_inc_size = 16

      basic_data_size = 0
      idx = 0

      ! the basic_data array is already allocated
      IF (ALLOCATED(basic_data)) THEN

        basic_data_size = SIZE(basic_data)

        ! find a matching entry or if there are not matches,
        ! the last empty entry
        DO i = 1, SIZE(basic_data)

          IF (basic_data(i)%n_g == n_g) THEN ! if matching entry -> early exit
            idx = i
            EXIT
          ELSE IF (basic_data(i)%n_g == -1) THEN ! if unused entry
            idx = i
          END IF

        END DO
      END IF

      ! if not matching entry was found
      IF (idx == 0) THEN

        ! resize array and get index of first new element
        CALL MOVE_ALLOC(basic_data, temp_basic_data)
        ALLOCATE(basic_data(basic_data_size + basic_data_alloc_inc_size))
        IF (basic_data_size > 0) &
          basic_data(1:basic_data_size) = temp_basic_data(1:basic_data_size)
        idx = basic_data_size + 1

      END IF

      ! if the current element is empty -> initialise structure
      IF (basic_data(idx)%n_ref == 0) CALL init_data(idx)

      ! increase reference counter
      basic_data(idx)%n_ref = basic_data(idx)%n_ref + 1
    END FUNCTION get_data_idx

    SUBROUTINE init_data(idx)
      INTEGER, INTENT(IN) :: idx

      basic_data(idx)%n_g = n_g
      ! if the process takes part in the reading
      IF (this_PE_does_IO) THEN
        basic_data(idx)%io_chunk = &
          uniform_partition( &
            extent(1, n_g), n_io_proc, p_pe_work/io_stride + 1)
      ELSE
        basic_data(idx)%io_chunk = extent(1, 0)
      END IF
    END SUBROUTINE init_data

  END SUBROUTINE setup_distrib_read_base

  ! -----------------------------------------------------------------
  ! Gets information of this PE of read decomposition associated with
  ! basic_data_index
  ! Return IO chunk assigned to this PE for basic_data(basic_data_index).
  !
  ! Input:
  !  - basic_data_index: handle to be basic_data entry used for reading data
  !                      of size n_g
  !
  ! Returns:
  !  - local IO chunkg
  FUNCTION distrib_read_get_decomp(basic_data_index)
    INTEGER, INTENT(IN) :: basic_data_index
    TYPE(extent) :: distrib_read_get_decomp

    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_get_decomp'

    CALL check_basic_data_index(basic_data_index, routine)

    distrib_read_get_decomp = basic_data(basic_data_index)%io_chunk
  END FUNCTION distrib_read_get_decomp

  ! -----------------------------------------------------------------
  ! Determine the owning I/O process for each global index in `glb_index`
  ! according to the read decomposition associated with `basic_data_index`.
  !
  ! Input:
  !  - basic_data_index: Handle identifying the basic_data entry that
  !                      defines the read decomposition of size n_g.
  !  - glb_index(:):     Global indices (1...n_g) of elements owned by
  !                      the local process in the user decomposition.
  !
  ! Return:
  !  - Rank indices of the I/O processes responsible for each element of
  !    glb_index according to the read decomposition.
  !
  ! Notes:
  !  - The mapping is computed using partidx_of_elem_uniform_deco,
  !    which distributes the global elements evenly across the set
  !    of I/O processes.
  !  - The returned process indices are shifted by rotate
  !    and scaled by io_stride to reflect the configured I/O layout.
  !  - If glb_index is not allocated or empty, a dummy array of length 1
  !    is allocated for owner.
  ! -----------------------------------------------------------------
  SUBROUTINE distrib_read_get_owner(basic_data_index, glb_index, owner)
    INTEGER, INTENT(IN) :: basic_data_index
    INTEGER, ALLOCATABLE, INTENT(in) :: glb_index(:) !< global ids stored on
                                                     !< the local process
    INTEGER, ALLOCATABLE, INTENT(OUT) :: owner(:)

    INTEGER :: n

    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_get_owner'

    IF (ALLOCATED(glb_index)) THEN
      IF (ANY(glb_index > basic_data(basic_data_index)%n_g)) &
        CALL finish(routine, "global index exceeds global size")
      IF (ANY(glb_index < 1)) &
        CALL finish(routine, "global index is smaller than 1")
      n = SIZE(glb_index)
    ELSE
      n = 0
    END IF
    ALLOCATE(owner(n))
    IF (n > 0) THEN
      owner(:) = &
        (partidx_of_elem_uniform_deco( &
           extent(1, basic_data(basic_data_index)%n_g), &
           n_io_proc,  glb_index(:)) - 1) * io_stride + rotate
    END IF

  END SUBROUTINE distrib_read_get_owner

  ! -----------------------------------------------------------------
  ! Determine count and displacement arguments for an MPI_Alltoallv call,
  ! that can be used to do redistribution of data from the read decomposition
  ! (associated with `basic_data_index`) to the user decomposition (described
  !  by `local_chunk).
  ! It is assumed that in the user decomposition, each process only has a
  ! single contiguous chunk of the global index space ([1..n_g]).
  !
  ! Input:
  !  - basic_data_index: Handle identifying the basic_data entry that
  !                      defines the read decomposition of size n_g.
  !  - local_chunk       Description of contiguous chunk of global indices of
  !                      elements owned by the local process in the user
  !                      decomposition.
  !
  ! Output:
  !  - sendcounts: Argument to be used in an MPI_Alltoallv call
  !  - sdispls:    Argument to be used in an MPI_Alltoallv call
  !  - recvcounts: Argument to be used in an MPI_Alltoallv call
  !  - rdispls:    Argument to be used in an MPI_Alltoallv call
  ! -----------------------------------------------------------------
  SUBROUTINE distrib_read_get_alltoallv_args( &
    basic_data_index, local_chunk, sendcounts, sdispls, recvcounts, rdispls)
    INTEGER, INTENT(IN) :: basic_data_index
    TYPE(extent), INTENT(in) :: local_chunk ! Description of local chunk of
                                            ! user decomposition
    INTEGER, ALLOCATABLE, INTENT(OUT) :: sendcounts(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: sdispls(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: recvcounts(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: rdispls(:)

    INTEGER :: io_rank
    INTEGER, ALLOCATABLE :: temp_sdispls(:)
    TYPE(extent) :: temp_local_chunk, intersect_chunk, curr_io_chunk

    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_get_alltoallv_args'

    CALL check_basic_data_index(basic_data_index, routine)

    ! ensure that the provided local chunk is valid
    IF ((local_chunk%first < 1) .OR. &
        (local_chunk%first > basic_data(basic_data_index)%n_g) .OR. &
        (local_chunk%size + local_chunk%first - 1 > &
         basic_data(basic_data_index)%n_g)) &
      CALL finish(routine, "invalid local_chunk")

    ! Allocate output arguments and initialize then with 0
    ALLOCATE(temp_sdispls(p_n_work), &
             sendcounts(p_n_work), sdispls(p_n_work), &
             recvcounts(p_n_work), rdispls(p_n_work), SOURCE = 0)

    ! Compute output arguments

    temp_local_chunk = local_chunk

    ! iterate over pieces of the local user chunk and map them to the IO ranks
    DO WHILE (temp_local_chunk%size > 0)

      ! determin which IO rank owns the first element of temp_local_chunk
      io_rank = &
        (partidx_of_elem_uniform_deco( &
           extent(1, basic_data(basic_data_index)%n_g), &
           n_io_proc,  temp_local_chunk%first) - 1) * io_stride + rotate

      ! reconstruct the current IO rank's chunk
      curr_io_chunk = &
        uniform_partition( &
          extent(1, basic_data(basic_data_index)%n_g), &
          n_io_proc, io_rank / io_stride + 1)

      ! intersection between the remaining part of local chunk and
      ! the current IO chunk
      intersect_chunk = extent_intersect(temp_local_chunk, curr_io_chunk)

      temp_sdispls(io_rank+1) = intersect_chunk%first - curr_io_chunk%first
      recvcounts(io_rank+1)   = intersect_chunk%size
      rdispls   (io_rank+1)   = intersect_chunk%first - local_chunk%first

      ! advance temp_local_chunk past the intersection
      temp_local_chunk%first = temp_local_chunk%first + intersect_chunk%size
      temp_local_chunk%size  = temp_local_chunk%size  - intersect_chunk%size
    END DO

    ! set sendcounts and sdispls by redistributing recvcounts and temp_sdispls
    CALL p_alltoall(recvcounts, sendcounts, p_comm_work)
    CALL p_alltoall(temp_sdispls, sdispls, p_comm_work)

    DEALLOCATE(temp_sdispls)

  END SUBROUTINE distrib_read_get_alltoallv_args

  ! -----------------------------------------------------------------
  ! Deregister and free resources associated with basic_data_index created by
  ! setup_distrib_read_base. Decrements the internal reference counter and
  ! frees index structures when no other user references the same n_g.
  !
  ! Input:
  !  - basic_data_index: handle to be basic_data entry used for reading data
  !                      of size n_g
  SUBROUTINE delete_distrib_read_base(basic_data_index)
    INTEGER, INTENT(IN) :: basic_data_index
    CHARACTER(*), PARAMETER :: &
      routine = modname//'::delete_distrib_read_base'

    CALL check_basic_data_index(basic_data_index, routine)

    basic_data(basic_data_index)%n_ref = basic_data(basic_data_index)%n_ref - 1

    IF (basic_data(basic_data_index)%n_ref .GT. 0) RETURN

    basic_data(basic_data_index)%n_g = -1
  END SUBROUTINE delete_distrib_read_base

  ! --------------------------------
  ! Open a NetCDF file for distributed reading.
  !
  ! Input:
  !  - path: filesystem path to NetCDF file
  !
  ! Returns:
  !  - ncid: NetCDF file id on IO ranks, -1 on non-IO ranks
  !
  ! Behavior:
  !  - If this PE is selected as an IO PE (this_PE_does_IO .TRUE.) and ICON is
  !    compiled with HAVE_PARALLEL_NETCDF, an attempt is made to open the file
  !    using nf90_open_par (parallel NetCDF, collective IO).
  !  - If the parallel open fails but the file exists, falls back to nf90_open
  !    (serial open) on IO ranks and emits a warning message.
  !  - Non-IO ranks return ncid = -1.
  !
  ! Notes:
  !  - Requires that the IO communicator (io_comm) has been configured by
  !    setup_distrib_read prior to calling if parallel NetCDF is used.
  INTEGER FUNCTION distrib_nf_open_base(path) RESULT(ncid)
    CHARACTER(*), INTENT(in) :: path

#if defined (HAVE_PARALLEL_NETCDF) && !defined (NOMPI)
    INTEGER              :: ierr, nvars, i
    INTEGER, ALLOCATABLE :: varids(:)
#endif
    LOGICAL              :: exists

    CHARACTER(*), PARAMETER :: routine = modname//'::distrib_nf_open_base'

    ! initialize module-level io variables lazily
    IF (n_io_proc == -1) CALL init_module_vars(routine)

    CALL message (routine, path)
    ncid = -1

    ! if this is an IO rank
    IF (this_PE_does_IO) THEN

      ! check whether the file existis
      INQUIRE(file=path, exist=exists)
      IF (.NOT. exists) &
        CALL finish(routine, "File "//TRIM(path)//" does not exist.")

#if defined (HAVE_PARALLEL_NETCDF) && !defined (NOMPI)

      ! try opening the file in parallel mode
      ierr = nf90_open_par( &
        path, IOR(nf90_nowrite, nf90_mpiio), io_comm, MPI_INFO_NULL, ncid)

      ! if the file was opened successfully in parallel mode
      IF (ierr == nf90_noerr) THEN
        ! Switch all vars to collective. Hopefully this is sufficient.
        CALL nf(nf90_inquire(ncid, nVariables = nvars), routine)
        ALLOCATE(varids(nvars))
        CALL nf(nf90_inq_varids(ncid, nvars, varids), routine)
        DO i = 1,nvars
          CALL nf( &
            nf90_var_par_access(ncid, varids(i), NF90_COLLECTIVE), routine)
        ENDDO
      ELSE
        CALL warning( &
          routine, 'falling back to serial semantics for opening netcdf file '//TRIM(path))
#endif

        ! open the file in serial mode if compiled without HAVE_PARALLEL_NETCDF,
        ! with NOMPI, or if opening the file in parallel mode failed
        CALL nf(nf90_open(path, nf90_nowrite, ncid), routine)

#if defined (HAVE_PARALLEL_NETCDF) && !defined (NOMPI)
      ENDIF
#endif
    END IF
  END FUNCTION distrib_nf_open_base

  ! ----------------------
  ! Close the NetCDF file descriptor on IO ranks. No-op on non-IO ranks.
  !
  ! Input:
  !  - ncid: NetCDF file id returned by distrib_nf_open_base
  SUBROUTINE distrib_nf_close_base(ncid)
    INTEGER, INTENT(in) :: ncid

    IF (this_PE_does_IO) THEN
      CALL nf(nf90_close(ncid), modname//'::distrib_nf_close_base')
    END IF
  END SUBROUTINE distrib_nf_close_base

  ! ------------------------------------------------
  ! Query whether var `vname` exists in open NetCDF file `ncid`.
  ! Only the process this_PE_does_IO performs the NetCDF inquiry; the
  ! integer error code is broadcast to all ranks and returned as logical.
  !
  ! Input:
  !  - ncid: NetCDF file id returned by distrib_nf_open_base
  !  - vname: Name of the variable to be checked in the file.
  !
  ! Return:
  !  - .TRUE. if the variable is defined in the file, .FALSE. otherwise
  FUNCTION distrib_nf_inq_varexists_base(ncid, vname)
    INTEGER, INTENT(in) :: ncid
    CHARACTER(*), INTENT(in) :: vname
    INTEGER :: err, vid

    LOGICAL :: distrib_nf_inq_varexists_base

    IF (this_PE_does_meta_IO) err = nf90_inq_varid(ncid, vname, vid)
    CALL p_bcast(err, parRootRank, p_comm_work)

    distrib_nf_inq_varexists_base = (err == nf90_noerr)

  END FUNCTION distrib_nf_inq_varexists_base

  ! ------------------------------------------------
  ! Query whether attribute `name` of variable `vname` exists in open
  ! NetCDF file `ncid` (`vname` == "global" references global attributes`).
  ! Only the process this_PE_does_IO performs the NetCDF inquiry; the
  ! integer error code is broadcast to all ranks and returned as logical.
  !
  ! Input:
  !  - ncid: NetCDF file id returned by distrib_nf_open_base
  !  - vname: Name of the variable to be checked in the file.
  !           ("global" reference a global attribute)
  !  - name: Name of the attribute to be checked in the file.
  !
  ! Return:
  !  - .TRUE. if the attribute is defined in the file, .FALSE. otherwise
  FUNCTION distrib_nf_inq_attexists_base(ncid, vname, name)
    INTEGER, INTENT(in) :: ncid
    CHARACTER(*), INTENT(in) :: vname
    CHARACTER(*), INTENT(in) :: name
    INTEGER :: err

    LOGICAL :: distrib_nf_inq_attexists_base

    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_nf_inq_attexists_base'

    IF (this_PE_does_meta_IO) THEN
      err = nf90_inquire_attribute(ncid, get_vid(routine, ncid, vname), name)
    END IF

    CALL p_bcast(err, parRootRank, p_comm_work)

    distrib_nf_inq_attexists_base = (err == nf90_noerr)

  END FUNCTION distrib_nf_inq_attexists_base

  ! -------------------------------------------------------------------
  ! Inquire the length of a NetCDF dimension collectively.
  !
  ! The root process reads the dimension length from the given NetCDF file
  ! and broadcasts the value to all other processes in the communicator.
  !
  ! Inputs:
  !  - ncid: NetCDF file id returned by distrib_nf_open_cont
  !  - dimname: name of the dimension
  !
  ! Outputs:
  !  - dimlen: length of the dimension
  SUBROUTINE distrib_nf_inq_dimlen_base(ncid, dimname, dimlen)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(len=*), INTENT(IN) :: dimname
    INTEGER, INTENT(OUT) :: dimlen

    CHARACTER(*), PARAMETER :: routine = modname//'::distrib_nf_inq_dimlen_base'
    INTEGER :: dimid

    !-----------------------------------------------------------------------
    ! Only the root rank reads the dimension length.
    !-----------------------------------------------------------------------
    IF (this_PE_does_meta_IO) THEN
      CALL nf(nf90_inq_dimid(ncid, dimname, dimid), routine)
      CALL nf(nf90_inquire_dimension(ncid, dimid, len=dimlen), routine)
    END IF

    !-----------------------------------------------------------------------
    ! Broadcast the dimension length to all ranks
    !-----------------------------------------------------------------------
    CALL p_bcast(dimlen, parRootRank, p_comm_work)

  END SUBROUTINE distrib_nf_inq_dimlen_base

  ! ------------------------------------------------
  ! Query the number and lengths of dimensions of variable `vname`
  ! in the open NetCDF file identified by `ncid`.
  ! Only the process with this_PE_does_IO performs the NetCDF
  ! inquiry; the resulting dimension count and lengths are broadcast
  ! to all ranks.
  !
  ! Input:
  !  - ncid: NetCDF file ID returned by distrib_nf_open_base
  !  - vname: Name of the variable whose dimensions are queried
  !
  ! Output:
  !  - var_ndims: Number of dimensions of the variable
  !  - var_dimlen(:): Lengths of the variable's dimensions (must be large
  !                   enough to hold all dimensions)
  !
  ! Notes:
  !  - On the root rank, the subroutine queries the NetCDF metadata using
  !    nf90_inquire_variable and nf90_inquire_dimension.
  !  - The results are then broadcast to all ranks via p_bcast.
  !  - An error is raised via finish() if var_dimlen is too small.
  !
  SUBROUTINE distrib_inq_var_dims_base(ncid, vname, var_ndims, var_dimlen)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: vname
    INTEGER, INTENT(OUT) :: var_ndims, var_dimlen(:)
    INTEGER :: i
    INTEGER :: temp_var_dimlen(NF90_MAX_VAR_DIMS), var_dimids(NF90_MAX_VAR_DIMS)
    CHARACTER(*), PARAMETER :: routine = modname//'::distrib_inq_var_dims_base'

    IF (this_PE_does_meta_IO) THEN
      CALL nf( &
        nf90_inquire_variable(&
          ncid, get_vid(routine, ncid, vname), &
          ndims = var_ndims, dimids = var_dimids), routine)
      DO i=1, var_ndims
        CALL nf( &
          nf90_inquire_dimension( &
            ncid, var_dimids(i), len = temp_var_dimlen(i)), routine)
      ENDDO
    END IF

    CALL p_bcast(var_ndims, parRootRank, p_comm_work)
    CALL p_bcast(temp_var_dimlen(1:var_ndims), parRootRank, p_comm_work)

    IF (SIZE(var_dimlen) < var_ndims) THEN
      CALL finish(routine, "array size of argument var_dimlen is too small")
    END IF

    var_dimlen(1:var_ndims) = temp_var_dimlen(1:var_ndims)

  END SUBROUTINE distrib_inq_var_dims_base

  ! -------------------------------------------------------------------
  ! Read an integer NetCDF attribute collectively.
  !
  ! The root process reads the integer attribute from the given NetCDF file
  ! and broadcasts the value to all other processes in the communicator.
  ! If `vname` is `'global'`, the attribute is assumed to be a file-global
  ! attribute and `nf90_global` is used as the variable ID.
  !
  ! Inputs:
  !  - ncid: NetCDF file id returned by distrib_nf_open_cont
  !  - vname: variable name ("global" references global attributes)
  !  - name: attribute name
  !
  ! Outputs:
  !  - value: attribute value
  SUBROUTINE distrib_nf_get_att_int_base(ncid, vname, name, value)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(len=*), INTENT(IN) :: vname
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: value

    CHARACTER(*), PARAMETER :: routine = &
      modname//'::distrib_nf_get_att_int_base'

    !-----------------------------------------------------------------------
    ! Only the root rank reads the attribute from file.
    !-----------------------------------------------------------------------
    IF (this_PE_does_meta_IO) THEN
      CALL nf( &
        nf90_get_att(ncid, get_vid(routine, ncid, vname), name, value), routine)
    END IF

    !-----------------------------------------------------------------------
    ! Broadcast the attribute value to all ranks in the communicator.
    !-----------------------------------------------------------------------
    CALL p_bcast(value, parRootRank, p_comm_work)

  END SUBROUTINE distrib_nf_get_att_int_base

  ! -------------------------------------------------------------------
  ! Read a string NetCDF attribute collectively.
  !
  ! The root process reads the string attribute from the given NetCDF file
  ! and broadcasts the value to all other processes in the communicator.
  ! If `vname` is 'global', the attribute is assumed to be a file-global
  ! attribute and `nf90_global` is used as the variable ID.
  !
  ! Inputs:
  !  - ncid: NetCDF file id returned by distrib_nf_open_cont
  !  - vname: variable name ("global" references global attributes)
  !  - name: attribute name
  !
  ! Outputs:
  !  - value: attribute value
  SUBROUTINE distrib_nf_get_att_string_base(ncid, vname, name, value)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(len=*), INTENT(IN) :: vname
    CHARACTER(len=*), INTENT(IN) :: name
    CHARACTER(len=*), INTENT(OUT) :: value

    CHARACTER(*), PARAMETER :: routine = &
        modname//'::distrib_nf_get_att_string_base'

    INTEGER :: att_len, varid
    CHARACTER(len=:), ALLOCATABLE :: tmp_value

    !-----------------------------------------------------------------------
    ! Only the root rank reads the attribute from file.
    !-----------------------------------------------------------------------
    IF (this_PE_does_meta_IO) THEN

      ! Get length of attribute
      varid = get_vid(routine, ncid, vname)
      CALL nf( &
        nf90_inquire_attribute(ncid, varid, name, len=att_len), routine)

    END IF

    !-----------------------------------------------------------------------
    ! Broadcast the attribute length and allocate the receive buffer
    !-----------------------------------------------------------------------
    CALL p_bcast(att_len, parRootRank, p_comm_work)
    ALLOCATE(CHARACTER(len=att_len) :: tmp_value)
    tmp_value(:) = " "

    !-----------------------------------------------------------------------
    ! Read attribute
    !-----------------------------------------------------------------------
    IF (this_PE_does_meta_IO) THEN
      CALL nf(nf90_get_att(ncid, varid, name, tmp_value), routine)
    END IF

    !-----------------------------------------------------------------------
    ! Broadcast the attribute length first, then the string itself
    !-----------------------------------------------------------------------
    CALL p_bcast(tmp_value, parRootRank, p_comm_work)

    ! Return value
    value = TRIM(tmp_value)

    ! Deallocate temporary variable
    DEALLOCATE(tmp_value)

  END SUBROUTINE distrib_nf_get_att_string_base

  !-----------------------------------------------------------------------------
  ! Basic read routines for 1D/2D INT/FLOAT/DOUBLE arrays
  !-----------------------------------------------------------------------------

  SUBROUTINE distrib_read_base_1D_int(ncid, vname, buffer, basic_data_index)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: vname
    INTEGER, INTENT(INOUT) :: buffer(:)
    INTEGER, INTENT(IN) :: basic_data_index
    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_base_1D_int'

    CALL check_distrib_read_args(basic_data_index, SIZE(buffer, 1), routine)

    IF (this_PE_does_IO) THEN
      CALL nf( &
        nf90_get_var(&
          ncid, get_vid(routine, ncid, vname, ref_vtype=NF90_INT), buffer(:), &
          (/basic_data(basic_data_index)%io_chunk%first/), &
          (/basic_data(basic_data_index)%io_chunk%size/)), routine)
    END IF

  END SUBROUTINE distrib_read_base_1D_int

  SUBROUTINE distrib_read_base_2D_int( &
    ncid, vname, buffer, basic_data_index, start_ext_dim)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: vname
    INTEGER, INTENT(INOUT) :: buffer(:,:)
    INTEGER, INTENT(IN) :: basic_data_index
    INTEGER, INTENT(IN), OPTIONAL :: start_ext_dim(1)
    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_base_2D_int'

    INTEGER :: start(2)

    CALL check_distrib_read_args(basic_data_index, SIZE(buffer, 1), routine)

    IF (this_PE_does_IO) THEN
      start(1) = basic_data(basic_data_index)%io_chunk%first
      IF (PRESENT(start_ext_dim)) THEN
        start(2:) = start_ext_dim(:)
      ELSE
        start(2:) = 1
      END IF
      CALL nf( &
        nf90_get_var(&
          ncid, get_vid(routine, ncid, vname, ref_vtype=NF90_INT), &
          buffer(:,:), start, &
          (/basic_data(basic_data_index)%io_chunk%size, SIZE(buffer,2)/)), &
          routine)
    END IF

  END SUBROUTINE distrib_read_base_2D_int

  SUBROUTINE distrib_read_base_3D_int( &
    ncid, vname, buffer, basic_data_index, start_ext_dim)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: vname
    INTEGER, INTENT(INOUT) :: buffer(:,:,:)
    INTEGER, INTENT(IN) :: basic_data_index
    INTEGER, INTENT(IN), OPTIONAL :: start_ext_dim(2)
    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_base_3D_int'

    INTEGER :: start(3)

    CALL check_distrib_read_args(basic_data_index, SIZE(buffer, 1), routine)

    IF (this_PE_does_IO) THEN
      start(1) = basic_data(basic_data_index)%io_chunk%first
      IF (PRESENT(start_ext_dim)) THEN
        start(2:) = start_ext_dim(:)
      ELSE
        start(2:) = 1
      END IF
      CALL nf( &
        nf90_get_var(&
          ncid, get_vid(routine, ncid, vname, ref_vtype=NF90_INT), &
          buffer, start, &
          (/basic_data(basic_data_index)%io_chunk%size, &
            SIZE(buffer,2), SIZE(buffer,3)/)), routine)
    END IF

  END SUBROUTINE distrib_read_base_3D_int

  SUBROUTINE distrib_read_base_1D_sp(ncid, vname, buffer, basic_data_index)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: vname
    REAL(sp), INTENT(INOUT) :: buffer(:)
    INTEGER, INTENT(IN) :: basic_data_index
    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_base_1D_sp'

    CALL check_distrib_read_args(basic_data_index, SIZE(buffer, 1), routine)

    IF (this_PE_does_IO) THEN
      CALL nf( &
        nf90_get_var(&
          ncid, get_vid(routine, ncid, vname), buffer(:), &
          (/basic_data(basic_data_index)%io_chunk%first/), &
          (/basic_data(basic_data_index)%io_chunk%size/)), routine)
    END IF

  END SUBROUTINE distrib_read_base_1D_sp

  SUBROUTINE distrib_read_base_2D_sp( &
    ncid, vname, buffer, basic_data_index, start_ext_dim)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: vname
    REAL(sp), INTENT(INOUT) :: buffer(:,:)
    INTEGER, INTENT(IN) :: basic_data_index
    INTEGER, INTENT(IN), OPTIONAL :: start_ext_dim(1)
    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_base_2D_sp'

    INTEGER :: start(2)

    CALL check_distrib_read_args(basic_data_index, SIZE(buffer, 1), routine)

    IF (this_PE_does_IO) THEN
      start(1) = basic_data(basic_data_index)%io_chunk%first
      IF (PRESENT(start_ext_dim)) THEN
        start(2:) = start_ext_dim(:)
      ELSE
        start(2:) = 1
      END IF
      CALL nf( &
        nf90_get_var(&
          ncid, get_vid(routine, ncid, vname), buffer, start, &
          (/basic_data(basic_data_index)%io_chunk%size, SIZE(buffer,2)/)), &
        routine)
    END IF

  END SUBROUTINE distrib_read_base_2D_sp

  SUBROUTINE distrib_read_base_3D_sp( &
    ncid, vname, buffer, basic_data_index, start_ext_dim)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: vname
    REAL(sp), INTENT(INOUT) :: buffer(:,:,:)
    INTEGER, INTENT(IN) :: basic_data_index
    INTEGER, INTENT(IN), OPTIONAL :: start_ext_dim(2)
    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_base_3D_sp'

    INTEGER :: start(3)

    CALL check_distrib_read_args(basic_data_index, SIZE(buffer, 1), routine)

    IF (this_PE_does_IO) THEN
      start(1) = basic_data(basic_data_index)%io_chunk%first
      IF (PRESENT(start_ext_dim)) THEN
        start(2:) = start_ext_dim(:)
      ELSE
        start(2:) = 1
      END IF
      CALL nf( &
        nf90_get_var(&
          ncid, get_vid(routine, ncid, vname), buffer, start, &
          (/basic_data(basic_data_index)%io_chunk%size, &
           SIZE(buffer,2), SIZE(buffer,3)/)), routine)
    END IF

  END SUBROUTINE distrib_read_base_3D_sp

  SUBROUTINE distrib_read_base_1D_dp(ncid, vname, buffer, basic_data_index)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: vname
    REAL(dp), INTENT(INOUT) :: buffer(:)
    INTEGER, INTENT(IN) :: basic_data_index
    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_base_1D_dp'

    CALL check_distrib_read_args(basic_data_index, SIZE(buffer, 1), routine)

    IF (this_PE_does_IO) THEN
      CALL nf( &
        nf90_get_var(&
          ncid, get_vid(routine, ncid, vname), buffer, &
          (/basic_data(basic_data_index)%io_chunk%first/), &
          (/basic_data(basic_data_index)%io_chunk%size/)), routine)
    END IF

  END SUBROUTINE distrib_read_base_1D_dp

  SUBROUTINE distrib_read_base_2D_dp( &
    ncid, vname, buffer, basic_data_index, start_ext_dim)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: vname
    REAL(dp), INTENT(INOUT) :: buffer(:,:)
    INTEGER, INTENT(IN) :: basic_data_index
    INTEGER, INTENT(IN), OPTIONAL :: start_ext_dim(1)
    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_base_2D_dp'

    INTEGER :: start(2)

    CALL check_distrib_read_args(basic_data_index, SIZE(buffer, 1), routine)

    IF (this_PE_does_IO) THEN
      start(1) = basic_data(basic_data_index)%io_chunk%first
      IF (PRESENT(start_ext_dim)) THEN
        start(2:) = start_ext_dim(:)
      ELSE
        start(2:) = 1
      END IF
      CALL nf( &
        nf90_get_var(&
          ncid, get_vid(routine, ncid, vname), buffer, start, &
          (/basic_data(basic_data_index)%io_chunk%size, SIZE(buffer,2)/)), &
          routine)
    END IF

  END SUBROUTINE distrib_read_base_2D_dp

  SUBROUTINE distrib_read_base_3D_dp( &
    ncid, vname, buffer, basic_data_index, start_ext_dim)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: vname
    REAL(dp), INTENT(INOUT) :: buffer(:,:,:)
    INTEGER, INTENT(IN) :: basic_data_index
    INTEGER, INTENT(IN), OPTIONAL :: start_ext_dim(2)
    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_base_3D_dp'

    INTEGER :: start(3)

    CALL check_distrib_read_args(basic_data_index, SIZE(buffer, 1), routine)

    IF (this_PE_does_IO) THEN
      start(1) = basic_data(basic_data_index)%io_chunk%first
      IF (PRESENT(start_ext_dim)) THEN
        start(2:) = start_ext_dim(:)
      ELSE
        start(2:) = 1
      END IF
      CALL nf( &
        nf90_get_var(&
          ncid, get_vid(routine, ncid, vname), buffer, start, &
          (/basic_data(basic_data_index)%io_chunk%size, &
            SIZE(buffer,2), SIZE(buffer,3)/)), routine)
    END IF

  END SUBROUTINE distrib_read_base_3D_dp

  !-----------------------------------------------------------------------------
  ! Basic read routines for root-based 2D INT arrays
  !-----------------------------------------------------------------------------

  SUBROUTINE distrib_read_root_base_1D_int(ncid, vname, buffer, start, count)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(len=*), INTENT(IN) :: vname
    INTEGER, INTENT(INOUT) :: buffer(:)
    INTEGER, INTENT(IN) :: start(:), count(:)

    CHARACTER(*), PARAMETER :: &
      routine = modname//'::distrib_read_root_base_1D_int'

    !-----------------------------------------------------------------------
    ! Only the root rank reads the array from the file.
    ! (even though the data is only broadcasted by root, we still have to
    !  access the data on all IO processes, because variable access may be
    !  configured for collective access)
    !-----------------------------------------------------------------------
    IF (this_PE_does_meta_IO) THEN
      CALL nf( &
        nf90_get_var( &
          ncid, get_vid(routine, ncid, vname, ref_vtype=NF90_INT), &
          buffer, start(:), count(:)), routine)
    END IF

    !-----------------------------------------------------------------------
    ! Broadcast the array to all ranks.
    !-----------------------------------------------------------------------
    CALL p_bcast(buffer, parRootRank, p_comm_work)

  END SUBROUTINE distrib_read_root_base_1D_int

  ! ----------------------------------------------------------------------------
  ! initialise internal static module variables
  SUBROUTINE init_module_vars(caller)
    CHARACTER(*), INTENT(IN) :: caller
    INTEGER :: temp_n, temp_stride
#if defined (HAVE_PARALLEL_NETCDF) && !defined (NOMPI)
    INTEGER :: ierr, myColor
#endif

    IF (n_io_proc /= -1) &
      CALL finish(caller, "internal variables have already been initialised")

    ! 1) compute io_stride (must be > 0)

    ! if user provided io stride
    IF (io_process_stride > 0) THEN
      io_stride = MAX(1, MODULO(io_process_stride, p_n_work))
    ELSE
      ! no user provided io stride
      !   --> do IO on around SQRT(p_n_work) processes
      temp_n = NINT(SQRT(REAL(p_n_work)))
      temp_stride = (p_n_work + temp_n - 1) / temp_n
      ! improve io process stride by rounding to the next power of two
      io_stride = 2**CEILING(LOG(REAL(temp_stride))/LOG(2.))
    END IF

    ! 2) rotation offset (in [0, io_stride-1])

    ! if user provided rank rotation
    IF (io_process_rotate > 0) THEN
      rotate = MODULO(io_process_rotate, io_stride)
    ELSE
      rotate = 0
    END IF

    ! 3) compute number of IO processes

    n_io_proc = (p_n_work - rotate + io_stride - 1) / io_stride

    ! 4) decide whether this rank is an IO rank

    this_PE_does_IO = (MOD(p_pe_work, io_stride) == rotate)

    ! 5) determine root rank (highest IO rank)

    parRootRank = rotate + (n_io_proc - 1) * io_stride

#if defined (HAVE_PARALLEL_NETCDF) && !defined (NOMPI)
    myColor = MERGE(1, MPI_UNDEFINED, this_PE_does_IO)
    IF (io_comm .EQ. MPI_COMM_NULL) &
      CALL MPI_Comm_split(p_comm_work, myColor, 0, io_comm, ierr)

    ! In case of parallel netcdf all process in io_comm take part in all
    ! file operations (including reading or meta data), because the operation
    ! can be implemented as collectives.
    this_PE_does_meta_IO = this_PE_does_IO
#else

    ! In case of serial netcdf only the root processes access the files, if
    ! meta data is requested.
    this_PE_does_meta_IO = (p_pe_work == parRootRank)
#endif
  END SUBROUTINE init_module_vars

  !-----------------------------------------------------------------------------
  ! Utility routines
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !   Validate that a given basic_data_index refers to a valid and
  !   initialized entry in the global basic_data array. This is a defensive
  !   runtime check used throughout the distributed I/O infrastructure to
  !   detect programming errors and uninitialized handles early.
  !
  ! Input:
  !  - basic_data_index: Index into the global basic_data(:) array.
  !  - routine:  Name of the calling routine, used for diagnostic messages.
  !
  ! Behavior:
  !   - Checks that the global array `basic_data` is allocated.
  !   - Verifies that `basic_data_index` lies within valid bounds.
  !   - Ensures that the reference counter `n_ref` of the pointed
  !     basic_data entry is at least one (i.e. the structure is in use).
  SUBROUTINE check_basic_data_index(basic_data_index, routine)
    INTEGER, INTENT(IN) :: basic_data_index
    CHARACTER(*), INTENT(IN) :: routine

    IF (.NOT. ALLOCATED(basic_data)) &
      CALL finish(routine, "basic_data not allocated")
    IF (basic_data_index .GT. SIZE(basic_data) &
        .OR. basic_data_index .LT. 1) &
      CALL finish(routine, "invalid basic_data_index")
    IF (basic_data(basic_data_index)%n_ref .LT. 1) &
      CALL finish(routine, "invalid basic_data reference counter")

  END SUBROUTINE check_basic_data_index

  !-----------------------------------------------------------------------
  !   Perform consistency checks on arguments passed to distributed read
  !   routines. Ensures that the given receive buffer is dimensioned in
  !   accordance with the basic read decomposition information.
  !
  ! Input:
  ! - basic_data_index: Index of the basic_data entry defining the
  !                     current read decomposition.
  ! - buffer_dim1_size: Size of the first dimension of the local receive
  !                     buffer that will hold data read from file.
  ! - routine: Name of the calling routine, used in diagnostic messages.
  !
  ! Behavior:
  !   - Invokes check_basic_data_index to ensure `basic_data_index` is valid.
  !   - Compares `buffer_dim1_size` with
  !       basic_data(basic_data_index)%io_chunk%size
  !     and aborts if the sizes do not match.
  SUBROUTINE check_distrib_read_args( &
    basic_data_index, buffer_dim1_size, routine)
    INTEGER, INTENT(IN) :: basic_data_index
    INTEGER, INTENT(IN) :: buffer_dim1_size
    CHARACTER(*), INTENT(IN) :: routine

    CALL check_basic_data_index(basic_data_index, routine)

    IF (buffer_dim1_size /= basic_data(basic_data_index)%io_chunk%size) &
      CALL finish( &
        routine, "mismatching size of first dimension of receive buffer")
  END SUBROUTINE check_distrib_read_args

  !-----------------------------------------------------------------------
  !   Obtain the NetCDF variable ID (varid) for a variable named `vname`
  !   in an open NetCDF file, and verify that its type matches the expected
  !   reference type `ref_vtype`.
  !
  ! Input:
  ! - caller: Name of the calling routine for error messages.
  ! - ncid: NetCDF file ID as returned by distrib_nf_open_base.
  ! - vname: Name of the variable to be located in the file.
  !
  ! Optional input:
  ! - ref_vtype: Expected NetCDF data type
  !              (e.g. NF90_REAL, NF90_DOUBLE, NF90_INT).
  !
  ! Return:
  ! - get_vid: The NetCDF variable ID corresponding to `vname`, if found (and of
  !            matching type). nf90_global is returned if `vname` = "global"
  FUNCTION get_vid(caller, ncid, vname, ref_vtype)
    CHARACTER(*), INTENT(IN) :: caller
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: vname
    INTEGER, OPTIONAL, INTENT(IN) :: ref_vtype

    INTEGER :: get_vid

    INTEGER :: vtype

    IF (TRIM(vname) == "global") THEN

      get_vid = nf90_global

    ELSE

      CALL nf(nf90_inq_varid(ncid, vname, get_vid), caller)

      ! if this routine is expected to check the type of the variable
      IF (PRESENT(ref_vtype)) THEN

        CALL nf(nf90_inquire_variable(ncid, get_vid, xtype = vtype), caller)
        IF (ref_vtype /= vtype) &
            CALL finish(caller, "invalid type of variable " // TRIM(vname))
      END IF

    END IF

  END FUNCTION get_vid

END MODULE mo_read_netcdf_distributed_base
