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

MODULE mo_reader_cams

  USE mo_kind,                    ONLY: wp, i8
  USE mo_parallel_config,         ONLY: get_nproma
  USE mo_exception,               ONLY: finish
  USE mo_reader_abstract,         ONLY: t_abstract_indexed_reader
  USE mo_reader_util,             ONLY: read_timestamps_from_netcdf
  USE mo_impl_constants,          ONLY: n_camsaermr
  USE mo_io_units,                ONLY: FILENAME_MAX
  USE mo_model_domain,            ONLY: t_patch
  USE mo_netcdf_errhandler,       ONLY: nf
  USE mo_netcdf
  USE mtime,                      ONLY: julianday, getJulianDayFromDatetime, datetime, OPERATOR(>)
  USE mo_mpi,                     ONLY: my_process_is_mpi_workroot, process_mpi_root_id, p_comm_work, p_bcast
  USE mo_read_netcdf_distributed, ONLY: distrib_nf_open, distrib_read, distrib_nf_close, idx_blk_time
  USE fortran_support,            ONLY: t_ptr_4d_wp
  USE mo_radiation_config,        ONLY: irad_aero, iRadAeroCAMSclim, iRadAeroCAMStd

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_cams_reader

  TYPE, EXTENDS(t_abstract_indexed_reader) :: t_cams_reader

    TYPE(t_patch), POINTER      :: p_patch => NULL()
    CHARACTER(len=NF90_MAX_NAME)  :: varnames(n_camsaermr+1)
    CHARACTER(len=FILENAME_MAX) :: filename
    INTEGER                     :: dist_fileid, nlev_cams
    LOGICAL                     :: lopened = .FALSE.

    TYPE(julianday), ALLOCATABLE :: times(:)
    INTEGER                     :: index = 1

  CONTAINS

    PROCEDURE :: init            => cams_init_reader
    PROCEDURE :: prev            => cams_prev
    PROCEDURE :: next            => cams_next
    PROCEDURE :: goto            => cams_goto
    PROCEDURE :: seek            => cams_seek
    PROCEDURE :: is_valid        => cams_is_valid

    PROCEDURE :: read            => cams_read
    PROCEDURE :: get_julian_day  => cams_get_julian_day
    PROCEDURE :: get_index       => cams_get_index

    PROCEDURE :: deinit          => cams_deinit_reader

    PROCEDURE :: get_nblks       => cams_get_nblks
    PROCEDURE :: get_npromz      => cams_get_npromz

  END TYPE t_cams_reader

  CHARACTER(len=*), PARAMETER :: modname = 'mo_reader_cams'

CONTAINS

  SUBROUTINE cams_init_reader(this, p_patch, filename)
    CLASS(t_cams_reader),    INTENT(inout) :: this
    TYPE(t_patch),      TARGET, INTENT(in   ) :: p_patch
    CHARACTER(len=*),           INTENT(in   ) :: filename

    CHARACTER(len=*), PARAMETER :: routine = 'cams_init_reader'

    INTEGER :: ntimes

    this%filename = TRIM(filename)

    IF (irad_aero == iRadAeroCAMSclim) THEN

     this%varnames(1)  = "Sea_Salt_bin1"              ! Sea salt bin1 aerosol mass mixing ratio kg kg-1
     this%varnames(2)  = "Sea_Salt_bin2"              ! Sea salt bin2 aerosol mass mixing ratio kg kg-1
     this%varnames(3)  = "Sea_Salt_bin3"              ! Sea salt bin3 aerosol mass mixing ratio kg kg-1
     this%varnames(4)  = "Mineral_Dust_bin1"          ! Mineral dust bin1 aerosol mass mixing ratio kg kg-1
     this%varnames(5)  = "Mineral_Dust_bin2"          ! Mineral dust bin2 aerosol mass mixing ratio kg kg-1
     this%varnames(6)  = "Mineral_Dust_bin3"          ! Mineral dust bin3 aerosol mass mixing ratio kg kg-1
     this%varnames(7)  = "Organic_Matter_hydrophilic" ! Hydrophilic organic matter aerosol mass mixing ratio  kg kg-1
     this%varnames(8)  = "Organic_Matter_hydrophobic" ! Hydrophobic organic matter aerosol mass mixing ratio kg kg-1
     this%varnames(9)  = "Black_Carbon_hydrophilic"   ! Hydrophilic black carbon aerosol mass mixing ratio kg kg-1
     this%varnames(10) = "Black_Carbon_hydrophobic"   ! Hydrophobic black carbon aerosol mass mixing ratio kg kg-1
     this%varnames(11) = "Sulfates"                   ! Sulfates aerosol mass mixing ratio kg kg-1
     this%varnames(12) = "pressure"                   ! air_pressure Pressure at layer centres (Pa)

     this%nlev_cams = 21

  ELSEIF (irad_aero == iRadAeroCAMStd) THEN

     this%varnames(1)  = "aermr01" ! Sea_Salt_bin1 mixing ratio (kg/kg)
     this%varnames(2)  = "aermr02" ! Sea_Salt_bin2 mixing ratio (kg/kg)
     this%varnames(3)  = "aermr03" ! Sea_Salt_bin3 mixing ratio (kg/kg)
     this%varnames(4)  = "aermr04" ! Mineral_Dust_bin1 mixing ratio (kg/kg)
     this%varnames(5)  = "aermr05" ! Mineral_Dust_bin2 mixing ratio (kg/kg)
     this%varnames(6)  = "aermr06" ! Mineral_Dust_bin3 mixing ratio (kg/kg)
     this%varnames(7)  = "aermr07" ! Organic_Matter_hydrophilic mixing ratio (kg/kg)
     this%varnames(8)  = "aermr08" ! Organic_Matter_hydrophobic mixing ratio (kg/kg)
     this%varnames(9)  = "aermr09" ! Black_Carbon_hydrophilic mixing ratio (kg/kg)
     this%varnames(10) = "aermr10" ! Black_Carbon_hydrophobic mixing ratio (kg/kg)
     this%varnames(11) = "aermr11" ! Sulfates mixing ratio (kg/kg)
     this%varnames(12) = "pres"    ! Pressure at base of layer (Pa)

     this%nlev_cams = 137

    ENDIF

    this%p_patch => p_patch

    IF (.NOT. this%lopened) THEN
      this%dist_fileid = distrib_nf_open(TRIM(this%filename))
      this%lopened = .TRUE.

      IF (my_process_is_mpi_workroot()) THEN
        CALL nf(read_timestamps_from_netcdf(TRIM(this%filename), this%times), routine)
        ntimes = SIZE(this%times)
      END IF

      CALL p_bcast(ntimes, process_mpi_root_id, p_comm_work)
      IF (.NOT. ALLOCATED(this%times)) THEN
        ALLOCATE(this%times(ntimes))
      END IF

      CALL p_bcast(this%times(:)%day, process_mpi_root_id, p_comm_work)
      CALL p_bcast(this%times(:)%ms, process_mpi_root_id, p_comm_work)
    ENDIF

  END SUBROUTINE cams_init_reader

  SUBROUTINE cams_prev (this)
    CLASS(t_cams_reader), INTENT(INOUT) :: this

    this%index = this%index - 1
  END SUBROUTINE cams_prev

  SUBROUTINE cams_next (this)
    CLASS(t_cams_reader), INTENT(INOUT) :: this

    this%index = this%index + 1
  END SUBROUTINE cams_next

  SUBROUTINE cams_goto (this, target_datetime)
    CLASS(t_cams_reader), INTENT(INOUT) :: this
    TYPE(datetime), INTENT(IN) :: target_datetime

    TYPE(julianday) :: jd

    INTEGER :: i

    CALL getJulianDayFromDatetime(target_datetime, jd)

    DO i = 1, SIZE(this%times)
      IF (this%times(i) > jd) EXIT
    END DO

    this%index = i - 1
  END SUBROUTINE cams_goto

  SUBROUTINE cams_seek (this, index)
    CLASS(t_cams_reader), INTENT(INOUT) :: this
    INTEGER(i8), INTENT(IN) :: index

    this%index = INT(index)
  END SUBROUTINE cams_seek

  FUNCTION cams_is_valid (this, msg) RESULT(valid)
    CLASS(t_cams_reader), INTENT(IN) :: this
    CHARACTER(len=*), INTENT(OUT), OPTIONAL :: msg
    LOGICAL :: valid

    valid = (this%index > 0 .AND. this%index <= SIZE(this%times))

    IF (PRESENT(msg)) msg = 'Index out of bounds'
  END FUNCTION cams_is_valid

  SUBROUTINE cams_read(this, varname, dat)
    CLASS(t_cams_reader), INTENT(INOUT)  :: this
    CHARACTER(len=*), INTENT(IN)         :: varname
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: dat(:,:,:,:)
    REAL(wp), ALLOCATABLE, TARGET        :: temp(:,:,:,:)
    TYPE(t_ptr_4d_wp)                    :: tmp(1)
    INTEGER                              :: var_dimlen(3),var_start(3), var_end(3), jt


    ALLOCATE(temp(get_nproma(), this%nlev_cams, this%p_patch%nblks_c, n_camsaermr+1))
    IF (ALLOCATED(dat)) DEALLOCATE(dat)
    ALLOCATE( dat(get_nproma(), this%nlev_cams, this%p_patch%nblks_c, n_camsaermr+1))

    var_dimlen(2) = SIZE(temp, 2) ! number of vertical levels
    var_dimlen(3) = 1             ! number of time steps = 1
    var_start(:)  = (/1, 1, 1/)
    var_end(:)    = var_dimlen(:)
    var_start(3)  = this%index
    var_end(3)    = this%index

    temp(:,:,:,:) = -1.0_wp

    IF (.NOT. this%lopened) &
      CALL finish(modname, 'CAMS climatology file not open!')
    IF (TRIM(varname) /= '') &
      CALL finish(modname, 'CAMS climatology: Only bulk reading of all variables implemented!')

    tmp(1)%p => temp

    DO jt = 1, n_camsaermr+1

       CALL distrib_read(this%dist_fileid, this%varnames(jt), tmp, &
         & (/this%p_patch%cells%dist_io_data/), edim=var_dimlen(2:3), dimo=idx_blk_time, &
         & start_ext_dim=var_start(2:3), end_ext_dim=var_end(2:3))

         dat(:,:,:,jt) = temp(:,:,:,1)
    ENDDO

    DEALLOCATE(temp)

  END SUBROUTINE cams_read

  FUNCTION cams_get_julian_day (this) RESULT(jd)
    CLASS(t_cams_reader), INTENT(IN) :: this
    TYPE(julianday) :: jd

    jd = this%times(this%index)
  END FUNCTION cams_get_julian_day

  FUNCTION cams_get_index (this) RESULT(index)
    CLASS(t_cams_reader), INTENT(IN) :: this
    INTEGER(i8) :: index

    index = this%index
  END FUNCTION cams_get_index

  FUNCTION cams_get_nblks (this) RESULT(nblks)
    CLASS(t_cams_reader), INTENT(in   ) :: this
    INTEGER                                :: nblks
    nblks = this%p_patch%nblks_c
  END FUNCTION cams_get_nblks

  FUNCTION cams_get_npromz (this) RESULT(npromz)
    CLASS(t_cams_reader), INTENT(in   ) :: this
    INTEGER                                :: npromz
    npromz = this%p_patch%npromz_c
  END FUNCTION cams_get_npromz

  SUBROUTINE cams_deinit_reader(this)
    CLASS(t_cams_reader), INTENT(inout) :: this

    IF (ASSOCIATED(this%p_patch)) NULLIFY(this%p_patch)
    IF (this%lopened) THEN
      CALL distrib_nf_close(this%dist_fileid)
    END IF
  END SUBROUTINE cams_deinit_reader

END MODULE mo_reader_cams
