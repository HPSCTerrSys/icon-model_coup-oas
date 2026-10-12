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

! Read anthropogenic emission of CO2

MODULE mo_bc_anthro_emission

  USE mo_kind,            ONLY: wp, i8
  USE mo_exception,       ONLY: finish, message, message_text
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights
  USE mo_run_config,      ONLY: msg_level
  USE mo_grid_config,     ONLY: n_dom
  USE mo_parallel_config, ONLY: nproma
  USE mo_impl_constants,  ONLY: MAX_CHAR_LENGTH,max_dom
  USE mo_time_config,     ONLY: time_config
  USE mo_model_domain,    ONLY: t_patch
  USE mo_cdi,             ONLY: streamOpenRead, streamInqVlist, streamClose, &
    & vlistInqTaxis, streamInqTimestep, taxisInqVdate
#ifdef __SINGLE_PRECISION
  USE mo_cdi,                ONLY: streamReadVarSliceF
#else
  USE mo_cdi,                ONLY: streamReadVarSlice
#endif
  USE mo_util_cdi,           ONLY: cdiGetStringError
  USE mo_mpi,             ONLY: my_process_is_mpi_workroot, p_bcast, &
    &                              process_mpi_root_id, p_comm_work
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
       &                               calculate_time_interpolation_weights


  IMPLICIT NONE
  PRIVATE

  TYPE t_ext_emis
    REAL(wp), CONTIGUOUS, POINTER :: co2ant(:,:,:) => NULL()
  END TYPE t_ext_emis

  TYPE(t_ext_emis), TARGET :: ext_emis(max_dom)

  PUBLIC :: read_bc_anthro_emission
  PUBLIC :: get_current_bc_anthro_emission_year
  PUBLIC :: bc_anthro_emission_time_interpolation

  INTEGER(i8), SAVE :: current_year = -1

  INTEGER :: nyears
  INTEGER :: imonth_beg, imonth_end
  LOGICAL :: lend_of_year

  TYPE(t_time_interpolation_weights) :: tiw_beg
  TYPE(t_time_interpolation_weights) :: tiw_end

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_bc_anthro_emission'

CONTAINS

  SUBROUTINE read_bc_anthro_emission(year,p_patch)
    INTEGER(i8),   INTENT(IN) :: year
    TYPE(t_patch), INTENT(IN) :: p_patch
    INTEGER :: jg
    CHARACTER(len=21) :: fn

    CHARACTER(len=*), PARAMETER :: routine = thismodule//':read_bc_anthro_emission'

    jg = p_patch%id                ! Only relevant if working with several inside nests
    lend_of_year = ( time_config%tc_stopdate%date%month  == 1  .AND. &
      &              time_config%tc_stopdate%date%day    == 1  .AND. &
      &              time_config%tc_stopdate%time%hour   == 0  .AND. &
      &              time_config%tc_stopdate%time%minute == 0  .AND. &
      &              time_config%tc_stopdate%time%second == 0 )

    nyears = time_config%tc_stopdate%date%year - time_config%tc_startdate%date%year + 1
    IF ( lend_of_year ) nyears = nyears - 1

    ! ----------------------------------------------------------------------
    tiw_beg = calculate_time_interpolation_weights(time_config%tc_startdate)
    tiw_end = calculate_time_interpolation_weights(time_config%tc_stopdate)

    IF ( nyears > 1 ) THEN
      imonth_beg = 0
      imonth_end = 13
      IF ( tiw_beg%month1_index == 0 ) imonth_end = tiw_end%month2_index
    ELSE
      imonth_beg = tiw_beg%month1_index
      imonth_end = tiw_end%month2_index
    ENDIF

    IF ( lend_of_year ) imonth_end = 13

    WRITE(message_text,'(a,i2,a,i2)') &
       & ' Allocating CO2_anthro for months ', imonth_beg, ' to ', imonth_end
    CALL message(routine, message_text)

    IF ( imonth_beg > imonth_end ) THEN
      WRITE (message_text, '(a)') 'imonth_beg < imonth_end'
      CALL finish(routine, message_text)
    ENDIF

    IF (n_dom > 1) THEN
      WRITE(fn, '(a,i2.2,a)') 'bc_anthro_emission_DOM', jg, '.nc'
    ELSE
      fn = 'bc_anthro_emission.nc'
    END IF
    IF (my_process_is_mpi_workroot()) THEN
      WRITE(message_text,'(3a,i0)') 'Read CO2 anthropogenic emission from ', TRIM(fn), ' for ', year
      CALL message('',message_text)
    ENDIF
    IF (.NOT.ASSOCIATED(ext_emis(jg)%co2ant)) THEN
      ALLOCATE (ext_emis(jg)%co2ant(nproma, p_patch%nblks_c, imonth_beg:imonth_end))
      !$ACC ENTER DATA CREATE(ext_emis(jg)%co2ant)
    ENDIF
    CALL read_anthro_emission_data(p_patch, ext_emis(jg)%co2ant, TRIM(fn), year)

    IF (jg==n_dom) current_year = year

    !$ACC UPDATE DEVICE(ext_emis(jg)%co2ant) ASYNC(1)

  END SUBROUTINE read_bc_anthro_emission

  SUBROUTINE read_anthro_emission_data(p_patch, dst, fn, y)
!TODO: switch to reading via mo_read_netcdf_distributed?
    TYPE(t_patch), INTENT(in) :: p_patch
    REAL(wp), CONTIGUOUS, INTENT(INOUT) :: dst(:,:,imonth_beg:)
    CHARACTER(len=*), INTENT(IN) :: fn
    INTEGER(i8), INTENT(in) :: y
    REAL(wp), ALLOCATABLE :: zin(:)
    REAL(wp) :: dummy(0)
    INTEGER :: vlID, taxID, tsID, ts_idx, strID, nmiss, vd, vy, vm, ts_found
    LOGICAL :: found_last_ts, lexist
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

    CHARACTER(len=*), PARAMETER :: routine = thismodule//':read_anthro_emission_data'

    IF (my_process_is_mpi_workroot()) THEN
      INQUIRE (file=fn, exist=lexist)
      IF (.NOT.lexist) THEN
        WRITE (message_text, '(3a)') 'Could not open file ', fn, ': run terminated.'
        CALL finish(routine, message_text)
      ENDIF
      strID = streamOpenRead(fn)

      IF ( strID < 0 ) THEN
        CALL cdiGetStringError(strID, cdiErrorText)
        WRITE (message_text, '(4a)') 'Could not open file ', fn, ': ', cdiErrorText
        CALL finish(routine, message_text)
      END IF
      vlID = streamInqVlist(strID)
      taxID = vlistInqTaxis(vlID)
      tsID = 0
      found_last_ts = .FALSE.
      ts_found = 0

      ALLOCATE(zin(p_patch%n_patch_cells_g))
      DO WHILE (.NOT. found_last_ts)
        IF (streamInqTimestep(strID, tsID) == 0) EXIT
        vd = taxisInqVdate(taxID)
        vy = vd/10000
        vm = (vd/100)-vy*100
        ts_idx = -1

        IF (INT(vy,i8) == y-1_i8 .AND. vm == 12) THEN
          IF ( imonth_beg == 0 ) THEN
              ts_idx = 0
              ts_found = ts_found + 1
          END IF
        ELSE IF (INT(vy,i8) == y) THEN
          IF ( vm >= imonth_beg .AND. vm <= imonth_end ) THEN
            ts_idx = vm
            ts_found = ts_found + 1
            IF ( vm == imonth_end ) found_last_ts = .TRUE.
          END IF
        ELSE IF (INT(vy,i8) == y+1_i8 .AND. vm == 1) THEN
          IF ( imonth_end == 13 ) THEN
            ts_idx = 13
            ts_found = ts_found + 1
            found_last_ts = .TRUE.
          END IF
        END IF
        IF (ts_idx /= -1) THEN
#ifdef __SINGLE_PRECISION
          CALL streamReadVarSliceF(strID, 0, 0, zin, nmiss)
#else
          CALL streamReadVarSlice (strID, 0, 0, zin, nmiss)
#endif
          CALL p_bcast(ts_idx, process_mpi_root_id, p_comm_work)
          dst(:,SIZE(dst,2),ts_idx) = 0._wp
          CALL p_patch%comm_pat_scatter_c%distribute(zin, dst(:,:,ts_idx), .FALSE.)
        ENDIF
        tsID = tsID+1
      END DO
      IF (ts_found < imonth_end - imonth_beg + 1) &
          & CALL finish (routine, 'could not read required data from input file')
      ts_idx = -1
      CALL p_bcast(ts_idx, process_mpi_root_id, p_comm_work)
      DEALLOCATE(zin)
      CALL streamClose(strID)
    ELSE
      ts_idx = 0
      DO
        CALL p_bcast(ts_idx, process_mpi_root_id, p_comm_work)
        IF(ts_idx .EQ. -1) EXIT
        dst(:,SIZE(dst,2),ts_idx) = 0._wp
        CALL p_patch%comm_pat_scatter_c%distribute(dummy, dst(:,:,ts_idx), .FALSE.)
      END DO
    END IF
  END SUBROUTINE read_anthro_emission_data

  SUBROUTINE bc_anthro_emission_time_interpolation(tiw, co2ant_out, p_patch)

    TYPE( t_time_interpolation_weights), INTENT(in) :: tiw
    REAL(wp)       , INTENT(inout) :: co2ant_out(:,:)
    TYPE(t_patch)  , INTENT(in)  :: p_patch

    REAL(wp), CONTIGUOUS, POINTER :: co2ant_in(:,:,:)

    INTEGER  :: jc, jb, jg, jce, nblk

    CHARACTER(len=*), PARAMETER :: routine = thismodule//':bc_anthro_emission_time_interpolation'

    jg = p_patch%id

    co2ant_in => ext_emis(jg)%co2ant

    jce  = SIZE(co2ant_out,1)
    nblk = SIZE(co2ant_out,2)

!$OMP PARALLEL DO PRIVATE(jb, jc)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) COPYIN(tiw) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO jb = 1, nblk
      DO jc = 1, jce
        co2ant_out(jc,jb) = tiw%weight1 * co2ant_in(jc,jb,tiw%month1_index) + tiw%weight2 * co2ant_in(jc,jb,tiw%month2_index)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO

  END SUBROUTINE bc_anthro_emission_time_interpolation


  FUNCTION get_current_bc_anthro_emission_year() RESULT(this_year)
    INTEGER(i8) :: this_year
    this_year = current_year
  END FUNCTION get_current_bc_anthro_emission_year


END MODULE mo_bc_anthro_emission
