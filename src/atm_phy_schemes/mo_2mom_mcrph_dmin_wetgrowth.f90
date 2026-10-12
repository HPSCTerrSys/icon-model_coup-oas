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

! Two-moment bulk microphysics after Seifert, Beheng and Blahak
!
! Description:
!
! Provides code to generate a lookup table for the critical wet growth diameter
! of graupel / hail by iteratively solving the heat balance equation
! of a particle undergoing riming, evaporation and collisions with
! other frozen particles. There is also code to write the table to a netcdf file.
!
! Method: See Appendix A of
!
! A. Khain, D. Rosenfeld, A. Pokrovsky, U. Blahak, A. Ryzhkov, 2011:
! The role of CCN in precipitation and hail in a mid-latitude storm as seen
! in simulations using a spectral (bin) microphysics model in a 2D dynamic frame,
! Atmospheric Research, 99, 129-146

!NEC$ options "-finline-max-depth=3 -finline-max-function-size=1000"

MODULE mo_2mom_mcrph_dmin_wetgrowth

  USE mo_kind,             ONLY : wp
  USE mo_2mom_mcrph_types, ONLY : particle
  USE mo_exception,        ONLY : message, finish, message_text

  USE netcdf, ONLY :  &
       NF90_NOERR, &
       NF90_CLOBBER, &
       NF90_NOWRITE, &
       nf90_create, &
       nf90_open, &
       nf90_def_dim, &
       nf90_close, &
       nf90_put_att, &
       nf90_strerror, &
       nf90_global, &
       nf90_enddef, &
       nf90_put_var, &
       nf90_def_var, &
       NF90_NETCDF4, &
       NF90_CLASSIC_MODEL, &
       NF90_DOUBLE

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dmin_wetgrowth_lookupcreate, write_dmin_wetgrowth_table, get_filebase_dmin_lut_file


  CHARACTER(len=*), PARAMETER :: modulename = 'mo_2mom_mcrph_dmin_wetgrowth'
  REAL(wp), PARAMETER :: &
       pi      = 3.14159265358_wp , &    ! pi
       c_p     = 1005._wp         , &    ! specific heat capacity (air) [J/(kg K)]
       Rgas    = 8.3145_wp        , &    ! universal gas constant [J/(mol K)]
       M_w     = 18.01528e-3_wp   , &    ! molecular mass of water [kg/mol]
       M_a     = 28.96e-3_wp      , &    ! molecular mass of air [kg/mol]
       R_w     = Rgas / M_w       , &    ! gas constant for water [J/(kg K)]
       R_a     = Rgas / M_a       , &    ! gas constant for air [J/(kg K)]
       T_0C    = 273.16_wp

  REAL(wp), PARAMETER :: miss_value = -999.99_wp
  REAL(wp), PARAMETER :: miss_thresh = -900.0_wp

CONTAINS

  !**********************************************************************
  !**********************************************************************

  !**********************************************************************
  !
  ! Subroutine for creating a lookup table for the critical wet growth diameter
  ! of graupel / hail by iteratively solving the heat balance equation
  ! of a particle undergoing riming, evaporation and collisions with
  ! other frozen particles.
  !
  ! Method: See Appendix A of
  !
  ! A. Khain, D. Rosenfeld, A. Pokrovsky, U. Blahak, A. Ryzhkov, 2011:
  ! The role of CCN in precipitation and hail in a mid-latitude storm as seen
  ! in simulations using a spectral (bin) microphysics model in a 2D dynamic frame,
  ! Atmospheric Research, 99, 129-146
  !
  !**********************************************************************

  SUBROUTINE dmin_wetgrowth_lookupcreate (phail, Dmin_lut, n_pres_o, n_temp_o, n_qw_o, n_qi_o, pres, temp, qw, qi)

    IMPLICIT NONE

    TYPE(PARTICLE), INTENT(in) :: phail
     !.. Lookup-table for Minimum diameter for wet-growth
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(out) :: Dmin_lut     ! m
    ! .. Table vectors for pressure, temperature, qw and qi
    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: pres, temp, qw, qi ! Pa, deg C, kg/m3, kg/m3

    !.. Array sizes for output as information to the calling routine:
    INTEGER, INTENT(out) :: &
         n_pres_o,     &  ! (8)  Number of Pressure values; range [300,1000]*1.e2 Pa
         n_temp_o,     &  ! (23) Number of Temperature values; range [-30,0] degrees Celsius
         n_qw_o,       &  ! (48) Number of bulk liquid water values; range [0.6,10]*1.e-3 kg/m3
         n_qi_o           ! (48) Number of bulk ice water values; range [0.6,10]*1.e-3 kg/m3

    INTEGER, PARAMETER :: maxiter = 100
    !.. Array sizes as internal fixed parameters:
    INTEGER, PARAMETER :: &
         n_pres = 8,   &  ! (8)  Number of Pressure values; range [300,1000]*1.e2 Pa
         n_temp = 23,  &  ! (23) Number of Temperature values; range [-30,0] degrees Celsius
         n_qw   = 48,  &  ! (48) Number of bulk liquid water values; range [0.6,10]*1.e-3 kg/m3
         n_qi   = 48      ! (48) Number of bulk ice water values; range [0.6,10]*1.e-3 kg/m3


    CHARACTER(len=*), PARAMETER :: routine = TRIM(modulename)//'::dmin_wetgrowth_lookupcreate'

    !..Initial and interval values in double precision:
    REAL(wp) :: p0, delta_p, T0, delta_T, &
         qw0, delta_qw, qi0, delta_qi

    INTEGER               :: i, j, k, l, ii, nerr
    INTEGER,  ALLOCATABLE :: iiv(:,:), rowcount(:), lfill(:)
    REAL(wp)              :: Dold, Dnew, D0
    REAL(wp), ALLOCATABLE :: Doldv(:,:), Dnewv(:,:), D0v(:,:), Dmin(:,:,:,:), rowmean(:)

    message_text(:) = ' '
    WRITE(message_text,'(a)') 'creating wet-growth Dmin table for: '//phail%name
    CALL message (TRIM(routine), TRIM(message_text))

    !=================================================================================
    ! Initializations:

    n_pres_o = n_pres
    n_temp_o = n_temp
    n_qw_o   = n_qw
    n_qi_o   = n_qi

    p0 = 300.e2_wp  ! Starting value for the table vector for pressure
    T0 = -30._wp    !    "-"  for temperature
    qw0 = 0.6e-3_wp !    "-"  for supercooled LWC
    qi0 = 0.6e-3_wp !    "-"  for IWC (ice, snow)

    delta_p  = 100.e2_wp  ! equi-distant increments for the table vectors
    delta_T = 1._wp       ! may change to non-equidistant below, if n_temp = 23
    delta_qw = 0.2e-3_wp
    delta_qi = 0.2e-3_wp

    ! .. Allocate variables, deallocate first for safety:
    IF (ALLOCATED(pres))     DEALLOCATE (pres)
    IF (ALLOCATED(temp))     DEALLOCATE (temp)
    IF (ALLOCATED(qw))       DEALLOCATE (qw)
    IF (ALLOCATED(qi))       DEALLOCATE (qi)
    IF (ALLOCATED(Dmin_lut)) DEALLOCATE (Dmin_lut)

    ALLOCATE( pres(n_pres) )                     ; pres = 0._wp
    ALLOCATE( temp(n_temp) )                     ; temp = 0._wp
    ALLOCATE( qw(n_qw) )                         ; qw = 0._wp
    ALLOCATE( qi(n_qi) )                         ; qi = 0._wp
    ALLOCATE( Dmin(n_qw,n_qi,n_pres,n_temp) )    ; Dmin = 0._wp     ! suitable index order for vectorized iteration

    ! .. Built table vectors before the iteration loop:
    IF (n_temp == 23) THEN ! special case on non-equi. distant T
      ! .. Non-equidistant temperature values:
      delta_T = 2._wp
      DO j = 1, n_temp
        IF (j == 1) THEN
          temp(j) = T0
        ELSE
          temp(j) = temp(j-1) + delta_T
        END IF
        IF (temp(j)+T_0C > 263.0_wp) delta_T = 1.0_wp
        IF (temp(j)+T_0C > 271.0_wp) delta_T = 0.5_wp
      END DO
    ELSE
      ! .. Equidistant temperature values [Celsius]:
      DO j = 1, n_temp
        temp(j) = T0 + (j-1)*delta_T
      END DO
    END IF

    ! .. Pressure values:
    DO i = 1, n_pres
      pres(i) = p0 + (i-1)*delta_p
    END DO

    ! .. Liquid water values
    DO k = 1, n_qw
      qw(k) = qw0 + (k-1)*delta_qw
    END DO

    ! .. Ice water values
    DO l = 1, n_qi
      qi(l) = qi0 + (l-1)*delta_qi
    END DO

    !=================================================================================
    ! .. Fixpoint iteration for the equilibrium diameter Dnew

#if defined (__SX__) || defined (__NEC_VH__) || defined (__NECSX__)

    !=================================================================================
    ! .. Vectorized version:

    ALLOCATE(iiv(n_qw,n_qi), Doldv(n_qw,n_qi), Dnewv(n_qw,n_qi), D0v(n_qw,n_qi))

!$OMP parallel do private(i,j,k,l,iiv,Doldv,D0v,Dnewv)
    DO j = 1, n_temp

      IF (temp(j) < 0.0_wp) THEN

        DO i = 1, n_pres

          ! .. 1) initialization of the iteration:
          iiv(:,:) = 0
          Doldv(:,:) = 5.e-3_wp
          D0v(:,:) = 1.0_wp

          ! .. Actual vectorizable iteration:
          iteration: DO WHILE (ANY(ABS(D0v) > 1.e-4_wp .AND. iiv < maxiter))

            DO l = 1, n_qi
!$NEC ivdep
              DO k = 1, n_qw

                IF ( ABS(D0v(k,l)) > 1.e-4_wp .AND. iiv(k,l) < maxiter) THEN
                  iiv(k,l) = iiv(k,l) + 1
                  Dnewv(k,l) = dwg_fpi(phail%a_geo,phail%b_geo,phail%a_vel,phail%b_vel,Doldv(k,l),temp(j),pres(i),qw(k),qi(l))
                  IF (Dnewv(k,l) <= 0.0_wp) THEN
                    ! This means that ice and supercooled liquid are so cold
                    ! that wet growth cannot occur.
                    iiv(k,l) = maxiter + 1;
                  END IF
                  D0v(k,l)   = ( Dnewv(k,l)-Doldv(k,l) ) / Doldv(k,l)
                  Doldv(k,l) = Dnewv(k,l)
                END IF

              END DO
            END DO

          END DO iteration

          DO l = 1, n_qi
!$NEC ivdep
            DO k = 1, n_qw
              IF (ABS(D0v(k,l)) > 1.e-4_wp) THEN
                Dmin(k,l,i,j) = miss_value
              ELSE
                Dmin(k,l,i,j) = MIN(Dnewv(k,l),999.99_wp)
              END IF
            END DO
          END DO

        END DO

      ELSE

        ! Enforce correct behaviour for Tvec(i) = 0 grad C, because
        ! iteration would lead to DIV0 in the first step otherwise:
        Dmin(:,:,:,j) = 0.0_wp

      END IF

    END DO
!$OMP end parallel do

    DEALLOCATE(iiv,Doldv,Dnewv,D0v)

#else

    !=================================================================================
    ! .. Unvectorized version:

!$OMP parallel do private(i,j,k,l,ii,Dold,D0,Dnew)
    DO j = 1, n_temp
      DO i = 1, n_pres
        DO l = 1, n_qi
          DO k = 1, n_qw

            ! .. Fixpoint iteration (vectorized) for the equilibrium diameter Dnew

            ! .. 1) initialization of the iteration:
            ii = 0
            Dold = 5.e-3_wp
            D0 = 1.0_wp

            IF (temp(j) < 0.0_wp) THEN

              ! .. Actual iteration:
             DO WHILE (ABS(D0) > 1.e-4_wp .AND. ii < maxiter)
                ii = ii +1
                Dnew = dwg_fpi(phail%a_geo,phail%b_geo,phail%a_vel,phail%b_vel,Dold,temp(j),pres(i),qw(k),qi(l))
                IF (Dnew <= 0.0_wp) THEN
                  ! This means that ice and supercooled liquid are so cold
                  ! that wet growth cannot occur.
                  ii = maxiter + 1;
                END IF
                D0 = ( Dnew - Dold ) / Dold
                Dold = Dnew
              END DO

              IF (ABS(D0) > 1.e-4_wp) THEN
                Dmin(k,l,i,j) = miss_value
              ELSE
                Dmin(k,l,i,j) = MIN(Dnew,999.99_wp)
              END IF

            ELSE

              ! Enforce correct behaviour for Tvec(i) = 0 grad C, because
              ! iteration would lead to DIV0 in the first step otherwise:
              Dmin(k,l,i,j) = 0.0_wp

            END IF


          END DO
        END DO
      END DO
    END DO
!$OMP end parallel do

#endif

    !=================================================================================
    ! .. Check table:

    IF (ANY(Dmin < miss_thresh)) THEN

      nerr = COUNT((Dmin < miss_thresh))
      WRITE(message_text,'(a,i0,a)') 'newly created wet-growth Dmin table for '//TRIM(phail%name)// &
           ' contains ',nerr,' non-converged values! Trying a correction!'
      CALL message ('WARNING: '//TRIM(routine), TRIM(message_text))

      !------------------------------------------------------------------------------
      ! .. 1st try: fill missing values by nearest neighbour values of the next
      !    smaller qi-value (here, Dmin varies the least):

      ALLOCATE(lfill(n_qw))
!$OMP parallel do private(i,j,k,l,lfill)
      DO j = 1, n_temp
        DO i = 1, n_pres

          ! .. First correction pass starting from l=1:
          DO l = 1, n_qi
!$NEC ivdep
            DO k = 1, n_qw
              IF (Dmin(k,l,i,j) < miss_thresh) THEN
                IF (l == 1) THEN
                  IF (Dmin(k,2,i,j) >= 0.0_wp) THEN
                    lfill(k) = 2
                  ELSE IF (Dmin(k,3,i,j) >= 0.0_wp) THEN
                    lfill(k) = 3
                  ELSE IF (Dmin(k,4,i,j) >= 0.0_wp) THEN
                    lfill(k) = 4
                  END IF
                ELSE
                  lfill(k) = l - 1
                END IF
                Dmin(k,l,i,j) = Dmin(k,lfill(k),i,j)
              END IF
            END DO
          END DO

          ! .. Second pass starting from l=n_qi in negative direction:
          DO l = n_qi, 1, -1
!$NEC ivdep
            DO k = 1, n_qw
              IF (Dmin(k,l,i,j) < miss_thresh) THEN
                IF (l == n_qi) THEN
                  IF (Dmin(k,n_qi-1,i,j) >= 0.0_wp) THEN
                    lfill(k) = n_qi-1
                  ELSE IF (Dmin(k,n_qi-2,i,j) >= 0.0_wp) THEN
                    lfill(k) = n_qi-2
                  ELSE IF (Dmin(k,n_qi-3,i,j) >= 0.0_wp) THEN
                    lfill(k) = n_qi-3
                  END IF
                ELSE
                  lfill(k) = l - 1
                END IF
                Dmin(k,l,i,j) = Dmin(k,lfill(k),i,j)
              END IF
            END DO
          END DO
        END DO
      END DO
!$OMP end parallel do
      DEALLOCATE(lfill)

      IF (ANY(Dmin < miss_thresh)) THEN

        nerr = COUNT((Dmin < miss_thresh))
        WRITE(message_text,'(a,i0,a)') 'first try to correction of wet-growth Dmin table for '//TRIM(phail%name)// &
             ' was not successful at ',nerr,' points! Doing hard correction!'
        CALL message ('WARNING: '//TRIM(routine), TRIM(message_text))

        !------------------------------------------------------------------------------
        ! .. Fill missing values by mean value over qi-direction where possible
        !    (here, Dmin varies the least), or by hard 999.99 if all else fails:

        ALLOCATE (rowmean(n_qw), rowcount(n_qw))
!$OMP parallel do private(i,j,k,l,rowmean,rowcount)
        DO j = 1, n_temp
          DO i = 1, n_pres
            DO l = 1, n_qi
              rowmean(:)  = SUM  (Dmin(:,:,i,j), dim=2, mask=(Dmin(:,:,i,j) > miss_thresh))
              rowcount(:) = COUNT(               dim=2, mask=(Dmin(:,:,i,j) > miss_thresh))
              WHERE (rowcount > 0)
                rowmean = rowmean / rowcount
              ELSEWHERE
                rowmean = miss_value
              END WHERE
              DO k = 1, n_qw
                IF (Dmin(k,l,i,j) < miss_thresh) THEN
                  Dmin(k,l,i,j) = rowmean(k)
                END IF
              END DO
            END DO
          END DO
        END DO
!$OMP end parallel do
        DEALLOCATE (rowmean, rowcount)

!$OMP parallel
        WHERE (Dmin < miss_thresh) Dmin = 999.99_wp
!$OMP end parallel

      ELSE

        WRITE(message_text,'(a,i0,a)') 'Correction for newly created wet-growth Dmin table for '//TRIM(phail%name)// &
             ' was successful!'
        CALL message ('INFO: '//TRIM(routine), TRIM(message_text))

      END IF

    END IF

    !=================================================================================
    ! .. Permute dimensions of Dmin_lut to Dmin_nc for output:

    ALLOCATE( Dmin_lut(n_pres,n_temp,n_qw,n_qi) ); Dmin_lut = 0._wp ! index order for final LUT

!$OMP parallel do private(i,j,k,l)
    DO j = 1, n_temp
      DO i = 1, n_pres
        DO l = 1, n_qi
!$NEC ivdep
          DO k = 1, n_qw
            Dmin_lut(i,j,k,l) = Dmin(k,l,i,j)
          END DO
        END DO
      END DO
    END DO
!$OMP end parallel do

    DEALLOCATE(Dmin)

  END SUBROUTINE dmin_wetgrowth_lookupcreate

  !**********************************************************************
  !**********************************************************************

  !**********************************************************************
  !
  ! Subroutine for writing a lookup table to a netcdf file
  !
  !**********************************************************************

  SUBROUTINE write_dmin_wetgrowth_table(basename, phail, Dmin_lut, pres, temp, qw, qi, ierr, filename)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: basename
    TYPE(PARTICLE),   INTENT(in) :: phail
    REAL(wp),         INTENT(in) :: Dmin_lut(:,:,:,:), pres(:), temp(:), qw(:), qi(:) ! in SI units!
    INTEGER,          INTENT(out):: ierr
    CHARACTER(len=*), INTENT(in), OPTIONAL :: filename

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modulename)//'::write_dmin_wetgrowth_table'

    CHARACTER(len=300) :: ncfile

    INTEGER            :: i, j, k, l
    INTEGER            :: n_pres, n_temp, n_qw, n_qi
    LOGICAL            :: fileexist
    INTEGER            :: cmode, fid, status
    INTEGER            :: pres_dimid, temp_dimid, qw_dimid, qi_dimid
    INTEGER            :: pres_varid, temp_varid, qw_varid, qi_varid, dmin_varid

    CHARACTER(8)  :: date
    CHARACTER(10) :: time
    CHARACTER(5)  :: zone

    ierr = 0

    n_pres = SIZE(pres)
    n_temp = SIZE(temp)
    n_qw   = SIZE(qw)
    n_qi   = SIZE(qi)

    ncfile(:) = ' '
    IF (PRESENT(filename)) THEN
      ncfile = TRIM(filename)
    ELSE
      ncfile = TRIM(get_filebase_dmin_lut_file(TRIM(basename), phail)) // '.nc'
    END IF

    ! Write table to file, if the file does not yet exist:

    ! Inquire, if the table file still does not exist. Maybe another Ensemble member
    ! has produced the table inbetween in case of ...
    INQUIRE(file=TRIM(ncfile), exist=fileexist)

    IF (.NOT. fileexist) THEN

      ! 1) Open and create Netcdf File
      ! ------------------------------

      message_text(:) = ' '
      WRITE(message_text,'(a)') 'writing NetCDF file '//TRIM(ncfile)//' ...'
      CALL message (TRIM(routine), TRIM(message_text))

      cmode = NF90_CLOBBER
!!$ UB: We go back to netcdf-3, because it is faster on input. For this, the special netcdf-4 flags are commented out:
!!$      cmode = IOR(cmode, NF90_NETCDF4)
!!$      cmode = IOR(cmode, NF90_CLASSIC_MODEL)  ! no fancy netcdf4 features (new data types, compounds, multiple unlimited dims, etc.), therefore smaller files

      fid = -1
      status = check_nc( nf90_create (TRIM(ncfile), cmode, fid), 'nf90_create '//TRIM(ncfile) )

      IF (status /= nf90_noerr) THEN

        message_text(:) = ' '
        WRITE(message_text,'(a)') 'ERROR when creating NetCDF lookup table file '//TRIM(ncfile)// &
             TRIM(nf90_strerror(status))
        fid = -1
        CALL message(TRIM(routine), TRIM(message_text))
        ierr = 1

      ELSE

        status = check_nc( nf90_def_dim(fid, 'npres',  n_pres, pres_dimid), 'nf90_def_dim pres')
        status = check_nc( nf90_def_dim(fid, 'ntemp',  n_temp, temp_dimid), 'nf90_def_dim temp')
        status = check_nc( nf90_def_dim(fid, 'nqw',  n_qw, qw_dimid), 'nf90_def_dim qw')
        status = check_nc( nf90_def_dim(fid, 'nqi',  n_qi, qi_dimid), 'nf90_def_dim qi')

        status = check_nc( nf90_def_var(fid, 'qw',    NF90_DOUBLE, (/qw_dimid/), qw_varid), 'nf90_def_var qw')
        status = check_nc( nf90_put_att(fid, qw_varid , 'units', 'g m-3')                 , 'nf90_put_att units for qw')
        status = check_nc( nf90_put_att(fid, qw_varid , 'long_name', 'Supercooled liquid water content'),  &
             'nf90_put_att long_name for qw')
        status = check_nc( nf90_put_att(fid, qw_varid , 'short_name', 'qw'),  'nf90_put_att short_name for qw')

        status = check_nc( nf90_def_var(fid, 'qi',    NF90_DOUBLE, (/qi_dimid/), qi_varid), 'nf90_def_var qi')
        status = check_nc( nf90_put_att(fid, qi_varid , 'units', 'g m-3')                 , 'nf90_put_att units for qi')
        status = check_nc( nf90_put_att(fid, qi_varid , 'long_name', 'Ice+snow water content'),  &
             'nf90_put_att long_name for T')
        status = check_nc( nf90_put_att(fid, qi_varid , 'short_name', 'qi'),  'nf90_put_att short_name for qi')

        status = check_nc( nf90_def_var(fid, 'T',    NF90_DOUBLE, (/temp_dimid/), temp_varid), 'nf90_def_var T')
        status = check_nc( nf90_put_att(fid, temp_varid , 'units', 'Celsius')               , 'nf90_put_att units for T')
        status = check_nc( nf90_put_att(fid, temp_varid , 'long_name', 'Temperature'),  &
             'nf90_put_att long_name for T')
        status = check_nc( nf90_put_att(fid, temp_varid , 'short_name', 'T'),  'nf90_put_att short_name for T')

        status = check_nc( nf90_def_var(fid, 'p',    NF90_DOUBLE, (/pres_dimid/), pres_varid), 'nf90_def_var p')
        status = check_nc( nf90_put_att(fid, pres_varid , 'units', 'hPa')                   , 'nf90_put_att units for p')
        status = check_nc( nf90_put_att(fid, pres_varid , 'long_name', 'Pressure'),  &
             'nf90_put_att long_name for p')
        status = check_nc( nf90_put_att(fid, pres_varid , 'short_name', 'p'),  'nf90_put_att short_name for p')

        status = check_nc( nf90_def_var(fid, 'Dmin_wetgrowth_table', NF90_DOUBLE, &
             (/pres_dimid, temp_dimid, qw_dimid, qi_dimid/), dmin_varid), 'nf90_def_var Dmin_wetgrowth_table')
        status = check_nc( nf90_put_att(fid, dmin_varid , 'units', 'mm'),  'nf90_put_att units for Dmin_wetgrowth_table')
        status = check_nc( nf90_put_att(fid, dmin_varid , 'long_name', 'Threshold diameter for wet-growth regime in mm'),  &
             'nf90_put_att long_name for Dmin_wetgrowth_table')
        status = check_nc( nf90_put_att(fid, dmin_varid , '_FillValue', miss_value),  &
             'nf90_put_att _FillValue for Dmin_wetgrowth_table')


        status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'title', 'Wetgrowth table for ICON 2mom scheme'), &
             'nf90_put_att title')
        status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'creator', 'ICON'), &
             'nf90_put_att creator')
        CALL date_and_TIME(date=date, time=time, zone=zone)
        status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'creationDateTime', date//time//' UCT'//zone), &
             'nf90_put_att creationDateTime')
        status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'ncdfFile', TRIM(ncfile)), &
             'nf90_put_att ncdf-file')
        status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'hydrometeorType', TRIM(phail%name)), &
             'nf90_put_att hydrometeorType')
        status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'a_geo', REAL(phail%a_geo)), &
             'nf90_put_att a_geo')
        status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'b_geo', REAL(phail%b_geo)), &
             'nf90_put_att b_geo')
        status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'a_vel', REAL(phail%a_vel)), &
             'nf90_put_att a_vel')
        status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'b_vel', REAL(phail%b_vel)), &
             'nf90_put_att b_vel')

        status = check_nc( nf90_enddef(fid) , 'nf90_enddef')

        ! 2) Write data to Netcdf File
        ! ----------------------------

        status = check_nc( nf90_put_var(fid, pres_varid, pres*1d-2, start=(/1/)), 'nf90_put_var pres')
        status = check_nc( nf90_put_var(fid, temp_varid, temp, start=(/1/)),     'nf90_put_var T'   )
        status = check_nc( nf90_put_var(fid, qw_varid, qw*1d3, start=(/1/)),     'nf90_put_var qw'  )
        status = check_nc( nf90_put_var(fid, qi_varid, qi*1d3, start=(/1/)),     'nf90_put_var qi'  )

        status = check_nc( nf90_put_var(fid, dmin_varid, Dmin_lut*1d3, start=(/1,1,1,1/)), &
             'nf90_put_var Dmin_wetgrowth_table')

        ! 3) Close Netcdf File
        ! --------------------

        status = check_nc( nf90_close(fid), 'nf90_close '//TRIM(ncfile) )
        IF (status /= nf90_noerr) THEN
          message_text(:) = ' '
          WRITE (message_text,'(a)') 'Error in closing file '//TRIM(ncfile)// &
               TRIM(nf90_strerror(status))
          CALL message(TRIM(routine), TRIM(message_text))
          ierr = 2
        END IF
      END IF

    ELSE

      WRITE(*,'(a)') 'INFO '//TRIM(routine)//': NetCDF file '//TRIM(ncfile)//' already exists! Not written!'

    END IF

  CONTAINS

    FUNCTION check_nc ( istat, rinfo ) RESULT (ostat)
      INTEGER, INTENT(in)          :: istat
      CHARACTER(len=*), INTENT(in) :: rinfo
      INTEGER                      :: ostat

      ! .. return the error status, so that the below error handling might be
      !     relaxed in the future (e.g., WARNING instead of abort_run())
      !     and the calling routine can decide what to do:
      ostat = istat
      IF (istat /= NF90_NOERR) THEN
        message_text(:) = ' '
        WRITE(message_text,'(a)') &
            'writing NetCDF: '//TRIM(rinfo)//': '//TRIM(NF90_strerror(istat))
        CALL message('ERROR '//TRIM(routine), TRIM(message_text))
        ierr = 3
      END IF

    END FUNCTION check_nc

  END SUBROUTINE write_dmin_wetgrowth_table
!
!
  FUNCTION get_filebase_dmin_lut_file(basename, phail) RESULT(filebase)

    CHARACTER(len=*), INTENT(in) :: basename
    TYPE(PARTICLE),   INTENT(in) :: phail

    CHARACTER(len=300)           :: filebase

    REAL(wp)           :: particle_parameters(4)
    INTEGER            :: i, j, ppar_int(4), ppar_dec(4)
    CHARACTER(len=100) :: cparams, ppar_c, ncint, ncdec

    particle_parameters(1) = phail%a_geo
    particle_parameters(2) = phail%b_geo
    particle_parameters(3) = phail%a_vel
    particle_parameters(4) = phail%b_vel

    ! .. Construct filename from basename and phail parameters. Make sure that
    !    the parameters are in decimal format with at least one digit before the komma and 4 decimals after the Komma.
    !    This hack is necessary, because different compilers implement different behaviour of Format "i0"
    ppar_int = INT(particle_parameters)
    ppar_dec = NINT(MODULO(particle_parameters, 1.0_wp) * 1e4)

    cparams(:) = ' '
    DO i=1, SIZE(particle_parameters)
      ncint(:) = ' '
      ncdec(:) = ' '
      IF (ppar_int(i) >= 0) THEN
        j = MAX(ppar_int(i), 1)
      ELSE
        j = MAX(ABS(ppar_int(i)), 1) + 1  ! 1 add char for minus sign
      END IF
      WRITE(ncint,'(i3)') MAX(INT(LOG10(REAL(j))+1), 1)
      IF (ppar_dec(i) >= 0) THEN
        j = MAX(ppar_dec(i), 1)
      ELSE
        j = MAX(ABS(ppar_dec(i)), 1) + 1  ! 1 add char for minus sign
      END IF
      WRITE(ncdec,'(i3)') MAX(INT(LOG10(REAL(j))+1), 1)
      ppar_c(:) = ' '
      WRITE (ppar_c,'("_",i'//TRIM(ADJUSTL(ncint))//',".",i'//TRIM(ADJUSTL(ncdec))//')') ppar_int(i), ppar_dec(i)
      cparams = TRIM(cparams)//TRIM(ADJUSTL(ppar_c))
    END DO

    filebase = TRIM(basename) // TRIM(cparams)

  END FUNCTION get_filebase_dmin_lut_file

  !
!**********************************************************************
!**********************************************************************
!
  REAL(wp) FUNCTION dwg_fpi(ageo_x,bgeo_x,avel_x,bvel_x,D,T,p,qw,qi)
!
! Implizite Funktion fuer den Gleichgewichtsdurchmesser bei Ts = 0 Grad C eines fallenden
! (trockenen) Hagelkorns unter der (konstanten) Wirkung von Riming, Akkretion von
! Wolken- und Regenwasser und Eis-Sublimation. Bezueglich letzterem
! wird angenommen, dass sich die Wolkenluft relativ zu Wasser in
! Saettigung befindet. Der Gleichgewichtsdurchmesser gibt die Trennung
! zwischen Dry- und Wet Growth an.

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: ageo_x, bgeo_x, avel_x, bvel_x, D, T, p, qw, qi

    REAL(wp) :: T_K, Cwat, Ci, Ci2, Dv, rho_l, vf, fv, K, tmp_dwg

    T_K = T + T_0C

    Cwat = Cwater(T_K)
    Ci = Cice(T_K)
    Ci2 = Cice(T_0C)
    Dv = diff_vap_air(T_K,p);

    rho_l = p / (R_a * T_K * (1._wp+0.61_wp*(0.622_wp*esat_w(T_K)/p)))

    vf = vfall_graupel(ageo_x,bgeo_x,avel_x,bvel_x,D,T,p)

    fv = vent_v(D,rho_l,T_K,vf)       ! Rasmussen & Pruppacher

    K = conduct_air(T_K)

    tmp_dwg = (2._wp*pi*fv*(Lh_s(T_0C)*Dv* &
         (esat_i(T_0C)/T_0C - esat_w(T_K)/T_K)/R_w - K*T) ) / &
         (pi/4._wp*vf*(e_collw(T_K)*qw*(Lh_m(T_K)+Cwat*T)+e_colli(T_0C)*qi*Ci*T))

    dwg_fpi = tmp_dwg

    RETURN

  END FUNCTION dwg_fpi
!
!**********************************************************************
!**********************************************************************
!
! Ventilation factor for vapor diffusion / heat transfer from a
! water droplet, see Pruppacher and Klett, 1997, p. 541 of Beard
! and Pruppacher (1971):
!
  REAL(wp) FUNCTION vent_v(D,rhol,T,v)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: D, rhol, & ! D in m, rhol in kg/m3
                            T, v       ! T in K and v in m/s

    REAL(wp), PARAMETER :: Nsc = 0.71_wp  ! for water vapor in air; depends on the fluid mixture!!!
    ! Note that the value of 0.71 is wrong, because back in 1971,
    ! they used wrong values of Dv as function of T, and n_sc = kinemat. visc.(T) / Dv(T).
    ! However, since in the experiments leading to the below formula,
    ! fv has been measured and related to computed values of fakt,
    ! one has to use 0.71 to retrieve the correct values for fv!

    REAL(wp), PARAMETER :: chi_lim = 1.4_wp
    REAL(wp) :: nu_l, Re, chi, tmp

    nu_l = dyn_visc_air(T) / rhol
    Re = v * D / nu_l

    chi = Nsc**0.33333_wp * SQRT(Re)

    IF (chi < chi_lim) THEN
       tmp = 1._wp + 0.108_wp * chi**2
    ELSE
       tmp = 0.78_wp + 0.308_wp * chi
    END IF

    vent_v = tmp

    RETURN

  END FUNCTION vent_v
!
!**********************************************************************
!**********************************************************************
!
! Sticking efficiency nach Lin et al. (1983) fuer die Kollision
! von Hagel mit Eis- und Schneepartikeln:
!
  REAL(wp) FUNCTION e_colli(T)
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: T

    e_colli = MIN(EXP(0.09_wp*(T-T_0C)),1._wp)

    RETURN
  END FUNCTION e_colli
!
!**********************************************************************
!**********************************************************************
!
! Riming eff. of cloud and rain water on graupel/hail, assumed as 1
!
  REAL(wp) FUNCTION e_collw(T)
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: T

    e_collw = 1._wp

    RETURN
  END FUNCTION e_collw
!
!**********************************************************************
!**********************************************************************
!
! Fall velocity of graupel/hail
!
  REAL(wp) FUNCTION vfall_graupel(ageo_x,bgeo_x,avel_x,bvel_x,D,T,p)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: ageo_x, bgeo_x, avel_x, bvel_x, D, T, p
    REAL(wp) :: rho, T_K, vel, x

    T_K = T + T_0C

    rho = p / (R_a * T_K * (1_wp+0.61_wp*(0.622_wp*esat_w(T_K)/p)))

    x = (D/ageo_x)**(1._wp/bgeo_x)

    vel = avel_x * x**bvel_x * (1.21_wp/rho)**0.5

    vfall_graupel = vel

    RETURN

  END FUNCTION vfall_graupel
!
!**********************************************************************
!**********************************************************************

!**********************************************************************
!
!   THERMODYNAMIC FUNCTIONS:
!
!   DESCRIPTION: All relevant thermodynamic functions
!                for liquid and ice water.
!
!   SOURCE:      From MATLAB codes by U. Blahak
!
!**********************************************************************
!
!
!**********************************************************************
!
!  DESCRIPTION:
!       Specific Heat capacity C of ice as function of temperature,
!       after Landolt-Boernstein, New Series V/4b (Laube, Hoeller)
!
!  Note: Cice in J/kg/K
!
  REAL(wp) FUNCTION Cice(T)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: T     ! T in K

    INTEGER,  PARAMETER :: nt = 11
    REAL(wp), PARAMETER :: temp0(nt) = (/ &
         173.16_wp, 183.16_wp, 193.16_wp, 203.16_wp, &
         213.16_wp, 223.16_wp, 233.16_wp, 243.16_wp, &
         253.16_wp, 263.16_wp, 273.16_wp /)
    REAL(wp), PARAMETER :: C0(nt) = (/ &
         1382._wp, 1449._wp, 1520._wp, 1591._wp, &
         1662._wp, 1738._wp, 1813._wp, 1884._wp, &
         1959._wp, 2031._wp, 2106._wp /)
    REAL(wp), PARAMETER :: dtemp = temp0(2) - temp0(1)

    INTEGER  :: jt
    REAL(wp) :: temp, C

    temp = MAX(MIN(T,temp0(nt)),temp0(1))

!    CALL locate(temp0,nt,temp,jt)
    jt = MAX(MIN(FLOOR((temp-temp0(1))/dtemp) + 1, nt-1), 1)

    Cice = linint(C0(jt),C0(jt+1),temp0(jt),temp0(jt+1),temp)

    RETURN

  END FUNCTION Cice
!
!**********************************************************************
!**********************************************************************
!
!  DESCRIPTION:
!       Specific Heat capacity C of water as function of temperature,
!       after Landolt-Boernstein, New Series V/4b (Laube, Hoeller)
!
!  Note: Cwater in J/kg/K
!
  REAL(wp) FUNCTION Cwater(T)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: T     ! T in K

    INTEGER,  PARAMETER :: nt = 11
    REAL(wp), PARAMETER :: temp0(nt) = (/ &
         223.16_wp, 233.16_wp, 243.16_wp, 253.16_wp, &
         263.16_wp, 273.16_wp, 283.16_wp, 293.16_wp, &
         303.16_wp, 313.16_wp, 323.16_wp /)
    REAL(wp), PARAMETER :: C0(nt) = (/ &
         5400._wp, 4770._wp, 4520._wp, 4350._wp, &
         4270._wp, 4217.8_wp, 4192.3_wp, 4181.8_wp, &
         4178.5_wp, 4178.5_wp, 4180.6_wp /)
    REAL(wp), PARAMETER :: dtemp = temp0(2) - temp0(1)

    INTEGER  :: jt
    REAL(wp) :: temp, C

    temp = MAX(MIN(T,temp0(nt)),temp0(1))

    jt = MAX(MIN(FLOOR((temp-temp0(1))/dtemp) + 1, nt-1), 1)

    Cwater = linint(C0(jt),C0(jt+1),temp0(jt),temp0(jt+1),temp)

    RETURN

  END FUNCTION Cwater
!
!**********************************************************************
!**********************************************************************
!
!  DESCRIPTION:
!       Calculates the diffusivity of water vapour in air as a fct. of
!       temperature and pressure in m2/s
!       Following Montgomery (1947)
!
  REAL(wp) FUNCTION diff_vap_air(T,p)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: T, p   ! Temp. in K and pressure in Pa

    Diff_vap_air = 2.26e-5_wp * (T/T_0C)**1.81 * (1000.e2_wp/p)

  END FUNCTION diff_vap_air
!
!**********************************************************************
!**********************************************************************
!
!  DESCRIPTION:
!        Berechnet die dynamische Zaehigkeit eta von Luft in [Pa s]
!        nach der Sutherland-Formel (Sutherland 1893, Landolt-Boernstein)
!        als Funktion der Temperatur T in K. Gueltigkeitsbereich nach
!        Landolt-Boernstein: T element aus [-70.0, 40.0] grad C
!
  REAL(wp) FUNCTION dyn_visc_air(T)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: T     ! Temp. in K
    REAL(wp) :: temp

    temp = MAX(MIN(T,313.16_wp),203.16_wp)

    dyn_visc_air = 1.8325e-5_wp * (temp/296.16_wp)**1.5 * &
         (296.16_wp+120._wp)/(temp+120._wp)

  END FUNCTION dyn_visc_air
!
!**********************************************************************
!**********************************************************************
!
!  DESCRIPTION:
!        Saturation vapour pressure in Pa over ice as a funct. of
!        temperature in Kelvin. As in LM.
!
  REAL(wp) FUNCTION esat_i(T)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: T     ! Temp. in K
    REAL(wp), PARAMETER  :: Tlim = 7.66_wp + 1.e-6_wp
    REAL(wp) :: tmp

    IF (T > Tlim) THEN
       tmp = 610.78_wp * EXP(21.8745584_wp*(T-T_0C)/(T-7.66_wp))
    ELSE
       tmp = 0._wp
    END IF

    esat_i = tmp

    RETURN

  END FUNCTION esat_i
!
!**********************************************************************
!**********************************************************************
!
!  DESCRIPTION:
!        Saturation vapour pressure in Pa over water as a funct. of
!        temperature in Kelvin. As in LM.
!
  REAL(wp) FUNCTION esat_w(T)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: T     ! Temp. in K
    REAL(wp), PARAMETER  :: Tlim = 35.86_wp + 1.e-6_wp
    REAL(wp) :: tmp

    IF (T > Tlim) THEN
       tmp = 610.78_wp * EXP(17.2693882_wp*(T-T_0C)/(T-35.86_wp))
    ELSE
       tmp = 0._wp
    END IF

    esat_w = tmp

    RETURN

  END FUNCTION esat_w
!
!**********************************************************************
!**********************************************************************
!
!  DESCRIPTION:
!        Saturation vapour pressure in Pa over water as a funct. of
!        temperature in Kelvin. Following Bolton.
!
  REAL(wp) FUNCTION esat_w_bolton(T)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: T     ! Temp. in K
    REAL(wp), PARAMETER  :: Tlim = 29.66_wp + 1.e-6_wp
    REAL(wp) :: tmp

    IF (T > Tlim) THEN
       tmp = 611.2_wp * EXP(17.67_wp*(T-T_0C)/(T-29.66_wp))
    ELSE
       tmp = 0._wp
    END IF

    esat_w_bolton = tmp

    RETURN

  END FUNCTION esat_w_bolton
!
!**********************************************************************
!**********************************************************************
!
!  DESCRIPTION:
!       Latent heat of melting l of water as function of temperature,
!       after Landolt-Boernstein, New Series V/4b (Laube, Hoeller)
!
!  Note: Lh_m in J/kg
!
  REAL(wp) FUNCTION Lh_m(T)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: T     ! T in K

    INTEGER,  PARAMETER :: nt = 6
    REAL(wp), PARAMETER :: temp0(nt) = (/ &
         223.16_wp, 233.16_wp, 243.16_wp, &
         253.16_wp, 263.16_wp, 273.16_wp /)
    REAL(wp), PARAMETER :: L0(nt) = 1.e6 * (/ &
         0.2035_wp, 0.2357_wp, 0.2638_wp, &
         0.2889_wp, 0.3119_wp, 0.3337_wp /)
    REAL(wp), PARAMETER :: dtemp = temp0(2) - temp0(1)

    INTEGER  :: jt
    REAL(wp) :: temp

    temp = MAX(MIN(T,temp0(nt)),temp0(1))

    jt = MAX(MIN(FLOOR((temp-temp0(1))/dtemp) + 1, nt-1), 1)

    Lh_m = linint(L0(jt),L0(jt+1),temp0(jt),temp0(jt+1),temp)

    RETURN

  END FUNCTION Lh_m
!
!**********************************************************************
!**********************************************************************
!
!  DESCRIPTION:
!       Latent heat of sublimation l of water as function of temperature,
!       after Landolt-Boernstein, New Series V/4b (Laube, Hoeller)
!
!  Note: Lh_s in J/kg
!
  REAL(wp) FUNCTION Lh_s(T)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: T     ! T in K

    INTEGER,  PARAMETER :: nt = 11
    REAL(wp), PARAMETER :: temp0(nt) = (/ &
         173.16_wp, 183.16_wp, 193.16_wp, 203.16_wp, &
         213.16_wp, 223.16_wp, 233.16_wp, 243.16_wp, &
         253.16_wp, 263.16_wp, 273.16_wp /)
    REAL(wp), PARAMETER :: L0(nt) = 1.e6 * (/ &
         2.8236_wp, 2.8278_wp, 2.8316_wp, 2.8345_wp, &
         2.8366_wp, 2.8383_wp, 2.8387_wp, 2.8387_wp, &
         2.8383_wp, 2.8366_wp, 2.8345_wp /)
    REAL(wp), PARAMETER :: dtemp = temp0(2) - temp0(1)

    INTEGER  :: jt
    REAL(wp) :: temp

    temp = MAX(MIN(T,temp0(nt)),temp0(1))

    jt = MAX(MIN(FLOOR((temp-temp0(1))/dtemp) + 1, nt-1), 1)

    Lh_s = linint(L0(jt),L0(jt+1),temp0(jt),temp0(jt+1),temp)

    RETURN

  END FUNCTION Lh_s
!
!**********************************************************************
!**********************************************************************
!
!  DESCRIPTION:
!       Latent heat of evaporation l of water as function of temperature,
!       after Landolt-Boernstein, New Series V/4b (Laube, Hoeller)
!
!  Note: Lh_e in J/kg
!
  REAL(wp) FUNCTION Lh_e(T)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: T     ! T in K

    INTEGER,  PARAMETER :: nt = 11
    REAL(wp), PARAMETER :: temp0(nt) = (/ &
         223.16_wp, 233.16_wp, 243.16_wp, 253.16_wp, &
         263.16_wp, 273.16_wp, 283.16_wp, 293.16_wp, &
         303.16_wp, 313.16_wp, 323.16_wp /)
    REAL(wp), PARAMETER :: L0(nt) = 1.e6 * (/ &
         2.6348_wp, 2.6030_wp, 2.5749_wp, 2.5494_wp, &
         2.5247_wp, 2.50084_wp, 2.4774_wp, 2.4535_wp, &
         2.4300_wp, 2.4062_wp, 2.3823_wp /)
    REAL(wp), PARAMETER :: dtemp = temp0(2) - temp0(1)

    INTEGER  :: jt
    REAL(wp) :: temp

    temp = MAX(MIN(T,temp0(nt)),temp0(1))

    jt = MAX(MIN(FLOOR((temp-temp0(1))/dtemp) + 1, nt-1), 1)

    Lh_e = linint(L0(jt),L0(jt+1),temp0(jt),temp0(jt+1),temp)

    ! Naeherung von Bolton:
!!$    Lh_e = (2.501_wp - 0.00237_wp*(T-T_0C) ) * 1.e6_wp

    RETURN

  END FUNCTION Lh_e
!
!**********************************************************************
!**********************************************************************
!
!  DESCRIPTION:
!     Berechnet die Waermeleitfaehigkeit von Luft in [W/m/K] aus dem
!     Wert bei 0 Grad C und derselben Temperaturabhaenigkeit
!     wie bei der dyn. viskositaet. Verwendet wird die
!     Sutherland-Formel (Sutherland 1893, Landolt-Boernstein)
!     als Funktion der Temperatur T in K. Siehe auch Montgomery (1947),
!     speziell fuer den Referenzwert bei 0 Grad C.
!     Liefert sehr aehnliche Werte wie die einfache lineare Regression aus
!     Landolt/Boernstein, so wie sie in der Funktion heatconduct_air(T)
!     programmiert ist.
!
!  Note: conduct_air in W/m/K
!
  REAL(wp) FUNCTION conduct_air(T)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: T    ! Temp. in K

    conduct_air = 0.0242_wp*(T/T_0C)**1.5 * &
         (T_0C+120._wp)/(T+120._wp)

    RETURN

  END FUNCTION conduct_air
!
!**********************************************************************
!**********************************************************************
!
!  DESCRIPTION:
!     Waermeleitfaehigkeit von Wasser in W/m/K als Funktion der Temperatur
!     T in K, Fit gueltig fuer 0 - 60 Grad C
!     Unbekannte Herkunft!!!
!
!  Note: conduct_water in W/m/K
!
  REAL(wp) FUNCTION conduct_water(T)

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: T    ! Temp. in K
    REAL(wp) :: temp

    temp = MIN(MAX(T-273.16_wp,0._wp),60._wp)

    conduct_water = 0.56905_wp + 0.0019025_wp*temp - &
         0.000008125_wp*(temp)**2

    RETURN

  END FUNCTION conduct_water
!
!**********************************************************************
!**********************************************************************
!
  !
  !**********************************************************************
  !
  !    DESCRIPTION:
  !       Locate value x in arrax xx and return index j, such that
  !       x is located between xx(j) and xx(j+1).
  !       IMPORTANT: xx has to be monotonic (increasing or decreasing).
  !
  !    METHOD: Bi-section
  !
  !**********************************************************************
  !
  SUBROUTINE locate(xx,n,x,j)

    IMPLICIT NONE

    INTEGER,  INTENT(in) :: n
    REAL(wp), INTENT(in) :: xx(n), x
    INTEGER, INTENT(out) :: j

    INTEGER              :: jl, jm, ju
    REAL(wp)             :: xi

    jl=0
    ju=n+1

    IF (xx(1) <= xx(n)) THEN
      ! monotonically increasing
      xi = MIN(MAX(x,xx(1)),xx(n))
    ELSE
      ! monotonically decreasing
      xi = MAX(MIN(x,xx(1)),xx(n))
    END IF

    searchloop: DO
      jm = (ju+jl) / 2
      IF ((xx(n) >= xx(1)) .EQV. (x >= xx(jm))) THEN
        jl=jm
      ELSE
        ju=jm
      END IF
      IF (ju-jl <= 1) EXIT searchloop
    END DO searchloop

    IF (x == xx(1)) THEN
       j = 1
    ELSE IF (x == xx(n)) THEN
       j = n-1
    ELSE
       j = jl
    END IF

    RETURN

  END SUBROUTINE locate

  ! ********************************************************************
  !
  ! Linear interpolation between two points (x1,y1) and (x2,y2)
  ! for a given y.
  !
  ! ********************************************************************

  REAL(wp) FUNCTION linint(x1,x2,y1,y2,y)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: x1, x2, y1, y2, y
    REAL(wp)             :: x
    REAL(wp), PARAMETER  :: eps = 1e-20_wp

    IF (ABS(y2-y1) < eps) THEN
      linint = x1
!      WRITE (*, '(a)') 'Error in thermo::linint(): table nodes y1 and y2 are equal!'
    ELSE
      !..interpolated value
      linint = x1 + (x2-x1)/(y2-y1) * (y-y1)
    END IF

    RETURN

  END FUNCTION linint
!
!**********************************************************************
!**********************************************************************
!

END MODULE mo_2mom_mcrph_dmin_wetgrowth
