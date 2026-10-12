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

! Read and apply monthly aerosol optical properties of S. Kinne
! from yearly files.

! ---------------------------
#include "consistent_fma.inc"
! ---------------------------

MODULE mo_bc_aeropt_splumes_opt

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish
  USE mo_read_interface,       ONLY: openInputFile, read_1D, &
                                   & read_bcast_real_2D, read_bcast_real_3D, &
                                   & closeFile
  USE mo_model_domain,         ONLY: p_patch
  USE mo_math_constants,       ONLY: rad2deg
  USE mtime,                   ONLY: datetime, getDayOfYearFromDateTime, &
       &                             getNoOfSecondsElapsedInDayDateTime, &
       &                             getNoOfDaysInYearDateTime
  USE mo_mpi,                  ONLY: p_io, p_comm_work, p_bcast
  USE mo_bc_aeropt_splumes_memory, ONLY: bc_spl_field, &
       &                                 nplumes, nfeatures, ntimes, nyears

  IMPLICIT NONE

  PRIVATE
  PUBLIC                  :: add_bc_aeropt_splumes_opt, setup_bc_aeropt_splumes_opt

  REAL(wp), ALLOCATABLE, TARGET ::                    &
       plume_lat   (:)  ,& !< (nplumes) latitude where plume maximizes
       plume_lon   (:)  ,& !< (nplumes) longitude where plume maximizes
       beta_a      (:)  ,& !< (nplumes) parameter a for beta function
                           !< vertical profile
       beta_b      (:)  ,& !< (nplumes) parameter b for beta function
                           !< vertical profile
       aod_spmx    (:)  ,& !< (nplumes) aod at 550 for simple plume (maximum)
       aod_fmbg    (:)  ,& !< (nplumes) aod at 550 for fine mode
                           !< natural background (for twomey effect)
       asy550      (:)  ,& !< (nplumes) asymmetry parameter for plume at 550nm
       ssa550      (:)  ,& !< (nplumes) single scattering albedo for
                           !< plume at 550nm
       angstrom    (:)  ,& !< (nplumes) angstrom parameter for plume
       sig_lon_E   (:,:),& !< (nfeatures,nplumes) Eastward extent of
                           !< plume feature
       sig_lon_W   (:,:),& !< (nfeatures,nplumes) Westward extent of
                           !< plume feature
       sig_lat_E   (:,:),& !< (nfeatures,nplumes) Southward extent of
                           !< plume feature
       sig_lat_W   (:,:),& !< (nfeatures,nplumes) Northward extent of
                           !< plume feature
       theta       (:,:),& !< (nfeatures,nplumes) Rotation angle of feature
       ftr_weight  (:,:),& !< (nfeatures,nplumes) Feature weights =
                           !< (nfeatures + 1) to account for BB background
       year_weight (:,:)    ,& !< (nyear,nplumes) Yearly weight for plume
       ann_cycle   (:,:,:)     !< (nfeatures,ntimes,nplumes) annual cycle for feature
  REAL(wp)                 :: &
       time_weight (nfeatures,nplumes), &    !< Time-weights to account for BB background
     & time_weight_bg (nfeatures,nplumes)    !< as time_wight but for natural background in Twomey effect

  CHARACTER(LEN=256)       :: cfname

  CONTAINS

  ! -----------------------------------------------------------------
  ! SETUP_BC_AEROPT_SPLUMES:  This subroutine should be called at initialization to
  !            read the netcdf data that describes the simple plume
  !            climatology.  The information needs to be either read
  !            by each processor or distributed to processors.
  !
  SUBROUTINE setup_bc_aeropt_splumes_opt
    !
    ! ----------
    !
    INTEGER           :: ifile_id

      cfname='MACv2.0-SP_v1.nc'
    CALL openInputFile(ifile_id, cfname)

    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='plume_lat',&
                       & alloc_array=plume_lat,    file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='plume_lon',&
                       & alloc_array=plume_lon,    file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='beta_a',   &
                       & alloc_array=beta_a,       file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='beta_b',   &
                       & alloc_array=beta_b,       file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='aod_spmx', &
                       & alloc_array=aod_spmx,     file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='aod_fmbg', &
                       & alloc_array=aod_fmbg,     file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='ssa550',   &
                       & alloc_array=ssa550,       file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='asy550',   &
                       & alloc_array=asy550,       file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='angstrom', &
                       & alloc_array=angstrom,     file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lat_W',&
                       & alloc_array=sig_lat_W,    file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lat_E',&
                       & alloc_array=sig_lat_E,    file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lon_W',&
                       & alloc_array=sig_lon_W,    file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lon_E',&
                       & alloc_array=sig_lon_E, file_name=cfname,            &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='theta',    &
                       & alloc_array=theta,        file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes_opt',             &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='ftr_weight',&
                       & alloc_array=ftr_weight,   file_name=cfname,          &
                       & variable_dimls=(/nfeatures,nplumes/),                &
                       & module_name='mo_bc_aeropt_splumes_opt',              &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='year_weight',&
                       & alloc_array=year_weight, file_name=cfname,            &
                       & variable_dimls=(/nyears,nplumes/),                    &
                       & module_name='mo_bc_aeropt_splumes_opt',               &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL read_3d_wrapper(ifile_id=ifile_id,        variable_name='ann_cycle', &
                       & alloc_array=ann_cycle,    file_name=cfname,          &
                       & variable_dimls=(/nfeatures,ntimes,nplumes/),         &
                       & module_name='mo_bc_aeropt_splumes_opt',              &
                       & sub_prog_name='setup_bc_aeropt_splumes_opt'         )

    CALL closeFile(ifile_id)

    !$ACC ENTER DATA COPYIN(beta_a, beta_b, plume_lat, plume_lon, sig_lon_W, sig_lon_E, sig_lat_W, sig_lat_E, ftr_weight) &
    !$ACC   COPYIN(aod_spmx, aod_fmbg, ssa550, asy550, year_weight, ann_cycle, angstrom, theta) &
    !$ACC   CREATE(time_weight, time_weight_bg) ASYNC(1)
    !Host has these pointer already and since we just read them, copyin is correct.

    RETURN
  END SUBROUTINE SETUP_BC_AEROPT_SPLUMES_OPT
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SET_TIME_WEIGHT:  The simple plume model assumes that meteorology constrains plume shape and that only source strength
  ! influences the amplitude of a plume associated with a given source region.   This routine retrieves the temporal weights
  ! for the plumes.  Each plume feature has its own temporal weights which varies yearly.  The annual cycle is indexed by
  ! week in the year and superimposed on the yearly mean value of the weight.
  !
  SUBROUTINE set_time_weight(year_fr, opt_use_acc)
    !
    ! ----------
    !
    REAL(wp), INTENT(IN) ::  &
         year_fr           !< Fractional Year (1850.0 - 2100.99)

    INTEGER          ::  &
         iyear          ,& !< Integer year values between 1 and 156 (1850-2100)
         iweek          ,& !< Integer index (between 1 and ntimes); for ntimes=52 this corresponds to weeks (roughly)
         iplume            ! plume number

    LOGICAL, INTENT(in), OPTIONAL   :: opt_use_acc
    LOGICAL                         :: use_acc
    IF (PRESENT(opt_use_acc)) use_acc = opt_use_acc
    !
    ! ----------
    !
    iyear = FLOOR(year_fr) - 1849
    iweek = FLOOR((year_fr - FLOOR(year_fr)) * ntimes) + 1

    IF ((iweek > ntimes) .OR. (iweek < 1) .OR. (iyear > nyears) .OR. (iyear < 1)) STOP 'Time out of bounds in set_time_weight'
    !$ACC PARALLEL LOOP GANG(STATIC: 1) VECTOR DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
    DO iplume=1,nplumes
      time_weight(1,iplume) = year_weight(iyear,iplume) * ann_cycle(1,iweek,iplume)
      time_weight(2,iplume) = year_weight(iyear,iplume) * ann_cycle(2,iweek,iplume)
      time_weight_bg(1,iplume) = ann_cycle(1,iweek,iplume)
      time_weight_bg(2,iplume) = ann_cycle(2,iweek,iplume)
    END DO
    !$ACC END PARALLEL LOOP
    RETURN
  END SUBROUTINE set_time_weight
  !
  ! ---------------------------------------------------------------------------------------------
  ! SP_AOP_PROFILE:  This subroutine calculates the simple plume aerosol and cloud active optical
  !                  properites based on the the simple plume fit to the MPI Aerosol Climatology
  !                  (Version 2).  It sums over nplumes to provide a profile of aerosol
  !                  optical properties on a host models vertical grid.
  !
  SUBROUTINE sp_aop_profile(           nlevels        ,&
     & jcs            ,ncol           ,ncol_max       ,lambda         ,oro            , &
     & year_fr        ,f1             ,f2             ,f3             ,f4             , &
     & z              ,dz             ,dNovrN         ,aod_prof       ,ssa_prof       ,asy_prof       , &
     & opt_use_acc                                                                                      )
    !
    ! ----------
    !
    INTEGER, INTENT(IN)        :: &
       & nlevels,                 & !< number of levels
       & jcs,                     & !< start index in block
       & ncol,                    & !< number of columns (end index)
       & ncol_max                   !< first dimension of 2d-vars as declared in calling (sub)program [nproma]

    REAL(wp), INTENT(IN)       :: &
       & lambda,                  & !< wavelength
       & year_fr,                 & !< Fractional Year (1903.0 is the 0Z on the first of January 1903, Gregorian)
       & oro(ncol),               & !< orographic height (m)
       & z (ncol_max,nlevels),    & !< height above sea-level (m)
       & dz(ncol_max,nlevels),    & !< level thickness (difference between half levels)
       & f1(ncol_max,nplumes),    &
       & f2(ncol_max,nplumes),    &
       & f3(ncol_max,nplumes),    &
       & f4(ncol_max,nplumes)

    REAL(wp), INTENT(OUT)      ::     &
       & dNovrN(ncol)               , & !< anthropogenic increment to cloud drop number concentration
       & aod_prof(ncol_max,nlevels) , & !< profile of aerosol optical depth
       & ssa_prof(ncol_max,nlevels) , & !< profile of single scattering albedo
       & asy_prof(ncol_max,nlevels)     !< profile of asymmetry parameter

    INTEGER                    :: iplume, icol, k

    REAL(wp)                   ::  &
       & eta(ncol_max,nlevels),    & !< normalized height (by 15 km)
       & z_beta(ncol_max,nlevels), & !< profile for scaling column optical depth
       & prof(ncol_max,nlevels),   & !< scaled profile (by beta function)
       & beta_sum(ncol),           & !< vertical sum of beta function
       & ssa(ncol),                & !< aerosol optical depth
       & asy(ncol),                & !< aerosol optical depth
       & cw_an(ncol),              & !< column weight for simple plume (anthropogenic) aod at 550 nm
       & cw_bg(ncol),              & !< column weight for fine-mode indurstrial background aod at 550 nm
       & caod_sp(ncol),            & !< column simple plume (anthropogenic) aod at 550 nm
       & caod_bg(ncol),            & !< column fine-mode natural background aod at 550 nm
       & f1s,                      & !< contribution from feature 1
       & f2s,                      & !< contribution from feature 2
       & f3s,                      & !< contribution from feature 1 in natural background of Twomey effect
       & f4s,                      & !< contribution from feature 2 in natural background of Twomey effect
       & aod_550,                  & !< aerosol optical depth at 550nm
       & aod_lmd,                  & !< aerosol optical depth at input wavelength
       & lfactor                     !< factor to compute wavelength dependence of optical properties


    LOGICAL, INTENT(IN), OPTIONAL               :: opt_use_acc
    LOGICAL                                     :: use_acc
    IF (PRESENT(opt_use_acc)) use_acc = opt_use_acc

    !$ACC DATA CREATE(eta, z_beta, prof, beta_sum, ssa, asy, cw_an, cw_bg, caod_sp, caod_bg) ASYNC(1) IF(use_acc)
    !
    ! ----------
    !
    ! get time weights
    CALL set_time_weight(year_fr, opt_use_acc=use_acc)
    !
    !
    ! initialize variables, including output
    !

    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
    DO k=1,nlevels
      DO icol=jcs,ncol
        aod_prof(icol,k) = 0.0_wp
        ssa_prof(icol,k) = 0.0_wp
        asy_prof(icol,k) = 0.0_wp
        z_beta(icol,k)   = MERGE(1.0_wp, 0.0_wp, z(icol,k) >= oro(icol))
        eta(icol,k)      = MAX(0.0_wp,MIN(1.0_wp,z(icol,k)/15000._wp))
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
    DO icol=jcs,ncol
      dNovrN(icol)   = 1.0_wp
      caod_sp(icol)  = 0.00_wp
      caod_bg(icol)  = 0.02_wp
    END DO
    !$ACC END PARALLEL LOOP
    !
    ! sum contribution from plumes to construct composite profiles of aerosol otpical properties
    !
    DO iplume=1,nplumes
      !
      ! calculate vertical distribution function from parameters of beta distribution
      !
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
      DO icol=jcs,ncol
        beta_sum(icol) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
      DO icol=jcs,ncol
        !$ACC LOOP SEQ
        DO k=1,nlevels
          prof(icol,k)   = (eta(icol,k)**(beta_a(iplume)-1._wp) * (1._wp-eta(icol,k))**(beta_b(iplume)-1._wp))*dz(icol,k)
          beta_sum(icol) = beta_sum(icol) + prof(icol,k)
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
      DO k=1,nlevels
        DO icol=jcs,ncol
          prof(icol,k)   = prof(icol,k) / beta_sum(icol) * z_beta(icol,k)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !
      ! calculate plume weights
      !
      !PREVENT_INCONSISTENT_IFORT_FMA
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
      DO icol=jcs,ncol
        ! calculate contribution to plume from its different features, to get a column weight for the anthropogenic
        ! (cw_an) and the fine-mode background aerosol (cw_bg)
        !
        f1s = time_weight(1,iplume) * f1(icol,iplume)
        f2s = time_weight(2,iplume) * f2(icol,iplume)
        f3s = time_weight_bg(1,iplume) * f3(icol,iplume)
        f4s = time_weight_bg(2,iplume) * f4(icol,iplume)


        cw_an(icol) = f1s * aod_spmx(iplume) + f2s * aod_spmx(iplume)
        cw_bg(icol) = f3s * aod_fmbg(iplume) + f4s * aod_fmbg(iplume)
        !
        ! calculate wavelength-dependent scattering properties
        !
        lfactor   = MIN(1.0_wp,700.0_wp/lambda)
        ssa(icol) = (ssa550(iplume) * lfactor**4) / ((ssa550(iplume) * lfactor**4) + ((1-ssa550(iplume)) * lfactor))
        asy(icol) =  asy550(iplume) * SQRT(lfactor)
      END DO
      !$ACC END PARALLEL LOOP
      !
      ! distribute plume optical properties across its vertical profile weighting by optical depth and scaling for
      ! wavelength using the anstrom parameter.
      !
      lfactor = EXP(-angstrom(iplume) * LOG(lambda/550.0_wp))
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) ASYNC(1) IF(use_acc)
      DO k=1,nlevels
        DO icol = jcs,ncol
          aod_550          = prof(icol,k)     * cw_an(icol)
          aod_lmd          = aod_550          * lfactor
          caod_sp(icol)    = caod_sp(icol)    + aod_550
          caod_bg(icol)    = caod_bg(icol)    + prof(icol,k) * cw_bg(icol)
          asy_prof(icol,k) = asy_prof(icol,k) + aod_lmd * ssa(icol) * asy(icol)
          ssa_prof(icol,k) = ssa_prof(icol,k) + aod_lmd * ssa(icol)
          aod_prof(icol,k) = aod_prof(icol,k) + aod_lmd
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !
    ! complete optical depth weighting
    !
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
    DO k=1,nlevels
      DO icol = jcs,ncol
        ! VM: deleted MERGE which causes problems on cray if used in this manner
        !asy_prof(icol,k) = MERGE(asy_prof(icol,k)/ssa_prof(icol,k), 0.0_wp, ssa_prof(icol,k) > TINY(1._wp))
        !ssa_prof(icol,k) = MERGE(ssa_prof(icol,k)/aod_prof(icol,k), 1.0_wp, aod_prof(icol,k) > TINY(1._wp))
        IF (ssa_prof(icol,k) > TINY(1._wp)) THEN
          asy_prof(icol,k) = asy_prof(icol,k)/ssa_prof(icol,k)
        ELSE
          asy_prof(icol,k) = 0.0_wp
        END IF

        !vm test:
        IF (aod_prof(icol,k) > TINY(1._wp)) THEN
          ssa_prof(icol,k) = ssa_prof(icol,k)/aod_prof(icol,k)
        ELSE
          ssa_prof(icol,k) = 1.0_wp
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !
    ! calcuate effective radius normalization (divisor) factor
    !
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
    DO icol=jcs,ncol
      dNovrN(icol) = LOG((1000.0_wp * (caod_sp(icol) + caod_bg(icol))) + 1.0_wp)/LOG((1000.0_wp * caod_bg(icol)) + 1.0_wp)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
    !$ACC END DATA
    RETURN
  END SUBROUTINE sp_aop_profile
  !
  !PRECALC_SP_AOP_PROFILE: introfuced to avoid repetitive heavy computions in sp_aop_profile.
  !
  SUBROUTINE precalc_sp_aop_profile(  &
       jcs            ,ncol           ,ncol_max       ,lon            , &
     & lat            ,f1             ,f2             ,f3             , &
     & f4             ,opt_use_acc                                                  )
        INTEGER, INTENT(IN)        :: &
       & jcs,                     & !< start index in block
       & ncol,                    & !< number of columns (end index)
       & ncol_max                   !< first dimension of 2d-vars as declared in calling (sub)program [nproma]

    REAL(wp), INTENT(IN)       :: &
       & lon(ncol),               & !< longitude in degrees E
       & lat(ncol)                  !< latitude in degrees N

    REAL(wp), INTENT(INOUT)    :: &
       & f1(ncol_max, nplumes), f2(ncol_max, nplumes), f3(ncol_max, nplumes), f4(ncol_max, nplumes)

    INTEGER                     :: iplume, icol
    REAL(wp)                    :: &
       & delta_lat,                & !< latitude offset
       & delta_lon,                & !< longitude offset
       & delta_lon_t,              & !< threshold for maximum longitudinal plume extent used in transition from 360 to 0 degrees
       & a_plume1,                 & !< gaussian longitude factor for feature 1
       & a_plume2,                 & !< gaussian longitude factor for feature 2
       & b_plume1,                 & !< gaussian latitude factor for feature 1
       & b_plume2,                 & !< gaussian latitude factor for feature 2
       & lon1,                     & !< rotated longitude for feature 1
       & lat1,                     & !< rotated latitude for feature 2
       & lon2,                     & !< rotated longitude for feature 1
       & lat2                        !< rotated latitude for feature 2

    LOGICAL, INTENT(IN), OPTIONAL               :: opt_use_acc
    LOGICAL                                     :: use_acc
    IF (PRESENT(opt_use_acc)) use_acc = opt_use_acc

    IF (f1(jcs,1) < 0._wp) THEN
    DO iplume=1,nplumes
      !
      ! calculate plume weights
      !
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
      DO icol=jcs,ncol
        !
        ! get plume-center relative spatial parameters for specifying amplitude of plume at given lat and lon
        !
        delta_lat = lat(icol) - plume_lat(iplume)
        delta_lon = lon(icol) - plume_lon(iplume)
        delta_lon_t = MERGE (260._wp, 180._wp, iplume == 1)
        delta_lon = MERGE ( delta_lon-SIGN(360._wp,delta_lon) , delta_lon , ABS(delta_lon) > delta_lon_t)

        a_plume1  = 0.5_wp / (MERGE(sig_lon_E(1,iplume), sig_lon_W(1,iplume), delta_lon > 0.0_wp)**2)
        b_plume1  = 0.5_wp / (MERGE(sig_lat_E(1,iplume), sig_lat_W(1,iplume), delta_lon > 0.0_wp)**2)
        a_plume2  = 0.5_wp / (MERGE(sig_lon_E(2,iplume), sig_lon_W(2,iplume), delta_lon > 0.0_wp)**2)
        b_plume2  = 0.5_wp / (MERGE(sig_lat_E(2,iplume), sig_lat_W(2,iplume), delta_lon > 0.0_wp)**2)
        !
        ! adjust for a plume specific rotation which helps match plume state to climatology.
        !
        lon1 =   COS(theta(1,iplume))*(delta_lon) + SIN(theta(1,iplume))*(delta_lat)
        lat1 = - SIN(theta(1,iplume))*(delta_lon) + COS(theta(1,iplume))*(delta_lat)
        lon2 =   COS(theta(2,iplume))*(delta_lon) + SIN(theta(2,iplume))*(delta_lat)
        lat2 = - SIN(theta(2,iplume))*(delta_lon) + COS(theta(2,iplume))*(delta_lat)
        !
        ! calculate contribution to plume from its different features, to get a column weight for the anthropogenic
        ! (cw_an) and the fine-mode background aerosol (cw_bg)
        !
        f1(icol,iplume) = ftr_weight(1,iplume) * EXP(-1._wp* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2))))
        f2(icol,iplume) = ftr_weight(2,iplume) * EXP(-1._wp* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2))))
        f3(icol,iplume) = ftr_weight(1,iplume) * EXP(-1.* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2))))
        f4(icol,iplume) = ftr_weight(2,iplume) * EXP(-1.* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2))))
      END DO
   END DO
  END IF
  END SUBROUTINE precalc_sp_aop_profile

  ! -----------------------------------------------------------------------------------------------
  ! ADD_BC_AEROPT_SPLUMES:  This subroutine provides the interface to simple plume (sp) fit to the
  !                         MPI Aerosol Climatology (Version 2). It does so by collecting or
  !                         deriving spatio-temporal information and calling the simple plume
  !                         aerosol subroutine and incrementing the background aerosol properties
  !                         (and effective radius) with the anthropogenic plumes.
  !
  SUBROUTINE add_bc_aeropt_splumes_opt( jg                                             ,&
     & jcs            ,jce            ,nproma         ,klev           ,jb          ,&
     & nb_sw          ,this_datetime  ,zf             ,dz             ,z_sfc       ,&
     & sw_wv1         ,sw_wv2         ,aod_sw_vr      ,ssa_sw_vr      ,asy_sw_vr   ,&
     & x_cdnc         ,opt_use_acc                                                 )
    !
    ! --- 0.1 Variables passed through argument list
    INTEGER, INTENT(IN) ::            &
         jg                          ,& !< domain index
         jcs                         ,& !< start index in current block
         jce                         ,& !< end index in current block
         nproma                      ,& !< block dimension
         klev                        ,& !< number of full levels
         jb                          ,& !< index for current block
         nb_sw                          !< number of bands in short wave

    TYPE(datetime), INTENT(IN)   :: this_datetime

    REAL(wp), INTENT (IN)        :: &
         zf(nproma,klev),            & !< geometric height at full level [m]
         dz(nproma,klev),            & !< geometric height thickness     [m]
         z_sfc(nproma),              & !< geometric height of surface    [m]
         sw_wv1(nb_sw),              & !< smallest wave number in each of the sw bands
         sw_wv2(nb_sw)                !< largest  wave number in each of the sw bands

    REAL(wp), INTENT (INOUT) ::       &
         aod_sw_vr(nproma,klev,nb_sw) ,& !< Aerosol shortwave optical depth
         ssa_sw_vr(nproma,klev,nb_sw) ,& !< Aerosol single scattering albedo
         asy_sw_vr(nproma,klev,nb_sw)    !< Aerosol asymmetry parameter
    REAL(wp), INTENT(OUT), OPTIONAL:: &
         x_cdnc(nproma)                  !< Scale factor for Cloud Droplet Number Concentration

    !
    ! --- 0.2 Dummy variables
    !
    INTEGER ::                        &
         jk                          ,& !< index for looping over vertical dimension
         jki                         ,& !< index for looping over vertical dimension for reversing
         jl                          ,& !< index for looping over block
         jwl                            !< index for looping over wavelengths

    REAL(wp) ::                       &
         year_fr                     ,& !< time in year fraction (1989.0 is 0Z on Jan 1 1989)
         lambda                      ,& !< wavelength at central band wavenumber [nm]
         lon_sp(nproma)              ,& !< longitude passed to sp
         lat_sp(nproma)              ,& !< latitude passed to sp
         z_fl_vr(nproma,klev)        ,& !< level height [m], vertically reversed indexing (1=lowest level)
         dz_vr(nproma,klev)          ,& !< level thickness [m], vertically reversed
         sp_aod_vr(nproma,klev)      ,& !< simple plume aerosol optical depth, vertically reversed
         sp_ssa_vr(nproma,klev)      ,& !< simple plume single scattering albedo, vertically reversed
         sp_asy_vr(nproma,klev)      ,& !< simple plume asymmetry factor, vertically reversed indexing
         sp_xcdnc(nproma)               !< drop number scale factor
    REAL(wp), POINTER  ::             &
         aer_sp_f1(:,:),   &  !< geographical weighting factor
         aer_sp_f2(:,:),   &  !< geographical weighting factor
         aer_sp_f3(:,:),   &  !< geographical weighting factor
         aer_sp_f4(:,:)       !< geographical weighting factor

    LOGICAL, INTENT(IN), OPTIONAL           :: opt_use_acc
    LOGICAL                                 :: use_acc
    use_acc = .False.
    IF (PRESENT(opt_use_acc)) use_acc = opt_use_acc

    year_fr = REAL(this_datetime%date%year,wp) &
         +((REAL(getDayOfYearFromDateTime(this_datetime),wp) &
         +REAL(getNoOfSecondsElapsedInDayDateTime(this_datetime),wp)/86400.0_wp) &
         /REAL(getNoOfDaysInYearDateTime(this_datetime),wp))
    IF (this_datetime%date%year > 1850) THEN
    !$ACC DATA CREATE(lon_sp, lat_sp, z_fl_vr, dz_vr, sp_aod_vr, sp_ssa_vr, sp_asy_vr, sp_xcdnc) ASYNC(1) IF(use_acc)
      !
      ! --- 1.1 geographic information
      !

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
      DO jk=1,klev
        DO jl=jcs,jce
          dz_vr(jl,jk) = dz(jl,klev-jk+1)
          z_fl_vr(jl,jk) = zf(jl,klev-jk+1)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
      DO jl= jcs, jce
        lon_sp(jl) = p_patch(jg)%cells%center(jl, jb)%lon * rad2deg
        lat_sp(jl) = p_patch(jg)%cells%center(jl, jb)%lat * rad2deg
      END DO
      !$ACC END PARALLEL LOOP
      !
      ! --- 1.2 Aerosol Shortwave properties
      !
      ! get aerosol optical properties in each band, and adjust effective radius
      !
      aer_sp_f1=>bc_spl_field(jg)%aer_sp_f1(:,:,jb)
      aer_sp_f2=>bc_spl_field(jg)%aer_sp_f2(:,:,jb)
      aer_sp_f3=>bc_spl_field(jg)%aer_sp_f3(:,:,jb)
      aer_sp_f4=>bc_spl_field(jg)%aer_sp_f4(:,:,jb)

      CALL precalc_sp_aop_profile(  &
       jcs            ,jce            ,nproma         ,lon_sp         , &
     & lat_sp         ,aer_sp_f1      ,aer_sp_f2      ,aer_sp_f3      , &
     & aer_sp_f4      ,opt_use_acc=use_acc                              )


      DO jwl = 1,nb_sw
        lambda = 1.e7_wp/ (0.5_wp * (sw_wv1(jwl) + sw_wv2(jwl)))
        CALL sp_aop_profile(                                              klev                , &
           & jcs                ,jce                ,nproma              ,lambda              , &
           & z_sfc(:)           ,year_fr            ,aer_sp_f1           ,aer_sp_f2           , &
           & aer_sp_f3          ,aer_sp_f4                                                    , &
           & z_fl_vr(:,:)       ,dz_vr(:,:)         ,sp_xcdnc(:)         ,sp_aod_vr(:,:)      , &
           & sp_ssa_vr(:,:)     ,sp_asy_vr(:,:)     ,opt_use_acc=use_acc             )
        !$ACC PARALLEL LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
        DO jk=1,klev
          DO jl=jcs,jce
            asy_sw_vr(jl,jk,jwl) = asy_sw_vr(jl,jk,jwl) * ssa_sw_vr(jl,jk,jwl) * aod_sw_vr(jl,jk,jwl)    &
                 + sp_asy_vr(jl,jk)   * sp_ssa_vr(jl,jk)    * sp_aod_vr(jl,jk)
            ssa_sw_vr(jl,jk,jwl) = ssa_sw_vr(jl,jk,jwl) * aod_sw_vr(jl,jk,jwl)                           &
                 + sp_ssa_vr(jl,jk)   * sp_aod_vr(jl,jk)
            aod_sw_vr(jl,jk,jwl) = aod_sw_vr(jl,jk,jwl) + sp_aod_vr(jl,jk)
            IF (ssa_sw_vr(jl,jk,jwl) > TINY(1.0_wp)) THEN
              asy_sw_vr(jl,jk,jwl) = asy_sw_vr(jl,jk,jwl)/ssa_sw_vr(jl,jk,jwl)
            ELSE
              asy_sw_vr(jl,jk,jwl) = asy_sw_vr(jl,jk,jwl)
            END IF

            IF (aod_sw_vr(jl,jk,jwl) > TINY(1.0_wp)) THEN
              ssa_sw_vr(jl,jk,jwl) = ssa_sw_vr(jl,jk,jwl)/aod_sw_vr(jl,jk,jwl)
            ELSE
              ssa_sw_vr(jl,jk,jwl) = ssa_sw_vr(jl,jk,jwl)
            END IF
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      IF (PRESENT (x_cdnc)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
        DO jl=jcs,jce
          x_cdnc(jl) = sp_xcdnc(jl)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
      !$ACC WAIT(1)
      !$ACC END DATA
      RETURN
    END IF
  END SUBROUTINE add_bc_aeropt_splumes_opt

  SUBROUTINE read_1d_wrapper(ifile_id,                 variable_name,        &
                           & alloc_array,              file_name,            &
                           & variable_dimls,           module_name,          &
                           & sub_prog_name                                   )
    INTEGER, INTENT(in)            :: ifile_id      !< file id from which
                                                    !< variable is read
    CHARACTER(LEN=*),INTENT(in)    :: variable_name !< name of variable
                                                    !< to be read
    REAL(wp), ALLOCATABLE, INTENT(out), TARGET  :: alloc_array(:) !< values of variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: file_name     !< file name of file
                                                    !< contain respective var.
    INTEGER, INTENT(in),OPTIONAL   :: variable_dimls(1)!< dimension length of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: module_name!< name of module
                                           !< containing calling subprogr.
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: sub_prog_name!< name of calling
                                                        !< subprogram
    CHARACTER(LEN=32)                   :: ci_length, cj_length
    CHARACTER(LEN=1024)                 :: message1, message2

    CALL read_1D(file_id=ifile_id,         variable_name=variable_name,      &
                 alloc_array=alloc_array                                     )
    IF (PRESENT(variable_dimls)) THEN
       IF (SIZE(alloc_array,1)/=variable_dimls(1)) THEN
         WRITE(ci_length,*) SIZE(alloc_array,1)
         WRITE(cj_length,*) variable_dimls(1)
         IF (PRESENT(sub_prog_name)) THEN
           message1=TRIM(ADJUSTL(sub_prog_name))//' of'
         ELSE
           message1='Unknown subprogram of '
         END IF
         IF (PRESENT(module_name)) THEN
           message1=TRIM(ADJUSTL(message1))//' '//TRIM(ADJUSTL(module_name))
         ELSE
           message1=TRIM(ADJUSTL(message1))//' unknown module'
         END IF
         message2=TRIM(ADJUSTL(variable_name))//'('// &
                & TRIM(ADJUSTL(cj_length))//') has wrong dimension length '// &
                & TRIM(ADJUSTL(ci_length))//' in'
         IF (PRESENT(file_name)) THEN
           message2=TRIM(ADJUSTL(message2))//' '//TRIM(ADJUSTL(file_name))
         ELSE
           message2=TRIM(ADJUSTL(message2))//' unknown file'
         END IF
         WRITE(0,*) TRIM(ADJUSTL(message1))
         WRITE(0,*) TRIM(ADJUSTL(message2))
         CALL finish(TRIM(ADJUSTL(message1)),TRIM(ADJUSTL(message2)))
       END IF
    END IF
  END SUBROUTINE read_1d_wrapper
  SUBROUTINE read_2d_wrapper(ifile_id,                 variable_name,        &
                           & alloc_array,              file_name,            &
                           & variable_dimls,           module_name,          &
                           & sub_prog_name                                   )
    INTEGER, INTENT(in)            :: ifile_id      !< file id from which
                                                    !< variable is read
    CHARACTER(LEN=*),INTENT(in)    :: variable_name !< name of variable
                                                    !< to be read
    REAL(wp), ALLOCATABLE, INTENT(out), TARGET  :: alloc_array(:,:) !< values of variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: file_name     !< file name of file
                                                    !< contain respective var.
    INTEGER, INTENT(in),OPTIONAL   :: variable_dimls(2)!< dimension length of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: module_name!< name of module
                                           !< containing calling subprogr.
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: sub_prog_name!< name of calling
                                                        !< subprogram
    CHARACTER(LEN=32)                   :: ci_length(2), cj_length(2)
    CHARACTER(LEN=1024)                 :: message1, message2

    CALL read_bcast_REAL_2D(file_id=ifile_id,                                &
                         &  variable_name=variable_name,                     &
                         &  alloc_array=alloc_array                          )
    IF (PRESENT(variable_dimls)) THEN
       IF (SIZE(alloc_array,1)/=variable_dimls(1) .OR. &
         & SIZE(alloc_array,2)/=variable_dimls(2)) THEN
         WRITE(ci_length(1),*) SIZE(alloc_array,1)
         WRITE(cj_length(1),*) variable_dimls(1)
         WRITE(ci_length(2),*) SIZE(alloc_array,2)
         WRITE(cj_length(2),*) variable_dimls(2)
         IF (PRESENT(sub_prog_name)) THEN
           message1=TRIM(ADJUSTL(sub_prog_name))//' of'
         ELSE
           message1='Unknown subprogram of '
         END IF
         IF (PRESENT(module_name)) THEN
           message1=TRIM(ADJUSTL(message1))//' '//TRIM(ADJUSTL(module_name))
         ELSE
           message1=TRIM(ADJUSTL(message1))//' unknown module'
         END IF
         message2=TRIM(ADJUSTL(variable_name))//'('//&
                & TRIM(ADJUSTL(cj_length(1)))//','//&
                & TRIM(ADJUSTL(cj_length(2)))//&
                & ') has wrong dimension length ('//&
                & TRIM(ADJUSTL(ci_length(1)))//','//&
                & TRIM(ADJUSTL(ci_length(2)))//&
                & ') in'
         IF (PRESENT(file_name)) THEN
           message2=TRIM(ADJUSTL(message2))//' '//TRIM(ADJUSTL(file_name))
         ELSE
           message2=TRIM(ADJUSTL(message2))//' unknown file'
         END IF
         WRITE(0,*) TRIM(ADJUSTL(message1))
         WRITE(0,*) TRIM(ADJUSTL(message2))
         CALL finish(TRIM(ADJUSTL(message1)),TRIM(ADJUSTL(message2)))
       END IF
    END IF
  END SUBROUTINE read_2d_wrapper
  SUBROUTINE read_3d_wrapper(ifile_id,                 variable_name,        &
                           & alloc_array,              file_name,            &
                           & variable_dimls,           module_name,          &
                           & sub_prog_name                                   )
    INTEGER, INTENT(in)            :: ifile_id      !< file id from which
                                                    !< variable is read
    CHARACTER(LEN=*),INTENT(in)    :: variable_name !< name of variable
                                                    !< to be read
    REAL(wp), ALLOCATABLE, INTENT(out), TARGET  :: alloc_array(:,:,:) !< values of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: file_name     !< file name of file
                                                    !< contain respective var.
    INTEGER, INTENT(in),OPTIONAL   :: variable_dimls(3)!< dimension length of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: module_name!< name of module
                                           !< containing calling subprogr.
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: sub_prog_name!< name of calling
                                                        !< subprogram
    CHARACTER(LEN=32)                   :: ci_length(3), cj_length(3)
    CHARACTER(LEN=1024)                 :: message1, message2

    CALL read_bcast_REAL_3D(file_id=ifile_id,                                &
                         &  variable_name=variable_name,                     &
                         &  alloc_array=alloc_array                          )
    IF (PRESENT(variable_dimls)) THEN
       IF (SIZE(alloc_array,1)/=variable_dimls(1) .OR. &
         & SIZE(alloc_array,2)/=variable_dimls(2) .OR. &
         & SIZE(alloc_array,3)/=variable_dimls(3)) THEN
         WRITE(ci_length(1),*) SIZE(alloc_array,1)
         WRITE(cj_length(1),*) variable_dimls(1)
         WRITE(ci_length(2),*) SIZE(alloc_array,2)
         WRITE(cj_length(2),*) variable_dimls(2)
         WRITE(ci_length(3),*) SIZE(alloc_array,3)
         WRITE(cj_length(3),*) variable_dimls(3)
         IF (PRESENT(sub_prog_name)) THEN
           message1=TRIM(ADJUSTL(sub_prog_name))//' of'
         ELSE
           message1='Unknown subprogram of '
         END IF
         IF (PRESENT(module_name)) THEN
           message1=TRIM(ADJUSTL(message1))//' '//TRIM(ADJUSTL(module_name))
         ELSE
           message1=TRIM(ADJUSTL(message1))//' unknown module'
         END IF
         message2=TRIM(ADJUSTL(variable_name))//'('//&
                & TRIM(ADJUSTL(cj_length(1)))//','//&
                & TRIM(ADJUSTL(cj_length(2)))//','//&
                & TRIM(ADJUSTL(cj_length(3)))//&
                & ') has wrong dimension length ('//&
                & TRIM(ADJUSTL(ci_length(1)))//','//&
                & TRIM(ADJUSTL(ci_length(2)))//','//&
                & TRIM(ADJUSTL(ci_length(3)))//&
                & ') in'
         IF (PRESENT(file_name)) THEN
           message2=TRIM(ADJUSTL(message2))//' '//TRIM(ADJUSTL(file_name))
         ELSE
           message2=TRIM(ADJUSTL(message2))//' unknown file'
         END IF
         CALL finish(TRIM(ADJUSTL(message1)),TRIM(ADJUSTL(message2)))
       END IF
    END IF
  END SUBROUTINE read_3d_wrapper

END MODULE mo_bc_aeropt_splumes_opt
