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

! Routines for the initialization of the wave const state.
! The const state includes time-constant index arrays and weights
! for source term evaluation and integration.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_phy_init

  USE mo_kind,                 ONLY: wp
  USE mo_mpi,                  ONLY: my_process_is_stdio
  USE mo_exception,            ONLY: message, finish
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_math_constants,       ONLY: rad2deg
  USE mo_wave_types,           ONLY: t_wave_const
  USE mo_wave_config,          ONLY: t_wave_config
  USE mo_wave_phy_util,        ONLY: jonswap_nonblk, fetch_law_nonblk

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_wave_phy

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_phy_init'

CONTAINS

  !>
  !! Initialize time-constant index arrays and weights
  !! for the computation of wave energy source terms and
  !! wave source time integration.
  !!
  SUBROUTINE init_wave_phy(wave_config, const)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':init_wave_phy'

    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wave_const),          INTENT(INOUT) :: const

    TYPE(t_wave_config), POINTER :: wc => NULL()

    ! Calculate the permissible minimum values of energy for a
    ! a given range of wind speed bins for each frequency
    ! from 1 to wave_config%jmax.
    !
    CALL min_energy(wave_config = wave_config,       & !IN
      &             flminfr_tab = const%flminfr_tab)   !INOUT

    ! initialisation of the nonlinear transfer computations
    ! computes time-constant index arrays and weights
    !
    CALL init_wave_nonlinear(wave_config = wave_config,  & !IN
      &                      const       = const)          !INOUT

  END SUBROUTINE init_wave_phy


  !>
  !! Calculation of index arrays and weights for the computation of
  !! the nonlinear transfer rate for shallow water.
  !!
  !! Computation of parameters used in discrete interaction
  !! parameterization of nonlinear transfer.
  !!
  !! Adoptation of NLWEIGT from WAM 4.5.
  !!
  !! SUSANNE HASSELMANN JUNE 86.
  !! H. GUNTHER   ECMWF/GKSS  DECEMBER 90 - CYCLE_4 MODIFICATIONS.
  !! P. Janssen   ECMWF June 2005                                         !
  !! H. Gunther   HZG   January 2015  cycle_4.5.4
  !!
  !! Reference
  !! S. Hasselmann and K. Hasselmann, JPO, 1985
  !!
  SUBROUTINE init_wave_nonlinear(wave_config, const)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':init_wave_nonlinear'

    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wave_const),          INTENT(INOUT) :: const

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: jf,jd
    INTEGER :: error

    INTEGER :: nfreqs, ndirs
    INTEGER :: klp1, ic, kh, klh, k, ks, icl1, icl2, isg, k1, k11, k2, k21
    INTEGER :: m, ikn, i, ie

    REAL(wp) :: deltha, cl1, cl2, al11, al12, ch, cl1h, cl2h
    REAL(wp) :: f1p1, frg, flp, flm, fkp, fkm

    REAL(wp), ALLOCATABLE, DIMENSION(:) :: frlon

    ! Parameters for discrete approximation of nonlinear transfer
    REAL(wp), PARAMETER :: alamd   = 0.25_wp   ! lambda
    REAL(wp), PARAMETER :: con     = 3000.0_wp ! weight for discrete approximation of nonlinear transfer
    REAL(wp), PARAMETER :: delphi1 = -11.48_wp
    REAL(wp), PARAMETER :: delphi2 = 33.56_wp


    wc => wave_config
    nfreqs = wc%nfreqs
    ndirs = wc%ndirs

    ALLOCATE(frlon(2*nfreqs+2), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    ! 1. Computation for angular grid
    deltha = wc%delth * rad2deg

    cl1 = delphi1/deltha
    cl2 = delphi2/deltha

    ! 1.1 computation of indices of angular cell.
    klp1 = ndirs+1
    ic = 1

    DO kh = 1,2
       klh = ndirs
       IF (kh.eq.2) klh=klp1
       DO k = 1,klh
          ks = k
          IF (kh.gt.1) ks=klp1-k+1
          IF (ks.gt.ndirs) CYCLE
          ch = ic*cl1
          const%ja1(ks,kh) = jafu(ch,k,ndirs)
          ch = ic*cl2
          const%ja2(ks,kh) = jafu(ch,k,ndirs)
       END DO
       ic = -1
    END DO

    ! 1.2 computation of angular weights
    icl1 = cl1
    cl1  = cl1 - icl1
    icl2 = cl2
    cl2  = cl2 - icl2
    wc%acl1 = ABS(cl1)
    wc%acl2 = ABS(cl2)
    wc%cl11 = 1._wp - wc%acl1
    wc%cl21 = 1._wp - wc%acl2
    al11 = (1._wp + alamd)**4
    al12 = (1._wp - alamd)**4
    wc%dal1 = 1._wp / al11
    wc%dal2 = 1._wp / al12

    ! 1.3 computation of angular indices
    isg = 1
    DO kh = 1,2
       cl1h = isg*cl1
       cl2h = isg*cl2
       DO k = 1,ndirs
          ks = k
          IF (kh.eq.2) ks = ndirs-k+2
          IF(k.eq.1) ks = 1
          k1 = const%ja1(k,kh)
          const%k1w(ks,kh) = k1
          IF (cl1h.lt.0.) THEN
             k11 = k1-1
             IF (k11.lt.1) k11 = ndirs
          ELSE
             k11 = k1+1
             IF (k11.gt.ndirs) k11 = 1
          END IF
          const%k11w(ks,kh) = k11
          k2 = const%ja2(k,kh)
          const%k2w(ks,kh) = k2
          IF (cl2h.lt.0) THEN
             k21 = k2-1
             IF(k21.lt.1) k21 = ndirs
          ELSE
             k21 = k2+1
             IF (k21.gt.ndirs) k21 = 1
          END IF
          const%k21w(ks,kh) = k21
       END DO
       isg = -1
    END DO

    ! 2. computation for frequency grid
    frlon(1:nfreqs) = wc%freqs(1:nfreqs)

    DO m = nfreqs+1,2*nfreqs+2
       frlon(m) = wc%co*frlon(m-1)
    END DO

    f1p1 = LOG10(wc%co)
    DO m = 1,nfreqs+4
       frg = frlon(m)
       const%af11(m) = con * frg**11
       flp = frg*(1.+alamd)
       flm = frg*(1.-alamd)
       ikn = INT(LOG10(1._wp+alamd)/f1p1+.000001_wp)
       ikn = m+ikn
       const%ikp(m) = ikn
       fkp = frlon(const%ikp(m))
       const%ikp1(m) = const%ikp(m)+1
       const%fklap(m) = (flp-fkp)/(frlon(const%ikp1(m))-fkp)

       const%fklap1(m) = 1._wp-const%fklap(m)
       IF (frlon(1).ge.flm) THEN
          const%ikm(m) = 1
          const%ikm1(m) = 1
          const%fklam(m) = 0._wp
          const%fklam1(m) = 0._wp
       ELSE
          ikn = INT(LOG10(1._wp-alamd)/f1p1+.0000001_wp)
          ikn = m+ikn-1
          IF (ikn.lt.1) ikn = 1
          const%ikm(m) = ikn
          fkm = frlon(const%ikm(m))
          const%ikm1(m) = const%ikm(m)+1
          const%fklam(m) = (flm-fkm)/(frlon(const%ikm1(m))-fkm)

          const%fklam1(m) = 1._wp-const%fklam(m)
       END IF
    END DO


    ! 3. compute tail frequency ratios
    wc%frh(:) = 0._wp ! initialisation

    ie = MIN(30,nfreqs+3)
    DO i = 1,ie
       m = nfreqs+i-1
       wc%frh(i) = (frlon(nfreqs)/frlon(m))**5
    END DO

    IF (ALLOCATED(frlon))          DEALLOCATE(frlon)

    !print nonlinear status
    IF (my_process_is_stdio()) THEN
       WRITE(0,'(/,'' -------------------------------------------------'')')
       WRITE(0,*)'        non linear interaction parameters'
       WRITE(0,'(  '' -------------------------------------------------'')')
       WRITE(0,'(/,''  frequency arrays'')')
       WRITE(0,'(''     acl1       acl2       cl11       cl21   '',                &
            &            ''    dal1       dal2'')')
       WRITE(0,'(1x,6f11.8)') wc%acl1, wc%acl2, wc%cl11, wc%cl21, wc%dal1, wc%dal2
       WRITE(0,*) ' '
       WRITE(0,'(''  m   ikp ikp1  ikm ikm1   fklap       fklap1 '',               &
            &            ''   fklam       fklam1     af11'')')

       DO jf = 1,size(const%ikp)
          WRITE(0,'(1x,i2,4i5,4f11.8,e11.3)') jf, const%ikp(jf), const%ikp1(jf), const%ikm(jf), const%ikm1(jf), &
               &            const%fklap(jf), const%fklap1(jf), const%fklam(jf), const%fklam1(jf), const%af11(jf)
       END DO

       WRITE(0,'(/,''  angular arrays'')')
       WRITE(0,'(''   |--------kh = 1----------||--------kh = 2----------|'')')
       WRITE(0,'(''  k   k1w   k2w  k11w  k21w   k1w   k2w  k11w  k21w'')')
       DO jd = 1,size(const%k1w,1)
          WRITE(0,'(1x,i2,8i6)') jd,(const%k1w(jd,kh), const%k2w(jd,kh), const%k11w(jd,kh),              &
               &                            const%k21w(jd,kh),kh=1,2)
       END DO

       WRITE(0,'(/,''  tail array frh'')')
       WRITE(0,'(1x,8f10.7)') wc%frh(1:30)
    END IF


    CALL message(routine,'finished')

  END SUBROUTINE init_wave_nonlinear


  !>
  !! Function to compute the index array for the angles of the
  !! interacting wavenumbers.
  !!
  !! Adopted from WAM 4.5 JAFU
  !!
  !!  S. Hasselmann        MPIFM        01/12/1985
  !!
  !! Indices defining bins in frequency and direction plane into
  !! which nonlinear energy transfer increments are stored. Needed
  !! for computation of the nonlinear energy transfer.
  !!
  !! Reference
  !! S. Hasselmann and K. Hasselmann, JPO, 1985 B
  INTEGER FUNCTION JAFU (CL, J, IAN)

    REAL(wp),    INTENT(IN) :: CL !! weights.
    INTEGER, INTENT(IN) :: J      !! index in angular array.
    INTEGER, INTENT(IN) :: IAN    !! number of angles in array.

    JAFU = J + INT(CL)
    IF (JAFU.LE.0)   JAFU = JAFU+IAN
    IF (JAFU.GT.IAN) JAFU = JAFU-IAN

  END FUNCTION JAFU


  !>
  !! Computation of minimum energy in spectral bins.
  !!
  !! Adopted from WAM 4.5 MIN_ENERGY
  !!
  !! J. Bidlot    ECMWF
  !!
  !! Computate of a table for minimum energy in spectral bins.
  !!
  !! Method:
  !!
  !! For each windspeed u10 as listed in the stress table the JONSWAP
  !! spectrum is computed from fetch laws and 1% of its value is stored for
  !! each frequency bin. Too small values are replaced by a given minimum.
  !!
  SUBROUTINE min_energy(wave_config, flminfr_tab)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':min_energy'

    TYPE(t_wave_config), INTENT(IN)    :: wave_config
    REAL(wp),            INTENT(INOUT) :: flminfr_tab(:,:) ! minimum energy for a given frequency and wind spped bin (jmax,nfreqs)

    INTEGER  :: j, nfreqs, jmax

    REAL(wp) :: flmin       !! absolute minimum energy in spectral bins
    REAL(wp) :: delu        !! wind increment

    REAL(wp) :: u10(SIZE(flminfr_tab,1))
    REAL(wp) :: fpk(SIZE(flminfr_tab,1))      ! JONSWAP peak frequency
    REAL(wp) :: alphaj0(SIZE(flminfr_tab,1))  ! JONSWAP alpha

    nfreqs = wave_config%nfreqs
    jmax = wave_config%jmax
    flmin = wave_config%flmin
    delu = wave_config%delu

    ! sanity check
    IF (SIZE(flminfr_tab,1) /= jmax) THEN
      CALL finish(routine, "Size of field flminfr_tab must be equal to jmax.")
    ENDIF

    ! compute windspeeds from table parameters
    DO j = 1,jmax
      u10(j) = REAL(j,wp)*delu
    END DO

    ! peak frequencies and alpha parameter from fetch law
    CALL fetch_law_nonblk(                              &
          & i_startidx  = 1,                            & !in
          & i_endidx    = jmax,                         & !in
          & fetch       = wave_config%fetch_min_energy, & !in
          & fpmax       = wave_config%freqs(nfreqs),    & !in
          & sp10m       = u10(:),                       & !in
          & fp          = fpk(:),                       & !out
          & alphaj      = alphaj0(:))                     !out

    ! JONSWAP spectra
    CALL jonswap_nonblk(                         &
          & i_startidx = 1,                      & !in
          & i_endidx   = jmax,                   & !in
          & freqs      = wave_config%freqs,      & !in
          & gamma      = wave_config%GAMMA_wave, & !in
          & sa         = wave_config%SIGMA_A,    & !in
          & sb         = wave_config%SIGMA_B,    & !in
          & flmin      = flmin,                  & !in
          & alphaj     = alphaj0(:)*0.01_wp,     & !in
          & FP         = fpk(:),                 & !in
          & ET         = flminfr_tab(:,:))         !out
  END SUBROUTINE min_energy

END MODULE mo_wave_phy_init
