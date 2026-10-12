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

! -------------------------------------------------------------------------------
! on DWD's NEC it is recommended to use the Advanced Scientific Library which
! provides vectorized spherical harmonics and a random number generator.
! For this -lasl_sequential has to be added to the VE configure script (adding
! it to LAPACK_LIBS works fine) and the pre-processor flag __ASL__
! can be defined here or in the VE configure script using -D__ASL__.
! -------------------------------------------------------------------------------
! #ifdef __NEC__
! #define __ASL__
! #endif

MODULE mo_stoch_pattern_generator

  !------------------------------------------------------------------------------
  !
  ! Description:
  !  This module contains a stochastic pattern generator using spherical harmonics
  !  based on ECMWF's SPPT scheme as described in Palmer et al. (2009)
  !
  !==============================================================================

  USE mo_kind,               ONLY: dp, sp, wp, i4, i8
  USE mo_math_constants,     ONLY: pi, rpi_2, pi_2, rad2deg
  USE mo_math_legendre,      ONLY: PlmON, PlmIndex
  USE mo_physical_constants, ONLY: earth_radius
  USE mtime,                 ONLY: datetime
  USE mo_exception,          ONLY: finish, message, txt=>message_text
  USE mo_mpi,                ONLY: get_my_global_mpi_id
  USE mo_gribout_config,     ONLY: gribout_config
  USE mo_run_config,         ONLY: msg_level
  USE mo_grid_config,        ONLY: l_limited_area

  ! all those are only needed for the stochastic_pattern_boundaries subroutine
  USE mo_model_domain,       ONLY: t_patch
  USE mo_sync,               ONLY: global_max, global_min
  USE mo_impl_constants,     ONLY: min_rlcell_int, SUCCESS
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_loopindices,        ONLY: get_indices_c

  ! for the NEC Auroa VEs is it convenient to use the vectorized Advanced Scientific Library (ASL)
#ifdef __ASL__
  USE asl_unified, ONLY : asl_library_initialize,        &
                          asl_library_is_initialized,    &
                          asl_random_create,             &
                          asl_random_initialize,         &
                          asl_random_generate_s,         &
                          asl_random_distribute_normal,  &
                          asl_random_distribute_normal_box_muller
#endif

  !==============================================================================

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: stochastic_pattern_init, stochastic_pattern_step, stochastic_pattern_generator, &
       &    stochastic_pattern_boundaries, stochastic_pattern_destruct

  !==============================================================================

  INTEGER, PARAMETER :: nmaxwave = 250      ! largest zonal wave number
  INTEGER, PARAMETER :: smax = 5            ! maximum number of scales
                                            ! note: smax=5 is hardcoded in the ICON namelist modules
  INTEGER :: nmax, ntime

  COMPLEX(KIND=sp), ALLOCATABLE, DIMENSION(:,:,:) ::  rcoeff   ! spectral coefficient in single precision
  REAL(KIND=sp),    ALLOCATABLE, DIMENSION(:,:)   ::  sigma    ! spectral variance for global (spherical harmonics)
  REAL(KIND=sp),    ALLOCATABLE, DIMENSION(:,:,:) ::  gcoeff   ! spectral variance for limited area (Fourier modes)

  REAL(KIND=sp), DIMENSION(smax) :: tau_spg
  INTEGER,       DIMENSION(smax) :: nstart, nend, nlen

  INTEGER :: itype_random_normals = 0            ! 0: built-in, 1: NEC's ASL
  INTEGER :: itype_legendre_polys = 0            ! 0: SHTOOLS,  1: NEC's ASL

  INTEGER :: randomhandle       ! unique handle for ASL random number generator
  INTEGER :: nstep              ! multiple of model timestep
  INTEGER :: nspg               ! number of spatial patterns
  INTEGER :: nscales            ! number of spatio-temporal scales
  LOGICAL :: lfourier

  REAL(KIND=wp) :: lat_min, lat_max, lon_min, lon_max, lat_ctr, lon_ctr, lat_len, lon_len

  REAL (sp), PARAMETER ::  pi1 = REAL(pi,kind=sp)
  REAL (sp), PARAMETER ::  pi2 = REAL(pi,kind=sp)*2.0_sp
  REAL (sp), PARAMETER ::  pi4 = REAL(pi,kind=sp)*4.0_sp
  REAL (sp), PARAMETER ::  pi8 = REAL(pi,kind=sp)*8.0_sp

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_stoch_pattern_generator'

CONTAINS

  !==============================================================================
  !  Module procedure to allocate and initialize some coefficients.
  !  This is called from mo_nwp_phy_init.
  !------------------------------------------------------------------------------

  SUBROUTINE stochastic_pattern_init(dtime,mtime_current,pspg,plam,plength,ptime,pvar,pasl,pmodes)
    REAL(KIND=wp),  INTENT(IN)               :: dtime          ! time step
    INTEGER,        INTENT(IN)               :: pspg           ! number of patterns in grid point space
    REAL(KIND=wp),  INTENT(IN), DIMENSION(:) :: plength        ! pattern length scale
    REAL(KIND=wp),  INTENT(IN), DIMENSION(:) :: ptime          ! pattern time scale
    REAL(KIND=wp),  INTENT(IN), DIMENSION(:) :: pvar           ! pattern variance
    INTEGER,        INTENT(IN), DIMENSION(:) :: pmodes         ! number of wave modes
    LOGICAL,        INTENT(IN)               :: plam           ! Fourier mode
    LOGICAL,        INTENT(IN)               :: pasl           ! ASL library on NEC
    TYPE(datetime), INTENT(IN)               :: mtime_current  ! current datetime

    ! local variables
    REAL(KIND=wp) :: dlength
    REAL(KIND=sp) :: fsum, sigm
    INTEGER :: iseed, nstart0, nend0, ndim
    INTEGER :: n,m,i,j,jj,k,nloc,kmax,istat
    COMPLEX(KIND=sp), DIMENSION(nmaxwave*nmaxwave/2+nmaxwave+1) :: random_normals
    REAL(KIND=sp),    DIMENSION(nmaxwave,smax)                  :: sigmacoeff

    REAL(KIND=sp), DIMENSION(smax) ::    &
       phi,         & ! one-step correlation of AR1 process
       fzero,       & ! normalization factor of spectral coefficients
       var_r,       & ! variance in grid point space
       kappaT         ! determines the horizontal length scale (but see pattern length below)

    nspg = pspg
    nmax = nmaxwave

    ! set number of scales based on length scale entries and check for consistency
    nscales  = COUNT(plength(1:nspg) > 0.0_wp)
    IF (nscales /= COUNT(ptime > 0.0_wp)) THEN
      IF(istat /= SUCCESS) CALL finish(modname, 'number of non-zero length and time scales has to be the same')
    END IF
    IF (nscales /= COUNT(pvar > 0.0_wp)) THEN
      IF(istat /= SUCCESS) CALL finish(modname, 'number of non-zero variances has to match length scales')
    END IF
    IF (nscales /= COUNT(pmodes > 0.0_wp)) THEN
      IF(istat /= SUCCESS) CALL finish(modname, 'number of non-zero modes has to match length scales')
    END IF
    IF (nscales > smax) THEN
      IF(istat /= SUCCESS) CALL finish(modname, 'number of length scales is larger than smax=6')
    END IF

    IF (pasl) THEN
      itype_random_normals = 1     ! use NEC's ASL
      itype_legendre_polys = 1     ! use NEC's ASL
    END If

    IF (.NOT.l_limited_area) THEN
      ! For global we set lfourier always to false and ignore actual namelist
      lfourier = .FALSE.
    ELSE
      ! Pattern generator using spherical harmonics or limited-area with Fourier modes.
      ! Note that spherical harmonics can be used for LAM as well because the
      ! spectral pattern generator is completely independent from the ICON grid
      ! and the transformation back into spatial space in only done for the
      ! actual grid points. The choice is mostly a question of computational efficiency.
      lfourier = plam
    END IF

    IF (lfourier) THEN
      lat_ctr = lat_min                ! set center to lower left corner
      lon_ctr = lon_min                ! and make spectral domain big enough
      lat_len = 2.0*(lat_max-lat_min)  ! to avoid any symmetry and periodicity
      lon_len = 2.0*(lon_max-lon_min)  ! in the limited-area domain
      dlength = earth_radius * MIN(lat_len,lon_len)
    ELSE
      dlength = earth_radius
    END IF

    ! time step of pattern generator has to be smaller than tau_spg but can be
    ! larger than dtime. Here we set nstep to run the pattern generator at an
    ! integer multiple of dtime with an upper limit of 10
    nstep = MAX(MIN(INT(MINVAL(ptime)/(4.0_wp*dtime)),10),1)

    ! initialize random number generator
    iseed = mtime_current%date%year  + mtime_current%time%minute * 13  &
          + mtime_current%date%month * mtime_current%date%day          &
          + mtime_current%time%hour  * 42                              &
          + gribout_config(1)%perturbationNumber * 3 + 1

    ! get size of random seed
    CALL random_seed(size=k)

    ! initialize random number generator
    CALL random_initialize(k, iseed)

    DO j=1,nscales   ! loop over scales

      ! core parameters of stochastic pattern generator set from namelist
      kappaT(j)  = REAL( (plength(j)/dlength)**2, kind=sp) * 0.5_sp
      tau_spg(j) = REAL( ptime(j), kind=sp)
      var_r(j)   = REAL( pvar(j),  kind=sp)

      ! one-step correlation of AR1 process
      phi(j) = EXP( -REAL(dtime,kind=sp)/tau_spg(j)*nstep)

      IF (lfourier) THEN

        ! use full wave number spectrum up to pmodes in limited area
        nstart(j) = 1
        nend(j)   = pmodes(j)
        nlen(j)   = (nend(j)-nstart(j)+1)*(nend(j)-nstart(j)+1)
        nloc      = 0

        ! Equation (A4) of appendix of Berner et al. (2015, MWR)
        ! (see also Thompson et al. 2021, MWR)
        fsum = 0.0_sp
        DO n=nstart(j),nend(j)
          DO m=nstart(j),nend(j)
            fsum = fsum + EXP(-pi8*kappaT(j)*(n*n+m*m))
          END DO
        END DO
        ! additional factor 4 in denominator because we use only 1/4 of the spectral domain
        ! and factor (2*pi)**2 as normalization of spectral (normalized) domain size
        fzero(j) = SQRT( var_r(j) * (1.0_sp - phi(j)**2) / (8.0_sp*fsum*pi2*pi2) )

      ELSE

        ! calculate sigma_n using equations (17) and (18) of Palmer et al.
        fsum = 0.0_sp
        DO n=1,nmax
          fsum = fsum + (2*n+1) * EXP(-kappaT(j)*n*(n+1))
        END DO
        fzero(j) = SQRT( var_r(j) * (1.0_sp - phi(j)**2) / (2.0_sp*fsum) )
        DO n=1,nmax
          sigmacoeff(n,j) = fzero(j) * EXP( -kappaT(j) * n*(n+1)/2.0_sp )
        END DO

        ! determine dominant wave number
        nloc = 0
        sigm = 0.0_sp
        DO n=1,nmax
          IF ( n*sigmacoeff(n,j).GT.sigm ) THEN
            sigm = n*sigmacoeff(n,j)
            nloc = n
          END IF
        END DO

        ! limit to most relevant wave modes
        IF (pmodes(j).GT.0) THEN
          nend(j)   = nloc+pmodes(j)/2
          nstart(j) = nloc-pmodes(j)/2
        ELSE
          nend(j)   = nmax
          nstart(j) = nloc/3
        ENDIF

        IF (msg_level > 0 .AND. nend(j) > nmax) &
             CALL message(modname,'WARNING: nend(j) > nmax, increasing nmaxwave might be necessary')

        nend(j)   = MIN(nend(j),nmax)
        nstart(j) = MAX(nstart(j),1)
        nlen(j)   = nend(j)*(nend(j)+3)/2 - (nstart(j)-3)*nstart(j)/2 + 1

      END IF
      !
    END DO

    ! allocate array for spherical harmonics coefficients for random pattern generator
    ! The if(allocated) statement is needed because the IAU calls phy_init twice
    IF (.NOT.ALLOCATED(rcoeff)) THEN

      nstart0 = MINVAL(nstart(1:nspg))
      nend0   = MAXVAL(nend(1:nspg))
      ndim    = nspg*nscales

      IF (lfourier) THEN
        ALLOCATE(rcoeff(0:nend0,nstart0:nend0,ndim), STAT=istat)
        IF(istat /= SUCCESS) CALL finish(modname, 'Allocation of rcoeff failed')
        ALLOCATE(gcoeff(nstart0:nend0,nstart0:nend0,ndim), STAT=istat)
        IF(istat /= SUCCESS) CALL finish(modname, 'Allocation of gcoeff failed')
        rcoeff(:,:,:) = 0.0_wp
        gcoeff(:,:,:) = 0.0_wp
      ELSE
        ALLOCATE(rcoeff(0:nend0,nstart0:nend0,ndim), STAT=istat)
        IF(istat /= SUCCESS) CALL finish(modname, 'Allocation of rcoeff failed')
        ALLOCATE(sigma(nstart0:nend0,nscales), STAT=istat)
        IF(istat /= SUCCESS) CALL finish(modname, 'Allocation of gcoeff failed')
        rcoeff(:,:,:) = 0.0_wp
        sigma(:,:) = 0.0_wp
      ENDIF
      IF (msg_level > 0) CALL message(modname,'spectral coefficients allocated')
    END IF


    DO i=1,nspg      ! loop over patterns in grid point space
      DO j=1,nscales   ! loop over scales per pattern

        jj = j + (i-1)*nscales

        ! calculate random normals
        CALL get_complex_random_normals(nlen(j),random_normals)

        ! initialization of spectral coefficients
        k = 1
        IF (lfourier) THEN
          DO n=nstart(j),nend(j)
            DO m=nstart(j),nend(j)

              ! g(n,m) as given by Eq. (A4) of Berner et al (2015, MWR)
              ! see also Thompson et al (2021, MWR), their Equations (1)-(3)
              gcoeff(m,n,jj) = fzero(j)*EXP(-pi4*kappaT(j)*(n*n+m*m))

              ! here we adopt Eq. (19) of Palmer et al.
              rcoeff(m,n,jj) = 1.0_sp/SQRT(1.0_sp - phi(j)**2) * gcoeff(m,n,j) * random_normals(k)

              k = k+1
            END DO
          END DO
          kmax = k-1
        ELSE
          DO n=nstart(j),nend(j)
            sigma(n,j) = sigmacoeff(n,j)
            DO m=0,n

              ! here we use Eq. (19) of Palmer et al.
              rcoeff(m,n,jj) = 1.0_sp/SQRT(1.0_sp - phi(j)**2) * sigma(n,j) * random_normals(k)

              k = k+1
            END DO
          END DO
          kmax = k-1
        END IF

      END DO
    END DO

    WRITE (txt,'(A,L1)') 'initialization complete, spg_fourier_mode = ',lfourier
    CALL message(modname,txt)
    IF (msg_level > 5) THEN
      IF (lfourier) THEN
        WRITE (txt,'(A,f10.2)') '   lat_max = ',lat_max*rad2deg ; CALL message('   ',txt)
        WRITE (txt,'(A,f10.2)') '   lat_min = ',lat_min*rad2deg ; CALL message('   ',txt)
        WRITE (txt,'(A,f10.2)') '   lon_max = ',lon_max*rad2deg ; CALL message('   ',txt)
        WRITE (txt,'(A,f10.2)') '   lon_min = ',lon_min*rad2deg ; CALL message('   ',txt)
      END IF
      WRITE (txt,'(A,i10)')   '   nspg    = ',nspg    ; CALL message('   ',txt)
      WRITE (txt,'(A,i10)')   '   nscales = ',nscales ; CALL message('   ',txt)
      DO j=1,nscales   ! loop over scales
        WRITE (txt,'(A,I4)')    '   scale j = ',j         ; CALL message('   ',txt)
        WRITE (txt,'(A,f10.1)') '   tau_spg = ',tau_spg(j); CALL message('   ',txt)
        WRITE (txt,'(A,e10.3)') '   kappaT  = ',kappaT(j) ; CALL message('   ',txt)
        WRITE (txt,'(A,e10.3)') '   plength = ',plength(j); CALL message('   ',txt)
        WRITE (txt,'(A,f10.2)') '   var_r   = ',var_r(j)  ; CALL message('   ',txt)
        WRITE (txt,'(A,i10)')   '   nstart  = ',nstart(j) ; CALL message('   ',txt)
        WRITE (txt,'(A,i10)')   '   nend    = ',nend(j)   ; CALL message('   ',txt)
        WRITE (txt,'(A,i10)')   '   nlen    = ',nlen(j)   ; CALL message('   ',txt)
      END DO
      WRITE (txt,'(A,i10)')   '   iseed   = ',iseed  ; CALL message('   ',txt)
      WRITE (txt,'(A,i10)')   '   nmax    = ',nmax   ; CALL message('   ',txt)
      WRITE (txt,'(A,i10)')   '   kmax    = ',kmax   ; CALL message('   ',txt)
      WRITE (txt,'(A,i10)')   '   nstep   = ',nstep  ; CALL message('   ',txt)
      WRITE (txt,'(A,f10.1)') '   dtime   = ',dtime  ; CALL message('   ',txt)
      WRITE (txt,'(A,e10.3)') '   dlength = ',dlength ; CALL message('   ',txt)
    END IF

  END SUBROUTINE stochastic_pattern_init

  !==============================================================================
  !  AR1 process for the time evolution of the spectral coefficients
  !------------------------------------------------------------------------------

  SUBROUTINE stochastic_pattern_step(dtime)
    REAL(wp), INTENT(IN) ::  dtime          ! time step
    REAL(KIND=sp) :: phi
    INTEGER       :: n,m,i,j,jj,k
    COMPLEX(KIND=sp), DIMENSION(nmaxwave*nmaxwave) ::  random_normals

    ntime = ntime+1

    IF (MOD(ntime,nstep).EQ.0) THEN

      IF (msg_level > 15) THEN
        WRITE (txt,'(A,i8)')  'step, ntime = ',ntime  ; CALL message(modname,txt)
      END IF

      DO i=1,nspg
        DO j=1,nscales

          jj = j + (i-1)*nscales

          ! one-step correlation of AR1 process
          phi = EXP( -dtime/tau_spg(j)*nstep )

          ! calculate new random normals
          CALL get_complex_random_normals(nlen(j),random_normals)

          ! time stepping of AR1 process for all spectral coefficients
          k = 1
          IF (lfourier) THEN
            DO n=nstart(j),nend(j)
              DO m=nstart(j),nend(j)
                rcoeff(m,n,jj) = phi*rcoeff(m,n,jj) + gcoeff(m,n,jj) * random_normals(k)
                k = k+1
              END DO
            END DO
          ELSE
            DO n=nstart(j),nend(j)
              DO m=0,n
                rcoeff(m,n,jj) = phi*rcoeff(m,n,jj) + sigma(n,j) * random_normals(k)
                k = k+1
              END DO
            END DO
          END IF
        END DO
      END DO

    END IF

  END SUBROUTINE stochastic_pattern_step

  !==============================================================================

  ! Stochastic pattern generator using spherical harmonics (Palmer et al. 2009)
  SUBROUTINE stochastic_pattern_generator( &
       nproma, istart, iend,       & ! indices
       spg,                        & ! spatial random patter
       clat,                       & ! latitude
       clon)                         ! longitude
    INTEGER,  INTENT(IN)        ::  nproma, istart, iend
    REAL(wp), INTENT(INOUT)     ::  spg(:,:)          ! spatial random pattern
    REAL(wp), INTENT(IN)        ::  clat(:)           ! center latitude
    REAL(wp), INTENT(IN)        ::  clon(:)           ! center longitude

    ! local variables (only on VE to avoid warning from VH compiler)
    INTEGER  :: js, jk, jj, jc, nn, mm, ierr

    ! Legendre polynomials need double precision
    REAL(sp), DIMENSION(nproma)              :: xlat, xlon
    REAL(wp), DIMENSION(nproma)              :: xvec
    REAL(wp), DIMENSION(nproma,nmaxwave+1)   :: plg

#ifdef __ASL__
    ! Work array for WINPLG library function
    REAL(wp), DIMENSION(3*nproma+nmaxwave+1) :: pwork
#else
    ! Work array for SHTOOLS legendre polynomials
    REAL(wp), DIMENSION((nmaxwave+1)*(nmaxwave+2)/2) :: plm
#endif

#ifndef __NVCOMPILER_FORTRAN__
    IF (MOD(ntime,nstep).EQ.0) THEN

      IF (msg_level > 15) THEN
        WRITE (txt,'(A,i8)')  'spg,  ntime = ',ntime
        CALL message(modname,txt)
      END IF

      IF (lfourier) THEN

        ! ICON domain is mapped into a double-periodic spectral space
        xlat(:) = pi2*REAL((clat(:)-lat_ctr)/lat_len, kind=sp)
        xlon(:) = pi2*REAL((clon(:)-lon_ctr)/lon_len, kind=sp)

        spg(:,:) = 0.0_wp
        DO jk=1,nspg
          DO js=1,nscales
            jj = js + (jk-1)*nscales
            DO nn=nstart(js),nend(js)
              DO mm=nstart(js),nend(js)
                IF ( rcoeff(mm,nn,jj)%re > 1e-16_sp .OR. rcoeff(mm,nn,jj)%im > 1e-16_sp ) THEN
                  DO jc=istart,iend
                    spg(jc,jk) = spg(jc,jk) + REAL( rcoeff(mm,nn,jj)%re * COS(mm*xlon(jc)) * COS(nn*xlat(jc)) &
                         &                        - rcoeff(mm,nn,jj)%im * SIN(mm*xlon(jc)) * SIN(nn*xlat(jc)), kind=wp)
                  END DO
                END IF
              END DO
            END DO
          END DO
        END DO

      ELSE

        ! global with spherical harmonics and clon is [-pi,pi] and clat is [-pi/2,pi/2],
        ! Legendre polynomials are defined on [-1,1]. Hence, we simply divide clat by pi/2
        xvec(:) = clat(:)*rpi_2
        xlon(:) = REAL( clon(:), kind=sp)

        spg(:,:) = 0.0_wp

        IF ( itype_legendre_polys == 1 ) THEN
          IF (msg_level > 15 .or. ntime < 2) &
               CALL message(modname,'spectral pattern with winplg')
#ifdef __ASL__
          DO jk=1,nspg
            DO js=1,nscales
              jj = js + (jk-1)*nscales
              DO nn=nstart(js),nend(js)
                ! calculates Legendre polynomials for this nn and all mm
                CALL winplg(nproma,xvec,nn,plg,nproma,pwork,ierr)
                DO mm=0,nn
                  DO jc=istart,iend
                    ! we need only the real part of spherical harmonics
                    spg(jc,jk) = spg(jc,jk) + plg(jc,mm+1) * REAL( ( rcoeff(mm,nn,jj)%re * COS(mm*xlon(jc)) &
                         &                                         - rcoeff(mm,nn,jj)%im * SIN(mm*xlon(jc)) ), kind=wp)
                  END DO
                END DO
              END DO
            END DO
          END DO
#endif
        ELSE
          IF (msg_level > 25) &  ! within OpenMP loop
             CALL message(modname,'spectral pattern with PlmON')
#ifndef __ASL__
          DO jk=1,nspg
            DO js=1,nscales
              jj = js + (jk-1)*nscales
              DO jc=istart,iend
                CALL PlmON(plm,nend(js),xvec(jc),1,0,ierr)
                DO nn=nstart(js),nend(js)
                  DO mm=0,nn
                    ! we need only the real part of spherical harmonics
                    spg(jc,jk) = spg(jc,jk) + plm(PlmIndex(nn,mm)) * REAL( ( rcoeff(mm,nn,jj)%re * COS(mm*xlon(jc)) &
                         &                                                 - rcoeff(mm,nn,jj)%im * SIN(mm*xlon(jc)) ), kind=wp)
                  END DO
                END DO
              END DO
            END DO
          END DO
#endif
        END IF

      END IF
    END IF
#endif

  END SUBROUTINE stochastic_pattern_generator

  !==============================================================================

  ! This subroutine fills crnd(:) with new complex random normals
  SUBROUTINE get_complex_random_normals(ndim,crnd)
    INTEGER,     INTENT(in)                   :: ndim
    COMPLEX(sp), INTENT(out), DIMENSION(ndim) :: crnd
    REAL(sp),    DIMENSION(2*ndim)            :: rnd

    CALL get_random_normals(2*ndim,rnd)

    crnd = CMPLX( rnd(1:ndim), rnd(ndim+1:2*ndim), kind=sp)

  END SUBROUTINE get_complex_random_normals

  ! This subroutine fills random_normals(:) with new random normals
  SUBROUTINE get_random_normals(ndim,random_normals)
    INTEGER,  INTENT(in)                   :: ndim
    REAL(sp), INTENT(out), DIMENSION(ndim) :: random_normals

    IF (itype_random_normals == 0) THEN  ! Box-Muller algorithm

      ! get ndim random numbers
      CALL random_NUMBER(random_normals)

      ! transform them to random normals
      CALL random_normal_bm(random_normals)

    ELSEIF (itype_random_normals == 1) THEN

#ifdef __ASL__
      CALL asl_random_generate_s(randomhandle, ndim, random_normals)
#endif
    END IF

  END SUBROUTINE get_random_normals

  SUBROUTINE random_normal_bm(rnd)
    ! Box & Muller, 1958. "A Note on the Generation of Random Normal Deviates".
    ! Annals of Mathematical Statistics, 29, 610-611, doi:10.1214/aoms/1177706645
    REAL(sp), INTENT(inout), DIMENSION(:) :: rnd
    REAL(sp) :: u(SIZE(rnd))
    INTEGER  :: m

    DO m=1,SIZE(rnd)-1,2
       u(m)   = SQRT(-2.0_sp*LOG(rnd(m)))*COS(2.0_sp*pi1*rnd(m+1))
       u(m+1) = SQRT(-2.0_sp*LOG(rnd(m)))*SIN(2.0_sp*pi1*rnd(m+1))
    END DO
    rnd(:) = u(:)

  END SUBROUTINE random_normal_bm

  ! This subroutine will initialize the random number generator
  SUBROUTINE random_initialize(nseed,number)
    INTEGER, INTENT(in)  :: nseed, number
    INTEGER     :: seed(nseed)
    INTEGER     :: n, ierr
    INTEGER(i8) :: i

    DO n=1,nseed
       i = INT(n*13*number,i8)
       seed(n) = myseed(i)
    END DO

    IF (itype_random_normals == 1) THEN

      IF (msg_level > 0) CALL message(modname,'Initalize ASL random number generator')

#ifdef __ASL__

      IF (.not.asl_library_is_initialized()) THEN
        ! initialize ASL libary
        CALL asl_library_initialize()
      END IF

      ! initialize ASL random number generator
      CALL asl_random_create(randomhandle, 0, ierr)
      CALL asl_random_initialize(randomhandle, nseed, seed, ierr)

      ! initalize random normals with mean of zero and stddev of one
      CALL asl_random_distribute_normal(randomhandle, 0.0_dp, 1.0_dp, ierr)
#endif

    ELSE

      IF (msg_level > 0) CALL message(modname,'initalize built-in random number generator')

      ! initialize built-in random number generator
      CALL random_seed(put=seed)

    END IF

  END SUBROUTINE random_initialize

  !==============================================================================

  ELEMENTAL INTEGER FUNCTION myseed(seed)
    INTEGER(i8), INTENT(in) :: seed
    INTEGER(i8) :: s
    IF (seed == 0) THEN
       s = 713105
    ELSE
       s = MOD(seed, 4159265358_i8)
    END IF
    s = MOD(s * 238462643_i8, 4197169399_i8)
    myseed = INT(MOD(s, INT(HUGE(0), i8)), KIND(0))
  END FUNCTION myseed

  !==============================================================================
  ! This subroutine determines the boundaries in lat/lon coordinates of the
  ! outermost ICON domain. This is needed for the stochastic pattern generator
  ! in limited-area mode.

  SUBROUTINE stochastic_pattern_boundaries(p_patch)
    TYPE(t_patch),  INTENT(IN)      :: p_patch        ! patches

    ! Local variables
    REAL(wp), DIMENSION(p_patch%nblks_c) :: clatmax,  clonmax,  clatmin,  clonmin
    REAL(wp)                             :: clatmaxi, clonmaxi, clatmini, clonmini

    ! loop indices
    INTEGER :: jc,jb,jg
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblock, i_endblock  !> blocks
    INTEGER :: i_startidx, i_endidx      !< slices

    IF (msg_level > 0) CALL message(modname,'domain boundaries')

    jg   = p_patch%id

    ! here we want only the outermost domain boundaries
    IF ( jg > 1 ) RETURN

    ! Exclude the nest boundary zone (is this actually necessary for jg=1?)
    rl_start     = grf_bdywidth_c+1
    rl_end       = min_rlcell_int
    i_startblock = p_patch%cells%start_block(rl_start)
    i_endblock   = p_patch%cells%end_block(rl_end)

    ! Find local min/max
    clonmax = -10.0_wp
    clonmin =  10.0_wp
    clatmax = -10.0_wp
    clatmin =  10.0_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblock, i_endblock
      CALL get_indices_c(p_patch, jb, i_startblock, i_endblock, &
                         i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        clatmax(jb)  = MAX(clatmax(jb), p_patch%cells%center(jc,jb)%lat)
        clatmin(jb)  = MIN(clatmin(jb), p_patch%cells%center(jc,jb)%lat)
        clonmax(jb)  = MAX(clonmax(jb), p_patch%cells%center(jc,jb)%lon)
        clonmin(jb)  = MIN(clonmin(jb), p_patch%cells%center(jc,jb)%lon)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! maximum/minimum over blocks
    clatmaxi = MAXVAL(clatmax(i_startblock:i_endblock))
    clatmini = MINVAL(clatmin(i_startblock:i_endblock))
    clonmaxi = MAXVAL(clonmax(i_startblock:i_endblock))
    clonmini = MINVAL(clonmin(i_startblock:i_endblock))

    ! maximum/minimum over all PEs
    lat_max = global_max(clatmaxi)
    lat_min = global_min(clatmini)
    lon_max = global_max(clonmaxi)
    lon_min = global_min(clonmini)

  END SUBROUTINE stochastic_pattern_boundaries

  !==============================================================================

  SUBROUTINE stochastic_pattern_destruct

    IF (ALLOCATED(rcoeff)) DEALLOCATE(rcoeff)
    IF (ALLOCATED(gcoeff)) DEALLOCATE(gcoeff)
    IF (ALLOCATED(sigma)) DEALLOCATE(sigma)

  END SUBROUTINE stochastic_pattern_destruct

  !==============================================================================
END MODULE mo_stoch_pattern_generator
