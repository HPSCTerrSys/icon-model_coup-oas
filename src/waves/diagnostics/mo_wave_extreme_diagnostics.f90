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

! Computes diagnostic extreme wave parameters in the wave model:
! kurtosis, bfis, qp_goda, steepness, thp_adj, sigma_f, sigma_th, hmaxn, hmax
!
! This version of extreme diagnostics is based on WAM 4.5 + ECMWF 2016
!
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_extreme_diagnostics
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_wave_config,         ONLY: t_wave_config
  USE mo_wave_types,          ONLY: t_wave_diag_dyn, t_wave_diag_out, t_wesd
  USE mo_impl_constants,      ONLY: min_rlcell
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_math_constants,      ONLY: pi, pi2
  USE mo_parallel_config,     ONLY: nproma
  USE mo_kind,                ONLY: wp
  USE mo_wave_constants,      ONLY: EMIN
  USE mo_wave_common_diagnostics,    ONLY: significant_wave_height

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_extreme_diagnostics'

  PUBLIC :: calculate_extreme_diagnostics

CONTAINS
  !>
  !! Calculation of extreme-wave diagnostic parameters
  !!
  SUBROUTINE calculate_extreme_diagnostics(p_patch, wave_config, depth, wesd, diag_dyn, diag_out)

    TYPE(t_patch),                 INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET,   INTENT(IN)    :: wave_config
    REAL(wp),                      INTENT(IN)    :: depth(:,:)    ! water depth
    TYPE(t_wesd),                  INTENT(IN)    :: wesd(:)       ! energy spectral bins
    TYPE(t_wave_diag_dyn),         INTENT(IN)    :: diag_dyn
    TYPE(t_wave_diag_out),         INTENT(INOUT) :: diag_out

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':calculate_extreme_parameters'

    ! calculate peak direction and width of 1D freq. and dir. spectral peaks
    !
    CALL peak_approx(p_patch = p_patch,           &
      &          wave_config = wave_config,       &
      &                 wesd = wesd,              &
      &              thp_adj = diag_out%thp_adj,  & ! OUT
      &              sigma_f = diag_out%sigma_f,  & ! OUT
      &              qp_goda = diag_out%qp_goda,  & ! OUT
      &             sigma_th = diag_out%sigma_th)   ! OUT

    ! calculate steepness, BFI, LH broadbandedness, relative width, kurtosis
    !
    CALL wave_kurtosis(p_patch = p_patch,              &
      &                  depth = depth,                &
      &                  emean = diag_dyn%emean,       &
      &               sigma_th = diag_out%sigma_th,    & ! IN
      &                    tm1 = diag_dyn%tm1,         & ! IN
      &                    tm2 = diag_dyn%tm2,         & ! IN
      &                sigma_f = diag_out%sigma_f,     & ! IN
      &                     kp = diag_out%kp,          & ! IN
      &              steepness = diag_out%steepness,   & ! OUT
      &                nu_f_LH = diag_out%nu_f_LH,     & ! OUT
      &                   bfis = diag_out%bfis,        & ! OUT
      &               relw_fth = diag_out%relw_fth,    & ! OUT
      &               kurtosis = diag_out%kurtosis)      ! OUT

    ! calculate maximum wave height and maximum wave period
    !
    CALL max_wave_height(p_patch = p_patch,       &
      &        wave_config = wave_config,         &
      &          steepness = diag_out%steepness,  &
      &           kurtosis = diag_out%kurtosis,   &
      &            nu_f_LH = diag_out%nu_f_LH,    &
      &                tpp = diag_out%tpp,        &
      &              emean = diag_dyn%emean,      &
      &                tm1 = diag_dyn%tm1,        &
      &              hmaxn = diag_out%hmaxn,      & ! OUT
      &               hmax = diag_out%hmax,       & ! OUT
      &               Tmax = diag_out%Tmax)         ! OUT

  END SUBROUTINE calculate_extreme_diagnostics


  !> Applies a parabolic approximation to estimate spectral peak width in f and theta
  !
  ! Adaptation of WAM 4.5 code
  ! SUBROUTINE PEAK_FREQ
  !
  SUBROUTINE peak_approx(p_patch, wave_config, wesd, thp_adj, qp_goda, sigma_f, sigma_th)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//'approx_peak'

    TYPE(t_patch),       INTENT(IN)         :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN) :: wave_config
    TYPE(t_wesd),        INTENT(IN)         :: wesd(:)           !energy spectral bins
    REAL(wp),            INTENT(INOUT)      :: thp_adj(:,:)      !adjusted peak direction
    REAL(wp),            INTENT(INOUT)      :: qp_goda(:,:)      !Goda peakedness parameter
    REAL(wp),            INTENT(INOUT)      :: sigma_f(:,:)      !width of frequency peak
    REAL(wp),            INTENT(INOUT)      :: sigma_th(:,:)     !width of directional peak

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd,jdm,jdp
    INTEGER :: peak_indf(nproma), peak_indd(nproma)

    REAL(wp):: Amax, Bmax, Cmax, pplus, pminus, Splus, Sminus, sigma_th2, fp_adj
    REAL(wp):: specf(nproma,wave_config%nfreqs), specd(nproma,wave_config%ndirs)
    REAL(wp):: specf_max(nproma), specd_max(nproma), sigma_th1(nproma)
    REAL(wp):: sum0(nproma), fesq(nproma), sum1(nproma), sum2(nproma)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! save some paperwork
    wc => wave_config

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,jd,jdm,jdp,i_startidx,i_endidx,specf,specd,specf_max,sum0,fesq,fp_adj, &
!$OMP & Amax,Bmax,Cmax,pplus,pminus,Splus,Sminus,specd_max,sigma_th1,peak_indf,peak_indd,sum1,sum2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      ! Initialisation of temporary local fields
      !----------------------------------------------------------------------------------------
      DO jc = i_startidx, i_endidx
        specf_max(jc) = 0._wp
        specd_max(jc) = 0._wp
        peak_indd(jc) = 1
        peak_indf(jc) = 1
        sum0(jc)      = 0._wp
        fesq(jc)      = 0._wp
        sum1(jc)      = 0._wp
        sum2(jc)      = 0._wp
        sigma_th1(jc) = 1._wp
      END DO

      ! Estimate Goda parameter and spectral width in frequency
      !----------------------------------------------------------------------------------------
      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          specf(jc,jf)   = 0._wp
        END DO

        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx
            specf(jc,jf) = specf(jc,jf) + wesd(jf)%ptr(jc,jd,jb)*wc%DELTH
          END DO
        END DO  ! jd

        DO jc = i_startidx, i_endidx
          IF (specf(jc,jf) > specf_max(jc)) THEN
            specf_max(jc) = specf(jc,jf)
            peak_indf(jc) = jf
          END IF
        END DO
      END DO  ! jf

      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          IF (specf(jc,jf) > 0.4_wp*specf_max(jc)) THEN
            sum0(jc) = sum0(jc) + specf(jc,jf)*wc%dfreqs(jf)
            fesq(jc) = fesq(jc) + 2._wp*specf(jc,jf)**2 * wc%dfreqs_freqs(jf)
          END IF
        END DO
      END DO  ! jf

      DO jc = i_startidx, i_endidx
        IF (sum0(jc) .ge. 10._wp*EMIN) THEN
          qp_goda(jc,jb) = MAX(MIN(fesq(jc)/sum0(jc)**2, 15._wp), 0.5_wp)
          IF (peak_indf(jc) .GE. wc%nfreqs .or. peak_indf(jc) .LE. 1) CYCLE
          fp_adj = wc%freqs(peak_indf(jc))
          Amax   = specf_max(jc)
          pplus  = wc%freqs(peak_indf(jc)+1)-fp_adj
          pminus = wc%freqs(peak_indf(jc)-1)-fp_adj
          Splus  = ( specf(jc,peak_indf(jc)+1)-Amax )/pplus
          Sminus = ( specf(jc,peak_indf(jc)-1)-Amax )/pminus
          Cmax   = (Sminus-Splus)/(pminus-pplus)
          IF (Cmax .GE. 0._wp) CYCLE
          Bmax = (pminus*Splus-Splus*pminus)/(pminus-pplus)
          fp_adj  = fp_adj - Bmax/(2._wp*Cmax)
          specf_max(jc)  = Amax - Bmax**2/(4._wp*Cmax)
          sigma_f(jc,jb) = MIN( SQRT(-specf_max(jc)/(2._wp*Cmax))/fp_adj, 1._wp/(SQRT(pi)*qp_goda(jc,jb)) )
        ELSE
          qp_goda(jc,jb) = 0._wp
          sigma_f(jc,jb) = 0._wp
        END IF
      END DO

      ! Estimate peak direction and spectral width in angle
      !----------------------------------------------------------------------------------------
      DO jd = 1,wc%ndirs
        DO jc = i_startidx, i_endidx
          specd(jc,jd)   = 0._wp
        END DO

        DO jf = 1,wc%nfreqs
          DO jc = i_startidx, i_endidx
            specd(jc,jd) = specd(jc,jd) + wesd(jf)%ptr(jc,jd,jb)*wc%dfreqs(jf)
          END DO
        END DO  ! jf

        DO jc = i_startidx, i_endidx
          IF (specd(jc,jd) > specd_max(jc)) THEN
            specd_max(jc) = specd(jc,jd)
            peak_indd(jc) = jd
          END IF
        END DO
      END DO  ! jd

      DO jc = i_startidx, i_endidx
        thp_adj(jc,jb) = wc%dirs(peak_indd(jc))
        Amax    = specd_max(jc)
        pplus   =  wc%DELTH
        pminus  = -wc%DELTH
        jdp = peak_indd(jc)+1
        IF (jdp > wc%ndirs) jdp = jdp - wc%ndirs
        jdm = peak_indd(jc)-1
        IF (jdm < 1)        jdm = jdm + wc%ndirs
        Amax    = specd_max(jc)
        Splus   = ( specd(jc,jdp)-Amax )/pplus
        Sminus  = ( specd(jc,jdm)-Amax )/pminus
        Cmax    = (Sminus-Splus)/(pminus-pplus)
        IF (Cmax .GE. 0._wp) CYCLE
        Bmax  = (pminus*Splus-pplus*Sminus)/(pminus-pplus)
        thp_adj(jc,jb) = thp_adj(jc,jb) - Bmax/(2._wp*Cmax)
        specd_max(jc)  = Amax - Bmax**2/(4._wp*Cmax)
        sigma_th1(jc)  = SQRT(-specd_max(jc)/(2._wp*Cmax))/thp_adj(jc,jb)
      END DO !jc

      DO jd = 1,wc%ndirs
        DO jc = i_startidx,i_endidx
        IF (COS(wc%dirs(jd)-thp_adj(jc,jb)) < 0.5_wp) CYCLE
          sum1(jc) = sum1(jc) + specd(jc,jd)
          sum2(jc) = sum2(jc) + COS(wc%dirs(jd)-thp_adj(jc,jb))*specd(jc,jd)
        END DO !jc
      END DO  !jd

      DO jc = i_startidx,i_endidx
        sigma_th2 = MERGE(1._wp, SQRT(2._wp*(1._wp-sum2(jc)/sum1(jc))), sum1(jc) .le. 0._wp)
        sigma_th(jc,jb) = MIN(sigma_th1(jc),sigma_th2)
      END DO

    END DO !jb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE peak_approx

  !<
  ! calculation of extreme wave parameters: slope, spectral widths, Goda peakedness, BFI, kurtosis
  !
  ! Adapted from WAM 4.5 model
  ! KURTOSIS, TRANSF, TRANSF2
  ! References:
  ! P. Janssen & J.-R. Bidlot (2008), P. Janssen & M. Onorato (2007), ECMWF IFS documentation (2016)
  !
  SUBROUTINE wave_kurtosis(p_patch, depth, emean, sigma_th, tm1, tm2, sigma_f, &
                          & kp, steepness, nu_f_LH, bfis, relw_fth, kurtosis)

  CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':wave_kurtosis'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    REAL(wp),                    INTENT(IN)    :: depth(:,:)      !< local ocean depth [m]
    REAL(wp),                    INTENT(IN)    :: emean(:,:)      !< total energy [m^2]
    REAL(wp),                    INTENT(IN)    :: sigma_th(:,:)   !< directional bandwidth [s-1]
    REAL(wp),                    INTENT(IN)    :: tm1(:,:)        !< m1 wave period [s]
    REAL(wp),                    INTENT(IN)    :: tm2(:,:)        !< m2 wave period [s]
    REAL(wp),                    INTENT(IN)    :: sigma_f(:,:)    !< frequency bandwidth [s-1]
    REAL(wp),                    INTENT(IN)    :: kp(:,:)         !< peak wavenumber [m-1]
    REAL(wp),                    INTENT(INOUT) :: steepness(:,:)  !< wave steepness [-]
    REAL(wp),                    INTENT(INOUT) :: nu_f_LH(:,:)    !< broadbandness parameter [-]
    REAL(wp),                    INTENT(INOUT) :: bfis(:,:)       !< Benjamin-Feir index [-]
    REAL(wp),                    INTENT(INOUT) :: relw_fth(:,:)   !< ratio of dir. peak to freq. peak [-]
    REAL(wp),                    INTENT(INOUT) :: kurtosis(:,:)   !< spectral kurtosis [-]

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb

    REAL(wp), PARAMETER :: pi3s3 = pi/(3._wp*SQRT(3._wp))
    REAL(wp), PARAMETER :: sfmin = 1._wp/(SQRT(pi)*15._wp)
    REAL(wp):: akpd, Th, Crel, Xnlfac, Xnl, c4_tot, bfidw

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,bfidw,akpd,Crel,Th,Xnlfac,Xnl,c4_tot) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx,i_endidx
        steepness(jc,jb) = kp(jc,jb)*SQRT(emean(jc,jb))
        nu_f_LH(jc,jb)   = SQRT((tm1(jc,jb)/tm2(jc,jb))**2-1)
        bfidw            = SQRT(2._wp)*steepness(jc,jb)/MAX(sigma_f(jc,jb),sfmin)

        !----------------------------------------------------------------------------
        ! Parameters for calculating the intermediate-depth BFI and the kurtosis
        !----------------------------------------------------------------------------
        ! Xnl = (9Th^4-10Th^2+9)/(8Th^3) - (1 + (2cg-cp/2)^2/(cs^2-cg^2))/(kd)
        ! BFIs^2 = -BFIdw^2 (cg/cp)^2 g X_nl/(k \omega \omega'')
        ! C4 = kappa_4/8 - pi/(3\sqrt(3)) BFIs^2/sqrt(1 + 3.5(\sigma_th/\sigma_f)^2)
        !----------------------------------------------------------------------------
        !
        ! constants depending on kd
        !
        akpd   = kp(jc,jb)*depth(jc,jb)
        Th     = TANH(akpd)
        Crel   = (1._wp + akpd*(1._wp-Th**2)/Th)/2._wp
        Xnlfac = 4._wp*Th*Crel**2/( (Th-akpd*(1._wp-Th**2))**2 +(1._wp-Th**2)*(2._wp*akpd*Th)**2 )
        Xnl    = ( 9._wp*(1._wp+Th**4)-10._wp*Th**2 )/(8._wp*Th**3) -  &
               & ( 1._wp + (4._wp*Crel-1._wp)**2/(4._wp*akpd/Th-4._wp*Crel**2) )/akpd
        !
        ! output variables
        !
        bfis(jc,jb)     = SQRT(MAX(0._wp,Xnlfac*Xnl))*bfidw
        relw_fth(jc,jb) = MERGE(0._wp, 0.5_wp*(sigma_th(jc,jb)/nu_f_LH(jc,jb))**2, bfidw<0.0001_wp)
        c4_tot          = 3._wp*steepness(jc,jb)**2 +pi3s3*bfis(jc,jb)**2/SQRT(1+7._wp*relw_fth(jc,jb))
        kurtosis(jc,jb) = MAX(-0.33_wp, MIN(1._wp, c4_tot) )
      END DO
    END DO

!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE wave_kurtosis

  !>
  !! Calculation of normalised and bare maximum wave height
  !! based on WAM 4.5 formulation
  !!
  !! Value for skewness and bound-wave kurtosis are from ECMWF IFS 2016
  !!
  SUBROUTINE max_wave_height(p_patch, wave_config, steepness, kurtosis,    &
           &  nu_f_LH, tpp, emean, tm1, hmaxn, hmax, Tmax)

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: steepness(:,:)
    REAL(wp),                    INTENT(IN)    :: kurtosis(:,:)
    REAL(wp),                    INTENT(IN)    :: nu_f_LH(:,:)
    REAL(wp),                    INTENT(IN)    :: tpp(:,:)
    REAL(wp),                    INTENT(IN)    :: emean(:,:)
    REAL(wp),                    INTENT(IN)    :: tm1(:,:)
    REAL(wp),                    INTENT(INOUT) :: hmaxn(:,:)    ! OUT
    REAL(wp),                    INTENT(INOUT) :: hmax(:,:)     ! OUT
    REAL(wp),                    INTENT(INOUT) :: Tmax(:,:)     ! OUT

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':max_wave_height'

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb
    REAL(wp):: alpha, beta, z0h, argh, zh, c3, zeps
    REAL(wp):: hs(SIZE(emean,1),SIZE(emean,2))   ! significant wave height

    REAL(wp), PARAMETER :: sq_2pi = 2._wp/SQRT(pi2)
    REAL(wp), PARAMETER :: gammaE = 0.57721566_wp ! Euler gamma constant
    REAL(wp), PARAMETER :: zeta_3 = 1.20205690_wp ! Riemann Zeta function at z=3
    ! G1, G2, G3 are 1st, 2nd and 3rd order derivatives of Gamma function at z=1.
    REAL(wp), PARAMETER :: G1 = -gammaE
    REAL(wp), PARAMETER :: G2 = (G1**2+pi**2/6._wp)
    REAL(wp), PARAMETER :: G3 = -2._wp*zeta_3+G1**3+0.5_wp*G1*pi**2


    CALL significant_wave_height(p_patch = p_patch,    &
      &                          emean   = emean(:,:), &
      &                          hs      = hs(:,:))    ! OUT

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,c3,z0h,alpha,beta,argh,zh,zeps) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        c3    = steepness(jc,jb)*SQRT(0.75_wp)
        z0h   = 0.5_wp*LOG(sq_2pi*nu_f_LH(jc,jb)*wave_config%Tlength/tpp(jc,jb))
        alpha = 2._wp*z0h*(z0h-1._wp) + (1._wp-2._wp*z0h)*G1 + 0.5_wp*G2
        beta  = z0h*(4._wp*z0h**2-12._wp*z0h+6._wp) - (6._wp*z0h**2-12._wp*z0h+3._wp)*G1 &
              & + 3._wp*(z0h-1._wp)*G2 - 0.5_wp*G3
        argh  = MAX(0.1_wp, 1._wp+alpha*kurtosis(jc,jb)+beta*c3**2)
        zh    = MERGE(1._wp, z0h + 0.5_wp*(-G1+LOG(argh)), c3<SQRT(EMIN))
        hmaxn(jc,jb) = SQRT(zh)
        hmax(jc,jb)  = hmaxn(jc,jb)*hs(jc,jb)
        zeps = nu_f_LH(jc,jb)/( SQRT(2.0_wp)*hmaxn(jc,jb) )
        Tmax(jc,jb) = ( 1._wp + 0.5_wp*zeps**2 + 0.75_wp*zeps**4 )*tm1(jc,jb)
      END DO
    END DO

!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE max_wave_height


END MODULE mo_wave_extreme_diagnostics
