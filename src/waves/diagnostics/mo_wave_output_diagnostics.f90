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

! Computes output diagnostics for scientific analysis and visualization.
! Fields have no impact on the wave solution.
! updated at output time steps

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_output_diagnostics
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_wave_config,         ONLY: t_wave_config
  USE mo_wave_types,          ONLY: t_wave_diag_dyn, t_wave_diag_out, t_wesd
  USE mo_impl_constants,      ONLY: min_rlcell
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_physical_constants,  ONLY: grav
  USE mo_math_constants,      ONLY: pi2, rad2deg
  USE mo_parallel_config,     ONLY: nproma
  USE mo_fortran_tools,       ONLY: init
  USE mo_wave_constants,      ONLY: EMIN
  USE mo_wave_common_diagnostics, ONLY: significant_wave_height

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_output_diagnostics'

  PUBLIC :: calculate_output_diagnostics

CONTAINS
  !>
  !! Calculation of purely diagnostic parameters
  !!
  SUBROUTINE calculate_output_diagnostics(p_patch, wave_config, sp10m, dir10m, wesd, diag_dyn, diag_out)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_config),   INTENT(IN)    :: wave_config
    REAL(wp),              INTENT(IN)    :: sp10m(:,:)
    REAL(wp),              INTENT(IN)    :: dir10m(:,:)   ! wind direction in 10m [rad]
    TYPE(t_wesd),          INTENT(IN)    :: wesd(:)       ! energy spectral bins
    TYPE(t_wave_diag_dyn), INTENT(IN)    :: diag_dyn
    TYPE(t_wave_diag_out), INTENT(INOUT) :: diag_out

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':calculate_output_diagnostics'

    ! calculate swell separation mask
    !
    CALL swell_separation(p_patch = p_patch, &
      &               wave_config = wave_config, &
      &                    dir10m = dir10m, &
      &                     ustar = diag_dyn%ustar, &
      &                swell_mask = diag_out%swell_mask) ! OUT

    ! separate total energy according to swell mask
    !
    CALL total_energy_sep(p_patch = p_patch, &
      &               wave_config = wave_config, &
      &                      wesd = wesd, &
      &                      mask = diag_out%swell_mask, &
      &                   emeanws = diag_out%emean_sea,  & ! OUT
      &                    emeans = diag_out%emean_swell)  ! OUT

    ! separate mean frequency energy according to swell mask
    !
    CALL mean_frequency_energy_sep(p_patch = p_patch, &
      &               wave_config = wave_config, &
      &                      wesd = wesd, &
      &                      mask = diag_out%swell_mask, &
      &                   emeanws = diag_out%emean_sea, &
      &                    emeans = diag_out%emean_swell, &
      &                  femeanws = diag_out%femean_sea, & ! OUT
      &                   femeans = diag_out%femean_swell) ! OUT


    ! calculate significant wave height from sea wave energy
    !
    CALL significant_wave_height(p_patch = p_patch, &
      &                          emean   = diag_out%emean_sea(:,:), &
      &                          hs      = diag_out%hs_sea(:,:)) ! OUT

    ! calculate significant wave height from swell wave energy
    !
    CALL significant_wave_height(p_patch = p_patch, &
      &                          emean   = diag_out%emean_swell(:,:), &
      &                          hs      = diag_out%hs_swell(:,:)) ! OUT

    ! calculate mean wave direction and directional wave spread
    !
    CALL mean_wave_direction_spread(p_patch = p_patch, &
      &                         wave_config = wave_config, &
      &                                wesd = wesd, &
      &                            mean_dir = diag_out%hs_dir, & ! OUT
      &                         mean_spread = diag_out%ds)       ! OUT

    ! calculate mean wave direction and directional wave spread
    ! separately for sea and swell according to given swell mask.
    !
    CALL mean_wave_direction_spread_sep(p_patch = p_patch, &
      &                             wave_config = wave_config, &
      &                                    wesd = wesd, &
      &                                    mask = diag_out%swell_mask, &
      &                                  md_sea = diag_out%hs_sea_dir, &   ! OUT
      &                                  ms_sea = diag_out%ds_sea, &       ! OUT
      &                                md_swell = diag_out%hs_swell_dir, & ! OUT
      &                                ms_swell = diag_out%ds_swell)       ! OUT

    ! calculate mean wave period from mean frequency wave energy
    !
    CALL mean_wave_period(p_patch = p_patch, &
      &                    femean = diag_dyn%femean(:,:), &
      &                        mp = diag_out%tmp(:,:)) ! OUT

    ! calculate mean wave period from wind sea mean frequency wave energy
    !
    CALL mean_wave_period(p_patch = p_patch, &
      &                    femean = diag_out%femean_sea(:,:), &
      &                        mp = diag_out%mp_sea(:,:)) ! OUT

    ! calculate mean wave period from wind swell mean frequency wave energy
    !
    CALL mean_wave_period(p_patch = p_patch, &
      &                    femean = diag_out%femean_swell(:,:), &
      &                        mp = diag_out%mp_swell(:,:)) ! OUT

    ! separate M1 and M2 periods according to swell mask
    !
    CALL m1_m2_periods_sep(p_patch = p_patch,              &
      &                wave_config = wave_config,          &
      &                       wesd = wesd,                 &
      &                       mask = diag_out%swell_mask,  &
      &                    emeanws = diag_out%emean_sea,   &
      &                     emeans = diag_out%emean_swell, &
      &                       m1ws = diag_out%m1_sea,      & ! OUT
      &                        m1s = diag_out%m1_swell,    & ! OUT
      &                       m2ws = diag_out%m2_sea,      & ! OUT
      &                        m2s = diag_out%m2_swell,    & ! OUT
      &                   f1meanws = diag_out%f1mean_sea,  & ! OUT
      &                    f1means = diag_out%f1mean_swell)  ! OUT

    ! calculate peak wave period and wavenumber
    !
    CALL peak_wave_period_wavenumber(p_patch = p_patch, &
      &               wave_config = wave_config,        &
      &                      wesd = wesd,               &
      &                wave_num_c = diag_dyn%wave_num_c,  &
      &                        pp = diag_out%tpp,       & ! OUT
      &                        kp = diag_out%kp)          ! OUT

    ! Calculate peak wind sea and swell wave period and wavenumber according to swell mask
    !
    CALL peak_wave_period_wavenumber_sep(p_patch = p_patch,  &
      &                   wave_config = wave_config,         &
      &                          wesd = wesd,                &
      &                    wave_num_c = diag_dyn%wave_num_c,   &
      &                          mask = diag_out%swell_mask, &
      &                        pp_sea = diag_out%pp_sea,     & ! OUT
      &                      pp_swell = diag_out%pp_swell,   & ! OUT
      &                        kp_sea = diag_out%kp_sea,     & ! OUT
      &                      kp_swell = diag_out%kp_swell)     ! OUT

    ! calculate wave drag coefficient and normalised stress
    !
    CALL wave_drag_stress_ch_par(p_patch = p_patch, &
      &                     sp10m = sp10m,          &
      &                     ustar = diag_dyn%ustar,   &
      &                      tauw = diag_dyn%tauw,    &
      &                        z0 = diag_dyn%z0,      &
      &                      drag = diag_out%drag,  & ! OUT
      &                     tauwn = diag_out%tauwn, & ! OUT
      &                      beta = diag_out%beta)    ! OUT

  END SUBROUTINE calculate_output_diagnostics

  !>
  !! Separation of M1 and M2 periods according to swell mask
  !!
  !! Adaptation of WAM 4.5 code.
  !! TM1_TM2_PERIODS_B
  !! Integration of spectra and adding of tail factors.
  !!
  SUBROUTINE m1_m2_periods_sep(p_patch, wave_config, wesd, mask, emeanws, emeans, &
    &                          m1ws, m1s, m2ws, m2s, f1meanws, f1means)
    CHARACTER(len=*), PARAMETER ::  &
         &  routine = modname//'m1_m2_periods_sep'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wesd),                INTENT(IN)    :: wesd(:)         !< energy spectral bins (nproma,ndirs,nblks_c,nfreqs)
    INTEGER,                     INTENT(IN)    :: mask(:,:,:,:)   !< swell mask
    REAL(wp),                    INTENT(IN)    :: emeanws(:,:)    !< wind sea energy (nproma,nblks_c)
    REAL(wp),                    INTENT(IN)    :: emeans(:,:)     !< swell energy (nproma,nblks_c)
    REAL(wp),                    INTENT(INOUT) :: m1ws(:,:)       !< wind sea m1 period (nproma,nblks_c)
    REAL(wp),                    INTENT(INOUT) :: m1s(:,:)        !< swell m1 period (nproma,nblks_c)
    REAL(wp),                    INTENT(INOUT) :: m2ws(:,:)       !< wind sea m2 period (nproma,nblks_c)
    REAL(wp),                    INTENT(INOUT) :: m2s(:,:)        !< swell m2 period (nproma,nblks_c)
    REAL(wp),                    INTENT(INOUT) :: f1meanws(:,:)   !< wind sea m1 frequency (nproma,nblks_c)
    REAL(wp),                    INTENT(INOUT) :: f1means(:,:)    !< swell m1 frequency (nproma,nblks_c)

    ! local
    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd
    REAL(wp):: temp(nproma,wave_config%nfreqs), temp1(nproma,wave_config%nfreqs)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! save some paperwork
    wc => wave_config

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,jd,i_startidx,i_endidx,temp,temp1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      ! compute sum of all wesd components that match a specific frequency
      DO jf = 1,wc%nfreqs

        ! initialization
        DO jc = i_startidx, i_endidx
          temp(jc,jf) = 0._wp
          temp1(jc,jf) = 0._wp
        END DO

        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx
            IF (mask(jc,jd,jf,jb)==1) THEN ! belongs to swell
              temp(jc,jf) = temp(jc,jf) + wesd(jf)%ptr(jc,jd,jb)
            ELSE
              temp1(jc,jf) = temp1(jc,jf) + wesd(jf)%ptr(jc,jd,jb)
            END IF
          END DO
        END DO  ! jd

      END DO  ! jf

      ! tail part
      DO jc = i_startidx, i_endidx
        m1s(jc,jb)  = wc%MP1_TAIL * temp(jc,wc%nfreqs)
        m2s(jc,jb)  = wc%MP2_TAIL * temp(jc,wc%nfreqs)
        m1ws(jc,jb) = wc%MP1_TAIL * temp1(jc,wc%nfreqs)
        m2ws(jc,jb) = wc%MP2_TAIL * temp1(jc,wc%nfreqs)
      END DO

      ! add all other frequencies
      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          m1s(jc,jb)  = m1s(jc,jb) + temp(jc,jf) * wc%dfim_fr(jf)
          m2s(jc,jb)  = m2s(jc,jb) + temp(jc,jf) * wc%dfim_fr2(jf)
          m1ws(jc,jb) = m1ws(jc,jb) + temp1(jc,jf) * wc%dfim_fr(jf)
          m2ws(jc,jb) = m2ws(jc,jb) + temp1(jc,jf) * wc%dfim_fr2(jf)
        END DO
      END DO

      ! clipping
      DO jc = i_startidx, i_endidx
        IF (emeans(jc,jb).gt.EMIN) THEN
          m1s(jc,jb) = emeans(jc,jb) / m1s(jc,jb)
          m2s(jc,jb) = SQRT(emeans(jc,jb) / m2s(jc,jb))
        ELSE
          m1s(jc,jb) =  1.0_wp
          m2s(jc,jb) =  1.0_wp
        END IF
        f1means(jc,jb) = 1.0_wp / m1s(jc,jb)

        IF (emeanws(jc,jb).gt.EMIN) THEN
          m1ws(jc,jb) = emeanws(jc,jb) / m1ws(jc,jb)
          m2ws(jc,jb) = SQRT(emeanws(jc,jb) / m2ws(jc,jb))
        ELSE
          m1ws(jc,jb) =  1.0_wp
          m2ws(jc,jb) =  1.0_wp
        END IF
        f1meanws(jc,jb) = 1.0_wp / m1ws(jc,jb)
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE m1_m2_periods_sep


  !>
  !! Separation of mean frequency energy according to mask
  !!
  !! Integration over frequencies for calculation of
  !! of mean frequency energy. Adaptation of WAM 4.5 code
  !! of the subroutine FEMEAN developed by S.D. HASSELMANN,
  !! optimized by L. Zambresky and H. Guenther, GKSS, 2001                              !
  !!
  SUBROUTINE mean_frequency_energy_sep(p_patch, wave_config, wesd, mask, emeanws, emeans, femeanws, femeans)
    CHARACTER(len=*), PARAMETER :: &
         & routine =  modname//'mean_frequency_energy_sep'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wesd), INTENT(IN)    :: wesd(:)         !energy spectral bins
                                                   ! wesd(nfreqs)%(nproma,ndirs,nblks_c)
    INTEGER,      INTENT(IN)    :: mask(:,:,:,:)   !=1 - swell           (nproma,ndirs,nfreqs,nblks_c)
    REAL(wp),     INTENT(IN)    :: emeanws(:,:)    !wind sea energy      (nproma,nblks_c)
    REAL(wp),     INTENT(IN)    :: emeans(:,:)     !swell energy         (nproma,nblks_c)
    REAL(wp),     INTENT(INOUT) :: femeans(:,:)    !swell mean frequency energy    (nproma,nblks_c)
    REAL(wp),     INTENT(INOUT) :: femeanws(:,:)   !wind sea mean frequency energy (nproma,nblks_c)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd

    REAL(wp) :: temp(nproma,wave_config%nfreqs), temp_1(nproma,wave_config%nfreqs)
    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! save some paperwork
    wc => wave_config


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,jd,i_startidx,i_endidx,temp,temp_1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs

        DO jc = i_startidx, i_endidx
          temp(jc,jf)   = 0._wp
          temp_1(jc,jf) = 0._wp
        ENDDO

        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx
            IF (mask(jc,jd,jf,jb)==1) THEN ! belongs to swell
              temp(jc,jf) = temp(jc,jf) + wesd(jf)%ptr(jc,jd,jb)
            ELSE ! belongs to wind sea
              temp_1(jc,jf) = temp_1(jc,jf) + wesd(jf)%ptr(jc,jd,jb)
            ENDIF
          ENDDO
        ENDDO  ! jd

      END DO  ! jf

      DO jc = i_startidx, i_endidx
        femeans(jc,jb)  = wc%MM1_TAIL * temp(jc,wc%nfreqs)
        femeanws(jc,jb) = wc%MM1_TAIL * temp_1(jc,wc%nfreqs)
      END DO

      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          femeans(jc,jb)  = femeans(jc,jb) + temp(jc,jf) * wc%DFIMOFR(jf)
          femeanws(jc,jb) = femeanws(jc,jb) + temp_1(jc,jf) * wc%DFIMOFR(jf)
        END DO
      END DO

      DO jc = i_startidx, i_endidx
        femeans(jc,jb)  = emeans(jc,jb) / MAX(femeans(jc,jb),EMIN)
        femeanws(jc,jb) = emeanws(jc,jb) / MAX(femeanws(jc,jb),EMIN)
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE mean_frequency_energy_sep

  !>
  !! Calculation of separation of total energy according to mask
  !!
  !! Calculation of total energy by integration over directions and frequencies.
  !! A tail correction is added.
  !! Adaptation of WAM 4.5 code of the subroutine TOTAL_ENERGY
  !! developed by S.D. HASSELMANN, optimized by L. Zambresky
  !! and H. Guenther, GKSS, 2001
  !!
  SUBROUTINE total_energy_sep(p_patch, wave_config, wesd, mask, emeanws, emeans)
    CHARACTER(len=*), PARAMETER ::  &
         &  routine = modname//'total_energy_sep'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wesd), INTENT(IN)    :: wesd(:)         !energy spectral bins (nproma,ndirs,nblks_c,nfreqs)
    INTEGER,      INTENT(IN)    :: mask(:,:,:,:)   !=1 where energy belongs to swell emeans, =0 - belongs to wind sea emeanws
    REAL(wp),     INTENT(INOUT) :: emeanws(:,:)    !wind sea energy (nproma,nblks_c)
    REAL(wp),     INTENT(INOUT) :: emeans(:,:)     !swell energy (nproma,nblks_c)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd

    REAL(wp):: sum1(nproma,wave_config%nfreqs), sum2(nproma,wave_config%nfreqs)
    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! save some paperwork
    wc => wave_config


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,jd,i_startidx,i_endidx,sum1,sum2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)


      ! compute sum of all wesd components that match a specific frequency
      DO jf = 1,wc%nfreqs

        ! initialization
        DO jc = i_startidx, i_endidx
          sum1(jc,jf) = 0._wp
          sum2(jc,jf) = 0._wp
        ENDDO

        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx
            IF (mask(jc,jd,jf,jb) == 1) THEN
              sum1(jc,jf) = sum1(jc,jf) + wesd(jf)%ptr(jc,jd,jb)
            ELSE
              sum2(jc,jf) = sum2(jc,jf) + wesd(jf)%ptr(jc,jd,jb)
            END IF
          END DO
        END DO  ! jd

      ENDDO  ! jf

      ! initialization
      DO jc = i_startidx, i_endidx
        emeans(jc,jb)  = wc%MO_TAIL * sum1(jc,wc%nfreqs)
        emeanws(jc,jb) = wc%MO_TAIL * sum2(jc,wc%nfreqs)
      ENDDO

      ! sum over all frequencies
      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          emeans(jc,jb)  = emeans(jc,jb)  + sum1(jc,jf) * wc%DFIM(jf)
          emeanws(jc,jb) = emeanws(jc,jb) + sum2(jc,jf) * wc%DFIM(jf)
        END DO
      ENDDO  ! jf

      ! clipping
      DO jc = i_startidx, i_endidx
        emeans(jc,jb)  = MAX(emeans(jc,jb),EMIN)
        emeanws(jc,jb) = MAX(emeanws(jc,jb),EMIN)
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE total_energy_sep


  !>
  !! Calculation of the swell separation mask
  !!
  !! Based on WAM 4.5 subroutine SWELL_SEPARATION
  !!     P.Lionello    February   87
  !!     L.Zambresky   November   87  GKSS/ECMWF
  !!     J.Bidlot      Febrary  1996       ECMWF
  !!     J.Bidlot      May      1999       ECMWF
  !!     J.Bidlot      April    2000       ECMWF
  !!     J.Bidlot      December 2003       ECMWF
  !!     L. Aouf                2013          MF
  !!
  !! The waves which do not interact with the wind are
  !! considered as swell.
  !!
  SUBROUTINE swell_separation(p_patch, wave_config, dir10m, ustar, swell_mask)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':swell_separation'

    TYPE(t_patch),               INTENT(IN) :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN) :: wave_config
    REAL(wp),                    INTENT(IN) :: dir10m(:,:)
    REAL(wp),                    INTENT(IN) :: ustar(:,:)
    INTEGER,                  INTENT(INOUT) :: swell_mask(:,:,:,:)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd


    TYPE(t_wave_config), POINTER :: wc => NULL()

    REAL(wp) :: fric, dw_phase_vel, trhld

    fric = 28._wp

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! save some paperwork
    wc => wave_config

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jc,i_startidx,i_endidx,dw_phase_vel,trhld) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jf = 1,wc%nfreqs
        dw_phase_vel = grav/(pi2*wc%freqs(jf))
        DO jd = 1, wc%ndirs
          DO jc = i_startidx, i_endidx

            trhld = fric/dw_phase_vel * 1.2_wp*ustar(jc,jb)*COS(wc%dirs(jd) - dir10m(jc,jb))
            swell_mask(jc,jd,jf,jb) = MERGE(1, 0, trhld < 1._wp)

          END DO
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE swell_separation


  !>
  !! Calculation of mean wave period 1/femean
  !! based on WAM 4.5 formulation
  !!
  SUBROUTINE mean_wave_period(p_patch, femean, mp)

    TYPE(t_patch),     INTENT(IN)    :: p_patch
    REAL(wp),          INTENT(IN)    :: femean(:,:) !< mean frequency energy [m^2]
    REAL(wp),          INTENT(INOUT) :: mp(:,:)     !< mean wave period [s]

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':mean_wave_period'

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb


    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx
        mp(jc,jb) = 1.0_wp / femean(jc,jb)
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE mean_wave_period

  !>
  !! Calculation of peak wave period and peak wavenumber
  !! based on WAM 4.5 formulation
  !!
  SUBROUTINE peak_wave_period_wavenumber(p_patch, wave_config, wesd, wave_num_c, pp, kp)
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':peak_wave_period_wavenumber'

    TYPE(t_patch),       INTENT(IN)         :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN) :: wave_config
    TYPE(t_wesd),        INTENT(IN)         :: wesd(:)           !energy spectral bins
    REAL(wp),            INTENT(IN)         :: wave_num_c(:,:,:) !wavenumber (1/m)
    REAL(wp),            INTENT(INOUT)      :: pp(:,:)           !peak period
    REAL(wp),            INTENT(INOUT)      :: kp(:,:)           !peak wavenumber

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd
    INTEGER :: peak_ind(nproma)

    REAL(wp):: temp(nproma,wave_config%nfreqs)
    REAL(wp):: temp_max(nproma)        !< maximum of temp over all frequencies

    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! save some paperwork
    wc => wave_config

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,jd,i_startidx,i_endidx,temp,temp_max,peak_ind) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs

        DO jc = i_startidx, i_endidx
          temp(jc,jf)   = 0._wp
        END DO

        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx
            temp(jc,jf) = temp(jc,jf) + wesd(jf)%ptr(jc,jd,jb)
          END DO
        END DO  ! jd

      END DO  ! jf

      DO jc = i_startidx, i_endidx
        temp_max(jc) = 0._wp
        peak_ind(jc) = 1
      ENDDO

      ! get frequency index with maximum energy
      !
      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          IF (temp(jc,jf) > temp_max(jc)) THEN
            temp_max(jc) = temp(jc,jf)
            peak_ind(jc) = jf
          ENDIF
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
        pp(jc,jb) = MERGE(1.0_wp, 1.0_wp /wc%freqs(peak_ind(jc)), peak_ind(jc) == 1)
        kp(jc,jb) = MERGE(1.0_wp, wave_num_c(jc,peak_ind(jc),jb), peak_ind(jc) == 1)
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE peak_wave_period_wavenumber


  !>
  !! Calculation of peak period and wavenumber for wind sea and swell according to mask
  !! based on WAM 4.5 formulation
  !!
  SUBROUTINE peak_wave_period_wavenumber_sep(p_patch, wave_config, wesd, wave_num_c, mask, pp_sea, pp_swell, kp_sea, kp_swell)
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':peak_wave_period_wavenumber_sep'

    TYPE(t_patch),       INTENT(IN)         :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN) :: wave_config
    TYPE(t_wesd),        INTENT(IN)         :: wesd(:)           !< energy spectral bins
    REAL(wp),            INTENT(IN)         :: wave_num_c(:,:,:) !< wavenumber (1/m)
    INTEGER,             INTENT(IN)         :: mask(:,:,:,:)     !< =1 where energy belongs to swell, =0 - belongs to wind sea
    REAL(wp),            INTENT(INOUT)      :: pp_sea(:,:)       !< wind sea peak wave period
    REAL(wp),            INTENT(INOUT)      :: pp_swell(:,:)     !< swell peak wave period
    REAL(wp),            INTENT(INOUT)      :: kp_sea(:,:)       !< wind sea peak wavenumber
    REAL(wp),            INTENT(INOUT)      :: kp_swell(:,:)     !< swell peak wavenumber

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx

    INTEGER :: jc,jb,jf,jd
    INTEGER :: peak_ind1(nproma)
    INTEGER :: peak_ind2(nproma)

    REAL(wp):: temp1(nproma,wave_config%nfreqs)
    REAL(wp):: temp2(nproma,wave_config%nfreqs)
    REAL(wp):: temp1_max(nproma)        !< maximum of temp1 over all frequencies
    REAL(wp):: temp2_max(nproma)        !< maximum of temp2 over all frequencies

    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! save some paperwork
    wc => wave_config
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,jd,i_startidx,i_endidx,temp1,temp2,temp1_max,temp2_max, &
!$OMP            peak_ind1,peak_ind2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          temp1(jc,jf)   = 0._wp
          temp2(jc,jf)   = 0._wp
        END DO
        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx
            IF (mask(jc,jd,jf,jb) == 1) THEN !swell
              temp1(jc,jf) = temp1(jc,jf) + wesd(jf)%ptr(jc,jd,jb)
            ELSE
              temp2(jc,jf) = temp2(jc,jf) + wesd(jf)%ptr(jc,jd,jb)
            END IF
          END DO
        END DO  ! jd
      END DO  ! jf
      DO jc = i_startidx, i_endidx
        temp1_max(jc) = 0._wp
        temp2_max(jc) = 0._wp
        peak_ind1(jc) = 1
        peak_ind2(jc) = 1
      ENDDO
      ! get frequency index with maximum energy
      !
      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          IF (temp1(jc,jf) > temp1_max(jc)) THEN
            temp1_max(jc) = temp1(jc,jf)
            peak_ind1(jc) = jf
          ENDIF
          IF (temp2(jc,jf) > temp2_max(jc)) THEN
            temp2_max(jc) = temp2(jc,jf)
            peak_ind2(jc) = jf
          ENDIF
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
        pp_swell(jc,jb) = MERGE(1.0_wp, 1.0_wp/wc%freqs(peak_ind1(jc)), peak_ind1(jc) == 1)
        pp_sea(jc,jb)   = MERGE(1.0_wp, 1.0_wp/wc%freqs(peak_ind2(jc)), peak_ind2(jc) == 1)
        kp_swell(jc,jb) = MERGE(1.0_wp, wave_num_c(jc,peak_ind1(jc),jb), peak_ind1(jc) == 1)
        kp_sea(jc,jb)   = MERGE(1.0_wp, wave_num_c(jc,peak_ind2(jc),jb), peak_ind2(jc) == 1)
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE peak_wave_period_wavenumber_sep

  !>
  !! Calculation of wave drag coefficient and normalised stress
  !! based on WAM 4.5 formulation
  !!
  SUBROUTINE wave_drag_stress_ch_par(p_patch, sp10m, ustar, tauw, z0, drag, tauwn, beta)

    TYPE(t_patch),     INTENT(IN)    :: p_patch
    REAL(wp),          INTENT(IN)    :: sp10m(:,:)
    REAL(wp),          INTENT(IN)    :: ustar(:,:)
    REAL(wp),          INTENT(IN)    :: tauw(:,:)
    REAL(wp),          INTENT(IN)    :: z0(:,:)
    REAL(wp),          INTENT(INOUT) :: drag(:,:)  ! drag coefficient
    REAL(wp),          INTENT(INOUT) :: tauwn(:,:) ! normalised wave stress
    REAL(wp),          INTENT(INOUT) :: beta(:,:)  ! Charnock parameter

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':wave_drag_stress_ch_par'

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb


    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx

        drag(jc,jb) = (ustar(jc,jb)*ustar(jc,jb) + 0.0001_wp) &
          &         / (sp10m(jc,jb)*sp10m(jc,jb) + 0.01_wp)

        tauwn(jc,jb) = tauw(jc,jb) / (ustar(jc,jb)*ustar(jc,jb) + 0.0001_wp)

        beta(jc,jb) = grav * z0(jc,jb) / MAX(ustar(jc,jb)*ustar(jc,jb), 1.0E-6_wp)

      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE wave_drag_stress_ch_par

  !>
  !! Calculation of mean wave direction and spread
  !! based on WAM 4.5 formulation
  !!
  SUBROUTINE mean_wave_direction_spread(p_patch, wave_config, wesd, mean_dir, mean_spread)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':mean_wave_direction_spread'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wesd),                INTENT(IN)    :: wesd(:)    !energy spectral bins
    REAL(wp),                    INTENT(INOUT) :: mean_dir(:,:)
    REAL(wp),                    INTENT(INOUT) :: mean_spread(:,:)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    REAL(wp) :: si(nproma), ci(nproma), temp_dsum(nproma)
    REAL(wp) :: temp(nproma,wave_config%ndirs)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd

    wc => wave_config

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,jd,i_startidx,i_endidx,si,ci,temp,temp_dsum) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      ! initialisation of si, ci
      DO jc = i_startidx, i_endidx
        si(jc) = 0._wp
        ci(jc) = 0._wp
        temp_dsum(jc) = 0._wp
      END DO

      ! initialisation of temp
      DO jd = 1, wc%ndirs
        DO jc = i_startidx, i_endidx
          temp(jc,jd) = 0._wp
        END DO
      END DO

      DO jf = 1,wc%nfreqs
        DO jd = 1, wc%ndirs
          DO jc = i_startidx, i_endidx
            temp(jc,jd) = temp(jc,jd) + wesd(jf)%ptr(jc,jd,jb) * wc%DFIM(jf)
          END DO
        END DO
      END DO

      DO jd = 1, wc%ndirs
        DO jc = i_startidx, i_endidx
          si(jc) = si(jc) + temp(jc,jd) * wc%sin_dir(jd)
          ci(jc) = ci(jc) + temp(jc,jd) * wc%cos_dir(jd)
          temp_dsum(jc) = temp_dsum(jc) + temp(jc,jd)
        END DO
      END DO

      ! compute mean direction and spread
      DO jc = i_startidx, i_endidx
        ! mean direction
        IF (ci(jc) .EQ. 0.0_wp) ci(jc) = 0.1E-30_wp
        mean_dir(jc,jb) = ATAN2(si(jc),ci(jc))
        IF (mean_dir(jc,jb) .LT. 0.0_wp) mean_dir(jc,jb) = mean_dir(jc,jb) + pi2
        IF (mean_dir(jc,jb) .GT. pi2-0.0001_wp) mean_dir(jc,jb) = 0.0_wp

        mean_dir(jc,jb) = mean_dir(jc,jb) * rad2deg

        ! mean spread
        mean_spread(jc,jb) = temp_dsum(jc) !SUM(temp(jc,:))
        IF (ABS(ci(jc)) .LT. 0.1E-15_wp) ci(jc) = SIGN(0.1E-15_wp, ci(jc))
        IF (ABS(si(jc)) .LT. 0.1E-15_wp) si(jc) = SIGN(0.1E-15_wp, si(jc))
        IF (ABS(mean_spread(jc,jb)) .LT. 0.1E-15_wp) mean_spread(jc,jb) = &
            SIGN(0.1E-15_wp, mean_spread(jc,jb))

        mean_spread(jc,jb) = 2.0_wp * &
            (1.0_wp - SQRT(si(jc)*si(jc) + ci(jc)*ci(jc)) / mean_spread(jc,jb))

        IF (mean_spread(jc,jb) .LE. 0.0_wp) THEN
          mean_spread(jc,jb) = TINY(1.0_wp)
        ELSE
          mean_spread(jc,jb) = SQRT(mean_spread(jc,jb))
        END IF

        mean_spread(jc,jb) = mean_spread(jc,jb) * rad2deg
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE mean_wave_direction_spread


  !>
  !! Calculation of mean wind sea and swell wave direction and spread
  !! according to mask
  !! based on WAM 4.5 formulation
  !!
  SUBROUTINE mean_wave_direction_spread_sep(p_patch, wave_config, wesd, mask, md_sea, ms_sea, md_swell, ms_swell)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':mean_wave_direction_spread'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wesd),                INTENT(IN)    :: wesd(:)         !energy spectral bins
    INTEGER,                     INTENT(IN)    :: mask(:,:,:,:)   !=1 where energy belongs to swell,
                                                                  !=0 - belongs to wind sea
    REAL(wp),                    INTENT(INOUT) :: md_sea(:,:) ! wind sea direction
    REAL(wp),                    INTENT(INOUT) :: ms_sea(:,:) ! wind sea spread
    REAL(wp),                    INTENT(INOUT) :: md_swell(:,:) ! swell direction
    REAL(wp),                    INTENT(INOUT) :: ms_swell(:,:) ! swell spread

    TYPE(t_wave_config), POINTER :: wc => NULL()

    REAL(wp) :: si1(nproma), ci1(nproma), si2(nproma), ci2(nproma), temp1_dsum(nproma), temp2_dsum(nproma)
    REAL(wp) :: temp1(nproma,wave_config%ndirs), temp2(nproma,wave_config%ndirs)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd

    wc => wave_config

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,jd,i_startidx,i_endidx,si1,ci1,temp1,si2,ci2,temp2,temp1_dsum,temp2_dsum) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      ! initialisation of si, ci
      DO jc = i_startidx, i_endidx
        si1(jc) = 0._wp
        ci1(jc) = 0._wp
        si2(jc) = 0._wp
        ci2(jc) = 0._wp
        temp1_dsum(jc) = 0._wp
        temp2_dsum(jc) = 0._wp
      END DO

      ! initialisation of temps
      DO jd = 1, wc%ndirs
        DO jc = i_startidx, i_endidx
          temp1(jc,jd) = 0._wp
          temp2(jc,jd) = 0._wp
        END DO
      END DO

      DO jf = 1,wc%nfreqs
        DO jd = 1, wc%ndirs
          DO jc = i_startidx, i_endidx
            IF (mask(jc,jd,jf,jb).EQ.1) THEN
              temp1(jc,jd) = temp1(jc,jd) + wesd(jf)%ptr(jc,jd,jb) * wc%DFIM(jf) ! swell
            ELSE
              temp2(jc,jd) = temp2(jc,jd) + wesd(jf)%ptr(jc,jd,jb) * wc%DFIM(jf) ! wind sea
            END IF
          END DO
        END DO
      END DO

      DO jd = 1, wc%ndirs
        DO jc = i_startidx, i_endidx
          si1(jc) = si1(jc) + temp1(jc,jd) * wc%sin_dir(jd)
          ci1(jc) = ci1(jc) + temp1(jc,jd) * wc%cos_dir(jd)
          si2(jc) = si2(jc) + temp2(jc,jd) * wc%sin_dir(jd)
          ci2(jc) = ci2(jc) + temp2(jc,jd) * wc%cos_dir(jd)
          temp1_dsum(jc) = temp1_dsum(jc) + temp1(jc,jd)
          temp2_dsum(jc) = temp2_dsum(jc) + temp2(jc,jd)
        END DO
      END DO

      ! compute mean direction and spread
      DO jc = i_startidx, i_endidx
        ! swell part
        ! mean direction
        IF (ci1(jc) .EQ. 0.0_wp) ci1(jc) = 0.1E-30_wp
        md_swell(jc,jb) = ATAN2(si1(jc),ci1(jc))
        IF (md_swell(jc,jb) .LT. 0.0_wp) md_swell(jc,jb) = md_swell(jc,jb) + pi2
        IF (md_swell(jc,jb) .GT. pi2-0.0001_wp) md_swell(jc,jb) = 0.0_wp

        md_swell(jc,jb) = md_swell(jc,jb) * rad2deg

        ! mean spread
        ms_swell(jc,jb) = temp1_dsum(jc) !SUM(temp1(jc,:))
        IF (ABS(ci1(jc)) .LT. 0.1E-15_wp) ci1(jc) = SIGN(0.1E-15_wp, ci1(jc))
        IF (ABS(si1(jc)) .LT. 0.1E-15_wp) si1(jc) = SIGN(0.1E-15_wp, si1(jc))
        IF (ABS(ms_swell(jc,jb)) .LT. 0.1E-15_wp) ms_swell(jc,jb) = &
            SIGN(0.1E-15_wp, ms_swell(jc,jb))

        ms_swell(jc,jb) = 2.0_wp * &
            (1.0_wp - SQRT(si1(jc)*si1(jc) + ci1(jc)*ci1(jc)) / ms_swell(jc,jb))

        IF (ms_swell(jc,jb) .LE. 0.0_wp) THEN
          ms_swell(jc,jb) = TINY(1.0_wp)
        ELSE
          ms_swell(jc,jb) = SQRT(ms_swell(jc,jb))
        END IF

        ms_swell(jc,jb) = ms_swell(jc,jb) * rad2deg

        ! wind sea part
        ! mean direction
        IF (ci2(jc) .EQ. 0.0_wp) ci2(jc) = 0.1E-30_wp
        md_sea(jc,jb) = ATAN2(si2(jc),ci2(jc))
        IF (md_sea(jc,jb) .LT. 0.0_wp) md_sea(jc,jb) = md_sea(jc,jb) + pi2
        IF (md_sea(jc,jb) .GT. pi2-0.0001_wp) md_sea(jc,jb) = 0.0_wp

        md_sea(jc,jb) = md_sea(jc,jb) * rad2deg

        ! mean spread
        ms_sea(jc,jb) = temp2_dsum(jc) !SUM(temp2(jc,:))
        IF (ABS(ci2(jc)) .LT. 0.1E-15_wp) ci2(jc) = SIGN(0.1E-15_wp, ci2(jc))
        IF (ABS(si2(jc)) .LT. 0.1E-15_wp) si2(jc) = SIGN(0.1E-15_wp, si2(jc))
        IF (ABS(ms_sea(jc,jb)) .LT. 0.1E-15_wp) ms_sea(jc,jb) = &
            SIGN(0.1E-15_wp, ms_sea(jc,jb))

        ms_sea(jc,jb) = 2.0_wp * &
            (1.0_wp - SQRT(si2(jc)*si2(jc) + ci2(jc)*ci2(jc)) / ms_sea(jc,jb))

        IF (ms_sea(jc,jb) .LE. 0.0_wp) THEN
          ms_sea(jc,jb) = TINY(1.0_wp)
        ELSE
          ms_sea(jc,jb) = SQRT(ms_sea(jc,jb))
        END IF

        ms_sea(jc,jb) = ms_sea(jc,jb) * rad2deg

      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE mean_wave_direction_spread_sep

END MODULE mo_wave_output_diagnostics
