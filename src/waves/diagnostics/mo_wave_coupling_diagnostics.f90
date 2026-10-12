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

! Computes diagnostic variables for coupling and output, if requested.
! Fields have no impact on the wave solution, but may impact coupled components.
! In case of coupling, diagnostics are updated every time step.
! In standalone runs, diagnostics are updated at output time steps.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_coupling_diagnostics
  USE mo_kind,                ONLY: wp, vp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_wave_config,         ONLY: t_wave_config
  USE mo_wave_types,          ONLY: t_wave_diag_dyn, t_wave_diag_cpl, t_wesd
  USE mo_impl_constants,      ONLY: min_rlcell, MAX_CHAR_LENGTH
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_physical_constants,  ONLY: grav
  USE mo_math_constants,      ONLY: pi2, rad2deg
  USE mo_parallel_config,     ONLY: nproma
  USE mo_fortran_tools,       ONLY: init
  USE mo_wave_constants,      ONLY: EMIN
  USE mo_wave_stokes,         ONLY: stokes_profile_spectrum, stokes_profile_breivik, &
    &                               stokes_drift
  USE mo_wave_common_diagnostics, ONLY: significant_wave_height

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_coupling_diagnostics'

  PUBLIC :: calculate_coupling_diagnostics
  PUBLIC :: wave_stress_ocean

CONTAINS
  !>
  !! Calculation of purely diagnostic parameters
  !!
  SUBROUTINE calculate_coupling_diagnostics(p_patch, wave_config, depth, wesd, diag_dyn, diag_cpl)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_config),   INTENT(IN)    :: wave_config
    REAL(wp),              INTENT(IN)    :: depth(:,:)    ! water depth
    TYPE(t_wesd),          INTENT(IN)    :: wesd(:)       ! energy spectral bins
    TYPE(t_wave_diag_dyn), INTENT(IN)    :: diag_dyn
    TYPE(t_wave_diag_cpl), INTENT(INOUT) :: diag_cpl

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':calculate_coupling_diagnostics'


    IF (ASSOCIATED(diag_cpl%hs)) THEN
      ! calculate significant wave height from total wave energy
      !
      CALL significant_wave_height(p_patch = p_patch, &
        &                          emean   = diag_dyn%emean(:,:), &
        &                          hs      = diag_cpl%hs(:,:)) ! OUT
    ENDIF


    ! calculate stokes drift velocities
    !
    ! surface values
    ! Note that the Breivik-type 3D stokes profile below requires
    ! Stokes surface values. Hence, we compute u_stokes and v_stokes
    ! unconditionally, for simplicity.
    !
    CALL stokes_drift(p_patch = p_patch,             &
      &           wave_config = wave_config,         &
      &            wave_num_c = diag_dyn%wave_num_c, &
      &                 depth = depth,               &
      &                  wesd = wesd,                &
      &              u_stokes = diag_cpl%u_stokes,   & ! OUT
      &              v_stokes = diag_cpl%v_stokes)     ! OUT

    ! vertical profile
    IF (ASSOCIATED(diag_cpl%last_idx_depth) .AND. &
      & ASSOCIATED(diag_cpl%u3d_stokes)     .AND. &
      & ASSOCIATED(diag_cpl%v3d_stokes)) THEN

      IF (wave_config%stokes_method == 1) THEN

        CALL stokes_profile_spectrum(p_patch = p_patch,      &
          &           wave_config = wave_config,             &
          &            wave_num_c = diag_dyn%wave_num_c,     &
          &                 depth = depth,                   &
          &        last_idx_depth = diag_cpl%last_idx_depth, &
          &                  wesd = wesd,                    &
          &            u3d_stokes = diag_cpl%u3d_stokes,     & ! OUT
          &            v3d_stokes = diag_cpl%v3d_stokes)       ! OUT

      ELSE

        CALL stokes_profile_breivik(p_patch = p_patch,       & ! IN
          &           wave_config = wave_config,             & ! IN
          &            wave_num_c = diag_dyn%wave_num_c,     & ! IN
          &                 depth = depth,                   & ! IN
          &        last_idx_depth = diag_cpl%last_idx_depth, & ! IN
          &                  wesd = wesd,                    & ! IN
          &              u_stokes = diag_cpl%u_stokes,       & ! IN
          &              v_stokes = diag_cpl%v_stokes,       & ! IN
          &            u3d_stokes = diag_cpl%u3d_stokes,     & ! OUT
          &            v3d_stokes = diag_cpl%v3d_stokes)       ! OUT
      END IF
    END IF

  END SUBROUTINE calculate_coupling_diagnostics


  !>
  !! Calculation of wave-to-ocean stress and energy flux
  !!
  !! Adaptation of WAM 4.5 code STRESSO.f90
  !! When coupled to waves, the ocean experiences the 'wave-to-ocean' stress
  !! \tau_oc = \tau_a - (\tau_w + \tau_ds)
  !! i.e. atmospheric stress minus stress which contributes to wave growth, plus stress
  !! imparted to the ocean by wave breaking. The wave-to-ocean stress is the difference
  !! between \tau_a and the integral of sl/c, with sl the sum of all the source terms.
  !!
  SUBROUTINE wave_stress_ocean(p_patch, wave_config, taua_x, taua_y, sl, diag_dyn, diag_cpl)
     CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
          &  routine = modname//'wave_stress_ocean'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: taua_x(:,:)
    REAL(wp),                    INTENT(IN)    :: taua_y(:,:)
    REAL(vp),                    INTENT(IN)    :: sl(:,:,:,:)
    TYPE(t_wave_diag_dyn),       INTENT(IN)    :: diag_dyn
    TYPE(t_wave_diag_cpl),       INTENT(INOUT) :: diag_cpl

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd

    REAL(wp) :: stotplus, cmrhowgdfth_rhoa
    REAL(wp) :: rhowgdfth(nproma,wave_config%nfreqs)
    REAL(wp) :: cm(nproma,wave_config%nfreqs)
    REAL(wp) :: sumt(nproma), sumx(nproma), sumy(nproma)
    REAL(wp) :: xstress(nproma), ystress(nproma)
    REAL(wp) :: roair    ! air density

    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! save some paperwork
    wc => wave_config

    roair = MAX(wc%roair,1._wp)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,jd,i_startidx,i_endidx,cm,sumt,sumx,sumy,   &
!$OMP            stotplus,rhowgdfth,cmrhowgdfth_rhoa,xstress,ystress) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,          &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          cm(jc,jf) = diag_dyn%wave_num_c(jc,jf,jb) * 1.0_wp/(pi2*wc%freqs(jf))
          rhowgdfth(jc,jf) = MERGE(wc%rhowg_dfim(jf), 0.0_wp, jf <= diag_dyn%last_prog_freq_ind(jc,jb))
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
        jf = diag_dyn%last_prog_freq_ind(jc,jb)
        IF (jf /= wc%nfreqs) rhowgdfth(jc,jf) = 0.5_wp * rhowgdfth(jc,jf)

        !initialisation
        xstress(jc) = taua_x(jc,jb)
        ystress(jc) = taua_y(jc,jb)
        diag_cpl%phioc(jc,jb) = diag_dyn%phiaw(jc,jb)
      END DO

      !sum
      DO jf = 1, MAXVAL(diag_dyn%last_prog_freq_ind(i_startidx:i_endidx,jb))
        DO jc = i_startidx, i_endidx
          sumt(jc) = 0._wp
          sumx(jc) = 0._wp
          sumy(jc) = 0._wp
        END DO

        DO jd = 1, wc%ndirs
          DO jc = i_startidx, i_endidx
            stotplus = MAX(sl(jc,jd,jf,jb),0._wp)
            sumt(jc) = sumt(jc) + stotplus
            sumx(jc) = sumx(jc) + stotplus * wc%sin_dir(jd)
            sumy(jc) = sumy(jc) + stotplus * wc%cos_dir(jd)
          END DO
        END DO

        DO jc = i_startidx, i_endidx
          diag_cpl%phioc(jc,jb) =  diag_cpl%phioc(jc,jb) - sumt(jc)*rhowgdfth(jc,jf)
          cmrhowgdfth_rhoa = cm(jc,jf) * rhowgdfth(jc,jf)/roair
          xstress(jc) = xstress(jc) - sumx(jc)*cmrhowgdfth_rhoa
          ystress(jc) = ystress(jc) - sumy(jc)*cmrhowgdfth_rhoa
        END DO
      END DO  ! jf

      DO jc = i_startidx, i_endidx
        diag_cpl%tauoc_x(jc,jb) = xstress(jc)
        diag_cpl%tauoc_y(jc,jb) = ystress(jc)
        diag_cpl%tauoc(jc,jb)   = SQRT(xstress(jc)**2 + ystress(jc)**2)
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE wave_stress_ocean

END MODULE mo_wave_coupling_diagnostics
