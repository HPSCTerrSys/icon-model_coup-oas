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
! Computes the Stokes depth profile using either full spectral information
! or a parametrisation

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_stokes
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_wave_config,         ONLY: t_wave_config
  USE mo_wave_types,          ONLY: t_wesd
  USE mo_impl_constants,      ONLY: min_rlcell
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_physical_constants,  ONLY: grav
  USE mo_math_constants,      ONLY: pi2, rad2deg
  USE mo_parallel_config,     ONLY: nproma
  USE mo_kind,                ONLY: wp
  USE mo_fortran_tools,       ONLY: init
  USE mo_wave_constants,      ONLY: EMIN, EPS1


  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_stokes'

  PUBLIC :: stokes_drift
  PUBLIC :: stokes_profile_spectrum
  PUBLIC :: stokes_profile_breivik

CONTAINS

  !>
  !! Calculation of Stokes drift components
  !!
  !! Adaptation of WAM 4.5 code of the subroutine STOKES_DRIFT
  !! developed by M.REISTAD, O.SAETRA, and H.GUNTHER
  !!
  !! References:
  !! Kern E. Kenyon, JGR, Vol 74 NO 28, 1969
  !! O. Breivik, J.-R. Bidlot & P. Janssen, 2016 (high-frequency tail)
  !!
  SUBROUTINE stokes_drift(p_patch, wave_config, wave_num_c, depth, wesd, u_stokes, v_stokes)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':stokes_drift'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: wave_num_c(:,:,:)  !< wave number (1/m)
    REAL(wp),                    INTENT(IN)    :: depth(:,:)
    TYPE(t_wesd),                INTENT(IN)    :: wesd(:)            !< energy spectral bins
    REAL(wp),                    INTENT(INOUT) :: u_stokes(:,:)
    REAL(wp),                    INTENT(INOUT) :: v_stokes(:,:)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    REAL(wp) :: ak, akd, fact
    REAL(wp) :: si(nproma), ci(nproma)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd

    wc => wave_config

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
    CALL init(u_stokes, lacc=.FALSE.)
    CALL init(v_stokes, lacc=.FALSE.)
!$OMP BARRIER
!$OMP DO PRIVATE(jb,jc,jf,jd,i_startidx,i_endidx,ak,akd,si,ci,fact) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      ! initialisation of si, ci
      DO jc = i_startidx, i_endidx
        si(jc) = 0._wp
        ci(jc) = 0._wp
      END DO

      freqs:DO jf = 1,wc%nfreqs
        DO jd = 1, wc%ndirs
          DO jc = i_startidx, i_endidx
            si(jc) = si(jc) + wesd(jf)%ptr(jc,jd,jb) * wc%sin_dir(jd)
            ci(jc) = ci(jc) + wesd(jf)%ptr(jc,jd,jb) * wc%cos_dir(jd)
          END DO
        END DO

        DO jc = i_startidx, i_endidx
          ak = wave_num_c(jc,jf,jb)
          akd = ak * depth(jc,jb)
          fact = 2._wp*grav*ak**2/(pi2*wc%freqs(jf)*TANH(2._wp*akd)) * wc%DFIM(jf)
          si(jc) = fact * si(jc)
          ci(jc) = fact * ci(jc)
          u_stokes(jc,jb) = u_stokes(jc,jb) + si(jc)
          v_stokes(jc,jb) = v_stokes(jc,jb) + ci(jc)
        END DO

      END DO freqs

      ! Addition of HF tail following Breivik (2016)
      DO jc = i_startidx, i_endidx
        si(jc)  = 0._wp
        ci(jc)  = 0._wp
      ENDDO

      DO jd = 1, wc%ndirs
        DO jc = i_startidx, i_endidx
          si(jc) = si(jc) + 2._wp*wesd(wc%nfreqs)%ptr(jc,jd,jb) * wc%sin_dir(jd) *  &
                          &  wave_num_c(jc,wc%nfreqs,jb) * pi2*wc%freqs(wc%nfreqs)**2
          ci(jc) = ci(jc) + 2._wp*wesd(wc%nfreqs)%ptr(jc,jd,jb) * wc%cos_dir(jd) *  &
                          &  wave_num_c(jc,wc%nfreqs,jb) * pi2*wc%freqs(wc%nfreqs)**2
        END DO
      END DO

      DO jc = i_startidx, i_endidx
        u_stokes(jc,jb) = u_stokes(jc,jb) + si(jc)
        v_stokes(jc,jb) = v_stokes(jc,jb) + ci(jc)
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE stokes_drift


  !>
  !! Calculation of Stokes drift profile from the full spectral integral
  !! Reference:
  !! Kern E. Kenyon, JGR, Vol 74 NO 28, 1969
  !! O. Breivik, J.-R. Bidlot & P. Janssen (2016) (high-frequency tail)
  !!
  SUBROUTINE stokes_profile_spectrum(p_patch, wave_config, wave_num_c, depth, last_idx_depth, wesd, u3d_stokes, v3d_stokes)

    CHARACTER(*), PARAMETER :: routine = modname//'::stokes_profile'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: wave_num_c(:,:,:)  !< wave number (1/m)
    REAL(wp),                    INTENT(IN)    :: depth(:,:)
    INTEGER,                     INTENT(IN)    :: last_idx_depth(:,:)
    TYPE(t_wesd),                INTENT(IN)    :: wesd(:)            !< energy spectral bins
    REAL(wp),                    INTENT(INOUT) :: u3d_stokes(:,:,:)
    REAL(wp),                    INTENT(INOUT) :: v3d_stokes(:,:,:)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    REAL(wp) :: ak, akd, fact, akcz
    REAL(wp) :: si(nproma), ci(nproma), kc(nproma)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd,jk

    wc => wave_config

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
    CALL init(u3d_stokes, lacc=.FALSE.)
    CALL init(v3d_stokes, lacc=.FALSE.)
!$OMP BARRIER
!$OMP DO PRIVATE(jb,jc,jf,jd,jk,i_startidx,i_endidx,ak,akd,si,ci,fact,kc,akcz) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1,wc%oce_stokes_nlev
        ! initialisation of si, ci
        DO jc = i_startidx, i_endidx
          si(jc) = 0._wp
          ci(jc) = 0._wp
        END DO

        freqs:DO jf = 1,wc%nfreqs
          DO jd = 1, wc%ndirs
            DO jc = i_startidx, i_endidx
              si(jc) = si(jc) + wesd(jf)%ptr(jc,jd,jb) * wc%sin_dir(jd)
              ci(jc) = ci(jc) + wesd(jf)%ptr(jc,jd,jb) * wc%cos_dir(jd)
            END DO
          END DO

          DO jc = i_startidx, i_endidx
            ak = wave_num_c(jc,jf,jb)
            akd = ak * MIN(depth(jc,jb),wc%stokes_depth)

            ! Stokes drift integrand factor as function of depth
            fact = 2._wp*grav*ak**2  * COSH(2._wp*akd - 2._wp*ak*MIN(wc%oce_stokes_mc(jk),wc%oce_stokes_mc(last_idx_depth(jc,jb)))) &
              &    / ( pi2*wc%freqs(jf) * SINH(2._wp*akd) ) * wc%DFIM(jf)

            si(jc) = fact * si(jc)
            ci(jc) = fact * ci(jc)

            u3d_stokes(jc,jk,jb) = u3d_stokes(jc,jk,jb) + si(jc)
            v3d_stokes(jc,jk,jb) = v3d_stokes(jc,jk,jb) + ci(jc)

          END DO
        END DO freqs

      END DO


    !------------------------------------------------------------------------------------------------
    !! Calculate Stokes drift high-frequency tail as function of z: Philips spectrum assumed for f>fc
    ! Reference:
    ! O. Breivik, J.-R. Bidlot & P. Janssen, 2016
    !------------------------------------------------------------------------------------------------

      DO jc = i_startidx, i_endidx
        si(jc) = 0._wp
        ci(jc) = 0._wp
        kc(jc) = wave_num_c(jc,wc%nfreqs,jb)
      END DO

      DO jd = 1, wc%ndirs
       DO jc = i_startidx, i_endidx
          si(jc) = si(jc) + wesd(wc%nfreqs)%ptr(jc,jd,jb) * wc%sin_dir(jd) * pi2*wc%freqs(wc%nfreqs) * kc(jc)
          ci(jc) = ci(jc) + wesd(wc%nfreqs)%ptr(jc,jd,jb) * wc%cos_dir(jd) * pi2*wc%freqs(wc%nfreqs) * kc(jc)
       END DO
      END DO

      DO jk = 1,MAXVAL(last_idx_depth(i_startidx:i_endidx,jb))
        DO jc = i_startidx, i_endidx
          IF (jk <= last_idx_depth(jc,jb)) THEN
            akcz = -kc(jc)*wc%oce_stokes_mc(jk)
            u3d_stokes(jc,jk,jb) = u3d_stokes(jc,jk,jb) + si(jc)*( EXP(2._wp*akcz) - SQRT(-pi2*akcz)*ERFC(SQRT(-2._wp*akcz)) )
            v3d_stokes(jc,jk,jb) = v3d_stokes(jc,jk,jb) + ci(jc)*( EXP(2._wp*akcz) - SQRT(-pi2*akcz)*ERFC(SQRT(-2._wp*akcz)) )
          END IF
        END DO
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE stokes_profile_spectrum



  !>
  !! Parametrisation of Stokes profile from Philips spectrum using surface value and Stokes transport
  !! Reference:
  !! O. Breivik, J.-R. Bidlot & P. Janssen (2016)
  !!
  SUBROUTINE stokes_profile_breivik(p_patch, wave_config, wave_num_c, depth, last_idx_depth, wesd, &
                                  & u_stokes, v_stokes, u3d_stokes, v3d_stokes)

    CHARACTER(*), PARAMETER :: routine = modname//'::stokes_profile'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: wave_num_c(:,:,:)  !< wave number (1/m)
    REAL(wp),                    INTENT(IN)    :: depth(:,:)
    INTEGER,                     INTENT(IN)    :: last_idx_depth(:,:)
    TYPE(t_wesd),                INTENT(IN)    :: wesd(:)            !< energy spectral bins
    REAL(wp),                    INTENT(IN)    :: u_stokes(:,:)
    REAL(wp),                    INTENT(IN)    :: v_stokes(:,:)
    REAL(wp),                    INTENT(INOUT) :: u3d_stokes(:,:,:)
    REAL(wp),                    INTENT(INOUT) :: v3d_stokes(:,:,:)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    REAL(wp) :: ak, akbz, akcz
    REAL(wp) :: temp(nproma), si(nproma), ci(nproma), kc(nproma), ust(nproma), vst(nproma)
    REAL(wp) :: kbar(nproma)          ! Breivik wavenumber (1/m)
    REAL(wp) :: Tstokes(nproma)       ! Magnitude of Stokes transport (m2s-1)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd,jk

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    wc => wave_config


!$OMP PARALLEL
    CALL init(u3d_stokes, lacc=.FALSE.)
    CALL init(v3d_stokes, lacc=.FALSE.)
!$OMP BARRIER
!$OMP DO PRIVATE(jb,jc,jf,jd,jk,i_startidx,i_endidx,si,ci,kc,ust,vst,temp,ak,Tstokes,kbar,akbz,akcz) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      !---------------------------------------------------------------------------------------------
      ! Subtraction of HF tail from surface Stokes drift to calculate vertical profile

      DO jc = i_startidx, i_endidx
        si(jc)  = 0._wp
        ci(jc)  = 0._wp
        kc(jc) = wave_num_c(jc,wc%nfreqs,jb)
        kbar(jc) = 0._wp
        Tstokes(jc) = 0._wp
      END DO

      DO jd = 1, wc%ndirs
        DO jc = i_startidx, i_endidx
          si(jc) = si(jc) + 2._wp*wesd(wc%nfreqs)%ptr(jc,jd,jb) * wc%sin_dir(jd)*kc(jc)*pi2*wc%freqs(wc%nfreqs)**2
          ci(jc) = ci(jc) + 2._wp*wesd(wc%nfreqs)%ptr(jc,jd,jb) * wc%cos_dir(jd)*kc(jc)*pi2*wc%freqs(wc%nfreqs)**2
        END DO
      END DO

      DO jc = i_startidx, i_endidx
        ust(jc) = u_stokes(jc,jb) - si(jc)
        vst(jc) = v_stokes(jc,jb) - ci(jc)
      END DO

      !---------------------------------------------------------------------------------------------


      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          temp(jc) = 0._wp
        END DO

        DO jd = 1, wc%ndirs
          DO jc = i_startidx, i_endidx
            temp(jc) = temp(jc) + wesd(jf)%ptr(jc,jd,jb)
          END DO
        END DO

        DO jc = i_startidx, i_endidx
          ak = wave_num_c(jc,jf,jb)
          Tstokes(jc) = Tstokes(jc)+temp(jc)*grav*ak/(pi2*wc%freqs(jf))*wc%DFIM(jf)
        END DO
      END DO !jf

      DO jc = i_startidx, i_endidx
        kbar(jc) = SQRT(ust(jc)**2+vst(jc)**2)/MAX(6._wp*Tstokes(jc),EPS1)
      END DO

      DO jk = 1,MAXVAL(last_idx_depth(i_startidx:i_endidx,jb))
        DO jc = i_startidx, i_endidx
          IF (jk <= last_idx_depth(jc,jb)) THEN
            akbz = -kbar(jc)*wc%oce_stokes_mc(jk)
            u3d_stokes(jc,jk,jb) = ust(jc)*( EXP(2._wp*akbz) - SQRT(-pi2*akbz)*ERFC(SQRT(-2._wp*akbz)) )
            v3d_stokes(jc,jk,jb) = vst(jc)*( EXP(2._wp*akbz) - SQRT(-pi2*akbz)*ERFC(SQRT(-2._wp*akbz)) )
          END IF
        END DO
      END DO

      !------------------------------------------------------------------------------------------------
      ! Calculate Stokes drift high-frequency tail as function of z: Philips spectrum assumed for f>fc
      ! Reference:
      ! O. Breivik, J.-R. Bidlot & P. Janssen, 2016
      !------------------------------------------------------------------------------------------------

      DO jk = 1,MAXVAL(last_idx_depth(i_startidx:i_endidx,jb))
        DO jc = i_startidx, i_endidx
          IF (jk <= last_idx_depth(jc,jb)) THEN
            akcz = -kc(jc)*wc%oce_stokes_mc(jk)
            u3d_stokes(jc,jk,jb) = u3d_stokes(jc,jk,jb) + si(jc)*( EXP(2._wp*akcz) - SQRT(-pi2*akcz)*ERFC(SQRT(-2._wp*akcz)) )
            v3d_stokes(jc,jk,jb) = v3d_stokes(jc,jk,jb) + ci(jc)*( EXP(2._wp*akcz) - SQRT(-pi2*akcz)*ERFC(SQRT(-2._wp*akcz)) )
          END IF
        END DO
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE stokes_profile_breivik

END MODULE mo_wave_stokes
