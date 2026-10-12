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

MODULE mo_hamocc_swr_absorption

  USE mo_param1_bgc, ONLY : iphy, icya
  USE mo_hamocc_nml, ONLY : l_cyadyn
  USE mo_kind,    ONLY: wp
  USE mo_control_bgc, ONLY: bgc_zlevs, bgc_nproma
  USE mo_bgc_memory_types, ONLY  : t_bgc_memory
  USE mo_fortran_tools, ONLY     : set_acc_host_or_device
  USE mo_exception, ONLY      : finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: swr_absorption


CONTAINS

SUBROUTINE swr_absorption(local_bgc_mem, start_idx,end_idx, klevs, pfswr, psicomo, dzw, lacc)

    TYPE(t_bgc_memory), POINTER :: local_bgc_mem
    INTEGER, INTENT(in):: start_idx
    INTEGER, INTENT(IN):: end_idx
    INTEGER :: klevs(bgc_nproma)
    REAL(wp), INTENT(in):: pfswr(bgc_nproma)
    REAL(wp), INTENT(in):: psicomo(bgc_nproma)
    REAL(wp), INTENT(in):: dzw(bgc_nproma,bgc_zlevs)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    !! Analogue to Zielinski et al., Deep-Sea Research II 49 (2002), 3529-3542

    REAL(wp), PARAMETER :: redfrac=0.4_wp !< red fraction of the spectral domain (> 580nm)

    REAL(wp), PARAMETER :: c_to_chl=12.0_wp/60.0_wp   !< ration Carbon to Chlorophyll
    REAL(wp), PARAMETER :: r_car=122.0_wp   !< Redfield ratio
    REAL(wp), PARAMETER :: pho_to_chl=r_car*c_to_chl*1.e6_wp !< 1 kmolP = (122*12/60)*10^6 mg[Chlorophyll]

    REAL(wp), PARAMETER :: atten_r=0.35_wp !< attenuation of red light [m-1]
    REAL(wp), PARAMETER :: atten_w=0.03_wp !< attenuation of blue/green light
                                           !! in clear water between 400nm and 580nm [m-1]
    REAL(wp), PARAMETER :: atten_c=0.04_wp !< attenuation of blue/green light
                                           !! by chlorophyll [m-1]

#ifndef __LVECTOR__
    REAL(wp) :: swr_r
    REAL(wp) :: swr_b
#else
    REAL(wp) :: swr_r(start_idx:end_idx)
    REAL(wp) :: swr_b(start_idx:end_idx)
#endif

    REAL(wp) :: rcyano

    INTEGER :: k, kpke, j, max_klevs
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

#if defined(__LVECTOR__) && defined(_OPENACC)
    IF (lzacc) CALL finish("", "LVECTOR variant after reworking not properly ported/tested on GPUs")
#endif

    ! if prognostic cyanobacteria are calculated
    ! use them in absorption (rcyano=1)
    rcyano=merge(1._wp,0._wp,l_cyadyn)

#ifndef __LVECTOR__

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO j = start_idx, end_idx

      local_bgc_mem%strahl(j) = pfswr(j) * (1._wp - psicomo(j))

      local_bgc_mem%swr_frac(j,1) = 1.0_wp


      kpke = klevs(j)

      IF(kpke > 0) then

      swr_r = redfrac
      swr_b = (1._wp-redfrac)
      !$ACC LOOP SEQ
      DO k=2,kpke

           swr_r = swr_r * EXP(-dzw(j,k-1) *  atten_r)
           swr_b = swr_b * EXP(-dzw(j,k-1) * (atten_w +&
        &    atten_c*pho_to_chl*MAX(0.0_wp,(local_bgc_mem%bgctra(j,k-1,iphy)+rcyano*local_bgc_mem%bgctra(j,k-1,icya)))))
           local_bgc_mem%swr_frac(j,k) = swr_r + swr_b

      END DO
      !$ACC LOOP SEQ
      DO k=1,kpke-1
           local_bgc_mem%meanswr(j,k) = (local_bgc_mem%swr_frac(j,k) + local_bgc_mem%swr_frac(j,k+1))/2._wp
      END DO
      local_bgc_mem%meanswr(j,kpke) = local_bgc_mem%swr_frac(j,k)

      ENDIF
    ENDDO
    !$ACC END PARALLEL

#else

    max_klevs = MAXVAL(klevs(start_idx:end_idx))

    !NEC$ nomove
    DO j = start_idx, end_idx
        local_bgc_mem%strahl(j) = pfswr(j) * (1._wp - psicomo(j))
        local_bgc_mem%swr_frac(j,1) = 1.0_wp

        IF(klevs(j) > 0) THEN
            swr_r(j) = redfrac
            swr_b(j) = (1._wp-redfrac)
        END IF
    END DO

    DO k = 2, max_klevs
        !NEC$ nomove
        DO j = start_idx, end_idx

            IF(k <= klevs(j)) THEN
                swr_r(j) = swr_r(j) * EXP(-dzw(j,k-1) * &
                       &   atten_r)

                swr_b(j) = swr_b(j) * EXP(-dzw(j,k-1) * &
                       &  (atten_w + atten_c * pho_to_chl * MAX(0.0_wp,(local_bgc_mem%bgctra(j,k-1,iphy) + rcyano * local_bgc_mem%bgctra(j,k-1,icya)))))

                local_bgc_mem%swr_frac(j,k) = swr_r(j) + swr_b(j)
            END IF

        END DO
    END DO

    DO k = 1, max_klevs - 1
        !NEC$ nomove
        DO j = start_idx, end_idx

            IF(k <= klevs(j) - 1) THEN
                local_bgc_mem%meanswr(j,k) = (local_bgc_mem%swr_frac(j,k) + local_bgc_mem%swr_frac(j,k+1))/2._wp
            END IF

        END DO
    END DO

    !NEC$ nomove
    DO j = start_idx, end_idx
        IF(klevs(j) > 0) THEN
            local_bgc_mem%meanswr(j,klevs(j)) = local_bgc_mem%swr_frac(j,klevs(j))
        END IF
    END DO

#endif

END SUBROUTINE swr_absorption
END MODULE
