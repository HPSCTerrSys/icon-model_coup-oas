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

! This module contains parameters for the solar constant and the
! routine for computing the cosine of the solar zenith angle, both
! together needed to define the incident solar flux.
!
! The underlying data are part of the standard CMPI6 input, referenced in
!
! Matthes, K., Funke, B., Anderson, M. E., Barnard, L., Beer, J.,
! Charbonneau, P., Clilverd, M. A., Dudok de Wit, T., Haberreiter, M.,
! Hendry, A., Jackman, C. H., Kretschmar, M., Kruschke, T., Kunze, M.,
! Langematz, U., Marsh, D. R., Maycock, A., Misios, S., Rodger, C. J.,
! Scaife, A. A., Seppala, A., Shangguan, M., Sinnhuber, M., Tourpali,
! K., Usoskin, I., van de Kamp, M., Verronen, P. T., and S. Versick,
! 2017: Solar Forcing for CMIP6 (v3.2). Geosci. Model Dev., 10,
! doi:10.5194/gmd-10-2247-2017.
!
! and provided by HEPPA/SOLARIS Geomar.

MODULE mo_solar_parameters

  USE mo_kind,   ONLY: wp
  USE mo_model_domain,    ONLY: t_patch
  USE mo_parallel_config, ONLY: nproma
  USE mo_impl_constants,  ONLY: max_dom
  USE mo_math_constants,  ONLY: pi, pi2, pi_2, pi_4 ! pi, pi*2, pi/2, pi/4

IMPLICIT NONE

  PRIVATE

  PUBLIC :: solar_parameters

CONTAINS


  !---------------------------------------------------------------------------
  !>
  !! @brief Scans a block and fills with solar parameters
  !!
  !! @remarks: This routine calculates the solar zenith angle for each
  !! point in a block of data.  For simulations with no diurnal cycle
  !! the cosine of the zenith angle is set to its average value (assuming
  !! negatives to be zero and for a day divided into nds intervals).
  !! Additionally a field is set indicating the fraction of the day over
  !! which the solar zenith angle is greater than zero.  Otherwise the field
  !! is set to 1 or 0 depending on whether the zenith angle is greater or
  !! less than 1.
  !
  SUBROUTINE solar_parameters(decl_sun,    time_of_day,     &
       &                      icosmu0,     dt_ext,          &
       &                      ldiur,       l_sph_symm_irr,  &
       &                      p_patch,                      &
       &                      cos_mu0,     daylight_frc,    &
       &                      lacc                          )

    REAL(wp), INTENT(in)  :: &
         decl_sun,           & !< delination of the sun
         time_of_day,        & !< time_of_day (in radians)
         dt_ext                !< time interval overfor which the insolated area is extended
    INTEGER , INTENT(in)  :: &
         icosmu0               !< defines if and how cosmu0 is defined for extended sunlit areas
    LOGICAL               :: &
         ldiur,              & !< diurnal cycle ON (ldiur=.TRUE.) or OFF (ldiur=.FALSE.)
         l_sph_symm_irr        !< spherical symmetric irradiation ON (l_sph_symm_irr=.TRUE.)
                               !< or OFF (l_sph_symm_irr=.FALSE.)
    TYPE(t_patch), INTENT(in) ::      p_patch
    REAL(wp), INTENT(out) :: &
         cos_mu0(:,:),       & !< cos_mu_0, cosine of the solar zenith angle
         daylight_frc(:,:)     !< daylight fraction (0 or 1) with diurnal cycle
    LOGICAL, INTENT(in)   :: lacc


    INTEGER     :: i, j
    REAL(wp)    :: zen1, zen2, zen3, xx

    INTEGER, PARAMETER :: nds = 128 !< number of diurnal samples

    LOGICAL  , SAVE :: initialized_sincos = .FALSE.
    REAL (wp), SAVE :: cosrad(nds), sinrad(nds)
    REAL (wp)       :: xsmpl(nds), xsmpl_day(nds), xnmbr_day(nds)
    REAL (wp), ALLOCATABLE :: sinlon(:,:), sinlat(:,:), coslon(:,:), coslat(:,:)
    REAL (wp), ALLOCATABLE :: mu0(:,:)

    REAL (wp), PARAMETER   :: eps=1.e-9_wp ! to converge in ca. 30 iterations
    LOGICAL  , SAVE        :: initialized_mu0s = .FALSE.
    REAL (wp), SAVE        :: mu0s, sin_mu0s
    REAL (wp)              :: dmu0, dcos_mu0, mu0min, mu0max
    INTEGER  , SAVE        :: nmu0

    INTEGER :: nprom, npromz, nblks

    !$ACC DECLARE CREATE(cosrad, sinrad)

    nprom=nproma
    npromz=p_patch%npromz_c
    nblks=p_patch%nblks_c

    IF (.NOT.ldiur.AND..NOT.initialized_sincos) THEN
       !
       ! sin and cos arrays for zonal mean radiation
       !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(xx) ASYNC(1) IF(lacc)
       DO i = 1, nds
          xx = pi2*(i-1.0_wp)/nds
          sinrad(i) = SIN(xx)
          cosrad(i) = COS(xx)
       END DO
       !$ACC END PARALLEL LOOP
       !
       initialized_sincos = .TRUE.
    END IF
    !
    IF (ldiur.AND.(icosmu0==4).AND.(.NOT.initialized_mu0s).AND.(dt_ext/=0.0_wp)) THEN
       !
       ! find mu0s for radiation with diurnal cycle
       dmu0   = dt_ext/86400._wp*pi
       mu0min = 0.0_wp
       mu0max = pi_2
       mu0s   = mu0max
       nmu0   = 0
       DO
          xx = (pi_2+dmu0-mu0s)*SIN(mu0s) - COS(mu0s)
          IF      (xx >  eps ) THEN
             mu0max = mu0s
          ELSE IF (xx < -eps ) THEN
             mu0min = mu0s
          ELSE
             EXIT  ! ABS(xx)<=eps, mu0s is found
          END IF
          mu0s = 0.5_wp*(mu0min+mu0max)
          nmu0 = nmu0+1
          !
       END DO
       !
       ! slope of the linear section of cos_mu0, see below, where used
       sin_mu0s = SIN(mu0s)
       !
       initialized_mu0s = .TRUE.
    END IF
    !
    zen1 = SIN(decl_sun)
    zen2 = COS(decl_sun)*COS(time_of_day)
    zen3 = COS(decl_sun)*SIN(time_of_day)

    IF (.NOT.l_sph_symm_irr) THEN       ! - spherically variable irradiation
       !
       ALLOCATE(sinlon(nprom,nblks))
       ALLOCATE(sinlat(nprom,nblks))
       ALLOCATE(coslon(nprom,nblks))
       ALLOCATE(coslat(nprom,nblks))
       ALLOCATE(mu0(nprom,nblks))

       !$ACC DATA CREATE(sinlon, sinlat, coslon, coslat, mu0) IF(lacc)

       !$ACC PARALLEL LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lacc)
       do j=1,nblks-1
         do i=1,nprom
           sinlon(i,j)=SIN(p_patch%cells%center(i,j)%lon)
           sinlat(i,j)=SIN(p_patch%cells%center(i,j)%lat)
           coslon(i,j)=COS(p_patch%cells%center(i,j)%lon)
           coslat(i,j)=COS(p_patch%cells%center(i,j)%lat)
         end do
       end do
       !$ACC END PARALLEL LOOP

       !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
       !$ACC LOOP GANG(STATIC: 1) VECTOR
       do i=1,npromz
         sinlon(i,nblks)=SIN(p_patch%cells%center(i,nblks)%lon)
         sinlat(i,nblks)=SIN(p_patch%cells%center(i,nblks)%lat)
         coslon(i,nblks)=COS(p_patch%cells%center(i,nblks)%lon)
         coslat(i,nblks)=COS(p_patch%cells%center(i,nblks)%lat)
       end do
       !$ACC END PARALLEL

       !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
       !$ACC LOOP GANG(STATIC: 1) VECTOR
       do i=npromz+1,nprom
         sinlon(i,nblks)=0.0_wp
         sinlat(i,nblks)=0.0_wp
         coslon(i,nblks)=0.0_wp
         coslat(i,nblks)=0.0_wp
       end do
       !$ACC END PARALLEL

       !
       IF (ldiur) THEN                  ! - with local diurnal cycle
          !
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          ! cos(zenith angle), positive for sunlit hemisphere
          do j=1,nblks
            do i=1,nprom
              cos_mu0(i,j)     =  zen1*sinlat(i,j)                &
                 &             -zen2*coslat(i,j)*coslon(i,j)    &
                 &             +zen3*coslat(i,j)*sinlon(i,j)
              !
              ! zenith angle
              mu0(i,j) = ACOS(cos_mu0(i,j))
            end do
          end do
          !$ACC END PARALLEL LOOP
          !
          ! increment of mu0 to include a rim of width dmu0= (dt_ext/2)*2pi/1day around
          ! the sunlit hemisphere at the radiation time so that the extended area includes
          ! the sunlit areas of all time step within the radiation interval of length dt_ext
          ! centered at the radiation time
          dmu0     = dt_ext/86400._wp*pi
          dcos_mu0 = SIN(dmu0)
          !
          ! add increment dcos_mu0 for the definition of the extended daylight area
          ! set day/night indicator to 1/0
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          do j=1,nblks
            do i=1,nprom
              daylight_frc(i,j) = 1.0_wp
              if (cos_mu0(i,j)+dcos_mu0 < 0.0_wp) then
                daylight_frc(i,j) = 0.0_wp
              end if
            end do
          end do
          !$ACC END PARALLEL LOOP
          !
          IF (dt_ext/=0.0_wp) THEN
             !

             SELECT CASE (icosmu0)
             CASE (0)
                !
                ! no modification -> nothing to do
                !
             CASE (1)
                !
                ! minimum value
               !
                !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lacc)
                do j=1,nblks
                  do i=1,nprom
                    if (daylight_frc(i,j) == 1.0_wp) then
                      cos_mu0(i,j) = MAX(0.1_wp,cos_mu0(i,j))
                    end if
                  end do
                end do
                !$ACC END PARALLEL LOOP
                !
             CASE (2)
                !
                ! shift and rescale
                !

                !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lacc)
                do j=1,nblks
                  do i=1,nprom
                    if (daylight_frc(i,j) == 1.0_wp) then
                      cos_mu0(i,j) = (cos_mu0(i,j)+dcos_mu0)/(1._wp+dcos_mu0)
                    end if
                  end do
                end do
                !$ACC END PARALLEL LOOP
                !
             CASE (3)
                !
                ! slope in band [pi/2-dmu0,pi/2+dmu0]
                !
                ! redefine cos_mu0 in [pi/2-dmu0,pi/2+dmu0] using a linear function of mu0 such that:
                !   inner edge                       : mu0=pi/2-dmu0 --> cos_mu0 = cos(pi/2-dmu0)
                !   center     = original terminator : mu0=pi/2      --> cos_mu0 = cos(pi/2-dmu0)/2
                !   outer edge = extended terminator : mu0=pi/2+dmu0 --> cos_mu0 = 0
                !
                !$ACC PARALLEL LOOP DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(lacc)
                do j=1,nblks
                  do i=1,nproma
                    if (abs(mu0(i,j)-pi_2)<dmu0) then
                      cos_mu0(i,j) = 0.5_wp*sin(dmu0)*(1._wp-(mu0(i,j)-pi_2)/dmu0)
                    end if
                  end do
                end do
                !$ACC END PARALLEL LOOP
                !
             CASE (4)
                !
                ! tangent slope in [mu0s,pi/2+dmu0]
                !
                ! redefine cos_mu0 in [mu0s,pi/2+dmu0] using a linear function of mu0 such that:
                !   inner edge                       : mu0=mu0s      --> cos_mu0 = cos(mu0s)
                !   in between                       : mu0           --> cos_mu0 = sin(mu0s)*(pi/2+dmu0-mu0)
                !   outer edge = extended terminator : mu0=pi/2+dmu0 --> cos_mu0 = 0
                !
                ! mu0s is the solution of : cos(mu0s) = sin(mu0s)*(pi/2+dmu0-mu0s)
                ! so that cos_mu0 is a C1 function in [0,pi/2+dmu0]
                !
                !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lacc)

                do j=1,nblks
                  do i=1,nproma
                    if ((mu0s<mu0(i,j)).AND.(mu0(i,j)<(pi_2+dmu0))) then
                      cos_mu0(i,j) = sin_mu0s*(pi_2+dmu0-mu0(i,j))
                    end if
                  end do
                end do
                !$ACC END PARALLEL LOOP
                !
             END SELECT
          END IF
          !
       ELSE                             ! - with zonally symmetric irradiation
          !
          ! For each cell (i,j), compute first cos(mu0) for nds regularly spaced longitudes on the
          ! latitude circle of the cell. Then compute the zonal mean of cos(mu0) for the latitude
          ! of this cell as the average over all longitudes, where cos(mu0)>=epsilon.
          !
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(lacc)
          DO j = 1, SIZE(cos_mu0,2)
            DO i = 1, SIZE(cos_mu0,1)
                !
                ! cos(zenith angle) for nds longitudes on the latitude circle of cell (i,j)
                xsmpl(:) =  zen1*sinlat(i,j)              &
                     &     -zen2*coslat(i,j)*cosrad(:)    &
                     &     +zen3*coslat(i,j)*sinrad(:)
                !
                ! set day/night indicator of sampled longitudes to 1/0
                ! and mask out cosmu0 values in day-only and night-only cells
                !
                xsmpl_day(:) = 0.0_wp
                xnmbr_day(:) = 0.0_wp
                WHERE (xsmpl(:) >= EPSILON(1.0_wp))     ! night side
                   xsmpl_day(:) = xsmpl(:)
                   xnmbr_day(:) = 1.0_wp
                END WHERE
                !
                cos_mu0(i,j)      = SUM(xsmpl_day(:))
                daylight_frc(i,j) = SUM(xnmbr_day(:))
                !
                IF (daylight_frc(i,j) > EPSILON(1.0_wp)) THEN         ! at least one point on the latitude
                   cos_mu0(i,j)      = cos_mu0(i,j)/daylight_frc(i,j) ! circle of the cell has daylight
                   daylight_frc(i,j) = daylight_frc(i,j)/REAL(nds,wp) ! and a zonal mean can be defined
                ELSE
                   cos_mu0(i,j)      = 0.0_wp                         ! set cos(mu0) = 0 on polar night
                   daylight_frc(i,j) = 0.0_wp                         ! latitude circle
                END IF
                !
             END DO
          END DO
          !$ACC END PARALLEL LOOP
          !
       END IF
       !
       !$ACC WAIT(1)
       !$ACC END DATA
       DEALLOCATE(sinlon)
       DEALLOCATE(sinlat)
       DEALLOCATE(coslon)
       DEALLOCATE(coslat)
       DEALLOCATE(mu0)
       !
    ELSE                                ! - spherically symmetric irradiation (for RCE),
       !                                    all grid points have the same solar incoming
       !                                    flux at TOA
       !
       IF (ldiur) THEN                  ! - with diurnal cycle of (0degE,0degN), i.e.
          !                                 local noon is at 12:00 UTC in all points
          !
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          do j=1,nblks
            do i=1,nproma
              cos_mu0(i,j) = -zen2          !   = cos_mu0 at (0degE,0degN)
              IF (-zen2 < 0.0_wp) THEN
                daylight_frc(i,j) = 0.0_wp
              ELSE
                daylight_frc(i,j) = 1.0_wp
              END IF
            end do
          end do
          !$ACC END PARALLEL LOOP
          !
       ELSE                             ! - without diurnal cycle
          !                                 all grid points have the same constant cos_mu0
          !
          ! cos_mu0(:,:) = pi_4           !  = pi/4 (why this choice?)
          ! quickhack for RCEMIP_analytical
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          do j=1,nblks
            do i=1,nproma
              cos_mu0(i,j) = COS(42.05_wp*pi/180._wp)           !  = pi/4 (why this choice?)
              daylight_frc(i,j) = 1.0_wp
            end do
          end do
          !$ACC END PARALLEL LOOP
          !
       END IF
       !
       !
    END IF

    !$ACC WAIT(1)
  END SUBROUTINE solar_parameters

END MODULE mo_solar_parameters
