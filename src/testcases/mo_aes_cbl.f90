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

! Module containing subroutines for the initialization of a
! convective boundary layer (CBL) on a torus

MODULE mo_aes_cbl
#ifndef __NO_AES__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish, print_value
  USE mo_model_domain,        ONLY: t_patch
  USE mo_physical_constants,  ONLY: rd, cvd, cpd, p0ref, cvd_o_rd, rd_o_cpd, &
     &                              grav, alv, vtmpc1, lh_v=>alv
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics, t_nh_ref
  USE mo_parallel_config,     ONLY: nproma
  USE mo_nh_testcases_nml,    ONLY: p_sfc, th_sfc, gamma, &
                                    nlev_pert, th_perturb, w_perturb, &
                                    itheta_init, t_cbl_sol,           &
                                    isrfc_type, shflx, lhflx
  USE mo_run_config,          ONLY: iqv, iqc
  USE mo_aes_thermo,          ONLY: specific_humidity, sat_pres_water
  USE mo_aes_vdf_config      ,ONLY: aes_vdf_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_aes_cbl_dry, print_aes_cbl_testcase_config

CONTAINS
  !
  !============================================================================
  !
  SUBROUTINE init_aes_cbl_dry( ptr_patch, ptr_nh_prog, ptr_nh_ref, ptr_nh_diag, ptr_metrics)

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,  INTENT(IN)   :: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_nh_prog),       INTENT(INOUT):: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT):: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   :: &  !< NH metrics state
      &  ptr_metrics
    TYPE(t_nh_ref),        INTENT(INOUT):: &  !< reference state vector
      &  ptr_nh_ref

    REAL(wp) :: z_exner_h(1:nproma,ptr_patch%nlev+1)

    INTEGER  :: jb,jk,jl,nblks_c,npromz_c,nlen,nlev,nlevp1  !<loop indices and control

    REAL(wp), ALLOCATABLE :: temp(:,:,:), rh(:,:,:)
    REAL(wp)              :: theta(1:nproma), theta_init
    REAL(wp)              :: x, y, z, dz, Kh, theta_sfc, theta_flx

    REAL(wp), PARAMETER :: zh0     = 0._wp      !< height (m) above which temperature increases

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    ALLOCATE(temp(nproma,nlev,nblks_c))
    ALLOCATE(rh(nproma,nlev,nblks_c))

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = p_sfc

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      END IF

      IF (itheta_init == 1) THEN
        ! Initialize linear potential temperature profile
        DO jl = 1, nlen
          DO jk = nlev, 1, -1
            z = ptr_metrics%z_mc(jl,jk,jb)-zh0

            CALL linear_theta_v_profile(z, th_sfc, gamma, theta_init)

            ptr_nh_prog%theta_v(jl,jk,jb) = theta_init
          END DO
        END DO

      ELSE IF (itheta_init == 2) THEN
        ! Initialize potential temperature profile with analytical
        ! solution for diffusion problem of theta
        DO jl = 1, nlen
          DO jk = nlev, 1, -1
            z = ptr_metrics%z_mc(jl,jk,jb)-zh0
            Kh = aes_vdf_config(1)%km_const * aes_vdf_config(1)%rturb_prandtl
            theta_flx = shflx * cvd / cpd

            CALL analytical_solution_cbl_const_sflx(z, t_cbl_sol, th_sfc, theta_flx, gamma, Kh, theta_init)

            ptr_nh_prog%theta_v(jl,jk,jb) = theta_init
          END DO
        END DO
      ELSE
        CALL finish('testcases/mo_aes_cbl.f90', 'Initial theta profile for CBL test case only accepts itheta_init == 1 or 2')
      END IF

      !Get hydrostatic exner at the surface using surface pressure
      z_exner_h(1:nlen,nlevp1) = (p_sfc/p0ref)**rd_o_cpd

      !Get exner at full levels starting from exner at surface
      ! from formula dPi/dz = -g/cpd/theta_v
      DO jk = nlev, 1, -1
        !exner at next half level after surface
        z_exner_h(1:nlen,jk) = z_exner_h(1:nlen,jk+1) - grav/cpd *    &
                               ptr_metrics%ddqz_z_full(1:nlen,jk,jb)/ &
                               ptr_nh_prog%theta_v(1:nlen,jk,jb)

        !exner at main levels
        ptr_nh_prog%exner(1:nlen,jk,jb) = 0.5_wp * (z_exner_h(1:nlen,jk)+z_exner_h(1:nlen,jk+1))
        IF ( ptr_nh_prog%exner(1,jk,jb)<0.0_wp ) THEN
          CALL finish('testcases/mo_aes_cbl.f90', 'Exner function is negative')
        ENDIF
      END DO ! end of jk

      DO jk = 1 , nlev
        ptr_nh_prog%rho(1:nlen,jk,jb) = (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)
        ! Mistake in following formula???
        ! Shouldn't it be p = rho*Rd*Tv with Tv = Exner * theta_v ???
        !ptr_nh_diag%pres(1:nlen,jk,jb) = ptr_nh_prog%rho(1:nlen,jk,jb)*rd*th_cbl(1)
        ptr_nh_diag%pres(1:nlen,jk,jb) = rd * ptr_nh_prog%rho(1:nlen,jk,jb) &
                                        * ptr_nh_prog%exner(1:nlen,jk,jb)   &
                                        * ptr_nh_prog%theta_v(1:nlen,jk,jb)
      END DO !jk

    END DO ! end of jb


    !meridional and zonal wind
    ptr_nh_prog%vn = 0._wp
    ptr_nh_ref%vn_ref = ptr_nh_prog%vn

    !vertical wind
    ptr_nh_prog%w = 0._wp
    ptr_nh_ref%w_ref = ptr_nh_prog%w

    ! Setting water vapour to zero
    ptr_nh_prog%tracer(:, :, :, iqv) = 0.0_wp

  END SUBROUTINE init_aes_cbl_dry
  !
  !============================================================================
  !
  ! This subroutine imposes  a linear profile for potential temperature based
  ! on the surface potential temperature (th_sfc) and lapse rate (gamma).
  !
  SUBROUTINE linear_theta_v_profile(z, theta_sfc, gamma, theta_init)

    REAL(wp), INTENT(IN) :: z, theta_sfc, gamma
    REAL(wp), INTENT(OUT) :: theta_init

    theta_init = theta_sfc + gamma*z

  END SUBROUTINE linear_theta_v_profile
  !
  !============================================================================
  !
  ! This subroutine computes the analytical solution for the diffusion problem
  ! of potential temperature (d theta/dt = Kh d^2 theta / d z^2), with constant
  ! eddy diffusivity, for a convective boundary layer setting. The boundary
  ! condition at the bottom is a constant surface flux.
  !
  SUBROUTINE analytical_solution_cbl_const_sflx(z, t_start, theta_sfc, shflx, gamma, Kh, theta_sol)

    REAL(wp), INTENT(IN) :: z, t_start, theta_sfc, shflx, gamma, Kh
    REAL(wp), INTENT(OUT) :: theta_sol

    REAL(wp) :: h, eta, g, PI

    PI=4._wp*ATAN(1._wp)

    h = sqrt(4._wp*Kh*t_start)
    eta = z/h

    g = exp(-eta*eta)/sqrt(PI) - eta * ( 1._wp - erf(eta) )

    theta_sol = theta_sfc + gamma*z + (shflx/Kh + gamma) * h * g

  END SUBROUTINE analytical_solution_cbl_const_sflx
  !
  !============================================================================
  !
  SUBROUTINE print_aes_cbl_testcase_config()

    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','aes_cbl testcase configuration')
    CALL message    ('','======================================')
    CALL message    ('','')

    CALL print_value('    p_sfc:          ',p_sfc         )
    CALL print_value('    th_sfc:         ',th_sfc        )
    CALL print_value('    gamma:          ',gamma         )

    CALL print_value('    isrfc_type:     ',isrfc_type    )
    CALL print_value('    shflx:          ',shflx         )
    CALL print_value('    lhflx:          ',lhflx         )

    CALL print_value('    nlev_pert:      ',nlev_pert     )
    CALL print_value('    th_perturb:     ',th_perturb    )
    CALL print_value('    w_perturb:      ',w_perturb     )

    CALL print_value('    itheta_init:    ',itheta_init   )
    CALL print_value('    t_cbl_sol:      ',t_cbl_sol     )

  END SUBROUTINE print_aes_cbl_testcase_config
!!!=============================================================================================

#endif
END MODULE mo_aes_cbl
