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
!
! Classes and functions for the turbulent mixing package (tmx)
!


!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_tmx_smagorinsky

  USE mo_kind,                ONLY: wp, vp
  USE mo_exception,           ONLY: message
  USE mo_tmx_field_class,     ONLY: t_domain
  USE mo_vdf_atmo_memory,     ONLY: t_vdf_atmo_config, t_vdf_atmo_inputs, t_vdf_atmo_diags
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_sync,                ONLY: SYNC_C, sync_patch_array
  USE mo_physical_constants,  ONLY: grav,rgrav
  USE mo_fortran_tools,       ONLY: init

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Smagorinsky_init, Smagorinsky_model

  REAL(wp) :: eps_louis = 1.0E-28_wp
  !$ACC DECLARE CREATE(eps_louis)

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_smagorinsky'

  CONTAINS
    !============================================================================
    SUBROUTINE Smagorinsky_init(domain,config,inputs,diagnostics)

      TYPE(t_domain),           INTENT(in)    :: domain
      CLASS(t_vdf_atmo_config), INTENT(in)    :: config
      CLASS(t_vdf_atmo_inputs), INTENT(in)    :: inputs
      CLASS(t_vdf_atmo_diags),  INTENT(inout) :: diagnostics

      REAL(wp), POINTER, DIMENSION(:,:,:) :: mixing_length_sq, gepot_agl_ic
      REAL(vp), POINTER, DIMENSION(:,:,:) :: dzh
      REAL(wp), POINTER, DIMENSION(:,:)   :: scaling_factor_louis
      REAL(wp), POINTER :: smag_constant, max_turb_scale
      LOGICAL,  POINTER :: use_louis

      CHARACTER(len=*), PARAMETER :: routine = modname//':Smagorinsky_init'

      CALL message(routine, '')

      mixing_length_sq  => diagnostics%mix_len_sq%Get_ptr_r3d()

      use_louis         => config%use_louis%Get_ptr_l0d()
      smag_constant     => config%smag_constant%Get_ptr_r0d()
      max_turb_scale    => config%max_turb_scale%Get_ptr_r0d()

      dzh               => inputs%dz_ic%get_ptr_v3d()
      gepot_agl_ic      => inputs%geopot_agl_ic%get_ptr_r3d()

      ! Compute mixing length for Smagorinsky model
      CALL compute_mixing_length(domain, dzh, gepot_agl_ic, smag_constant, max_turb_scale, mixing_length_sq)

      ! Computation of variables for stability correction function
      IF (use_louis) THEN

        scaling_factor_louis  => diagnostics%louis_factor%Get_ptr_r2d()

        ! compute scaling_factor_louis in init!
        CALL compute_scaling_factor_louis(domain,scaling_factor_louis)

      END IF

    END SUBROUTINE Smagorinsky_init
    !============================================================================
    !
    ! Computes the square of the SGS mixing length for the Smagorinsky model.
    !
    ! lambda^2 = (Cs * Delta)^2 *(kappa*x_3)^2 / ((Cs * Delta)^2 + (kappa*x_3)^2)
    !          = (Cs * Delta * x_3)^2 / ( (Cs * Delta / kappa)^2 + x_3^2 )
    !
    ! with   Cs    : Smagorinsky constant
    !        Delta : filter/grid width
    !        x_3   : vertical coordinate
    !        kappa : von Karman constant (kappa=0.40)
    !
    ! reference see Dipankar et al. (2015)
    !
    SUBROUTINE compute_mixing_length(domain,dzh, gepot_agl_ic,smag_constant,max_turb_scale,mixing_length_sq)

      TYPE(t_domain), INTENT(in)    :: domain
      REAL(wp),       INTENT(in)    :: smag_constant, max_turb_scale
      REAL(wp), DIMENSION(:,:,:), INTENT(in)    :: gepot_agl_ic
      REAL(vp), DIMENSION(:,:,:), INTENT(in)    :: dzh
      REAL(wp), DIMENSION(:,:,:), INTENT(out)   :: mixing_length_sq

      REAL(wp)  :: kappa, les_filter, z_mc
      INTEGER :: jg, jk, jb, jc
      INTEGER :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
      INTEGER :: nlevp1

      CHARACTER(len=*), PARAMETER :: routine = modname//':compute_mixing_length'

      ! CALL message(routine, '')

      ! von Karman constant
      kappa = 0.4_wp

      jg     = domain%patch%id
      nlevp1 = domain%nlev + 1

      rl_start   = 3
      rl_end     = min_rlcell_int
      i_startblk = domain%patch%cells%start_block(rl_start)
      i_endblk   = domain%patch%cells%end_block(rl_end)

!$OMP PARALLEL
      CALL init(mixing_length_sq, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx, i_endidx, z_mc,les_filter) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(domain%patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) PRIVATE(z_mc, les_filter) ASYNC(1)
        DO jk = 1 , nlevp1
          DO jc = i_startidx, i_endidx
            z_mc  = gepot_agl_ic(jc,jk,jb) * rgrav

            les_filter =  smag_constant           &
                          * MIN( max_turb_scale,  &
                                (dzh(jc,jk,jb) * domain%area(jc,jb))**0.33333_wp &
                               )
            !
            mixing_length_sq(jc,jk,jb) = (les_filter*z_mc)**2._wp    &
                                         / ((les_filter/kappa)**2._wp+z_mc**2._wp)

          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END SUBROUTINE compute_mixing_length
    !============================================================================
    !
    ! Computes the scaling factor for Louis constant b.
    ! Scaling factor for Louis constant b is designed to be 1 with R2B8 setup.
    !
    SUBROUTINE compute_scaling_factor_louis(domain,scaling_factor_louis)
      TYPE(t_domain),           INTENT(in)    :: domain
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: scaling_factor_louis

      INTEGER :: jb, jc
      INTEGER :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx

      ! Global mean of cell area for R2B8 [m]
      REAL(wp), PARAMETER :: mean_area_R2B8 = 97294071.23714285_wp

      CHARACTER(len=*), PARAMETER :: routine = modname//':compute_scaling_factor_louis'

      ! CALL message(routine, '')

      rl_start   = 3
      rl_end     = min_rlcell_int
      i_startblk = domain%patch%cells%start_block(rl_start)
      i_endblk   = domain%patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(domain%patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO jc = i_startidx, i_endidx
          scaling_factor_louis(jc,jb) = mean_area_R2B8 / domain%area(jc,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO

    END SUBROUTINE compute_scaling_factor_louis

    !============================================================================
    !
    ! This subroutine calls the models to compute the eddy viscosity and
    ! diffusivity based on the Smagorinsky-Lilly eddy viscosity model.
    ! Depending on the configuration, the classical version (Lilly 1962) or the
    ! Louis formulation (Louis 1979) for the stability correction function is used
    ! with the option to exclude land cells from the Louis formulation.
    !
    SUBROUTINE Smagorinsky_model( &
      domain,                     &
      mech_prod,                  &
      bruvais,                    &
      rho_ic,                     &
      mixing_length_sq,           &
      rturb_prandtl,              &
      use_louis,                  &
      use_louis_land,             &
      use_louis_ice,              &
      louis_constant_b,           &
      scaling_factor_louis,       &
      fract_land,                 &
      fract_ice,                  &
      km_ic,                      &
      kh_ic                       &
      )

      TYPE(t_domain), INTENT(in)    :: domain
      REAL(wp), INTENT(in), DIMENSION(:,:,:)  :: mech_prod, bruvais, rho_ic, mixing_length_sq
      REAL(wp), INTENT(in), DIMENSION(:,:)    :: scaling_factor_louis, fract_land, fract_ice
      REAL(wp), INTENT(in) :: rturb_prandtl, louis_constant_b
      LOGICAL,  INTENT(in) :: use_louis, use_louis_land, use_louis_ice

      REAL(wp), POINTER, INTENT(inout), DIMENSION(:,:,:) :: km_ic, kh_ic

      REAL(wp) :: stability_term
      INTEGER :: jb, jc, jk, nlev, nlevp1
      INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
      INTEGER :: rl_start, rl_end

      nlev = domain%nlev
      nlevp1 = nlev+1

      rl_start   = 3
      rl_end     = min_rlcell_int
      i_startblk = domain%patch%cells%start_block(rl_start)
      i_endblk   = domain%patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, stability_term) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(domain%patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) &
        !$ACC   PRIVATE(stability_term)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2 , nlev
#else
        DO jk = 2 , nlev
          DO jc = i_startidx, i_endidx
#endif
            IF (use_louis) THEN
              IF  (      (.NOT. use_louis_land .AND. fract_land(jc,jb) > 0.5_wp) &
                &   .OR. (.NOT. use_louis_ice  .AND. fract_ice(jc,jb)  > 0.5_wp) &
                & ) THEN
                ! If Louis formula is used but not over land and/or sea ice, use classic formulation for cells
                ! with more than 50% land fraction or more than 50% ice fraction.
                stability_term = stability_term_classic(mech_prod(jc,jk,jb), bruvais(jc,jk,jb), rturb_prandtl)
              ELSE
                stability_term = stability_term_louis(mech_prod(jc,jk,jb), bruvais(jc,jk,jb), rturb_prandtl, &
                  &  louis_constant_b, scaling_factor_louis(jc,jb))
              END IF
            ELSE
              stability_term = stability_term_classic(mech_prod(jc,jk,jb), bruvais(jc,jk,jb), rturb_prandtl)
            END IF
            km_ic(jc,jk,jb) = rho_ic(jc,jk,jb)               &
                              * mixing_length_sq(jc,jk,jb)   &
                              * stability_term

            kh_ic(jc,jk,jb) = km_ic(jc,jk,jb) * rturb_prandtl
          END DO
        END DO
        !$ACC END PARALLEL LOOP

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO jc = i_startidx, i_endidx
          kh_ic(jc,1,jb)      = kh_ic(jc,2,jb)
          kh_ic(jc,nlevp1,jb) = kh_ic(jc,nlev,jb)
          km_ic(jc,1,jb)      = km_ic(jc,2,jb)
          km_ic(jc,nlevp1,jb) = km_ic(jc,nlev,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO

     !$ACC WAIT(1)

      CALL sync_patch_array(SYNC_C, domain%patch, kh_ic, lacc=.TRUE.)
      CALL sync_patch_array(SYNC_C, domain%patch, km_ic, lacc=.TRUE.)

    END SUBROUTINE Smagorinsky_model
    !
    !============================================================================
    !
    ! This function computes the stability correction term for the eddy viscosity.
    ! The stability term includes the strain rate into the stability
    ! correction function, e.g.
    !
    !   stability term = sqrt(|S|^2 - N^2 / Pr_t )
    !
    !   in  Km = rho * lambda^2 * stability_term
    !
    ! Km     : eddy viscosity
    ! |S|    : magnitude of strain rate  (-> |S|^2 = 0.5*mech_prod)
    ! N^2    : bruvais
    ! Pr_t   : turbulent Prandtl number
    ! lambda : mixing length
    !
#ifndef _OPENACC
    ELEMENTAL &
#endif
    PURE FUNCTION stability_term_classic(mech_prod, bruvais, rturb_prandtl) RESULT(stability_term)

      REAL(wp), INTENT(in) :: rturb_prandtl,  mech_prod, bruvais
      REAL(wp)             :: stability_term

      !$ACC ROUTINE SEQ

      stability_term = SQRT(MAX( 0._wp, 0.5_wp * mech_prod - rturb_prandtl * bruvais ))

    END FUNCTION stability_term_classic
    !
    !============================================================================
    !
    ! This function computes the stability term for the eddy viscosity based on
    ! the stability correction function of Louis (1979).
    ! The stability term includes the strain rate into the stability
    ! correction function, e.g.
    !     Km = rho * lambda^2 * stability_term
    !
    !     - stability_term = sqrt(|S|^2 * stability_factor_louis )
    !
    !     - stability_factor_louis = 1 / (1 + b * Ri )**n
    !
    ! Km        : eddy viscosity
    ! |S|^2     : square of strain rate  (-> |S|^2 = 0.5*mech_prod)
    ! N^2       : bruvais
    ! Pr_t      : turbulent Prandtl number
    ! lambda    : mixing length
    ! Ri        : Richardson number
    ! b         : Louis constant
    !
#ifndef _OPENACC
    ELEMENTAL &
#endif
    PURE FUNCTION stability_term_louis(mech_prod, bruvais, rturb_prandtl, &
      &  louis_constant_b, scaling_factor_louis) RESULT(stability_term)

      REAL(wp), INTENT(in) :: mech_prod, bruvais, rturb_prandtl, louis_constant_b, scaling_factor_louis
      REAL(wp)             :: stability_term

      REAL(wp) :: Ri, stability_function

      !$ACC ROUTINE SEQ

      Ri  = 2._wp * bruvais / MAX(eps_louis, mech_prod)

      stability_function = MAX( 1.0_wp - Ri * rturb_prandtl,                   &
                                MIN(1._wp, (1._wp / (1._wp + louis_constant_b  &
                                                      * scaling_factor_louis   &
                                                      * ABS(Ri)                &
                                            ))**4._wp                          &
                             ))

      stability_term = SQRT( 0.5_wp * mech_prod * stability_function )

    END FUNCTION stability_term_louis

END MODULE mo_tmx_smagorinsky
