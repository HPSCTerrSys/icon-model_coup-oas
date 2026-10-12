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

! Two-moment bulk microphysics by Axel Seifert, Klaus Beheng and Uli Blahak
!
! Description:
! Provides various derived types for two-moment bulk microphysics

MODULE mo_2mom_mcrph_types

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message, txt => message_text

  IMPLICIT NONE

  ! The following parameter is supposed to raise an error, if mistakenly used.
  ! It will be overwritten if initialised.
  REAL(wp), PARAMETER      :: missing =  -HUGE(1.0_wp)

  PUBLIC

  !==================================================================================
  ! Type declarations:
  !==================================================================================

  ! Derived type for atmospheric variables
  TYPE ATMOSPHERE
    REAL(wp), POINTER, DIMENSION(:,:) :: w, p, t, rho, qv, zh, tke
  END TYPE ATMOSPHERE

  ! Derived type for hydrometeor species including pointers to data
  TYPE particle
    CHARACTER(20) :: name       !..name of particle class
    REAL(wp)      :: nu         !..first shape parameter of size distribution
    REAL(wp)      :: mu         !..2nd shape parameter
    REAL(wp)      :: x_max      !..max mean particle mass
    REAL(wp)      :: x_min      !..min mean particle mass
    REAL(wp)      :: a_geo      !..pre-factor in diameter-mass relation
    REAL(wp)      :: b_geo      !..exponent in diameter-mass relation
    REAL(wp)      :: a_vel      !..pre-factor in power law fall speed (all particles have a power law fall speed,
    REAL(wp)      :: b_vel      !..exponent in power law fall speed    some have also an Atlas-type relation)
    REAL(wp)      :: a_ven      !..first parameter in ventilation coefficient
    REAL(wp)      :: b_ven      !..2nd parameter in ventilation coefficient
    REAL(wp)      :: cap        !..coefficient for capacity of particle
    REAL(wp)      :: vsedi_max  !..max bulk sedimentation velocity
    REAL(wp)      :: vsedi_min  !..min bulk sedimentation velocity
    REAL(wp), POINTER, DIMENSION(:,:) :: n     !..number density
    REAL(wp), POINTER, DIMENSION(:,:) :: q     !..mass density
    REAL(wp), POINTER, DIMENSION(:,:) :: rho_v !..density correction of terminal fall velocity

    REAL(wp)      :: ecoll_c    !..maximum collision efficiency with cloud droplets
    REAL(wp)      :: D_crit_c   !..D-threshold for cloud riming
    REAL(wp)      :: q_crit_c   !..q-threshold for cloud riming
    REAL(wp)      :: s_vel      !..dispersion of fall velocity for collection kernel (see SB2006, Eqs 60-63)
    ! The following are needed only needed for lwf
    REAL(wp)      :: lwf_cnorm1  !..1st parameter for normalized diameter
    REAL(wp)      :: lwf_cnorm2  !..2nd parameter for normalized diameter
    REAL(wp)      :: lwf_cnorm3  !..3rd parameter for normalized diameter
    REAL(wp)      :: lwf_cmelt1  !..1st parameter for melting integral
    REAL(wp)      :: lwf_cmelt2  !..2nd parameter for melting integral
    REAL(wp), pointer, dimension(:,:) :: l  !..mass density of liquid water on ice (per unit volume of air)
  END TYPE particle

  ! .. Because of OpenMP we have to separate the data pointers from the run-time-invariant coefficients.
  !    Therefore we carry 2 data structures for each particle species, e.g. graupel and graupel_coeff.
  !    The following derived types are for the run-time coefficients

  TYPE particle_coeffs
    REAL(wp)      :: a_f  ! ventilation coefficient, vent_coeff_a(particle,1)
    REAL(wp)      :: b_f  ! ventilation coefficient, vent_coeff_b(particle,1) * N_sc**n_f / SQRT(nu_l)
    REAL(wp)      :: c_i  ! 1.0/particle%cap
    REAL(wp)      :: c_z  ! coefficient for 2nd mass moment

  ! .. for spherical particles we need to store the coefficients for the
  !    power law bulk sedimentation velocity
    REAL(wp)      :: coeff_alfa_n = missing
    REAL(wp)      :: coeff_alfa_q = missing
    REAL(wp)      :: coeff_lambda = missing

  ! .. non-spherical particles have an Atlas-type terminal fall velocity relation
    REAL(wp)      :: alfa = missing   !..1st parameter in Atlas-type fall speed
    REAL(wp)      :: beta = missing   !..2nd parameter in Atlas-type fall speed
    REAL(wp)      :: gama = missing   !..3rd parameter in Atlas-type fall speed

  ! .. raindrops have an Atlas-type terminal fall velocity relation
  !    and a mu-D-relation which is used in sedimentation and evaporation
  !    (see Seifert 2008, J. Atmos. Sci.)
    REAL(wp)      :: cmu0 = missing   !..Parameters for mu-D-relation of rain: max of left branch
    REAL(wp)      :: cmu1 = missing   !     max of right branch
    REAL(wp)      :: cmu2 = missing   !     scaling factor
    REAL(wp)      :: cmu3 = missing   !     location of min value = breakup equilibrium diameter
    REAL(wp)      :: cmu4 = missing   !     min value of relation
    INTEGER       :: cmu5             !     exponent

    REAL(wp)      :: k_au = missing   !..Parameters for autoconversion
    REAL(wp)      :: k_sc = missing   !    and selfcollection

    REAL(wp)      :: sc_coll_n  = missing !..Parameters for self-collection
    REAL(wp)      :: sc_delta_n = missing !..Parameters for self-collection
    REAL(wp)      :: sc_theta_n = missing !   of snow
    REAL(wp)      :: sc_delta_q = missing !   of cloud ice
    REAL(wp)      :: sc_theta_q = missing
  END TYPE particle_coeffs

  TYPE aerosol_ccn
     REAL(wp)      :: Ncn0      ! CN concentration at ground
     REAL(wp)      :: Nmin      ! minimum value for CCN
     REAL(wp)      :: lsigs     ! log(sigma_s)
     REAL(wp)      :: R2        ! in mum
     REAL(wp)      :: etas      ! soluble fraction
     REAL(wp)      :: wcb_min   ! min updraft speed for Segal&Khain nucleation
     REAL(wp)      :: z0        ! parameter for height-dependency, constant up to z0_nccn
     REAL(wp)      :: z1e       ! 1/e scale height of N_ccn profile
  END TYPE aerosol_ccn

  TYPE aerosol_in
     REAL(wp)      :: N0     ! CN concentration at ground
     REAL(wp)      :: z0     ! parameter for height-dependency, constant up to z0_nccn
     REAL(wp)      :: z1e    ! 1/e scale height of N_ccn profile
  END TYPE aerosol_in

    !..these are coefficients for collection processes of the type a+b->a
  TYPE collection_coeffs
     REAL(wp) :: delta_n_aa, delta_n_ab, delta_n_bb, &
          &      delta_q_aa, delta_q_ab, delta_q_bb, &
          &      theta_n_aa, theta_n_ab, theta_n_bb, &
          &      theta_q_aa, theta_q_ab, theta_q_bb
  END TYPE collection_coeffs

  !..these are coefficients for collection processes of the type a+b->c
  TYPE rain_riming_coeffs
    REAL(wp) :: delta_n_aa,delta_n_ab,delta_n_bb, &
         &      delta_q_aa,delta_q_ab,delta_q_ba,delta_q_bb, &
         &      theta_n_aa,theta_n_ab,theta_n_bb, &
         &      theta_q_aa,theta_q_ab,theta_q_ba,theta_q_bb
  END TYPE rain_riming_coeffs

  !..these are coefficients for partial moments collection of type a+b->a
  TYPE coll_coeffs_ir_pm
    REAL(wp) :: moma(0:1,1:5), momb(0:1,1:5)
    REAL(wp) :: delta_aa(0:1,0:1), delta_ab(0:1,0:1), delta_bb(0:1,0:1)
    REAL(wp) :: theta_aa(0:1,0:1), theta_ab(0:1,0:1), theta_bb(0:1,0:1)
    REAL(wp) :: lamfakt_a, lamfakt_b
  END TYPE coll_coeffs_ir_pm

  TYPE dep_imm_coeffs
    REAL(wp) :: alf_dep, bet_dep, nin_dep, &
                alf_imm, bet_imm, nin_imm
  END TYPE dep_imm_coeffs

  !..these are coefficients for diagnostics, which should be pre-computed at
  !   the initial stage:
  TYPE t_diag_coeffs_2mom
    ! for hail kinetic energy flux including horizontal transport:
    REAL(wp) :: coeff_v3_kef_hail, coeff_v2_kef_hail, coeff_v1_kef_hail
  END TYPE t_diag_coeffs_2mom

  ! Type declaration for a general 1D equidistant lookup table:
  TYPE lookupt_1D
    LOGICAL :: is_initialized = .FALSE.
    INTEGER :: iflag = -HUGE(1) ! general-purpose flag
    CHARACTER(len=40) :: name   ! general-purpose name
    INTEGER :: n1  ! number of grid points in x1-direction
    REAL(wp), DIMENSION(:), POINTER :: x1 => NULL()   ! grid vector in x1-direction
    REAL(wp)                        :: dx1   ! dx1   (grid distance w.r.t. x1)
    REAL(wp)                        :: odx1  ! one over dx1
    REAL(wp), DIMENSION(:), POINTER :: ltable => NULL()
  END TYPE lookupt_1D

  ! Type declaration for a general 2D equidistant lookup table:
  TYPE lookupt_2D
    LOGICAL :: is_initialized = .FALSE.
    INTEGER :: iflag = -HUGE(1) ! general-purpose flag
    CHARACTER(len=40) :: name   ! general-purpose name
    INTEGER :: n1  ! number of grid points in x1-direction
    INTEGER :: n2  ! number of grid points in x2-direction
    REAL(wp), DIMENSION(:), POINTER :: x1 => NULL()   ! grid vector in x1-direction
    REAL(wp), DIMENSION(:), POINTER :: x2 => NULL()   ! grid vector in x1-direction
    REAL(wp)                        :: dx1   ! dx1   (grid distance w.r.t. x1)
    REAL(wp)                        :: dx2   ! dx2   (grid distance w.r.t. x2)
    REAL(wp)                        :: odx1  ! one over dx1
    REAL(wp)                        :: odx2  ! one over dx2
    REAL(wp), DIMENSION(:,:), POINTER :: ltable => NULL()
  END TYPE lookupt_2D

  ! Type declaration for a general 4D equidistant lookup table:
  TYPE lookupt_4D
    LOGICAL :: is_initialized = .FALSE.
    INTEGER :: iflag = -HUGE(1) ! general-purpose flag
    CHARACTER(len=40) :: name   ! general-purpose name
    INTEGER :: n1  ! number of grid points in x1-direction
    INTEGER :: n2  ! number of grid points in x2-direction
    INTEGER :: n3  ! number of grid points in x3-direction
    INTEGER :: n4  ! number of grid points in x4-direction
    REAL(wp), DIMENSION(:), POINTER :: x1 => NULL()  ! grid vector in x1-direction
    REAL(wp), DIMENSION(:), POINTER :: x2 => NULL()  ! grid vector in x1-direction
    REAL(wp), DIMENSION(:), POINTER :: x3 => NULL()  ! grid vector in x1-direction
    REAL(wp), DIMENSION(:), POINTER :: x4 => NULL()  ! grid vector in x1-direction
    REAL(wp)                     :: dx1          ! dx1   (grid distance w.r.t. x1)
    REAL(wp)                     :: dx2          ! dx2   (grid distance w.r.t. x2)
    REAL(wp)                     :: dx3          ! dx3   (grid distance w.r.t. x3)
    REAL(wp)                     :: dx4          ! dx4   (grid distance w.r.t. x4)
    REAL(wp)                     :: odx1         ! one over dx1
    REAL(wp)                     :: odx2         ! one over dx2
    REAL(wp)                     :: odx3         ! one over dx3
    REAL(wp)                     :: odx4         ! one over dx4
    REAL(wp), DIMENSION(:,:,:,:), POINTER :: ltable => NULL()
  END TYPE lookupt_4D

  CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_types'

  !==================================================================================
  ! Actual instances of types:
  !==================================================================================

  ! Types to hold the equidistant lookup table for graupel wetgrowth diameter:
  TYPE(lookupt_4d), TARGET :: ltabdminwgg
  ! Types to hold the equidistant lookup table for hail wetgrowth diameter:
  TYPE(lookupt_4d), TARGET :: ltabdminwgh

  ! Types to hold the equidistant lookup table for sticking efficiencies:
  TYPE(lookupt_1d), TARGET :: ltab_estick_ice
  TYPE(lookupt_1d), TARGET :: ltab_estick_snow
  TYPE(lookupt_1d), TARGET :: ltab_estick_parti

  !==================================================================================
  !==================================================================================

  PRIVATE :: routine

END MODULE mo_2mom_mcrph_types
