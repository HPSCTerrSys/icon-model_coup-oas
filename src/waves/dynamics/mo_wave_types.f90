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

! Contains the types to set up the wave model.

MODULE mo_wave_types

  USE mo_kind,                ONLY: wp, vp
  USE mo_var_list,            ONLY: t_var_list_ptr
  USE mo_fortran_tools,       ONLY: t_ptr_2d3d, t_ptr_2d3d_vp, t_ptr_3d_int

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_wave_prog
  PUBLIC :: t_wave_source
  PUBLIC :: t_wave_const
  PUBLIC :: t_wave_diag_dyn
  PUBLIC :: t_wave_diag_out
  PUBLIC :: t_wave_diag_cpl
  PUBLIC :: t_wave_state
  PUBLIC :: t_wave_state_lists
  PUBLIC :: t_wesd

  ! derived type for wave energy spectral density
  !
  TYPE t_wesd
    REAL(wp), POINTER, CONTIGUOUS :: & ! wave energy spectral density for fixed frequency [m^2 s]
    ptr(:,:,:) => NULL()               ! over all directions and cells (nproma,ndirs,nblks_c)

    TYPE(t_ptr_2d3d), ALLOCATABLE :: dir(:) !< pointer array for direction bins
  END TYPE t_wesd


  TYPE t_wave_prog
    !
    TYPE(t_wesd), ALLOCATABLE :: wesd(:)  ! wave energy spectral density array [m^2 s]
                                          ! wesd(nfreqs)%ptr(nproma,ndirs,nblks_c)
  END TYPE t_wave_prog


  ! source function state vector object
  !
  TYPE t_wave_source
    REAL(vp), POINTER, CONTIGUOUS :: &
      &  sl(:,:,:,:),         & ! total source function                    (nproma,ndirs,nfreqs,nblks_c) (-)
      &  fl(:,:,:,:)            ! diagonal matrix of functional derivative (nproma,ndirs,nfreqs,nblks_c) (-)

    INTEGER, POINTER, CONTIGUOUS ::  &
      &  llws(:,:,:,:)          ! 1 - where sinput is positive (nproma,ndirs,nfreqs,nblks_c) (-)

    TYPE(t_ptr_2d3d_vp), ALLOCATABLE :: sl_ptr(:)   !< pointer array: one pointer for each frequency
    TYPE(t_ptr_2d3d_vp), ALLOCATABLE :: fl_ptr(:)   !< pointer array: one pointer for each frequency
    TYPE(t_ptr_3d_int),  ALLOCATABLE :: llws_ptr(:) !< pointer array: one pointer for each frequency
  END TYPE t_wave_source


  ! diagnostic variables required for time integration (affect solution)
  ! Validity: every time step
  !
  TYPE t_wave_diag_dyn
    REAL(wp), POINTER, CONTIGUOUS :: &
      &  gv_c(:,:,:),         & ! group velocity                             (nproma,nfreqs,nblks_c)  (m/s)
      &  gv_e(:,:,:),         & ! group velocity                             (nproma,nfreqs,nblks_e)  (m/s)
      &  wave_num_c(:,:,:),   & ! wave number at cell centers as a function
                                ! of circular frequency and water depth      (nproma,nfreqs,nblks_c) (1/m)
      &  wave_num_e(:,:,:),   & ! wave number at cell edges as a function of
                                ! circular frequency and water depth         (nproma,nfreqs,nblks_e) (1/m)
      &  akmean(:,:),         & ! mean wavenumber based on sqrt(1/k)-moment  (nproma,nblks_c) (1/m)
      &  xkmean(:,:),         & ! mean wavenumber based on sqrt(k)-moment    (nproma,nblks_c) (1/m)
      &  ustar(:,:),          & ! friction velocity                          (nproma,nblks_c) (m/s)
      &  z0(:,:),             & ! roughness length                           (nproma,nblks_c) (m)
      &  tauhf(:,:),          & ! high-frequency stress divided by air density (nproma,nblks_c) (m/s)^2
      &  phihf(:,:),          & ! high-frequency energy flux into waves      (nproma,nblks_c) (W m-2)
      &  tauw(:,:),           & ! wave stress divided by air density         (nproma,nblks_c) (m/s)^2
      &  phiaw(:,:),          & ! energy flux from wind into waves           (nproma,nblks_c) (W m-2)
                                ! integrated over the full frequency range
      ! total waves
      &  emean(:,:),          & ! total energy                               (nproma,nblks_c) (m^2)
      &  emeanws(:,:),        & ! total wind sea input energy                (nproma,nblks_c) (m^2)
      &  femean(:,:),         & ! mean frequency (m0/m-1)                    (nproma,nblks_c) (s-1)
      &  f1mean(:,:),         & ! mean frequency based on f-moment (m1/m0)   (nproma,nblks_c) (s-1)
      &  femeanws(:,:),       & ! wind sea mean frequency (m0/m-1)           (nproma,nblks_c) (s-1)
      &  tm1(:,:),            & ! total wave m1 period                       (nproma,nblks_c) (s)
      &  tm2(:,:),            & ! total wave m2 period                       (nproma,nblks_c) (s)
      &  hrms_frac(:,:),      & ! square ratio (Hrms / Hmax)**2              (nproma,nblks_c) (-)
      &  wbr_frac(:,:)          ! fraction of breaking waves                 (nproma,nblks_c) (-)

    INTEGER, POINTER, CONTIGUOUS ::  &
      &  last_prog_freq_ind(:,:)! last frequency index of the prognostic range (nproma,nblks_c) (-)

  END type t_wave_diag_dyn


  ! diagnostic variables for output only (no impact on solution)
  ! Validity: at output times
  !
  TYPE t_wave_diag_out
    REAL(wp), POINTER, CONTIGUOUS :: &
      &  beta(:,:),           & ! Charnock parameter               (nproma,nblks_c) (-)
      &  drag(:,:),           & ! drag coefficient                 (nproma,nblks_c) (-)
      &  ds(:,:),             & ! total directional wave spread    (nproma,nblks_c) (deg)
      &  hs_dir(:,:),         & ! total mean wave direction        (nproma,nblks_c) (deg)
      &  kp(:,:),             & ! total peak wavenumber            (nproma,nblks_c) (m-1)
      &  tauwn(:,:),          & ! normalised wave stress           (nproma,nblks_c) (-)
      &  tmp(:,:),            & ! total mean wave period           (nproma,nblks_c) (s)
      &  tpp(:,:),            & ! total peak wave period           (nproma,nblks_c) (s)
      !
      ! wind sea
      &  ds_sea(:,:),         & ! sea directional spreed           (nproma,nblks_c) (deg)
      &  emean_sea(:,:),      & ! wind sea energy                  (nproma,nblks_c) (m^2)
      &  femean_sea(:,:),     & ! wind sea mean frequency (m0/m-1) (nproma,nblks_c) (s-1)
      &  f1mean_sea(:,:),     & ! wind sea mean frequency (m1/m0)  (nproma,nblks_c) (s-1)
      &  hs_sea(:,:),         & ! sea significant wave height      (nproma,nblks_c) (m)
      &  hs_sea_dir(:,:),     & ! sea mean wave direction          (nproma,nblks_c) (deg)
      &  kp_sea(:,:),         & ! sea peak wavenumber              (nproma,nblks_c) (m-1)
      &  mp_sea(:,:),         & ! sea mean period                  (nproma,nblks_c) (s)
      &  m1_sea(:,:),         & ! sea m1-period                    (nproma,nblks_c) (s)
      &  m2_sea(:,:),         & ! sea m2-period                    (nproma,nblks_c) (s)
      &  pp_sea(:,:),         & ! sea peak period                  (nproma,nblks_c) (s)
      !
      ! swell
      &  ds_swell(:,:),       & ! swell directional spreed         (nproma,nblks_c) (deg)
      &  emean_swell(:,:),    & ! swell energy                     (nproma,nblks_c) (m^2)
      &  femean_swell(:,:),   & ! swell mean frequency (m0/m-1)    (nproma,nblks_c) (s-1)
      &  f1mean_swell(:,:),   & ! swell mean frequency (m1/m0)     (nproma,nblks_c) (s-1)
      &  hs_swell(:,:),       & ! swell significant wave height    (nproma,nblks_c) (m)
      &  hs_swell_dir(:,:),   & ! swell mean wave direction        (nproma,nblks_c) (deg)
      &  kp_swell(:,:),       & ! swell peak wavenumber            (nproma,nblks_c) (m-1)
      &  mp_swell(:,:),       & ! swell mean period                (nproma,nblks_c) (s)
      &  m1_swell(:,:),       & ! swell m1-period                  (nproma,nblks_c) (s)
      &  m2_swell(:,:),       & ! swell m2-period                  (nproma,nblks_c) (s)
      &  pp_swell(:,:),       & ! swell peak period                (nproma,nblks_c) (s)
      !
      ! extreme wave parameters
      &  steepness(:,:),      & ! integral wave steepness        (nproma,nblks_c) (-)
      &  thp_adj(:,:),        & ! adjusted peak wave direction   (nproma,nblks_c) (deg)
      &  sigma_f(:,:),        & ! frequency peak width           (nproma,nblks_c) (-)
      &  sigma_th(:,:),       & ! directional peak width         (nproma,nblks_c) (-)
      &  qp_goda(:,:),        & ! Goda peakedness parameter      (nproma,nblks_c) (-)
      &  nu_f_LH(:,:),        & ! Longuet-Higgins width parameter(nproma,nblks_c) (-)
      &  bfis(:,:),           & ! intermediate-depth BFI         (nproma,nblks_c) (-)
      &  relw_fth(:,:),       & ! relative direction-freq width  (nproma,nblks_c) (-)
      &  kurtosis(:,:),       & ! spectral wave kurtosis         (nproma,nblks_c) (-)
      &  hmaxn(:,:),          & ! normalised max wave height     (nproma,nblks_c) (-)
      &  hmax(:,:),           & ! maximum wave height            (nproma,nblks_c) (m)
      &  Tmax(:,:)              ! maximum wave period            (nproma,nblks_c) (s)

    INTEGER, POINTER, CONTIGUOUS ::  &
      &  swell_mask(:,:,:,:)       ! swell separation mask (nproma,ndirs,nfreqs,nblks_c) (-)

    TYPE(t_ptr_3d_int), ALLOCATABLE :: swmask_ptr(:)     !< pointer array: one pointer for each variance bin
  END TYPE t_wave_diag_out

  ! diagnostic variables for coupling and output, if requested
  ! (no impact on wave solution, but may impact coupled components)
  ! Validity: every time step in coupled runs
  !           at output times (if requested) in standalone runs
  !
  TYPE t_wave_diag_cpl
    REAL(wp), POINTER, CONTIGUOUS :: &
      &  hs(:,:),             & ! total significant wave height       (nproma,nblks_c)         (m)
      &  u_stokes(:,:),       & ! U-component of surface Stokes drift (nproma,nblks_c)         (m/s)
      &  v_stokes(:,:),       & ! V-component of surface Stokes drift (nproma,nblks_c)         (m/s)
      &  u3d_stokes(:,:,:),   & ! U-component of 3d Stokes drift      (nproma,ndepths,nblks_c) (m/s)
      &  v3d_stokes(:,:,:),   & ! V-component of 3d Stokes drift      (nproma,ndepths,nblks_c) (m/s)
      &  phioc(:,:),          & ! wave-to-ocean energy flux           (nproma,nbkls_c)         (kg/s^3)
      &  tauoc_x(:,:),        & ! zonal component of wave-to-ocean stress (nproma,nblks_c)     (m/s)^2
      &  tauoc_y(:,:),        & ! meridional component of wave-to-ocean stress (nproma,nblks_c) (m/s)^2
      &  tauoc(:,:)             ! magnitude of wave-to-ocean stress   (nproma,nblks_c)         (m/s)^2

    INTEGER, POINTER, CONTIGUOUS ::  &
      &  last_idx_depth(:,:)    ! last index of Stokes depth layer    (nproma,nblks_c)         (-)
  END TYPE t_wave_diag_cpl


  ! time independent fields and coefficients which are calculated once
  ! during initialization.
  !
  TYPE t_wave_const
    REAL(wp), POINTER, CONTIGUOUS :: &
      &  flminfr_tab(:,:)       ! minimum value of energy for a given
                                ! frequency and wind speed bin                     (jmax,nfreqs)

    REAL(vp), POINTER, CONTIGUOUS :: &
      &  AF11(:),                 & ! for discrete approximation of nonlinear transfer (nfreqs+4) (-)
      &  FKLAP(:), FKLAP1(:),     & ! --//-- (nfreqs+4) (-)
      &  FKLAM(:), FKLAM1(:)        ! --//-- (nfreqs+4) (-)

    INTEGER, POINTER, CONTIGUOUS ::  &
      &  ikp(:), ikp1(:),         & ! for discrete approximation of nonlinear transfer (nfreqs+4) (-)
      &  ikm(:), ikm1(:),         & ! --//-- (nfreqs+4) (-)
      &  k1w(:,:), k2w(:,:),      & ! --//-- (ndirs, 2) (-)
      &  k11w(:,:), k21w(:,:),    & ! --//-- (ndirs, 2) (-)
      &  ja1(:,:), ja2(:,:)         ! --//-- (ndirs, 2) (-)
  END TYPE t_wave_const

  !
  ! Grouping of diagnostic fields
  !

  ! List     |   Validity      | Used internally | Used by coupled components | Computation trigger
  !========================================================================================================
  !          |                 |                 |                            |
  ! diag_dyn |  always         |       yes       |       maybe                | every timestep
  ! diag_cpl |  at coupling    |       no        |       yes                  | coupling schedule
  !          |  optional: at   |                 |                            | optional: output schedule
  !          |  output         |                 |                            |
  ! diag_out |  at output      |       no        |       no                   | output schedule
  ! const    |  always         |       yes       |       no                   | at initialization step

  TYPE t_wave_state
    !
    !< array of prognostic states at different timelevels
    TYPE(t_wave_prog), ALLOCATABLE    :: prog(:)
    !
    ! < source function state vector
    TYPE(t_wave_source)               :: source
    !
    !< coefficients for model integration
    !  computed once during model initialization
    TYPE(t_wave_const)                :: const
    !
    !< diagnostics required for time integration (affect solution)
    TYPE(t_wave_diag_dyn)             :: diag_dyn
    !
    !< diagnostics required for coupling and optionally for output.
    !  Computed at coupling times and at output times, if requested.
    !  No impact on solution, but may affect other components.
    TYPE(t_wave_diag_cpl)             :: diag_cpl
    !
    !< diagnostics for output only (no impact on solution)
    TYPE(t_wave_diag_out)             :: diag_out
  END TYPE t_wave_state


  TYPE t_wave_state_lists
    ! array of prognostic state lists at different timelevels
    TYPE(t_var_list_ptr), ALLOCATABLE :: prog_list(:)  !< shape: (timelevels)
    TYPE(t_var_list_ptr)              :: source_list
    TYPE(t_var_list_ptr)              :: const_list
    TYPE(t_var_list_ptr)              :: diag_dyn_list
    TYPE(t_var_list_ptr)              :: diag_cpl_list
    TYPE(t_var_list_ptr)              :: diag_out_list
  END TYPE t_wave_state_lists

END MODULE mo_wave_types
