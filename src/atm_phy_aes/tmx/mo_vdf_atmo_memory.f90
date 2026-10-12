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

! Memory management for the vertical diffusion atmospheric component
!
! This module provides the memory management structures and routines for
! the vertical diffusion component of the atmospheric model in TMX. It defines
! the input, configuration, and diagnostic data structures required for
! vertical diffusion calculations and provides routines for their initialization.
!
! The module contains three main data structures:
! - t_vdf_atmo_config: Configuration parameters for vertical diffusion
! - t_vdf_atmo_inputs: Input meteorological fields required by the diffusion scheme
! - t_vdf_atmo_diags: Diagnostic fields produced by the diffusion scheme

MODULE mo_vdf_atmo_memory

  USE mo_kind,            ONLY: wp
  USE mo_tmx_var,         ONLY: t_tmx_var
  USE mo_tmx_field_class, ONLY: t_domain
  USE mo_run_config,      ONLY: ntracer

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_vdf_atmo_config, build_vdf_atmo_config
  PUBLIC :: t_vdf_atmo_inputs, build_vdf_atmo_inputs
  PUBLIC :: t_vdf_atmo_diags, build_vdf_atmo_diags, t_vel_grad_tensor

  ! Configuration parameters for atmospheric vertical diffusion
  !
  ! This type contains all the physical constants, turbulence parameters,
  ! and configuration settings needed for vertical diffusion calculations.
  TYPE t_vdf_atmo_config
    ! Physics constants and configuration parameters
    TYPE(t_tmx_var) :: &
      & cpd              , & !< Specific heat capacity at constant pressure for dry air [J/(kg K)]
      & cvd              , & !< Specific heat capacity at constant volume for dry air [J/(kg K)]
      & smag_constant    , & !< Smagorinsky constant []
      & max_turb_scale   , & !< Maximum turbulence length scale [m]
      & rturb_prandtl    , & !< Reverse Prandtl number (1/Pr) []
      & turb_prandtl     , & !< Prandtl number (Pr) []
      & km_min           , & !< Minimum exchange coefficient for momentum [m2/s]
      & k_s              , & !< Surface exchange coefficient []
      & louis_constant_b , & !< Louis constant b []
      & km_const         , & !< Constant turbulent diffusivity value [m2/s]
      & scale_turb_energy_flux, & !< Scaling factor for turbulent energy flux []
      & dtime            , & !< Physics time step [s]
      & dissipation_factor    !< Factor for kinetic energy dissipation calculation []
    TYPE(t_tmx_var) :: &
      & use_louis        , & !< Switch to activate Louis stability formula []
      & use_louis_land   , & !< Switch to activate Louis formula over land (default: true) []
      & use_louis_ice    , & !< Switch to activate Louis formula over sea ice (default: true) []
      & use_km_const     , & !< Switch to use constant turbulent diffusivity []
      & use_scale_turb_energy_flux, & !< Switch to scale turbulent energy flux []
      & solver_type      , & !< Type of solver (1=explicit, 2=implicit) []
      & energy_type      , & !< Type of energy calculation []
      & l_co2                !< Use CO2 tracer for vertical diffusion []
  END TYPE t_vdf_atmo_config

  ! Input fields required for vertical diffusion calculations
  !
  ! This type contains all the meteorological fields needed as inputs for
  ! the vertical diffusion calculations. Fields are organized by their
  ! vertical position (full or half levels) and horizontal location.
  TYPE t_vdf_atmo_inputs
    ! Full levels on cells
    TYPE(t_tmx_var) ::  &
      & tracer_c     ,  & !< Tracer concentrations [kg/kg]
      & temp_c       ,  & !< Air temperature [K]
      & u_wind_c     ,  & !< Zonal wind component [m/s]
      & v_wind_c     ,  & !< Meridional wind component [m/s]
      & geo_height_c ,  & !< Geopotential height [m]
      & rho_c        ,  & !< Air density [kg/m3]
      & temp_virt_c  ,  & !< Virtual temperature [K]
      & pres_c       ,  & !< Atmospheric pressure [Pa]
      & moist_mass_c ,  & !< Mass of moist air in grid cell [kg/m2]
      & cv_air_c     ,  & !< Specific heat capacity at constant volume [J/(kg K)]
      & z_c          ,  & !< Height of full levels [m]
      & dz_c         ,  & !< Vertical grid spacing at full levels [m]
      & inv_dz_c     ,  & !< Inverse of vertical grid spacing (1/dz) [1/m]
      & fract_land   ,  & !< Fraction of land in grid cell [] (for use_louis_land=false only)
      & fract_ice         !< Fraction of sea ice in grid cell [] (for use_louis_ice=false only)

    ! Half levels on cells (interfaces)
    TYPE(t_tmx_var) ::  &
      w_wind_ic      ,  & !< Vertical wind component at half levels [m/s]
      pres_ic        ,  & !< Pressure at half levels [Pa]
      geo_height_ic  ,  & !< Geopotential height at half levels [m]
      geopot_agl_ic  ,  & !< Geopotential above ground level at half levels [m2/s2]
      z_ic           ,  & !< Height of half levels [m]
      dz_ic          ,  & !< Vertical grid spacing at half levels [m]
      inv_dz_ic           !< Inverse of vertical grid spacing at half levels [1/m]

    ! Full levels on edges
    ! TYPE(t_tmx_var) ::  &
    !   vn_e                !< Normal component of wind at cell edges [m/s]
  END TYPE t_vdf_atmo_inputs

  ! Velocity gradient tensor structure for storing 3D gradient components
  !
  ! This type provides pointers to the components of the velocity gradient
  ! tensor in the normal, tangential and vertical directions.
  TYPE t_vel_grad_tensor
    REAL(wp), DIMENSION(:,:,:), POINTER :: ptr => NULL() ! Pointer to a specific gradient component
  END TYPE t_vel_grad_tensor

  ! Diagnostic fields produced by vertical diffusion calculations
  !
  ! This type contains all diagnostic fields produced by the
  ! vertical diffusion process, including exchange coefficients,
  ! stability parameters, velocity components, and energy terms.
  ! Fields are organized by location (cells, edges, vertices) and
  ! vertical position (full or half levels).
  TYPE t_vdf_atmo_diags
    LOGICAL :: is_initialized = .FALSE. !< Flag indicating if the structure is initialized

    ! 3D diagnostics on cells (full levels)
    TYPE(t_tmx_var) ::  &
      & ghf          ,  & !< Full level geopotential height above ground [m]
      & ctgz         ,  & !< Static energy [m2/s2]
      & div_c        ,  & !< Divergence on cell centers []
      & theta_v     ,  & !< Virtual potential temperature [K]
      & pprfac       ,  & !< Air density divided by layer thickness [kg/m4]
      & km           ,  & !< Exchange coefficient for momentum [m2/s]
      & kh           ,  & !< Exchange coefficient for heat [m2/s]
      & km_c         ,  & !< Exchange coefficient for momentum at cell centers [m2/s]
      & heating      ,  & !< Layer heating rate [W/m2]
      & dissip_ke         !< Dissipation of kinetic energy [W/m2]

    ! 3D diagnostics on cell interfaces (half levels)
    TYPE(t_tmx_var) ::  &
      & rho_ic       ,  & !< Air density at interfaces [kg/m3]
      & bruvais      ,  & !< Brunt-Vaisala frequency [1/s]
      & stab_func    ,  & !< Stability function []
      & mech_prod    ,  & !< Mechanical production of turbulence [m2/s3]
      & km_ic        ,  & !< Exchange coefficient for momentum at interfaces [m2/s]
      & kh_ic        ,  & !< Exchange coefficient for heat at interfaces [m2/s]
      & mix_len_sq        !< Square of turbulent mixing length [m2]

    ! 3D diagnostics on edges (full levels)
    TYPE(t_tmx_var) ::  &
      & vn           ,  & !< Normal wind component [m/s]
      & shear        ,  & !< Shear rate [1/s]
      & div_stress   ,  & !< Stress divergence [m/s2]
      & vn_ie        ,  & !< Normal wind at edge interfaces [m/s]
      & vt_ie        ,  & !< Tangential wind at edge interfaces [m/s]
      & w_ie         ,  & !< Vertical wind at edge interfaces [m/s]
      & km_ie             !< Exchange coefficient for momentum at edge interfaces [m2/s]

    ! 3D diagnostics at vertices
    TYPE(t_tmx_var) ::  &
      & u_vert       ,  & !< Zonal wind at vertices [m/s]
      & v_vert       ,  & !< Meridional wind at vertices [m/s]
      & w_vert       ,  & !< Vertical wind at vertices [m/s]
      & km_iv             !< Exchange coefficient for momentum at vertex interfaces [m2/s]

    ! 2D vertically integrated diagnostics
    TYPE(t_tmx_var) ::  &
      & ctgzvi       ,  & !< Vertically integrated static energy [J/m2]
      & dissip_ke_vi ,  & !< Vertically integrated dissipation of kinetic energy [W/m2]
      & int_energy_vi,  & !< Vertically integrated moist internal energy [J/m2]
      & int_energy_vi_tend   !< Tendency of vertically integrated moist internal energy [W/m2]

    ! Surface flux diagnostics
    TYPE(t_tmx_var) ::  &
      & louis_factor ,  & !< Scaling factor for Louis constant b []
      & lhflx        ,  & !< Latent heat flux at surface [W/m**2]
      & shflx             !< Sensible heat flux at surface [W/m**2]

    ! Velocity gradient tensor components at edges
    TYPE(t_tmx_var) ::  &
      & dvn_dn       ,  & !< Normal gradient of normal wind [1/s]
      & dvt_dn       ,  & !< Normal gradient of tangential wind [1/s]
      & dw_dn        ,  & !< Normal gradient of vertical wind [1/s]
      & dvn_dt       ,  & !< Tangential gradient of normal wind [1/s]
      & dvt_dt       ,  & !< Tangential gradient of tangential wind [1/s]
      & dw_dt        ,  & !< Tangential gradient of vertical wind [1/s]
      & dvn_dz       ,  & !< Vertical gradient of normal wind [1/s]
      & dvt_dz       ,  & !< Vertical gradient of tangential wind [1/s]
      & dw_dz             !< Vertical gradient of vertical wind [1/s]

    ! Velocity gradient tensor on edges
    !
    ! This 3x3 tensor stores the components of the velocity gradient
    ! on edges. The components are arranged as follows:
    ! (1,1): dvn_dn - normal gradient of normal wind
    ! (1,2): dvn_dt - tangential gradient of normal wind
    ! (1,3): dvn_dz - vertical gradient of normal wind
    ! (2,1): dvt_dn - normal gradient of tangential wind
    ! (2,2): dvt_dt - tangential gradient of tangential wind
    ! (2,3): dvt_dz - vertical gradient of tangential wind
    ! (3,1): dw_dn - normal gradient of vertical wind
    ! (3,2): dw_dt - tangential gradient of vertical wind
    ! (3,3): dw_dz - vertical gradient of vertical wind
    TYPE(t_vel_grad_tensor) :: vel_grad_e(3,3)

  END TYPE t_vdf_atmo_diags

  ! Module name for error messages and debugging
  CHARACTER(LEN=*), PARAMETER :: modname="mo_vdf_atmo_memory"

  CONTAINS

  ! Builds and initializes the vdf_atmo_config data structure
  !
  ! This subroutine initializes all configuration parameters in the t_vdf_atmo_config
  ! structure with appropriate dimensions and attributes. The structure includes
  ! physical constants like specific heat capacities, turbulence parameters like
  ! the Smagorinsky constant and mixing length, as well as solver control parameters.
  SUBROUTINE build_vdf_atmo_config(this, domain)

    TYPE(t_vdf_atmo_config), INTENT(inout) :: this    !< The vdf_atmo_config structure to be initialized
    TYPE(t_domain),                POINTER :: domain  !< Domain specification containing grid dimensions

    INTEGER :: patch_id
    INTEGER :: shape_0d(0) ! Zero-sized array for scalar configuration parameters

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"build_vdf_atmo_config"

    ! Extract information from domain
    patch_id = domain%id      ! Patch ID for the domain

    ! Initialize configuration fields
    CALL this%cpd              %Init('cpd',                'double',  dims=shape_0d, patch_id=patch_id)
    CALL this%cvd              %Init('cvd',                'double',  dims=shape_0d, patch_id=patch_id)
    CALL this%smag_constant    %Init('smag_constant',      'double',  dims=shape_0d, patch_id=patch_id)
    CALL this%max_turb_scale   %Init('max_turb_scale',     'double',  dims=shape_0d, patch_id=patch_id)
    CALL this%rturb_prandtl    %Init('rturb_prandtl',      'double',  dims=shape_0d, patch_id=patch_id)
    CALL this%turb_prandtl     %Init('turb_prandtl',       'double',  dims=shape_0d, patch_id=patch_id)
    CALL this%km_min           %Init('km_min',             'double',  dims=shape_0d, patch_id=patch_id)
    CALL this%k_s              %Init('k_s',                'double',  dims=shape_0d, patch_id=patch_id)
    CALL this%use_louis        %Init('use_louis',          'logical', dims=shape_0d, patch_id=patch_id)
    CALL this%use_louis_land   %Init('use_louis_land',     'logical', dims=shape_0d, patch_id=patch_id)
    CALL this%use_louis_ice    %Init('use_louis_ice',      'logical', dims=shape_0d, patch_id=patch_id)
    CALL this%louis_constant_b %Init('louis_constant_b',   'double',  dims=shape_0d, patch_id=patch_id)
    CALL this%use_km_const     %Init('use_km_const',       'logical', dims=shape_0d, patch_id=patch_id)
    CALL this%km_const         %Init('km_const',           'double',  dims=shape_0d, patch_id=patch_id)
    CALL this%use_scale_turb_energy_flux%Init('use_scale_turb_energy_flux', 'logical',  dims=shape_0d, patch_id=patch_id)
    CALL this%scale_turb_energy_flux%Init('scale_turb_energy_flux', 'double',  dims=shape_0d, patch_id=patch_id)
    CALL this%dtime            %Init('dtime',              'double',  dims=shape_0d, patch_id=patch_id)
    CALL this%dissipation_factor%Init('dissipation_factor','double',  dims=shape_0d, patch_id=patch_id)
    CALL this%solver_type      %Init('solver_type',        'integer', dims=shape_0d, patch_id=patch_id)
    CALL this%energy_type      %Init('energy_type',        'integer', dims=shape_0d, patch_id=patch_id)
    CALL this%l_co2            %Init('l_co2',              'logical', dims=shape_0d, patch_id=patch_id)

  END SUBROUTINE build_vdf_atmo_config

  ! Builds and initializes the vdf_atmo_inputs data structure
  !
  ! This subroutine initializes all input fields required for vertical diffusion
  ! calculations. The fields include atmospheric state variables like temperature,
  ! wind components, density, and pressure, as well as grid metrics like height
  ! and layer thickness at both full and half levels.
  SUBROUTINE build_vdf_atmo_inputs(this, domain)

    TYPE(t_vdf_atmo_inputs), INTENT(inout) :: this    !< The vdf_atmo_inputs structure to be initialized
    TYPE(t_domain),                POINTER :: domain  !< Domain specification containing grid dimensions

    INTEGER :: nproma, nblks_c, nblks_e, nlev, nlevp1, patch_id
    INTEGER :: shape_2d(2), shape_3d(3), shape_4d(4)

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"build_vdf_atmo_inputs"

    ! Extract information from domain
    nproma   = domain%nproma  ! Number of grid points per block
    nblks_c  = domain%nblks_c ! Number of blocks for cells
    nblks_e  = domain%nblks_e ! Number of blocks for edges
    nlev     = domain%nlev    ! Number of vertical levels
    nlevp1   = nlev + 1       ! Number of half levels
    patch_id = domain%id      ! Patch ID for the domain

    ! Initialize tracer field (4D)
    shape_4d = [nproma,nlev,nblks_c,ntracer]
    CALL this%tracer_c     %Init('tracer_c',      'double', dims=shape_4d, patch_id=patch_id)

    ! Initialize 3D fields at full levels on cells
    shape_3d = [nproma,nlev,nblks_c]
    CALL this%temp_c       %Init('temp_c',        'double', dims=shape_3d, patch_id=patch_id)
    CALL this%u_wind_c     %Init('u_wind_c',      'double', dims=shape_3d, patch_id=patch_id)
    CALL this%v_wind_c     %Init('v_wind_c',      'double', dims=shape_3d, patch_id=patch_id)
    CALL this%geo_height_c %Init('geo_height_c',  'double', dims=shape_3d, patch_id=patch_id)
    CALL this%rho_c        %Init('rho_c',         'double', dims=shape_3d, patch_id=patch_id)
    CALL this%temp_virt_c  %Init('temp_virt_c',   'double', dims=shape_3d, patch_id=patch_id)
    CALL this%pres_c       %Init('pres_c',        'double', dims=shape_3d, patch_id=patch_id)
    CALL this%moist_mass_c %Init('moist_mass_c',  'double', dims=shape_3d, patch_id=patch_id)
    CALL this%cv_air_c     %Init('cv_air_c',      'double', dims=shape_3d, patch_id=patch_id)
    CALL this%z_c          %Init('z_c',           'double', dims=shape_3d, patch_id=patch_id)
    CALL this%dz_c         %Init('dz_c',          'double', dims=shape_3d, patch_id=patch_id)
    CALL this%inv_dz_c     %Init('inv_dz_c',      'mixed' , dims=shape_3d, patch_id=patch_id)

    ! Initialize 3D fields at half levels on cells
    shape_3d = [nproma,nlevp1,nblks_c]
    CALL this%w_wind_ic    %Init('w_wind_ic',     'double', dims=shape_3d, patch_id=patch_id)
    CALL this%pres_ic      %Init('pres_ic',       'double', dims=shape_3d, patch_id=patch_id)
    CALL this%geo_height_ic%Init('geo_height_ic', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%geopot_agl_ic%Init('geopot_agl_ic', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%z_ic         %Init('z_ic',          'double', dims=shape_3d, patch_id=patch_id)
    CALL this%dz_ic        %Init('dz_ic',         'mixed' , dims=shape_3d, patch_id=patch_id)
    CALL this%inv_dz_ic    %Init('inv_dz_ic',     'double', dims=shape_3d, patch_id=patch_id)

    ! Initialize 3D fields at full levels on edges
    shape_3d = [nproma,nlev,nblks_e]
    ! CALL this%vn_e         %Init('vn_e',          'double', dims=shape_3d, patch_id=patch_id)

    shape_2d = [nproma,nblks_c]
    CALL this%fract_land   %Init('fract_land',    'double', dims=shape_2d, patch_id=patch_id)
    CALL this%fract_ice    %Init('fract_ice',     'double', dims=shape_2d, patch_id=patch_id)

  END SUBROUTINE build_vdf_atmo_inputs

  ! Builds and initializes the vdf_atmo_diags data structure
  !
  ! This subroutine initializes all diagnostic fields generated during the
  ! vertical diffusion process. These include exchange coefficients, stability
  ! parameters, divergence, shear, velocity components at various locations,
  ! and the components of the velocity gradient tensor. Additionally, it
  ! initializes the 3x3 velocity gradient tensor structure by setting its
  ! pointers to the corresponding gradient components.
  SUBROUTINE build_vdf_atmo_diags(this, domain)

    TYPE(t_vdf_atmo_diags), INTENT(inout) :: this    !< The vdf_atmo_diags structure to be initialized
    TYPE(t_domain),               POINTER :: domain  !< Domain specification containing grid dimensions

    INTEGER :: nproma, nblks_c, nblks_e, nblks_v, nlev, nlevp1, patch_id
    INTEGER :: shape_2d(2), shape_3d(3)

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"build_vdf_atmo_diags"

    ! Extract information from domain
    nproma   = domain%nproma  ! Number of grid points per block
    nblks_c  = domain%nblks_c ! Number of blocks for cells
    nblks_e  = domain%nblks_e ! Number of blocks for edges
    nblks_v  = domain%nblks_v ! Number of blocks for vertices
    nlev     = domain%nlev    ! Number of vertical levels
    nlevp1   = nlev + 1       ! Number of half levels
    patch_id = domain%id      ! Patch ID for the domain

    ! Initialize 3D fields at full levels on cells
    shape_3d = [nproma,nlev,nblks_c]
    CALL this%ghf         %Init('ghf',          'double', dims=shape_3d, patch_id=patch_id)
    CALL this%ctgz        %Init('ctgz',         'double', dims=shape_3d, patch_id=patch_id)
    CALL this%div_c       %Init('div_c',        'double', dims=shape_3d, patch_id=patch_id)
    CALL this%theta_v     %Init('theta_v',      'double', dims=shape_3d, patch_id=patch_id)
    CALL this%pprfac      %Init('pprfac',       'double', dims=shape_3d, patch_id=patch_id)
    CALL this%km          %Init('km',           'double', dims=shape_3d, patch_id=patch_id)
    CALL this%kh          %Init('kh',           'double', dims=shape_3d, patch_id=patch_id)
    CALL this%km_c        %Init('km_c',         'double', dims=shape_3d, patch_id=patch_id)
    CALL this%heating     %Init('heating',      'double', dims=shape_3d, patch_id=patch_id)
    CALL this%dissip_ke   %Init('dissip_ke',    'double', dims=shape_3d, patch_id=patch_id)

    ! Initialize 3D fields at half levels on cells
    shape_3d = [nproma,nlevp1,nblks_c]
    CALL this%rho_ic      %Init('rho_ic',       'double', dims=shape_3d, patch_id=patch_id)
    CALL this%bruvais     %Init('bruvais',      'double', dims=shape_3d, patch_id=patch_id)
    CALL this%stab_func   %Init('stab_func',    'double', dims=shape_3d, patch_id=patch_id)
    CALL this%mech_prod   %Init('mech_prod',    'double', dims=shape_3d, patch_id=patch_id)
    CALL this%km_ic       %Init('km_ic',        'double', dims=shape_3d, patch_id=patch_id)
    CALL this%kh_ic       %Init('kh_ic',        'double', dims=shape_3d, patch_id=patch_id)
    CALL this%mix_len_sq  %Init('mix_len_sq',   'double', dims=shape_3d, patch_id=patch_id)

    ! Initialize 3D fields at full levels on edges
    shape_3d = [nproma,nlev,nblks_e]
    CALL this%vn          %Init('vn',           'double', dims=shape_3d, patch_id=patch_id)
    CALL this%shear       %Init('shear',        'double', dims=shape_3d, patch_id=patch_id)
    CALL this%div_stress  %Init('div_stress',   'double', dims=shape_3d, patch_id=patch_id)

    ! Initialize 3D fields at half levels on edges
    shape_3d = [nproma,nlevp1,nblks_e]
    CALL this%vn_ie       %Init('vn_ie',        'double', dims=shape_3d, patch_id=patch_id)
    CALL this%vt_ie       %Init('vt_ie',        'double', dims=shape_3d, patch_id=patch_id)
    CALL this%w_ie        %Init('w_ie',         'double', dims=shape_3d, patch_id=patch_id)
    CALL this%km_ie       %Init('km_ie',        'double', dims=shape_3d, patch_id=patch_id)

    ! Initialize 3D fields at full levels on vertices
    shape_3d = [nproma,nlev,nblks_v]
    CALL this%u_vert      %Init('u_vert',       'double', dims=shape_3d, patch_id=patch_id)
    CALL this%v_vert      %Init('v_vert',       'double', dims=shape_3d, patch_id=patch_id)

    ! Initialize 3D fields at half levels on vertices
    shape_3d = [nproma,nlevp1,nblks_v]
    CALL this%w_vert      %Init('w_vert',       'double', dims=shape_3d, patch_id=patch_id)
    CALL this%km_iv       %Init('km_iv',        'double', dims=shape_3d, patch_id=patch_id)

    ! Initialize 2D (vertically integrated) diagnostics
    shape_2d = [nproma,nblks_c]
    CALL this%ctgzvi      %Init('ctgzvi',       'double', dims=shape_2d, patch_id=patch_id)
    CALL this%dissip_ke_vi%Init('dissip_ke_vi', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%int_energy_vi%Init('int_energy_vi','double', dims=shape_2d, patch_id=patch_id)
    CALL this%int_energy_vi_tend%Init('int_energy_vi_tend','double', dims=shape_2d, patch_id=patch_id)

    ! Louis scaling factor and surface fluxes
    CALL this%louis_factor%Init('louis_factor', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%lhflx       %Init('lhflx',        'double', dims=shape_2d, patch_id=patch_id)
    CALL this%shflx       %Init('shflx',        'double', dims=shape_2d, patch_id=patch_id)

    ! Initialize velocity gradient tensor components on edges
    shape_3d = [nproma,nlev,nblks_e]
    CALL this%dvn_dn      %Init('dvn_dn',       'double', dims=shape_3d, patch_id=patch_id)
    CALL this%dvt_dn      %Init('dvt_dn',       'double', dims=shape_3d, patch_id=patch_id)
    CALL this%dw_dn       %Init('dw_dn',        'double', dims=shape_3d, patch_id=patch_id)
    CALL this%dvn_dt      %Init('dvn_dt',       'double', dims=shape_3d, patch_id=patch_id)
    CALL this%dvt_dt      %Init('dvt_dt',       'double', dims=shape_3d, patch_id=patch_id)
    CALL this%dw_dt       %Init('dw_dt',        'double', dims=shape_3d, patch_id=patch_id)
    CALL this%dvn_dz      %Init('dvn_dz',       'double', dims=shape_3d, patch_id=patch_id)
    CALL this%dvt_dz      %Init('dvt_dz',       'double', dims=shape_3d, patch_id=patch_id)
    CALL this%dw_dz       %Init('dw_dz',        'double', dims=shape_3d, patch_id=patch_id)

    ! Initialize velocity gradient tensor on edges
    IF (this%is_initialized) THEN
      this%vel_grad_e(1,1)%ptr => this%dvn_dn%Get_ptr_r3d()
      this%vel_grad_e(1,2)%ptr => this%dvn_dt%Get_ptr_r3d()
      this%vel_grad_e(1,3)%ptr => this%dvn_dz%Get_ptr_r3d()
      this%vel_grad_e(2,1)%ptr => this%dvt_dn%Get_ptr_r3d()
      this%vel_grad_e(2,2)%ptr => this%dvt_dt%Get_ptr_r3d()
      this%vel_grad_e(2,3)%ptr => this%dvt_dz%Get_ptr_r3d()
      this%vel_grad_e(3,1)%ptr => this%dw_dn%Get_ptr_r3d()
      this%vel_grad_e(3,2)%ptr => this%dw_dt%Get_ptr_r3d()
      this%vel_grad_e(3,3)%ptr => this%dw_dz%Get_ptr_r3d()
      !$ACC ENTER DATA COPYIN(this)
      __acc_attach(this%vel_grad_e(1,1)%ptr)
      __acc_attach(this%vel_grad_e(1,2)%ptr)
      __acc_attach(this%vel_grad_e(1,3)%ptr)
      __acc_attach(this%vel_grad_e(2,1)%ptr)
      __acc_attach(this%vel_grad_e(2,2)%ptr)
      __acc_attach(this%vel_grad_e(2,3)%ptr)
      __acc_attach(this%vel_grad_e(3,1)%ptr)
      __acc_attach(this%vel_grad_e(3,2)%ptr)
      __acc_attach(this%vel_grad_e(3,3)%ptr)
    END IF

    this%is_initialized = .TRUE.

  END SUBROUTINE build_vdf_atmo_diags

END MODULE mo_vdf_atmo_memory
