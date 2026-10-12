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

! Memory management for the vertical diffusion surface component
!
! This module provides the memory management structures and routines for
! the vertical diffusion component of the surface model. It defines
! the input data structures required for the surface calculations
! and provides routines for their initialization.
!

MODULE mo_vdf_sfc_memory

  USE mo_kind,            ONLY: wp, vp, i1, i4
  USE mo_tmx_var,         ONLY: t_tmx_var
  USE mo_tmx_field_class, ONLY: t_domain

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_vdf_sfc_inputs, build_vdf_sfc_inputs
  PUBLIC :: t_vdf_sfc_diags, build_vdf_sfc_diags
  PUBLIC :: t_vdf_sfc_config, build_vdf_sfc_config

  ! Data structure containing configuration parameters for surface VDF
  !
  ! This type contains all the configuration parameters needed for the
  ! vertical diffusion calculations at the surface.
  TYPE t_vdf_sfc_config
    ! Physics constants and configuration parameters
    TYPE(t_tmx_var) :: &
      & dtime         , & !< Time step [s]
      & cpd           , & !< Specific heat capacity at constant pressure for dry air [J/(kg K)]
      & cvd           , & !< Specific heat capacity at constant volume for dry air [J/(kg K)]
      & cvv           , & !< Specific heat capacity at constant volume for water vapor [J/(kg K)]
      & min_sfc_wind  , & !< Minimum surface wind speed [m/s]
      & wind_gustiness, & !< Wind gustiness factor [-]
      & min_rough     , & !< Minimum roughness length [m]
      & rough_m_oce   , & !< Ocean roughness length for momentum [m]
      & rough_m_ice   , & !< Ice roughness length for momentum [m]
      & fsl               !< Weight for interpolation to surface layer mid level [-]
    TYPE(t_tmx_var) :: &
      & l_co2         , & !< Use CO2 in calculations (logical) [-]
      & nice_thickness_classes !< Number of sea ice thickness classes [-]
  END TYPE t_vdf_sfc_config

  ! Data structure containing all input fields required for surface diffusion
  !
  ! This type contains all the meteorological fields needed as inputs for
  ! the surface diffusion calculations.
  TYPE t_vdf_sfc_inputs

    ! Atmospheric values at lowest level
    TYPE(t_tmx_var) :: &
      & ta          , & !< Air temperature at lowest level [K]
      & tv          , & !< Virtual temperature at lowest level [K]
      & ua          , & !< Zonal wind component at lowest level [m/s]
      & va          , & !< Meridional wind component at lowest level [m/s]
      & qa          , & !< Specific humidity at lowest level [kg/kg]
      & rho_atm     , & !< Air density at lowest level [kg/m3]
      & pa          , & !< Pressure at lowest level [Pa]
      & psfc        , & !< Surface pressure [Pa]
      & zf          , & !< Height of lowest atmospheric full level [m]
      & zh              !< Height of surface [m]

    ! Surface precipitation fluxes
    TYPE(t_tmx_var) :: &
      & rsfl        , & !< Surface rain flux, large-scale [kg/(m2 s)]
      & ssfl            !< Surface snow flux, large-scale [kg/(m2 s)]

    ! Radiation fluxes
    TYPE(t_tmx_var) :: &
      & rlds        , & !< Surface downward longwave radiation [W/m2]
      & rsds        , & !< Surface downward shortwave radiation [W/m2]
      & rvds_dir    , & !< Direct visible radiation [W/m2]
      & rnds_dir    , & !< Direct near-IR radiation [W/m2]
      & rpds_dir    , & !< Direct PAR radiation [W/m2]
      & rvds_dif    , & !< Diffuse visible radiation [W/m2]
      & rnds_dif    , & !< Diffuse near-IR radiation [W/m2]
      & rpds_dif    , & !< Diffuse PAR radiation [W/m2]
      & emissivity  , & !< Longwave surface emissivity [-]
      & cosmu0          !< Cosine of zenith angle [-]

    ! Other atmospheric and surface properties
    TYPE(t_tmx_var) :: &
      & co2          , & !< Atmospheric CO2 concentration [ppmv]
      & co2flx_ant   , & !< CO2 flux from anthropogenic sources [kg/(m2 s)]
      & ocean_u,       & !< U component of ocean current [m/s]
      & ocean_v,       & !< V component of ocean current [m/s]
      & ice_u,         & !< U component of sea ice velocity [m/s]
      & ice_v,         & !< V component of sea ice velocity [m/s]
      & ice_thickness    !< Thickness of sea ice [m]

    ! Tile-based fields
    TYPE(t_tmx_var) :: &
      & tsfc_tile   , & !< Surface temperature per tile [K]
      & fract_tile      !< Fractional area per tile [-]

    ! Reference height in surface layer
    TYPE(t_tmx_var) :: dz
  END TYPE t_vdf_sfc_inputs

  ! Data structure containing all diagnostic fields for surface diffusion
  !
  ! This type contains all the diagnostic fields needed to analyze and understand
  ! the surface diffusion process, including surface fluxes and properties.
  TYPE t_vdf_sfc_diags
    LOGICAL :: is_initialized = .FALSE.

    ! Basic atmospheric and surface properties
    TYPE(t_tmx_var) :: &
      & wind_rel_tile, & !< Atmospheric wind speed at lowest level rel. to surface [m/s]
      & theta_atm   , & !< Atmospheric potential temperature at lowest level [K]
      & thetav_atm      !< Atmospheric virtual potential temperature at lowest level [K]

    ! Surface roughness
    TYPE(t_tmx_var) :: &
      & rough_h     , & !< Grid-mean roughness length for heat [m]
      & rough_m     , & !< Grid-mean roughness length for momentum [m]
      & rough_h_tile, & !< Roughness length for heat per tile [m]
      & rough_m_tile    !< Roughness length for momentum per tile [m]

    ! Surface temperature and humidity
    TYPE(t_tmx_var) :: &
      & tsfc        , & !< Grid-mean surface temperature [K]
      & tsfc_rad    , & !< Grid-mean radiative surface temperature [K]
      & qsat_tile       !< Saturation specific humidity per tile [kg/kg]

    ! Surface exchange coefficients
    TYPE(t_tmx_var) :: &
      & kh                  , & !< Grid-mean exchange coefficient for heat [-]
      & km                  , & !< Grid-mean exchange coefficient for momentum [-]
      & kh_tile             , & !< Exchange coefficient for heat per tile [-]
      & km_tile             , & !< Exchange coefficient for momentum per tile [-]
      & interp_fac_10m_tile , & !< Interpolation coefficient in 10m height per tile [-]
      & interp_fac_2m_tile  , & !< Interpolation coefficient in 2m height per tile [-]
      & interp_fac_tsfc_tile    !< Product of interpolation coefficient and surface temperature per tile [K]

    ! Surface fluxes - grid mean
    TYPE(t_tmx_var) :: &
      & evapotrans  , & !< Grid-mean surface evapotranspiration [W/m2]
      & lhfl        , & !< Grid-mean latent heat flux [W/m2]
      & shfl        , & !< Grid-mean sensible heat flux [W/m2]
      & ustress     , & !< Grid-mean zonal wind stress [N/m2]
      & vstress     , & !< Grid-mean meridional wind stress [N/m2]
      & ufts        , & !< Energy flux at surface from thermal exchange [W/m2]
      & ufvs        , & !< Energy flux at surface from vapor exchange [W/m2]
      & lwfl_up     , & !< Grid-mean upward longwave radiation [W/m2]
      & swfl_up     , & !< Grid-mean upward shortwave radiation [W/m2]
      & co2flx_nat  , & !< CO2 flux from natural sources [kg/(m2 s)]
      & co2flx          !< Total CO2 flux at surface [kg/(m2 s)]

    ! Surface fluxes - per tile
    TYPE(t_tmx_var) :: &
      & evapotrans_tile, & !< Surface evapotranspiration per tile [W/m2]
      & lhfl_tile      , & !< Latent heat flux per tile [W/m2]
      & shfl_tile      , & !< Sensible heat flux per tile [W/m2]
      & ustress_tile   , & !< Zonal wind stress per tile [N/m2]
      & vstress_tile   , & !< Meridional wind stress per tile [N/m2]
      & lwfl_net_tile  , & !< Net longwave flux per tile [W/m2]
      & swfl_net_tile  , & !< Net shortwave flux per tile [W/m2]
      & co2flx_nat_tile    !< CO2 flux from natural sources per tile [kg/(m2 s)]

    ! Surface properties per tile
    TYPE(t_tmx_var) :: &
      & theta_tile  , & !< Potential temperature per tile [K]
      & thetav_tile , & !< Virtual potential temperature per tile [K]
      & rho_tile        !< Air density at surface per tile [kg/m3]

    ! Stability parameters
    TYPE(t_tmx_var) :: &
      & moist_rich_tile !< Moist Richardson number per tile [-]

    ! Albedo - grid mean
    TYPE(t_tmx_var) :: &
      & albvisdir   , & !< Grid-mean albedo for visible, direct [-]
      & albvisdif   , & !< Grid-mean albedo for visible, diffuse [-]
      & albnirdir   , & !< Grid-mean albedo for near-IR, direct [-]
      & albnirdif   , & !< Grid-mean albedo for near-IR, diffuse [-]
      & albedo          !< Grid-mean broadband albedo [-]

    ! Albedo - per tile
    TYPE(t_tmx_var) :: &
      & albvisdir_tile, & !< Albedo for visible, direct per tile [-]
      & albvisdif_tile, & !< Albedo for visible, diffuse per tile [-]
      & albnirdir_tile, & !< Albedo for near-IR, direct per tile [-]
      & albnirdif_tile, & !< Albedo for near-IR, diffuse per tile [-]
      & albedo_tile       !< Broadband albedo per tile [-]

    ! Sea ice specific diagnostics
    TYPE(t_tmx_var) :: &
      & q_ice_top     , & !< Energy flux available for surface melting of sea ice [W/m2]
      & q_ice_bot     , & !< Energy flux at ice-ocean interface [W/m2]
      & snow_thickness    !< Thickness of snow on sea ice [m]

    ! Land specific diagnostics
    TYPE(t_tmx_var) :: &
      & q_snocpymlt_lnd !< Heating used to melt snow on canopy [W/m2]

    ! 10m and 2m diagnostics
    TYPE(t_tmx_var) :: &
      & t2m         , & !< 2m temperature [K]
      & t2m_tile    , & !< 2m temperature per tile [K]
      & hus2m       , & !< 2m specific humidity [kg/kg]
      & hus2m_tile  , & !< 2m specific humidity per tile [kg/kg]
      & dew2m       , & !< 2m dewpoint temperature [K]
      & dew2m_tile  , & !< 2m dewpoint temperature per tile [K]
      & wind10m     , & !< 10m wind speed [m/s]
      & u10m        , & !< 10m zonal wind [m/s]
      & v10m        , & !< 10m meridional wind [m/s]
      & wind10m_tile, & !< 10m wind speed per tile [m/s]
      & u10m_tile   , & !< 10m zonal wind per tile [m/s]
      & v10m_tile       !< 10m meridional wind per tile [m/s]

    ! Tile index handling
    TYPE(t_tmx_var) :: &
      & nvalid      , & !< Number of valid points per block and tile
      & indices         !< Indices of valid points on chunk per block and tile
  END TYPE t_vdf_sfc_diags

  ! Module name for error messages and debugging
  CHARACTER(LEN=*), PARAMETER :: modname="mo_vdf_sfc_memory"

  CONTAINS

  ! Builds and initializes the vdf_sfc_config data structure
  !
  ! This subroutine initializes all fields in the t_vdf_sfc_config structure
  ! with appropriate attributes.
  !
  ! this    The vdf_sfc_config structure to be initialized
  ! domain  Domain specification containing grid dimensions
  SUBROUTINE build_vdf_sfc_config(this, domain)

    TYPE(t_vdf_sfc_config), INTENT(inout) :: this
    TYPE(t_domain),              POINTER  :: domain

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"build_vdf_sfc_config"

    INTEGER :: patch_id
    INTEGER :: shape_0d(0) ! Zero-szied array for scalar configuration parameters

    ! Extract information from domain
    patch_id = domain%id  ! Patch ID for the domain

    ! Initialize scalar configuration parameters (0D)
    CALL this%dtime%Init('time step', 'double', dims=shape_0d, patch_id=patch_id)
    CALL this%cpd%Init('cpd', 'double', dims=shape_0d, patch_id=patch_id)
    CALL this%cvd%Init('cvd', 'double', dims=shape_0d, patch_id=patch_id)
    CALL this%cvv%Init('cvv', 'double', dims=shape_0d, patch_id=patch_id)
    CALL this%min_sfc_wind%Init('minimum surface wind speed', 'double', dims=shape_0d, patch_id=patch_id)
    CALL this%wind_gustiness%Init('wind gustiness factor', 'double', dims=shape_0d, patch_id=patch_id)
    CALL this%min_rough%Init('minimal roughness length', 'double', dims=shape_0d, patch_id=patch_id)
    CALL this%rough_m_oce%Init('ocean roughness length', 'double', dims=shape_0d, patch_id=patch_id)
    CALL this%rough_m_ice%Init('ice roughness length', 'double', dims=shape_0d, patch_id=patch_id)
    CALL this%fsl%Init('weight for interpolation to surface_layer mid level', 'double', dims=shape_0d, patch_id=patch_id)
    CALL this%nice_thickness_classes%Init('number of sea ice thickness classes', 'integer', dims=shape_0d, patch_id=patch_id)
    CALL this%l_co2%Init('use CO2', 'logical', dims=shape_0d, patch_id=patch_id)

  END SUBROUTINE build_vdf_sfc_config

  ! Builds and initializes the vdf_sfc_inputs data structure
  !
  ! This subroutine initializes all fields in the t_vdf_sfc_inputs structure
  ! with appropriate dimensions and attributes based on the domain specifications.
  !
  ! this    The vdf_sfc_inputs structure to be initialized
  ! domain  Domain specification containing grid dimensions
  SUBROUTINE build_vdf_sfc_inputs(this, domain)

    TYPE(t_vdf_sfc_inputs), INTENT(inout) :: this
    TYPE(t_domain),               POINTER :: domain

    INTEGER :: nproma, nblks_c, ntiles, patch_id
    INTEGER :: shape_2d(2), shape_3d(3)

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"build_vdf_sfc_inputs"

    ! Extract information from domain
    nproma   = domain%nproma   ! Number of grid points per block
    nblks_c  = domain%nblks_c  ! Number of blocks for cells
    ntiles   = domain%ntiles   ! Number of surface tiles
    patch_id = domain%id       ! Patch ID for the domain

    ! Initialize 2D fields (horizontal domain)
    shape_2d = [nproma, nblks_c]

    ! Atmospheric values at lowest level
    CALL this%ta%Init('atm temperature', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%tv%Init('atm virtual temperature', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%ua%Init('atm zonal wind', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%va%Init('atm meridional wind', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%qa%Init('atm total water', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%rho_atm%Init('atm density', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%pa%Init('atm full level pressure', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%psfc%Init('surface pressure', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%zf%Init('atm geometric height full', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%zh%Init('sfc geometric height half', 'double', dims=shape_2d, patch_id=patch_id)

    ! Surface precipitation fluxes
    CALL this%rsfl%Init('surface rain flux, large-scale', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%ssfl%Init('surface snow flux, large-scale', 'double', dims=shape_2d, patch_id=patch_id)

    ! Radiation fluxes
    CALL this%rlds%Init('surface downward longwave radiation', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%rsds%Init('surface downward shortwave radiation', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%rvds_dir%Init('all-sky surface downward direct visible radiation', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%rnds_dir%Init('all-sky surface downward direct near-IR radiation', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%rpds_dir%Init('all-sky surface downward direct PAR radiation', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%rvds_dif%Init('all-sky surface downward diffuse visible radiation', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%rnds_dif%Init('all-sky surface downward diffuse near-IR radiation', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%rpds_dif%Init('all-sky surface downward diffuse PAR radiation', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%emissivity%Init('longwave surface emissivity', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%cosmu0%Init('cosine of zenith angle', 'double', dims=shape_2d, patch_id=patch_id)

    ! CO2 fluxes
    CALL this%co2flx_ant%Init('CO2 flux from anthropogenic sources', 'double', dims=shape_2d, patch_id=patch_id)

    ! Other atmospheric and surface properties
    CALL this%co2%Init('atm CO2 concentration', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%ocean_u%Init('u-component of ocean current', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%ocean_v%Init('v-component of ocean current', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%ice_u%Init('u-component of sea ice velocity', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%ice_v%Init('v-component of sea ice velocity', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%ice_thickness%Init('thickness of sea ice', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%dz%Init('reference height in surface layer times 2', 'mixed', dims=shape_2d, patch_id=patch_id)

    ! Tile-based fields (3D: horizontal domain + tiles)
    shape_3d = [nproma, nblks_c, ntiles]
    CALL this%tsfc_tile%Init('surface temperature, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%fract_tile%Init('grid box fraction of tiles', 'double', dims=shape_3d, patch_id=patch_id)

  END SUBROUTINE build_vdf_sfc_inputs

  ! Builds and initializes the vdf_sfc_diags data structure
  !
  ! This subroutine initializes all diagnostic fields in the t_vdf_sfc_diags structure
  ! with appropriate dimensions and attributes based on the domain specifications.
  !
  ! this    The vdf_sfc_diags structure to be initialized
  ! domain  Domain specification containing grid dimensions
  SUBROUTINE build_vdf_sfc_diags(this, domain)

    TYPE(t_vdf_sfc_diags), INTENT(inout) :: this
    TYPE(t_domain),              POINTER :: domain

    INTEGER :: nproma, nblks_c, ntiles, patch_id
    INTEGER :: shape_2d(2), shape_3d(3), shape_nvalid(2), shape_indices(3)

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"build_vdf_sfc_diags"

    ! Extract information from domain
    nproma   = domain%nproma   ! Number of grid points per block
    nblks_c  = domain%nblks_c  ! Number of blocks for cells
    ntiles   = domain%ntiles   ! Number of surface tiles
    patch_id = domain%id       ! Patch ID for the domain

    ! Initialize 2D fields (horizontal domain)
    shape_2d = [nproma, nblks_c]

    ! Basic atmospheric and surface properties
    CALL this%theta_atm%Init('atm potential temperature', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%thetav_atm%Init('atm virtual potential temperature', 'double', dims=shape_2d, patch_id=patch_id)

    ! Surface roughness - grid mean
    CALL this%rough_h%Init('roughness length heat', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%rough_m%Init('roughness length momentum', 'double', dims=shape_2d, patch_id=patch_id)

    ! Surface temperature and surface fluxes - grid mean
    CALL this%tsfc%Init('sfc temperature', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%tsfc_rad%Init('sfc radiative temperature', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%evapotrans%Init('sfc evapotranspiration', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%lhfl%Init('sfc latent heat flux', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%shfl%Init('sfc sensible heat flux', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%ustress%Init('sfc zonal wind stress', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%vstress%Init('sfc mer. wind stress', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%ufts%Init('energy flux at surface from thermal exchange', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%ufvs%Init('energy flux at surface from vapor exchange', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%lwfl_up%Init('sfc longwave upward flux', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%swfl_up%Init('sfc shortwave upward flux', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%co2flx_nat%Init('CO2 flux from natural sources', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%co2flx%Init('total CO2 flux at surface', 'double', dims=shape_2d, patch_id=patch_id)

    ! Surface exchange coefficients - grid mean
    CALL this%kh%Init('exchange coefficient for scalar', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%km%Init('exchange coefficient for momentum', 'double', dims=shape_2d, patch_id=patch_id)

    ! Albedo - grid mean
    CALL this%albvisdir%Init('albedo VIS direct', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%albvisdif%Init('albedo VIS diffuse', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%albnirdir%Init('albedo NIR direct', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%albnirdif%Init('albedo NIR diffuse', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%albedo%Init('albedo', 'double', dims=shape_2d, patch_id=patch_id)

    ! Sea ice specific diagnostics
    CALL this%q_ice_top%Init('energy flux available for surface melting of sea ice', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%q_ice_bot%Init('energy flux at ice-ocean interface', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%snow_thickness%Init('thickness of snow on sea ice', 'double', dims=shape_2d, patch_id=patch_id)

    ! Land specific diagnostics
    CALL this%q_snocpymlt_lnd%Init('heating used to melt snow on canopy', 'double', dims=shape_2d, patch_id=patch_id)

    ! 10m and 2m diagnostics - grid mean
    CALL this%t2m%Init('2m temperature', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%hus2m%Init('2m specific humidity', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%dew2m%Init('2m dewpoint temperature', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%wind10m%Init('10m wind speed', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%u10m%Init('10m zonal wind', 'double', dims=shape_2d, patch_id=patch_id)
    CALL this%v10m%Init('10m meridional wind', 'double', dims=shape_2d, patch_id=patch_id)

    ! Tile-based fields (3D: horizontal domain + tiles)
    shape_3d = [nproma, nblks_c, ntiles]

    ! Surface properties per tile
    CALL this%wind_rel_tile%Init('atm rel. wind speed, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%rough_h_tile%Init('roughness length heat, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%rough_m_tile%Init('roughness length momentum, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%qsat_tile%Init('sfc saturation specific humidity, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%theta_tile%Init('sfc potential temperature, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%thetav_tile%Init('sfc virtual potential temperature, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%rho_tile%Init('sfc density, tile', 'double', dims=shape_3d, patch_id=patch_id)

    ! Surface exchange coefficients per tile
    CALL this%kh_tile%Init('exchange coefficient for scalar, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%km_tile%Init('exchange coefficient for momentum, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%interp_fac_2m_tile%Init('interpolation coef. at 2m height, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%interp_fac_10m_tile%Init('interpolation coef. at 10m heigh, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%interp_fac_tsfc_tile%Init('product of interp. coef. and sfc temp., tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%moist_rich_tile%Init('moist richardson number, tile', 'double', dims=shape_3d, patch_id=patch_id)

    ! Surface fluxes per tile
    CALL this%evapotrans_tile%Init('sfc evapotranspiration, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%lhfl_tile%Init('sfc latent heat flux, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%shfl_tile%Init('sfc sensible heat flux, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%ustress_tile%Init('sfc zonal wind stress, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%vstress_tile%Init('sfc mer. wind stress, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%lwfl_net_tile%Init('sfc longwave net flux, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%swfl_net_tile%Init('sfc shortwave net flux, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%co2flx_nat_tile%Init('CO2 flux from natural sources, tile', 'double', dims=shape_3d, patch_id=patch_id)

    ! Albedo per tile
    CALL this%albvisdir_tile%Init('albedo VIS direct, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%albvisdif_tile%Init('albedo VIS diffuse, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%albnirdir_tile%Init('albedo NIR direct, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%albnirdif_tile%Init('albedo NIR diffuse, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%albedo_tile%Init('albedo, tile', 'double', dims=shape_3d, patch_id=patch_id)

    ! 10m and 2m diagnostics per tile
    CALL this%t2m_tile%Init('2m temperature, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%hus2m_tile%Init('2m specific humidity, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%dew2m_tile%Init('2m dewpoint temperature, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%wind10m_tile%Init('10m wind speed, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%u10m_tile%Init('10m zonal wind, tile', 'double', dims=shape_3d, patch_id=patch_id)
    CALL this%v10m_tile%Init('10m meridional wind, tile', 'double', dims=shape_3d, patch_id=patch_id)

    ! Tile index handling
    shape_nvalid = [nblks_c, ntiles]
    shape_indices = [nproma, nblks_c, ntiles]
    CALL this%nvalid%Init('number of valid points per block, tile', 'integer', dims=shape_nvalid, patch_id=patch_id)
    CALL this%indices%Init('indices of valid points on chunk per block, tile', 'integer', dims=shape_indices, patch_id=patch_id)

    this%is_initialized = .TRUE.

  END SUBROUTINE build_vdf_sfc_diags

END MODULE mo_vdf_sfc_memory
