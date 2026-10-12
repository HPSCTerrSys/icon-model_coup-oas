# Release notes for icon-2026.04

## ICON-Atmo

### Dynamics

- Swap dimensions of coefficient fields for igradp_method from (sparse, nproma, ...) to (nproma, sparse, ...) for consistency with icon4py.
- The velocity tendencies are now accessed with same pattern as prognostics (using nnow and nnew). This assumes `itime_scheme>=4` for the dycore.
- Fixed computation of total integrals diagnostic with output namelist option "totint"
- Performance optimization on GPU for the computation of runtime maximum wind speed diagnostic

#### Assimilation Interfaces

- Update of `dace_icon` and changes in data assimilation interfaces.
- EMVORADO: Replace NetCDF77 procedures by NetCDF90 procedures.

### AES Physics

- Turbulent mixing (TMX)
  - Use Memory Manager library for data management
  - Use tile-specific relative wind speed for land/ocean/sea ice tiles in computations of surface roughness, exchange coefficients and fluxes
  - Added namelist options to not use Louis stability formula over land and/or sea ice
  - Fixes
    - Fixed GPU async queues in surface module to prevent possible race conditions
    - Implemented several OpenACC-related fixes to prevent model crashes in some circumstances
    - Fixed loop indexing in surface module to prevent acting on cells outside of valid domain
    - Fixed the vertical interpolation in the computation of near-surface diagnostics effecting 2m temperature and 10m wind

- Microphysics
  - Added C++/Kokkos implementation of AES microphysics (used with `--enable-ragnarok`)  
    **Development of the Fortran implementation of AES microphysics is frozen and it may be removed in a future release**
  - Added a `cloud_mig` tuning parameter `cia` to change the ice sticking efficiency
  - Changed the default maximum height up to which microphysics is computed from 33000m to 22500m
    (same as default height up to which transport of water tracers is computed)
  - Fixed double definition of cloud droplet number concentration
  - Corrections for rain microphysics (evaporation and accretion)

- RTE-RRTMGP radiation
  - Simplified the radiation interface by using only one unique criterion for the presence of clouds
  - Restructured the radiation interface for better readability

### NWP Physics

- Implementation of gravity waves parametrization MS-GWaM
- Revision of diagnostics for ceiling and visibility to provide a better agreement with SYNOP measurements.
- Tuning options to improve 10m-winds in mountainous regions.
- Modified cloud fraction for two-moment cloud ice microphysics.
- Change default setting for `init_latbc_from_fg` to `TRUE`, to prevent mismatches between initial conditions and the initial set of lateral boundary data in LAM setups.
- Synchronization of ecRad with ECMWF code.
- Provide interactive CO2 to ecRad radiation.
- Finish when using Intel+OpenMP+ecRad: Using OpenMP for ecRad with Intel compiler leads to non-reproducible results due to various race conditions.
- Coupling CAMS aerosols with Segal and Khain cloud droplets activation scheme.
- Upgrades to Single Column and LES Mode.
- Added diagnostic output fields: accumulated freezing rain, optional turbulence diagnostics, hail diagnostics for the 2-moment microphysics.
- Remove init mode `MODE_IAU_OLD`.
- Deprecate RRTM radiation and associated ozone option (`irad_o3 = 6`).
- Fix uninitialized variables in ecRad interface.
- Fix OpenMP directive in Kinne aerosol preparation for ecRad.

#### Upper atmosphere physics

- Further encapsulation and consolidation of the code base

## ICON-Ocean

- Implementation of symmetric incremental analysis update (IAU) in the ocean
- Implement the core features for limited area ocean setups
- Ocean ERA5 forcing via python provider script and YAC coupling
- Solver option lhs_direct=.true. now also available for zstar
- Optimizations
  - Port solver global sums on GPU
  - OpenMP improvements
- Diagnostics
  - Ocean vorticity on cells for output with hiopy
  - Surface-level-type for Sea Ice Vars changed to ice_class
  - Seaice principal stress diagnostic
- Fixes
  - Fix write_initial_state (variable assignments and transfer to GPU)
  - Include jstep_shift to fix timing of ocean output in case of enabled timeshift
  - Bugfix for HAMOCC binary identical restarts
  - Correct unit of amoc26n diagnostic
- Add shallow water as supported configuration

## ICON-Waves

- Support NetCDF file input for model initialization
- Restructure analytic wave spectrum initialization
- Revise wave propagation testcase
- Improve memory layout of 4D wave fields
- Support grib2 encoding of wave energy spectra
- Add coupling to the ocean component
- Code improvements regarding active and passive diagnostics
- Add diagnostics: wave steepness, kurtosis, Benjamin-Feir index, peak frequency and direction
- Fixes:
  - Boundary condition for wave group velocity at ocean-land boundary
  - Units and descriptions for wave-specific fields

## Soil and Surface

### Climate: ICON-Land

- Enabled running ICON-Land standalone for one or more single sites
- Enabled JSBACH simulations with natural and anthropogenic land cover change
- Refactored the calculation of global diagnostics and ported to GPU
- QUINCY
  - Improved vegetation Nitrogen uptake in high latitudes
  - Added interactive nitrogen and phosphorus flags in soil biogeochemistry
  - Fixed the effects of ice and supercooled water on vegetation growth and nutrient leaching
  - Fixed a bug in the vegetation litter decomposition in the carbon-only model
  - Fixed an issue in restart file comparisons caused by a purely diagnostic variable only used with optimality code
  - Fixed a bug in the soil nutrient vertical transport
  - Fixed multiple issues in soil biogeochemistry Langmuir kinetics
  - Refactored the execution control of the slow pool spinup accelerator in the soil biogeochemistry
  - Removed code parts not required in ICON-Land
  - Removed deprecated code from the SPP1685 project
- Hydrology
  - Fix concerning the effect of supercooled water on plant water stress
  - Removed checks needed with deprecated initial files
  - Simplified the reading of soil boundary condition data
- Code cleaning and technical fixes
  - Cleaned up the use of string functions from `fortran-support` library
  - Fixed freeing of some leftover memory when reading initial data
  - Fixed error message in case of invalid variable names defined in
    the output namelist
- Pre- and postprocessing scripts
  - New script to process ERA5 data for ICON-Land standalone forcing
  - Updated the script to equilibrate JSBACH humus carbon pools

### NWP: TERRA

- Fix to avoid soil ice content exceeding total soil water content.

## ICON-ART

- GPU optimization of OEM and VPRM
- Introduction of MieAI
- enable overwriting of GFAS data at selected time intervals
- Emissions of natural aerosols and trace gases from fully coupled and interactive ocean and land modules:
  - Dust emission parameterization with JSBACH and vdiff coupling
  - Dry deposition for aerosols with JSBACH and VDIFF
  - JSBACH PFTs mapped to MEGAN
- Updated submodule CloudJ from version 7.4 to version 8.0.2
- Daily varying solar irradiances / photon fluxes for CloudJ
- Rename diagnostic variables of MECCA and CLOUDJ; completed/corrected photolysis indices
- DMS online emissions for MECCA tracers
- Adjustment for allowing more than one time interval for latbc
- Reintroduce calculation of effective radii for cloud-j to fixe an issue where droplet radii were evaluated as zero, effectively ignoring clouds
- Added TLEV_UPDATE_ADVECTION group to all ART tracers to allow ComIn to properly reference them when forwarded to plugins

## Externals

- Added Memory Manager library 1.0.0 as external in externals/memman
- Added Gravity waves parametrization MS-GWaM as external in externals/msgwam
- Switch to fortran-support 2.2.1
- Switch to iconmath 1.3.0
- Switch to probtest v3.1
- Switch to YAC v3.14.0_p1
- Switch to CDI v2.6.0
- Switch to ComIn 0.5
- Switch to YAXT 0.11.5
- Switch to mkexp 1.5.0

## Documentation

- Rework of the documentation frontpage
- Include contribution, reviewing, merging and documentation guidelines
- Convert the namelist overview PDF to markdown and include it in the documentation
- Add section to announce deprecated features
- Various updates to the Getting Started and Documentation parts

## Infrastructure

- Support file input for RBF coefficients
- Fix timeshift crosscheck for single precision runs
- Add single precision tests for buildbot-CI
- Refactor output coupling to allow masking and half-levels
- Enable coupling between atmosphere on GPU and ocean on CPU
- Coupling of icon-waves and ocean
- Avoid read-in of deprecated `output_jfile` parameter from restart for model initialization
- Improve distributed patch generation from grid files
- Added support for new Memory Manager library in ICON variable lists infrastructure
  - Logical switch `lmemman` in run config, currently only set to true if `--enable-aes` and TMX turbulence in AES physics is used
  - Added `lmemman` to `add_var` parameters in order to register and allocate variables in Memory Manager
    (no effect on using such variables in traditional way)
- Fix computation of `lon_dim/lat_dim` for regular output grids. In rare cases one line of lon/lat was missing.

### Scripting and testing

- mkexp
  - Simplify CI using mkexp 1.5.0 extended command line features
  - Update documentation/examples (input file, jsbach offline config)
  - Update model coupling (fix restart, output coupling integration)
  - Fix config handling (name handling "R02B10", lists in job environment)
  - Improve scripting (run-time helpers, name handling "R02B10", progress indicator, non-SLURM systems, log timestamps, Python override at run-time)
  - Extend/update platform support (host settings, OpenMP support, vector blocking)
  - Extend/update supported configurations (AES AMIP: fix defaults, add R2B7, output config; OES: shallow water setup)
  - hiopy
    - Rewrite hiopy integration
    - Handle unimplemented time intervals and calendars
    - Fix interaction with mkexp Python settings
    - Added hiopy configuration to store simulated data directly to S3
    - Extended hiopy scripting to allow creating Zarr 2 datasets
- make_runscripts and buildbot
  - Added parameter changes for the Sapphire-2.0 configuration and basic settings and links for more resolutions to the AES AMIP template
  - Removed unused and deprecated hiopy content from the templates
  - Updated several run scripts: use recent land data to allow running the current ICON-Land version
  - Increased stacksize for default targets
  - Implemented experiment template for Jablonowski Williamson test case `jabw` and added it to single precision builders
  - Cleanup of experiment lists for single precision and related bugfixes
  - Avoid pre-fetching for ocean Williamson test
  - Deleted outdated templates
- Test cases
  - Fixed the RCE test case by making explicit which surface temperatures it can be run with and by adding the correct initial profiles for each
  - Added constant u-velocity from namelist to some ocean torus tests

### Building

- BLAS is made an optional dependency
- Report to the relevant `config.log` when calling `config.status` of the bundled packages
- Fix detection of an external ComIn
- Configure YAC and YAXT in parallel when possible
- Fix usage of configure option `--disable-async-io-rma`
- Fetch DACE from the private repository at build time when `--enable-dace`
- Fix inconsistent usage of the `GRIBAPI` macro
- Added option `--enable-lvector` to control optimizations for vector engines in the ocean component of ICON
- Refactored ICON-Land model version and release handling following ICON OSR naming scheme
- Add mkexp setup config to `make install` target
- Make the generic GCC wrapper more flexible
- Allow for CMake definitions as configure script command-line arguments
- Optionally enable tests of CMake-based bundled packages
- Disentangle dependency of NWP physics on LES configuration on GPUs (make `--disable-les` work with `--enable-gpu`)

### Miscellaneous

- Replace frettify with Codee for Fortran formatting
- Ragnarok:
  - Add contributing guidelines
  - Refactoring and cleanup
- Fix a CPU/GPU divergence in snow fraction
- Implemented workarounds to support Cray compiler 20 for AMD GPUs
- Support filter-based NetCDF-4 compression for output
- Introduction of in-place variant for `mpi_allreduce`
- Clean up `USE` declarations: removed unnecessary includes and enforced the usage of the ONLY keyword
- Reorganize grid-related modules
- Improve error handling in case of memory tracing
- Fix mtime-related memory leaks
- Replace deprecated Intel SIMD directives with OMP SIMD


# Release notes for icon-2025.10-2

## ICON-Atmo

- Improved consistency of surface roughness, drag and momentum flux over ocean/sea ice in VDIFF
- Added viscous term in computation of surface roughness over ocean in VDIFF
- Implement a new version of the Mironov scheme that accounts for the thermal effects of snow on seaice

### NWP Physics

- Correction of the roughness length computation over ocean in TMX
- Add multiscale stochastic pattern generator and iSPPT
- Bugfixes for restart issues with NWP sea ice
- Added diagnostics: pressure at convective cloud top/base
- Added C-cycle flags for C4MIP in XPP configuration

## ICON-Ocean

- Vectorization of HAMOCC subroutines
- Use consistent freezing temperature for sea ice in coupled XPP and AES
- Bugfixes for extended N-cycle in HAMOCC

## Soil and Surface

### Climate: ICON-Land

- QUINCY development
  - Enable anthropogenic land cover change with QUINCY biogeochemistry
  - Updated the Jena-Soil-Model for use with QUINCY as an alternative to the simple soil-biogeochemistry model
  - Improve paramaterization of vegetation phosphorus (P) uptake to avoid crops and natural vegetation dying by P limitation
  - Improve the calculation of plant water stress with frozen soil
  - Minor refactoring: unification of namelist names
  - Correction of aggregation of area-dependent variables in the QUINCY model
  - Update default PI control values for atmospheric 13CO2 and 14CO2 composition according to CMIP7 forcing (Graven, H. (2025))
  - Revisit QUINCY output ensuring that all output variables required for diagnostics and global budget calculations are present
  - Update nitrogen cycle parameters
  - Add script to generate CMIP7 based deposition data
  - Refactor vegetation memory structure to reduce code complexity
  - Assure restart identity
- Refactoring of anthropogenic land cover change process
- Hydrology: Added ford inline documentation
- Land initial files
  - Major update: 12 and 13 tile setups for jsbach and quincy
  - Fix for the skin layer conductivity
  - Estimation of initial soil moisture from vegetation fraction
  - Automized HD parameter file generation also for high resolution grids (internal HD)
- Land initialization: Initial soil moisture is turned to ice at temperatures below zero degrees.
- Revised 'basic' output list
- Anthropogenic emission files: Created anthropogenic emission data including aviation sources from the CMIP7 dataset
- Update BSD-3C licence year to 2026
- More flexible handling of the number of PFTs

## Infrastructure

### Building

- Fix building of the YAC Python interface when the Python interface of MTIME is disabled
- Consistently install ICON and the relevant bundled packages to the specified installation prefixes and subdirectories


# Release notes for icon-2025.10-1

## ICON-Atmo

### NWP Physics

- Bugfix in vdiff interface concerning restart reproducibility
- Introducing CO2 prognostic tracer to ecRad to enable emission-driven simulations with interactive land-ocean-atmosphere carbon cycle
- Allow more than one time interval for lateral boundary conditions
- Set taudecay in cover_koe separately for shallow, mid-level and deep convection
- Modified call to SR set_cdnc_from_extdata to get correct update for cloud_num in radiation after restarts
- Clean-up of mo_nwp_gscp_interface and 2mom microphysics scheme
- Fixes in Emvorado:
  - determination of nearest ICON cell
  - fix interface to eccodes routine codes_open_file

## ICON-Ocean

- Use ice class as vertical axis attribute for 3d sea ice variables

## ICON-ART

- Update the testsuite setup and scripts

## Externals

- Switch to ComIn 0.4.0

## Infrastructure

- Add distributed GRIB decoding
- Merge 2-Moment Microphysics Type Extensions into Parent Types

### Building

- Clean up and clarify single-precision options
- Fixes:
  - Handle apostrophes in the hostname


# Release notes for icon-2025.10

## ICON-Atmo

### AES Physics

- Implemented a new convective boundary layer (CBL) test case for TMX validation
- Optimization of the read-in of ozone data
- Updated the atm_qubicc_test experiment and its checksuite
- Re-introduced optimized simple plumes
- Implemented full (interactive) carbon cycle with hamocc and land coupling, anthropogenic
  emissions and CO2 tracer transport
- Implemented new rain microphysics
- Added gustiness to surface turbulent exchange
- Decreased albedo of ocean water from 0.07 to 0.05
- Implemented riming of cloud water and snow and faster snow fall velocity
- Refactored implementation of inhomx factors
- Removed snow from the radiation
- Ported AES Thermodynamics to Kokkos
- Fixes for:
  - GPU port of solar_parameters function
  - reading ozone input data
  - atmospheric energy diagnostics for GPU
  - tropopause height calculation when using vertical nesting
  - initialization of TMX with nested mesh dt_loc
  - surface exchange coefficients in TMX
  - Implemented additional checks and fixes for using 6-hourly SST/SIC data

### NWP Physics

- Allow for external specification of trace gases in ecrad
- New mixed-phase Spectral Bin Microphysics
- enable reading and interpolating ozone and Kinne aerosol data using YAC
- Tuning changes for better prediction of fog / low stratus
- Extension of adaptive parameter tuning to reduce cold/moist bias around sunset
- Improved tuning of filtering time scales for adaptive parameter tuning
- Modified shallow cloud cover as function of cloud droplet number
- New deep convection options to improve precipitation in the tropics
- Use cdnc scaling factor of the year 1850 in picontrol mode
- Extensions for prognostic 2D-Aerosol Prog2DAero
- Include a namelist switch to enable linear dependency of gwd momentum flux on precipitation
- Regime-dependent FSD parameterization
- More flexible options for entrainment/detrainment tuning in convection scheme
- Implementing interactive carbon-cycle
- Fixes:
  - Inconsistent initialization time step length in vdiff interface and jsbach
  - Use functions for sat. vapor pressure consistently between NWP physics and VDIFF/ICON-Land (ICON-XPP)
  - Make ecRad compilable with Cray 17.0.1 for AMD GPUs
  - Fix echotop and echotopinm diagnostics for clouds reaching model top
  - fix P-E lake double counting in JSBACH
  - LHN: Correction of some smaller bugs in observation processing
  - LHN: event bugfix
  - Bug fix for nudging interpolation
  - Bug fix in pressure bias correction for IFS lateral boundary data
  - Bugfix in mo_nonhydro_state for turbdiff, affecting SCM

## ICON-Ocean

- Add optimised variant of the CG-solver for ocean sea surface height
- Improve sea ice evp solver convergence
- Routines to read in First Guess and Analysis Files in the Ocean
- GPU, OMP, vectorisation improvements of ocean and sea ice code
- Optional fillValues on dry ocean cells
- Output of vertical integrals in HAMOCC
- Compilation of sea ice thermodynamics for atmosphere-only builds
- Fixes:
  - Fix dimension mismatch in ocean analytic forcing
  - Prevent Langmuir cell depth from becoming zero
  - Bugfix for the mixed layer depth diagnostics on vector/gpu
  - Fix concurrent HAMOCC
  - Correct vertical axis attributes for some HAMOCC and sea ice variables

## ICON-Waves

- Implementation of Stokes depth diagnostic
- new diagnostic: peak wavenumber
- Introducing the ocean layers table in Icon-Waves
- Restructure wave model initialization

## Soil and Surface

### Climate: ICON-Land

- QUINCY development
  - Integrated QUINCY processes with JSBACH physics processes (HYDRO, SEB, SSE, TURB, RAD)
  - Implementation of land management processes incl. reading of land-use data
    - Forest management (harvest)
    - Agriculture (8 crop functional types, cropland phenology, harvest and fertilisation)
  - First implementation allowing anthropogenic land cover change with QUINCY biogeophysics processes (canopy mode)
  - Improved calibration of nitrogen biogeochemistry
  - Enabled reading of biogeochemical vegetation and soil states from input file (external spin-up)
  - Enabled setting of selected calibration vegetation and soil-biogeochemistry parameters via namelist
  - Reduced number of restart variables
  - Refactored canopy configuration of QUINCY to permit calculation on GPU for standalone and coupled runs
  - Fixes:
    - Fixes in turnover, mortality and growth to prevent negative vegetation pools with extremely small fluxes
    - Fixed unit conversion for land CO2 flux in QUINCY
- Implemented a common interface for QUINCY and JSBACH
- Updated the anthropogenic land cover change process, including reading next years land cover map as target in case of daily land cover changes
- Added JSBACH usecase with 12 PFTs now including C4 crops
- Updated inline documentation of hydrology process
- Fixes:
  - Verification check for some JSBACH lctlib parameters when compiling with single precision
  - Account for proc0_shift from the parallel namelist when running ICON-Land standalone
  - OpenACC fixes for lumi
  - Fix needed with jsbach standalone simulations on GPUs
  - Fixed sequence of variables in surface temperature routine calls
  - Removed unnecessary mpi-all-reduce calls to reduce run time
  - Fixed uninitialized value for function get_time_dt in standalone model
- Code cleaning:
  - Only use one variable representing time step length
  - Revised surface temperature computation wrt variable names and comments
- Updated JSBACH usecase with TMX and PFTs
- Added support for CUDA graphs with AES physics (VDIFF and TMX)
- Implemented fix for using older restart files when not using skin temperature scheme (TMX or standalone)
- Surface water ponds (if enabled) now modify the top layer soil heat capacity
- Updated soil and root depths in land initial files
- Revised effect of phase change in snow, surface water and top layer soil storages
- Updated snow density parametrization to include a dependency on mean snow temperature

### NWP: TERRA and other surface issues

- Fixes:
  - rime formation term for interception storage
  - w_i nonconservation
  - bug fix for bare-soil evaporation

## ICON-ART

- Implement heat emission from wildfire as sensible heat flux from the surface
- Update of OEM and VPRM code
- no limit for number of point sources
- Fixes:
  - check if numbers in xml are numbers
  - removed clipping of tracers in diagnostics interface
  - optical properties
  - colum calculation in full chemistry

## Coupling

- Optimization of synchronization between atmosphere and ocean in coupled configurations
- Expose valid_masks in output_coupling
- Implement component-specific finish.status files
- Fixes for using vertex-based 3d fields in the output coupling
- Add support for AES model start date being shifted back in time by timeshift
- Bugfix to ensure restartability of coupled ocean-atmosphere configurations
- Fix overwriting of Tskin after ocean coupling

## Externals

- Switch to fortran-support 2.2.0
- Switch to CDI 2.5.2.1
- Switch to HD v5.2.4
- Switch to probtest v1.1
- Switch to YAC v3.9.2_p2
- Switch to mkexp 1.4.3

## Infrastructure

- Add option to build ICON in single precision (experimental)
- single precision for ecRad
- Add infrastructure for the Kokkos C++ development (Ragnarok)
- Add infrastructure for unit testing
- Stochastic pattern generator based on spherical harmonics
- refactor SST(SIC reader and time interpolation classes
- Read/Write RBF coefficients from/to file
- Port CDNC Scaling to GPU
- CUDA graphs support for NWP seaice
- Empa GPU developments
- Fix for the CUDA 13 release coming with NVHPC 25.9

### Scripting and testing

- Added AMIP test for TMX with PFTs in ICON-Land
- Added 158, 49, 5 and 2.1 km AMIP template with AES physics for mkexp (cpu + gpu, distributed IO)
- Added 49km AMIP BuildBot test with AES physics on one Levante GPU
- Support more ocean setups and output options in omip run script template
- Allow mkexp experiments to use CUDA graph utilization when configured for current host environment
- Refactor QUBICC test scripts
- Add low-resolution mkexp test setups for AES based bubble, Radiative Convective Equilibrium, Aquaplanet, AMIP, nested, A/O coupled and OMIP
- New AMIP setup in exp.aes_amip
- Re-enable restart mechanism in AMIP-runscripts, generated with make_runscripts

### Building

- Find and check MPI_LAUNCH command in the configure script
- Detect GPU architecture in the configure script
- Configure CUDA/HIP C++ compiler and flags
- Pass ICON_LDFLAGS to the Fortran compiler only
- Fixes for disable switches --disable-aes and --disable-jsbach
- Drop the HIP event handling suppression of the Cray OpenACC runtime
- Fix cached configuration

### Miscellaneous

- Some OpenACC optimizations and fixes
- Cleaned up obsolete workarounds for NVHPC/PGI compiler (OpenACC and OpenMP)
- Further extension of skin temperature over the ocean in Vdiff
- Add option to disable GPU memory usage output
- Fixes for the quad-precision handling


# Release notes for icon-2025.04-2

## Infrastructure

- Increase buffer size for auxiliary field during nest start


# Release notes for icon-2025.04-1

## ICON-Atmo

### NWP Physics

- Fixes:
  - Enables using spun-up FLake variables from first guess in MODE_COMBINED

## ICON-ART

Fixes:
  - changed used height array (z_mc -> z_ifc) for column computation in chemistry
  - changed usage of tracer names to mode names in diagnostics routine regarding optical properties
  - fixed setup information for standard configurations

## Externals

- Switch to the latest version of ICON-TIXI
- Switch to YAXT 0.11.4


# Release notes for icon-2025.04

## ICON-Atmo

### AES Physics

- Add new radiation fields:
    - rsntcs: toa net clear-sky shortwave
    - rsnscs: surface net clear sky shortwave
    - rlntcs: toa net clear-sky longwave
    - rlnscs: surface net longwave
- Added support for time steps with fractional seconds in AES/ICON-Land
- Add diagnostics for cloud condensed water and total water mass fraction
- OpenACC port of solar_parameters() and pre_rte_rrtmgp_radiation()
- Fixes for:
    - TMX: GPU and OpenMP, \_wp inconsistency, gcc-14
    - Bubble test
    - Vertical integrals for cloud condensed water + ice and related diagnostics
    - Total ice in cloud microphysics
    - Phillips nucleation von GPU
    - Initialization for nesting and latbc
    - Nesting and restarting with sea-ice on GPU

### NWP Physics

- 1-moment Microphysics
  - Option to use more accurate coefficients for saturation pressure from IFS: itype_satpres_coeffs
  - Three new namelist variables for tuning of graupel microphysics to reduce overprediction of high precip intensities:
     - tune_supsat_limfac: allows for supersaturation in updrafts in the saturation adjustment
     - lvariable_rain_n0: if .TRUE., the variable intercept parameter is activated. The multiplicative factor rain_n0_factor is used for drizzle (small qr) while the default value is approached for heavy rain (large qr).
     - tune_box_ice: scales the ice cloud cover in cloud cover scheme
  - GPU port for icpl_aero_gscp = 3  (MODIS climatology for cloud-droplet number)
- 2-moment Microphysics
  - Aerosol-cloud interaction for cloudice2mom using ART dust for ice nucleation
     - icpl_aero_ice=3: simplified coupling with fixed modal diameters
     - icpl_aero_ice=4: full coupling with aerosol size information from ART
  - Removed a confusion between specific mass and mean mass in the parameterization of homogeneous ice nucleation in cloudice2mom (gscp=3)
  - Added consistent treatment of effective radius for the two-moment cloud ice scheme (gscp=3).
    In addition, the cloud ice number source for the snow drift term is introduced.
  - New ensemble perturbation parameters are added (range_ccn_Ncn0, range_in_fact, range_avel_i, range_avel_g, range_cap_snow, range_cap_ice).
- Radiation
  - CAMS climatology: The format of the CAMS aerosol climatologies has changed, and aerosol is now provided as mixing ratio (kg/kg) instead of layer integrated mass.
  - ecRad: new logical namelist variable ecrad_check_input: if .TRUE., several input fields are checked for physical consistency and the verbosity of ecRad is increased.
- Convection: new namelist tune_grzdc_offset: Scaling factor for offset in CAPE closure for grayzone deep convection. Positive values reduce the activity of the convection scheme and suppress convective drizzle (recommendation: 0.1--0.2)
- Turbulence: The complexity of the operational turbulence scheme for NWP has been reduced in several steps.
- SSO: Fix of array access for inwp_sso=2 in a diagnostic computation.
- Tuning: New NWP physics options for tuning and improved physical consistency
  - itype_dissip_heat=2: new option to take into account dissipative heating from turbulent momentum dissipation in the turbulence interface
  - icpl_o3_tp=2: improved option for ozone-tropopause coupling to avoid excessive additional ozone for low tropopauses (in combination with itune_o3=3)
  - shift_ratsea: option to shift ensemble mean of rat_sea w.r.t. the deterministic value
  - shift_boxliq_asy: option to shift ensemble mean of tune_box_liq_asy w.r.t. the deterministic value
- Diagnostics
  - Update visibility diagnostic to be consistent with icpl_rad_reff=1 (typical RUC settings)
  - GPU port of the diagnosis of planetary boundary layer height
  - Add accumulated and max/min diagnostic fields to nest start interpolation
  - EMVORADO: Bugfix call to polarimetric dbz diagnostic in case of Tmatrix

### Modifications by the CLM Community

- Fixes in soil-moisture dependent albedo tuning to enable reproducibilty of ICON-CLM on GPU

## ICON-Ocean

- Add diagnostics for all terms of the temperature and salinity budget
- Add diagnostics for upper ocean heat content (hc300m and hc700m)
- Add new output variables (tos, sos, sivol, snvol)
- Add GRIB codes for ocean variables (mld, mlotst, normal_velocity, stretch_c, hctm, hc300m, hc700m, snhc, sihc)
- Improved GPU code and performance optimizations and validation tests
- More consistent precision of literal floating point values used in TKE and IDEMIX
- Bugfix correct unit conversion for tos
- Bugfix for dynamic short wave radiation absorption
- Bugfix when initializing the ocean from restart files

## ICON-Waves

- Memory layout (blocks as last dimension) and runtime improvements (precomputation of expensive operations where possible)
- Revision of wave-atmosphere coupling, specifically for the usage of z0 (roughness length) passed from the wave model to the atmosphere
- Wave initialization and coupling:
  - Similar to the atmosphere model, the wave model start date can be shifted backwards in time (new Namelist parameter dt_shift in initwave_nml)
  - When shifting the YAC startdate backwards in time as well (using src_lag/tgt_lag to match with dt_shift of the coupled models) it is possible to start the coupling at the very first integration step.
- Implement depth-limited level of energy
- Add timers for ICON-waves
- Several fixes for:
  - Asynchronous output writing in coupled mode
  - Calculation of 10m wind direction for wind forcing
  - Implementation of a minimum allowed level of wave energy for each frequency
  - Setting of the output source time level for prognostic wave energy
  - Usage of bathymetry versus water depth field
  - Formula for Charnock output parameter

## Soil and Surface

### Climate: ICON-Land

- Added per-process namelist option lrestart_cont to allow restarting from other experiments run
  without that process
- QUINCY development
  - Added experiment file for ICON-Land standalone runs using QUINCY in CANOPY mode (no biogeochemistry)
  - Initial GPU port of QUINCY - running in CANOPY mode with the ICON-Land standalone driver
  - Minor scientific updates
    - Improvement in the first soil-layer hydrology
    - Bugfix: snow melt calculations
    - Bugfix: calculation of grassland phenology
    - Bugfix: added minimum level of C limitation on nitrification and denitrification
    - Bugfix: diffusion water flux limitation in QUINCY soil physics
    - Fixed calculation of saturated water content from input data
    - Improvements in the computation of several rate modifiers used in soil biogeochemistry calculations
    - Clean-up of calculation of stand-replacing harvest
    - Runtime optimisation: reduced number of aggregated variables
    - Reduced number of variables in the restart file
    - Included forcing and output of carbon isotopes
    - Read elevation external parameter for QUINCY from file
    - Improved handling of n and p deposition reading from forcing data
    - Inclusion of self-thinning and herbivory in grasslands and pastures (but not crops)
    - Improvements regarding C:N and N:P ratios in leaves and soil organic matter
  - Merged the radiation process of QUINCY into the radiation process of JSBACH
  - Use JSBACH4 canopy, soil and snow albedos with QUINCY albedo calculations
  - Use the turbulence process of JSBACH in QUINCY replacing QUINCY-specific turbulence-code
  - Preparations for using further JSBACH physics processes with QUINCY
  - Consolidated and cleaned up namelist handling and physical parameters between QUINCY and JSBACH
  - Implementation of a harvest process for QUINCY (for now using a global constant)
  - Technical implementation of an agriculture process for QUINCY (scientifically not yet ready for use)
  - Bugfix: static reals were missing decimal
  - Bugfix: some local REAL variables were missing kind statement
- Made PFT parameters available in memory init functions
- New optional tag for the memory usage report
- New functions for time control: get_previous_month_length and get_previous_year_length
- Memory reduction: array allocation only if needed with the specific setup
- The script suite to generate ICON-Land input data now also includes scripts to generate
  HD parameter files (for internal HD) and HD receive masks (for external HD), besides
  high-resolution Merit-Rema topography data is used.
- Added diagnostic variable for volumetric soil moisture content for soil layers
- Represent soil ice as ice volume, not as water equivalent anymore, thereby fixing
  an energy balance inconsistency during soil ice melt
- Changed handling of excess soil moisture
- Changes and fixes for inline documentation
- Improved vectorization on NEC machines
- Clean-up of ICON-Land code
- Introduction of an output group for jsbach monitoring variables
- Implemented soil hydrology parametrization for uniform scale
- Updated thaw depth diagnostics
- Enabled JSBACH usecase with PFTs when using TMX
- Interface: New switch to supress YAC call during initialization phase
- Updated HD-YAC coupling interface to use mo_coupling_utils
- Implemented experimental formulation of skin temperature
- Fixes:
  - Bare soil evaporation and modification of roughness (heat) and photosynthetic efficiency parameters
  - Calculation of snow aging
  - Computation of fast drainage within the ARNO scheme of JSBACH hydrology
  - For supporting single precision
  - Initialization of carbon pools from file (read_cpools)
  - ICON-Land standalone concerning nproma
  - Simulations with JSBACH assimilation and LAI prescribed from climatology
  - HD global water conservation test

### NWP: TERRA and other surface issues

- Modularization: The one long TERRA routine has been split into smaller subroutines containing one task each. Some issues for the water budget have been fixed
- nwp_sfc_interface: move accumulation of runoff_[sg] to nwp_statistics
- itype_ahf=3: Option for time-dependent specification of anthropogenic heat flux based on time-filtered T2M
- Improve coupled-model water conservation by adding TERRA's water nonconservation diagnostic to the subsurface runoff that gets passed to the HD model.
- Implementation of ocean surface layer parameterization: warm layer, cold skin
- Fix for extpar data and the usage of snow analysis increments on subgrid-scale glacier points

## ICON-ART

- Implemented LinozV3 ozone chemistry parameterization
- Added Upper Boundary NOy to stratospheric SimNOy scheme
- Enhanced simplified OH chemistry with O(1D) calculation and consideration of CFC
- Fixed bugs in LinozV2 polar chemistry
- Added a dedicated module for chemical constants
- Improved flexibility and configuration for prescribed aerosol optical properties
- Introduced Subpollen Particles (SPP) parameterization
- Updated FPlume module
- Improved robustness of dust radiation calculations
- Implemented a new dry deposition scheme for gases
- Added wet deposition processes for trace gases

### Coupling

- First implementation of coupling the nested AES atmosphere to the ocean model
- Support coupled setups with one component starting from IAU
- Support component wise model initialization from restarts
- Improvements for coupling timers, metadata handling and GPU setups
- Output coupling improvements for prognostic variables, variables on pressure levels and more vertical grids
- Add coupling infrastructure for one-way coupling with Super-Droplet cloud micropphysics model CLEO

## Externals

- Replace math-support and math-interpolation with iconmath 1.2.0
- Make use of the math-horizontal component of iconmath 1.2.0
- Switch to fortran-support 2.1.0
- Switch to mtime 1.3.0
- Switch to YAXT 0.11.2
- Switch to CDI 2.5.2
- Switch to Comin 0.3.0
- Switch to YAC v3.6.2_p2
- Switch to the latest version of HD

## Infrastructure

- Several improvements for single precision: communication, mpi, io
- Fix the generation of index lists on GPUs for LAM and nested simulation runs

### Scripting and testing

- Update land data for simulations with ICON-Land on R02B04 grid 0049
- Add check for memory consumption on vector machines
- Add tests for LAM/Nest and AMIP
- Add offline check for global water conservation of coupled ICON XPP
- Add 40km mesh size (R02B06) for AMIP development test
- Add mkexp setup for ICON-Land standalone and AMIP-style (NWP ATM)
- Add 10km AMIP template for mkexp (cpu + gpu, distributed IO)
- Restructure mkexp templates and add bubble test
- Add log file monitoring to mkexp

### Building

- Improve building on macOS
- Update the generic configure wrapper and its documentation

### Miscellaneous

- Improve support for Cray compiler 17+ for AMD GPUs
- Single precision improvements for TMX
- Add contribution guidelines
- Add various issue and merge request templates
- Replace Poisson routine in stochastic NWP convection code
- Add logging of GPU memory usage
- Clean up documentation
- Integration of online documentation


# Release notes for icon-2024.10

The following lists give an overview on the main changes since the last release icon-2024.07.
Note that this release now also contains the external model HD-couple, which is ready for
open source now.

## ICON-Atmo

### DyCore

- Revise projection to tangent plane for the FFSL scheme
- Algorithmic optimization of the MIURA3 transport scheme
- Optimize quadrature routines for tracer transport
- Bug fix for linear advection quadrature
- Namelist option for CFL monitoring frequency
- Fixed accumulation of small epsilon on contravariant mass flux (in PPM)
- AES: Fix faulty usage of NWP variable `prm_diag` with nested domains
- OpenACC bugfix in interpolation of ozone from pressure levels to model levels
- Cleaned up some time-related constants (wrong place / doubled definitions)

### NWP Physics

- Extension of adaptive parameter tuning
- Preparing a major revision of the NWP turbulence code including the integration of some not yet
  considered effects of surface roughness
- Tuning for prognostic 2D aerosol scheme
- VDIFF turbulence: deep-atmosphere fixes
  - Consistently use full geopotential when converting between dry static energy and temperature
  - Relax limits on vapor pressure table lookups
- Cleanup in radiation and aerosol code parts
- Initialise aerosol fields in case of Kinne/CAMS Aerosol

### AES Physics

- VDIFF turbulence: deep-atmosphere fixes (similar to NWP Physics)
  - Consistently use full geopotential when converting between dry static energy and temperature
  - Relax limits on vapor pressure table lookups
- TMX turbulence
  - Refactoring for increased modularization
  - Add 2m dewpoint temperature diagnostic
  - Fixes for OpenACC
- Add diagnostics to trace atmospheric energy
- Re-activate output of aerosol optical properties with RTE-RRTMGP
- Revise clear sky radiation computations

## ICON-Ocean

- Add ocean isopycnal transport diagnostic
- Bug fix: calculate sea water density always on fixed depth levels
- Bug fix and refactoring of the ocean age tracer
- Optimize ocean surface solver
- Continue GPU porting of ocean
- ICON-Waves: prepare coupling of surface waves to the ocean
- ICON-Waves: Add restart and checkpointing functionality

## Soil and Surface

### Climate: ICON-Land

- Fixes for using older restart files from before the JSBACH pond scheme was implemented
- Improvements in JSBACH soil hydrology
  - Change lower and upper limits of soil moisture
  - Add option to use uniform distribution of soil moisture for infiltration and drainage as alternative
    to semi-distributed parameterization that accounts for sub-grid variability (Arno scheme)
  - New option to force initialization of soil moisture from a file instead of from IFS analysis
  - Fixes for OpenACC loops in JSBACH hydrology
  - Bug fixes for JSBACH pond scheme
- Update of the scripts to generate ICON-Land initial (ic) and boundary condition (bc) files
- Implement daily execution of anthropogenic land cover change by interpolation of annual maps
- Allow running both: anthropogenic and natural land cover change
- QUINCY development
  - Refactoring of the quincy soil physics process
  - Updates incl. first implementation of coupling with ICON-Atmo
  - Implementation of a spin-up accelerator for the slow biogeochemical soil pools
  - Implementation of wood product pools
  - Implementation of a carbon conservation test
  - Minor code fixes towards usability and style recommendations
  - Minor scientific updates
    - First step to include stem area (SAI, stem area index) into radiation scheme
    - Fix of slow growth in early season in cold grassland sites
    - Calibration of self-thinning for trees
- Switch from deprecated YAC interface in HD model
- Several fixes for DSL pre-processor script `dsl4jsb.py`
- Port JSBACH carbon and disturbance modules to GPU
- New options to reduce diagnostic output in log file from water balance checks
- Fix for the land cover fraction diagnostics of simulations with natural or anthropogenic
  land cover change
- Fix too cold soil temperatures for partially snow-covered grid cells
- Fixes for natural land cover change

### NWP: TERRA

- Encapsulate initialization of land use-related parameters for NWP
- Add option for ICON-internal soil moisture adjustment

## Externals

- Update HD model
- ComIn 0.2.0
- Updates for YAC
  - Switch to version 3.4.0_p2
  - Fix output coupling and enable python interfaces (yac,mtime,comin) for testing
- Introduce the math-support and math-interpolation libraries
- Update to MTIME 1.2.2

## Infrastructure

- MPI: check worker architecture during communicator creation
- Restructured Subroutine initicon_inverse_post_op
- Update the mechanism for source provenance collection

### Coupling

- Add detailed timers for output_coupling
- Fix OpenACC bugs that affect coupled het jobs
- Do runoff diagnostic only when new data are received from YAC
- Interface aes/ocean: move ocean coupling call from `mo_interface_iconam_aes.f90` to `mo_nh_stepping.f90`

### Scripting and testing

- Fix component names in mkexp experiment to prevent wrong coupling setups
- New option -P for mars4icon_smi to use surface pressure instead of lnsp
- Activate test for lgrayzone_deepconv
- Several improvements to the experiment setup with the `mkexp` run script generation system
- Adjust LUMI-G defaults in create_target_header
- Additional experiment runscript template `run/exp.atm_nwp_jsbach-C` to test carbon cycle with ICON-XPP
   (NWP atmosphere simulations with jsbach)
- DCMIP Tropical Cyclone experiments: removed `dcmip_tc_51` test case and activated
   `dcmip_tc_52` also for AES physics
- Clean up and consolidate several run script templates with AES physics (NextGEMS, AMIP, nested, land)
- Update of JSC run scripts
- Fix buildbot test scripts for Juwels and Booster

### Building

- Several minor fixes for the configure script
- Fix `USE mtime`, `INCLUDE netcdf.inc` and check for NetCDF Fortran 77 API
- New makefile target 'env' to retrieve the build environment (`BUILD_ENV`) set in a configure wrapper
- Expose BUILD_ENV to the runscript generators
- Several minor fixes and improvements for the build system
- Introduce configure option `--enable-bundled-python` to build the Python
   interfaces of `MTIME`, `YAC` and `COMIN`

### GPU port, technical developments and optimizations

- GPU port of radiation namelist switch irad_o3=5
- Prepare for eccodes versions >= 2.32.0
- GPU port for stratocumulus tuning parameters `tune_sc_*`
- Disentangle ext_data state construction, separates the state construction from its initialization
- Adjust nproma to 256 B alignment to improve GPU performance
- Optionally suppress HIP event handling of the Cray OpenACC runtime
- Move and split CUDA/HIP source files
- Several bug fixes for OpenMP and OpenACC


# Release notes for icon-2024.07

These are the release notes of the ICON model.
Below the main changes as compared to the previous ICON release.

## ICON-Atmo

### DyCore

- Improve adaptive CFL reduction at model start
- Update and extend diagnostics in supervise_total_integrals_nh
- Fix OpenACC loops in vertical advection

### NWP Physics

- Improvements to adaptive parameter tuning
- Implemented option for Cloud droplet number from MODIS climatology
- RUC (Rapid Update Cycle) diagnostics, fixes and testcases
- DeMott ice nucleation parametrization
- New output variables: diagnostic tot_pr_max; surface radiative fluxes without islope_rad corrections
- changed random normal values generation in Stochastically Perturbed Physics Tendencies (SPPT)
- Gust diagnosis option for large-eddy permitting configurations
- Option for moisture diffusion (water vapor and cloud water)
- 2-moment-microphysics: new hail/graupel shedding and limitation of graupel production
- Tuning parameters for TERRA-URB (Soil parameterization) and wind gusts
- fix atmospheric water budget issues
- Optimize emvorado compile time, bugfix and new namelist parameters
- Aerosol optical depth (AOD) output from Kinne and enabling fr_glac
- Computation of instantaneous grid scale precipitation rate at every physics time step
- Include parametrization of wave breaking into the ICON-waves physics
- GPU Port of the `ldass_lhn` switch for the two-moment microphysics scheme
- Fix for atmo-wave roughness length coupling
- Wave-dependent sea surface roughness in icon-nwp
- Fix initialisation of aerosol fields in case of Kinne/CAMS Aerosol

### Modifications by the CLM Community

- CDNC scaling for climate projections
- Removed the alteration of Modis Cdnc
- Add namelist parameter (rat_lam,rsmin_fac) to enhance tuning capabilities
- Added correction for albedo dependent on soil moisture for soil types 3 to 6

### Various bug fixes

- Several bug fixes for OpenACC and Cray compiler
- Make limited-area model (LAM) and nested grids work on CPU
- Add fixes to graupel microphysics (synchronize with granule)
- Surface net longwave radiation flux
- init call with OpenMP in aes-ocean coupling
- nextGEMS prefinal
- Forgot loop indices in OpenMP PRIVATE statement for one loop in TMX turbulence package
- Removed inconsistency for turbulent diffusion coefficients between TTE and Smagorinsky schemes for VDIFF and TMX turbulence packages
- Replace two ACC KERNELS constructs with explicit loops to improve performance in TMX turbulence package
- Make sensible heat flux in surface energy balance consistent with atmosphere when using TMX turbulence package
- Remove ACC TILE constructs from single (not nested) loops in TMX turbulence package
- Performance fixes for OpenACC and fixes for OpenMP in TMX turbulence package
- Fixed some potentially non-contiguous OpenACC transfers in coupling code plus removed some OpenACC KERNELS constructs in VDIFF

### New features and other modifications

- Replace `lrad_yac` for coupling O3 and aerosols via python processes
- Add namelist parameter for lower tropospheric stability correction
- Add output diagnostic for relative humidity
- Expose `output_nml` variables via YAC
- Add diagnostic for 2m specific humidity in TMX turbulence package
- Add tuning options in TMX turbulence package (kinetic energy dissipation and stability correction) and change cloud droplet number concentration from 200 to 50
- Consistently apply the physical tendencies at constant mass and volume in the cloud microphysics (mig), the radiant energy transfer (rad) and the  turbulent mixing (vdf/tmx)
- Updates for energetics with TMX turbulence package - consistently use internal energy and account for the energetic contributions of temperature-varying phase change energy
- Various fixes and updates for TMX turbulence package (diffusion of internal energy, diffusion of u/v velocities)
- Added a cloud inhomogeneity factor for snow and clean-up of factors for cloud liquid water
- Remove unused code in the RTE-RRTMGP interface
- Change setup for Sapphire base/nextGEMS prefinal

## ICON-ART

- Changed initialization of chemtracers for init_gas=0 in ART (now sets tracers to 0)
- INAS ice nucleation scheme preparation, dusty cirrus, meteogram output
- Repair ICON-ART standard cases
- New plume rise model implementation as option...
- Clean-up of constants and adding emissions into wet air
- Remove usage of RRTM in ART
- Avoid double execution of the phenology update...
- Fix exponent sign bug in Hande et al 2016 CCN activation
- Fix for double allocation of some variables
- Bugfix for wrong pressure unit in PSC scheme in ART
- Water content computation of sea salt aerosol
- Refactored and GPU-ported OEM (Online Emission Module) code

## ICON-Ocean

- New feature: Add diagnostic for ocean bottom pressure
- New Feature: Add ocean potential temperature in Kelvin as optional output
- New feature: Add diagnostic age tracer and age tracer squared
- Bugfix: Enable sea-ice drag on ocean as default
- Bugfix: Reintroduce missing coupling fields in ocean restart files
- Refactoring: Improve initialization of 2D variables
- Call interface_aes_ocean only if coupled_to_ocean

## ICON-Land

- Adapted ICON-Land to the interface change of init and copy functions of mo_fortran_tools
- Bugfix: Fixed restartability of offline model
- New feature: Ported offline model to GPU
- Bugfix: Fixed an indexing bug in phenology
- Remove hard-coded number of soil layers in JSBACH and get it from input data instead
- New feature: Implementation of QUINCY biogeochemistry model
- Update of inline documentation and code cleaning
- Added various pre-processing scripts for ICON-Land/JSBACH
- New feature: Implemented ponds in JSBACH
- Replaced copyrighted code for solver of tridiagonal linear equation system
- New feature: Updates for ICON-Land offline simulations (filter for surface drag computation, refactor of driver)
- New feature: Ported (internal) Hydrologic Discharge (HD) model to GPU
- Bugfix: Fixed wrong vector assignment inside loop
- Use new versions of JSBACH input data for NextGEMS atmosphere-only experiments
- New feature: Improvements for JSBACH soil hydrology (soil organic matter, soil freezing and thawing, soil moisture, new solver for vertical water transport) plus cleanup of code and variable names
- Updates to the bookkeeping of how processes in ICON-Land run on different tiles
- New feature: Implementation of forest age classes
- Refactoring and fixes for the LCC (land cover change) framework including the processes for disturbances (wind/fire), fuel/carbon and anthropogenic land cover change
- OpenACC updates and fixes: use of ASYNCs and other changes to follow ICON OpenACC programming guidelines; remove some ACC WAITs for optimization

## Externals

- Updates to ecRad: gas model split & spectral weights; fixes from ecrad master
- CAMS forecasted aerosol as an option for ecRad radiation
- Update ComIn to version 0.1.1
- Update libfortran-support to version 1.2.0
- Update YAXT to version 0.10.2
- Update to CDI version 2.4.0
- Update to RTE+RRTMGP version 1.7

## Infrastructure

Several bug and technical fixes, below the most relevant ones:

- 1mom-microphysics granule: Build and Testing Infrastructure
- Extensions of the action feature, which is used in ICON to reinitialize variables at regular intervals
- Configure: try to find python3 before python
- Fixes bug in mo_communication_yaxt
- Fix for array bound violation in src/atm_phy_nwp/mo_nwp_phy_init.f90
- Improve configure for use with external comin
- Fix the time and date exposure when ComIn is enabled
- Fix abort_mpi: exit with a non-zero exit code on error
- Switch to NetCDF Fortran 90 interface

And some modifications to coupling:

- Moved coupler initialisation
- Refactoring of coupling
- Hd (Hydrological Discharge) updates
- Coupled TERRA - HD through YAC
- Online diagnose of global sum runoffs before and after the YAC coupling
- Correct calculation of Qtop and re-enable LSM patching in NWP-Ocean coupling

## GPU

- GPU performance fixes and workarounds for NVIDIA and AMD GPUs
- Removing majority of i_am_accel_node instances (deprecated debugging feature)
- Adding nblocks_e namelist option + updated logic for nproma, nblocks_c
- Cuda graph : activation through logicals + testing infrastructure
- Added WAITS for bit reproducibility on GPU


# Release notes for icon-2024.01-1

A workaround for the NVIDIA compiler has been implemented to be able to compile yac.

Necessary backports from the main ICON repository for buildbot testing are also included.


# Release notes for icon-2024.01

It took a lot of work, sweat and tears: ICON's first Open Source
(BSD 3 Clause) release is out!

Nevertheless, a lot of changes have been applied to icon-2024.01 since
the last released version icon-2.6.6. The given list might not contain
all changes (eg.bug fixes, ...), but is a good overview.

## Buildbot

- Consolidate probtest
- Read/write access for experiment directory on NEC
- MCH: Add mch_kenda-ch1_small experiment
- Reorganization of nwp buildbot data
- Refactor START_MODEL_function
- Migrate data pool on Balfrin
- Add qubicc tests to levante_gpu/cpu_nvhpc builders
- Add bubble update test
- Add mixed precision builders on balfrin with performance test
- Introduces the concept of one tolerance hash per experiment and builder
- Add build-only tests on LUMI
- New `aes_bubble_land` test case with land (JSBACH)
- Move the configuration of the bundled libraries to the build stage
- GPU port of tracer transport and HAMOCC functions
- Introduction of new experiment list **merge2rc** (== icon-dev + dwd)
- Introduction of new experiment lists: art, oce for specific component tests
- Several improvements and fixes for buildbot
- Add latest GPU build wrapper for levante
- Add the ocean-gpu Buildbot list for variants of ocean_omip experiments with NVHPC on Levante

## NWP

- Redesign of bias corretion of reference precipitation.
- Modified time dependence of EPS perturbations for LHN tuning parameters
- Revision of adaptive parameter tuning for seaice scheme
- Anthropogenic aerosol for 2D-aerosol scheme
- Coupling with hydrological model HD
- First set of tuning changes / new tuning options for ICON-D05
- adapt HZEROCL diagnostic for MCH
- WAVES:
  - new diagnostic output fields
  - Horizontal transport of ocean surface wave energy (part I)
- Bugfix in vertical interpolation in cams climatology
- WAVES: Coupling yac3 ICON-NWP <-> ICON-waves
- Modified z_pbl output for NWP
- CAMS climatological aerosol as an option for ecRad radiation
- Modification to rh based vis diagnostic and clean up of vis diagnostic.
- sfc_seaice: heat transfer from the sea-ice slab is modified
- EMVORADO: New option for dynamic wet growth T-limits in radar forward operator
- Additional namelist options for improved SSO/gravity-wave tuning
- JSBACH/VDIFF: add interface between NWP physics package and vdiff
- Remove optional mass fixer from (incremental) feedback routine
- Remove optional divergence averaging
- Remove the turbulence scheme EDMF from ICON
- Rewrite gamma functions to avoid license issues
- Remove RTTOV coefficients from /data directory due to license isssues
- EMVORADO:
  - Move two radar interface modules from ICON source code to submodule emvorado
  - update emvorado headers and gamma functions
- Integrate separate deep-atmosphere dycore into standard dycore
- Two-moment scheme changes for RUC and CLC
- Implementation of warm-rain spectral bin microphysics (SBM)
- Generic hydrometors for ecrad
- latitude-dependent decorrelation length scale for cloud overlap
- set of changes needed to allow numerically stable integrations at mesh sizes below about 100 m
- waves: implement full suite of parameterizations
- waves standalone: Read-in of forcing data
- new namelist parameter `tune_capethresh`
- ECRAD: updated version by ECMWF
- ECRAD: activate ecckd gas optics
- TERRA-URB: anthropogenic heat flux added
- Latent heat nudging: fixed several incorrect loop boundaries (nproma -> i_endidx)
- Enable reduced radiation grid for ART+ecrad
- wave model: additional model components added
  - ext_data state, read in, output
  - wave energy initialization
  - time loop
- Removed `__ICON__`, `__COSMO__` and `HAVE_FC_ATTRIBUTE_CONTIGUOUS` macros
- Remove optional mass fixer from (incremental) feedback routine
- Remove optional divergence averaging
- Remove the turbulence scheme EDMF from ICON
- Remove RTTOV coefficients from /data directory due to license isssues
- EMVORADO:
  - Move two radar interface modules from ICON source code to submodule emvorado
  - Update emvorado headers and gamma functions
- EMVORADO: Fix serious problem with the new dt_obs(1:3) notation in automatic reading from radar meta data from obs file
- ICON-seamless prototype 2
- Revised NWP ocean interface
- Updates to the two-moment microphysics scheme
- New option for horizontal Smagorinsky diffusion on vertical wind speed
- TERRA-URB: prevent evaporation from bare soil
- Option to modify diagnostic cloud scheme to enhance cloud cover in stratocumulus regions
- Add option for slope-dependent radiation with shading but no sky-view factor effects
- waves: add complete set of wave physics
- MVSTREAM: Horizontal mean for local grids (LAM)
- New namelist switch for ozone tuning
- New namelist parameter for specifying the CFL-W threshold for adaptive time step reduction
- New configure option for disabling/enabling the NWP physics package
- iterative IAU: revised implementation by looping over perform_nh_stepping
- Cleanup advection schemes used for density and potential temperature
- Remove optional open upper boundary condition l_open_ubc=.TRUE. and remnants of the hexagonal code
- DACE: Changes to allow DACE to be used on daint on GPU
- securing that cloud water is not being diffused above the level 'kstart_moist'
- Bug fix for nest initialization during runtime: add filtered wind speed increment

## GPU

- GPU port of the visibility diagnostic
- GPU port of the explicit two-moment microphysics scheme
- Fix uninitialized data reads in kenda-ch1_small
- port of Hailcast
- initial ART port and port of the pollen-related processes to GPU
- add missing ASYNC(1) clause to PARALLEL regions
- Fix out-of-bound memory access in vertical advection for GPU
- CUDA Graphs for Terra and turbtran and GPU optimizations
- OpenACC port of totint
- Optimized OpenACC port of two-moment microphysics scheme and LHN
- Make ICON buildable on LUMI
- Add GPU testing for levante to the default test suites
- OpenACC port of
  - vertical output interpolation
  - 3D turbulence scheme
  - diagnostics
  - SPPT
  - sstice_mode=6 (update with user-defined interval)
- fix for iterative IAU test on GPU
- Workaround for cray mpi bug using async/IO

## MPIM Atmosphere

- Rewrite and corrections of one-moment microphysics (graupel)
- ICON-Land update for land offline simulations
- Separate orbit and solar calculations from radiation in `src/atm_phy_rte_rrtmgp` to enable any combination of `--en/disable-jsbach` and `--en/disable-rte-rrtmgp`
- Fix potential array bounds error in ocean
- Revision of the interface structure for the AES physics
- clean up the Hadley test runscript
- remove graupel (mig) namelist parameter that are no longer in use
- collect changes for nextGEMS cycle 3 except hiopy and compresm
- Refactoring of coupling via YAC

## MPIM Ocean

- refactor TKE (ocean)
- port the ocean time loop to GPU
- OpenACC port of the ocean tracer transport
- OpenACC port of LVECTOR code used in ocean tracer transport
- enable concurrent hamocc to work with the ICON-O z-star levels
- Bugfixes and refactoring for MOC diagnostic

## JSBACH

- Replaced CLAW directives and preprocessing with pure OpenACC
- JSBACH workarounds for Cray compiler on lumi GPU
- Prepare ICON-Land for `mo_exception` API change and new math library in ICON-C
- mprove NEC vectorization in JSBACH for ICON-Seamless
- fix a bug in carbon conservation test
- harmonize scripts for ICON-Land ic and bc file generation

## Update yac2 -> yac3

- XML config file is replaced by yaml
- all coupled experiments in buildbot are converted to yaml
- many non-buildbot tests are converted to yaml

## ART

- Ensemble and VPRM functionality of the Online Emission Module
- Extension of the pollen treatment with hazel (CORY)
- Tuning, additional diagnostics and introduction of SI units for dust
- Possibility of choosing between 32- and 64-bit float output applied to ART
- XML reading accepts comments in XML-file
- Moving of _src/art_interface_ from the ICON-Code to _interface_ in the ART-Code
- Diagnostics for emissions and washout and reintroducing accidentally lost emission call
- Coagulation and mode shifting to mixed modes
- Aerosol climatology naming

# Release notes for icon-2.6.6

Later than expected: this is another ICON release

- OpenACC port for Nvidia GPUs and enabling HPE-Cray OpenACC use with AMD GPUs
  - almost all of the ICON atmospheric components are ready for use on
    GPU equiped machines. This has been a big effort by Nvidia,
    HPE/Cray, MeteoSwiss, DWD, CSCS, CSC, and MPIM (the order is
    arbitrary).
  - there is further progress coming up
- CSCS
  - preparation of daint/dom transition for Nvidia compiler bumb to 22.5.
  - added Docker files for EXCLAIM project in scripts/docker/exclaim
  - more test cases
- DWD
  - still persisting problem with OpenMP (non-reproducability) in buildbot tests.
  - bufixes and improvements in the data assimilation part and physics
  - two-moment cloud microphysics has been replaced by a new version
  - initial steps for wave model implementation
  - added jsbach/vdiff from former ICON echam-physics
  - enable/fix more problems with respect to run icon-seamless
  - some code moves and consolidation
  - OpenACC beautification application running to allow passing a required CI test
  - more test cases
- MPIM
  - bugfixes in the ocean
  - added additional diagnostics
  - improvements in JSBACH (HD mask, LCC, former dynveg)
  - updating of build on levante and macos
  - fixes for claw
  - fixes for building ICON with some data assimilation components
  - Nvidia 22.9 is needed for JSBACH
  - ocean updated to zstar vertical coordinate
  - bump rte-rrtmgp to new version
- DKRZ
  - fix many errors and enabled new features with respect to coupling external components
  - new and updated config-wrappers
  - improve buildbot testing

For details take a look at:

https://gitlab.dkrz.de/icon/wiki/-/wikis/Protocol-of-Release-Commits

# Release notes for icon-2.6.5

It is time for yet another ICON release.

After 10 month we managed to provide a new ICON release with new
machine configurations, a new buildbot, significantly extended,
improved GPU capabilities, first steps into ICON seamless and a
cleaned up version of the MPIM AES physics (now available as aes
physics and not echam anymore.

- CSCS: Introduction of more and optimized OpenACC code parts and
        introducing work to allow later more modularization for
        improving and consolidating of the ICON code.

- DKRZ: Bugfixes and more bugfixes, cdi-1.8.4 incorporates nom (joint
        effort of DKRZ and CIMD) cdi-pio, new versions of all tools
        maintained by DKRZ, a new buildbot (Ralf Müller) and strong
        support for getting ICON running on levante.

- DWD: Introduction of many improvements in the parameterizations, in
       tuning and adding bugfixes. Setting up on the new ECMWF system
       and supporting the data assimilation system. First steps for
       ICON-seamless are done, especially the coupling to ICON ocean
       is technically available and most of the time varying external
       data for climate simulations have been added. However, the code
       variant has to be tuned and extensively tested and some parts are
       still missing. This activity is a joint effort of DWD and MPIM
       colleagues.

- MCH: Is actively providing more, and more model components of the
       NWP physics for the use on GPUs. Furthermore they are
       supporting all activities related to GPUs and testing of large
       parts of the ICON code.

- MPIM AES: To simplify and streamline further high resolution
       modelling only, the echam physics has been renamed as aes
       physics and some former components have been removed as are the
       SSO gravity drag parameterization, as well the gravity drag
       parameterization for the middle atmosphere, and most prominent
       the convection parameterization. As cloud paramterizations the
       DWD graupel and two-moment scheme are available. A few more,
       minor, changes have been applied. Newly added is the former
       large eddy model component the Smagorinsky 'vertical' diffusion
       scheme.

- MPIM OES: Added a new C-grid seaice model (C. Mehlmann) and extensive
       improvements of OpenMP for hamocc.

- MPIM LES: JSBACH got more model components and has been adapted for
       to all changes above its top-level.  - KIT: The aerosoly
       microphysics and plume components have been added as well as a
       large number of bugfixes.

- MPIM CIMD: just a lot of things in the area of coupling, configuring
       and building and model infrastrucuture, new machine setups
       for lumi (CSC), juwels/booster (JSC) and levante (DKRZ)

- All developers: many of the work has been done in cross
  institutional working groups.


A special thanks is going to Sergey Kosukhin (MPIM-CIMD) for his work
on the configuration and build system of ICON making life of
developers much easier.

The work on the GPU code version, implemented based on OpenACC has
been supported by NVIDIA. Thanks a lot, Dmitry Alexeev (and his
colleagues).

For detailed information, please, have a look at

https://gitlab.dkrz.de/icon/wiki/-/wikis/Protocol-of-Release-Commits

# Release notes for icon-2.6.4

Another new  release of icon is available.

It consists of

- fix for the soil water budget and snow-tiles
- integrate YAC2, corresponding changes in jsbach and the yaxt version are included
- bundled libraries are no checked at DWD and DKRZ

- fix current CDI-PIO setup for coupled setups
- use plain netcdf intead of CDI for restart IO

- fixes for the building system regarding new compilers, CLAW compatibility, changes in RTTOV and building against eccodes/grib_api

- bugfix for hfbasin diagnostic

- changes for the DWD NEC Aurora

- adjustments for ruby-0 setup (clean up of lsm masks, coupling, add irad_aero=12, buildbot test, scripts)

- New ICON-Land/JSBACH version: jsbach:master@39758f03
  This merge contains many changes, bug fixes, and improvements/new features in ICON-Land/JSBACH.
  The most important new features are:
    Time handling in JSBACH now supports nested domains (with ECHAM physics), plus fix for the time albedo is calculated for radiation
    New formulation for computation of roughness length (ported from JSBACH3) Implementation of anthropogenic land cover change
    Diagnostic 1d global mean JSBACH output for monitoring
    Implementation of the standalone JSBACH model as an ICON model component


- Implementation of the two-moment microphysics scheme by Seifert and Beheng (2006). The original NWP routines are used 'as is'. Use encapsulation.
- Unify variable long names.
- Add new testcases and/or bug fixes: RCEMIP_analytical, RCEtorus, Tconst, bubble.
- Add buildbot tests (S. Rast).
- Make iqneg output distinguishable. Fix typos.
- Add templates for 2 moment in NWP and SP physics.
- Add minmax diagnostics for microphysics inside SP.

- make dom_start_time and dom_end_time relative to experiment start/end time (Hauke Schulz)

- Simple fix of integer range in interface function for CDI-PIO

- update to icon-oes-1.3.05
- optional seaice initialization
- optional calculation of windstress in uncoupled ocean runs
- coupled ruby/dyamond: combine ocean and ice velocities in the coupling interface to the atmosphere

- extended_N-cycle (c8d10814 to 421c584e)

- update of mkexp scripts and mtime

- add perp. month/day

- For coupled configurations, the land and runoff are created using a
  fractional mask which is generated from a selected pair of ocean and
  atmosphere grids and mask. For clarity, these data are stored in a
  directories with names made up of both - atmosphere and ocean - grid
  IDs.

- With the current physics coupling only one set of pressure variables
  is needed. Pressure variables in mo_echam_phy_memory are therefore
  replaced. Use add_ref instead of add_var for several more fields.
  Reintroduce an optimised reading of Kinne aerosols, only those month
  are read in which are need for the current job rather than reading
  always the full year.

- introduce standard emissivity to all test cases
- Fix for Cariolle scheme
- Fix a problem with 10 m wind diagnostics in the coupled model
- In icon-dev.checksuite: Change LC-CTYPE to LC-ALL. With this correction icon-dev.checksuite works again on MacOS Big Sur.
- Bugfix for graupel initialization in Sapphire physics.
- Bug fix for ocean surface coupling in VDIFF:
     Missing A_klev+1/2 term is added, which induce ocean surface momentum stress to the opposite direction.
     Missing alpha term in wind stress diagnostic is added.
     Timing of wind stress diagnostic calculation is fixed.

- new test script exp.atm_ape_mlo_test for mixed layer ocean together with aquaplanet
- adjust experiment test scripts, remove several outdated and unused templates
- changes paths in templates from /pool/data/ICON/grids/private/rene/mpim to /pool/data/ICON/grids/public/mpim
- remove unused and completely outdated icon-authors.txt file
- update simple plume input file for amip reference experiment.

 - New reference data, tolerance ranges and AMIP reference experiment.

- New option for adaptive pressure bias correction at lateral boundaries (limited-area mode only)
- Revised upper boundary condition for vertically nested grids: reduces spurious reflection of vertically propagating sound and gravity waves at the uppermost nest interface level.
- Improved process description for gravity wave emission in SSO scheme (NWP)
- Update of Emvorado radar forward operator, including new options for dual polarization radars
- Update of effective radius coupling with radiation (NWP)
- OpenACC port for LHN code (NWP)
- optional diagnostic of lightning potential index and lightning flash density (NWP)
     new output fields: lfd_con, lfd_con_max, lpi_con, lpi_con_max, mlpi_con, mlpi_con_max, koi
- close output stream if open: fixes rare and random model crashes at simulation end on NEC

- removed hydrostatic DyCore
- removed interface to PSRAD radiation scheme (for NWP)
- removed COSMO ifdefs
- removed unused options for snow-cover fraction diagnosis (for NWP)

# Release notes for icon-2.6.3

The new release of icon is available.

It consists of

- a large number of bugfixes, refactorings, and optimizations for the models infrastructure
- bugfixes and improvements in the build environment
- consolidation of all QUBICC based enhancements and bugfixes (for the sapphire physics)
- many add-ons for the sapphire physics port to GPUs based on OpenACC and many further steps on porting the NWP physics to GPUs
- ART has made its first step to an git submodule external (draft implementation not to be used yet- no warranty)
- fixes and improvements in the ocean code including hamocc
- improvements in the data assimilation NWP physics coupling
- tuning of data assimilation
- added rrtm-gp as radiation scheme for the sapphire physics on GPU
- much progress on cdi-pio use in many icon components
- refactoring of mpi communication library (with a focus on GPU to GPU communication)

For many more details visit:

https://gitlab.dkrz.de/icon/wiki/-/wikis/Protocol-of-Release-Commits

June 7th, 2021
