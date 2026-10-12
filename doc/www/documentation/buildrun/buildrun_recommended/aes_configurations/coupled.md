```{eval-rst}
:orphan:
```

(ref_buildrun_aes_coupled)=
# Coupled AES Configuration

Input dependencies:
: (for Levante)
```
root = /pool/data/ICON/grids/public/mpim
model_dir = <checkout/working-directory>

# Grids
$root/0036/icon_grid_0036_R02B04_O.nc
$root/0043-0036/land/r0002/hdpara_icon_r2b4_013_using_grid_lsmask.nc
$root/0043/icon_grid_0043_R02B04_G.nc

# Restart files
$root/0036/ocean/restart/r0001/ler1166_restart_oce_21000101T000000Z.nc
$root/0043-0036/land/r0002/hdrestart_R02B04_013_G_210120_1334_with_grid_file_mask.nc

# Initial conditions
$root/0043-0036/land/r0002/ic_land_soil.nc
$root/0043/initial_condition/r0001/ifs2icon_1979010100_R02B04_G.nc

# Boundary conditions
$root/0043-0036/land/r0002/bc_land_frac.nc
$root/0043-0036/land/r0002/bc_land_phys.nc
$root/0043-0036/land/r0002/bc_land_soil.nc
$root/0043-0036/land/r0002/bc_land_sso.nc
$root/0043/aerosol_kinne/r0001/bc_aeropt_kinne_lw_b16_coa.nc
$root/0043/aerosol_kinne/r0001/bc_aeropt_kinne_sw_b14_coa.nc
$root/0043/aerosol_kinne/r0001/bc_aeropt_kinne_sw_b14_fin_2014.nc
$root/0043/ozone/r0001/bc_ozone_historical_2014.nc
$root/0043/sst_and_seaice/r0001/bc_sic_1979_2016.nc
$root/0043/sst_and_seaice/r0001/bc_sst_1979_2016.nc
$root/independent/greenhouse_gases/greenhouse_ssp245.nc
$root/independent/solar_radiation/3.2/swflux_14band_cmip6_1850-2299-v3.2.nc

# Model parameters
$model_dir/data/ECHAM6_CldOptProps_rrtmgp_lw.nc
$model_dir/data/ECHAM6_CldOptProps_rrtmgp_sw.nc
$model_dir/data/rrtmgp-gas-lw-g128.nc
$model_dir/data/rrtmgp-gas-sw-g112.nc
$model_dir/externals/jsbach/data/lctlib_nlct21.def
```

Compatible machines and compilers:
:   Levante CPU (add `--enable-openmp` for optimized performance)

Recommended resources:
: 4 CPU nodes

Estimated runtime (for resources indicated above):
: _The simplified test configuration (see [Description](#description)) is numerically stable only for a couple of days. It runs 6 model hours in about 40 s_

Scripting:
:   &nbsp;
    - {{ '[Coupled AES mkexp config]({}/run/mkexp/types/AES/coupled-R02B04L90.config)'.format(base_url) }} (add your slurm account `ACCOUNT` in the mkexp script to run in Levante)

Analysis/postprocessing:
: For all models participating in the [nextGEMS](https://nextgems-h2020.eu/) project, examples, software and recipes for data processing have been collected on the [easyGEMS](https://easy.gems.dkrz.de/) site. See the ICON specific entries there for more information.


## Description

Unified coupled atmosphere/ocean configuration based on the AES atmosphere package as used for the [nextGEMS](https://nextgems-h2020.eu/) production experiments.
The AES package features a reduced set of sub-grid scale parametrizations due to the targeted km-scale model resolutions.
At the time of this writing, only a simple test model is implemented, using 160 km (R2B4) horizontal resolution with a comparably short timestep to compensate instabilities due to the reduced parametrization set. For scientific production, definitions for horizontal resolutions of at least 10 km atmospheric and 5 km oceanic (R2B8_R2B9) are currently under development.
