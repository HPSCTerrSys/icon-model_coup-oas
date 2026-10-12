```{eval-rst}
:orphan:
```

(ref_buildrun_amip)=
# AMIP Configuration

AMIP is a standard experimental protocol for global Atmospheric General Circulation Models (AGCMs) and prescribes realistic Sea Surface Temperature and Sea Ice from 1979 till near present. Virtually the entire international climate modeling community has participated in this project since its inception in 1990. This framework enables a diverse community of scientists to analyze AGCMs in a systematic fashion, a process which serves to facilitate model improvement. It is managed by PCMDI (Earth System Model Evaluation Project). See [AMIP](https://pcmdi.llnl.gov/mips/amip/home/overview.html).

Input dependencies:
: (for Levante)
```
root = /pool/data/ICON/grids/public/mpim/0052  # Any grid spacing works, e.g. 0052, 0054, etc

# Grid
$root/icon_grid_0052_R02B06_G.nc  # icon_grid_0054_R02B08_G.nc if grid spacing is 0054, etc

# Initial conditions
$root/initial_conditions/*

# Boundary conditions
$root/ozone/*
$root/sst_and_ice/*
$root/aerosol_kinne/*
$root/land/*
```

Compatible machines and compilers:
: Levante CPU (Intel,GNU,NAG) and GPU (NVHPC).

Recommended resources:
: Depends of the horizontal resolution and if CPUs or GPUs are used. For CPUs, R02B06 run with 4 nodes, R02B08 with 32 nodes, R02B09 with ~500 nodes and R02B10 with ~600 nodes. For GPUs, R02B04 and R02B06 use one node with nblocks_c = 1 and nblocks_sub=4, R02B08 can run with 8, 16 and more than 32 nodes. However, the configuration of the nblocks_c and nblocks_sub differs. With 8 nodes, nblocks_c = 2 and nblocks_sub = 5. With 16 nodes, nblocks_c = 1 and nblocks_sub = 4. With 32 nodes and more, nblocks_c = 1 nblocks_sub = 1. R02B09 run with 24 nodes with nblocks_c = 1 and nproma_sub = 800.

Estimated runtime  (for resources indicated above):
: It takes 0.8 node hour for one hour simulation in R02B08. A similar performance is expected for GPUs.

Scripting:
: {{ '[AMIP mkexp config for CPU]({}/run/examples/amip-aes.config)'.format(base_url) }}, {{ '[AMIP mkexp config for GPU]({}/run/examples/amip-aes_gpu.config)'.format(base_url) }} (add your slurm account `ACCOUNT` in the mkexp script to run in Levante).

## Description
ICON can run in an AMIP configuration. This means solving the fluid dynamics equations in the atmopshere with three main paramaterization (microphysics, turbulence and radiation) on the entire globe. In this configuration, the atmosphere is coupled to a dynamical 1D land module. This means that land does not transport energy nor water horizontally. The ocean and the sea-ice are not dynamically active (uncoupled). In other words, sea surface temperature and sea-ice extension are prescribed.

In the default configuration, the atmosphere is divided in 90 vertical levels and the land uses 5 soil layers. The configuration of the grid is global and it can be used with horizontal grid spacing of 160 km and finer.

The default configuration has the initial conditions of January 1st 2020. Other initial conditions can be found `/pool/data/ICON/grids/public/mpim/XXXX/initial_conditions/r0100`. Two dates are available to be used as initial conditions across the different grid configurations: 1979-01-01 and 2020-01-01. But certain grid configurations have more available dates, e.g., R02B08 has four: 1979-01-01, 1990-01-01, 2020-01-01, and 2020-08-01 (`/pool/data/ICON/grids/public/mpim/0054/initial_conditions/r0100`). If another date want to be used, use the script `/pool/data/ICON/grids/public/mpim/XXX/initial_conditions/r0100/11u-make-initial-data-ifs2icon-from-era5.sh` (only one script) to generate a new initial conditions from a pool of files in `/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/ecmwf_initial_data/initial_conditions)`.

Aside from the initial conditions, the AMIP configuration needs time-dependent input data from aerosols, ozone, land properties, sea surface temperature, and sea ice extent. These inputs are automatically updated with the simulation time.
