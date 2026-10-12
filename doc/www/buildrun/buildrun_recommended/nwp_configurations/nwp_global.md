```{eval-rst}
:orphan:
```

(ref_buildrun_nwp_global)=
# Use Case for Global ICON-NWP (R02B06)

Data Download:
: [Click here](https://swift.dkrz.de/v1/dkrz_4d992e1b-f237-4258-a2bc-138ca6a1cf59/icon-model-use-cases/nwp-global-R02B06.tar.bz2) or use

```shell
wget https://swift.dkrz.de/v1/dkrz_4d992e1b-f237-4258-a2bc-138ca6a1cf59/icon-model-use-cases/nwp-global-R02B06.tar.bz2
```

Input files:
: [**Grid files**](ref_buildrun_grids)
  : `nwp-global-R02B06/icon_grid_0057_R02B05_R.nc`: R02B05 grid  (80 km effective mesh size), [reduced radiation grid](ref_atmosphere_ecrad_redgrid)
  : `nwp-global-R02B06/icon_grid_0058_R02B06_G.nc`: R02B06 grid  (40 km effective mesh size), global, domain 1
  : `nwp-global-R02B06/icon_grid_0059_R02B07_N02.nc`: R02B07 grid (20 km effective mesh size), nested domain over Europe, domain 2
: [**External parameter files**](ref_buildrun_external_param)
  : `nwp-global-R02B06/extpar/icon_extpar_0058_R02B06_G_20220825_tiles.nc`: External parameter file for domain 1
  : `nwp-global-R02B06/extpar/icon_extpar_0059_R02B07_N02_20220825_tiles.nc`: External parameter file for domain 2
: [**Initial conditions**](ref_buildrun_icbc)
  : `nwp-global-R02B06/data/igfff00000000_2024020700`: Initial data for February 7, 2024, 00 UTC (for domain 1, domain 2 is initialized internally due to a time-shift in the start time)
: **Dictionaries** `nwp-global-R02B06/mapfiles`: Files to map the ICON internal variable names to I/O names used for ECCODES (GRIB2 shortnames)

Running the test:
: **Run script** `nwp-global-R02B06/run_R02B06N07`: run-script for ECMWF Atos computer in Bologna
: **Log output** `nwp-global-R02B06/run_R02B06.39880836.out`: job-log file of a run on Atos with icon-model-2025.10

Recommended resources:
: 256 cores

Estimated runtime (for resources indicated above):
: about 5 minutes on 256 cores of the Atos Computer of ECMWF in Bologna

## Description

This use case runs a global ICON application with a resolution of about 40 km (R02B06). It includes a nest (R02B07) over Europe.

The run-script `run_R02B06N07` defines the batch-parameters (for SLURM on Atos), sets I/O directories
and all the necessary namelist input for ICON. It can easily be modified to run on different
platforms.

`run_R02B06N07_2025.04` is an older version of the run-script. Some namelist input has been adapted now.
