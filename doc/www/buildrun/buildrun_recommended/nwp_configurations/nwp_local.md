```{eval-rst}
:orphan:
```

(ref_buildrun_nwp_local)=
# Use Case for Limited-Area ICON-NWP (R19B07)

Data Download:
: [Click here](https://swift.dkrz.de/v1/dkrz_4d992e1b-f237-4258-a2bc-138ca6a1cf59/icon-model-use-cases/nwp-local-R19B07.tar.bz2) or use

```shell
wget https://swift.dkrz.de/v1/dkrz_4d992e1b-f237-4258-a2bc-138ca6a1cf59/icon-model-use-cases/nwp-local-R19B07.tar.bz2
```

Input files:
: [**Grid files**](ref_buildrun_grids)
  : `nwp-local-R19B07/grids/icon_grid_0046_R19B06_LR.nc`: R19B06 grid  (4 km effective mesh size), [reduced radiation grid](ref_atmosphere_ecrad_redgrid)
  : `nwp-local-R19B07/grids/icon_grid_0047_R19B07_L.nc`: R19B07 grid  (2 km effective mesh size), limited-area domain over Germany
  : `nwp-local-R19B07/grids/icon_grid_0047_R19B07_L_lbc.nc`: R19B07 grid (2 km effective mesh size), frame grid for sparse lateral boundary reading
: [**External parameter file**](ref_buildrun_external_param)
  : `nwp-local-R19B07/extpar/icon_extpar_0047_R19B07_L_20220601_tiles.nc`: External parameter file
: [**Initial conditions**](ref_buildrun_icbc)
  : `nwp-local-R19B07/data/igfff00000000`: Initial data for August 31, 2023, 00 UTC
: [**Lateral boundary conditions**](ref_buildrun_icbc)
  : `nwp-local-R19B07/data/*_lbc`: Lateral boundary conditions (3 hourly) for up to +24h
: **Dictionaries** `nwp-local-R19B07/mapfiles`: Files to map the ICON internal variable names to I/O names used for ECCODES (GRIB2 shortnames)

Running the test:
: **Run script** `nwp-local-R19B07/run_icon_d2`: run-script for ECMWF Atos computer in Bologna
: **Log output** `nwp-local-R19B07/run_R19B07.39865769.out`: job-log file of a run on Atos with icon-model-2025.10

Recommended resources:
: 1024 cores

Estimated runtime (for resources indicated above):
: 24h forecast about 42 minutes on 1024 cores of the Atos Computer of ECMWF in Bologna

## Description

This use case runs a limited-area (local) ICON application over Germany with a resolution of about 2 km (R19B07). This use case is comparable to DWD's operational application ICON-D2.

The run-script `run_icon_d2` defines the batch-parameters (for SLURM on Atos), sets I/O directories
and all the necessary namelist input for ICON. It can easily be modified to run on different
platforms.

The file `run_R19B07.39865769.out` contains the job log file and ASCII output from the ICON run.

`run_icon_d2_2025.04` is an older version of the run-script. Some namelist input has been adapted now.
