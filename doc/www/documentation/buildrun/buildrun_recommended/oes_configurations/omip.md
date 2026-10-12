```{eval-rst}
:orphan:
```

(ref_buildrun_omip)=
# OMIP Configuration

Input dependencies:
: (for Levante)
```
root = /pool/data/ICON/oes/input/r0002/OceanOnly_Icos_0158km_etopo40

# R2B4 L40 grid
$root/OceanOnly_Icos_0158km_etopo40.nc

# OMIP forcing (linked to `ocean-flux.nc`)
$root/omipForcing-mpiomDaily-OceanOnly_Icos_0158km_etopo40.nc

# PHC initial conditions (linked to `initial_state.nc`)
$root/omipInitialState-AnnualAverage-OceanOnly_Icos_0158km_etopo40-40levels.nc

# PHC salinity restoring (linked to `ocean-relax.nc`)
$root/omipRelaxSurface-OceanOnly_Icos_0158km_etopo40.nc
```

Compatible machines and compilers:
: Levante CPU (nag, gfortran, ifort), and GPU (nvhpc)

Recommended resources:
: Single node

Estimated runtime (for resources indicated above):
: few minutes

Scripting:
: {{ '[OMIP mkexp config]({}/run/checksuite.ocean_internal/test_oes_omip.config)'.format(base_url) }} (remove the `buildbot` option from `EXP_OPTIONS` and add your slurm account `ACCOUNT` in the mkexp script to run in Levante).

Analysis/postprocessing:
: (under development).


## Description

This is a coarse resolution configuration (r2b4 l40 which corresponds to 160km horizontal resolution) of the ICON ocean model component including sea ice. The configuration is driven with OMIP (Ocean Model Intercomparison Project) forcing data and uses a realistic bathymetry. It applies the `zstar` vertical coordinate.
