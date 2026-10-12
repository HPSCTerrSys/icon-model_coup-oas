```{eval-rst}
:orphan:
```

(ref_buildrun_shallow_water)=
# Shallow Water Configuration

Input dependencies:
: (for Levante)
```
root = /pool/data/ICON/oes/grids/AquaPlanets/

# Shallow water grid
$root/AquaPlanet_IcosDual_0158km_springOpt.nc
```

Compatible machines and compilers:
: Levante CPU (nag, gfortran, ifort), and GPU (nvhpc)

Recommended resources:
: Single node

Estimated runtime (for resources indicated above):
: few minutes

Scripting:
: {{ '[Shallow Water mkexp config]({}/run/checksuite.ocean_internal/test_ocean_WilliamsonTestCase2_Hex.config)'.format(base_url) }} (remove the `buildbot` option from `EXP_OPTIONS` and add your slurm account `ACCOUNT` in the mkexp script to run in Levante).

Analysis/postprocessing:
: Results should be compared to the initial conditions.


## Description

This is a coarse ICON ocean shallow water configuration based on Williamson et al., 1992 (Test 2,
Global Steady State Nonlinear Zonal Geostrophic Flow). There is only one vertical layer,
no surface forcing and a flat bottom. The dynamics results from the specified initial conditions.
This configuration tests the capability of ICON to run with shallow water equations.
