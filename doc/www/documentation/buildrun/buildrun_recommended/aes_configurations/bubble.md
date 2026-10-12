```{eval-rst}
:orphan:
```

(ref_buildrun_bubble)=
# Bubble Configuration

Input dependencies:
: (for Levante)
```
root = /pool/data/ICON/grids/public/mpim

# Grid (any torus grid works)
$root/Torus_Triangles_20x4_5000m/Torus_Triangles_20x4_5000m.nc  # (Example)
```

Compatible machines and compilers:
: Levante CPU, MacOS

Recommended resources:
: Single Processor.

Estimated runtime (for resources indicated above):
: Tens of seconds

Scripting:
: {{ '[Bubble mkexp config]({}/run/examples/bubble.config)'.format(base_url) }} (note the initialization is performed {{ '[here]({}/src/testcases/mo_aes_bubble.f90)'.format(base_url) }}, using {{ '[default parameter settings]({}/src/configure_model/mo_aes_bubble_config.f90)'.format(base_url) }}, which, along with model parameters, can be modified within Namelists) (add your slurm account `ACCOUNT` in the mkexp script to run in Levante).

Analysis/postprocessing:
: Checkout the examples on [easyGEMS](https://easy.gems.dkrz.de/simulations/ICON/bubble.html).


## Description

The basic set-up of the bubble experiment consists of an atmosphere over a [bi-periodic domain](https://easy.gems.dkrz.de/Processing/playing_with_triangles/toroidal_grids.html). A roughly 20 km deep atmosphere is discretized over 70 layers and initialized to form two uniformally stratified layers upon which a perturbation is imposed. The first layer extending upwards with a temperature profile given by

```{math}
T(z) = T_0 + z\Gamma_0
```

to a height z_0 above which the temperature follows the profile

```{math}
T(z) = T(z_0) + (z-z_0)\Gamma_1
```

where T_0, z_0, Gamma_0, and Gamma_1 are parameters that can be specified. The humidity is set to a constant prescribed value. The perturbation, or "bubble", is specified as a latitudinally slab symmetric gaussian perturbation in relative humidity and temperature that is centered at the surface and the central longitude. The surface is given the properties of saturated water at a fixed temperature, and the mean wind is initialized to 0 m/s.
In terms of the physical set up, it is often most useful to vary T_0 (`aes_bubble_config%t0` which is set in `aes_bubble_nml`) and the surface temperature (`ape_sst_val` which is set in `nh_testcase_nml`) as this determines the phase of the condensate that forms in the bubble.
