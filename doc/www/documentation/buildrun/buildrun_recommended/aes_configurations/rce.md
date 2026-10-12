```{eval-rst}
:orphan:
```

(ref_buildrun_rce)=
# Radiative Convective Equilibrium (RCE) Configuration

Input dependencies:
: (for Levante)
```
root = /pool/data/ICON/grids/public/mpim/Torus_Triangles_100000m_x_100000m_res1000m

# Grids
$root/Torus_Triangles_100000m_x_100000m_res1000m.nc

# Initial conditions
$root/rcemip_analytical_o3_100000m_x_100000m_res1000m.nc
```

Compatible machines and compilers:
: Levante CPU (Intel)

Recommended resources:
: Single node.

Estimated runtime (for resources indicated above):
: It takes 5 minutes for 4 simulated hours on one Levante node (2xAMD 7763, 256 Gb main memory). It takes 15 hours for 200 simulated days on two Levante nodes.

Scripting:
: {{ '[RCE mkexp config]({}/run/checksuite.atm/test_aes_rce.config)'.format(base_url) }} (remove the `buildbot` option from `EXP_OPTIONS` and add your slurm account `ACCOUNT` in the mkexp script to run in Levante).

## Description
This model configuration follows the small domain configuration described in [Wing et al., 2018](https://doi.org/10.1029/2020MS002138). On the horizontal, the setup consist of a doubly periodic domain [bi-periodic domain](https://easy.gems.dkrz.de/Processing/playing_with_triangles/toroidal_grids.html) on a 100km by 100km square domain, with a resolution of 1 km. The vertical dimension is discretized by means of 75 levels extending to 33947 meters.

The radiative forcing at the top of the atmosphere is set to 551.58 $W/m^2$, with a fixed zenith angle of 42.05&deg;. Rotation is turned off. The surface temperature is set to 300K. The profiles of trace gases have no spatial dependence. The initial profile of specific humidity is set by means of an analytic formula.
