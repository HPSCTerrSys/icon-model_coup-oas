```{eval-rst}
:orphan:
```

(ref_buildrun_aquaplanet)=
# Aquaplanet Configuration

Input dependencies:
: (for Levante)
```
root = /pool/data/ICON/grids/public/mpim/

# Grid (any grid spacing works, e.g. 0052, 0054, etc )
$root/0052/icon_grid_0052_R02B06_G.nc  # icon_grid_0054_R02B08_G.nc if grid spacing is 0054, etc

# Boundary conditions
$root/0052/ozone/*  # Any
$root/common/greenhouse_historical_plus.nc
```

Compatible machines and compilers:
: Levante

Recommended resources:
: Single Processor for development, for scientific purpose it should be four nodes for R02B04.

Estimated runtime (for resources indicated above):
: For R02B04, $\approx$ 640 SDPD for a single node, and $\approx$ 2500 SDPD for four nodes. For development, it is only required a few time steps; however, for scientific purposes, it is required six month of spin-up time and a minimum of four months and ideally one year of simulated time.

Scripting:
: {{ '[Aquaplanet mkexp config]({}/run/checksuite.atm/test_aes_ape.config)'.format(base_url) }} (the aquaplanet is configured and set in {{ '[here]({}/src/testcases/mo_nh_testcases.f90'.format(base_url) }}, using the time constant sea surface temperature {{ '[here]({}/src/testcases/mo_ape_params.f90)'.format(base_url) }}, which, along with model parameters, can be modified within Namelists) (remove the `buildbot` option from `EXP_OPTIONS` and add your slurm account `ACCOUNT` in the mkexp script to run in Levante).

Analysis/postprocessing:
: (under development)


## Description

The aquaplanet experiment follows the [APE protocol](https://doi.org/10.1006/asle.2000.0022). This configuration tests the atmospheric component of ICON and its interaction with an idealized ocean, i.e. surface fluxed modeled from a prescribed constant zonally symmetric sea surface temperature.

In more detail, an aquaplanet experiment consists of an Earth-sized planet with an Earth-like atmosphere whose lower boundary condition is consistent with a water-covered surface (no sea-ice) with a prescribed and zonally and temporally constant sea surface temperature. The [APE protocol](https://doi.org/10.1006/asle.2000.0022) proposed different sea surface temperature profiles; however, we use the _SST-QOBS_ configuration  since it resembles the zonal and time average of Earth's surface temperature. The surface temperature peaks at the equator with a temperature of 27 °C and drops to zero at the 60° latitude. In addition, it maintains a perpetual equinox with symmetric-constant radiation about the equator by setting the eccentricity and obliquity of the Earth to zero, and using a solar constant of 1361 W/m$\rm{^2}$. Ozone follows a constant zonally symmetric profile about the equator, the greenhouse gases are well mixed, and the interaction between aerosol and radiation is neglected. By using this configuration, the aquaplanet's forcing is symmetrical with respect to the equator, and hence its statistical climate state is symmetrical for long time integration.

The vertical grid setup of the aquaplanet experiment consist of a stretched vertical grid with 90 levels, where levels are more finely spaced close to the surface than at the model top at 75km, where a damping is also employed, increasing in strength from 44km upwards.
