```{eval-rst}
:orphan:
```

(ref_sbm_overview)=
# SBM Overview

The Spectral Bin Microphysics (SBM) parametrization is a detailed cloud microphysics scheme used in atmospheric models that explicitly represents the size (or mass) distribution of cloud particles by dividing them into discrete size (or mass) bins. Unlike bulk microphysics schemes that use only a few prognostic variables to represent entire particle populations, SBM tracks the number concentration of particles in each size bin, allowing for a more realistic representation of microphysical processes such as condensation, evaporation, collision-coalescence, and ice processes. This approach provides higher fidelity in simulating hydrometeors size distributions, which is crucial for accurately predicting precipitation formation, cloud optical properties, and aerosol-cloud interactions. However, SBM schemes are computationally expensive due to the large number of prognostic variables required (typically tens of size bins for different hydrometeor types), making them primarily suitable for high-resolution simulations or research applications where detailed microphysical representation is essential for understanding cloud processes and their climate impacts.
The SBM scheme was implemented in ICON following a similar implementation in WRF ({term}`Shpund et al. 2019`). It is a "light" version (so called Fast-SBM) of the comprehensive SBM scheme described in {term}`Khain et al. 2004` and {term}`Khain et al. 2008`. A partial "warm-phase" SBM version implemented in ICON several years ago is described in {term}`Khain et al. 2022`.

This page aims to describe the configuration of SBM implemented in ICON.

(ref_sbm_implementation)=
# SBM implementation in ICON

SBM scheme is activated using the namelist parameter **{term}`itype_gscp`=8**.
In contrast to the bulk schemes, the SBM scheme does not use the saturation adjustment procedure before and after the call to the microphysics scheme. To solve the equation of diffusional growth/evaporation SBM uses the supersaturation change during the last time step, which occurs due to advection with contribution from other physics processes. The coupled equation system for diffusional growth and supersaturation is solved, so that supersaturations over water and ice are non-zero to the end of microphysical time sub-steps.  Therefore, it is required to turn off saturation adjustment, i.e. set the namelist parameter **{term}`inwp_satad`=0**. In order to calculate the supersaturation change during the time step, a special "sbm storage" (see mo_sbm_storage.f90) was implemented, which allows to save and restore the temperature and water vapor of the previous time step.
The main "fast-physics" parametrizations called in ICON every time step are the turbulent mixing "turbdiff" and microphysics. It is reasonable for SBM to use the supersaturation right after the advection, rather than the smoothed supersaturation after "turbdiff". Therefore it is required to set the corresponding namelist parameter **{term}`lmicrophysicsFirst`=.true.**.
In order to fully couple SBM to ICON's dynamics, it is required to set the corresponding namelist parameter **{term}`lsbm_coupled`=.true.**. For testing purposes, there is an option to run ICON with 2-moment dynamics and microphysics, and on top of it, every time step run SBM with an option to analyze its output (e.g. using extra_3d variables). To choose this option, where SBM does not have any influence on the model dynamics or microphysics, one needs to select the namelist parameter **{term}`lsbm_coupled`=.false.**.
In addition to already existing tracers in ICON, SBM includes 132 advected tracers. 3 hydromoteors size distributions for drops, ice/snow and graupel/hail are defined on logarithmic doubling mass grid containing 33 bins (tracers). The minimum particles radius is 2 µm, maximum radius is 3.4 mm. Drops with radii exceeding 50 µm (bin number krdrop=15) are attributed to raindrops, while smaller drops are considered as cloud droplets. Similarly, ice/snow particles with bin number exceeding krice=18 are attributed to snow flakes and smaller particles are considered as ice crystals. In contrast, graupel/hail 33 bins are defined as pure graupel or pure hail distributions by hail_opt parameter.
Size distributions of dry CCN are defined on a logarithmic mass grid also containing 33 bins (tracers). The maximum CCN radius is taken equal to 2 µm. Initial size distribution of aerosols playing the role of cloud condensational nuclei (CCN) applied in SBM are described as a sum of three log-normal distributions, representing nuclei mode (ultrafine aerosols) centered at radius of 0.005 µm, accumulation mode centered at 0.035 µm and coarse mode centered at 0.31 µm ({term}`Ghan et al. 2011`). Supersaturation and CCN properties are used to calculate the critical size using the Kohler theory. The CCN with sizes exceeding the critical value are converted to droplets and corresponding bins in the aerosol size distribution become empty. Non-activated aerosols are advected similarly to hydrometeors, so wind transports both hydrometeors and non-activated CCN. Settling velocity of CCN is neglected.
There are several options in SBM to define an initial vertical profile of CCN. Among others, the code includes the definitions of theoretical maritime and continental options. By default, the SBM uses a continental initial profile. A namelist scaling factor **{term}`tune_sbmccn`** [0-1] allows to reduce the CCN concentration initial profile with respect to the continental case. **{term}`tune_sbmccn`=1** means the use of pure continental profile and **{term}`tune_sbmccn`=0** is an approximation of maritime case.

## Interaction with other parametrizations

ICON uses total water and ice contents (qc, qr, qi, qs, qg) in many places in the code. When using SBM, we keep advecting these total water contents, but overwrite them by the integrals over the corresponding size distributions, right after the call to the microphysics. The vertical turbulent diffusion "turbdiff" also uses cloud water, ice and snow (qc, qi, qs). For consistency, the corresponding SBM bins (tracers) were added to the vertical turbulent diffusion as well. Convection scheme might also detrain qc, qi, qr, qs. For mass conservation, we "evaporate" this water, transfering it to water wapor qv, and subtract the corresponding latent heat release. In appropreate conditions, this water vapor may condensate to hydrometeors as part of SBM. Total water contents are used in ICON initial and boundary conditions. In case of SBM, similar transformation to water vapor is applied.

## Description of SBM microphysical processes used in ICON

The equation for diffusional growth/evaporation is solved together with the equations for supersaturation, using microphysical sub-steps, which are three times smaller than the dynamical time step, so that the local and advective tendencies of supersaturation during microphysical time steps are taken into account. In this way, the water vapor condensates on the growing droplets and the supersaturation is reduced. In strong updrafts, where the supersaturation growth rapidly or in cases of low drop concentration, this reduction does not reach zero, leading to lower latent heat release than with saturation adjustment used in bulk schemes.
The hydrometeors size distributions evolution due to collisions between drops is calculated by solving the stochastic equation for collisions using the low diffusional method by {term}`Bott 1998`. Collision kernels depend on height as described in {term}`Khain & Pinsky 2018`. Probabilities for spontaneous breakup were taken from laboratory measurements ({term}`Kamra et al. 1991`). The changes of size distribution function by the collisional breakup of drops are described as in {term}`Seifert et. al. 2005`. Drop sedimentation is performed separately for each mass bin, so new droplet size distribution forms automatically after diffusion growth, collisions and sedimentation.
Ice particles are formed either by drop freezing or primary ice nucleation. Heterogeneous drop freezing is calculated using the method close to that of Bigg ({term}`Khain et al. 2004`), where probability of freezing increases with drop mass and with the decrease of temperature. Freezing of small droplets with radii lower than 100 um is assumed to lead to formation of ice crystals (smallest bins in the snow size distribution). Freezing of larger drops leads to formation of graupel or hail (in the hail version). Homogeneous drop freezing takes place at temperatures lower than -39C. Droplet nucleation at lower temperatures lead to their immediate freezing.
Rate of primary ice nucleation calculated using the parameterization of {term}`Meyers et al. 1992`, depends on supersaturation over ice and temperature.
Ice-ice and ice-water collisions are calculated using a modified Bott method. Snow-snow collisions lead to formation of larger snow particles. Water-graupel/ hail collisions (riming process) lead to the formation of larger graupel or hail, respectively. Water-snow collisions lead to formation of graupel (or hail) in case the mass of liquid drop is larger than that of ice particle. Collisions of large snow with small drops lead either to the formation of larger snow if the liquid water content in the environment is less than a critical value, or to the formation of graupel/hail otherwise. Fall velocities depend on type of ice particles, their size and altitude (air density) and are taken from empirical relationships ({term}`Khain & Pinsky 2018`). Properties of ice particles (density, aspect ratios, fall velocities) are expressed via their mass. It means that the changes of ice properties by deposition/sublimation) are automatically taken into account. Dependencies of particle shape (aspect rations) of water drops and ice particles allow one to use the output for calculation of polarimetric radar signatures. A breakup of the largest snowflakes is calculated with probability increasing with snowflake size.
The decrease in the masses of bins in ice particles by melting and the increase in mass of the corresponding drop bins by melting is calculated by prescribing the characteristic melting time scales, which depend on the ice particles type and bin number (i.e. particle size). The smallest melting time scale is assumed for snow and the largest time scale is assumed for hail. No shading is assumed in this SBM version.

## Limitations:

SBM is not supported on GPU (no option to compile with OPENACC).


# Glossary of Namelist Parameters

:::{glossary}
itype_gscp
  (`$nwp_phy_nml`) Choose the microphyics scheme. Set **8** to choose SBM

inwp_satad
  (`$nwp_phy_nml`) Turn on/off saturation adjustment. Set **0** when using SBM

lmicrophysicsFirst
  (`$nwp_phy_nml`) TRUE: run microphysics before turbdiff, FALSE: after turbdiff. Set **TRUE** when using SBM

lsbm_coupled
  (`$nwp_phy_nml`) FALSE: use two-moment scheme for feedback and run uncoupled SBM, TRUE: use SBM feedback. Recommended to set **TRUE** when using SBM

tune_sbmccn
  (`nwp_tuning_nml`) [0-1] scaling factor to reduce the ccn concentration initial profile with respect to the continental (polluted) case when using SBM. Recommended to set **1** for polluted case and **0.1** for maritime (clean) case.
:::
