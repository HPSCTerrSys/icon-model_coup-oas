```{eval-rst}
:orphan:
```

(ref_waves_overview)=
# Waves

🌊 The ocean surface gravity wave model ICON-waves (Dobrynin et al., GMD, in preparation) is designed to explicitly model ocean surface gravity waves and their feedbacks on the atmosphere and ocean within the Earth system. ICON-waves is a joint effort led by the Deutscher Wetterdienst (DWD) with contributions from the Max Planck Institute for Meteorology (MPI-M), the German Climate Computing Center (DKRZ) and the Helmholz research center HEREON. Surface waves influence the sea surface state, generate turbulence, modify ocean currents, and affect air-sea exchanges of heat, matter, and momentum. ICON-waves addresses these processes by providing a wave-spectrum-dependent interface within the ICON framework. ICON-waves introduces two-way coupling through the coupler YAC (Yet Another Coupler), providing dynamic feedback of wave-induced processes to both the atmosphere and ocean (under development). The integration of ICON-waves into the ICON framework represents a significant advancement in modeling the complexity of atmosphere-ocean feedbacks, enabling more realistic simulations of atmosphere-ocean interactions and offers potential benefits for weather and climate prediction.

(ref_waves_model)=
## Spectral approach of modeling of surface waves
ICON-waves solves the spectral wave energy balance equation in geographical dimensions with latitude (𝜙) and longitude (𝜆), and spectral dimensions directions (𝜃) and frequencies (ω). The model computes the evolution of the wave spectrum 𝐹(𝜔,𝜃,𝜙,𝜆,𝑡), simulating the energy propogation distributed across directions and frequencies over space and time. The model employs advanced numerical techniques, such as the Flux Form Semi-Lagrangian (FFSL) finite-volume scheme for spatial transport, which is balanced by the source function, describing the wave physics. The source function includes the wind energy input, energy dissipation (including whitecapping, bottom friction, and depth-induced breaking), and nonlinear wave-wave interactions terms. The physical parameterizations of the source function is adapted from the well-established WAM model ({term}`The Wamdi Group 1988`; {term}`Komen et al 1996`) and is fully integrated into the ICON code infrastructure, allowing for seamless operation on ICON's flexible icosahedral-triangular grid. Temporal integration uses a splitting method, treating advection and source terms separately. The source term integration allows full implicit schemes for stability.

(ref_waves_config)=
## Configuration and simulation modes
Conceptually, the ICON-waves component is an integral part of the ICON framework, similar to the ICON Atmosphere (ICON-A) and Ocean (ICON-O) components. In general, within the ICON framework, a specific component can be run by combining configuration switches and setting the model type and model name in the master model namelist file. A set of provided configuration wrappers, located in the **config** folder for different computational systems, can be used. For example, for DKRZ's Levante HPC, one of the configuration wrappers located in **config/dkrz**, along with the additional command line argument `--enable-waves`, will configure the ICON system to run the ICON-waves component as well (assuming the starting configuration is from the user-created **build** folder):

```bash
../config/dkrz/levante.intel --enable-waves
```

For a specific experiment with ICON-waves, the settings in the master model namelist must be:

```fortran
&master_model_nml
 model_type                  = 98
 model_name                  = "wave"
 ...
```

ICON-wave can be run in two modes: standalone and coupled. In both modes, the meridional and zonal components of the 10-meter wind speed, as well as ice concentration fields, are required (except for test case experiment, see below). In standalone mode, the forcing is provided through external files. In coupled mode with ICON-A, the forcing is exchanged via the YAC coupler. Following fields are exchanging in the coupled ICON-A - ICON-waves mode:

| send from ICON-A to ICON-waves      | send from ICON-waves to ICON-A |
|-------------------------------------|--------------------------------|
| Zonal wind at 10 m                  | Sea surface roughness length   |
| Meridional wind at 10 m             |                                |
| Fraction of ocean covered by sea ice|                                |
---

ICON-waves test numerical experiments can be run using the provided experiment templates located in the **run/checksuite.nwp** folder. Following experements with ICON-waves are available:

| Experiment                                      | Description                   |
|-------------------------------------------------|-----------------------------------------------------|
| nwpexp.run_ICON_18_R2B6_waves                   | Test case simulation on grid R2B6 (no external files requeried)|
| nwpexp.run_ICON_18_R2B4_waves_adv_nophys        | Waves with advection, no physics, on grid R2B4      |
| nwpexp.run_ICON_21_R2B4_waves_standalone_restart| Standalone wave run with restart on grid R2B4       |
| nwpexp.run_ICON_23_R2B4_atmo_waves_coupled      | Coupled atmosphere-waves run on grid R2B4          |
| nwpexp.run_ICON_28_R2B6_IAU_atmo_waves_coupled  | Coupled atmosphere-waves run with IAU on grid R2B6 |

 Run scripts for different experiments can be generated using `make_runscripts`, for example on DKRZ's Levante HPC, assuming the starting point is the user-created **build/run/checksuite.nwp** folder for the experiment `nwpexp.run_ICON_18_R2B6_waves`:
```bash
../../make_runscripts run_ICON_18_R2B6_waves -r run/checksuite.nwp
```
This creates a run script **nwpexp.run_ICON_18_R2B6_waves.run**, which can be submitted by:
```bash
sbatch nwpexp.run_ICON_18_R2B6_waves.run
```
The experiment output will be stored in **build/experiments/run_ICON_18_R2B6_waves**

(ref_waves_output)=
## Output parameters
ICON-waves supports the output of numerous prognostic and diagnostic parameters, including the full wave spectrum. The most commonly used outputs are the spectrum-integrated parameters, such as

| Name           | Units  | Description                          |
|----------------|--------|--------------------------------------|
| hs             | m      | Total significant wave height        |
| hs_dir         | deg    | Total mean wave direction            |
| tpp            | s      | Total wave peak period               |
| tmp            | s      | Total wave mean period               |
| tm1            | s      | Total m1 wave period                 |
| tm2            | s      | Total m2 wave period                 |
| ds             | deg    | Total directional wave spread        |
| hs_sea         | m      | Sea significant wave height          |
| hs_sea_dir     | deg    | Sea mean wave direction              |
| pp_sea         | s      | Sea wave peak period                 |
| mp_sea         | s      | Sea wave mean period                 |
| m1_sea         | s      | Sea m1 wave period                   |
| m2_sea         | s      | Sea m2 wave period                   |
| ds_sea         | deg    | Sea directional wave spread          |
| hs_swell       | m      | Swell significant wave height        |
| hs_swell_dir   | deg    | Swell mean wave direction            |
| pp_swell       | s      | Swell wave peak period               |
| mp_swell       | s      | Swell wave mean period               |
| m1_swell       | s      | Swell m1 wave period                 |
| m2_swell       | s      | Swell m2 wave period                 |
| ds_swell       | deg    | Swell directional wave spread        |
|||**Output of wave spectrum**|
| tracer_xxx     | m²Hz⁻¹ | Spectral bin of wave energy*         |

*xxx, for example 001, is the frequency index, where the index ranges from 1 to the length of tracer_xxx, representing each direction of the wave spectrum.
