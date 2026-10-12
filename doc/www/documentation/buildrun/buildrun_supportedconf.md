```{eval-rst}
:orphan:
```

(ref_buildrun_supportedconf)=
# Supported Configurations

ICON partners provide a number of supported scientific configurations, which already provide everything necessary for running ICON. Please note that there is no support for configurations that are not listed below.


## AES Configurations

These configurations are supported by the Max Planck Institute for Meteorology. AES (Atmosphere in the Earth System) indicates the use of AES physics, which is designed for high-performance computing platforms and simulations with horizontal grids of 10 km or less on a variety of domains. The AES configurations can be generated using {{ '[mkexp]({}/doc/www/buildrun/buildrun_running.md)'.format(base_url) }}.

[Aquaplanet](ref_buildrun_aquaplanet)
: This configuration runs the ICON atomspheric component in interaction with an idealized ocean (e.g., SST-QOBS). This set-up allows a fast integration to achieve a climate mean state in comparisson to AMIP simulations.

[Bubble](ref_buildrun_bubble):
: This is a simple idealized case used for code development.  It simulates a buoyant, slab symmetric, bubble, on a small bi-periodic domain. The short simulation time (120 min, 240 timesteps) and the small domain (160 cells) allows output to be written at every timestep and grid point. The slab-symmetric setupt makes it easy to visualize.

[Nest](ref_buildrun_nest):
: This case tests ICON's capability of handling several meshes, nested within one another. The configuration uses a global R2B4 parent mesh, a first-level nested domain that covers the Atlantic ocean, and two second-level nested domains focused over the ITCZ and stratocumulus southern region. The latter two domains share the same parent mesh, which slightly complexifies the case. This allows to test both the parent-child nesting and most features pertaining to limited-area runs.

[AMIP](ref_buildrun_amip):
: This configuration follows the Atmospheric Model Intercomparison Project (AMIP) protocol. In this configuration, ICON solves the fluid dynamics equations in the atmosphere on the entire Globe using horizontal grid spacing of 40 km and finer. The atmosphere is coupled to a 1-D land module and to a non-dynamical ocean and sea ice. This means that sea surface temperature and sea ice area is prescribed. The atmosphere is vertically discretized in 90 levels, and the land is represented by 5 soil layers.

[RCE](ref_buildrun_rce):
: This configuration is intended to study the Radiative Convective Equilibrium (RCE) achieved in ICON. RCE is an important conceptualization for understanding climate. RCE inquires about the thermodynamical equilibrium that a moist atmosphere would attain under constant incoming solar radiation, and is the result of a balance between the latent heating by condensation and the long-wave radiative cooling of the atmosphere. Simulated times in the order of 200 days give insight into how RCE is achieved in ICON. For the purposes of code development this case is run for 6 simulated hours.


### Coupled configurations

[Coupled AES](ref_buildrun_aes_coupled):
: Unified coupled atmosphere/ocean configuration based on a simplified version of the [nextGEMS](https://nextgems-h2020.eu/) experiments. The full version is currenty under development.


## OES Configurations

The OES (Ocean in the Earth System) configurations are supported by the Max Planck Institute for Meteorology.

[OMIP](ref_buildrun_omip):
: This is a close-to-production configuration for ICON's ocean model componed used for code development. In this ocean configuration all major available features can be enabled, as far as they can be activated simultaneously.

[Shallow Water](ref_buildrun_shallow_water):
: This is a shallow water configuration for ICON's ocean component running on a hex grid. It follows Williamson et al., 1992 (Test 2) and can be regarded as a test of the ocean dynamical core.

## NWP Use Cases

The use cases for NWP (Numerical Weather Prediction) are supported by Deutscher Wetterdienst (DWD).

[NWP Global (R2B06)](ref_buildrun_nwp_global):
: This use case runs a global ICON application with a resolution of about 40 km (R02B06). It includes a nest (R02B07) over Europe.

[NWP Local (R19B07)](ref_buildrun_nwp_local):
: This use case runs a limited-area (local) ICON application over Germany with a resolution of about 2 km (R19B07).
  It is comparable to DWD's operational application ICON-D2.


## ART Use Cases

The use cases for ART (Aerosol and Reactive Trace gases interactions) are available on the [DKRZ Swiftbrowser](https://swiftbrowser.dkrz.de/public/dkrz_4d992e1b-f237-4258-a2bc-138ca6a1cf59/icon-model-use-cases/).
Each use case can be downloaded as a ZIP file containing input data and a README with setup and runtime instructions. Further information is provided in the [ART User guide](https://www.icon-art.kit.edu/userguide/index.php?title=Main_Page).

art\_global\_r02b04\_NWP\_GASPHASE:
: This use case simulates the global atmosphere with a resolution of approximately 160 km (R02B04) including ART. It includes detailed gas-phase chemistry for ozone, particularly the extended Chapman cycle.

art\_global\_r02b05\_NWP\_LIFETIME:
: This use case simulates the global atmosphere with a resolution of approximately 80 km (R02B05) including ART. It employs a simplified lifetime-based approach for most chemical tracers and a linearized ozone scheme. Regional tracers are also activated in this configuration.

art\_global\_r02b05\_NWP\_OH\_CHEMISTRY:
: This use case simulates the global atmosphere with a horizontal resolution of approximately 80 km (R02B05) including ART. It focuses on the representation of hydroxyl radical (OH) chemistry and its interactions with other reactive species.

art\_global\_r02b06\_ALLAERO\_NORAD:
: This use case simulates the global atmosphere with a resolution of approximately 40 km (R02B06) including ART. It incorporates various types of aerosols (dust, sea salt, wildfire) and their emissions, but does not consider their impact on radiation.

art\_global\_r02b06\_DUST\_RAD:
: This use case simulates the global atmosphere with a resolution of approximately 40 km (R02B06), utilizing a nested high-resolution region (20 km, R02B07) for specific areas. ART is activated. It focuses on dust, including its emissions and radiative effects.

art\_global\_r02b06\_VOLAERO\_RAD:
: This use case simulates the global atmosphere with a resolution of approximately 40 km (R02B06) including ART. It incorporates both chemical and aerosol tracers, enabling the simulation of volcanic eruptions and their atmospheric impact.

art\_local\_LAM\_OEM:
: This testcase demonstrates the ICON-ART Online Emission Module (OEM) coupled with the regional atmospheric model. The model calculates biogenic and anthropogenic CO2 fluxes using the Vegetation Photosynthesis and Respiration Model (VPRM), driven by dynamic emission factors, regional vegetation parameters, and meteorological input. The domain covers central Europe (focused on Switzerland) at a horizontal resolution of approximately 2.5 km (R19B09). Emissions are modulated in time using diurnal, weekly, and annual scaling factors defined by OEM input files. Boundary conditions are updated hourly to simulate realistic atmospheric transport of trace gases.

art\_local\_POLLEN\_SPP:
: This use case employs high-resolution local area mode (LAM) over Europe with a resolution of approximately 6.5 km (R03B08) including ART. It focuses on simulating different types of pollen and the radioactive isotope Cesium-137 (Cs-137).


## ICON SCM and LES

[ICON-SCM/LES documentation](ref_buildrun_icon_scm_les):
: ICON can be run in a single column model (SCM) and large-eddy simulation (LES) mode.
