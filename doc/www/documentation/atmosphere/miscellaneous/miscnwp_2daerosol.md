```{eval-rst}
:orphan:
```

(ref_miscnwp_2daero)=
# Simplified Prognostic Aerosol Module _Prog2DAero_
{{author}}D. Rieger, DWD{{endauthor}}

ICON contains a simplified, computationally cheap aerosol module _Prog2DAero_. The additional cost in a global simulation is typically below 1%. The idea is to account for the main causes of variability in atmospheric aerosol concentrations not reflected in climatological values. It is based on the [Tegen climatological aerosol](ref_atmosphere_ecrad_aerosol_tegen) which is used for operational numerical weather prediction with ICON. Severeal design choices for _Prog2DAero_ were made based on characteristics of this climatology ([see here](ref_atmosphere_ecrad_aerosol_tegen)). _Prog2DAero_ makes use of the vertically integrated optical depth of the five species Sea Salt, Soil Dust, Sulfate, Organic Carbon and Black Carbon as prognostic variable and adds the following process descriptions:

- Horizontal transport using vertically averaged wind fields
- Artificial diffusion to account for vertical wind shear
- Emissions based on physical parameterizations for natural aerosol species (sea salt, soil dust)
- Emissions based on datasets for anthropogenic and wildfire emissions
- Source term for Sulfate to account for nucleation from the gas phase
- Wet deposition due to grid-scale and convective processes
- Relaxation approach to account for other sink terms (sedimentation, dry deposition)

The prognostic equation for aerosol species `j` can be written as

```{math}
\frac{\partial \psi_j}{\partial t} =
\overline{v_{H,j}}\nabla\psi_j
+ c_{diff,j} \nabla^2 \psi_j
+ S_{e,j}
+ S_{w,j}
+\frac{\psi_{clim,j}-\psi_j}{\tau_{clim,j}}
```

where the right hand side terms are advection, artificial diffusion, source terms, sink terms and relaxation.
{math}`\psi_j` denotes the vertically integrated optical depth, {math}`\overline{v_{H,j}}` the vertically averaged horizontal wind, {math}`c_{diff,j}` a diffusion coefficient, {math}`S_{e,j}` source terms (emission), {math}`S_{w,j}` sink terms (washout), {math}`\psi_{clim,j}` the climatological value of vertically integrated optical depth and {math}`\tau_{clim,j}` a relaxation time scale.

Parts of _Prog2DAero_ can be activated individually. By chosing the namelist parameter {term}`i2daero_dust` `=1`, the simplified prognostic description for soil dust can be activated. Likewise, {term}`i2daero_seas` `=1` activates sea salt. The three species Sulfate, Organic Carbon and Black Carbon can be activated by {term}`i2daero_anthro` `=1`. In addition, emissions from wild fires can be added to these species using the namelist option {term}`i2daero_fire`.

## Combination with ICON-ART

Since the more sophisticated aerosol from [ICON-ART](ref_atmosphere_art) uses [Tegen climatological aerosol](ref_atmosphere_ecrad_aerosol_tegen) complementary for species that are not considered by the [ICON-ART](ref_atmosphere_art) setup, _Prog2DAero_ can be used complementary to [ICON-ART](ref_atmosphere_art) in a similar way. For example by using natural aerosol (Mineral Dust, Sea Salt) from [ICON-ART](ref_atmosphere_art) and anthropogenic and wildfire species (Sulfate, Organic Carbon, Black Carbon) from _Prog2DAero_.

# Glossary of Namelist Parameters

_Operational NWP setting marked by {material-regular}`settings;1em;pst-color-secondary`_

:::{glossary}
i2daero_dust
  (`nwp_phy_nml`) Activate soil dust from _Prog2DAero_ 0:{material-regular}`settings;1em;pst-color-secondary` deactivated, 1: activated

i2daero_seas
  (`nwp_phy_nml`) Activate sea salt from _Prog2DAero_ 0:{material-regular}`settings;1em;pst-color-secondary` deactivated, 1: activated

i2daero_anthro
  (`nwp_phy_nml`) Activate anthropogenic aerosol (Sulfate, Organic Carbon and Black Carbon) from _Prog2DAero_ 0:{material-regular}`settings;1em;pst-color-secondary` deactivated, 1: activated _requires `enable_edgar` during [ExtPar dataset generation with Zonda](ref_tools_gridextpargui)_

i2daero_fire
  (`nwp_phy_nml`) Activate wild fire aerosol (Sulfate, Organic Carbon and Black Carbon) from _Prog2DAero_ 0:{material-regular}`settings;1em;pst-color-secondary` deactivated, 1: Read from file 2: Seasonal climatology
:::
