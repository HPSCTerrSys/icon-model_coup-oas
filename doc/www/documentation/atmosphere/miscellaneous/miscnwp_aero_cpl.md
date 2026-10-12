```{eval-rst}
:orphan:
```

(ref_miscnwp_aero_cpl)=
# Aerosol-Cloud Coupling in ICON

## Overview

The interaction between aerosols and clouds influences both radiation (via indirect radiative effects) and microphysical processes. The ICON 1-moment microphysics scheme requires Cloud Droplet Number Concentration (CDNC) for calculating both effective radius and auto-conversion rates. The coupling method is controlled by the namelist parameter **[`icpl_aero_gscp`](nwp_phy_nml-icpl_aero_gscp)**.
The implementation is based on the {term}`Segal and Khain 2006` (SK2006) droplets activation parametrization, which in its complete form uses a 4D lookup table with four parameters: CCN concentration (range: 50-6400 cm⁻³), CCN mode radius (range: 0.02-0.04 µm), geometric standard deviation (range: 0.1-0.5), and updraft velocity at cloud base (range: 0.5-5.0 m/s).
Different configuration options provide varying levels of complexity, from constant CDNC to full 4D parametrization. When using [CAMS aerosols](ref_atmosphere_ecrad_aerosol_cams), the system processes 3D aerosol fields containing 11 tracers: three dust bins (DU), three sea salt bins (SS), sulphate (SU), hydrophobic and hydrophilic black carbon (BC), and hydrophobic and hydrophilic organic matter (OM).

## Configuration Options

**[`icpl_aero_gscp`](nwp_phy_nml-icpl_aero_gscp)`=0`**: Constant CDNC - no aerosol-cloud interactions (Default)
:   Method: Initializes CDNC with a fixed value of 200 cm⁻³ that remains constant throughout the model run

**[`icpl_aero_gscp`](nwp_phy_nml-icpl_aero_gscp)`=1`**: Simplified Tegen-SK2006 Coupling
:   Requirements: **[`irad_aero`](radiation_nml-irad_aero)`=6, 7, or 9`**
:   Method: Implements the SK2006 parametrization with simplifications:
    -	Derives CCN concentration from [Tegen's](ref_atmosphere_ecrad_aerosol_tegen) 2D AOD using prescribed vertical profiles
    -	Combines 5 aerosol species into a single concentration using solubility-based weighting (0.1 for dust, 0.9 for organics)
    -	Uses fixed parameters: mode radius = 0.03 µm, geometric standard deviation = 0.3, vertical wind speed = 0.25 m/s
    -	Reduces the full 4D lookup table to 1D, making CDNC dependent solely on CCN concentration
    -	Calculates CDNC at cloud base and extends to 3D based on pressure levels

**[`icpl_aero_gscp`](nwp_phy_nml-icpl_aero_gscp)`=2`**: Full CAMS-SK2006 Coupling
:   Requirements: **[`irad_aero`](radiation_nml-irad_aero)`=7 or 8`**
:   Method: Implements full 4D-LUT SK2006 parametrization
    -	Species-specific number concentrations, mode radii, and geometric standard deviations
    -	Hygroscopicity factors control each aerosol type's contribution
    -	Currently uses grid-scale updraft velocity
    -	Direct 3D CDNC calculation without vertical profile assumptions

**[`icpl_aero_gscp`](nwp_phy_nml-icpl_aero_gscp)`=3`**: External CDNC Climatology
:   Data source: External parameter file containing cloud-droplet number climatology (can be generated using [Zonda](ref_tools_gridextpargui))
:   Method: Reads 2D climatology and extends to 3D using the same approach as Option 1. For further information see [here.](ref_atmosphere_ecrad_cdnc)
