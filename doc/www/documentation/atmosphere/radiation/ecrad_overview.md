```{eval-rst}
:orphan:
```

(ref_atmosphere_ecrad)=
# ecRad Overview

The NWP physics package of ICON uses the ecRad radiation scheme developed by [ECMWF](https://ecmwf.int).
This page aims to described the configuration and input options for ecRad that are available in ICON.
For a more detailed and complete description of the radiation scheme itself is provided in the [ECMWF Confluence Wiki](https://confluence.ecmwf.int/display/ECRAD).
The original source code is publicly available on [GitHub](https://github.com/ecmwf-ifs/ecrad), please note that ICON uses a modified version of this code.

(ref_atmosphere_ecrad_input)=
# Radiation Input Options

(ref_atmosphere_ecrad_gases)=
## Gas Input Options

There are multiple options for the specification of several components of the gaseous composition of the atmosphere available.
The corresponding namelist parameters are [`irad_h2o`](radiation_nml-irad_) for water vapor, [`irad_o3`](radiation_nml-irad_) for ozone, [`irad_co2`](radiation_nml-irad_) for carbon dioxide, [`irad_n2o`](radiation_nml-irad_) for nitrous oxide, [`irad_ch4`](radiation_nml-irad_) for methane, [`irad_o2`](radiation_nml-irad_) for oxygen, [`irad_cfc11`](radiation_nml-irad_) for trichlorofluoromethane and [`irad_cfc12`](radiation_nml-irad_) for dichlorodifluoromethane.

### External specification

For all of the above described gases, the option `-1` (e.g., [`irad_h2o`](radiation_nml-irad_)`=-1`) allows for an external specification of the gaseous concentrations. A variable `<gas>rad_ext` (e.g., `h2orad_ext`) is created which can be filled with mass mixing ratios ({math}`kg\,kg^{-1}`) from an external source, for example via the [Community Interface **ComIn**](ref_tools_comin). There is no cross-check that the arrays contain meaningful values. This is left to the user.

(ref_atmosphere_ecrad_aerosol)=
## Aerosol Input Options

(ref_atmosphere_ecrad_aerosol_tegen)=
### Tegen climatology

Climatological aerosol based on the {term}`Tegen et al. 1997` climatology can be selected by choosing **[`irad_aero`](radiation_nml-irad_aero)`=6`**.
This options has the following characteristics:

- Optical thicknesses at the wavelength 550 nm of the 5 species **Sea Salt**, **Soil Dust**, **Sulfate**, **Organic Carbon** and **Black Carbon** are provided in the [external parameter file](ref_buildrun_external_param).
- The annual cycle is considered by providing monthly data which is linearly interpolated inside ICON to the target date.
- The original data is vertically integrated optical thickness. For the use in ecRad, an exponentially decaying, normalized vertical profile is added by ICON.
- The target variables optical thickness (SW/LW), single scattering albedo (SW) and asymmetry parameter (SW) at the radiation wavelength bands are derived based on lookup tables in the ICON code.

### Simplified Prognostic Aerosol Module _Prog2DAero_

See [here](ref_miscnwp_2daero) for further information.

(ref_atmosphere_ecrad_aerosol_cams)=
### CAMS climatology or CAMS forecast aerosol

ICON currently supports the use of either the CAMS 49R2 aerosol climatology, or CAMS forecast aerosol for direction aerosol-radiation interactions. Aerosol-cloud interactions using CAMS is explained [here](ref_miscnwp_aero_cpl).

#### Using the 49R2 CAMS climatology

ICON supports use of the recent (December 2024) CAMS climatology, version 49R2. This option is activated with the namelist parameter **[`irad_aero`](radiation_nml-irad_aero)`=7`**.

Older versions (43R3) are no longer supported. Aerosol mixing ratios are supplied monthly on 21 pressure surfaces.
Aerosol species affected by human activity have an additional dimension `epoch` covering thirteen 5-year long periods from 1955 to 2015.
The original climatology file can be downloaded from [this ECMWF webpage](https://aux.ecmwf.int/ecpds/home/radiation/aerosol_climatology/aerosol_cams_climatology_49r2_1951-2019_4D.nc).

At this point in time, ICON does **not** support the new `epoch` dimension of the climatology (i.e. anthropogenic change of aerosol over 5+ year periods).
The ICON repository contains the script {{ '[`make_camsclim_onICONgrid.sh`]({}/scripts/preprocessing/make_camsclim_onICONgrid.sh)'.format(base_url) }} to extract the latest (2015) epoch from the original data file and interpolate the climatology onto an ICON grid of the user’s choice.
Installation of CDO, NCO and python3 (numpy, xarray) tools is required to run this script.
The CAMS climatology is based on reanalysis data from 2003-2014 ({term}`Bozzo et al. 2020`), providing 3D monthly aerosol profile climatologies.

#### Using CAMS forecasts

ICON can also use CAMS forecast aerosol fields on 137 model levels. This option is activated with the namelist parameter **[`irad_aero`](radiation_nml-irad_aero)`=8`**.

CAMS forecast aerosol fields can be retrieved from ECMWF via MARS request. The ICON repository contains the script {{ '[`make_camsforc_onICONgrid.sh`]({}/scripts/preprocessing/make_camsforc_onICONgrid.sh)'.format(base_url) }} which then interpolates the CAMS forecast aerosol onto an ICON grid of the user’s choice.
The script header contains more information on how to retrieve CAMS forecast aerosol from MARS.
CDO, python3 and ecmwf-toolbox are required to run this script.

#### Information relevant to both CAMS climatology and forecast


:::{admonition} A note on remapping to the ICON grid
:class: admonition-icontheme
CDO remapping is used for the interpolation. None of the available remapping options are perfect:
- Conservative remapping (**remapcon**) leads to visible 'squares' corresponding to the original, coarser CAMS climatology in the interpolated fields, which translates into visible 'square' shapes in the clear-sky radiation also. This is obviously undesireable.
- Bicubic remapping (**remapbic**) leads to 'overshooting' features around high orography (Himalayas, Andes). Also undesireable.
- The option considered to be best at the moment (and **implemented by default**) is bilinear remapping (**remapbil**), which produces a reasonably smooth field without obvious overshooting.

However, the 'best' option may depend on the application: for regional simulations away from steep orography, **remapbic** may be more advantageous (smoother).
:::

**Vertical interpolation onto ICON model levels**

ICON interpolates from the original CAMS vertical levels (21 pressure surfaces for climatology, 137 IFS model levels for forecast) to the ICON vertical model levels. Because the horizontal resolution of the original CAMS files can be much lower than the ICON resolution, the difference in surface pressure between ICON and CAMS can be in excess of 200hPa around steep orography. A straight interpolation between pressure levels would therefore neglect a significant amount of aerosol mass if the lowest 200hPa of the CAMS profile were ignored.
To avoid this, the CAMS profile is re-distributed between the ICON surface pressure and top of the atmosphere before interpolation, conserving total aerosol mass in the column.

**Differences to previous CAMS versions**

The previous CAMS climatology 43R3 provided aerosol as **layer mass**. This climatology now also exists (v2) in a format providing **mixing ratios** instead. Since layer mass will no longer be used from the CAMS side, reading in layer mass is no longer supported by ICON.

The 43R3 CAMS climatology also had some unrealistic features, such as high dust accumulations in the stratosphere. These artefacts have now been removed in the 49R2 climatology. Since it is unlikely that the older, flawed climatology will continue to be used, the option to read in the 43R3 climatology is also no longer supported. However, it is easy to make ICON compatible with this older version again, if so desired, by changing the number of CAMS levels (`nlev_cams`) in `mo_reader_cams.f90` from 21 back to 60, and using the more recent mixing ratio version of the [43R4 CAMS climatology which can be obtained from ECMWF](https://aux.ecmwf.int/ecpds/home/radiation/aerosol_climatology/aerosol_cams_climatology_43r3_v2_3D.nc).

**Appropriate aerosol optical properties**

The new 49R2 CAMS climatology was created using an updated version (CY48R1) of the CAMS aerosol model, which contains some significant changes relative to the older version used to create the 43R3 climatology. Therefore, the **new climatology should be used with an appropriate set of optical properties!** A tabulated list of the appropriate optical properties for each version of the climatology/forecasts can be found **[here](https://confluence.ecmwf.int/display/ECRAD/Aerosol-radiation+interactions+in+the+IFS)**.

For this implementation, the appropriate aerosol optical properties have been pre-selected for use with the CAMS 49R2 climatology. The pre-selected default for the CAMS forecasts is set to the most recent CAMS IFS cycle 49R1. If working with **older CAMS forecasts**, the **user has to adapt the selected properties** in `mo_nwp_ecrad_init.f90` (line 342 and following) according to the table linked above.

:::{admonition} Known limitations
:class: admonition-icontheme
Known limitations of the new 49R2 climatology include that the "far field" aerosol such as in the Arctic is too low.
The IFS is run with an additional artificial small constant background term to get the best results.
:::

### Aerosol options for simulations with transient climatologies

For climate simulations, the aerosol interface can also handle optical properties prescribed by transient aerosol climatologies, which are used for CMIP-type simulations. Input data are available on `levante` (e.g., `/pool/data/icon-xpp/common` or `/pool/data/CLMcom/ICON-CLM/data/rcm_new`) and comprise
- **tropospheric aerosol**: the full Kinne climatology ([MACv2](https://doi.org/10.5194/acp-19-10919-2019)) providing monthly fields for 1845--2100 or the simple plume version, which is modelling the full fields on the basis of nine spatial plumes associated with different major anthropogenic source regions  ([MACv2-SP](https://doi.org/10.5194/gmd-10-433-2017)). An update of the simple plume scheme is available for CMIP7 ([SPv2](https://zenodo.org/records/15283189); [Fiedler et al., 2025](https://doi.org/10.1029/2025MS005067)).
- **stratospheric aerosol**: the version used for CMIP6, based on an extension of the Pinatubo aerosol data set ([Stenchikov et al., 1998](https://doi.org/10.1029/98JD00693)) to cover the entire period between 1850 and 2014 ([Stevens et al., 2019](https://doi.org/10.1002/jame.20015)). For CMIP7, input data was retrieved from the [input4MIPs project](https://input4mips-cvs.readthedocs.io/en/latest/dataset-overviews/stratospheric-volcanic-so2-emissions-aod/) and converted into a similar format. For the historical period, the years 1850--2023 are available, as well as a "background" climatology for 1850.

The Kinne input data have to be interpolated offline onto the ICON grid used for the current simulation (an example script is provided at {{ '[`interpolate_KinneAero.sh`]({}/scripts/preprocessing/interpolate_KinneAero.sh)'.format(base_url) }}). An option for online interpolation using YAC is implemented as well, a corresponding description will be provided in an ICON report. The simple plume data are independent of the ICON grid so that while using these, only the Kinne background fields for the year 1850 have to be interpolated.\
The stratospheric aerosols are provided as zonal averages and do not have to be regridded.

The namelist switch **[`irad_aero`](radiation_nml-irad_aero)`=12,13,14,15,18,19`** provides different combinations of the transient tropospheric and stratospheric climatologies. Most commonly used options are [`irad_aero`](radiation_nml-irad_aero)` = 18` for historical simulations and [`irad_aero`](radiation_nml-irad_aero)` = 19` for climate projection scenarios.

(ref_atmosphere_ecrad_cdnc)=
## Cloud droplet number concentration (cdnc)

MODIS retrievals of cloud droplet number concentration can be used in ICON when the namelist parameter [`icpl_aero_gscp`](nwp_phy_nml-icpl_aero_gscp) is set to 3. This option was implemented and tested for the Kinne aerosols (i.e., the full climatology or the simple-plume version, see **[`irad_aero`](radiation_nml-irad_aero)**), as an alternative to the coupling between aerosols and clouds when using the Tegen aerosol climatology.
The external parameter file must then contain the field **cdnc** and can be generated with **[Extpar](https://docs.icon-model.org/tools/tools.html#ref-tools-gridextpargui)** (option `enable_cdnc` in the `Expert` mode of Zonda).
When using the external cdnc, it is advisable to set:
- [`icpl_aero_conv`](nwp_phy_nml-icpl_aero_conv)` = 1`: simple coupling between auto-conversion (in convection scheme) and aerosol

The external climatological cdnc can be scaled when setting namelist parameter [`scale_cdnc_mode`](nwp_phy_nml-scale_cdnc_mode) to 1 or 2, and providing the necessary input file for the Simple Plume scheme in the run directory:

- [`scale_cdnc_mode`](nwp_phy_nml-scale_cdnc_mode)` = 0` (default): no scaling; the cdnc used by ICON will be the same regardless of the simulation year.
- [`scale_cdnc_mode`](nwp_phy_nml-scale_cdnc_mode)` = 1`: apply a yearly varying scaling of cdnc, which is derived from scaling factor of the Simple Plume scheme (e.g. for climate simulations for the historical period after 1850 or for climate projections).
- [`scale_cdnc_mode`](nwp_phy_nml-scale_cdnc_mode)` = 2`: apply a constant scaling of cdnc representative for the year 1850 (e.g. for pre-industrial experiment)

Note that [`scale_cdnc_mode`](nwp_phy_nml-scale_cdnc_mode)` = 1,2` is only possible in combination with [`icpl_aero_gscp`](nwp_phy_nml-icpl_aero_gscp)` = 3`.

(ref_atmosphere_ecrad_fsd)=
## Condensate heterogeneity - the FSD parameter

ICON predicts one value for the condensate mixing ratios of liquid and ice for each grid box. This is supplemented by a cloud fraction from the diagnostic cloud scheme. Because the amount of radiation reflected or absorbed depends non-linearly on condensate amount, it matters how this condensate is distributed within each grid box. The default assumption is that the condensate amount is distributed within the cloudy part of the grid box with the functional shape of a Gamma distribution. The distribution average corresponds to the predicted grid box condensate amount, while the width of the distribution is given by the "fractional standard deviation" (FSD) parameter, which is defined as


```{math}
FSD=\frac{standard\, deviation}{ mean}
```

By default, ICON assumes that `FSD=1` everywhere, i.e. the normalised width of the condensate distribution is the same everywhere. However, observations show that this is not the case. Condensate is distributed more homogeneously in stratiform clouds compared to cumuliform clouds. Also, grid boxes containing cloud edges (i.e. not overcast, with a `cloud fraction < 1`) have wider condensate distributions than overcast grid boxes, because cloud edges naturally contain less condensate than the cloud interior.

To account for this effect, a regime-dependent parameterization for the FSD parameter can be used by setting the namelist parameter in the `radiation_nml`:

[`lcalculate_fsd`](radiation_nml-lcalculate_fsd)` = .true.`

This parameterization is based on the publications {term}`Ahlgrimm et al. 2016` and {term}`Ahlgrimm et al. 2017`, with some minor modifications documented in the ICON code. Broadly, the effect of using this parameterization is to make clouds appear more reflective in areas dominated by stratiform clouds (e.g. stratocumulus regions, extratropics) and less reflective in regions dominated by more convective cloud (e.g. trade cumulus regions, tropics).

Two additional parameters may be set in the `radiation_nml`:

[`fsd_background`](radiation_nml-fsd_background)` = 1` is the default value used when the parameterization is switched off entirely, or in cloud-free regions of the model atmosphere when the parameterization is active. When radiation is calculated on the reduced grid, the interpolation from the full grid to the reduced grid may interpolate between cloudy and cloud-free grid points, which therefore must be assigned a valid FSD value.

[`fsd_gridlen`](radiation_nml-fsd_gridlen(max_dom))` = 80` is the assumed horizontal grid spacing of the ICON grid (by default set to 80km). Observations show that the unresolved condensate heterogeneity that must be parameterized should reduce as the resolution of the model increases. This means the parameterized FSD parameter becomes smaller (clouds become more homogeneous) for a finer grid resolution. However, ICON has been operating with a fixed FSD value of 1 at all resolutions for years, and produces a resolution-independent top-of-the-atmosphere radiation balance with this fixed value. Replacing this fixed value with a resolution-dependent FSD value (in the absence of compensating changes elsewhere) would produce a resolution-dependent TOA radiation balance, which is not desirable. Therefore, for the time being, it is recommended to use the fixed gridlength value of 80km at all resolutions, as this produces FSD values that average out to approximately 1 globally, maintaining the usual TOA radiation balance.

Lastly, it should be mentioned that the FSD calculation for liquid clouds depends on the cloud fraction: High cloud fraction is a proxy for stratiform clouds, which are assigned lower FSD values. In cases where the model predicts an incorrect cloud fraction (e.g. prediction cloud fractions <50% in stratocumulus regions), the error in the cloud radiative effect may be enhanced when using the FSD parameterization. In this example, the parameterization would make clouds with fraction <50% less reflective in an area where cloud cover (and therefore cloud radiative effect) is already too low.


(ref_atmosphere_ecrad_outputVars)=
# Radiation output variables

  These lists should provide a quick overview on output variables for radiation.
  The variables are defined in {{ '[`mo_nwp_phy_state`]({}/src/atm_phy_nwp/mo_nwp_phy_state.f90)'.format(base_url) }}.
  Diagnostic ouput variables are often computed in {{ '[`mo_nwp_diagnosis`]({}/src/atm_phy_nwp/mo_nwp_diagnosis.f90)'.format(base_url) }}.
  The diagnostic derived type is defined in {{ '[`mo_nwp_phy_types`]({}/src/atm_phy_nwp/mo_nwp_phy_types.f90)'.format(base_url) }}.

  Net fluxes are defined as downward positive.

## Shortwave

   Some shortwave output variables exist without and with orographic shading (1) and with slope-dependent and orographic shading (2).
   The variables ending on either `_os` or `_tan_os` are defined depending on the namelist parameter `islope_rad(dom)` in the `radiation_nml`.

<a name="fradswback">

   | icon variable name | grib2 name<sup><a href="#frad1">1</a></sup>  | direction | specifics      | acc.<sup><a href="#frad2">2</a></sup>        | location      | Description |
   |--------------------|-------------|-----------|----------------|-------------|---------------|-------------------|
   | asob_s             | asob_s      | net       |                | x           | surface       | Surface net solar radiation since model start |
   | asob_t             | asob_t      | net       |                | x           | TOA           | TOA net solar radiation since model start |
   | asobclr_s          | asob_s_cs   | net       | clear-sky      | x           | surface       | Clear-sky surface net solar radiation since model start |
   | asod_s             | asod_s      | down      |                | x           | surface       | Surface down solar rad. since model start |
   | asod_s_os          |             | down      |                | x           | surface (1)   | Surface down solar rad. incl. orographic shading since model start |
   | asod_s_tan_os      |             | down      |                | x           | surface (2)   | Surface down solar rad. incl. slope-dependent and orographic shading since model start |
   | asod_t             | asod_t      | down      |                | x           | TOA           | Top down solar radiation |
   | asodifd_s          | aswdif_s    | down      | diffuse        | x           | surface       | Surface down solar diff. rad. since model start |
   | asodifu_s          | asdifu_s    | up        | diffuse        | x           | surface       | Surface up solar diff. rad.  since model start|
   | asodifu_s_os       |             | up        | diffuse        | x           | surface (1)   | Surface up solar diff. incl. orographic shading since model start |
   | asodifu_s_tan_os   |             | up        | diffuse        | x           | surface (2)   | Surface up solar diff. incl. slope-dependent and orographic shading since model start |
   | asodird_s          | aswdir_s    | down      | direct         | x           | surface       | Surface down solar direct rad. since model start |
   | asodird_s_os       |             | down      | direct         | x           | surface (1)   | Surface down solar direct rad. incl. orographic shading since model start |
   | asodird_s_tan_os   |             | down      | direct         | x           | surface (2)   | Surface down solar direct rad. incl. slope-dependent and orographic shading since model start |
   | asou_t             | uswrf       | up        |                | x           | TOA           | Top up solar radiation since model start |
   | sob_s_t_*<sup><a href="#frad3">3</a></sup>| sobs_rad  | net | |             | surface, tile | tile-based shortwave net flux at surface |
   | sob_s              | sobs_rad    | net       |                |             | surface       | shortwave net flux at surface |
   | sob_s_os           |             | net       |                |             | surface (1)   | shortwave net flux at surface incl. orographic shading |
   | sob_s_tan_os       |             | net       |                |             | surface (2)   | shortwave net flux at surface incl. slope-dependent and orgraphic shading |
   | sob_t              | sobt_rad    | net       |                |             | TOA           | shortwave net flux at TOA |
   | sobclr_s           | sobs_rad    | net       | clear-sky      |             | surface       | net shortwave clear-sky flux at surface |
   | sod_t              | sodt_rad    | down      |                |             | surface       | downward shortwave flux at TOA |
   | sodifd_s           | swdifds_rad | down      | diffuse        |             | surface       | shortwave diffuse downward flux at surface |
   | sou_s              | swdifus_rad | up        |                |             | surface       | shortwave upward flux at surface |
   | sou_s_os           |             | up        |                |             | surface (1)   | shortwave upward flux at surface incl. orographic shading |
   | sou_s_tan_os       |             | up        |                |             | surface (2)   | shortwave upward flux at surface incl. slope-dependent and orographic shading |
   | sou_t              | uswrf       | up        |                |             | TOA           | shortwave upward flux at TOA |
   | swflx_dn_clr       | -           | down      | clear-sky      |             | 3d            | shortwave downward clear-sky flux |
   | swflx_dn           | -           | down      |                |             | 3d            | shortwave downward flux |
   | swflx_up_clr       | -           | up        | clear_sky      |             | 3d            | shortwave upward clear-sky flux |
   | swflx_up           | -           | up        |                |             | 3d            | shortwave upward flux |
   | trsolall           | -           | net       | transmissivity |             | 3d            | shortwave net transmissivity |


## Thermal

   | icon variable name | grib2 name<sup><a href="#frad1">1</a></sup>  | direction | specifics      | acc.<sup><a href="#frad2">2</a></sup>        | location | Description |
   |--------------------|-------------|-----------|----------------|-------------|----------|-------------------|
   | athb_s             | athb_s      | net       |                | x           | surface  | surface net thermal radiation since model start |
   | athb_t             | athb_t      | net       |                | x           | TOA      | TOA net thermal radiation since model start |
   | athbclr_s          | athb_s_cs   | net       | clear-sky      | x           | surface  | clear-sky surface net thermal radiation since model start |
   | athd_s             | athd_s      | down      |                | x           | surface  | Surface down thermal radiation since model start |
   | athu_s             | athu_s      | up        |                | x           | surface  | Surface up thermal radiation since model start |
   | lwflx_dn_clr       | -           | down      | clear-sky      |             | 3d       | longwave downward clear-sky flux |
   | lwflx_dn           | -           | down      |                |             | 3d       | longwave downward flux |
   | lwflx_up_clr       | -           | up        | clear-sky      |             | 3d       | longwave upward clear-sky flux |
   | lwflx_up           | -           | up        |                |             | 3d       | longwave upward flux |
   | lwflxall           | nlwrf       | net       |                |             | 3d       | longwave net flux |
   | thb_s_t_*          | thbs_rad    | net       |                |             | surface  | tile-based longwave net flux at surface |
   | thb_s              | thbs_rad    | net       |                |             | surface  | longwave net flux at surface |
   | thb_t              | thbt_rad    | net       |                |             | TOA      | thermal net flux at TOA |
   | thbclr_s           | thbt_rad_cs | net       | clear-sky      |             | surface  | net longwave clear-sky flux at surface |
   | thu_s              | thus_rad    | up        |                |             | surface  | longwave upward flux at surface |

## Diagnostic for bands PAR, VIS, NIR
   The bands for diagnostic output are:
   - PAR : photosynthetically active flux ({math}`400 - 700 nm`)
   - VIS : visible ({math}`300 - 700 nm`)
   - NIR : near-infrared ({math}`0.7 - 5 \mu m`)

   | icon variable name    | grib2 name<sup><a href="#frad1">1</a></sup> | direction | specifics   | acc.<sup><a href="#frad2">2</a></sup> | location    | Description |
   |-----------------------|-----------|-----------|------------------|-------------|-------------|-------------------|
   | aswflx_par_sfc        | apab_s    | down      |                  | x           | surface     | Downward PAR flux |
   | aswflx_par_sfc_tan_os |           | down      |                  | x           | surface (2) | Downward PAR flux incl. slope-dependent and orographic shading |
   | fr_nir_sfc_diff       | -         | down      | diffuse fraction |             | surface     | diffuse fraction of downward near-infrared flux at surface |
   | fr_par_sfc_diff       | -         | down      | diffuse fraction |             | surface     | diffuse fraction of downward photosynthetically active flux at surface |
   | fr_vis_sfc_diff       | -         | down      | diffuse fraction |             | surface     | diffuse fraction of downward visible flux at surface |
   | swflx_nir_sfc         |           | down      |                  |             | surface     | downward near-infrared flux at surface |
   | swflx_par_sfc         |           | down      |                  |             | surface     | downward photosynthetically active flux at surface |
   | swflx_par_sfc_tan_os  |           | down      |                  |             | surface (2) | downward photosynthetically active flux at surface incl. slope-dependent and orographic shading |
   | swflx_vis_sfc         |           | down      |                  |             | surface     | downward visible flux at surface |

1. <a name="frad1"/>The grib2 names in the table refer to the short names resulting from the DWD grib definition files.<a href="#fradswback">{octicon}`undo;1em;pst-color-secondary`</a>
2. <a name="frad2"/>Output variables starting with the letter `a` are likely "accumulated" since model start. Depending on the switch [`lflux_avg`](io_nml-lflux_avg), they contain either averages since model start ([`lflux_avg`](io_nml-lflux_avg)`=.true.`) or accumulated values.<a href="#fradswback">{octicon}`undo;1em;pst-color-secondary`</a>
3. <a name="frad3"/>The asterisk stands for the number of the surface tile.<a href="#fradswback">{octicon}`undo;1em;pst-color-secondary`</a>

(ref_ecrad_implementation)=
# Implementation Details

(ref_atmosphere_ecrad_redgrid)=
## Reduced Radiation Grid

Radiation is one of the computationally expensive physical parameterizations.
There are several possibilities to decrease the computational cost of radiation, e.g. by reducing the temporal, spatial or spectral resolution.
By activating **[`lredgrid_phys`](grid_nml-lredgrid_phys)`=.true.`** and specifying the corresponding grid file with **[`radiation_grid_filename`](grid_nml-radiation_grid_filename)**, the radiation is calculated on a one grid level coarser domain to reduce the computational cost of the radiation by about a factor of 4.

**We highly recommend activating the reduced radiation grid option for the following reasons:**

- It is computationally cheaper by a factor of 4, usually without a degradation of the results. Since the radiation is treated as a slow physics process, it is not called every time step anyways. A coarser horizontal grid thus fits better to the advective time scale.
- Only for the reduced radiation grid, there is an additional option named **[`latm_above_top`](nwp_phy_nml-latm_above_top)** (see [here](ref_ecrad_atmabovetop)). This option adds an extra layer at the top to account for the incoming long-wave radiation. This reduces the biases at the model top significantly.
- For global domains, there is a load balancing for sunlit and shadowed parts of the earth for the reduced grid.

For a more detailed description of the reduced radiation grid implementation, see {term}`ICON Tutorial`.

(ref_ecrad_atmabovetop)=
## Considering Atmosphere above the Model Top

The namelist switch [`latm_above_top`](nwp_phy_nml-latm_above_top) is the only way in the ecRad implementation of ICON to take into account the layers above the top of the model ([`top_height`](sleve_nml-top_height)). In this case, an additional layer is added for the radiation calculations above the top of the model, which has 1.5 times the layer thickness of the layer below it (i.e., the top layer of the model).

This additional layer primarily influences long-wave back radiation, which can lead to differences up to 20K in quasi-equilibrium. Not using [`latm_above_top`](nwp_phy_nml-latm_above_top) can lead to large biases. In addition, the values for ozone, aerosols, and clouds from the top model layer are copied to the additional layer (i.e., zero gradient condition). In practice, however, these are usually zero for aerosols and clouds.

:::{admonition} Note on limited-area simulations
:class: admonition-icontheme
Obviously, there is a significant amount of ozone missing from the stratosphere, even when [`top_height`](sleve_nml-top_height) is taken into account. This manifests itself in increased heating rates in the top layer, which are partially compensated for by a corresponding lack of contributions from higher layers in the longwave back radiation.
:::

Vertically nested domains  can have significantly more additional layers via namelist setting. For ICON-EU, for example, [`nexlevs_rrg_vnest`](grid_nml-nexlevs_rrg_vnest)`=14`. In this case, atmospheric variables including trace gases are copied from the parent domain to the (equal-resolution) radiation grid of the nested domain. Note that in case of multiple vertical nesting, this mechanism does not work recursively. That is, the number of layers that can be added above a vertical nest interface is limited by the additional number of model levels the direct parent domain has.

{material-regular}`warning;2em;pst-color-secondary` _If not using a reduced radiation grid, i.e. calculating radiation on the dynamics grid, neither longwave back radiation nor ozone absorption are taken into account above the model top. Especially for limited-area configurations with a model top in the lower stratosphere, this can lead to large biases._

(ref_ecrad_debugging)=
## Debugging & known issues

The namelist option [`ecrad_check_input`](radiation_nml-ecrad_check_input) can be used to perform a series of checks and adjust some settings:

- ecRad verbosity is increased significantly
- Sanity checks on pressure input fields
- Abort if the combination Intel+OpenMP+ecRad is used (see below)
- Make use of `out_of_physical_bounds` checks provided by ecRad

{material-regular}`warning;2em;pst-color-secondary` _Two fields (`zrg_clc` and `opt_ptrs_sw(jw)%ptr_ssa`) contained values slightly larger than 1 after upscaling (rounding error). These fields are clipped to 1 if [`ecrad_check_input`](radiation_nml-ecrad_check_input)` = .TRUE.`._

:::{admonition} Combining Intel + OpenMP + ecRad
:class: admonition-icontheme
There is a known issue for the combination Intel + OpenMP + ecRad that leads to non-reproducible results.
Since the cause of this error could not be identified yet, we highly recommend to avoid this combination.
Note that ecRad with OpenMP in ICON works for other compilers.
:::

(ref_ecrad_single_precision)=
## Single precision

The ecRad radiation scheme supports single-precision computation, which improves
performance and reduces memory usage.
This approach is accurate enough for most applications.
It is especially beneficial for high-resolution simulations or large ensemble
runs, where computational efficiency is crucial.
You can enable single-precision computation via the appropriate configure
option.
Refer to the output of `./configure --help` for details.

Please note that the quality of this feature has not been evaluated for global
applications (i.e., altitudes higher than 25 km).
Users are therefore advised to run a benchmark simulation before using this
combination.
