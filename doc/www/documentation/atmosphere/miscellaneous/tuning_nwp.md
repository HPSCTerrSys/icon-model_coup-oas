```{eval-rst}
:orphan:
```

(ref_atmosphere_tuning)=
# Tuning parameters

_Please note that this description applies to the [](ref_atmosphere_nwp_physics)_.

To get the best performance out of a model for a specific application, or for a particular area on the globe, it can be useful to change certain parameters of the model specification within reasonable bounds.

**ICON model parameters suitable for model tuning**

The tables summarize the most important tuning variables for the ICON model, and are largely based on {term}`Reinert et al. 2025`, chapter 12.2 and {term}`Avgoustoglou et al. 2020`. Yet, the document and the list of variables should be handled with care. Purely varying some of the listed parameters blindly will most likely not give satisfactory results. A physical understanding of the identified model shortcomings/biases should be built up first, followed by a choice of the associated model parameters and a systematic variation and evaluation of simulations. The parameters of interest may strongly vary for the region of interest, the model resolution and the specific purpose.
Please also keep in mind that the list is neither exhaustive, nor complete. There may well be further model parameters that are more suitable for individual applications.

## SSO tuning
{material-regular}`warning;2em;pst-color-secondary` Tuning of the SSO and GWD parameters is dependent on the employed external parameters.

|Parameter     |  Description            |  Meaningful Range                                    |    Comment|
|-----------|--------------|--------|-------|
|  gkwake        |  low level wake drag constant {math}`C_d` for blocking |    {math}`1.5 \pm 0.5`   |Very strong dependency on raw data resolution: for ICON-D2 with ASTER data, we use 0.25|
|  gkdrag        |  gravity wave drag constant {math}`G`, a function of mountain sharpness  |       {math}`0.075 \pm 0.04`     | Should be zero (turned off) at convection-permitting resolutions|
|  gfrcrit       |  critical Froude number determining depth of blocked layer {math}`H_{n_{crit}}`| {math}`0.4 \pm 0.1`||
|  grcrit        |  critical Richardson number|     0.25||
|  tune_minsso   |  minimal value of SSO-STDH (m) where  SSO-effects are being considered | default 10 | must also be adapted in extpar!|
|  tune_blockred |  multiples of the  SSO-STDH, above which the SSO-blocking tendency is being reduced proportionally to STDH/z_AGL| 15 | default 100 = deactivated |


## GWD tuning

|Parameter     |  Description            |  Meaningful Range                                    |    Comment|
|-----------|--------------|--------|-------|
|     gfluxlaun  | variability range for non-orographic gravity wave launch momentum flux |  {math}`2.50 \cdot 10^{-3}` {math}`\pm  0.75\cdot 10^{-3}` [Pa] | relevant for global applications only |

## grid scale  microphysics

|  Parameter  |      Description |             Meaningful  Range |  Comment |
|  --------|--------|--------------|-----|
| zvz0i    |    terminal fall velocity  of ice |  {math}`0.85 \pm 0.25` [m/s] |   allows temperature bias tuning in the upper tropical troposphere as well as TOA long-wave fluxes |
| zceff_min|       minimum value for sticking efficiency |       0.01 - 0.075 | tropics|
| zcsg     | efficiency for cloud-graupel riming | 0.5  | |
|  v0snow  |       factor in the terminal velocity for snow |  10.0 - 30.0  | depending on microphysics scheme, see gscp_data.f90 |
|  icesedi_exp   |   exponent for density correction of cloud ice sedimentation | 0.3 - 0.33  |  no perturbation recommended |
|  rain_n0fac   |    multiplicative change of intercept parameter of raindrop size distribution | 0.25 - 4. | multiplicative perturbation |

## cloud cover

|  Parameter  |      Description |             Meaningful  Range |  Comment |
|  --------|--------|--------------|-----|
| box_liq  |    Box width for liquid clouds assumed in the cloud cover scheme | {math}`0.05 \pm 0.02` ||
| box_liq_asy  |    Asymmetry factor for liquid cloud cover diagnostic | 2.0 - 4.0 (def. 3.0) |  sensitive to TOA solar fluxes and to a lesser degree long-wave fluxes |
| box_liq_sfc_fac | Tuning factor for box_liq reduction near the surface | 1.0 ||
| box_ice  |    Box width for ice clouds assumed in the cloud cover scheme | 0.05 ||
| thicklayfac   |   factor for increasing the box width for layer thicknesses exceeding 150 m |   {math}`0.005 \pm 0.005` [1/m] |  accounting for vertical sub-grid overlap |
| sgsclifac  |     Scaling factor for turbulence-induced subgrid-scale contribution to diagnosed cloud ice | 0.0 - 1.0 |  0.0 turns this contribution off |
| supsat_limfac | Limiting factor for allowed supersaturation in satad | 0. ||
| allow_overcast |  Tuning factor for steeper dependence CLC (RH) |      {math}`\leq 1.0` |   setting allow_overcast {math}`\leq 1` together with reduction of tune_box_liq_asy causes steeper CLC(RH) dependence. {material-regular}`warning;2em;pst-color-secondary` recommendation: allow overcast<1 should not be used in combination with lsgs_cond=.TRUE. |

## turbulence

|  Parameter  |      Description |             Meaningful  Range |  Comment |
|  --------|--------|--------------|-----|
| q_crit |      critical value for normalised super-saturation |   {math}`1.6 \pm 1.0` ||
| rlam_heat | scaling factor of the laminar boundary layer for heat (scalars), the change in rlam_heat is accompanied by an inverse change of rat_sea in order to keep the evaporation over water (controlled by rlam_heat{math}`\cdot`rat_sea) the same. {material-regular}`warning;2em;pst-color-secondary` recommendation: the product of rlam heat and rat sea should not be significantly larger than 10. Otherwise, there will be too little evaporation over the oceans. | {math}`10.0 \pm 8.0` |  additive perturbation ||
| rat_sea   | controls latent and heat fluxes over water | 0.8 - 10.0 | lower values increase latent and sensible fluxes over water |
| a_hshr    | length scale factor for the separated horizontal shear mode | {math}`1.0 \pm 1.0` ||
| a_stab    | factor for stability correction of turbulent length scale |   {math}`0.0 \pm 1.0` |  turned off by default because it degrades global skill scores |
| c_diff    | length scale factor for vertical diffusion of TKE | {math}`0.2 \pm 2.0` ||
| alpha0 | lower bound of velocity-dependent Charnock parameter | 0.0123-0.0335 | additive ensemble perturbation of Charnock-parameter |
| alpha1 | parameter scaling the molecular roughness of water waves | 0.1-1.0 | lower values increase latent and sensible fluxes over water, particularly at low wind speeds.|
| tur_len | asymptotic maximal turbulent distance | {math}`500. \pm 150.` [m] | default is 150 m |
| tkhmin | scaling factor for minimum vertical diffusion coefficient for heat and moisture | {math}`0.75 \pm 0.2` | 0.75 |
| tkmmin | scaling factor for minimum vertical diffusion coefficient for momentum | {math}`0.75 \pm 0.2`   ||
| tkred_sfc |       multiplicative change of reduction of minimum diffusion coefficients near the surface | 0.25 - 4.0 |   multiplicative perturbation |

## TERRA

|  Parameter  |      Description |             Meaningful  Range |  Comment |
|  --------|--------|--------------|-----|
| c_soil  |  evaporating fraction of soil | {math}`1.0 \pm 0.25` ||
| cwimax_ml | scaling parameter for maximum interception storage |  {math}`5.\cdot 10^{-7} -  5.\cdot 10^{-4}` |  multiplicative perturbation, low values ({math}`< 10^{-6}`) turn off interception layer |
| minsnowfrac | Lower limit of snow cover fraction to which melting snow is artificially reduced in the context of the snow-tile approach | {math}`0.2 \pm 0.1` ||
| dust_abs | Tuning factor for enhanced LW absorption of mineral dust in the Saharan region |     0.0  | Reduces bias over Sahara for the RRTM scheme but not necessary and implemented with ecRad and itype_lwemiss=2 |

## convection

|  Parameter  |      Description |           Meaningful  Range    |  Comment |
|  --------|--------|--------------|-----|
| entrorg  | entrainment parameter in convection scheme valid for dx=20km | {math}`1.95\cdot 10^{-3}\pm 0.2\cdot 10^{-3}` |  corresponds to entr_sc in the shallow convection part of COSMO Tiedtke scheme |
| rdepths  | maximum allowed shallow  convection depth | {math}`2.0 \cdot 10^{4}` {math}`\pm 5.0\cdot 10^{3}` Pa | |
| rprcon   | coefficient for conversion of cloud water into precipitation | {math}`1.4\cdot 10^{-3}\pm 0.25\cdot 10^{-3}` ||
| capdcfac_et | fraction of CAPE diurnal cycle correction applied in the extratropics | {math}`0.5 \pm 0.7` ||
| capdcfac_tr | fraction of CAPE diurnal cycle correction applied in the tropics | {math}`0.5 \pm 0.75` ||
| lowcapefac  | tuning parameter for diurnal-cycle correction in convection scheme: reduction factor for low-cape situations | {math}`1.0 \pm 0.5` ||
| negpblcape  | tuning parameter for diurnal-cycle correction in convection scheme: maximum negative PBL CAPE allowed in the modified CAPE closure | -500.- 0.||
| rhebc_land  | RH threshold for onset of evaporation below cloud base over land | {math}`0.825 \pm 0.05` | 0.75 as default in code |
| rhebc_ocean | RH threshold for onset of evaporation below cloud base over sea | {math}`0.85 \pm 0.05` ||
| rhebc_land_trop | RH threshold ... over tropical land |   {math}`0.70 \pm 0.05`  | tropics ||
| rhebc_ocean_trop | RH threshold ...over tropical sea  |  {math}`0.76 \pm 0.05`   | tropics ||
| rcucov | convective area fraction used for computing evaporation below cloud base | 0.075 | 0.05 coded as default ||
| rcucov_trop | convective area fraction used for computing evaporation below cloud base, tropics | 0.03 | tropics ||
| texc | Excess value for temperature used in test parcel ascent | {math}`0.125 \pm 0.05` [K] ||
| qexc | Excess fraction of grid-scale QV used in test parcel ascent | {math}`0.0125 \pm 0.005` [kg/kg] ||
