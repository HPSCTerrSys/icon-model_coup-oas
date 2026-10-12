```{eval-rst}
:orphan:
```

(ref_infrastructure_deprecated)=
# List of Deprecated Features

<a href="#table_deprecated">Table 1</a> provides an overview of features to be removed in future releases.
Further details on these deprecated features are provided below.

If you want to activate a deprecated feature anyways, add `--enable-deprecated` to your [configuration](ref_buildrun_configuration).
Note, that you have to re-configure in this case.

{material-regular}`warning;2em;pst-color-secondary` _Please note that deprecated features are only minimally tested to keep them basically running until final removal._
<a name="table_deprecated">

:::{table} Table 1: List of deprecated features
:width: 65
:widths: auto
:align: center

| Option                                                     | Scheduled for Removal    |
| :-----------                                               | :------------            |
| [RRTM radiation (NWP)](ref_infrastructure_deprecated_rrtm) | icon-2027.04             |
| [GME Ozone (NWP)](ref_infrastructure_deprecated_ozone6)    | icon-2027.04             |
:::

(ref_infrastructure_deprecated_rrtm)=
## RRTM radiation scheme (NWP physics)

Description:
: Old operational radiation scheme derived from ECHAM/IFS

Associated Namelist Options:
: `inwp_radiation=1`  (`&nwp_phy_nml`)
: `lrtm_filename` (`&nwp_phy_nml`)
: `cldopt_filename` (`&nwp_phy_nml`)

Reasons:
: Unmaintained legacy code
: Known, unresolved bugs (zigzag patterns in long wave flux divergences)
: Many compiler warnings
: Out of operational use for nearly 5 years, not recommended
: Recent radiation input options are only available to alternative schemes (reduced maintenance overhead)

Recommended Alternatives:
: [ecRad](ref_atmosphere_ecrad): `inwp_radiation=4`

Scheduled for Removal
: icon-2027.04

(ref_infrastructure_deprecated_ozone6)=
## Ozone Option 6 (GME ozone) (NWP physics)

Description:
: Ozone climatology with T5 geographical distribution and Fourier series for seasonal cycle

Associated Namelist Options:
: `irad_o3=6`  (`&radiation_nml`)

Reasons:
: Was only available for `inwp_radiation=1`
: Not maintained for more than 10 years

Recommended Alternatives:
: Blending between GEMS and MACC ozone climatologies, `irad_o3=79` (`&radiation_nml`)

Scheduled for Removal
: icon-2027.04
