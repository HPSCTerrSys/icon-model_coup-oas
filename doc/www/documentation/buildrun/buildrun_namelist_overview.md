```{eval-rst}
:orphan:
```

(ref_buildrun_nml)=
# ICON Namelists

All available Fortran namelist parameters are tabulated by name, type, default value, unit, description, and scope:

 - _Type_ refers to the type of the Fortran variable, in which the
value is stored: I=INTEGER, L=LOGICAL, R=REAL, C=character string
 - _Default_ is the preset value, if defined, that is assigned to
this parameter within the programs.
 - _Unit_ shows the unit of the control parameter, where applicable.
 - _Description_ explains in a few words the purpose of the parameter.
 - _Scope_ explains under which conditions the namelist parameter
has any effect, if its scope is restricted to specific settings of
other namelist parameters.

Information on the file, where the namelist is defined and used, is
given at the end of each table.



# Namelist parameters defining the atmospheric model

Namelist parameters for the ICON models are organized in several thematic
Fortran namelists controling the experiment, and the properties of
dynamics, transport, physics etc.



(ref_buildrun_nml_aes_bubble_nml)=
## aes_bubble_nml

The following namelist controls the parameter setting for the testcase
'aes_bubble'. In the framework of this testcase, particular initial
conditions can be set by the parameters described in the table of the
namelist variables hereafter:


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (aes_bubble_nml-psfc)=
    **psfc**
  - R
  - 101325.0
  - Pa
  - Initial value of surface pressure.
  -

* - (aes_bubble_nml-t_am)=
    **t_am**
  - R
  - 180.0
  - K
  - Absolute minimum of atmospheric temperature in initial state.
  -

* - (aes_bubble_nml-t0)=
    **t0**
  - R
  - 303.5
  - K
  - Temperature at bottom of atmosphere in initial state.
  -

* - (aes_bubble_nml-gamma0)=
    **gamma0**
  - R
  - 0.009
  - K/m
  - Lapse rate in lowest atmospheric part in initial temperature profile.
  -

* - (aes_bubble_nml-z0)=
    **z0**
  - R
  - 3000.0
  - m
  - Below z0, the lapse rate gamma0 is applied in the initial temperature profile, above z0, the lapse rate is gamma1.
  -

* - (aes_bubble_nml-gamma1)=
    **gamma1**
  - R
  - 0.00001
  - K/m
  - Lapse rate above z0 in the initial temperature profile. However, temperature cannot fall below t_am in the initial temperature profile.
  -

* - (aes_bubble_nml-t_perturb)=
    **t_perturb**
  - R
  - 3.0
  - K
  - Maximum temperature perturbation in center of Gaussians in initial state.
  -

* - (aes_bubble_nml-relhum_bg)=
    **relhum_bg**
  - R
  - 0.7
  -
  - Background relative humidity in initial state.
  -

* - (aes_bubble_nml-relhum_mx)=
    **relhum_mx**
  - R
  - 0.95
  -
  - Maximum relative humidity in initial state.
  -

* - (aes_bubble_nml-hw_x)=
    **hw_x**
  - R
  - 12500.0
  - m
  - Half width in x-direction in meters of the bubble in initial state.
  -

* - (aes_bubble_nml-hw_z)=
    **hw_z**
  - R
  - 500.0
  - m
  - Half width in z-direction in meters of the bubble in initial state.
  -

* - (aes_bubble_nml-x_center)=
    **x_center**
  - R
  - 0.0
  - m
  - Placement of maximum of Gaussian relative to the origin in x-direction (if Gaussian is applied into x-direction only, lgaussxy=.FALSE.) or relative to the origin in x- and y-direction (if Gaussian is applied into x- and y- direction, lgaussxy=.TRUE.) in initial state.
  -

* - (aes_bubble_nml-lgaussxy)=
    **lgaussxy**
  - L
  - .FALSE.
  - K
  - .TRUE., if half width calculated for x-direction and x_center is applied also to y direction in initial state.
  -
```




(ref_buildrun_nml_aes_cop_nml)=
## aes_cop_nml

The parameterization of cloud optical properties for the AES physics is configured
by a data structure _aes_cop_config(jg=1:ndom)%\<param>_, which is a
1-dimensional array extending over all  domains. The structure contains parameters `i`
providing control over the parametrized effects:


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (aes_cop_nml-cn1lnd)=
    **cn1lnd**
  - R
  - 20.0
  - 1e6/m3
  - cloud droplet number concentration over land
  -

* - (aes_cop_nml-cn2lnd)=
    **cn2lnd**
  - R
  - 180.0
  - 1e6/m3
  - cloud droplet number concentration over land
  -

* - (aes_cop_nml-cn1sea)=
    **cn1sea**
  - R
  - 20.0
  - 1e6/m3
  - cloud droplet number concentration over sea
  -

* - (aes_cop_nml-cn2sea)=
    **cn2sea**
  - R
  - 80.0
  - 1e6/m3
  - cloud droplet number concentration over sea
  -

* - (aes_cop_nml-cinhomi)=
    **cinhomi**
  - R
  - 0.8
  -
  - ice cloud inhomogeneity factor
  -

* - (aes_cop_nml-cinhoml_cfl)=
    **cinhoml_cfl**
  - R
  - 0.4
  -
  - liquid cloud inhomogeneity factor for cumuliform clouds over land
  -

* - (aes_cop_nml-cinhoml_cfo)=
    **cinhoml_cfo**
  - R
  - 0.4
  -
  - liquid cloud inhomogeneity factor for cumuliform clouds over ocean
  -

* - (aes_cop_nml-cinhoml_sf)=
    **cinhoml_sf**
  - R
  - 0.8
  -
  - liquid cloud inhomogeneity factor for stratiform clouds
  -

* - (aes_cop_nml-cinhoml_del1)=
    **cinhoml_del1**
  - R
  - 2.0
  -
  - del1 is the ordinate value of the atan2 function blending between cumuliform and stratiform liquid cloud inhomogeneity factors depending on the lower tropospheric stability (lts)
  -

* - (aes_cop_nml-cinhoml_del2)=
    **cinhoml_del2**
  - R
  - 20.0
  -
  - lts - del2 is the abscissa value of the atan2 function blending between cumuliform and stratiform liquid cloud inhomogeneity factors, depending on the lower tropospheric stability (lts)
  -

* - (aes_cop_nml-cinhoml_lts_height)=
    **cinhoml_lts_height**
  - R
  - 3200.0
  - m
  - height over ocean used to determine the index of the highest level below this height. This model level is used to compute the lower tropospheric stability as a difference of the potential temperature at this level and at the surface to identify stratocumulus areas
  -

* - (aes_cop_nml-cthomi)=
    **cthomi**
  - R
  - tmelt-35.
  - K
  - maximum temperature for homogeneous freezing
  -

* - (aes_cop_nml-csecfrl)=
    **csecfrl**
  - R
  - 1.5E-5
  - kg/kg
  - minimum in-cloud water mass mixing ratio in mixed phase clouds
  -
```




(ref_buildrun_nml_aes_cov_nml)=
## aes_cov_nml

The parameterization of cloud cover for the AES physics is configured
by a data structure _aes_cov_config(jg=1:ndom)%\<param>_, which is a
1-dimensional array extending over all domains. The structure contains the
following control parameters:


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (aes_cov_nml-cqx)=
    **cqx**
  - R
  - 1.0e-8
  - kg/kg
  - critical mass fraction of cloud water + cloud ice in air, if exceeded cloud cover in cell = 1, otherwise = 0
  -
```




(ref_buildrun_nml_aes_phy_nml)=
## aes_phy_nml

The AES physics is configured by a data structure _aes_phy_config(jg=1:ndom)%\<param>_,
which is a 1-dimensional array extending over all domains. The structure contains
several parameters providing time control for the atmospheric forcing by the different
parameterizations. Further logical switches control how the atmospheric boundary conditions
for the AES physics are determined. Time control parameters are available for the atmospheric
processes tabulated below.


```{list-table}
:header-rows: 1

* - _prc_
  - _parameterized process_

* - rad
  - LW and SW radiation

* - vdf
  - vertical diffusion

* - mig
  - graupel microphysics

* - two
  - two moment microphysics

* - car
  - Cariolle's linearized ozone chemistry

* - art
  - ART chemistry
```



The time control for an atmospheric forcing by a process _prc_ consists of
three components, the time interval _dt_prc_ for re-computing the forcing,
and the start and end dates and times defining the interval _[sd_prc,ed_prc[_,
in which the forcing is either computed, if the date/time coincides with the interval
_dt_prc_, or recycled. Recycling means that the forcing stored from the last
computation is used again. Outside of the interval the forcing is set to zero.

If _dt_prc_ is not specified, or an empty string or a string of blanks or an
interval of length 0s, e.g. _"PT0S"_ is given, then the forcing is switched
off for the entire experiment and the start and end dates and times are irrelevant.

If _sd_prc_ or _ed_prc_ are not specified, or an empty string or a
string of blanks are given, then the experiment start date and the experiment
stop date are used, respectively.

Further the forcing control switch _fc_prc_ can be used to decide if an
active process (_dt_prc_ > 0) is used for the integration
(_fc_prc_ = 1) or only computed for diagnostic purposes (_fc_prc_ = 0).


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (aes_phy_nml-dt_prc)=
    **dt_prc**
  - C
  - ""
  -
  - This is the time interval in ISO 8601-2004 format at which the forcing by the process _prc_ is computed.
  - [iforcing](run_nml-iforcing) = 2

* - (aes_phy_nml-sd_prc)=
    **sd_prc**
  - C
  - ""
  -
  - Defines the start date/time in ISO 8601-2004 format of the interval _[sd_prc,ed_prc[_, in which the forcing by the process _prc_ is computed in intervals _dt_prc_.
  - [iforcing](run_nml-iforcing) = 2 and _dt_prc_  > 0.000s

* - (aes_phy_nml-ed_prc)=
    **ed_prc**
  - C
  - ""
  -
  - Defines the end date/time in ISO 8601-2004 format of the interval _[sd_prc,ed_prc[_, in which the forcing by the process _prc_ is computed in intervals _dt_prc_.
  - [iforcing](run_nml-iforcing) = 2 and _dt_prc_  > 0.000s

* - (aes_phy_nml-fc_prc)=
    **fc_prc**
  - I
  - 1
  -
  - Forcing control for process  _prc_.
    - fc_prc = 0: the forcing of the process is not used in the integration.
    - fc_prc = 1: the forcing of the process is used in the integration.
  - [iforcing](run_nml-iforcing) = 2 and _dt_prc_  > 0.000s

* - (aes_phy_nml-ljsb)=
    **ljsb**
  - L
  - .FALSE.
  -
  - .TRUE. for using the JSBACH land surface model
  - [iforcing](run_nml-iforcing) = 2

* - (aes_phy_nml-llake)=
    **llake**
  - L
  - .FALSE.
  -
  - .TRUE. for using lakes in JSBACH
  - [iforcing](run_nml-iforcing) = 2

* - (aes_phy_nml-lamip)=
    **lamip**
  - L
  - .FALSE.
  -
  - .TRUE. for AMIP boundary conditions
  - [iforcing](run_nml-iforcing) = 2

* - (aes_phy_nml-l2moment)=
    **l2moment**
  - L
  - .FALSE.
  -
  - .TRUE. for the 2-moment microphysics scheme
  - [iforcing](run_nml-iforcing) = 2

* - (aes_phy_nml-lmlo)=
    **lmlo**
  - L
  - .FALSE.
  -
  - .TRUE. for mixed layer ocean
  - [iforcing](run_nml-iforcing) = 2

* - (aes_phy_nml-lice)=
    **lice**
  - L
  - .FALSE.
  -
  - .TRUE. for sea-ice temperature calculation
  - [iforcing](run_nml-iforcing) = 2

* - (aes_phy_nml-lsstice)=
    **lsstice**
  - L
  - .FALSE.
  -
  - .TRUE. for inst. 6hourly sst and ice (prelim)
  - [iforcing](run_nml-iforcing) = 2

* - (aes_phy_nml-iqneg_d2p)=
    **iqneg_d2p**
  - I
  - 0
  -
  - If negative tracer mass fractions are found in the dynamics to physics interface, then:
    - 1,3: they are reported;
    - 2,3: they are replaced with zero
  - [iforcing](run_nml-iforcing) = 2

* - (aes_phy_nml-iqneg_p2d)=
    **iqneg_p2d**
  - I
  - 0
  -
  - If negative tracer mass fractions are found in the physics to dynamics interface, then:
    1,3: they are reported;
    2,3: they are replaced with zero
  - [iforcing](run_nml-iforcing) = 2

* - (aes_phy_nml-zmaxcloudy)=
    **zmaxcloudy**
  - R
  - 22500.0
  - m
  - maximum height (m) for cloud related computations
  -
```




(ref_buildrun_nml_aes_rad_nml)=
## aes_rad_nml

The input from AES physics to the rte_rrtmgp scheme is configured by
a data structure _aes_rad_config(jg=1:ndom)%\<param>_, which
is a 1-dimensional array extending over all  domains. The structure
contains parameters providing control over the Earth orbit, the
computation of the SW incoming flux at the top of the atmosphere and
the atmospheric composition:


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (aes_rad_nml-isolrad)=
    **isolrad**
  - I
  - 0
  -
  - Selects the spectral solar irradiation (SSI) at 1 AU distance from the sun
    - 0: SRTM default solar spectrum, TSI = 1368.222 W/m2.
    - 1: Time dependent solar sprectrum from file
    - 2: Average 1844--1856 of transient CMIP5 solar, TSI = 1360.875 W/m2
    - 3: Average 1979--1988 of transient CMIP5 solar spectrum, TSI = 1361.371 W/m2
    - 4: Solar flux for RCE simulations with diurnal cycle, TSI = 1069.315 W/m2
    - 5: Solar flux for RCE simulations without diurnal cycle, TSI = 433.3371 W/m2
    - 6: Average 1850-1873 of transient CMIP6 solar, TSI = 1360.744 W/m2
    - 7: Solar flux for RCEmip analytical simulations without diurnal cycle, TSI = 551.58 W/m2
  - dt_rad > 0.000s

* - (aes_rad_nml-fsolrad)=
    **fsolrad**
  - R
  - 1
  -
  - Scaling factor for solar irradiance
  - dt_rad > 0.000s

* - (aes_rad_nml-l_orbvsop87)=
    **l_orbvsop87**
  - L
  - .TRUE.
  -
  - .TRUE. for the realistic VSOP87 Earth orbit
    .FALSE. for the Kepler orbit
  - dt_rad > 0.000s

* - (aes_rad_nml-cecc)=
    **cecc**
  - R
  - 0.016715
  -
  - eccentricity of the Kepler orbit
  - dt_rad > 0.000s and l_orbvsop87 = .FALSE.

* - (aes_rad_nml-cobld)=
    **cobld**
  - R
  - 23.44100
  - deg
  - obliquity of the Earth rotation axis on the Kepler orbit
  - dt_rad > 0.000s  and l_orbvsop87 = .FALSE.

* - (aes_rad_nml-clonp)=
    **clonp**
  - R
  - 282.7000
  - deg
  - longitude of perihelion with respect to vernal equinox on the Kepler orbit
  - dt_rad > 0.000s  and l_orbvsop87 = .FALSE.

* - (aes_rad_nml-lyr_perp)=
    **lyr_perp**
  - L
  - .FALSE.
  -
  - .FALSE. for transient VSOP87 Earth orbit
    .TRUE.: VSOP87 Earth orbit of year yr_perp is perpertuated
  - dt_rad > 0.000s  and l_orbvsop87 = .TRUE.

* - (aes_rad_nml-yr_perp)=
    **yr_perp**
  - L
  - -99999
  -
  - year of vsop87 orbit to  be perpetuated for lyr_perp = .TRUE.
  - dt_rad > 0.000s  and l_orbvsop87 = .TRUE.

* - (aes_rad_nml-ldiur)=
    **ldiur**
  - L
  - .TRUE.
  -
  - .TRUE. for diurnal cycle in solar irradiation
    .FALSE. for zonally averaged solar irradiation
  - dt_rad > 0.000s

* - (aes_rad_nml-l_sph_symm_irr)=
    **l_sph_symm_irr**
  - L
  - .FALSE.
  -
  - .TRUE. for globally averaged irradiation; .FALSE. for lat (lon) dependent irradiation
  -

* - (aes_rad_nml-icosmu0)=
    **icosmu0**
  - I
  - 3
  -
  - PROVISIONAL - ONLY BEST METHODS WILL BE KEPT ("0" or "3")
    - 0: no adjustment, the original cosmu0 is used
    - 3: 0.5*SIN(dmu0)*(1+(pi/2-mu0)/dmu0), dmu0=pi*dt_rad/86400s
    Has small effects on the MA temp. and wind and the land surface temp.
  - dt_rad > 0.000s

* - (aes_rad_nml-irad_h2o)=
    **irad_h2o**
  - I
  - 1
  -
  - Selects source for concentration of water vapor, cloud water and cloud ice
    - 0: No H2O(gas,liq,ice) in radiation
    - 1: H2O(gas,liq,ice) mass mixing ratios from tracer fields
  - dt_rad > 0.000s

* - (aes_rad_nml-irad_co2)=
    **irad_co2**
  - I
  - 2
  -
  - Selects source for concentration of CO2
    - 0: No CO2 in radiation
    - 1: CO2 mass mixing ratio from tracer field
    - 2: CO2 volume mixing ratio set by [vmr_co2](aes_rad_nml-vmr_co2)
    - 3: CO2 volume mixing ratio from ghg scenario file
  - dt_rad > 0.000s and CO2 tracer is defined

* - (aes_rad_nml-irad_ch4)=
    **irad_ch4**
  - I
  - 2
  -
  - Selects source for concentration of CH4
    - 0: No CH4 in radiation
    - 2: CH4 volume mixing ratio set by 'vmr _ch4'
    - 3: CH4 vertically constant volume mixing ratio from ghg scenario file
    - 12: CH4 tanh-profile with surface volume mixing ratio set by 'vmr _ch4'
    - 13: CH4 tanh-profile with surface volume mixing ratio from ghg scenario file
  - dt_rad > 0.000s

* - (aes_rad_nml-irad_n2o)=
    **irad_n2o**
  - I
  - 2
  -
  - Selects source for concentration of N2O
    - 0: No N2O in radiation
    - 2: N2O volume mixing ratio set by [vmr_n2o](aes_rad_nml-vmr_ch4)
    - 3: N2O vertically constant volume mixing ratio from ghg scenario file
    - 12: N2O tanh-profile with surface volume mixing ratio set by [vmr_n2o](aes_rad_nml-vmr_ch4)
    - 13: N2O tanh-profile with surface volume mixing ratio from ghg scenario file
  - dt_rad > 0.000s

* - (aes_rad_nml-irad_cfc11)=
    **irad_cfc11**
  - I
  - 2
  -
  - Selects source for concentration of CFC11
    - 0: No CFC11 in radiation
    - 2: CFC11 volume mixing ratio set by 'vmr _cfc11'
    - 3: CFC11 volume mixing ration from ghg scenario file
  - dt_rad > 0.000s

* - (aes_rad_nml-irad_cfc12)=
    **irad_cfc12**
  - I
  - 2
  -
  - Selects source for concentration of CFC12
    - 0: No CFC12 in radiation
    - 2: CFC12 volume mixing ratio set by 'vmr _cfc12'
    - 3: CFC12 volume mixing ration from ghg scenario file
  - dt_rad > 0.000s

* - (aes_rad_nml-irad_o3)=
    **irad_o3**
  - I
  - 0
  -
  - Selects source for concentration of O3
    - 0: No O3 in radiation
    - 1: O3 mass mixing ratio from tracer field
    - 4: O3 constant-in-time 3-dim. volume mixing ratio from file
    - 5: O3 transient 3-dim. volume mixing ratio from file
    - 6: O3 clim. annual cycle 3-dim. volume mixing ratio from file
  - dt_rad > 0.000s

* - (aes_rad_nml-irad_o2)=
    **irad_o2**
  - I
  - 2
  -
  - Selects source for concentration of O2
    - 0: No O2 in radiation
    - 2: O2 volume mixing ratio set by 'vmr _o2'
  - dt_rad > 0.000s

* - (aes_rad_nml-vmr_co2)=
    **vmr_co2**
  - R
  - 348.0e-06
  - m3/m3
  - Volume mixing ratio of CO2
  - dt_rad > 0.000s

* - (aes_rad_nml-vmr_ch4)=
    **vmr_ch4**
  - R
  - 1650.0e-09
  - m3/m3
  - Volume mixing ratio of CH4
  - dt_rad > 0.000s

* - (aes_rad_nml-vmr_n2o)=
    **vmr_n2o**
  - R
  - 306.0e-09
  - m3/m3
  - Volume mixing ratio of N2O
  - dt_rad > 0.000s

* - (aes_rad_nml-vmr_o2)=
    **vmr_o2**
  - R
  - 0.20946
  - m3/m3
  - Volume mixing ratio of O2
  - dt_rad > 0.000s

* - (aes_rad_nml-vmr_cfc11)=
    **vmr_cfc11**
  - R
  - 214.5e-12
  - m3/m3
  - Volume mixing ratio of CFC11
  - dt_rad > 0.000s

* - (aes_rad_nml-vmr_cfc12)=
    **vmr_cfc12**
  - R
  - 371.1e-12
  - m3/m3
  - Volume mixing ratio of CFC11
  - dt_rad > 0.000s

* - (aes_rad_nml-frad_h2o)=
    **frad_h2o**
  - R
  - 1.0
  -
  - Scaling factor for concentration of water vapor, cloud water and cloud ice
  - dt_rad > 0.000s

* - (aes_rad_nml-frad_co2)=
    **frad_co2**
  - R
  - 1.0
  -
  - Scaling factor for concentration of  CO2
  - dt_rad > 0.000s

* - (aes_rad_nml-frad_ch4)=
    **frad_ch4**
  - R
  - 1.0
  -
  - Scaling factor for concentration of CH4
  - dt_rad > 0.000s

* - (aes_rad_nml-frad_n2o)=
    **frad_n2o**
  - R
  - 1.0
  -
  - Scaling factor for concentration of N2O
  - dt_rad > 0.000s

* - (aes_rad_nml-frad_o3)=
    **frad_o3**
  - R
  - 1.0
  -
  - Scaling factor for concentration of O3
  - dt_rad > 0.000s

* - (aes_rad_nml-frad_o2)=
    **frad_o2**
  - R
  - 1.0
  -
  - Scaling factor for concentration of O2
  - dt_rad > 0.000s

* - (aes_rad_nml-frad_cfc11)=
    **frad_cfc11**
  - R
  - 1.0
  -
  - Scaling factor for concentration of CFC11
  - dt_rad > 0.000s

* - (aes_rad_nml-frad_cfc12)=
    **frad_cfc12**
  - R
  - 1.0
  -
  - Scaling factor for concentration of CFC12
  - dt_rad > 0.000s

* - (aes_rad_nml-irad_aero)=
    **irad_aero**
  - I
  - 2
  -
  - Selects source of aerosol types
    - 0: No aerosols in radiation
    - 12: tropospheric background aerosol (Kinne). The filename linked to the file containing the background aerosols data must not contain a year.
    - 13: transient tropospheric aerosol  (Kinne) including anthropogenic aerosols
    - 18: tropospheric background aerosol (Kinne, filename without year) + stratospheric aerosols + simple plumes (parameterized time dependent anthropogenic aerosols)
    - 19: tropospheric background aerosol (Kinne, filename without year), no stratospheric aerosols + simple plumes (parameterized time dependent anthropogenic aerosols)
  - dt_rad > 0.000s
```






(ref_buildrun_nml_aes_vdf_nml)=
## aes_vdf_nml


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (aes_vdf_nml-lsfc_mom_flux)=
    **lsfc_mom_flux**
  - L
  - .TRUE.
  -
  - switch on/off surface momentum flux
  - dt_vdf > 0.000s

* - (aes_vdf_nml-lsfc_heat_flux)=
    **lsfc_heat_flux**
  - L
  - .TRUE.
  -
  - switch on/off surface heat flux
  - dt_vdf > 0.000s

* - (aes_vdf_nml-pr0)=
    **pr0**
  - R
  - 1.0
  -
  - neutral limit Prandtl number, can be varied from about 0.6 to 1.0
  - dt_vdf > 0.000s

* - (aes_vdf_nml-f_tau0)=
    **f_tau0**
  - R
  - 0.17
  -
  - neutral non-dimensional stress factor
  - dt_vdf > 0.000s

* - (aes_vdf_nml-c_f)=
    **c_f**
  - R
  - 0.185
  -
  - mixing length: coriolis term tuning parameter
  - dt_vdf > 0.000s

* - (aes_vdf_nml-c_n)=
    **c_n**
  - R
  - 2.0
  -
  - mixing length: stability term tuning parameter
  - dt_vdf > 0.000s

* - (aes_vdf_nml-wmc)=
    **wmc**
  - R
  - 0.5
  -
  - ratio of typical horizontal velocity to wstar at free convection
  - dt_vdf > 0.000s

* - (aes_vdf_nml-fsl)=
    **fsl**
  - R
  - 0.4
  -
  - fraction of first-level height at which surface fluxes are nominally evaluated, tuning param for sfc stress
  - dt_vdf > 0.000s

* - (aes_vdf_nml-fbl)=
    **fbl**
  - R
  - 3.0
  -
  - 1/fbl: fraction of BL height at which lmix hat its max
  - dt_vdf > 0.000s

* - (aes_vdf_nml-lmix_max)=
    **lmix_max**
  - R
  - 150.0
  - m
  - maximum mixing length
  - dt_vdf > 0.000s

* - (aes_vdf_nml-z0m_min)=
    **z0m_min**
  - R
  - 0.000015
  - m
  - minimum roughness length
  - dt_vdf > 0.000s

* - (aes_vdf_nml-z0m_ice)=
    **z0m_ice**
  - R
  - 0.001
  - m
  - roughness length for sea ice surfaces
  - dt_vdf > 0.000s

* - (aes_vdf_nml-z0m_oce)=
    **z0m_oce**
  - R
  - 0.001
  - m
  - roughness length for sea water surfaces
  - dt_vdf > 0.000s

* - (aes_vdf_nml-turb)=
    **turb**
  - I
  - 1
  -
  - 1: TTE scheme, 2: 3D Smagorinsky
  - dt_vdf > 0.000s

* - (aes_vdf_nml-smag_constant)=
    **smag_constant**
  - R
  - 0.23
  -
  -
  - dt_vdf > 0.000s

* - (aes_vdf_nml-max_turb_scale)=
    **max_turb_scale**
  - R
  - 300.0
  -
  - max turbulence length scale
  - dt_vdf > 0.000s

* - (aes_vdf_nml-turb_prandtl)=
    **turb_prandtl**
  - R
  - 0.33333333333
  -
  - turbulent prandtl number
  - dt_vdf > 0.000s

* - (aes_vdf_nml-km_min)=
    **km_min**
  - R
  - 0.001
  -
  - min mass weighted turbulent viscosity
  - dt_vdf > 0.000s

* - (aes_vdf_nml-min_sfc_wind)=
    **min_sfc_wind**
  - R
  - 0.3
  -
  - min sfc wind in free convection limit
  - dt_vdf > 0.000s

* - (aes_vdf_nml-wind_g)=
    **wind_g**
  - R
  - 3.0
  -
  - wind gust parameter
  - dt_vdf > 0.000s

* - (aes_vdf_nml-lcuda_graph_vdf)=
    **lcuda_graph_vdf**
  - L
  - .FALSE.
  -
  - Activate cuda graphs for JSBACH land model. Automatically set to .FALSE. if not compiled with the ICON_USE_CUDA_GRAPH cpp key.
  - ICON_USE_CUDA_GRAPH activated
```






(ref_buildrun_nml_aes_wmo_nml)=
## aes_wmo_nml

The diagnostics of the tropopause pressure, following the WMO definition is configured by a data structure
_aes_wmo_config(jg=1:ndom)%\<param>_, which is a 1-dimensional array extending over all  domains:


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (aes_wmo_nml-zmaxwmo)=
    **zmaxwmo**
  - R
  - 38000.0
  - m
  - maximum height for tropopause search
  -

* - (aes_wmo_nml-zminwmo)=
    **zminwmo**
  - R
  - 5000.0
  - m
  - minimum height for tropopause search
  -
```




(ref_buildrun_nml_art_nml)=
## art_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (art_nml-cart_aero_emiss_xml)=
    **cart_aero_emiss_xml**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_aerosol_xml)=
    **cart_aerosol_xml**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_cheminit_coord)=
    **cart_cheminit_coord**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_cheminit_file)=
    **cart_cheminit_file**
  - CHARACTER
  -
  -
  -
  -

* - (art_nml-cart_cheminit_type)=
    **cart_cheminit_type**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_chemtracer_xml)=
    **cart_chemtracer_xml**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_coag_xml)=
    **cart_coag_xml**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_diagnostics_xml)=
    **cart_diagnostics_xml**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_emiss_xml_file)=
    **cart_emiss_xml_file**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_ext_data_xml)=
    **cart_ext_data_xml**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_fplume_inp)=
    **cart_fplume_inp**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_input_folder)=
    **cart_input_folder**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_io_suffix)=
    **cart_io_suffix**
  - CHARACTER
  -
  -
  -
  -

* - (art_nml-cart_mecca_xml)=
    **cart_mecca_xml**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_modes_xml)=
    **cart_modes_xml**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_opt_props_nc)=
    **cart_opt_props_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_pntSrc_xml)=
    **cart_pntSrc_xml**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_radioact_file)=
    **cart_radioact_file**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_type_sedim)=
    **cart_type_sedim**
  - CHARACTER
  - "expl"
  -
  -
  -

* - (art_nml-cart_volcano_file)=
    **cart_volcano_file**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-cart_vortex_init_date)=
    **cart_vortex_init_date**
  - CHARACTER
  - ''
  -
  -
  -

* - (art_nml-iart_aci_cold)=
    **iart_aci_cold**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_aci_warm)=
    **iart_aci_warm**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_aero_washout)=
    **iart_aero_washout**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_anthro)=
    **iart_anthro**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_ari)=
    **iart_ari**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_dust)=
    **iart_dust**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_fire)=
    **iart_fire**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_fplume)=
    **iart_fplume**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_init_aero)=
    **iart_init_aero**
  - INTEGER
  -
  -
  -
  -

* - (art_nml-iart_init_gas)=
    **iart_init_gas**
  - INTEGER
  -
  -
  -
  -

* - (art_nml-iart_isorropia)=
    **iart_isorropia**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_modeshift)=
    **iart_modeshift**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_nonsph)=
    **iart_nonsph**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_pollen)=
    **iart_pollen**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_radioact)=
    **iart_radioact**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_seas_water)=
    **iart_seas_water**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_seasalt)=
    **iart_seasalt**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_volc_numb)=
    **iart_volc_numb**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-iart_volcano)=
    **iart_volcano**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-irad_multicall)=
    **irad_multicall**
  - INTEGER
  - 0
  -
  -
  -

* - (art_nml-lart_aerosol)=
    **lart_aerosol**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (art_nml-lart_chem)=
    **lart_chem**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (art_nml-lart_chemtracer)=
    **lart_chemtracer**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (art_nml-lart_conv)=
    **lart_conv**
  - LOGICAL
  - .TRUE.
  -
  -
  -

* - (art_nml-lart_debugRestart)=
    **lart_debugRestart**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (art_nml-lart_diag_out)=
    **lart_diag_out**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (art_nml-lart_diag_xml)=
    **lart_diag_xml**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (art_nml-lart_dusty_cirrus)=
    **lart_dusty_cirrus**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (art_nml-lart_emiss_turbdiff)=
    **lart_emiss_turbdiff**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (art_nml-lart_excl_end_pntSrc)=
    **lart_excl_end_pntSrc**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (art_nml-lart_mecca)=
    **lart_mecca**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (art_nml-lart_pntSrc)=
    **lart_pntSrc**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (art_nml-lart_psc)=
    **lart_psc**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (art_nml-lart_turb)=
    **lart_turb**
  - LOGICAL
  - .TRUE.
  -
  -
  -

* - (art_nml-nart_substeps_sedi)=
    **nart_substeps_sedi**
  - INTEGER
  -
  -
  -
  -

* - (art_nml-radioact_maxtint)=
    **radioact_maxtint**
  - REAL
  -
  -
  -
  -

* - (art_nml-rart_dustyci_crit)=
    **rart_dustyci_crit**
  - REAL
  - 70.0_wp
  -
  -
  -

* - (art_nml-rart_dustyci_rhi)=
    **rart_dustyci_rhi**
  - REAL
  - 0.90_wp
  -
  -
  -

* - (art_nml-rart_qv_fire)=
    **rart_qv_fire**
  - REAL
  - 0.0_wp
  -
  -
  -

* - (art_nml-rart_shfl_fire)=
    **rart_shfl_fire**
  - REAL
  - 0.0_wp
  -
  -
  -
```


Defined and used in: {{ '[src/namelists/mo_art_nml.f90]({}/src/namelists/mo_art_nml.f90)'.format(base_url) }}

(ref_buildrun_nml_assimilation_nml)=
## assimilation_nml

The main switch for the Latent heat nudging scheme is ldass_lhn.


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (assimilation_nml-nlhn_start)=
    **nlhn_start**
  - R(n_dom)
  - -9999.
  - s
  - time in seconds when LHN is applied for the first time
  - ldass_lhn = .TRUE.

* - (assimilation_nml-nlhn_end)=
    **nlhn_end**
  - R(n_dom)
  - -9999.
  - s
  - time in seconds when LHN is applied for the last time
  - ldass_lhn = .TRUE.

* - (assimilation_nml-lhn_coef)=
    **lhn_coef**
  - R(n_dom)
  - 1.0
  -
  - Nudging coefficient of adding the increments
  -

* - (assimilation_nml-fac_lhn_up)=
    **fac_lhn_up**
  - R(n_dom)
  - 2.0
  -
  - Upper limit of the scaling factor of the temperature profile.
  -

* - (assimilation_nml-fac_lhn_down)=
    **fac_lhn_down**
  - R(n_dom)
  - 0.5
  -
  - Lower limit of the scaling factor of the temperature profile.
  -

* - (assimilation_nml-lhn_logscale)=
    **lhn_logscale**
  - L
  - .TRUE.
  -
  - Apply all scaling factors as logarithmic values
  - fac_lhn_down, fac_lhn_up, fac_lhn_artif

* - (assimilation_nml-lhn_updt_rule)=
    **lhn_updt_rule**
  - I(n_dom)
  - 0
  -
  - Rule for humidity update by LHN:
    - 0: LHN updates qv (standard).
    - 1: LHN updates qi if qi>0 and T<0; qv update otherwise.
  -

* - (assimilation_nml-thres_lhn)=
    **thres_lhn**
  - R
  - 0.1/3600.
  - mm/s
  - Minimal value of precipitation rate, either of model or radar. LHN will be applied first for precipitation above it.
  -

* - (assimilation_nml-start_fadeout)=
    **start_fadeout**
  - R
  - 1.0
  -
  - Value to determine, at which model time step a fading out of the increments might start.
  -

* - (assimilation_nml-lhn_qrs)=
    **lhn_qrs**
  - L
  - .TRUE.
  -
  - Use a vertical average of precipitation fluxes as reference to compare with radar observed precipitation, to avoid severe overestimation due to displacement of model surface precipitation.
    If set .FALSE. the model surface precipitation rate is used as reference.
  -

* - (assimilation_nml-rqrsgmax)=
    **rqrsgmax**
  - R(n_dom)
  - 1.0
  -
  - This value determines the height of the vertical averaging, to obtain the reference precipitation rate.
    It is the model layer where the quotion of the maximal precipitation flux occurred for the first time.
  - lhn_qrs = .TRUE.

* - (assimilation_nml-lhn_refbias)=
    **lhn_refbias**
  - L
  - .FALSE.
  -
  - Apply a bias correction between so called reference precipitation (lhn_qrs = .TRUE.) and modelled precipitation at ground. This option is recommended when both quantities shows a systematic bias which cannot be adjusted by changing rqrsgmax.
  - lhn_qrs = .TRUE.

* - (assimilation_nml-ref_bias0)=
    **ref_bias0**
  - R
  - 1.0
  -
  - In case of lhn_refbias = .TRUE. the bias correction starts with this factor. So far, there is no cycling of the factor foreseen, but could be implemented, when it seems to be beneficial.
  - lhn_refbias = .TRUE.

* - (assimilation_nml-dtrefbias)=
    **dtrefbias**
  - R
  - 1800.0
  - s
  - Relaxation time, which defines how fast the bias correction is done.
  - lhn_refbias = .TRUE.

* - (assimilation_nml-lhn_hum_adj)=
    **lhn_hum_adj**
  - L
  - .TRUE.
  -
  - Apply an increment of specific humidity with respect to the estimated temperature increment to maintain the relative humidty
  -

* - (assimilation_nml-lhn_no_ttend)=
    **lhn_no_ttend**
  - L
  - .FALSE.
  -
  - Only apply moisture increments. Temperature increments will only be used for calculation of moisture increments
  - lhn_hum_adj=.TRUE.

* - (assimilation_nml-lhn_incloud)=
    **lhn_incloud**
  - L
  - .TRUE.
  -
  - Apply increments only in model layers where the underlying latent heat release of the model is positive.
  - lhn_artif_only=.FALSE.

* - (assimilation_nml-lhn_limit)=
    **lhn_limit**
  - L
  - .TRUE.
  -
  - Limitation of temperature increments
  - abs_lhn_lim

* - (assimilation_nml-abs_lhn_lim)=
    **abs_lhn_lim**
  - R(n_dom)
  - 50./3600.
  - K/s
  - Lower and upper limit for temperature increments to be added.
  - lhn_limit = .TRUE.

* - (assimilation_nml-lhn_filt)=
    **lhn_filt**
  - L
  - .TRUE.
  -
  - Vertical smoothing of the profile of temperature increments
  -

* - (assimilation_nml-lhn_relax)=
    **lhn_relax**
  - L(n_dom)
  - .FALSE.
  -
  - Horizontal smoothing of radar data but also of incorporated model fields
  - nlhn_relax

* - (assimilation_nml-nlhn_relax)=
    **nlhn_relax**
  - I(n_dom)
  - 2
  - grid points
  - Number of horizontal grid point, where smoothing is applied.
  - lhn_relax = .TRUE.

* - (assimilation_nml-lhn_wweight)=
    **lhn_wweight**
  - L(n_dom)
  - .FALSE.
  -
  - Reduction of the LHN temperature increment in case of strong advection, messured by horizontal wind in 950, 850 and 700 hPa.
    The reduction is done linearly down to zero.
  -

* - (assimilation_nml-lhn_artif)=
    **lhn_artif**
  - L
  - .TRUE.
  -
  - Apply an artificial temperature profile to estimate increments at model grid points without significant precipitation (determined by fac_lhn_artif).
  - fac_lhn_artif, tt_artif_max, zlev_artif_max, std_artif_ma

* - (assimilation_nml-fac_lhn_artif)=
    **fac_lhn_artif**
  - R
  - 5.0
  -
  - Value of the ratio of radar to model precipitation rate, from which an artificial temperature profile is applied
  - lhn_artif=.TRUE.

* - (assimilation_nml-fac_lhn_artif_tune)=
    **fac_lhn_artif_tune**
  - R
  - 1.0
  -
  - Tuning factor to optimize the effectiveness of the artificial profile.
  - lhn_artif=.TRUE.

* - (assimilation_nml-lhn_artif_only)=
    **lhn_artif_only**
  - L
  - .FALSE.
  -
  - Scaling the artificial temperature profile instead of local model profile of latent heat release for calculation the increments at any model grid point.
    The scaling factor is still be determined by the ratio of observed to modelled precipitation rate.
  - tt_artif_max, zlev_artif_max, std_artif_max

* - (assimilation_nml-tt_artif_max)=
    **tt_artif_max**
  - R
  - 0.0015
  - K
  - Maximal temperature of Gaussian shaped function used a artificial temperature profile.
  - lhn_artif, lhn_artif_only

* - (assimilation_nml-zlev_artif_max)=
    **zlev_artif_max**
  - R
  - 1000.0
  - m
  - Height of maximum of Gaussian shaped function used a artificial temperature profile.
  - lhn_artif, lhn_artif_only

* - (assimilation_nml-std_artif_max)=
    **std_artif_max**
  - R
  - 4.0
  - m
  - Parameter defining width of Gaussian shaped function used a artificial temperature profile.
  - lhn_artif, lhn_artif_only

* - (assimilation_nml-nlhnverif_start)=
    **nlhnverif_start**
  - R
  - -9999.
  - s
  - time in seconds when online verification within LHN is active for the first time
  - ldass_lhn = .TRUE.

* - (assimilation_nml-nlhnverif_end)=
    **nlhnverif_end**
  - R
  - -9999.
  - s
  - time in seconds when online verification within LHN is active for the last time
  - ldass_lhn = .TRUE.

* - (assimilation_nml-lhn_diag)=
    **lhn_diag**
  - L
  - .FALSE.
  -
  - Enable a extensive diagnostic output, writing into file lhn.log.
    lhn_diag is set .TRUE. automatically, when online verification is active.
  -

* - (assimilation_nml-lhn_dt_obs)=
    **lhn_dt_obs**
  - R
  - 300.0
  - s
  - Frequency of the radar observations
  -

* - (assimilation_nml-radar_in)=
    **radar_in**
  - C
  - './'
  -
  - Path where the radar data file is expected.
  -

* - (assimilation_nml-radardata_file)=
    **radardata_file**
  - C(n_dom)
  -
  -
  - Name of the radar data file. This might be either in GRIB2 or in NetCDF (recommended).
  -

* - (assimilation_nml-lhn_black)=
    **lhn_black**
  - L(n_dom)
  - .FALSE.
  -
  - Apply a blacklist information in the radar data obtained by comparison against satelite cloud information
  -

* - (assimilation_nml-blacklist_file)=
    **blacklist_file**
  - C(n_dom)
  - 'radarblacklist.nc'
  -
  - Name of blacklist file, containing a mask concerning the quality of the radar data.
    Value 1: good quality
    Value 0: bad quality
    This might be either in GRIB2 or in NetCDF (recommended).
  - lhn_black=.TRUE.

* - (assimilation_nml-lhn_bright)=
    **lhn_bright**
  - L(n_dom)
  - .FALSE.
  -
  - Apply a model intern bright band detection to avoid strong overestimation due to uncertain radar observations.
  -

* - (assimilation_nml-height_file)=
    **height_file**
  - C(n_dom)
  - 'radarheight.nc'
  -
  - Name of file containing the height of the lowest scan for each possible radar station within the given radar composite.
    This file is required, when applying bright band detection.
    This might be either in GRIB2 or in NetCDF (recommended).
  - lhn_bright=.TRUE.

* - (assimilation_nml-nradar)=
    **nradar**
  - I(n_dom)
  - 20
  -
  - Maximal number of radar height layers contained within height_file
  - lhn_bright=.TRUE.

* - (assimilation_nml-lhn_spqual)=
    **lhn_spqual**
  - L(n_dom)
  - .FALSE.
  -
  - Use quality index to infer the horizontal spatial weight of the LHN increments. The quality index is read in as RAD_QUAL variable (besides the RAD_PRECIP variable) from the LHN input file.
  -

* - (assimilation_nml-dace_coupling)=
    **dace_coupling**
  - L
  - .FALSE.
  -
  - Invoke DACE for model equivalents of observations
  - Requires iterate_iau=.T. if init_mode == MODE_IAU (5)

* - (assimilation_nml-dace_time_ctrl)=
    **dace_time_ctrl**
  - I(3)
  - 0
  -
  - Steering parameters for DACE time control: start,end,step
  -

* - (assimilation_nml-dace_debug)=
    **dace_debug**
  - I
  - 0
  -
  - Debugging level for DACE interface
  -

* - (assimilation_nml-dace_output_file)=
    **dace_output_file**
  - C
  - ""
  -
  - Filename for redirection of DACE stdout
  -

* - (assimilation_nml-dace_namelist_file)=
    **dace_namelist_file**
  - C
  - 'namelist'
  -
  - Filename of the file containing the dace namelist
  -
```



Defined and used in: {{ '[src/namelists/mo_assimilation_nml.f90]({}/src/namelists/mo_assimilation_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_ccycle_nml)=
## ccycle_nml

The coupling of the carbon cycle between the atmosphere and land and ocean is configured by the data structure _ccycle_config(jg=1:ndom)%\<param>_, which is a 1-dimensional array extending over all domains.


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ccycle_nml-iccycle)=
    **iccycle**
  - I
  - 0
  -
  - controls the carbon cycle mode:
    - 0: no C-cycle
    - 1: C-cycle with interactive atmospheric {math}`CO_2` concentration
    - 2: C-cycle with prescribed atmospheric {math}`CO_2` concentration
  - dt_vdf > 0.000s and [ljsb](aes_phy_nml-ljsb) = .TRUE. (and atmosphere is coupled to ocean with biogeochemistry)

* - (ccycle_nml-ico2conc)=
    **ico2conc**
  - I
  - 2
  -
  - controls the {math}`CO_2` concentration provided to land/JSBACH and - if coupled to the ocean - to the ocean/HAMOCC
    - 2: constant concentration as defined by  [vmr_co2](ccycle_nml-vmr_co2)
    - 4: transient concentration scenario from file `bc_greenhouse_gases.nc`
  - iccycle = 2

* - (ccycle_nml-vmr_co2)=
    **vmr_co2**
  - R
  - 284.32
  - ppmv
  - constant {math}`CO_2` volume mixing ratio of 1850 (CMIP6)
  - ico2conc = 2
```




(ref_buildrun_nml_comin_nml)=
## comin_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (comin_nml-plugin_list)=
    **plugin_list**
  -
  -
  -
  -
  -
```


Defined and used in: {{ '[src/namelists/mo_comin_nml.f90]({}/src/namelists/mo_comin_nml.f90)'.format(base_url) }}

(ref_buildrun_nml_cloud_mig_nml)=
## cloud_mig_nml

The parameterization of cloud microphysics 'graupel' for the AES physics is configured by a data structure _cloud_mig_config(jg=1:ndom)%\<param>_, which is a 1-dimensional array extending over all  domains.

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (cloud_mig_nml-cia)=
    **cia**
  - R
  - 1.0
  -
  - parameter to change the ice sticking efficiency
  -
```


(ref_buildrun_nml_cloud_two_nml)=
## cloud_two_nml

cloud_two_nml

(ref_buildrun_nml_coupling_mode_nml)=
## coupling_mode_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (coupling_mode_nml-is_coupled_mode)=
    **is_coupled_mode**
  - L
  - ???
  -
  -
  -

* - (coupling_mode_nml-coupled_mode)=
    **coupled_mode**
  - L
  - ???
  -
  -
  -

* - (coupling_mode_nml-coupled_to_ocean)=
    **coupled_to_ocean**
  - L
  - .FALSE.
  -
  - .TRUE.: Required for coupled ocean-atmosphere or ocean-wave similations.
    Indicates the coupling of the model component at hand (e.g. atmo, wave) to the ocean model. Yac coupling routines have to be called.
  -

* - (coupling_mode_nml-coupled_to_waves)=
    **coupled_to_waves**
  - L
  - .FALSE.
  -
  - .TRUE.: Required for coupled wave-atmosphere or wave-ocean similations.
    Indicates the coupling of the model component at hand (e.g. atmo, ocean) to the wave model. Yac coupling routines have to be called.
  -

* - (coupling_mode_nml-coupled_to_atmo)=
    **coupled_to_atmo**
  - L
  - .FALSE.
  -
  - .TRUE.: Required for coupled atmosphere-ocean or atmosphere-wave similations.
    Indicates the coupling of the model at hand (e.g. ocean, wave) to the atmosphere model. Yac coupling routines have to be called.
  -

* - (coupling_mode_nml-coupled_to_aero)=
    **coupled_to_aero**
  - L
  - .FALSE.
  -
  - .TRUE.: Activates the coupling of the atmosphere ([iforcing](run_nml-iforcing)=2,3) to Kinne aerosole input files.
    Kinne aerosol input is taken from python processes rather than direct reading from pre-processed input files. rte-rrtmgp and ecRad radiation supported only. In this case yac coupling routines are called for read in and spatial and temporal interpolation to the chosen ICON grid.
  - For inwp_radiation = 4 (ecRad), irad_aero = 13 and irad_aero = 19 are supported.

* - (coupling_mode_nml-coupled_to_o3)=
    **coupled_to_o3**
  - L
  - .FALSE.
  -
  - .TRUE.: Activates the coupling of the atmosphere ([iforcing](run_nml-iforcing) = 2,3) to o3 input files.
    O3 input is taken from python processes rather than direct reading from pre-processed input files. In this case yac coupling routines are called for read in and spacial and temporal interpolation to the chosen ICON grid.
  - For inwp_radiation = 4 (ecRad), irad_o3 = 5 is supported.

* - (coupling_mode_nml-coupled_to_river)=
    **coupled_to_river**
  - L
  - .FALSE.
  -
  - .TRUE.: Required for coupled atmosphere-river-ocean similations.
    Indicates the coupling of the model at hand (e.g. atmo, ocean) to the river model. Yac coupling routines have to be called.
  -

* - (coupling_mode_nml-use_sens_heat_flux_hack)=
    **use_sens_heat_flux_hack**
  - L
  - .FALSE.
  -
  - .TRUE.: ??
  -

* - (coupling_mode_nml-suppress_sens_heat_flux_hack_over_ice)=
    **suppress_sens_heat_flux_hack_over_ice**
  - L
  - .FALSE.
  -
  - .TRUE.: ??
  -

* - (coupling_mode_nml-coupled_to_output)=
    **coupled_to_output**
  - L
  - .FALSE.
  -
  - enables the output coupling - All suitable variables in the varlist are defined in the coupler for coupling with external output components.
  -
```



Defined and used in: {{ '[src/namelists/mo_coupling_nml.f90]({}/src/namelists/mo_coupling_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_diffusion_nml)=
## diffusion_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (diffusion_nml-lhdiff_temp)=
    **lhdiff_temp**
  - L
  - .TRUE.
  -
  - Diffusion on the temperature field
  -

* - (diffusion_nml-lhdiff_vn)=
    **lhdiff_vn**
  - L
  - .TRUE.
  -
  - Diffusion on the horizontal wind field
  -

* - (diffusion_nml-lhdiff_w)=
    **lhdiff_w**
  - L
  - .TRUE.
  -
  - Diffusion on the vertical wind field
  -

* - (diffusion_nml-lhdiff_q)=
    **lhdiff_q**
  - L
  - .FALSE.
  -
  - Diffusion on QV and QC (water vapor and cloud water)
  -

* - (diffusion_nml-hdiff_order)=
    **hdiff_order**
  - I
  - 5
  -
  - Order of {math}`\nabla` operator for diffusion:
    - -1: no diffusion
    - 2: {math}`\nabla^{2}` diffusion
    - 3: (removed)
    - 4: {math}`\nabla^{4}` diffusion
    - 5: Smagorinsky {math}`\nabla^{2}` diffusion combined with {math}`\nabla^{4}` background diffusion as specified via hdiff_efdt_ratio. Set hdiff_efdt_ratio<=0 for switching off background diffusion.
  -

* - (diffusion_nml-lsmag_3d)=
    **lsmag_3d**
  - L(max_dom)
  - .FALSE.
  -
  - .TRUE.: Use 3D Smagorinsky formulation for computing the horizontal diffusion coefficient (recommended at mesh sizes finer than 1 km if the LES turbulence scheme is not used)
  - hdiff_order=5; itype_vn_diffu=1

* - (diffusion_nml-lhdiff_smag_w)=
    **lhdiff_smag_w**
  - L(max_dom)
  - .FALSE.
  -
  - .TRUE.: Use additional Smagorinsky diffusion for w (recommended at mesh sizes finer than 500 m if the LES turbulence scheme is not used)
  - hdiff_order=5; lhdiff_w=.TRUE.

* - (diffusion_nml-itype_vn_diffu)=
    **itype_vn_diffu**
  - I
  - 1
  -
  - Reconstruction method used for Smagorinsky diffusion:
    - 1: u/v reconstruction at vertices only
    - 2: u/v reconstruction at cells and vertices
  - hdiff_order=5

* - (diffusion_nml-itype_t_diffu)=
    **itype_t_diffu**
  - I
  - 2
  -
  - Discretization of temperature diffusion:
    - 1: {math}`K_h \nabla^2 T`
    - 2: {math}`\nabla \cdot (K_h \nabla T)`
  - hdiff_order=5

* - (diffusion_nml-hdiff_efdt_ratio)=
    **hdiff_efdt_ratio**
  - R
  - 36.0
  -
  - ratio of e-folding time to time step (or 2* time step when using a 3 time level time stepping scheme) (values above 30 are recommended when using hdiff_order=5)
  -

* - (diffusion_nml-hdiff_w_efdt_ratio)=
    **hdiff_w_efdt_ratio**
  - R
  - 15.0
  -
  - ratio of e-folding time to time step for diffusion on vertical wind speed
  -

* - (diffusion_nml-hdiff_min_efdt_ratio)=
    **hdiff_min_efdt_ratio**
  - R
  - 1.0
  -
  - minimum value of hdiff_efdt_ratio near model top
  - hdiff_order=4

* - (diffusion_nml-hdiff_tv_ratio)=
    **hdiff_tv_ratio**
  - R
  - 1.0
  -
  - Ratio of diffusion coefficients for temperature and normal wind: {math}`T:v_{n}`
  -

* - (diffusion_nml-hdiff_multfac)=
    **hdiff_multfac**
  - R
  - 1.0
  -
  - Multiplication factor of normalized diffusion coefficient for nested domains
  - n_dom > 1

* - (diffusion_nml-hdiff_smag_faci)=
    **hdiff_smag_faci**
  - R
  - 0.015
  -
  - Scaling factor for Smagorinsky diffusion at height hdiff_smag_z and below. hdiff_smag_fac {math}`\geq 0`.
  -

* - (diffusion_nml-hdiff_smag_fac2)=
    **hdiff_smag_fac2**
  - R
  - {math}`2\cdot10^{-6}\cdot(1600+25000+\sqrt(1600\cdot(1600+50000))) \approx 0.071`
  -
  - Scaling factor for Smagorinsky diffusion at height hdiff_smag_z2. hdiff_smag_fac2 {math}`\geq 0`. Between hdiff_smag_z and hdiff_smag_z2 the scaling factor changes linearly from hdiff_smag_fac to hdiff_smag_fac2.
  -

* - (diffusion_nml-hdiff_smag_fac3)=
    **hdiff_smag_fac3**
  - R
  - 0.0
  -
  - Scaling factor for Smagorinsky diffusion at height hdiff_smag_z3. hdiff_smag_fac3 {math}`\geq 0`. The three points (hdiff_smag_z2, hdiff_smag_fac2), (hdiff_smag_z3, hdiff_smag_fac3), and (hdiff_smag_z4, hdiff_smag_fac4) determine the quadratic function for the scaling factor between hdiff_smag_z2 and hdiff_smag_z4.
  -

* - (diffusion_nml-hdiff_smag_fac4)=
    **hdiff_smag_fac4**
  - R
  - 1.0
  -
  - Scaling factor for Smagorinsky diffusion at height hdiff_smag_z4 and higher. hdiff_smag_fac4 {math}`\geq 0`.
  -

* - (diffusion_nml-hdiff_smag_z)=
    **hdiff_smag_z**
  - R
  - 32500.0
  - m
  - Height up to which hdiff_smag_fac is used, and where the linear profile up to height hdiff_smag_z2 starts.
  -

* - (diffusion_nml-hdiff_smag_z2)=
    **hdiff_smag_z2**
  - R
  - {math}`1600+50000+\sqrt(1600\cdot(1600+50000)) \approx 60686`
  - m
  - Height with scaling factor hdiff_smag_fac2 where the linear profile starting at hdiff_smag_z ends, and where the quadratic profile up to hdiff_smag_z4 starts. hdiff_smag_z < hdiff_smag_z2 < hdiff_smag_z4.
  -

* - (diffusion_nml-hdiff_smag_z3)=
    **hdiff_smag_z3**
  - R
  - 50000.0
  - m
  - Height with scaling factor hdiff_smag_fac3. Needed to determine the quadratic function between hdiff_smag_z2 and hdiff_smag_z4. hdiff_smag_z3 {math}`\neq` hdiff_smag_z2 {math}`\land` hdiff_smag_z3 {math}`\neq` hdiff_smag_z4.
  -

* - (diffusion_nml-hdiff_smag_z4)=
    **hdiff_smag_z4**
  - R
  - 90000.0
  - m
  - Height from which hdiff_smag_fac4 is used. hdiff_smag_z4 > hdiff_smag_z2.
  -
```



Defined and used in: {{ '[src/namelists/mo_diffusion_nml.f90]({}/src/namelists/mo_diffusion_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_dbg_index_nml)=
## dbg_index_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (dbg_index_nml-dbg_lat_in)=
    **dbg_lat_in**
  -
  -
  -
  -
  -

* - (dbg_index_nml-dbg_lon_in)=
    **dbg_lon_in**
  -
  -
  -
  -
  -

* - (dbg_index_nml-idbg_blk)=
    **idbg_blk**
  -
  -
  -
  -
  -

* - (dbg_index_nml-idbg_elev)=
    **idbg_elev**
  -
  -
  -
  -
  -

* - (dbg_index_nml-idbg_idx)=
    **idbg_idx**
  -
  -
  -
  -
  -

* - (dbg_index_nml-idbg_mxmn)=
    **idbg_mxmn**
  -
  -
  -
  -
  -

* - (dbg_index_nml-idbg_slev)=
    **idbg_slev**
  -
  -
  -
  -
  -

* - (dbg_index_nml-idbg_val)=
    **idbg_val**
  -
  -
  -
  -
  -

* - (dbg_index_nml-str_mod_tst)=
    **str_mod_tst**
  - CHARACTER
  -
  -
  -
  -
```


Defined and used in: {{ '[src/namelists/mo_dbg_nml.f90]({}/src/namelists/mo_dbg_nml.f90)'.format(base_url) }}

(ref_buildrun_nml_dynamics_nml)=
## dynamics_nml
This namelist is relevant if ldynamics=.TRUE.


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (dynamics_nml-divavg_cntrwgt)=
    **divavg_cntrwgt**
  - R
  - 0.5
  -
  - Weight of central cell for divergence averaging
  -

* - (dynamics_nml-lcoriolis)=
    **lcoriolis**
  - L
  - .TRUE.
  -
  - Coriolis force
  -

* - (dynamics_nml-ldeepatmo)=
    **ldeepatmo**
  - L
  - .FALSE.
  -
  - Switch for deep-atmosphere modification.  
    [Note: only the reversible part of the dynamics (largely coincident with what is commonly referred to as "the dynamical core") is modified for the deep atmosphere.  
    Irreversible dynamics of any kind (largely coincident with what is commonly referred to as "the physics") are not explicitly modified.  
    Neither are artificial numerical measures for stabilizing, smoothing and the like modified explicitly.]
  - [iforcing](run_nml-iforcing) = 0, 2, 3  
    is_plane_torus = .FALSE.

* - (dynamics_nml-lmoist_thdyn)=
    **lmoist_thdyn**
  - L
  - .TRUE.
  -
  - Include moisture-dependence of atmospheric heat capacities in thermodynamic equation (automatically reset to .FALSE. in dry model configurations)
  -
```

Defined and used in: {{ '[src/namelists/mo_dynamics_nml.f90]({}/src/namelists/mo_dynamics_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_ensemble_pert_nml)=
## ensemble_pert_nml


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ensemble_pert_nml-use_ensemble_pert)=
    **use_ensemble_pert**
  - L
  - .FALSE.
  -
  - Main switch to activate physics parameter perturbations for ensemble forecasts / ensemble data assimilation; the perturbations are applied via random numbers depending on the perturbationNumber (ensemble member ID) specified in [gribout_nml](ref_buildrun_nml_gribout_nml). Perturbations are always turned off if perturbationNumber {math}`\le` 0
  - [iforcing](run_nml-iforcing) = inwp

* - (ensemble_pert_nml-itype_pert_gen)=
    **itype_pert_gen**
  - I
  - 1
  -
  - Mode of ensemble perturbation generation
    - 1: Equal distribution within perturbation range
    - 2: Discrete distribution with 50% probability for default value and 25% probability for upper and lower extrema
  -

* - (ensemble_pert_nml-timedep_pert)=
    **timedep_pert**
  - I
  - 0
  -
  - Time-dependence of ensemble physics perturbations (except tkred_sfc, which oscillates with a time scale of 20 days)
    - 0: None
    - 1: Random seed for perturbation generation depends on initial date
    - 2: Time-dependent perturbations varying sinusoidally within their range
  - Note: LHN perturbations always follow option 2 if the time dependence is not turned off.

* - (ensemble_pert_nml-fac_rng_spinup)=
    **fac_rng_spinup**
  - I
  - 1
  -
  - Factor for number of spinup calls for random number generator
  -

* - (ensemble_pert_nml-range_gkwake)=
    **range_gkwake**
  - R
  - 1.5
  -
  - Variability range (multiplicative) for low level wake drag constant
  -

* - (ensemble_pert_nml-range_gkdrag)=
    **range_gkdrag**
  - R
  - 0.04
  -
  - Variability range for orographic gravity wave drag constant
  -

* - (ensemble_pert_nml-range_gfrcrit)=
    **range_gfrcrit**
  - R
  - 0.1
  -
  - Variability range for critical Froude number in SSO scheme
  -

* - (ensemble_pert_nml-range_gfluxlaun)=
    **range_gfluxlaun**
  - R
  - 0.75e-3
  -
  - Variability range for non-orographic gravity wave launch momentum flux
  -

* - (ensemble_pert_nml-range_zvz0i)=
    **range_zvz0i**
  - R
  - 0.25
  - m/s
  - Variability range for terminal fall velocity of cloud ice
  - inwp_gscp = 1 or 2

* - (ensemble_pert_nml-range_rain_n0fac)=
    **range_rain_n0fac**
  - R
  - 4.0
  -
  - Multiplicative change of intercept parameter of raindrop size distribution
  - inwp_gscp = 1 or 2

* - (ensemble_pert_nml-range_ccn_Ncn0)=
    **range_ccn_Ncn0**
  - R
  - 1.0
  -
  - For 2-moment mircophysics: multiplicative change of CCN concentration for Segal & Khain cloud activation parameterization. The base value `ccn_Ncn0` may be explicitly specified, otherwise the respective base value from the aerosol scenario `ccn_type` is automatically taken. Not time dependent, regardless of `timedep_pert`.
  - inwp_gscp = 4,5,7

* - (ensemble_pert_nml-range_in_fact)=
    **range_in_fact**
  - R
  - 1.0
  -
  - For 2-moment mircophysics: multiplicative change of IN concentration for ice nucleation parameterization. Not time dependent, regardless of `timedep_pert`.
  - inwp_gscp = 4,5,7

* - (ensemble_pert_nml-range_avel_i)=
    **range_avel_i**
  - R
  - 1.0
  -
  - For 2-moment mircophysics: multiplicative change of cloud ice fall speed. The base value `avel_i` may be explicitly specified, otherwise the default `avel_i` is taken. Not time dependent, regardless of `timedep_pert`.
  - inwp_gscp = 4,5,7

* - (ensemble_pert_nml-range_avel_g)=
    **range_avel_g**
  - R
  - 1.0
  -
  - For 2-moment mircophysics: multiplicative change of graupel fall speed. The base value `avel_g` may be explicitly specified, otherwise the default `avel_g` is taken. Not time dependent, regardless of `timedep_pert`.
  - inwp_gscp = 4,5,7

* - (ensemble_pert_nml-range_cap_ice)=
    **range_cap_ice**
  - R
  - 1.0
  -
  - For 2-moment mircophysics: multiplicative change of capacitance of cloud ice for depositional growth. The base value `cap_ice` may be explicitly specified, otherwise the default `cap_ice` is taken.
  - inwp_gscp = 4,5,7

* - (ensemble_pert_nml-range_cap_snow)=
    **range_cap_snow**
  - R
  - 1.0
  -
  - For 2-moment mircophysics: multiplicative change of capacitance of snow for depositional growth. The base value `cap_snow` may be explicitly specified, otherwise the default `cap_snow` is taken.
  - inwp_gscp = 4,5,7

* - (ensemble_pert_nml-range_entrorg)=
    **range_entrorg**
  - R
  - 0.2e-3
  - 1/m
  - Variability range (additive) for entrainment parameter in convection scheme
  - inwp_convection = 1

* - (ensemble_pert_nml-range_entrorg_mult)=
    **range_entrorg_mult**
  - R
  - 1
  -
  - Asymmetric-multiplicative variation for entrainment parameter in convection scheme, combined with a quadratic reduction of the convective adjustment time scale for positive perturbations. Should be used alternatively to the additive perturbation described above, i.e. setting a factor above 1 shall be combined with range_entrorg = 0.
  - inwp_convection = 1

* - (ensemble_pert_nml-range_rdepths)=
    **range_rdepths**
  - R
  - 5.e3
  - Pa
  - Variability range for maximum allowed shallow convection depth
  - inwp_convection = 1

* - (ensemble_pert_nml-range_rmfdeps)=
    **range_rmfdeps**
  - R
  - 1
  -
  - Multiplicative variation of the rmfdeps parameter, i.e. the fraction of the updraft mass flux that is used as a start value for the downdraft calculation at the level of free sinking
  - inwp_convection = 1

* - (ensemble_pert_nml-range_rprcon)=
    **range_rprcon**
  - R
  - 0.25e-3
  -
  - Variability range for tuning parameter controlling conversion of cloud water into precipitation
  - inwp_convection = 1

* - (ensemble_pert_nml-range_capdcfac_et)=
    **range_capdcfac_et**
  - R
  - 0.75
  -
  - Maximum fraction of CAPE diurnal cycle correction applied in the extratropics
  - icapdcycl = 3

* - (ensemble_pert_nml-range_rhebc)=
    **range_rhebc**
  - R
  - 0.05
  -
  - Variability range for RH threshold for the onset of evaporation below cloud base
  - inwp_convection = 1

* - (ensemble_pert_nml-range_texc)=
    **range_texc**
  - R
  - 0.05
  - K
  - Variability range for temperature excess value in test parcel ascent
  - inwp_convection = 1

* - (ensemble_pert_nml-range_qexc)=
    **range_qexc**
  - R
  - 0.005
  -
  - Variability range for mixing ratio excess value in test parcel ascent
  - inwp_convection = 1

* - (ensemble_pert_nml-range_box_liq)=
    **range_box_liq**
  - R
  - 0.01
  -
  - Variability range for box width scale of liquid clouds in cloud cover scheme
  - inwp_cldcover = 1

* - (ensemble_pert_nml-range_box_liq_asy)=
    **range_box_liq_asy**
  - R
  - 0.25
  -
  - Variability range for asymmetry factor for sub-grid scale liquid cloud distribution
  - inwp_cldcover = 1

* - (ensemble_pert_nml-range_thicklayfac)=
    **range_thicklayfac**
  - R
  - 0.0025
  -
  - Variability range for thick-layer correction factor for sub-grid scale liquid cloud distribution
  - inwp_cldcover = 1

* - (ensemble_pert_nml-range_fac_ccqc)=
    **range_fac_ccqc**
  - R
  - 4
  -
  - Factor for latent-heat correction in CLC-QC relationship in cloud cover scheme
  - inwp_cldcover = 1

* - (ensemble_pert_nml-range_tkhmin)=
    **range_tkhmin**
  - R
  - 0.2
  - m{math}`^2`s{math}`^{-1}`
  - Variability range for minimum vertical diffusion for heat/moisture
  - inwp_turb = 1

* - (ensemble_pert_nml-range_tkmmin)=
    **range_tkmmin**
  - R
  - 0.2
  - m{math}`^2`s{math}`^{-1}`
  - Variability range for minimum vertical diffusion for momentum
  - inwp_turb = 1

* - (ensemble_pert_nml-range_turlen)=
    **range_turlen**
  - R
  - 150
  - m
  - Variability range for turbulent mixing length
  - inwp_turb = 1

* - (ensemble_pert_nml-range_a_hshr)=
    **range_a_hshr**
  - R
  - 1
  -
  - Variability range for scaling factor for extended horizontal shear term
  - inwp_turb = 1

* - (ensemble_pert_nml-range_a_stab)=
    **range_a_stab**
  - R
  - 1
  -
  - Variability range for stability correction
  - inwp_turb = 1

* - (ensemble_pert_nml-range_c_diff)=
    **range_c_diff**
  - R
  - 2.0
  -
  - Range for multiplicative change of length scale factor for vertical diffusion
  - inwp_turb = 1

* - (ensemble_pert_nml-range_q_crit)=
    **range_q_crit**
  - R
  - 1
  -
  - Variability range for critical value for normalized supersaturation in turbulent cloud scheme
  - inwp_turb = 1

* - (ensemble_pert_nml-range_tkred_sfc)=
    **range_tkred_sfc**
  - R
  - 4.0
  -
  - Range for multiplicative change of reduction of minimum diffusion coefficients near the surface
  - inwp_turb = 1

* - (ensemble_pert_nml-range_rlam_heat)=
    **range_rlam_heat**
  - R
  - 8.0
  -
  - Variability range (additive) of laminar transport resistance parameter
  - inwp_turb = 1

* - (ensemble_pert_nml-range_charnock)=
    **range_charnock**
  - R
  - 1.5
  -
  - Variability range (multiplicative!) of upper and lower bound of wind-speed dependent Charnock parameter
  - inwp_turb = 1

* - (ensemble_pert_nml-range_minsnowfrac)=
    **range_minsnowfrac**
  - R
  - 0.1
  -
  - Variability range for minimum value to which snow cover fraction is artificially reduced in case of melting snow
  - idiag_snowfrac = 20

* - (ensemble_pert_nml-range_c_soil)=
    **range_c_soil**
  - R
  - 0.25
  -
  - Variability range for evaporating fraction of soil
  -

* - (ensemble_pert_nml-range_cwimax_ml)=
    **range_cwimax_ml**
  - R
  - 2.0
  -
  - Variability range for capacity of interception storage (multiplicative)
  -

* - (ensemble_pert_nml-range_lhn_coef)=
    **range_lhn_coef**
  - R
  - 0.0
  -
  - Scaling factor for latent heat nudging increments
  - latent heat nudging; i.e. ldass_lhn = .TRUE.

* - (ensemble_pert_nml-range_lhn_artif_fac)=
    **range_lhn_artif_fac**
  - R
  - 0.0
  -
  - Scaling factor for artificial heating profile in latent heat nudging
  - latent heat nudging; i.e. ldass_lhn = .TRUE.

* - (ensemble_pert_nml-range_lhn_down)=
    **range_lhn_down**
  - R
  - 0.0
  -
  - Lower limit for reduction of pre-existing latent heating in LHN
  - latent heat nudging; i.e. ldass_lhn = .TRUE.

* - (ensemble_pert_nml-range_lhn_up)=
    **range_lhn_up**
  - R
  - 0.0
  -
  - Upper limit for increase of pre-existing latent heating in LHN
  - latent heat nudging; i.e. ldass_lhn = .TRUE.

* - (ensemble_pert_nml-range_z0_lcc)=
    **range_z0_lcc**
  - R
  - 0.25
  -
  - Variability range (relative change) of roughness length attributed to each landuse class
  -

* - (ensemble_pert_nml-range_rootdp)=
    **range_rootdp**
  - R
  - 0.2
  -
  - Variability range (relative change) of root depth attributed to each landuse class
  -

* - (ensemble_pert_nml-range_rsmin)=
    **range_rsmin**
  - R
  - 0.2
  -
  - Variability range (relative change) of minimum stomata resistance attributed to each landuse class
  -

* - (ensemble_pert_nml-range_laimax)=
    **range_laimax**
  - R
  - 0.15
  -
  - Variability range (relative change) of leaf area index (maximum of annual cycle) attributed to each landuse class
  -

* - (ensemble_pert_nml-stdev_sst_pert)=
    **stdev_sst_pert**
  - R
  - 0.0
  - K
  - Inserting the standard deviation of SST perturbations (present in the model input data) activates a correction factor for the saturation vapor pressure over oceans, which compensates the systematic increase of evaporation due to the SST perturbations.
  -

* - (ensemble_pert_nml-shift_boxliq_asy)=
    **shift_boxliq_asy**
  - R
  - 0.0
  -
  - Option to shift ensemble mean of tune_box_liq_asy w.r.t. the deterministic value.
  -

* - (ensemble_pert_nml-shift_ratsea)=
    **shift_ratsea**
  - R
  - 0.0
  -
  - Option to shift ensemble mean of rat_sea w.r.t. the deterministic value.
  -
```



Defined and used in: {{ '[src/namelists/mo_ensemble_pert_nml.f90]({}/src/namelists/mo_ensemble_pert_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_gribout_nml)=
## gribout_nml


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (gribout_nml-preset)=
    **preset**
  - C
  - "determ{math}`\dots`"
  -
  - Setting this different to "none" enables a couple of defaults for the other parameters of this namelist. If, additionally, the user tries to set any of these other parameters to a conflicting value, an error message is thrown. Possible values are:
    * "none"
    * "deterministic"(a)
    * "ensemble"(b)
    * "modcomp:deterministic"(c)
    * "modcomp:ensemble"(d)
    They correspond to:
    typeOfGeneratingProcess = 2(a)/4(b)/2(c)/4(d)
    localDefinitionNumber = 254(a)/253(b)/230(c)/230(d)
    typOfProcessedData = 1(a)/5(b)/1(c)/5(d)
    typeOfEnsembleForecast = 192(b)/192(d)
    (Note: "modcomp:{math}`\ldots`" require ecCodes version >= 2.31.0)
  - filetype=2

* - (gribout_nml-tablesVersion)=
    **tablesVersion**
  - I
  - 15
  -
  - Main switch for Table version
  - filetype=2

* - (gribout_nml-backgroundProcess)=
    **backgroundProcess**
  - I
  - 0
  -
  - Background process
      - GRIB2 code table backgroundProcess.table
  - filetype=2

* - (gribout_nml-generatingCenter)=
    **generatingCenter**
  - I
  - -1
  -
  - Output generating center. If this key is not set, center information is taken from the grid file
    * 78: DWD
    * 98: MPIMET{math}`{}^{+}`
    * 98: ECMWF
    (`{}^{+}` The official WMO code for the MIPMET is 252.)
  - filetype=2

* - (gribout_nml-generatingSubcenter)=
    **generatingSubcenter**
  - I
  - -1
  -
  - Output generating Subcenter. If this key is not set, subcenter information is taken from the grid file
    * 255:           DWD
    * 232:           MPIMET
    * 0: ECMWF
  - filetype=2

* - (gribout_nml-generatingProcessIdentifier)=
    **generatingProcessIdentifier**
  - I(n_dom)
  - 1
  -
  - generating Process Identifier
      - GRIB2 code table generatingProcessIdentifier.table
  - filetype=2

* - (gribout_nml-numberOfForecastsInEnsemble)=
    **numberOfForecastsInEnsemble**
  - I
  - -1
  -
  - Local definiton for ensemble products, (only set if value changed from default)
  - filetype=2

* - (gribout_nml-perturbationNumber)=
    **perturbationNumber**
  - I
  - -1
  -
  - Local definiton for ensemble products, (only set if value changed from default)
  - filetype=2

* - (gribout_nml-productionStatusOfProcessedData)=
    **productionStatusOfProcessedData**
  - I
  - 1
  -
  - Production status of data
      - GRIB2 code table 1.3
  - filetype=2

* - (gribout_nml-significanceOfReferenceTime)=
    **significanceOfReferenceTime**
  - I
  - 1
  -
  - Significance of reference time
      - GRIB2 code table 1.2
  - filetype=2

* - (gribout_nml-typeOfEnsembleForecast)=
    **typeOfEnsembleForecast**
  - I
  - -1
  -
  - Local definiton for ensemble products (only set if value changed from default)
  - filetype=2

* - (gribout_nml-typeOfGeneratingProcess)=
    **typeOfGeneratingProcess**
  - I
  - -1
  -
  - Type of generating process
      - GRIB2 code table 4.3
  - filetype=2

* - (gribout_nml-typeOfProcessedData)=
    **typeOfProcessedData**
  - I
  - -1
  -
  - Type of data
      - GRIB2 code table 1.4
  - filetype=2

* - (gribout_nml-localDefinitionNumber)=
    **localDefinitionNumber**
  - I
  - -1
  -
  - local Definition Number:
    * 254: Deterministic system
    * 253: Ensemble system:
    * 230: Model composition
      - GRIB2 code table grib2LocalSectionNumber.78.table.
    Note that in case of 230 ("Model composition") `preset="modcomp:deterministic/ensemble"` has to be used to choose between deterministic or ensemble systems!
    In addition, 230 requires ecCodes version >= 2.31.0.
  - filetype=2
    generatingCenter=78/80/215

* - (gribout_nml-localNumberOfExperiment)=
    **localNumberOfExperiment**
  - I
  - 1
  -
  - local Number of Experiment
  - filetype=2
    generatingCenter=78/80/215

* - (gribout_nml-localTypeOfEnsembleForecast)=
    **localTypeOfEnsembleForecast**
  - I
  - -1
  -
  - Local definiton for ensemble products (only set if value changed from default)
  - filetype=2
    generatingCenter=78/80/215

* - (gribout_nml-typeOfGrib2TileTemplate)=
    **typeOfGrib2TileTemplate**
  - C
  - "DWD"
  -
  - Type of GRIB2 templates which are used for decoding tiled surface fields
    * "WMO": official WMO templates (55, 59)
    * "DWD": local DWD templates (40455, 40456)
  - filetype = 2

* - (gribout_nml-lspecialdate_invar)=
    **lspecialdate_invar**
  - L
  - .FALSE.
  -
  - Special reference date for invariant and climatological fields
    * .TRUE.: set special reference date 0001-01-01, 00:00
    * .FASLE.: no special reference date
  - filetype = 2

* - (gribout_nml-ldate_grib_act)=
    **ldate_grib_act**
  - L
  - .TRUE.
  -
  - GRIB creation date
    * .TRUE.: add creation date
    * .FALSE.: add dummy date
  - filetype=2

* - (gribout_nml-lgribout_24bit)=
    **lgribout_24bit**
  - L
  - .FALSE.
  -
  - If TRUE, write thermodynamic fields {math}`\rho`, {math}`\theta_{v}`, {math}`T`, {math}`p` with 24bit precision instead of 16bit
  - filetype=2

* - (gribout_nml-lgribout_compress_ccsds)
    **lgribout_compress_ccsds**
  - L
  - .FALSE.
  -
  - Enable CCSDS second level compression
  - filetype=2


* - (gribout_nml-grib_lib_compat)=
    **grib_lib_compat**
  - C
  - "current"
  -
  - Type of GRIB library backward compatibility adjustment:
    * "current": No adjustment
    * "eccodes:2.31.0": Please see Section below.
  - filetype=2
    ecCodes version >= 2.32.0

* - (gribout_nml-model_components)=
    **model_components**
  - C(4)
  - " ", " ", " "
  -
  - Model components. Currently, the following options:
    * "icon-nwp" (1)
    * "art-nwp"
    * "ocean-nwp"
    * "waves-nwp"
    can be combined (without repetition) to represent a model configuration with max. 4 components.  
    For example, use:  
    - "ocean-nwp" for a stand-alone ocean-model application.
    - "icon-nwp", ''art-nwp'' for an ICON-ART application (lart = .TRUE.).  
    Note that it is the responsibility of the user, to choose a setting that is in line with the actual model composition!  
    (1) This stands for the model component: atmosphere + land & water surface (as these are no indpendent components in the NWP configuration). This is the original basic component for NWP applications and therefore habitually associated with the term "icon-nwp" by users. The other model components were introduced later.)
  - filetype=2
    generatingCenter=78/80/215
    localDefinitionNumber=230
    ecCodes version >= 2.31.0

* - (gribout_nml-local-Production-Context)=
    **local-Production-Context**
  - I(n_dom)
  - -1
  -
  - Local production context:  
    - GRIB2 code table 2.233.table (local DWD definitions)
    This key encodes supplementary information to the generating-Process-Identifier (if required).
  - filetype=2  
    generatingCenter=78/80/215  
    localDefinitionNumber=230  
```




### Notes on the GRIB library backward compatiblity adjustment:


The I/O library CDI uses the ecCodes library for GRIB file handling.
Sometimes, version updates of ecCodes come along with a change in behavior
that results in some GRIB metadata having different values in GRIB output files (which cannot be avoided without additional measures).
This is undesirable, at least in operational NWP.
In order to allow for maintaining continuity, we try to "overwrite" such new behavior with an explicit reproduction of the prior behavior, if possible.
The following Table describes the options of the associated namelist parameter `grib_lib_compat`
in more detail.


```{list-table}
:header-rows: 1

* - (Notes on the GRIB library backward compatiblity adjustment:-grib_lib_compat)=
    **grib_lib_compat**
  - **Description**
  - **Applies to GRIB lib version**
  - **Notes**

* - "current"
  - No adjustment applied.
  -
  -

* - "eccodes:2.31.0"
  - CDI uses the ecCodes sample GRIB file "GRIB2.tmpl" as a starting file for ecCodes. The `SecondFixedSurface` GRIB keys are assigned the following values in GRIB2.tmpl:
    * `typeOfSecondFixedSurface` = 255 (MISSING)
    * `scaleFactorOfSecondFixedSurface` = 255 (MISSING)
    * `scaledValueOfSecondFixedSurface` = 2,147,483,647 (MISSING)
    Now, if `typeOfSecondFixedSurface` is set to a value >= 10, 102 say ("Specific altitude above mean sea level"), but `scaleFactorOfSecondFixedSurface` and `scaledValueOfSecondFixedSurface` are not explicitly set, the result is as follows:
    **(1) For ecCodes version < 2.32.0**
    * `typeOfSecondFixedSurface`
    = 102
    * `scaleFactorOfSecondFixedSurface` = 0
    * `scaledValueOfSecondFixedSurface` = 0
    **(2) For ecCodes version >= 2.32.0**
    * `typeOfSecondFixedSurface` = 102
    * `scaleFactorOfSecondFixedSurface` = MISSING
    * `scaledValueOfSecondFixedSurface` = MISSING
    With `grib_lib_compat` = "eccodes:2.31.0", we try to reproduce behavior (1) even for ecCodes version >= 2.32.0.
  - ecCodes version >= 2.32.0
  - If the used ecCodes version is
    < 2.32.0, "eccodes:2.31.0" will be overwritten with "current".
```



Defined and used in: {{ '[src/namelists/mo_gribout_nml.f90]({}/src/namelists/mo_gribout_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_grid_nml)=
## grid_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (grid_nml-lplane)=
    **lplane**
  - L
  - .FALSE.
  -
  - planar option
  -

* - (grid_nml-is_plane_torus)=
    **is_plane_torus**
  - L
  - .FALSE.
  -
  - f-plane approximation on triangular grid
  -

* - (grid_nml-corio_lat)=
    **corio_lat**
  - R
  - 0.0
  - deg
  - Center of the f-plane is located at this geographical latitude
  - lplane=.TRUE. and is_plane_torus=.TRUE.

* - (grid_nml-grid_angular _velocity)=
    **grid_angular _velocity**
  - R
  - Earth's
  - rad/s
  - The angular velocity in rad per sec.
  -

* - (grid_nml-l_scm_mode)=
    **l_scm_mode**
  - L
  - .FALSE.
  -
  - Single Column Model (SCM) mode. Can be extended to equivalent LES and CRM setups by setting ldynamics=.TRUE. .
  - is_plane_torus=.TRUE.

* - (grid_nml-l_limited_area)=
    **l_limited_area**
  - L
  - .FALSE.
  -
  -
  -

* - (grid_nml-grid_rescale_factor)=
    **grid_rescale_factor**
  - R
  - 1.0
  -
  - Defined as the inverse of the reduced-size earth reduction factor {math}`X`. Choose `grid_rescale_factor` < 1 for a reduced-size earth.
  -

* - (grid_nml-lrescale_timestep)=
    **lrescale_timestep**
  - L
  - .FALSE.
  -
  - if .TRUE. then the timestep will be multiplied by `grid_rescale_factor`.
  -

* - (grid_nml-lrescale_ang_vel)=
    **lrescale_ang_vel**
  - L
  - .FALSE.
  -
  - if .TRUE. then the angular velocity will be divided by `grid_rescale_factor`.
  -

* - (grid_nml-lfeedback)=
    **lfeedback**
  - L(n_dom)
  - .TRUE.
  -
  - Specifies if feedback to parent grid is performed. Setting lfeedback(1)=.FALSE. turns off feedback for all nested domains; to turn off feedback for selected nested domains, set lfeedback(1)=.TRUE. and set `.FALSE.` for the desired model domains
  - n_dom > 1

* - (grid_nml-ifeedback_type)=
    **ifeedback_type**
  - I
  - 2
  -
  - 1: incremental feedback
    - 2: relaxation-based feedback
    Note: vertical nesting requires option 2 to run numerically stable over longer time periods
  - n_dom > 1

* - (grid_nml-start_time)=
    **start_time**
  - R(n_dom)
  - 0.0
  - s
  - Time when a nested domain starts to be active. Relative time w.r.t.
    experiment start date (`ini_datetime_string` / `experimentStartDate`).
    (namelist entry is ignored for the global domain)
  - n_dom > 1

* - (grid_nml-end_time)=
    **end_time**
  - R(n_dom)
  - 1.E30
  - s
  - Time when a nested domain terminates. Relative time w.r.t.
    experiment start date (`ini_datetime_string` / `experimentStartDate`).
    (namelist entry is ignored for the global domain)
  - n_dom > 1

* - (grid_nml-patch_weight)=
    **patch_weight**
  - R(n_dom)
  - 0.0
  -
  - If patch_weight is set to a value > 0 for any of the first level child patches, processor splitting will be performed, i.e. every of the first level child patches gets a subset of the total number or processors corresponding to its patch_weight. A value of 0. corresponds to exactly 1 processor for this patch, regardless of the total number of processors. For the root patch and higher level childs, patch_weight is not used. However, patch_weight must be set to 0 for these patches to avoid confusion.
  - n_dom > 1

* - (grid_nml-lredgrid_phys)=
    **lredgrid_phys**
  - L(n_dom)
  - .FALSE.
  -
  - If set to .TRUE. radiation is calculated on a reduced grid (= one grid level higher).
    Needs to be set for each model domain separately; for the global domain, the file containing the reduced grid must be specified in the variable `radiation_grid_filename`
  -

* - (grid_nml-nexlevs_rrg_vnest)=
    **nexlevs_rrg_vnest**
  - I
  - 8
  -
  - Maximum number of extra (additional) model layers used for calculating radiation if vertical nesting is combined with a reduced radiation grid. For these layers, temperature and pressure are copied from the parent domain (thus, the difference in the number of model levels constitutes another upper limit). Higher values improve the consistency of radiative flux divergences near the top of a vertically nested domain.
    lredgrid_phys = .TRUE., lvert_nest = .TRUE., latm_above_top = .TRUE.
  -

* - (grid_nml-dynamics_grid_filename)=
    **dynamics_grid_filename**
  - C
  -
  -
  - Array of the grid filenames to be used by the dycore. May contain the keyword `<path>` which will be substituted by `model_base_dir`.
  -

* - (grid_nml-dynamics_parent_ grid_id)=
    **dynamics_parent_ grid_id**
  - I(n_dom)
  - {math}`i-1`
  -
  - Array of the indexes of the parent grid filenames, as described by the dynamics_grid_filename array. Indexes start at 1, an index of 0 indicates no parent.  Specification of this namelist parameter is only required if more than one domain is in use _and_ the grid files are rather old s.t. they do not contain a `uuidOfParHGrid` global attribute.
  -

* - (grid_nml-radiation_grid_filename)=
    **radiation_grid_filename**
  - C
  -
  -
  - Grid filename to be used for the radiation model on the coarsest grid. Filled only if the radiation grid is different from the dycore grid. May contain the keyword `<path>` which will be substituted by `model_base_dir`.
  - lredgrid_phys=.TRUE.

* - (grid_nml-create_vgrid)=
    **create_vgrid**
  - L
  - .FALSE.
  -
  - .TRUE.: Write vertical grid files containing (`vct_a`, `vct_b`, `z_ifc`, and `z_ifv`.
  -

* - (grid_nml-vertical_grid_filename)=
    **vertical_grid_filename**
  - C(n_dom)
  -
  -
  - Array of filenames. These files contain the vertical grid definition (`vct_a`, `vct_b`, `z_ifc`). If empty, the vertical grid is created within ICON during the setup phase.
  -

* - (grid_nml-vct_filename)=
    **vct_filename**
  - C
  -
  -
  - Filename of ASCII file containing the 1D vertical coordinate tables `vct_a`, `vct_b`. See [information on vertical level distribution](ref_buildrun_nml_vlev_distr) for further information on the format. If empty, `vct_a`, `vct_b` are created within ICON during the setup phase.
  -

* - (grid_nml-use_duplicated_connectivity)=
    **use_duplicated_connectivity**
  - L
  - .TRUE.
  -
  - if .TRUE., the zero connectivity is replaced by the last non-zero value
  -

* - (grid_nml-use_dummy_cell_closure)=
    **use_dummy_cell_closure**
  - L
  - .FALSE.
  -
  - if .TRUE. then create a dummy cell and connect it to cells and edges with no neighbor
  -
```



Defined and used in: {{ '[src/namelists/mo_grid_nml.f90]({}/src/namelists/mo_grid_nml.f90)'.format(base_url) }}




(ref_buildrun_nml_gridref_nml)=
## gridref_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (gridref_nml-grf_intmethod_c)=
    **grf_intmethod_c**
  - I
  - 2
  -
  - Interpolation method for grid refinement (cell-based dynamical variables):
    - 1: parent-to-child copying
    - 2: gradient-based interpolation
  - n_dom > 1

* - (gridref_nml-grf_intmethod_ct)=
    **grf_intmethod_ct**
  - I
  - 2
  -
  - Interpolation method for grid refinement (cell-based tracer variables):
    - 1: parent-to-child copying
    - 2: gradient-based interpolation
  - n_dom > 1

* - (gridref_nml-grf_intmethod_e)=
    **grf_intmethod_e**
  - I
  - 6
  -
  - Interpolation method for grid refinement (edge-based variables):
    - 1: (removed)
    - 2: RBF interpolation
    - 3: (removed)
    - 4: combination gradient-based / RBF
    - 5: (removed)
    - 6: same as 4, but direct interpolation of mass fluxes along nest interface edges
  - n_dom > 1

* - (gridref_nml-grf_velfbk)=
    **grf_velfbk**
  - I
  - 1
  -
  - Method of velocity feedback:
    - 1: average of child edges 1 and 2
    - 2: 2nd-order method using RBF interpolation
  - n_dom > 1

* - (gridref_nml-grf_scalfbk)=
    **grf_scalfbk**
  - I
  - 2
  -
  - Feedback method for dynamical scalar variables ({math}`T, p_{sfc}`):
    - 1: area-weighted averaging
    - 2: bilinear interpolation
  - n_dom > 1

* - (gridref_nml-grf_tracfbk)=
    **grf_tracfbk**
  - I
  - 2
  -
  - Feedback method for tracer variables:
    - 1: area-weighted averaging
    - 2: bilinear interpolation
  - n_dom > 1

* - (gridref_nml-rbf_vec_kern_grf_e)=
    **rbf_vec_kern_grf_e**
  - I
  - 1
  -
  - RBF kernel for grid refinement (edges):
    - 1: Gaussian
    - 2: {math}`1/(1+r^{2})`
    - 3: inverse multiquadric
  - n_dom > 1

* - (gridref_nml-rbf_scale_grf_e)=
    **rbf_scale_grf_e**
  - R(n_dom)
  - 0.5
  -
  - RBF scale factor for grid refinement (lateral boundary interpolation to edges). Refers to the respective parent domain and thus does not need to be specified for the innermost nest. Lower values than the default of 0.5 are needed for child mesh sizes less than about 500 m.
  - n_dom > 1

* - (gridref_nml-denom_diffu_t)=
    **denom_diffu_t**
  - R
  - 135
  -
  - Denominator for lateral boundary diffusion of temperature
  - n_dom > 1

* - (gridref_nml-denom_diffu_v)=
    **denom_diffu_v**
  - R
  - 200
  -
  - Denominator for lateral boundary diffusion of velocity
  - n_dom > 1

* - (gridref_nml-l_density_nudging)=
    **l_density_nudging**
  - L
  - .FALSE.
  -
  - .TRUE.: Apply density nudging near lateral nest boundary if grf_intmethod_e {math}`\in`
    {2,4
    }
  - n_dom > 1 .AND. lfeedback = .TRUE.

* - (gridref_nml-fbk_relax_timescale)=
    **fbk_relax_timescale**
  - R
  - 10800
  -
  - Relaxation time scale for feedback
  - n_dom>1 .AND. lfeedback = .TRUE. .AND. ifeedback_type = 2
```



Defined and used in: {{ '[src/namelists/mo_gridref_nml.f90]({}/src/namelists/mo_gridref_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_hamocc_nml)=
## hamocc_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (hamocc_nml-atm_co2)=
    **atm_co2**
  -
  - 278._wp
  -
  -
  -

* - (hamocc_nml-atm_n2)=
    **atm_n2**
  -
  - 802000._wp
  -
  -
  -

* - (hamocc_nml-atm_o2)=
    **atm_o2**
  -
  - 196800._wp
  -
  -
  -

* - (hamocc_nml-bkcya_Fe)=
    **bkcya_Fe**
  -
  - 30.e-8_wp
  -
  -
  -

* - (hamocc_nml-bkcya_P)=
    **bkcya_P**
  -
  - 5.e-8_wp
  -
  -
  -

* - (hamocc_nml-calmax)=
    **calmax**
  -
  - 0.15_wp
  -
  -
  -

* - (hamocc_nml-cya_growth_max)=
    **cya_growth_max**
  - REAL
  - 0.2_wp
  -
  -
  -

* - (hamocc_nml-cycdec)=
    **cycdec**
  -
  - 0.1_wp
  -
  -
  -

* - (hamocc_nml-deltacalc)=
    **deltacalc**
  -
  - 0._wp
  -
  -
  -

* - (hamocc_nml-deltaorg)=
    **deltaorg**
  -
  - 0._wp
  -
  -
  -

* - (hamocc_nml-deltasil)=
    **deltasil**
  -
  - 0._wp
  -
  -
  -

* - (hamocc_nml-denitrification)=
    **denitrification**
  -
  - 1.82e-3_wp
  -
  -
  -

* - (hamocc_nml-disso_op_q10)=
    **disso_op_q10**
  -
  - 2.3_wp
  -
  -
  -

* - (hamocc_nml-disso_op_tref)=
    **disso_op_tref**
  -
  - 20._wp
  -
  -
  -

* - (hamocc_nml-doc_remin_q10)=
    **doc_remin_q10**
  -
  - 2._wp
  -
  -
  -

* - (hamocc_nml-doc_remin_tref)=
    **doc_remin_tref**
  -
  - 10._wp
  -
  -
  -

* - (hamocc_nml-dremcalc)=
    **dremcalc**
  -
  - 0.075_wp
  -
  -
  -

* - (hamocc_nml-dremopal)=
    **dremopal**
  -
  - 0.01_wp
  -
  -
  -

* - (hamocc_nml-drempoc)=
    **drempoc**
  -
  - 0.026_wp
  -
  -
  -

* - (hamocc_nml-dzsed)=
    **dzsed**
  -
  -
  -
  -
  -

* - (hamocc_nml-grazra)=
    **grazra**
  - REAL
  - 1.0_wp
  -
  -
  -

* - (hamocc_nml-hion_solver)=
    **hion_solver**
  -
  - 1
  -
  -
  -

* - (hamocc_nml-i_settling)=
    **i_settling**
  -
  - 1
  -
  -
  -

* - (hamocc_nml-inpw)=
    **inpw**
  -
  -
  -
  -
  -

* - (hamocc_nml-ks)=
    **ks**
  - INTEGER
  - 0.31_wp
  -
  -
  -

* - (hamocc_nml-l_N_cycle)=
    **l_N_cycle**
  -
  -
  -
  -
  -

* - (hamocc_nml-l_PDM_settling)=
    **l_PDM_settling**
  -
  -
  -
  -
  -

* - (hamocc_nml-l_bgc_check)=
    **l_bgc_check**
  -
  -
  -
  -
  -

* - (hamocc_nml-l_cpl_co2)=
    **l_cpl_co2**
  -
  -
  -
  -
  -

* - (hamocc_nml-l_cyadyn)=
    **l_cyadyn**
  -
  - .TRUE.
  -
  -
  -

* - (hamocc_nml-l_doc_q10)=
    **l_doc_q10**
  -
  -
  -
  -
  -

* - (hamocc_nml-l_dynamic_pi)=
    **l_dynamic_pi**
  -
  -
  -
  -
  -

* - (hamocc_nml-l_hamocc_vertint)=
    **l_hamocc_vertint**
  -
  -
  -
  -
  -

* - (hamocc_nml-l_implsed)=
    **l_implsed**
  -
  -
  -
  -
  -

* - (hamocc_nml-l_init_bgc)=
    **l_init_bgc**
  -
  -
  -
  -
  -

* - (hamocc_nml-l_limit_sal)=
    **l_limit_sal**
  -
  -
  -
  -
  -

* - (hamocc_nml-l_opal_q10)=
    **l_opal_q10**
  -
  - .TRUE.
  -
  -
  -

* - (hamocc_nml-l_opalsed_q10)=
    **l_opalsed_q10**
  -
  -
  -
  -
  -

* - (hamocc_nml-l_poc_q10)=
    **l_poc_q10**
  -
  - .TRUE.
  -
  -
  -

* - (hamocc_nml-l_re)=
    **l_re**
  -
  - .TRUE.
  -
  -
  -

* - (hamocc_nml-l_up_sedshi)=
    **l_up_sedshi**
  -
  -
  -
  -
  -

* - (hamocc_nml-l_virtual_tep)=
    **l_virtual_tep**
  -
  - .TRUE.
  -
  -
  -

* - (hamocc_nml-mc_depth)=
    **mc_depth**
  - REAL
  - 100._wp
  -
  -
  -

* - (hamocc_nml-mc_fac)=
    **mc_fac**
  - REAL
  - 2.0_wp
  -
  -
  -

* - (hamocc_nml-no3nh4red)=
    **no3nh4red**
  -
  - 0.002_wp
  -
  -
  -

* - (hamocc_nml-no3no2red)=
    **no3no2red**
  -
  - 0.002_wp
  -
  -
  -

* - (hamocc_nml-opal_remin_q10)=
    **opal_remin_q10**
  -
  - 2.6_wp
  -
  -
  -

* - (hamocc_nml-opal_remin_tref)=
    **opal_remin_tref**
  -
  - 10._wp
  -
  -
  -

* - (hamocc_nml-poc_remin_q10)=
    **poc_remin_q10**
  -
  - 2.1_wp
  -
  -
  -

* - (hamocc_nml-poc_remin_tref)=
    **poc_remin_tref**
  -
  - 10._wp
  -
  -
  -

* - (hamocc_nml-sinkspeed_calc)=
    **sinkspeed_calc**
  -
  - 30._wp
  -
  -
  -

* - (hamocc_nml-sinkspeed_martin_ez)=
    **sinkspeed_martin_ez**
  -
  - 3.5_wp
  -
  -
  -

* - (hamocc_nml-sinkspeed_opal)=
    **sinkspeed_opal**
  -
  - 30._wp
  -
  -
  -

* - (hamocc_nml-sinkspeed_poc)=
    **sinkspeed_poc**
  -
  - 5._wp
  -
  -
  -
```


Defined and used in: {{ '[src/hamocc/icon_specific/mo_hamocc_nml.f90]({}/src/hamocc/icon_specific/mo_hamocc_nml.f90)'.format(base_url) }}

(ref_buildrun_nml_initicon_nml)=
## initicon_nml


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (initicon_nml-init_mode)=
    **init_mode**
  - I
  - 2
  -
  - 1: MODE_DWDANA
    start from DWD analysis or FG
    - 2: MODE_IFSANA
    start from IFS analysis
    - 3: MODE_COMBINED
    IFS atm + ICON/GME soil
    - 4: MODE_COSMO
    start from prognostic set of variables as used by COSMO
    - 5: MODE_IAU
    start from DWD analysis with incremental analysis update.
    - 7: MODE_ICONVREMAP
    start from hor. interpolated DWD initialized analysis data with subsequent vertical remapping
  -

* - (initicon_nml-dt_ana)=
    **dt_ana**
  - R
  - 10800
  - s
  - Time interval of assimilation cycle.
  - icpl_da_sfcevap>= 2

* - (initicon_nml-dt_iau)=
    **dt_iau**
  - R
  - 10800
  - s
  - Duration of incremental analysis update (IAU) procedure. Start time for IAU is the actual model start time (see below).
  - init_mode=5

* - (initicon_nml-dt_shift)=
    **dt_shift**
  - R
  - 0
  - s
  - Time by which the actual model start time is shifted ahead of the nominal date. The latter is given by either `ini_datetime_string` or `experimentStartDate`. dt_shift must be NEGATIVE, usually {math}`- 0.5` dt_iau.
  - init_mode=5

* - (initicon_nml-iterate_iau)=
    **iterate_iau**
  - L
  - .FALSE.
  -
  - If .TRUE., the IAU phase is calculated twice with halved dt_iau in first cycle. This allows writing a fully initialized analysis at the nominal initialization date while using a centered IAU window for the forecast.
  - init_mode=5
    and [dt_shift](initicon_nml-dt_shift) < 0

* - (initicon_nml-rho_incr_filter_wgt)=
    **rho_incr_filter_wgt**
  - R
  - 0
  -
  - Vertical filtering weight on density increments
  - init_mode=5

* - (initicon_nml-niter_diffu)=
    **niter_diffu**
  - I
  - 10
  -
  - Number of diffusion iterations applied on wind increments
  - init_mode=5

* - (initicon_nml-niter_divdamp)=
    **niter_divdamp**
  - I
  - 25
  -
  - Number of divergence damping iterations applied on wind increments
  - init_mode=5

* - (initicon_nml-type_iau_wgt)=
    **type_iau_wgt**
  - I
  - 1
  -
  - Weighting function for performing IAU
    - 1: Top-Hat
    - 2: SIN2
  - init_mode=5

* - (initicon_nml-nlevsoil_in)=
    **nlevsoil_in**
  - I
  - 4
  -
  - number of soil levels of input data
  - init_mode=2

* - (initicon_nml-zpbl1)=
    **zpbl1**
  - R
  - 500.0
  - m
  - bottom height (AGL) of layer used for gradient computation
  -

* - (initicon_nml-zpbl2)=
    **zpbl2**
  - R
  - 1000.0
  - m
  - top height (AGL) of layer used for gradient computation
  -

* - (initicon_nml-lread_ana)=
    **lread_ana**
  - L
  - .TRUE.
  -
  - If .FALSE., ICON is started from first guess only. Analysis field is not required, and skipped if provided.
  - init_mode=1,3

* - (initicon_nml-use_lakeiceana)=
    **use_lakeiceana**
  - L
  - .FALSE.
  -
  - If .TRUE., analysis data for sea ice fraction are also used for freshwater lakes (for the time being restricted to the Great Lakes; extension to other lakes needs to be tested)
  - init_mode=5

* - (initicon_nml-qcana_mode)=
    **qcana_mode**
  - I
  - 0
  -
  - If > 0, analysis increments for cloud water concentration are read and processed.
    - 1: QC increments are added to QV increments
    - 2: QC increments are added to QC if clouds are present, otherwise to QV increments
  - init_mode=5

* - (initicon_nml-qiana_mode)=
    **qiana_mode**
  - I
  - 0
  -
  - 1: analysis increments for cloud ice concentration are read and processed.
  - init_mode=5

* - (initicon_nml-qrsgana_mode)=
    **qrsgana_mode**
  - I
  - 0
  -
  - 1: analysis increments for rain, snow and graupel mass concentrations are read and processed. In case of the 2-moment microphysics (inwp_gscp=4,5,6), also hail mass concentration increments are processed.
  - init_mode=5

* - (initicon_nml-qnxana_2mom_mode)=
    **qnxana_2mom_mode**
  - I
  - 0
  -
  - Only effective in case of 2-moment microphysics (inwp_gscp=3,4,5,6). Affects the analysis increments of the the number concentrations of those hydrometeors in IAU which have been selected by the settings of qcana_mode, qiana_mode and qrsgana_mode:
    - 0: analysis increments are not taken from analysis files but diagnosed based on the mass concentrations (from fg) and mass increments.
    - 1: analysis increments are taken from the analysis files. If missing for a specific hydrometeor type, they are diagnosed similar to option 0 as a fallback.
  - init_mode=5, inwp_gscp=3,4,5,6

* - (initicon_nml-icpl_da_sfcevap)=
    **icpl_da_sfcevap**
  - I
  - 0
  -
  - Coupling between data assimilation and model parameters controlling surface evaporation (bare soil and plants). Choosing values > 0 requires itype_vegetation_cycle=2 :
    - 0: off
    - 1: use time-filtered T2M bias provided by the soil moisture analysis
    - 2: use in addition a time-filtered RH increment at the lowest model level (requires assimilation of RH2M)
    - 3: as option 2, but use a time-filtered temperature increment at the lowest model level instead of the T2M bias provided by the SMA (requires assimilation of T2M and RH2M)
    - 4: as option 3, but uses the minimum evaporation resistance (default set by cr_bsmin) instead of c_soil for adaptive tuning of bare-soil evaporation
    - 5: as option 4, but additional adjustment of hydraulic diffusivity (capillary transport) and asymmetry factor for stomata resistance
    - 6: as option 5, but additional time-dependent adjustment of rlam_heat and bare-soil evaporation resistance depending on a combination of meteorological conditions and time-filtered data assimilation increments
  - init_mode=5

* - (initicon_nml-smi_relax_timescale)=
    **smi_relax_timescale**
  - R
  - 20.0
  - days
  - Relaxation time scale for ICON-internal soil moisture adjustment, referring to a filtered RH increment of 1%. Setting the time scale to zero turns off the soil moisture adjustment.
  - icpl_da_sfcevap {math}`\ge` 2

* - (initicon_nml-itype_sma)=
    **itype_sma**
  - I
  - 1
  -
  - Type of soil moisture analysis used
    - 1: use external soil moisture analysis provided by the data assimilation
    - 2: use ICON-internal SMA based on adaptive parameter tuning
  - init_mode=5; icpl_da_sfcevap {math}`\ge` 3

* - (initicon_nml-icpl_da_snowalb)=
    **icpl_da_snowalb**
  - I
  - 0
  -
  - Coupling between temperature bias inferred from data assimilation and snow albedo
    - 0: off
    - 1: on; requires assimilation of T2M and cycling of a time-filtered temperature increment at the lowest model level
    - 2: as option 1, but additional adaptation of sea-ice albedo
    - 3: as option 2, but additional adaptation of snow-cover fraction diagnosis
  - init_mode=5; icpl_da_sfcevap {math}`\ge` 3

* - (initicon_nml-icpl_da_landalb)=
    **icpl_da_landalb**
  - I
  - 0
  -
  - Coupling between temperature/humidity bias inferred from data assimilation and albedo of snow-free land
    - 0: off
    - 1: on; requires assimilation of T2M and RH2M and cycling of the full set of filtered assimilation increments coming along with icpl_da_sfcevap {math}`\ge` 5
  - init_mode=5; icpl_da_sfcevap {math}`\ge` 5

* - (initicon_nml-icpl_da_seaice)=
    **icpl_da_seaice**
  - I
  - 0
  -
  - Coupling between temperature bias inferred from data assimilation and seaice scheme
    - 0: off
    - 1: add filtered T increment to initial seaice temperature
    - 2: as above, and additional adaptive tuning of bottom heat flux if lbottom_hflux = true
  - init_mode=5; icpl_da_sfcevap {math}`\ge` 3

* - (initicon_nml-icpl_da_skinc)=
    **icpl_da_skinc**
  - I
  - 0
  -
  - Coupling between bias of diurnal temperature amplitude inferred from data assimilation and skin conductivity
    - 0: off
    - 1: on; requires assimilation of T2M and cycling of a time-filtered weighted (with cosine of local time) temperature increment at the lowest model level
    - 2: as option 1, but additional adaptation of soil heat conductivity and heat capacity
  - init_mode=5

* - (initicon_nml-icpl_da_sfcfric)=
    **icpl_da_sfcfric**
  - I
  - 0
  -
  - Coupling between data assimilation and model parameters controlling surface friction (roughness length and SSO blocking tendency at lowest level).
    - 0: off
    - 1: on; requires assimilation of 10m-winds and cycling a time-filtered assimilation increment of absolute wind speed at the lowest model level; moreover, it is strongly recommended to use extpar data with full SSO information (generated in Feb. 2022 or later). Coupling is masked in large parts of Russia where the assimilation of 10m winds is blacklisted in the operational settings of 2022
    - 2: on without masking over Russia, to be combined with 10m wind assimilation without blacklisting
    - 3: same as 2, but reduced weight of assimilation increments occurring at PBL wind speeds below 7.5 m/s
  - init_mode=5

* - (initicon_nml-scalfac_da_sfcfric)=
    **scalfac_da_sfcfric**
  - R
  - 2.5
  -
  - Scaling factor for adaptive surface friction (see eqns. 3 and 4 in <https://doi.org/10.1002/qj.4535>)
  - icpl_da_sfcfric > 0

* - (initicon_nml-icpl_da_tkhmin)=
    **icpl_da_tkhmin**
  - I
  - 0
  -
  - Coupling between data assimilation and near-surface reduction profile for minimum vertical diffusion of heat
    - 0: off
    - 1: on
  - init_mode=5, icpl_da_sfcevap > 2 and icpl_da_skinc > 0

* - (initicon_nml-dt_filt)=
    **dt_filt**
  - R(2)
  - 2.5
  - days
  - Filtering time scale for filtered assimilation increments used for adaptive parameter tuning.
    1st value: standard time scale
    2nd value: time scale for t_avginc and rh_avginc (recommended tuning 5 days)
  - all icpl_da options

* - (initicon_nml-adjust_tso_tsnow)=
    **adjust_tso_tsnow**
  - L
  - .FALSE.
  -
  - If .TRUE., apply T increments for lowest model level also to snow and upper soil layers (full to upper 3 cm, half to 3-9 cm layer). Requires assimilation of T2M to be meaningful
  - init_mode=5

* - (initicon_nml-lconsistency_checks)=
    **lconsistency_checks**
  - L
  - .TRUE.
  -
  - If .FALSE., consistency checks for Analysis and First Guess fields are skipped. On default, checks are performed for _uuidOfHGrid_ and _validity time_.
  - init_mode=1,3,4,5,7

* - (initicon_nml-l_coarse2fine_mode)=
    **l_coarse2fine_mode**
  - L(n_dom)
  - .FALSE.
  -
  - If true, apply corrections for coarse-to-fine mesh interpolation to wind and temperature
  -

* - (initicon_nml-lp2cintp_incr)=
    **lp2cintp_incr**
  - L(n_dom)
  - .FALSE.
  -
  - If true, interpolate atmospheric data assimilation increments from parent domain.
    Can be specified separately for each nested domain; setting the first (global) entry to true activates the interpolation for all nested domains.
  - init_mode=5

* - (initicon_nml-lp2cintp_sfcana)=
    **lp2cintp_sfcana**
  - L(n_dom)
  - .FALSE.
  -
  - If true, interpolate atmospheric surface analysis data from parent domain.
    Can be specified separately for each nested domain; setting the first (global) entry to true activates the interpolation for all nested domains.
  - init_mode=5

* - (initicon_nml-ltile_init)=
    **ltile_init**
  - L
  - .FALSE.
  -
  - True: initialize tiled surface fields from a first guess coming from a run without tiles.
    Along coastlines and lake shores, a neighbor search is executed to fill the variables on previously non-existing land or water points with reasonable values. Should be combined with ltile_coldstart = .TRUE.
  - init_mode=1,5,7

* - (initicon_nml-ltile_coldstart)=
    **ltile_coldstart**
  - L
  - .FALSE.
  -
  - If true, tiled surface fields are initialized with tile-averaged fields from a previous run with tiles.
    A neighbor search is applied to subgrid-scale ocean points for SST and sea-ice fraction.
  - init_mode=1,5,7

* - (initicon_nml-lcouple_ocean_coldstart)=
    **lcouple_ocean_coldstart**
  - L
  - .TRUE.
  -
  - If true, initialize newly defined land points from ICON-O with default T and Q profiles.
  - is_coupled_mode=T

* - (initicon_nml-lvert_remap_fg)=
    **lvert_remap_fg**
  - L
  - .FALSE.
  -
  - If true, vertical remapping is applied to the atmospheric first-guess fields, whereas the analysis increments remain unchanged. The number of model levels must be the same for input and output fields, and the z_ifc (alias HHL) field pertaining to the input fields must be appended to the first-guess file.
  - init_mode=5

* - (initicon_nml-ifs2icon_filename)=
    **ifs2icon_filename**
  - C
  -
  -
  - Filename of IFS2ICON input file, default `<path>ifs2icon_R<nroot>B<jlev>_DOM <idom>.nc`. May contain the keywords `<path>` which will be substituted by `model_base_dir`, as well as `nroot`, `nroot0`, `jlev`, and `idom` defining the current patch.
  - init_mode=2

* - (initicon_nml-dwdfg_filename)=
    **dwdfg_filename**
  - C
  -
  -
  - Filename of DWD first-guess input file, default `<path>dwdFG_R<nroot>B<jlev>_DOM <idom>.nc`. May contain the keywords `<path>` which will be substituted by `model_base_dir`, as well as `nroot`, `nroot0`, `jlev`, and `idom` defining the current patch.
  - init_mode=1,3,5,7

* - (initicon_nml-dwdana_filename)=
    **dwdana_filename**
  - C
  -
  -
  - Filename of DWD analysis input file, default `<path>dwdana_R<nroot>B<jlev>_DOM <idom>.nc`. May contain the keywords `<path>` which will be substituted by `model_base_dir`, as well as `nroot`, `nroot0`, `jlev`, and `idom` defining the current patch.
  - init_mode=1,3,5

* - (initicon_nml-filetype)=
    **filetype**
  - I
  - -1 (undef.)
  -
  - One of CDI's FILETYPE_XXX constants. Possible values:
    2 (=FILETYPE_GRB2)
    4 (=FILETYPE_NC2).
    If this parameter has not been set, we try to determine the file type by its extension "*.grb*" or ".nc".
  -

* - (initicon_nml-check_fg(jg)%list)=
    **check_fg(jg)%list**
  - C(:)
  -
  -
  - In ICON a small subset of first guess input fields is declared 'optional', meaning that they are read in if present, but they are not mandatory to start the model. By adding optional fields to this list, they become mandatory for domain `jg`, such that the model aborts if any of them is missing. This list may include a subset of the optional first guess fields, or even the entire set of first guess fields. On default this list is empty, such that optional fields experience a cold-start initialization if they are missing and the model does not abort.
  - init_mode=1,5,7

* - (initicon_nml-check_ana(jg)%list)=
    **check_ana(jg)%list**
  - C(:)
  -
  -
  - List of mandatory analysis fields for domain `jg` that must be present in the analysis file. If these fields are not found, the model aborts. For all other analysis fields, the FG-fields will serve as fallback position.
  - init_mode=1,5

* - (initicon_nml-ana_varnames_map_ file)=
    **ana_varnames_map_ file**
  - C
  -
  -
  - Dictionary file which maps internal variable names onto GRIB2 shortnames or NetCDF var names. This is a text file with two columns separated by whitespace, where left column: ICON variable name, right column: GRIB2 short name or NetCDF var name.
  -

* - (initicon_nml-itype_vert_expol)=
    **itype_vert_expol**
  - I
  - 1
  -
  - Type of vertical extrapolation of initial data:
    - 1: Linear extrapolation (standard)
    - 2: Blend of linear extrapolation and simple climatology. Intended for upper-atmosphere simulations and specific settings can be found in [upatmo_nml](ref_buildrun_nml_upatmo_nml).
  - Requires: ivctype = 2; l_limited_area = .FALSE.

* - (initicon_nml-fire2d_filename)=
    **fire2d_filename**
  - C
  - 'gfas2d_emi_
    \<species>_
    \<gridfile>_
    \<yyyymmdd>.nc'
  -
  - Wildfire emission data sets for the \<species> bc, oc and so2. Possible keywords: \<species>, \<gridfile>, \<nroot>, \<nroot0>, \<jlev>, \<idom>, \<yyyymmdd>
  - Requires: i2daero_fire = 1

* - (initicon_nml-parallel_grib_decoding)=
    **parallel_grib_decoding**
  - L
  - .FALSE.
  -
  - Parallel decoding of GRIB2 input data by Work PEs:
    * .FALSE.: Workroot PE reads and decodes the records of a GRIB input file sequentially.
    * .TRUE.: Workroot PE reads the records and distributes them to all Work PEs for decoding in parallel.
    The second case requires input files of format GRIB2 with file extension '.grb', '.grb2', '.grib' or '.grib2'. In addition, a dictionary is mandatory (ana_varnames_map_file).  
    The data values packed in a GRIB file may have a bit depth up to bitsPerValue = 32. Note, however, that the data values are internally buffered in single precision (for reasons of efficiency) irrespective of the encoding accuracy, which may not deliver the maximum possible precision.  
    (A few corresponding timers are implemented and subordinated to the init_icon timer. Activating them requires timers_level > 8.)
  - ICON compiled with Fortran90 interfaces of package `ecCodes`;
    init_mode=1, 5, 7
```



Defined and used in: {{ '[src/namelists/mo_initicon_nml.f90]({}/src/namelists/mo_initicon_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_interpol_nml)=
## interpol_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (interpol_nml-l_intp_c2l)=
    **l_intp_c2l**
  - L
  - .TRUE.
  -
  - DEPRECATED
  -

* - (interpol_nml-l_mono_c2l)=
    **l_mono_c2l**
  - L
  - .TRUE.
  -
  - Monotonicity can be enforced by demanding that the interpolated value is not higher or lower than the stencil point values.
  -

* - (interpol_nml-llsq_high_consv)=
    **llsq_high_consv**
  - L
  - .TRUE.
  -
  - conservative (T) or non-conservative (F) least-squares reconstruction for high order transport
  -

* - (interpol_nml-lsq_high_ord)=
    **lsq_high_ord**
  - I
  - 3
  -
  - polynomial order of high order least-squares reconstruction for tracer transport
    - 1: linear
    - 2: quadratic
    - 3: cubic
  - ihadv_tracer > 2

* - (interpol_nml-llsq_lin_consv)=
    **llsq_lin_consv**
  - L
  - .FALSE.
  -
  - conservative (T) or non-conservative (F) least-squares reconstruction for 2nd order (linear) transport
  -

* - (interpol_nml-nudge_efold_width)=
    **nudge_efold_width**
  - R
  - 2.0
  -
  - e-folding width (in units of cell rows) for lateral boundary nudging coefficient. This switch and the following two pertain to one-way nesting and limited-area mode
  -

* - (interpol_nml-nudge_max_coeff)=
    **nudge_max_coeff**
  - R
  - 0.02
  -
  - Maximum relaxation coefficient for lateral boundary nudging. Recommended range of values for limited-area mode is 0.06 -- 0.075. The range of validity is [0 -- 0.2].
    **Please note** that the user value is internally multiplied by 5.
  -

* - (interpol_nml-nudge_zone_width)=
    **nudge_zone_width**
  - I
  - 8
  -
  - Total width (in units of cell rows) for lateral boundary nudging zone. For the limited-area mode, a minimum of 10 is recommended. If < 0 the patch boundary_depth_index is used.
  -

* - (interpol_nml-rbf_dim_c2l)=
    **rbf_dim_c2l**
  - I
  - 10
  -
  - stencil size for direct lon-lat interpolation:
    - 4: nearest neighbor,
    - 13: vertex stencil,
    - 10: edge stencil.
  -

* - (interpol_nml-rbf_scale_mode_ll)=
    **rbf_scale_mode_ll**
  - I
  - 2
  -
  - Specifies, how the RBF shape parameter is determined for lon-lat interpolation.
    1 : lookup table based on grid level
    2 : determine automatically.
    So far, this routine only estimates the smallest value for the shape parameter for which the Cholesky is likely to succeed in floating point arithmetic. 3 : explicitly set shape parameter in each output namelist (namelist parameter [rbf_scale](output_nml-rbf_scale)).
  -

* - (interpol_nml-rbf_vec_kern_c)=
    **rbf_vec_kern_c**
  - I
  - 1
  -
  - Kernel type for reconstruction at cell centres:
    - 1: Gaussian
    - 3: inverse multiquadric
  -

* - (interpol_nml-rbf_vec_kern_e)=
    **rbf_vec_kern_e**
  - I
  - 3
  -
  - Kernel type for reconstruction at edges:
    - 1: Gaussian
    - 3: inverse multiquadric
  -

* - (interpol_nml-rbf_vec_kern_ll)=
    **rbf_vec_kern_ll**
  - I
  - 1
  -
  - Kernel type for reconstruction at lon-lat-points:
    - 1: Gaussian
    - 3: inverse multiquadric
  -

* - (interpol_nml-rbf_vec_kern_v)=
    **rbf_vec_kern_v**
  - I
  - 1
  -
  - Kernel type for reconstruction at vertices:
    - 1: Gaussian
    - 3: inverse multiquadric
  -

* - (interpol_nml-rbf_vec_scale_c)=
    **rbf_vec_scale_c**
  - R(n_dom)
  - resolution-dependent
  -
  - Scale factor for RBF reconstruction at cell centres
  -

* - (interpol_nml-rbf_vec_scale_e)=
    **rbf_vec_scale_e**
  - R(n_dom)
  - resolution-dependent
  -
  - Scale factor for RBF reconstruction at edges
  -

* - (interpol_nml-rbf_vec_scale_v)=
    **rbf_vec_scale_v**
  - R(n_dom)
  - resolution-dependent
  -
  - Scale factor for RBF reconstruction at vertices
  -

* - (interpol_nml-rbf_coeffs_filename)=
    **rbf_coeffs_filename**
  - C
  - rbf_coeffs_dom\<i_dom>.nc
  -
  - Array of filenames to be used for reading/writing RBF coefficients if `lrbf_read=.TRUE. or lrbf_write=.TRUE.`. Each entry corresponds to the RBF coefficients for the grids defined in `dynamics_grid_filename`.
  -

* - (interpol_nml-lrbf_read)=
    **lrbf_read**
  - L
  - .FALSE.
  -
  - Flag. If .TRUE. then the RBF coefficients are read from file `rbf_coeffs_filename` instead of reconstructing. Coefficients are only read from file if scaling parameters, and gridsize are equal to those in input file. Else the model fails.
  -

* - (interpol_nml-lrbf_write)=
    **lrbf_write**
  - L
  - .FALSE.
  -
  - Flag. If .TRUE. then the RBF coefficients are written to file `rbf_coeffs_filename`.
  -

* - (interpol_nml-support_baryctr_intp)=
    **support_baryctr_intp**
  - L
  - .FALSE.
  -
  - Flag. If .FALSE. barycentric interpolation is replaced by a fallback interpolation.
  -

* - (interpol_nml-lreduced_nestbdry_stencil)=
    **lreduced_nestbdry_stencil**
  - L
  - .FALSE.
  -
  - Flag. If .TRUE. then the nest boundary points are taken out from the lat-lon interpolation stencil.
  -
```



Defined and used in: {{ '[src/namelists/mo_interpol_nml.f90]({}/src/namelists/mo_interpol_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_io_nml)=
## io_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (io_nml-lkeep_in_sync)=
    **lkeep_in_sync**
  - L
  - .FALSE.
  -
  - Sync output stream with file on disk after each timestep
  -

* - (io_nml-dt_diag)=
    **dt_diag**
  - R
  - 86400.0
  - s
  - diagnostic integral output interval
  - [output](run_nml-output) = "totint"`

* - (io_nml-dt_checkpoint)=
    **dt_checkpoint**
  - R
  - 0
  - s
  - Time interval for writing restart files. Note that if the value of dt_checkpoint resulting from model default or user's specification is longer than dt_restart, it will be reset (by the model) to dt_restart so that at least one restart file is generated during the restart cycle.
  - [output](run_nml-output) /= "none"`

* - (io_nml-inextra_2d)=
    **inextra_2d**
  - I
  - 0
  -
  - Number of extra 2D Fields for diagnostic/debugging output.
  -

* - (io_nml-inextra_3d)=
    **inextra_3d**
  - I
  - 0
  -
  - Number of extra 3D Fields for diagnostic/debugging output.
  -

* - (io_nml-lflux_avg)=
    **lflux_avg**
  - L
  - .TRUE.
  -
  - if .FALSE. the output fluxes are accumulated
    from the beginning of the run
    if .TRUE. the output fluxes are average values
    from the beginning of the run, except of
    TOT_PREC that would be accumulated
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-itype_hzerocl)=
    **itype_hzerocl**
  - I
  - 1
  -
  - Specifies setting of hzerocl if no freezing level is found.
    - 1: Height of orography,
    - 2: -999.0_wp (undef),
    - 3: extrapolated value below ground (assuming -6.5 K/km).
  -

* - (io_nml-itype_pres_msl)=
    **itype_pres_msl**
  - I
  - 1
  -
  - Specifies method for computation of mean sea level pressure (and geopotential at pressure levels below the surface).
    - 1: GME-type extrapolation,
    - 2: stepwise analytical integration,
    - 3: current IFS method,
    - 4: IFS method with consistency correction
    - 5: New DWD method constituting a mixture between IFS and old GME method (departure level for downward extrapolation between 10 m and 150 m AGL depending on elevation)
  -

* - (io_nml-itype_rh)=
    **itype_rh**
  - I
  - 1
  -
  - Specifies method for computation of relative humidity
    - 1: WMO-type: water only (e_s=e_s_water),
    - 2: IFS-type: mixed phase (water and ice),
    - 3: IFS-type with clipping ({math}`\mathrm{rh}\leq 100`)
  -

* - (io_nml-gust_interval)=
    **gust_interval**
  - R(n_dom)
  - 3600.0
  - s
  - Interval over which wind gusts are maximized
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-ff10m_interval)=
    **ff10m_interval**
  - R(n_dom)
  - 600.0
  - s
  - Interval over which 10-m winds are averaged (and used as basis for the gust diagnosis)
  - itune_gust_diag=4

* - (io_nml-celltracks_interval)=
    **celltracks_interval**
  - R(n_dom)
  - 3600.0
  - s
  - Interval over which celltrack variables are maximized (lpi_max, uh_max, vorw_ctmax, w_ctmax, tcond_max, tcond10_max, dbz_ctmax, tot_pr_max)
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-dt_celltracks)=
    **dt_celltracks**
  - R(n_dom)
  - 120.0
  - s
  - Interval at which celltrack variables except lpi (uh, vorw, w_ct, tcond, tcond10) are calculated to determine uh_max, vorw_ctmax, w_ctmax, tcond_max, tcond10_max and dbz_ctmax
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-dt_lpi)=
    **dt_lpi**
  - R(n_dom)
  - 180.0
  - s
  - Interval at which lpi is calculated for determining lpi_max
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-dt_hailcast)=
    **dt_hailcast**
  - R(n_dom)
  - 180.0
  - s
  - Interval at which hailcast is called for determining dhail_mx, dhail_sd, dhail_av
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-wdur_min_hailcast)=
    **wdur_min_hailcast**
  - R(n_dom)
  - 900.0
  - s
  - Minimal updraft persistence per column for hailcast to be activated
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-dt_radar_dbz)=
    **dt_radar_dbz**
  - R(n_dom)
  - 120.0
  - s
  - Interval at which radar reflectivity is calculated for determining dbz_ctmax
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-force_calc_optvar)=
    **force_calc_optvar**
  - I(n_dom)
  - 0
  -
  - Allows to force the computation of optional diagnostics in domains where no output is written, e.g. to have valid fields for nest boundary interpolation. By default, the computations are triggered by the output namelists of a given model domain. Setting, for instance, the second entry to 3 means that the output namelists of domain 3 are used to trigger the optional diagnostics in domain 2 as well.
    **Caution: If the output fields written in a nested domain are a subset of the fields written in the parent domain, using this option will cause a failure!**
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-precip_interval)=
    **precip_interval**
  - C(n_dom)
  - "P01Y"
  -
  - Interval over which precipitation variables are accumulated (rain_gsp, snow_gsp, graupel_gsp, ice_gsp, hail_gsp, prec_gsp, rain_con, snow_con, prec_con, tot_prec, prec_con_rate_avg, prec_gsp_rate_avg, tot_prec_rate_avg)
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-totprec_d_interval)=
    **totprec_d_interval**
  - C(n_dom)
  - "PT01H"
  -
  - Interval over which the special precipitation variable tot_prec_d is accumulated, which can be output alongside or alternatively to tot_prec and enables a different accumulation time for this field than precip_interval.
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-maxt_interval)=
    **maxt_interval**
  - C(n_dom)
  - "PT06H"
  -
  - Interval over which max/min 2-m temperatures are calculated
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-runoff_interval)=
    **runoff_interval**
  - C(n_dom)
  - "P01Y"
  -
  - Interval over which surface and soil water runoff are accumulated
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-sunshine_interval)=
    **sunshine_interval**
  - C(n_dom)
  - "P01Y"
  -
  - Interval over which sunshine duration is accumulated
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-itype_dursun)=
    **itype_dursun**
  - I
  - 0
  -
  - Type of sunshine. 0 for WMO standard and for sunshine duration counted if >120W/m{math}`^2`. In the case of type 1 (this is the MeteoSwiss definition) the sunshine duration is counted only if >200W/m{math}`^2`
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-wshear_uv_heights)=
    **wshear_uv_heights**
  - R(max_wshear)
    max_wshear=10
  - 1000.0, 3000.0, 6000.0
  -
  - List of height levels (m AGL) for which the vertical windshear output variables "wshear_u" and "wshear_v" are to be output.
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-srh_heights)=
    **srh_heights**
  - R(max_srh)
    max_srh=10
  - 1000.0, 3000.0
  -
  - List of height levels (m AGL) for which the storm relative helicity "srh" is to be output. "srh" is a vertical integral from the ground to a certain height. The listed height levels denote different upper bounds for this integration.
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-echotop_meta_)=
    **echotop_meta%...**
  - TYPE(n_dom)
    R(1)
    R(max_echotop)
    max_echotop=10
  - 3600.0
    (/18.0,25.0,35.0/)
  - s
    dBZ
  - Derived type to define properties of radar reflectivity echotops for each domain. Two types of echotops are available: minimum pressure (`'echotop'`) and maximum height (`'echotopinm'`) during a given time interval where a given reflectivity threshold is exeeded. Takes effect if 'echotop' and/or 'echotopinm' is/are present in the ml_varlist of any domain-specific namelist [output_nml](ref_buildrun_nml_output_nml).
  The derived type contains the echotop properties which are listed to the left, along with their defaults and units.
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-echotop_meta_-time_interval)=
    ***...time_interval**
  - TYPE(n_dom)
    R(1)
    R(max_echotop)
    max_echotop=10
  - 3600.0
    (/18.0,25.0,35.0/)
  - s
    dBZ
  - Time interval [s] over which echotops are calculated
    You have to specify properties for each domain separately, e.g.
    echotop_meta(1)%time_interval=3600.0
    echotop_meta(2)%time_interval=1800.0
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-echotop_meta_-dbzthresh)=
    **...dbzthresh**
  - TYPE(n_dom)
    R(1)
    R(max_echotop)
    max_echotop=10
  - 3600.0
    (/18.0,25.0,35.0/)
  - s
    dBZ
  - List of reflectivity thresholds [dBZ] for which echotops shall be computed
    You have to specify properties for each domain separately, e.g.
    echotop_meta(1)%dbzthresh=19.0,25.0,35.0,46.0
    echotop_meta(2)%dbzthresh=27.0,36.0
  - [iforcing](run_nml-iforcing)=3

* - (io_nml-output_nml_dict)=
    **output_nml_dict**
  - C
  - ' '
  -
  - File containing the mapping of variable names to the internal ICON names. May contain the keyword `<path>` which will be substituted by `model_base_dir`.
    The format of this file:
    One mapping per line, first the name as given in the `ml_varlist`, `hl_varlist`, `pl_varlist` or `il_varlist` parameters, then the internal ICON name, separated by an arbitrary number of blanks. The line may also start and end with an arbitrary number of blanks. Empty lines or lines starting with `#` are treated as comments.
    Names not covered by the mapping are used as they are.
  - [output_nml](ref_buildrun_nml_output_nml) namelists

* - (io_nml-linvert_dict)=
    **linvert_dict**
  - L
  - .FALSE.
  -
  - If .TRUE., columns in dictionary file output_nml_dict are evaluated in inverse order.
    This allows using the same dictionary file as for input (ana_varnames_map_file from parallel_grib_decoding).
  -

* - (io_nml-netcdf_dict)=
    **netcdf_dict**
  - C
  - ' '
  -
  - File containing the mapping from internal names to names written to NetCDF. May contain the keyword `<path>` which will be substituted by `model_base_dir`.
    The format of this file:
    One mapping per line, first the name written to NetCDF, then the internal name, separated by an arbitrary number of blanks (_inverse to the definition of `output_nml_dict_`). The line may also start and end with an arbitrary number of blanks. Empty lines or lines starting with `#` are treated as comments.
    Names not covered by the mapping are output as they are.
    Note that the specification of output variables, e.g.
    in `ml_varlist`, is independent from this renaming, see the namelist parameter `output_nml_dict` for this.
  - [output_nml](ref_buildrun_nml_output_nml) namelists, NetCDF output

* - (io_nml-lnetcdf_flt64_output)=
    **lnetcdf_flt64_output**
  - L
  - .FALSE.
  -
  - If .TRUE. floating point variable output in NetCDF files is written in 64-bit instead of 32-bit accuracy.
  -

* - (io_nml-restart_file_type)=
    **restart_file_type**
  - I
  - 4
  -
  - Type of restart file. One of CDI's FILETYPE_XXX. So far, only 4 (=FILETYPE_NC2) is allowed
  -

* - (io_nml-restart_write_mode)=
    **restart_write_mode**
  - C
  - " "
  -
  - Restart read/write mode.
    Allowed settings (character strings!) are listed below.
  -

* - (io_nml-nrestart_streams)=
    **nrestart_streams**
  - I
  - 1
  -
  - When using the restart write mode "dedicated procs multifile", it is possible to split the restart output into several files, as if nrestart_streams * num_io_procs restart processes were involved. This speeds up the read-in process, since all the files may then be read in parallel.
  - `restart_write_mode` = `"dedicated procs multifile"`

* - (io_nml-checkpoint_on_demand)=
    **checkpoint_on_demand**
  - L
  - F
  -
  - `.TRUE.` allows checkpointing (followed by stopping) during runtime triggered by a file named 'stop_icon' in the working directory. In addition, a file named 'ready_for_checkpoint' is generated in the working directory once the model is ready for checkpointing, i.e. after the end of the setup phase, or, if applicable, the end of the IAU phase.
  - Combination with `restart_write_mode` = `"joint procs multifile"` is strongly recommended

* - (io_nml-lmask_boundary)=
    **lmask_boundary**
  - L(n_dom)
  - F
  -
  - Set to `.TRUE.`, if interpolation zone should be masked in triangular output.
  -

* - (io_nml-ldiagnose_tke)=
    **ldiagnose_tke**
  - L(n_dom)
  - F
  -
  - Set to .TRUE. to calculate grid-scale TKE based on temporal wind variance over interval specified with gstke_interval. Also calculates SGS TKE and EDR averaged over interval, and activates additional diagnostic variables for EDR and turbulent length scale.
  - [iforcing](run_nml-iforcing)=3 Time-average SGS TKE, EDR and length scale output currently implemented only for inwp_turb=1 (Turbdiff). Fields will be zero otherwise.

* - (io_nml-gstke_interval)=
    **gstke_interval**
  - R(n_dom)
  - 3600.0
  - s
  - Time interval in seconds over which to calculate temporal wind variance for grid-scale TKE.
  - [iforcing](run_nml-iforcing)=3
```



### Restart read/write mode:


Allowed settings for `restart_write_mode` are:

`sync`

 - 'Old' synchronous mode. PE <br># 0 reads and writes restart files. All other PEs have to wait.

`async`

 - Asynchronous restart writing:
   Dedicated PEs (`num_restart_proc > 0`) write restart files
   while the simulation continues.  Restart PEs can only parallelize
   over different patches. --- Read-in: PE # 0 reads while other PEs
   have to wait.

`joint procs multifile`

 - All worker PEs write restart files to a dedicated
   directory. Therefore, the directory itself is called the restart
   file.  The information is stored in a way that it can be read back
   into the model independent from the processor count and the domain
   decomposition. --- Read-in: All worker PEs read the data in parallel.

`dedicated procs multifile`

 - In this case, all the restart data is first transferred to memory
   buffers in dedicated restart writer PEs.  After that, the work
   processes carry on with their work immediately, while the restart
   writers perform the actual restart writing asynchronously.  Restart
   PEs can parallelize over patches and horizontal indices.
   --- Read-in: All worker PEs read the data in parallel..


`""`

 - Fallback mode. <br> If num_restart_proc ` == 0`,
   then this behaves like `sync`, otherwise like `async`.



### Some notes on the output of optional diagnostics


{math}`\blacksquare` _How can I switch on the output of one of the available diagnostics?_

Let us assume that you would like to output _potential vorticity_
(see table of available diagnostics below) on model levels.
Simply add the following element to the desired output namelist (see [output_nml](ref_buildrun_nml_output_nml) ) in your run script:

```
&output_nml
...
ml_varlist = ..., 'pv'
...
/
```

Please note that the output of some diagnostics is restricted to the NWP mode
([iforcing](run_nml-iforcing) = inwp = 3, see column "Scope" in the table below).

{math}`\blacksquare` _Which optional diagnostics are currently available for output?_

Here is a table of the available diagnostics and some additional information on them.


```{list-table}
:header-rows: 1

* - (Some notes on the output of optional diagnostics-Short name)=
    **Short name**
  - **Long name**
  - **Unit**
  - **Scope**
  - **Shape**
  - **Specifications in [io_nml](ref_buildrun_nml_io_nml)**
  - **Place of computation in source code** **

* - rh
  - relative humidity
  - %
  - [iforcing](run_nml-iforcing) = inwp = 3
  - 3d
  - itype_rh
  - [1]

* - pv
  - potential vorticity
  - K m2 kg-1 s-1
  - [iforcing](run_nml-iforcing) = inwp
  - 3d
  -
  - [2]

* - sdi2
  - supercell detection index (SDI2)
  - s-1
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - lpi
  - lightning potential index (LPI)
  - J kg-1
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - lpi_max
  - lightning potential index, maximum during prescribed time interval
  - J kg-1
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  - celltracks_interval
    dt_lpi
  - [2]

* - ceiling
  - ceiling height
  - m
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - vis
  - near-surface horizontal visibility
  - m
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - hbas_sc
  - cloud base above msl, shallow convection
  - m
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - htop_sc
  - cloud top above msl, shallow convection
  - m
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - twater
  - total column-integrated water
  - kg m-2
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - q_sedim
  - specific content of precipitation particles
  - kg kg-1
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - tcond_max
  - total column-integrated condensate, maximum during prescribed time interval
  - kg m-2
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  - celltracks_interval
    dt_celltracks
  - [2]

* - tcond10_max
  - total column-integrated condensate above z(T=-10 degC), maximum during prescribed time interval
  - kg m-2
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  - celltracks_interval
    dt_celltracks
  - [2]

* - uh_max
  - updraft helicity, maximum during prescribed time interval
  - m2 s-2
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  - celltracks_interval
    dt_celltracks
  - [2]

* - vorw_ctmax
  - maximum rotation amplitude during prescribed time interval
  - s-1
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  - celltracks_interval
    dt_celltracks
  - [2]

* - w_ctmax
  - maximum updraft track during prescribed time interval
  - m s-1
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  - celltracks_interval
    dt_celltracks
  - [2]

* - dbz
  - radar reflectivity
  - dBZ
  - [iforcing](run_nml-iforcing) = inwp
  - 3d
  -
  - [2]

* - dbz_cmax
  - column maximum reflectivity
  - dBZ
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - dbz_850
  - reflectivity in approx. 850 hPa
  - dBZ
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - dbz_ctmax
  - column and time maximum reflectivity during prescribed time interval
  - dBZ
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  - celltracks_interval
    dt_radar_dbz
  - [2]

* - echotop
  - minimum pressure of exceeding radar reflectivity threshold during prescribed time interval
  - Pa
  - [iforcing](run_nml-iforcing) = inwp
  - 3d
  - celltracks_interval
    echotop_meta
  - [2]

* - echotopinm
  - maximum height of exceeding radar reflectivity threshold during prescribed time interval
  - m
  - [iforcing](run_nml-iforcing) = inwp
  - 3d
  - celltracks_interval
    echotop_meta
  - [2]

* - pres_msl
  - mean sea level pressure
  - Pa
  -
  - 2d
  - itype_pres_msl
  - [3]

* - omega
  - vertical (pressure) velocity
  - Pa s-1
  -
  - 3d
  -
  - [2]

* - tot_prec_d
  - total accumulated precipitation during a different time interval compared to tot_prec
  - kg m-2
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  - totprec_d_interval
  - [1], [4], [5]

* - tot_pr_max
  - maximum total precipitation rate during prescribed time interval
  - kg m-2 s-1
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  - celltracks_interval
  - [4]

* - lapse_rate
  - temperature gradient between 500 and 850 hPa
  - K m-1
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - mconv
  - low level horizontal moisture convergence averaged over 0-1000 m AGL layer based on specific humidity, {math}`\frac{1}{1\,\text{km}}\int_{0}^{1\,\text{km AGL}}\nabla_h\cdot(q_v\vec v_h)\,dz`
  - s-1
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - wshear_u
  - difference of U component between certain heights ("wshear_uv_heights") AGL and the lowest model level
  - m s-1
  - [iforcing](run_nml-iforcing) = inwp
  - 3d
  - wshear_uv_heights
  - [2]

* - wshear_v
  - difference of V component between certain heights ("wshear_uv_heights") AGL and the lowest model level
  - m s-1
  - [iforcing](run_nml-iforcing) = inwp
  - 3d
  - wshear_uv_heights
  - [2]

* - srh
  - storm relative helicity considering storm motion estimate of Bunkers et al. (2000) for right-movers. srh is a vertical intergal up to a certain height AGL and may be output for different upper bounds ("srh_heights").
  - m2 s-2
  - [iforcing](run_nml-iforcing) = inwp
  - 3d
  - srh_heights
  - [2]

* - cape_mu
  - approximate value of the most unstable CAPE considering a test parcel from the height level with largest equivalent potential temperature between the ground and 3000 m AGL
  - J kg-1
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - cin_mu
  - approximate value of the most unstable CIN consistent to cape_mu
  - J kg-1
  - [iforcing](run_nml-iforcing) = inwp
  - 2d
  -
  - [2]

* - dm_hail_s
  - mean mass diameter of hail at the surface
  - m
  - [iforcing](run_nml-iforcing) = inwp, [inwp_gsp](nwp_phy_nml-inwp_gscp) = 4,5,6,7
  - 2d
  -
  - [2]

* - dm_hail_max_s
  - Time maximum of mean mass diameter of hail at the surface
  - m
  - [iforcing](run_nml-iforcing) = inwp, [inwp_gsp](nwp_phy_nml-inwp_gscp) = 4,5,6,7
  - 2d
  - celltracks_interval
  - [2]

* - demax_hail_s
  - estimated maximum hail diameter at the surface
  - m
  - [iforcing](run_nml-iforcing) = inwp, [inwp_gsp](nwp_phy_nml-inwp_gscp) = 4,5,6,7
  - 2d
  - `[nwp_tuning_nml](ref_buildrun_nml_nwp_phy_nml)`: [tune_demax_hail_s](nwp_tuning_nml-tune_demax_hail_s), [prhthresh_demax_hail_s](nwp_tuning_nml-prhthresh_demax_hail_s), [kefthresh_demax_hail_s](nwp_tuning_nml-kefthresh_demax_hail_s), [qnhthresh_demax_hail_s](nwp_tuning_nml-qnhthresh_demax_hail_s), [lwindeffect_kef_hail_s](nwp_tuning_nml-lwindeffect_kef_hail_s)
  - [2]

* - demax_hail_tmax_s
  - time maximum of estimated maximum hail diameter at the surface
  - m
  - [iforcing](run_nml-iforcing) = inwp, [inwp_gsp](nwp_phy_nml-inwp_gscp) = 4,5,6,7
  - 2d
  - celltracks_interval
    `[nwp_tuning_nml](ref_buildrun_nml_nwp_phy_nml)`: [tune_demax_hail_s](nwp_tuning_nml-tune_demax_hail_s), [prhthresh_demax_hail_s](nwp_tuning_nml-prhthresh_demax_hail_s), [kefthresh_demax_hail_s](nwp_tuning_nml-kefthresh_demax_hail_s), [qnhthresh_demax_hail_s](nwp_tuning_nml-qnhthresh_demax_hail_s), [lwindeffect_kef_hail_s](nwp_tuning_nml-lwindeffect_kef_hail_s)
  - [2]

* - kef_hail_s
  - hail kinetic energy flux at the surface
  - W m-2
  - [iforcing](run_nml-iforcing) = inwp, [inwp_gsp](nwp_phy_nml-inwp_gscp) = 4,5,6,7
  - 2d
  -
  - [2]

* - kef_hail_max_s
  - time maximum of hail kinetic energy flux at the surface
  - W m-2
  - [iforcing](run_nml-iforcing) = inwp, [inwp_gsp](nwp_phy_nml-inwp_gscp) = 4,5,6,7
  - 2d
  - celltracks_interval
  - [4]

* - ke_hail_s
  - accumulated hail kinetic energy at the surface
  - J m-2
  - [iforcing](run_nml-iforcing) = inwp, [inwp_gsp](nwp_phy_nml-inwp_gscp) = 4,5,6,7
  - 2d
  - precip_interval
  - [4]
```



*:  To be used in [output_nml](ref_buildrun_nml_output_nml).

**: The keys, {[1]}, {[2]}, etc., are itemized under the following point.


{math}`\blacksquare` _Where can I find more about the computation of the diagnostics in the source code?_

_As for the ICON model component of the non-hydrostatic atmosphere:_

Each optional diagnostic has its own switch in the source code of ICON
which is set to `.TRUE.` if the diagnostic is found in one of the [output_nml](ref_buildrun_nml_output_nml) in your run script.
This configuration can be found in the module: <br>
{{ '[src/configure_model/mo_io_config.f90]({}/src/configure_model/mo_io_config.f90)'.format(base_url) }}.

Further information on the metadata of the diagnostics can be found
in their allocation area.
For the diagnostics that are meant for the NWP mode of ICON
([iforcing](run_nml-iforcing) = inwp = 3, see column "Scope" in table above),
the allocation takes place in:<br>
{{ '[src/atm_phy_nwp/mo_nwp_phy_state.f90]({}/src/atm_phy_nwp/mo_nwp_phy_state.f90)'.format(base_url) }}.<br>
Optional diagnostics with unrestricted scope are allocated in:<br>
{{ '[src/atm_dyn_iconam/mo_nonhydro_state.f90]({}/src/atm_dyn_iconam/mo_nonhydro_state.f90)'.format(base_url) }}.

The job control of the computation and output of most of the optional diagnostics
is organized by the post-processing scheduler:<br>
{{ '[src/atm_dyn_iconam/mo_pp_scheduler.f90]({}/src/atm_dyn_iconam/mo_pp_scheduler.f90)'.format(base_url) }}, <br>
{{ '[src/atm_dyn_iconam/mo_pp_tasks.f90]({}/src/atm_dyn_iconam/mo_pp_tasks.f90)'.format(base_url) }}, <br>
and integrated into the main time loop in:<br>
{{ '[src/atm_dyn_iconam/mo_nh_stepping.f90]({}/src/atm_dyn_iconam/mo_nh_stepping.f90)'.format(base_url) }}.<br>
The job control of a small portion of the diagnostics
is organized in: <br>
{{ '[src/atm_phy_nwp/mo_nwp_diagnosis.f90]({}/src/atm_phy_nwp/mo_nwp_diagnosis.f90)'.format(base_url) }}.

Finally, the computation of the individual diagnostics can be found
in the following modules
(the assignment of the keys, {[1]}, {[2]}, etc., to the respective diagnostic
is found in the column "Place of computation in source code" of table above):

 - [1]: {{ '[src/atm_phy_nwp/mo_util_phys.f90]({}/src/atm_phy_nwp/mo_util_phys.f90)'.format(base_url) }}
 - [2]: {{ '[src/atm_phy_nwp/mo_opt_nwp_diagnostics.f90]({}/src/atm_phy_nwp/mo_opt_nwp_diagnostics.f90)'.format(base_url) }}
 - [3]: {{ '[src/atm_phy_nwp/mo_nh_diagnose_pmsl.f90]({}/src/atm_phy_nwp/mo_nh_diagnose_pmsl.f90)'.format(base_url) }}
 - [4]: {{ '[src/atm_phy_nwp/mo_nwp_gscp_interface.f90]({}/src/atm_phy_nwp/mo_nwp_gscp_interface.f90)'.format(base_url) }}
 - [5]: {{ '[src/atm_phy_nwp/mo_nwp_diagnosis.f90]({}/src/atm_phy_nwp/mo_nwp_diagnosis.f90)'.format(base_url) }}


Defined and used in: {{ '[src/namelists/mo_io_nml.f90]({}/src/namelists/mo_io_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_les_nml)=
## les_nml

Parameters for LES turbulence scheme, valid for inwp_turb=5.

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (les_nml-sst)=
    **sst**
  - R
  - 300
  - K
  - sea surface temperature for idealized LES simulations
  - isrfc_type=5,4

* - (les_nml-shflx)=
    **shflx**
  - R
  - 0.1
  - Km/s
  - Kinematic sensible heat flux at surface
  - isrfc_type = 2

* - (les_nml-lhflx)=
    **lhflx**
  - R
  - 0
  - m/s
  - Kinematic latent heat flux at surface
  - isrfc_type = 2

* - (les_nml-isrfc_type)=
    **isrfc_type**
  - I
  - 1
  -
  - surface type
    - 0: No fluxes and zero shear stress
    - 1: TERRA land physics
    - 2: fixed surface fluxes
    - 3: fixed buoyancy fluxes
    - 4: RICO test case
    - 5: fixed SST
    - 6: time varying SST and qv_s with prescribed roughness length for semi-idealized setups
    - 7,8,9: various fixed fluxed options
    - 10: surface conditions as used for SCM
  -

* - (les_nml-ufric)=
    **ufric**
  - R
  - -999
  - m/s
  - friction velocity for idealized LES simulations; if < 0 then it is automatically diagnosed
  -

* - (les_nml-psfc)=
    **psfc**
  - R
  - -999
  - Pa
  - surface pressure for idealized LES simulations; if < 0 then it uses the surface pressure from dynamics
  -

* - (les_nml-min_sfc_wind)=
    **min_sfc_wind**
  - R
  - 1.0
  - m/s
  - Minimum surface wind for surface layer useful in the limit of free convection
  -

* - (les_nml-is_dry_cbl)=
    **is_dry_cbl**
  - L
  - .FALSE.
  -
  - switch for dry convective boundary layer simulations
  -

* - (les_nml-smag_constant)=
    **smag_constant**
  - R
  - 0.23
  -
  - Smagorinsky constant
  -

* - (les_nml-km_min)=
    **km_min**
  - R
  - 0.0
  -
  - Minimum turbulent viscosity
  -

* - (les_nml-smag_coeff_type)=
    **smag_coeff_type**
  - I
  - 1
  -
  - choose type of coefficient setting:
    - 1: Smagorinsky model (default)
    - 2: set coeff. externally by Km_ext, Kh_ext (for testing purposes, e.g. Straka et al. (1993))
  -

* - (les_nml-Km_ext)=
    **Km_ext**
  - R
  - 75.0
  - m{math}`^2`/s
  - externally set constant kinematic viscosity
  - smag_coeff_type=2

* - (les_nml-Kh_ext)=
    **Kh_ext**
  - R
  - 75.0
  - m{math}`^2`/s
  - externally set constant diffusion coeff.
  - smag_coeff_type=2

* - (les_nml-max_turb_scale)=
    **max_turb_scale**
  - R
  - 300.0
  -
  - Asymtotic maximum turblence length scale (useful for coarse grid LES and when grid is vertically stretched)
  -

* - (les_nml-turb_prandtl)=
    **turb_prandtl**
  - R
  - 0.333333
  -
  - turbulent Prandtl number
  -

* - (les_nml-bflux)=
    **bflux**
  - R
  - 0.0007
  - m{math}`^2`/s{math}`^3`
  - buoyancy flux for idealized LES simulations (Stevens 2007)
  - isrfc_type=3

* - (les_nml-tran_coeff)=
    **tran_coeff**
  - R
  - 0.02
  - m/s
  - transfer coefficient near surface for idealized LES simulation (Stevens 2007)
  - isrfc_type=3

* - (les_nml-vert_scheme_type)=
    **vert_scheme_type**
  - I
  - 2
  -
  - type of time integration scheme in vertical diffusion
    - 1: explicit
    - 2: fully implicit
  -

* - (les_nml-sampl_freq_sec)=
    **sampl_freq_sec**
  - R
  - 60
  - s
  - sampling frequency in seconds for statistical (1D and 0D) output
  -

* - (les_nml-avg_interval_sec)=
    **avg_interval_sec**
  - R
  - 900
  - s
  - (time) averaging interval in seconds for 1D statistical output
  -

* - (les_nml-expname)=
    **expname**
  - C
  - ICOLES
  -
  - expname to name the statistical output file
  -

* - (les_nml-ldiag_les_out)=
    **ldiag_les_out**
  - L
  - .FALSE.
  -
  - Control for the statistical output in LES mode
  -

* - (les_nml-les_metric)=
    **les_metric**
  - L
  - .FALSE.
  -
  - Switch to turn on Smagorinsky diffusion with 3D metric terms to account for topography
  -
```



Defined and used in: {{ '[src/namelists/mo_les_nml.f90]({}/src/namelists/mo_les_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_limarea_nml)=
## limarea_nml

Scope: l_limited_area=.TRUE.

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (limarea_nml-itype_latbc)=
    **itype_latbc**
  - I
  - 0
  -
  - Type of lateral boundary nudging.
    - 0: constant lateral boundary conditions derived from the initial conditions,
    - 1: time-dependent lateral boundary conditions provided by an external source (IFS, COSMO or a coarser-resolution ICON run)
  -

* - (limarea_nml-dtime_latbc)=
    **dtime_latbc**
  - R(3)
  - -1.0
  - s
  - Time difference between two consecutive boundary data. (Upper bound for asynchronous read-in: 1 day = 86400 s.)
    Up to 3 values for bc intervals are allowed; specifying more than one value requires setting interval bounds via bcintv_endtime.
  - itype_latbc {math}`\ge` 1

* - (limarea_nml-bcintv_endtime)=
    **bcintv_endtime**
  - R(3)
  - -1.0
  - s
  - Upper interval bounds for dtime_latbc, relative time w.r.t. model start (required if more than one value is specified for dtime_latbc; in any case, one entry less is needed for bcintv_endtime than for dtime_latbc).
  - itype_latbc {math}`\ge` 1

* - (limarea_nml-init_latbc_from_fg)=
    **init_latbc_from_fg**
  - L
  - .TRUE.
  -
  - If .TRUE., take lateral boundary conditions for initial time from first guess (or analysis) field
  - itype_latbc {math}`\ge` 1

* - (limarea_nml-nudge_hydro_pres)=
    **nudge_hydro_pres**
  - L
  - .TRUE.
  -
  - If .TRUE., hydrostatic pressure is used to compute lateral boundary nudging (recommended if boundary conditions contain hydrostatic pressure, which is usually the case)
  - itype_latbc {math}`\ge` 1

* - (limarea_nml-fac_latbc_presbiascor)=
    **fac_latbc_presbiascor**
  - R
  - 0.0
  -
  - Scaling factor for pressure bias correction at lateral boundaries. Requires running in data assimilation cycle. Recommended value for activating the option is 1.
  - itype_latbc {math}`\ge` 1, init_mode=5

* - (limarea_nml-latbc_filename)=
    **latbc_filename**
  - C
  -
  -
  - Filename of boundary data input file, these files must be located in the `latbc_path` directory.  Default:
    `prepiconR<nroot>B<jlev>_<y><m><d><h>.nc`.  The filename may contain keyword tokens (day, hour, etc.) which will be automatically replaced during the run-time. See the table below for a list of allowed keywords.
  - itype_latbc = 1

* - (limarea_nml-latbc_path)=
    **latbc_path**
  - C
  - './'
  -
  - Absolute path to boundary data.
  - itype_latbc = 1

* - (limarea_nml-latbc_boundary_grid)=
    **latbc_boundary_grid**
  - C
  - ' '
  -
  - Grid file defining the lateral boundary. Empty string means: whole domain is read for the lateral boundary. This NetCDF grid file must contain two integer index arrays: `int global_cell_index(cell)`, `int global_edge_index(edge)`, both with attributes `nglobal` which contains the global size size of the non-sparse cells and edges.
  - itype_latbc = 1

* - (limarea_nml-latbc_varnames_map_ file)=
    **latbc_varnames_map_ file**
  - C
  -
  -
  - Dictionary file which maps internal variable names onto GRIB2 shortnames or NetCDF var names. This is a text file with two columns separated by whitespace, where left column: ICON variable name, right column: GRIB2 short name. This list contains variables that are to be read asynchronously for boundary data nudging in a HDCP2 simulation. All new boundary variables that in the future, would be read asynchronously. Need to be added to text file dict.latbc in run folder.
  - num_prefetch_proc=1

* - (limarea_nml-latbc_contains_qcqi)=
    **latbc_contains_qcqi**
  - L
  - .TRUE.
  -
  - Set to .FALSE. if there is no qc, qi in latbc data.
  -

* - (limarea_nml-nretries)=
    **nretries**
  - I
  - 0
  -
  - If LatBC data is unavailable: number of retries
  -

* - (limarea_nml-retry_wait_sec)=
    **retry_wait_sec**
  - I
  - 10
  -
  - If LatBC data is unavailable: idle wait seconds between retries
  -
```



Defined and used in: {{ '[src/namelists/mo_limarea_nml.f90]({}/src/namelists/mo_limarea_nml.f90)'.format(base_url) }}



Keyword substitution in boundary data filename (`latbc_filename`):

```{list-table}
:header-rows: 1

* - Keyword
  - Description

* - `<y>`
  - substituted by year (four digits)

* - `<m>`
  - substituted by month (two digits)

* - `<d>`
  - substituted by day (two digits)

* - `<h>`
  - substituted by hour (two digits)

* - `<min>`
  - substituted by minute (two digits)

* - `<sec>`
  - substituted by seconds (two digits)

* - `<ddhhmmss>`
  - substituted by a _relative_ day-hour-minute-second string.

* - `<dddhh>`
  - substituted by a _relative_ (three-digit) day-hour string.
```




(ref_buildrun_nml_lnd_nml)=
## lnd_nml


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (lnd_nml-nlev_snow)=
    **nlev_snow**
  - I
  - 2
  -
  - number of snow layers
  - lmulti_snow=.TRUE.

* - (lnd_nml-ntiles)=
    **ntiles**
  - I
  - 1
  -
  - number of tiles
  -

* - (lnd_nml-zml_soil)=
    **zml_soil**
  - R(:)
  - 0.005, 0.02, 0.06,
    0.18, 0.54, 1.62,
    4.86, 14.58
  - m
  - soil full layer depths
  - init_mode = 2, 3

* - (lnd_nml-czbot_w_so)=
    **czbot_w_so**
  - R
  - 2.5
  - m
  - thickness of the hydrological active soil layer
  -

* - (lnd_nml-lsnowtile)=
    **lsnowtile**
  - L
  - .FALSE.
  -
  - .TRUE.: consider snow-covered and snow-free tiles separately
  - ntiles > 1

* - (lnd_nml-frlnd_thrhld)=
    **frlnd_thrhld**
  - R
  - 0.05
  -
  - fraction threshold for creating a land grid point
  - ntiles > 1

* - (lnd_nml-frlake_thrhld)=
    **frlake_thrhld**
  - R
  - 0.05
  -
  - fraction threshold for creating a lake grid point
  - ntiles > 1

* - (lnd_nml-frsea_thrhld)=
    **frsea_thrhld**
  - R
  - 0.05
  -
  - fraction threshold for creating a sea grid point
  - ntiles > 1

* - (lnd_nml-frlndtile_thrhld)=
    **frlndtile_thrhld**
  - R
  - 0.05
  -
  - fraction threshold for retaining the respective tile for a grid point
  - ntiles > 1

* - (lnd_nml-lmelt)=
    **lmelt**
  - L
  - .TRUE.
  -
  - .TRUE. soil model with melting process
  -

* - (lnd_nml-lmelt_var)=
    **lmelt_var**
  - L
  - .TRUE.
  -
  - .TRUE. freezing temperature dependent on water content
  -

* - (lnd_nml-lana_rho_snow)=
    **lana_rho_snow**
  - L
  - .TRUE.
  -
  - .TRUE. take rho_snow-values from analysis file
  - init_mode=1

* - (lnd_nml-lmulti_snow)=
    **lmulti_snow**
  - L
  - .FALSE.
  -
  - .TRUE. for use of multi-layer snow model (default is single-layer scheme)
  -

* - (lnd_nml-l2lay_rho_snow)=
    **l2lay_rho_snow**
  - L
  - .FALSE.
  -
  - .TRUE. predict additional snow density for upper part of the snowpack, having a maximum depth of max_toplaydepth
  - lmulti_snow = .FALSE.

* - (lnd_nml-max_toplaydepth)=
    **max_toplaydepth**
  - R
  - 0.25
  - m
  - maximum depth of uppermost snow layer
  - lmulti_snow=.TRUE. or l2lay_rho_snow=.TRUE.

* - (lnd_nml-idiag_snowfrac)=
    **idiag_snowfrac**
  - I
  - 1
  -
  - Type of snow-fraction diagnosis:
    - 1: based on SWE only
    - 2: more advanced method used operationally
    - 20: same as "2", but with artificial reduction of snow fraction in case of melting snow (should be used only in combination with lsnowtile=.TRUE.)
  -

* - (lnd_nml-itype_snowevap)=
    **itype_snowevap**
  - I
  - 2
  -
  - Tuning of snow evaporation in vegetated areas:
    - 1: Tuning turned off
    - 2: First level of tuning without additional control variables
    - 3: Second level of tuning with additional I/O variables for snow age and maximum snow depth (should be used only if these additional variables are avaliable from the DWD assimilation cycle)
  - lsnowtile=.TRUE.

* - (lnd_nml-itype_lndtbl)=
    **itype_lndtbl**
  - I
  - 3
  -
  - Table values used for associating surface parameters to land-cover classes:
    - 1: defaults from extpar (GLC2000 and GLOBCOVER2009)
    - 2: Tuned version based on IFS values for globcover classes (GLOBCOVER2009 only)
    - 3: even more tuned operational version (GLOBCOVER2009 only)
    - 4: tuned version for new bare soil evaporation scheme (itype_evsl=4)
  -

* - (lnd_nml-itype_root)=
    **itype_root**
  - I
  - 2
  -
  - type of root density distribution
    - 1: constant
    - 2: exponential
  -

* - (lnd_nml-itype_evsl)=
    **itype_evsl**
  - I
  - 2
  -
  - type of bare soil evaporation parameterization
    - 2: BATS scheme, Dickinson (1984)
    - 3: ISBA scheme, Noilhan and Planton (1989)
    - 4: Resistance-based formulation by Schulz and Vogel (2020)
    - 5: same as 4, but uses the minimum evaporation resistance (default set by cr_bsmin) instead of c_soil for tuning; the namelist parameter c_soil is ignored in this case, and a value of 2 is used internally
  -

* - (lnd_nml-itype_trvg)=
    **itype_trvg**
  - I
  - 2
  -
  - type of vegetation transpiration parameterization
    - 2: BATS scheme, Dickinson (1984)
    - 3: Extended BATS scheme with additional prognostic variable for integrated plant transpiration since sunrise; should be used only with an appropriate first guess for this variable coming from the DWD assimilation cycle
  -

* - (lnd_nml-itype_canopy)=
    **itype_canopy**
  - I
  - 1
  -
  - Type of canopy parameterization with respect to surface energy balance
    - 1: Surface energy balance equation solved at the ground surface, canopy energetically not represented
    - 2: Skin temperature formulation by Schulz and Vogel (2020), based on Viterbo and Beljaars (1995)
  -

* - (lnd_nml-cskinc)=
    **cskinc**
  - R
  - {math}`-1.0`
  - {math}`\mathrm{W m^{-2} K^{-1}}`
  - Skin conductivity
    For cskinc < 0, an external parameter field SKC is read and used
    For cskinc > 0, this globally constant value is used in the whole model domain
    Reasonable range: 10.0 -- 1000.0
  - itype_canopy = 2

* - (lnd_nml-tau_skin)=
    **tau_skin**
  - R
  - 3600.0
  - s
  - Relaxation time scale for the computation of the skin temperature
  - itype_canopy = 2

* - (lnd_nml-lterra_urb)=
    **lterra_urb**
  - L
  - .FALSE.
  -
  - If .TRUE., activate urban model TERRA_URB by Wouters et al.
    (2016, 2017)
    (see Schulz et al.
    2023)
  -

* - (lnd_nml-lurbalb)=
    **lurbalb**
  - L
  - .TRUE.
  -
  - If .TRUE., use urban albedo and emissivity (Wouters et al.
    2016)
  - lterra_urb = .TRUE.

* - (lnd_nml-itype_ahf)=
    **itype_ahf**
  - I
  - 2
  -
  - If {math}`\ge` 1, use urban anthropogenic heat flux (Wouters et al.
    2016)
    - 1: constant value given by the first entry in [tune_urbahf](nwp_tuning_nml-tune_urbahf)
    - 2: variable value depending on climatological T2M as specified in [tune_urbahf](nwp_tuning_nml-tune_urbahf)
    - 3: as option 2, but using a time-filtered predicted T2M rather than climatology
  - lterra_urb = .TRUE.

* - (lnd_nml-itype_kbmo)=
    **itype_kbmo**
  - I
  - 2
  -
  - Type of bluff-body thermal roughness length parameterisation
    - 1: Standard SAI-based turbtran (Raschendorfer 2001)
    - 2: Brutsaert-Kanda parameterisation for bluff-body elements (kB-1) (Kanda et al.
    2007)
    - 3: Zilitinkevich (1970)
  - lterra_urb = .TRUE.

* - (lnd_nml-itype_eisa)=
    **itype_eisa**
  - I
  - 3
  -
  - Type of evaporation from impervious surface area
    - 1: Evaporation like bare soil (see Schulz and Vogel 2020)
    - 2: No evaporation
    - 3: PDF-based puddle evaporation (Wouters et al.
    2015)
  - lterra_urb = .TRUE.

* - (lnd_nml-itype_heatcond)=
    **itype_heatcond**
  - I
  - 2
  -
  - type of soil thermal conductivity
    - 1: constant soil thermal conductivity
    - 2: moisture dependent soil thermal conductivity (see Schulz et al.
    2016)
    - 3: variant of option 2 with reduced near-surface thermal conductivity in the presence of plant cover
  -

* - (lnd_nml-itype_interception)=
    **itype_interception**
  - I
  - 1
  -
  - type of plant interception
    - 1: standard scheme, effectively switched off by tiny value cwimax_ml
  -

* - (lnd_nml-cwimax_ml)=
    **cwimax_ml**
  - R
  - {math}`1.e^{-6}`
  - m
  - scaling parameter for maximum interception storage (almost switched off);
    use {math}`5.e-4` to activate interception storage
  - itype_interception = 1

* - (lnd_nml-c_soil)=
    **c_soil**
  - R
  - 1.0
  -
  - surface area density of the (evaporative) soil surface
    allowed range: 0 -- 2
  - itype_evsl = 2,3,4

* - (lnd_nml-c_soil_urb)=
    **c_soil_urb**
  - R
  - 1.0
  -
  - surface area density of the (evaporative) soil surface, urban areas
    allowed range: 0 -- 2
  - itype_evsl = 2,3,4

* - (lnd_nml-cr_bsmin)=
    **cr_bsmin**
  - R
  - 110.0
  - s/m
  - minimum bare soil evaporation resistance (see Schulz and Vogel 2020)
    Note: c_soil and c_soil_urb are ignored in this case
  - itype_evsl = 5 or icpl_da_sfcevap = 4

* - (lnd_nml-rsmin_fac)=
    **rsmin_fac**
  - R
  - 1.0
  -
  - scaling factor for rsmin.
  -

* - (lnd_nml-itype_hydbound)=
    **itype_hydbound**
  - I
  - 1
  -
  - type of hydraulic lower boundary condition
    - 1: none
    - 3: ground water as lower boundary of soil column
  -

* - (lnd_nml-lstomata)=
    **lstomata**
  - L
  - .TRUE.
  -
  - If .TRUE., use map of minimum stomatal resistance
    If .FALSE., use constant value of {math}`150\, \mathrm{s/m}`.
  -

* - (lnd_nml-l2tls)=
    **l2tls**
  - L
  - .TRUE.
  -
  - If .TRUE., forecast with 2-Time-Level integration scheme (mandatory in ICON)
  -

* - (lnd_nml-lseaice)=
    **lseaice**
  - L
  - .TRUE.
  -
  - .TRUE. for use of sea-ice model
  -

* - (lnd_nml-lprog_albsi)=
    **lprog_albsi**
  - L
  - .FALSE.
  -
  - If .TRUE., sea-ice albedo is computed prognostically
  - lseaice=.TRUE.

* - (lnd_nml-lbottom_hflux)=
    **lbottom_hflux**
  - L
  - .FALSE.
  -
  - If .TRUE., use parameterization for bottom heat flux in seaice scheme
  - lseaice=.TRUE.

* - (lnd_nml-llake)=
    **llake**
  - L
  - .TRUE.
  -
  - .TRUE. for use of lake model
  -

* - (lnd_nml-sstice_mode)=
    **sstice_mode**
  - I
  - 1
  -
  - 1: SST and sea ice fraction are read from the analysis. The SST is kept constant whereas the sea ice fraction can be modified by the seaice model. (This mode also applices to coupled atmo/ocean simulations.)
    - 2: SST and sea ice fraction are read from the analysis. The SST is updated by climatological increments on a daily basis. The sea ice fraction can be modified by the seaice model.
    - 3: SST and sea ice fraction are updated daily, based on climatological monthly means
    - 4: SST and sea ice fraction are updated daily, based on actual monthly means
    - 5: SST and sea ice fraction are updated daily, based on actual daily means (**not yet implemented**)
    - 6: SST and sea ice fraction are updated with user-defined interval
  - [iforcing](run_nml-iforcing)=3

* - (lnd_nml-itype_oskin_warm)=
    **itype_oskin_warm**
  - I
  - 0
  -
  - SST skin formulation: warm layer component
    - 0: off
    - 1: activated
  -

* - (lnd_nml-itype_oskin_cold)=
    **itype_oskin_cold**
  - I
  - 0
  -
  - SST skin formulation: cold skin component
    - 0: off
    - 1: activated
  -

* - (lnd_nml-hice_min)=
    **hice_min**
  - R
  - 0.05
  - m
  - Minimum sea-ice thickness
  - lseaice=.TRUE.

* - (lnd_nml-hice_max)=
    **hice_max**
  - R
  - 3.0
  - m
  - Maximum sea-ice thickness (for coupled runs assure consistency with seaice_limit)
  - lseaice=.TRUE.

* - (lnd_nml-albsi_snow_max)=
    **albsi_snow_max**
  - R
  - 0.8
  -
  - Maximum albedo of snow over sea ice
  - lseaice=.TRUE.

* - (lnd_nml-albsi_snow_min)=
    **albsi_snow_min**
  - R
  - 0.5
  -
  - Minimum albedo of snow over sea ice
  - lseaice=.TRUE.

* - (lnd_nml-albsi_max)=
    **albsi_max**
  - R
  - 0.7
  -
  - Maximum albedo of sea ice
  - lseaice=.TRUE.

* - (lnd_nml-albsi_min)=
    **albsi_min**
  - R
  - 0.48
  -
  - Minimum albedo of sea ice
  - lseaice=.TRUE.

* - (lnd_nml-sst_file_interval)=
    **sst_file_interval**
  - C
  - 'P1M'
  -
  - ISO 8601 interval after which new SST/CI files are opened.
  - sstice_mode=6

* - (lnd_nml-sst_td_filename)=
    **sst_td_filename**
  - C
  -
  -
  - Filename of SST input files for time dependent SST. Default is `<path>SST_<year>_<month>_<gridfile>`. May contain the keyword `<path>` which will be substituted by model_base_dir
    In case sstice_mode=6, keywords for extpar files and `<day>`, `<hh>`, `<mm>`, `<ss>` are also available. NetCDF variable should be named SST.
  - sstice_mode=3,4,5,6

* - (lnd_nml-ci_td_filename)=
    **ci_td_filename**
  - C
  -
  -
  - Filename of sea ice fraction input files for time dependent sea ice fraction. Default is `<path>CI_<year>_<month>_<gridfile>`. May contain the keyword `<path>` which will be substituted by model_base_dir
    In case sstice_mode=6, keywords for extpar files and `<day>`, `<hh>`, `<mm>`, `<ss>` are also available. NetCDF variable should be named SIC.
  - sstice_mode=3,4,5,6

* - (lnd_nml-lcuda_graph_lnd)=
    **lcuda_graph_lnd**
  - L
  - .FALSE.
  -
  - Activate cuda graphs for the land scheme. Automatically set to .FALSE. if not compiled with the ICON_USE_CUDA_GRAPH cpp key.
  - ICON_USE_CUDA_GRAPH activated
```



Defined and used in: {{ '[src/namelists/mo_lnd_nwp_nml.f90]({}/src/namelists/mo_lnd_nwp_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_ls_forcing_nml)=
## ls_forcing_nml

Parameters for large-scale forcing, valid for torus geometry, is_plane_torus=.TRUE.

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ls_forcing_nml-is_ls_forcing)=
    **is_ls_forcing**
  - L
  - .TRUE.
  -
  - switch for enabling LS forcing
  -

* - (ls_forcing_nml-is_subsidence_moment)=
    **is_subsidence_moment**
  - L
  - .FALSE.
  -
  - switch for enabling LS vertical advection due to subsidence for momentum equations
  -

* - (ls_forcing_nml-is_subsidence_heat)=
    **is_subsidence_heat**
  - L
  - .FALSE.
  -
  - switch for enabling LS vertical advection due to subsidence for thermal equations
  -

* - (ls_forcing_nml-is_advection)=
    **is_advection**
  - L
  - .FALSE.
  -
  - switch for enabling LS horizontal advection
  -

* - (ls_forcing_nml-is_advection_uv)=
    **is_advection_uv**
  - L
  - .TRUE.
  -
  - switch for enabling LS horizontal advection for u and v
  - is_advection=.TRUE.

* - (ls_forcing_nml-is_advection_t)=
    **is_advection_t**
  - L
  - .TRUE.
  -
  - switch for enabling LS horizontal advection for temperature
  - is_advection=.TRUE.

* - (ls_forcing_nml-is_advection_q)=
    **is_advection_q**
  - L
  - .TRUE.
  -
  - switch for enabling LS horizontal advection for moisture
  - is_advection=.TRUE.

* - (ls_forcing_nml-is_nudging)=
    **is_nudging**
  - L
  - .FALSE.
  -
  - switch for enabling LS Newtonian relaxation (nudging)
  -

* - (ls_forcing_nml-is_nudging_uv)=
    **is_nudging_uv**
  - L
  - .TRUE.
  -
  - switch for enabling LS Newtonian relaxation (nudging) for horizontal winds only
  - is_nudging=.TRUE.

* - (ls_forcing_nml-is_nudging_tq)=
    **is_nudging_tq**
  - L
  - .TRUE.
  -
  - switch for enabling LS Newtonian relaxation (nudging) for temperature and specific humidity only
  - is_nudging=.TRUE.

* - (ls_forcing_nml-nudge_start_height_uv)=
    **nudge_start_height_uv**
  - R
  - 1000.0
  - m
  - height where nudging starts for winds
  - is_nudging=.TRUE.

* - (ls_forcing_nml-nudge_start_height_t)=
    **nudge_start_height_t**
  - R
  - 1000.0
  - m
  - height where nudging starts for temperature
  - is_nudging=.TRUE.

* - (ls_forcing_nml-nudge_start_height_q)=
    **nudge_start_height_q**
  - R
  - 1000.0
  - m
  - height where nudging starts for moisture
  - is_nudging=.TRUE.

* - (ls_forcing_nml-nudge_full_height_uv)=
    **nudge_full_height_uv**
  - R
  - 2000.0
  - m
  - height where nudging reaches full strength for winds
  - is_nudging=.TRUE.

* - (ls_forcing_nml-nudge_full_height_t)=
    **nudge_full_height_t**
  - R
  - 2000.0
  - m
  - height where nudging reaches full strength for temperature
  - is_nudging=.TRUE.

* - (ls_forcing_nml-nudge_full_height_q)=
    **nudge_full_height_q**
  - R
  - 2000.0
  - m
  - height where nudging reaches full strength for moisture
  - is_nudging=.TRUE.

* - (ls_forcing_nml-dt_relax_uv)=
    **dt_relax_uv**
  - R
  - 3600.0
  - s
  - relaxation time scale for the nudging for winds
  - is_nudging=.TRUE.

* - (ls_forcing_nml-dt_relax_t)=
    **dt_relax_t**
  - R
  - 3600.0
  - s
  - relaxation time scale for the nudging for temperature
  - is_nudging=.TRUE.

* - (ls_forcing_nml-dt_relax_q)=
    **dt_relax_q**
  - R
  - 3600.0
  - s
  - relaxation time scale for the nudging for moisture
  - is_nudging=.TRUE.

* - (ls_forcing_nml-is_geowind)=
    **is_geowind**
  - L
  - .FALSE.
  -
  - switch for enabling geostrophic wind
  -

* - (ls_forcing_nml-is_ls_coriolis)=
    **is_ls_coriolis**
  - L
  - .FALSE.
  -
  - Coriolis term for SCM or LES calculated based on domain average winds as part of large-scale forcing
  -

* - (ls_forcing_nml-is_rad_forcing)=
    **is_rad_forcing**
  - L
  - .FALSE.
  -
  - switch for enabling radiative forcing
  - inwp_rad=.FALSE.

* - (ls_forcing_nml-is_sim_rad)=
    **is_sim_rad**
  - L
  - .FALSE.
  -
  - switch for enabling a simplified radiation scheme
  - inwp_rad=.FALSE.

* - (ls_forcing_nml-is_theta)=
    **is_theta**
  - L
  - .FALSE.
  -
  - switch to indicate that the prescribed radiative forcing is for potential temperature
  - is_rad_forcing=.TRUE.
```



Defined and used in: {{ '[src/namelists/mo_ls_forcing_nml.f90]({}/src/namelists/mo_ls_forcing_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_master_nml)=
## master_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (master_nml-institute)=
    **institute**
  - C
  - ''
  -
  - Acronym of the institute for which the full institute name is printed in the log file. Options are DWD, MPIM, KIT, or CSCS. Otherwise the full names of MPIM and DWD are printed.
  -

* - (master_nml-lrestart)=
    **lrestart**
  - L
  - .FALSE.
  -
  - If .TRUE.: Current experiment is started from a restart.
  -

* - (master_nml-read_restart_namelists)=
    **read_restart_namelists**
  - L
  - .TRUE.
  -
  - If .TRUE.: Namelists are read from the restart file to override the default namelist settings, before reading new namelists from the run script. Otherwise the namelists stored in the restart file are ignored.
  -

* - (master_nml-lrestart_write_last)=
    **lrestart_write_last**
  - L
  - .FALSE.
  -
  - If .TRUE.: model run should create restart at experiment end. This is independent from the settings of the restart interval.
  -

* - (master_nml-model_base_dir)=
    **model_base_dir**
  - C
  - ''
  -
  - General path which may be used in file names of other name lists: If a file name contains the keyword `<path>`, then this `model_base_dir` will be substituted.
  -
```




(ref_buildrun_nml_master_model_nml)=
## master_model_nml

Repeated for each model.

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (master_model_nml-model_name)=
    **model_name**
  - C
  -
  -
  - Character string for naming this component.
  -

* - (master_model_nml-model_namelist_filename)=
    **model_namelist_filename**
  - C
  -
  -
  - File name containing the model namelists.
  -

* - (master_model_nml-model_type)=
    **model_type**
  - I
  - -1
  -
  - Identifies which component to run.
    - 1: atmosphere
    - 2: ocean
    - 3: radiation
    - 4: hamocc
    - 5: jsbach
    - 98: wave
    - 99: dummy_model
  -

* - (master_model_nml-model_do_restart)=
    **model_do_restart**
  - C
  - 'yes'/'no' based on lrestart
  -
  - Allows to overwrite the main restart switch lrestart for a single model component.
    Available options:
    `no` : no restart
    `yes`: classic restart
    `for_init`: initialize model from restart file (skips reading of restart attributes and ignores missing variables in restart file)
  -

* - (master_model_nml-model_min_rank)=
    **model_min_rank**
  - I
  - 0
  -
  - Start MPI rank for this model.
  -

* - (master_model_nml-model_max_rank)=
    **model_max_rank**
  - I
  - -1
  -
  - End MPI rank for this model.
  -

* - (master_model_nml-model_inc_rank)=
    **model_inc_rank**
  - I
  - 1
  -
  - Stride of MPI ranks.
  -

* - (master_model_nml-model_rank_group_size)=
    **model_rank_group_size**
  - I
  - 1
  -
  - ???
  -
```




(ref_buildrun_nml_master_time_control_nml)=
## master_time_control_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (master_time_control_nml-calendar)=
    **calendar**
  - C
  - "proleptic gregorian"
  -
  - Selects the calendar type to use:
    "proleptic gregorian" = proleptic Gregorian calendar
    "365 day year" = 365 day year without leap years
    "360 day year" = 360 day year with 30 day months
  -

* - (master_time_control_nml-experimentReferenceDate)=
    **experimentReferenceDate**
  - C
  - " "
  - ISO8601 formatted string
  - This specifies the reference date for the [calendar](master_time_control_nml-calendar) in use. It is an anchor date for cycling of events on the time line. If this namelist parameter is unspecified, then the reference date is set to the experiment start date.
  -

* - (master_time_control_nml-experimentStartDate)=
    **experimentStartDate**
  - C
  - " "
  - ISO8601 formatted string
  - This is the start date of an experiment, which remains valid for the whole experiment. The start date is also the reference date of the experiment, which is the anchor point for cycling events. In special cases the reference date might be reset. Reasons might be debugging purposes or spinning off experiments from an existing restart of an other experiment.
  -

* - (master_time_control_nml-experimentStopDate)=
    **experimentStopDate**
  - C
  - " "
  - ISO8601 formatted string
  - This is the date an experiment is finished.
  -

* - (master_time_control_nml-forecastLeadTime)=
    **forecastLeadTime**
  - C
  - " "
  - ISO8601 formatted string
  - Specifies the time span for a numerical weather forecast. It is used to set the experiment stop time with respect to the experiment start date.
  -

* - (master_time_control_nml-checkpointTimeIntVal)=
    **checkpointTimeIntVal**
  - C
  - " "
  - ISO8601 formatted string
  - Time interval for writing checkpoints.
  -

* - (master_time_control_nml-restartTimeIntVal)=
    **restartTimeIntVal**
  - C
  - " "
  - ISO8601 formatted string
  - Time interval for writing a restart file and interrupt the current running job.
  -
```






(ref_buildrun_nml_meteogram_output_nml)=
## meteogram_output_nml

This namelist is relevant if [output](run_nml-output)="nml".
Nearest neighbour 'interpolation' is used for all variables.

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (meteogram_output_nml-lmeteogram_enabled)=
    **lmeteogram_enabled**
  - L(n_dom)
  - .FALSE.
  -
  - Flag. True, if meteogram of output variables is desired.
  -

* - (meteogram_output_nml-zprefix)=
    **zprefix**
  - C(n_dom)
  - "METEO GRAM_"
  -
  - string with file name prefix for output file
  -

* - (meteogram_output_nml-ldistributed)=
    **ldistributed**
  - L(n_dom)
  - .TRUE.
  -
  - Flag. Separate files for each PE.
  -

* - (meteogram_output_nml-loutput_tiles)=
    **loutput_tiles**
  - L
  - .FALSE.
  -
  - Write tile-specific output for some selected surface/soil fields
  -

* - (meteogram_output_nml-n0_mtgrm)=
    **n0_mtgrm**
  - I(n_dom)
  - 0
  -
  - initial time step for meteogram output.
  -

* - (meteogram_output_nml-ninc_mtgrm)=
    **ninc_mtgrm**
  - I(n_dom)
  - 1
  -
  - output interval (in time steps)
  -

* - (meteogram_output_nml-stationlist_tot)=
    **stationlist_tot**
  -
  - 53.633,  9.983, 'Hamburg'
  -
  - list of meteogram stations (triples with lat, lon, name string)
  -

* - (meteogram_output_nml-silent_flush)=
    **silent_flush**
  - L(n_dom)
  - 1
  -
  - do not warn about flushing to disk if .TRUE.
  -

* - (meteogram_output_nml-max_time_stamps)=
    **max_time_stamps**
  - I(n_dom)
  - 1
  -
  - number of output time steps to record in memory before flushing to disk
  -

* - (meteogram_output_nml-var_list)=
    **var_list**
  - C(:)
  - " "
  -
  - Positive-list of variables (optional). Only variables contained in this list are included in the meteogram. If the default list is not changed by user input, then all available variables are added to the meteogram
  -
```



Defined and used in: {{ '[src/namelists/mo_mtgrm_nml.f90]({}/src/namelists/mo_mtgrm_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_mpiom_phy_nml)=
## mpiom_phy_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (mpiom_phy_nml-lmpiom_convection)=
    **lmpiom_convection**
  - LOGICAL
  -
  -
  -
  -

* - (mpiom_phy_nml-lmpiom_gentmcwill)=
    **lmpiom_gentmcwill**
  - LOGICAL
  -
  -
  -
  -

* - (mpiom_phy_nml-lmpiom_radiation)=
    **lmpiom_radiation**
  - LOGICAL
  -
  -
  -
  -
```


Defined and used in: {{ '[src/ocean/config/mo_ocean_nml.f90]({}/src/ocean/config/mo_ocean_nml.f90)'.format(base_url) }}

(ref_buildrun_nml_nonhydrostatic_nml)=
## nonhydrostatic_nml


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nonhydrostatic_nml-itime_scheme)=
    **itime_scheme**
  - I
  - 4
  -
  - Options for predictor-corrector time-stepping scheme:
    - 4: Contravariant vertical velocity is computed in the predictor step only, velocity tendencies are computed in the corrector step only (most efficient option)
    - 5: Contravariant vertical velocity is computed in both substeps (beneficial for numerical stability in very-high resolution setups with extremely steep slops, otherwise no significant impact)
    - 6: As 5, but velocity tendencies are also computed in both substeps (no apparent benefit, but more expensive)
  -

* - (nonhydrostatic_nml-rayleigh_type)=
    **rayleigh_type**
  - I
  - 2
  -
  - Type of Rayleigh damping
    - 1: CLASSICAL (requires velocity reference state!)
    - 2: Klemp (2008) type
  -

* - (nonhydrostatic_nml-rayleigh_coeff)=
    **rayleigh_coeff**
  - R(n_dom)
  - 0.05
  -
  - Rayleigh damping coefficient {math}`1/\tau_{0}` (Klemp, Dudhia, Hassiotis: MWR136, pp.3987-4004); higher values are recommended for R2B6 or finer resolution
  -

* - (nonhydrostatic_nml-damp_height)=
    **damp_height**
  - R(n_dom)
  - 45000
  - m
  - Height at which Rayleigh damping of vertical wind starts (needs to be adjusted to model top height; the damping layer should have a depth of at least 20 km when the model top is above the stratopause)
  -

* - (nonhydrostatic_nml-htop_moist_proc)=
    **htop_moist_proc**
  - R
  - 22500.0
  - m
  - Height above which moist physics and advection of cloud and precipitation variables are turned off
  -

* - (nonhydrostatic_nml-hbot_qvsubstep)=
    **hbot_qvsubstep**
  - R
  - 22500.0
  - m
  - Height above which QV is advected with substepping scheme
  - ihadv_tracer=22, 32, 42 or 52

* - (nonhydrostatic_nml-htop_aero_proc)=
    **htop_aero_proc**
  - R
  - 22500.0
  - m
  - Height above which physical processes and advection of ART aerosol tracer variables are turned off; the default value is set to the same value as htop_moist_proc. This value is taken for all ART aerosol tracers, but not chemical tracers for which physical processes and advection are computed in all model levels per default; it may be overwritten for specific ART tracers (also chemical tracers) by the tag 'htop_proc' in the XML file when defining the individual ART tracers.
  - ART aerosol tracers (with an index {math}`\ge` iqt)

* - (nonhydrostatic_nml-vwind_offctr)=
    **vwind_offctr**
  - R
  - 0.15
  -
  - Off-centering in vertical wind solver. Higher values may be needed for R2B5 or coarser grids when the model top is above 50 km. Negative values are not allowed
  -

* - (nonhydrostatic_nml-rhotheta_offctr)=
    **rhotheta_offctr**
  - R
  - -0.1
  -
  - Off-centering of density and potential temperature at interface level (may be set to 0.0 for R2B6 or finer grids; positive values are not recommended)
  -

* - (nonhydrostatic_nml-veladv_offctr)=
    **veladv_offctr**
  - R
  - 0.25
  -
  - Off-centering of velocity advection in corrector step. Negative values are not recommended
  -

* - (nonhydrostatic_nml-ivctype)=
    **ivctype**
  - I
  - 2
  -
  - Type of vertical coordinate:
    - 1: Gal-Chen hybrid
    - 2: SLEVE (uses [sleve_nml](ref_buildrun_nml_sleve_nml))
  -

* - (nonhydrostatic_nml-ndyn_substeps)=
    **ndyn_substeps**
  - I
  - 5
  -
  - number of dynamics substeps per fast-physics / transport step
  -

* - (nonhydrostatic_nml-vcfl_threshold)=
    **vcfl_threshold**
  - R
  - 1.05
  -
  - Threshold for vertical advection CFL number at which the adaptive time step reduction (increase of ndyn_substeps w.r.t. the fixed fast-physics time step) is triggered.
  -

* - (nonhydrostatic_nml-nlev_hcfl)=
    **nlev_hcfl**
  - I(max_dom)
  - 0
  -
  - Number of model levels (counted from the top) for which the horizontal CFL number is evaluated in addition and used for an adaptive dynamics time step reduction. In practice, doing this for the upper 10-15 levels is sufficient with a model top of 75 km.
  -

* - (nonhydrostatic_nml-cfl_monitoring_freq)=
    **cfl_monitoring_freq**
  - I
  - 5
  -
  - Monitoring frequency for CFL number (in units of fast-physics time steps of domain 1)
  -

* - (nonhydrostatic_nml-lextra_diffu)=
    **lextra_diffu**
  - L
  - .TRUE.
  -
  - .TRUE.: Apply additional momentum diffusion at grid points close to the stability limit for vertical advection (becomes effective extremely rarely in practice; this is mostly an emergency fix for pathological cases with very large orographic gravity waves)
  -

* - (nonhydrostatic_nml-divdamp_fac)=
    **divdamp_fac**
  - R
  - 0.0025
  -
  - Scaling factor for divergence damping at height [divdamp_z](nonhydrostatic_nml-divdamp_z) and below. [divdamp_fac](nonhydrostatic_nml-divdamp_fac) {math}`\geq 0`.
  -

* - (nonhydrostatic_nml-divdamp_fac2)=
    **divdamp_fac2**
  - R
  - 0.004
  -
  - Scaling factor for divergence damping at height [divdamp_z2](nonhydrostatic_nml-divdamp_z2). [divdamp_fac2](nonhydrostatic_nml-divdamp_fac2) {math}`\geq 0`. Between [divdamp_z](nonhydrostatic_nml-divdamp_z) and [divdamp_z2](nonhydrostatic_nml-divdamp_z2) the scaling factor changes linearly from [divdamp_fac](nonhydrostatic_nml-divdamp_fac) to [divdamp_fac2](nonhydrostatic_nml-divdamp_fac2).
  -

* - (nonhydrostatic_nml-divdamp_fac3)=
    **divdamp_fac3**
  - R
  - 0.004
  -
  - Scaling factor for divergence damping at height [divdamp_z3](nonhydrostatic_nml-divdamp_z3). [divdamp_fac3](nonhydrostatic_nml-divdamp_fac3) {math}`\geq 0`. The three points ([divdamp_z2](nonhydrostatic_nml-divdamp_z2), [divdamp_fac2](nonhydrostatic_nml-divdamp_fac2)), ([divdamp_z3](nonhydrostatic_nml-divdamp_z3), [divdamp_fac3](nonhydrostatic_nml-divdamp_fac3)), and [divdamp_z4](nonhydrostatic_nml-divdamp_z4), *divdamp_fac4*) determine the quadratic function for the scaling factor between [divdamp_z2](nonhydrostatic_nml-divdamp_z2) and [divdamp_z4](nonhydrostatic_nml-divdamp_z4).
  -

* - (nonhydrostatic_nml-divdamp_fac4)=
    **divdamp_fac4**
  - R
  - 0.004
  -
  - Scaling factor for divergence damping at height [divdamp_z4](nonhydrostatic_nml-divdamp_z4) and higher. [divdamp_fac4](nonhydrostatic_nml-divdamp_fac4) {math}`\geq 0`.
  -

* - (nonhydrostatic_nml-divdamp_z)=
    **divdamp_z**
  - R
  - 32500.0
  - m
  - Height up to which [divdamp_fac](nonhydrostatic_nml-divdamp_fac) is used, and where the linear profile up to height [divdamp_z2](nonhydrostatic_nml-divdamp_z2) starts.
  -

* - (nonhydrostatic_nml-divdamp_z2)=
    **divdamp_z2**
  - R
  - 40000.0
  - m
  - Height with scaling factor [divdamp_fac2](nonhydrostatic_nml-divdamp_fac2) where the linear profile starting at [divdamp_z](nonhydrostatic_nml-divdamp_z) ends, and where the quadratic profile up to [divdamp_z4](nonhydrostatic_nml-divdamp_z4) starts. [divdamp_z](nonhydrostatic_nml-divdamp_z) < [divdamp_z2](nonhydrostatic_nml-divdamp_z2) < [divdamp_z4](nonhydrostatic_nml-divdamp_z4).
  -

* - (nonhydrostatic_nml-divdamp_z3)=
    **divdamp_z3**
  - R
  - 60000.0
  - m
  - Height with scaling factor [divdamp_fac3](nonhydrostatic_nml-divdamp_fac3). Needed to determine the quadratic function between [divdamp_z2](nonhydrostatic_nml-divdamp_z2) and [divdamp_z4](nonhydrostatic_nml-divdamp_z4). [divdamp_z3](nonhydrostatic_nml-divdamp_z3) {math}`\neq` [divdamp_z2](nonhydrostatic_nml-divdamp_z2) {math}`\land` [divdamp_z3](nonhydrostatic_nml-divdamp_z3) {math}`\neq` [divdamp_z4](nonhydrostatic_nml-divdamp_z4).
  -

* - (nonhydrostatic_nml-divdamp_z4)=
    **divdamp_z4**
  - R
  - 80000.0
  - m
  - Height from which [divdamp_fac4](nonhydrostatic_nml-divdamp_fac4) is used. [divdamp_z4](nonhydrostatic_nml-divdamp_z4) > [divdamp_z2](nonhydrostatic_nml-divdamp_z2).
  -

* - (nonhydrostatic_nml-divdamp_order)=
    **divdamp_order**
  - I
  - 4
  -
  - Order of divergence damping:
    - 2: second-order divergence damping
    - 4: fourth-order divergence damping
    - 24: combined second-order and fourth-order divergence damping and enhanced vertical wind off-centering during the initial spinup phase (does not allow checkpointing/restarting earlier than 2.5 hours of integration)
  -

* - (nonhydrostatic_nml-divdamp_type)=
    **divdamp_type**
  - I
  - 3
  -
  - Type of divergence damping:
    - 2: divergence damping acting on 2D divergence
    - 3: divergence damping acting on 3D divergence
    - 32: combination of 3D div. damping in the troposphere with transition to 2D div. damping in the stratosphere
  -

* - (nonhydrostatic_nml-divdamp_trans_start)=
    **divdamp_trans_start**
  - R
  - 12500.0
  -
  - Lower bound of transition zone between 2D and 3D divergence damping
  - *divdamp_type* = 32

* - (nonhydrostatic_nml-divdamp_trans_end)=
    **divdamp_trans_end**
  - R
  - 17500.0
  -
  - Upper bound of transition zone between 2D and 3D divergence damping
  - *divdamp_type* = 32

* - (nonhydrostatic_nml-iadv_rhotheta)=
    **iadv_rhotheta**
  - I
  - 2
  -
  - Advection method for rho and rhotheta:
    - 1: simple second-order upwind-biased scheme
    - 2: 2nd order Miura horizontal
  -

* - (nonhydrostatic_nml-igradp_method)=
    **igradp_method**
  - I
  - 3
  -
  - Discretization of horizontal pressure gradient:
    - 1: conventional discretization with metric correction term
    - 2: Taylor-expansion-based reconstruction of pressure (advantageous at very high resolution)
    - 3: Similar discretization as option 2, but uses hydrostatic approximation for downward extrapolation over steep slopes
    - 4: Cubic/quadratic polynomial interpolation for pressure reconstruction
    - 5: Same as 4, but hydrostatic approximation for downward extrapolation over steep slopes
  -

* - (nonhydrostatic_nml-l_zdiffu_t)=
    **l_zdiffu_t**
  - L
  - .TRUE.
  -
  - .TRUE.: Compute Smagorinsky temperature diffusion truly horizontally over steep slopes
  - hdiff_order=5 .AND. lhdiff_temp = .TRUE.

* - (nonhydrostatic_nml-thslp_zdiffu)=
    **thslp_zdiffu**
  - R
  - 0.025
  -
  - Slope threshold above which truly horizontal temperature diffusion is activated
  - hdiff_order=5 .AND. lhdiff_temp=.TRUE. .AND. l_zdiffu_t=.TRUE.

* - (nonhydrostatic_nml-thhgtd_zdiffu)=
    **thhgtd_zdiffu**
  - R
  - 200
  - m
  - Threshold of height difference between neighboring grid points above which truly horizontal temperature diffusion is activated (alternative criterion to thslp_zdiffu)
  - hdiff_order=5 .AND. lhdiff_temp=.TRUE. .AND. l_zdiffu_t=.TRUE.

* - (nonhydrostatic_nml-exner_expol)=
    **exner_expol**
  - R
  - 1./3.
  -
  - Temporal extrapolation (fraction of dt) of Exner function for computation of horizontal pressure gradient. This damps horizontally propagating sound waves. For R2B5 or coarser grids, values between 1/2 and 2/3 are recommended. Model will be numerically unstable for negative values.
  -
```



Defined and used in: {{ '[src/namelists/mo_nonhydrostatic_nml.f90]({}/src/namelists/mo_nonhydrostatic_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_nudging_nml)=
## nudging_nml

Parameters for the upper boundary nudging in the limited-area mode
(l_limited_area = .TRUE.) or global nudging.
For the lateral boundary nudging, please see [interpol_nml](ref_buildrun_nml_interpol_nml) and [limarea_nml](ref_buildrun_nml_limarea_nml).
The characteristics of the driving data for the nudging can be specified
in [limarea_nml](ref_buildrun_nml_limarea_nml).


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nudging_nml-nudge_type)=
    **nudge_type**
  - I(n_dom)
  - 0
  -
  - Nudging type:
    * 0: none
    * 1: upper boundary nudging
    * 2: global nudging
    Please note:
    * nudge_type = 1 requires l_limited_area = .TRUE.
    * nudge_type = 1 is also applicable to nested domains. Nudging is performed against the same forcing data set for all domains. If nudging is enabled for one or more nested domains, it needs to be enabled for the base domain, as well.
    * nudge_type = 2 (global nudging) is applied in primary domain only
    * for global nudging the following settings in [limarea_nml](ref_buildrun_nml_limarea_nml) are mandatory:
      - itype_latbc = 1 (time-dependent driving data)
      - dtime_latbc = {math}`\ldots`
      - latbc_path = "{math}`\ldots`"
      - latbc_boundary_grid = " " (no boundary grid: driving data have to be available on entire grid)
      - latbc_varnames_map_file = "{math}`\ldots`" (e.g., run/dict.latbc), if num_prefetch_proc = 1 (asynchronous read-in of driving data)
    * defaults and (additional) scopes for global nudging are marked by {math}`(\cdot)_\text{glbndg}`, if a parameter applies to both upper boundary and global nudging
  - [iforcing](run_nml-iforcing) = 3 (NWP)
    ivctype = 2 (SLEVE)

* - (nudging_nml-max_nudge_coeff_vn)=
    **max_nudge_coeff_vn**
  - R
  - 0.04 {math}`(0.016)_\text{glbndg}`
  - 1
  - Max. nudging coefficient for the horizontal wind (i.e.
    the edge-normal wind component {math}`v_n`). Given the wind update due to the nudging term on the rhs:
    {math}`v_n(t) = v_n^\star(t) + \text{nudge\_coeff\_vn}(z) * \text{ndyn\_substeps} * [\overline{v_n}(t) - v_n^\star(t)]`,
    where {math}`t` and {math}`z` denote time and height, respectively, {math}`\overline{v_n}(t)` is the target wind to nudge to, and {math}`v_n^\star` is the value before the nudging, the vertical profile of the coefficient for upper boundary nudging reads:
    {math}`\text{nudge\_coeff\_vn}(z) = \text{max\_nudge\_coeff\_vn} * [(z - \text{nudge\_start\_height})/(\text{top\_height} - \text{nudge\_start\_height})]^2`,
    for {math}`\text{nudge\_start\_height} \leq z \leq \text{top\_height}` (see nudge_start_height below), and is zero elsewhere. The range of validity is {math}`\text{max\_nudge\_coeff\_vn} \in [0, \sim 0.2]`, where the lower boundary is mandatory.
    **Please note** that the user value is internally multiplied by 5.
  - nudge_type > 0
    (nudge_var = "all" or "...,vn,..."){math}`{}_\text{glbndg}`

* - (nudging_nml-max_nudge_coeff_thermdyn)=
    **max_nudge_coeff_thermdyn**
  - R
  - 0.075 {math}`(0.03)_\text{glbndg}`
  - 1
  - Max. nudging coefficient for the thermodynamic variables selected by nudge_hydro_pres in case of upper boundary nudging and by thermdyn_type in case of global nudging. The range of validity is {math}`\text{max\_nudge\_coeff\_thermdyn} \in [0, \sim 0.2]`, where the lower boundary is mandatory.
    **Please note** that the user value is internally multiplied by 5.
  - nudge_type > 0
    (nudge_var = "all" or "...,thermdyn,..."){math}`{}_\text{glbndg}`

* - (nudging_nml-max_nudge_coeff_qv)=
    **max_nudge_coeff_qv**
  - R
  - 0.008
  - 1
  - Max. nudging coefficient for water vapor. The range of validity is {math}`\text{max\_nudge\_coeff\_qv} \in [0, \sim 0.2]`, where the lower boundary is mandatory. (For global nudging only.)
    **Please note** that the user value is internally multiplied by 5.
  - nudge_type = 2
    nudge_var = "{all}" { or} "...,qv,..."

* - (nudging_nml-nudge_start_height)=
    **nudge_start_height**
  - R
  - 12000 {math}`(2000)_\text{glbndg}`
  - m
  - Nudging is applied for:
    {math}`\text{nudge\_start\_height} \leq z \leq \text{top\_height}` in case of upper boundary nudging and for:
    {math}`\text{nudge\_start\_height} \leq z \leq \text{nudge\_end\_height}` in case of global nudging, where {math}`z` denotes the nominal height of the grid layer center, and top_height is the height of the model top (see [sleve_nml](ref_buildrun_nml_sleve_nml)). For upper boundary nudging the range of validity is {math}`\text{nudge\_start\_height} \in [0, \text{top\_height}]`, where both boundaries are mandatory. For global nudging a nudge_start_height in the range {math}`[0, \text{top\_height}]` has to satisfy nudge_start_height < nudge_end_height. Values outside {math}`[0, \text{top\_height}]` will be interpreted as nudge_start_height = 0.
  - nudge_type > 0

* - (nudging_nml-nudge_end_height)=
    **nudge_end_height**
  - R
  - 40000
  - m
  - Nudging is applied for:
    {math}`\text{nudge\_start\_height} \leq z \leq \text{nudge\_end\_height}`, where {math}`z` denotes the nominal height of the grid layer center. A nudge_end_height in the range {math}`[0, \text{top\_height}]` has to satisfy nudge_start_height < nudge_end_height. Values outside {math}`[0, \text{top\_height}]` will be interpreted as nudge_start_height = top_height. (For global nudging only.)
  - nudge_type = 2

* - (nudging_nml-nudge_profile)=
    **nudge_profile**
  - I
  - 4
  -
  - Vertical profile of the nudging coefficient (nudging strength) between nudge_start_height and nudge_end_height:
    * 1: squared scaled vertical distance from nudge_start_height (this is the profile used for upper boundary nudging)
    * 2: constant profile
    * 3: hyperbolic tangent profile
    * 4: trapezoidal profile
    The profile values range from 0 to 1. A multiplication with max_nudge_coeff_vn/thermdyn/qv and ndyn_substeps yields the final value of the nudging coefficient. (For global nudging only.)
  - nudge_type = 2

* - (nudging_nml-nudge_scale_height)=
    **nudge_scale_height**
  - R
  - 3000
  - m
  - Scale height of nudging profile. (For global nudging only.)
  - nudge_type = 2
    nudge_profile = 3 or 4

* - (nudging_nml-nudge_var)=
    **nudge_var**
  - C
  - "all"
  -
  - Select the variables that shall be nudged:
    * "vn": horizontal wind
    * "thermdyn": thermodynamic variables
    * "qv": water vapor
    * comma-separated list: e.g., "vn,thermdyn"
    * "all": all available variables (i.e. equivalent to "vn,thermdyn,qv")
    Please note that the nudging of water vapor requires ltransport = .TRUE. (For global nudging only.)
  - nudge_type = 2

* - (nudging_nml-thermdyn_type)=
    **thermdyn_type**
  - I
  - 1
  -
  - Set of variables used to compute the thermodynamic nudging increments:
    * 1: hydrostatic set (pressure and temperature)
    * 2: non-hydrostatic set (density and virtual potential temperature)
  - nudge_type = 2
    nudge_var = "all" or "...,thermdyn,..."
```



Defined and used in: {{ '[src/namelists/mo_nudging_nml.f90]({}/src/namelists/mo_nudging_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_nwp_phy_nml)=
## nwp_phy_nml

The switches for the physics schemes and the time steps can be set for each model domain individually.
If only one value is specified, it is copied to all child domains, implying that the same set
of parameterizations and time steps is used in all domains. If the number of values given
in the namelist is larger than 1 but less than the number of model domains, then the settings
from the highest domain ID are used for the remaining model domains.

If the time steps are not an integer multiple of the advective time step (dtime), then the time step of the
respective physics parameterization is automatically rounded to the next higher integer multiple
of the advective time step. If the radiation time step is not an integer multiple of the cloud-cover
time step it is automatically rounded to the next higher integer multiple of the cloud cover time step.


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nwp_phy_nml-inwp_gscp)=
    **inwp_gscp**
  - I (max_ dom)
  - 1
  -
  - cloud microphysics and precipitation
    - -1: Externally specified (e.g. ComIn)
    - 0: none
    - 1: cloud ice scheme (based on COSMO-EU microphysics, 2-cat ice: cloud ice, snow)
    - 2: graupel scheme (based on COSMO-DE microphysics, 3-cat ice: cloud ice, snow, graupel)
    - 3: Two-moment cloud ice scheme
    - 4,5,6,7: Two-moment microphysics of ICON-RUC, further configuration possible in namelist [twomom_mcrph_nml](ref_buildrun_nml_twomom_mcrph_nml)
    - 8: Spectral Bin Microphysics (SBM) by A. Khain
    - 9: Kessler scheme
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-qi0)=
    **qi0**
  - R
  - 0.0
  - kg/kg
  - cloud ice threshold for autoconversion
  - inwp_gscp=1

* - (nwp_phy_nml-qc0)=
    **qc0**
  - R
  - 0.0
  - kg/kg
  - cloud water threshold for autoconversion
  - inwp_gscp=1

* - (nwp_phy_nml-mu_rain)=
    **mu_rain**
  - R
  - 0.0
  -
  - shape parameter in gamma distribution for rain
  - inwp_gscp > 0

* - (nwp_phy_nml-rain_n0_factor)=
    **rain_n0_factor**
  - R
  - 1.0
  -
  - tuning factor for intercept parameter of raindrop size distribution
  - inwp_gscp > 0

* - (nwp_phy_nml-lvariable_rain_n0)=
    **lvariable_rain_n0**
  - L
  - .FALSE.
  -
  - variable intercept parameter of raindrop size distribution: the multiplicative factor rain_n0_factor is used for drizzle (small {math}`q_r`) while the default value is approached for heavy rain (large {math}`q_r`)
  - inwp_gscp=2

* - (nwp_phy_nml-mu_snow)=
    **mu_snow**
  - R
  - 0.0
  -
  - shape parameter in gamma distribution for snow
  - inwp_gscp > 0

* - (nwp_phy_nml-icpl_aero_gscp)=
    **icpl_aero_gscp**
  - I
  - 0
  -
  - 0: Constant CDNC - no aerosol-cloud interactions (default)
    - 1: Simple coupling between Segal and Khain (2006) activation scheme and Tegen aerosol climatology; requires irad_aero=6, 7 or 9
    - 2: Full coupling between CAMS climatology or forecasted aerosols with the Segal and Khain (2006) activation scheme; requires irad_aero= 7 or 8
    - 3: Use CDNC climatology from external parameter file. External parameter files containing cloud-droplet number climatology can be generated with Zonda
    More advanced options are in preparation
  - currently only for inwp_gscp = 1

* - (nwp_phy_nml-scale_cdnc_mode)=
    **scale_cdnc_mode**
  - I
  - 0
  -
  - Controls scaling of external climatological cloud droplet number concentration (CDNC) data. Only effective if external CDNC is used.
    - 0: No additional scaling applied, only scaling from tuning namelist is used.
    - 1: Apply time-varying, spatially resolved scaling factors derived from the simple plume scheme during the model run. The scaling field varies daily in time and 2D in space, and is applied according to the experiment’s actual (simulation) year.
    Note: For years 2000–2020, it is usually preferable to disable scaling (mode 0), since MODIS data are available for those years directly.
    - 2: Apply constant scaling of external CDNC to represent pre-industrial conditions. All years are scaled to 1850 values, using spatially varying scaling factors.
  - Climate projections using
    icpl_aero_gscp = 3
    irad_aero = 18,19
    Picontrol experiments using
    icpl_aero_gscp = 3
    irad_aero = 12

* - (nwp_phy_nml-icpl_aero_ice)=
    **icpl_aero_ice**
  - I
  - 0
  -
  - 0: ice crystal concentration defined by temperature using Cooper (1987) formula
    - 1: ice nucleating particles concentration defined by temperature and dust concentration using DeMott (2015) formula
    - 2: ice nucleating particles parameterization of Ullrich et al. (2017)
    - 3: prognostic mineral dust used as ice nucleating particles (dust number only)
    - 4: prognostic mineral dust used as ice nucleating particles (dust number and size)
  - 0: inwp_gscp= 1, 2
    - 1: inwp_gscp=1, 2 and irad_aero=6, 7 or 8
    - 2: inwp_gscp=3
    - 3: inwp_gscp=3 and ART
    - 4: inwp_gscp=3 and ART

* - (nwp_phy_nml-inwp_convection)=
    **inwp_convection**
  - I (max_ dom)
  - 1
  -
  - convection
    - 0: none
    - 1: Tiedtke/Bechtold convection
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-lshallowconv_only)=
    **lshallowconv_only**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: use shallow convection only
  - inwp_convection = 1; cannot be combined with lgrayzone_deepconv

* - (nwp_phy_nml-lgrayzone_deepconv)=
    **lgrayzone_deepconv**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: activates shallow and deep convection but not mid-level convection, together with some tuning measures targeted at grayzone (convection-permitting) model resolutions
  - inwp_convection = 1; cannot be combined with lshallowconv_only

* - (nwp_phy_nml-itype_parcel_ascent)=
    **itype_parcel_ascent**
  - I (max_ dom)
  - 1
  -
  - 1: ICON-NWP, 2: IFS/CY41r2
  - inwp_convection = 1

* - (nwp_phy_nml-ldetrain_conv_prec)=
    **ldetrain_conv_prec**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: Activate detrainment of convective rain and snow
  - inwp_convection = 1

* - (nwp_phy_nml-lconv_cdnc_interp)=
    **lconv_cdnc_interp**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: cloud_num based interpolation for rmfdeps and rhebc to replace land/sea mask. Should be combined with icpl_aero_conv=1.
  - inwp_convection = 1, icpl_aero_gscp=3

* - (nwp_phy_nml-icapdcycl)=
    **icapdcycl**
  - I
  - 0
  -
  - Type of CAPE correction to improve diurnal cycle for convection:
    - 0: none (IFS default prior to autumn 2013)
    - 1: intermediate testing option
    - 2: correctoins over land and water now operational at ECMWF
    - 3: correction over land as in 2 restricted to the tropics, no correction over water (this choice optimizes the NWP skill scores)
  - inwp_convection = 1

* - (nwp_phy_nml-lstoch_expl)=
    **lstoch_expl**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: activate explicit stochastic shallow convection scheme
    EXPERIMENTAL! will not produce clean restart
    to be used in conjunction with lrestune_off=.T. and lmflimiter_off=.T.
  - inwp_convection = 1

* - (nwp_phy_nml-lstoch_sde)=
    **lstoch_sde**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: activate stochastic differential equation (SDE) shallow convection scheme
    to be used in conjunction with lrestune_off=.T. and lmflimiter_off=.T.
  - inwp_convection = 1

* - (nwp_phy_nml-lstoch_deep)=
    **lstoch_deep**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: activate stochastic differential equation (SDE) deep convection scheme
  - inwp_convection = 1

* - (nwp_phy_nml-lrestune_off)=
    **lrestune_off**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: switches off resolution-dependent tuning of shallow convection parameters
  - inwp_convection = 1

* - (nwp_phy_nml-lmflimiter_off)=
    **lmflimiter_off**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: disables mass flux limiter by setting it to high values that are rarely reached by shallow convection
  - inwp_convection = 1

* - (nwp_phy_nml-lvvcouple)=
    **lvvcouple**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: use vertical velocity at 650hPa as criterion to couple shallow convection
    with resolved deep convection
  - inwp_convection = 1

* - (nwp_phy_nml-lvv_shallow_deep)=
    **lvv_shallow_deep**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: use vertical velocity at 650hPa to distinguish between shallow and
    deep convection within convection routines (instead of cloud depth)
  - inwp_convection = 1

* - (nwp_phy_nml-lstoch_spinup)=
    **lstoch_spinup**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: spin up cloud ensemble to equilibrium in stochastic shallow convection schemes,
    only takes effect when lstoch_expl=T or lstoch_sde=T
  - inwp_convection = 1

* - (nwp_phy_nml-nclds)=
    **nclds**
  - I (max_ dom)
  - 5000
  -
  - maximum possible number of shallow clouds per grid box in explicit stochastic cloud ensemble.
    only takes effect when lstoch_expl=T
  - inwp_convection = 1

* - (nwp_phy_nml-icpl_aero_conv)=
    **icpl_aero_conv**
  - I
  - 0
  -
  - 0: off
    - 1: simple coupling between autoconversion (in convection scheme) and aerosol climatology; requires irad_aero=6, 7, 12, 13, 14, 15, 18, 19
  - when irad_aero=12, 13, 14, 15, 18, 19 (Kinne aerosol), icpl_aero_conv=1 should be used with iclp_aero_gscp=3 (external cloud droplet number)

* - (nwp_phy_nml-i2daero_anthro)=
    **i2daero_anthro**
  - I
  - 0
  -
  - 0: off
    - 1: simple prognostic aerosol scheme for anthropogenic species
  - irad_aero=6

* - (nwp_phy_nml-i2daero_fire)=
    **i2daero_fire**
  - I
  - 0
  -
  - 0: off
    - 1: add wild fire emissions from file to black carbon, organic carbon and sulfate.
    Requires i2daero_anthro=1 and fire2d_filename. 2: add wild fire emissions from file to black carbon, organic carbon and sulfate from seasonal climatology. Requires the fields bcfire_clim, ocfire_clim and so2fire_clim in the external data.
  - irad_aero=6

* - (nwp_phy_nml-i2daero_dust)=
    **i2daero_dust**
  - I
  - 0
  -
  - 0: off
    - 1: simple prognostic aerosol scheme for mineral dust
  - irad_aero=6

* - (nwp_phy_nml-i2daero_seas)=
    **i2daero_seas**
  - I
  - 0
  -
  - 0: off
    - 1: simple prognostic aerosol scheme for sea salt
  - irad_aero=6

* - (nwp_phy_nml-icpl_gwd_prec)=
    **icpl_gwd_prec**
  - I
  - 0
  -
  - 0: launched momentum flux of GWD is proportional to gfluxlaun.
    - 1: launched momentum flux of GWD depends linearly from total precipitation
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-icpl_o3_tp)=
    **icpl_o3_tp**
  - I
  - 1
  -
  - 0: off
    - 1: simple coupling between the ozone mixing ratio and the thermal tropopause, restricted to the extratropics
    - 2: improved variant of 1 that avoids excessive additional ozone for low tropopauses
  - irad_o3 = 7 or 9

* - (nwp_phy_nml-itype_dissip_heat)=
    **itype_dissip_heat**
  - I
  - 1
  -
  - Options for calculating dissipative heating in NWP interface
    - 0: off; setting ltmpcor = .TRUE. forces reset to 0 because dissipative heating is then calculated in turbdiff
    - 1: dissipative heating from SSO, GWD and Rayleigh friction
    - 2: as 1, plus dissipative heating from momentum dissipation by turbulence
  - ltmpcor = .FALSE.

* - (nwp_phy_nml-inwp_cldcover)=
    **inwp_cldcover**
  - I (max_ dom)
  - 1
  -
  - cloud cover scheme for radiation
    - 0: no clouds (only QV)
    - 1: diagnostic cloud cover (by Martin Koehler)
    - 2: prognostic total water variance (not yet started)
    - 3: clouds from COSMO SGS cloud scheme
    - 4: clouds as in turbulence (turbdiff)
    - 5: grid scale clouds
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-lsgs_cond)=
    **lsgs_cond**
  - L (max_ dom)
  - .TRUE.
  -
  - Apply subgrid-scale condensational heating related to the non-convective part of diagnosed cloud water
  - inwp_cldcover = 1

* - (nwp_phy_nml-itype_icecloud_diag)=
    **itype_icecloud_diag**
  - I
  - 1
  -
  - 1: standard
    - 2: option for two-moment microphysics
  - inwp_cldcover = 1,
    option 2 requires inwp_gscp=3,4,5,6,7,8

* - (nwp_phy_nml-lsbm_coupled)=
    **lsbm_coupled**
  - L (max_ dom)
  - .TRUE.
  -
  - TRUE: use SBM feedback, FALSE: use 2M for feedback and run uncoupled SBM
  - inwp_gscp = 8

* - (nwp_phy_nml-lmicrophysicsFirst)=
    **lmicrophysicsFirst**
  - L (max_ dom)
  - .FALSE.
  -
  - TRUE: run microphysics before turbdiff, FALSE: after turbdiff
  -

* - (nwp_phy_nml-inwp_radiation)=
    **inwp_radiation**
  - I (max_ dom)
  - 1
  -
  - radiation
    - 0: none
    - 1: RRTM radiation
    - 2: (removed)
    - 3: (removed)
    - 4: ecRad radiation
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-inwp_satad)=
    **inwp_satad**
  - I
  - 1
  -
  - saturation adjustment
    - 0: none
    - 1: saturation adjustment at constant density
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-inwp_turb)=
    **inwp_turb**
  - I (max_ dom)
  - 1
  -
  - vertical diffusion and transfer
    - 0: none
    - 1: COSMO diffusion and transfer
    - 2: GME turbulence scheme
    - 5: Classical Smagorinsky diffusion
    - 6: VDIFF turbulence scheme (requires inwp_surface = 2)
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-inwp_sso)=
    **inwp_sso**
  - I (max_ dom)
  - 1
  -
  - subgrid scale orographic drag
    - 0: none
    - 1: Lott and Miller scheme (COSMO)
  - [iforcing](run_nml-iforcing) = inwp
    inwp_turb > 0

* - (nwp_phy_nml-inwp_gwd)=
    **inwp_gwd**
  - I (max_ dom)
  - 1
  -
  - non-orographic gravity wave drag
    - 0: none
    - 1: Orr-Ern-Bechtold-scheme (IFS)
  - [iforcing](run_nml-iforcing) = inwp
    inwp_turb > 0

* - (nwp_phy_nml-inwp_surface)=
    **inwp_surface**
  - I (max_ dom)
  - 1
  -
  - surface scheme
    - 0: none
    - 1: TERRA
    - 2: JSBACH (requires inwp_turb = 6)
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-ustart_raylfric)=
    **ustart_raylfric**
  - R
  - 160.0
  - m/s
  - wind speed at which extra Rayleigh friction starts
  - inwp_gwd > 0

* - (nwp_phy_nml-efdt_min_raylfric)=
    **efdt_min_raylfric**
  - R
  - 10800.0
  - s
  - minimum e-folding time of Rayleigh friction (effective for u > ustart_raylfric + 90 m/s)
  - inwp_gwd > 0

* - (nwp_phy_nml-latm_above_top)=
    **latm_above_top**
  - L (max_ dom)
  - .FALSE.
  -
  - .TRUE.: take into account atmosphere above model top for radiation computation
  - inwp_radiation > 0

* - (nwp_phy_nml-itype_z0)=
    **itype_z0**
  - I
  - 2
  -
  - Type of roughness length data used for turbulence scheme:
    - 1: land-cover-related roughness including contribution from sub-scale orography (does not account for tiles)
    - 2: land-cover-related roughness based on tile-specific landuse class
    - 3: land-cover-related roughness based on tile-specific landuse class including contribution from sub-scale orography
  - inwp_turb > 0

* - (nwp_phy_nml-itype_satpres_coeffs)=
    **itype_satpres_coeffs**
  - I
  - 1
  -
  - Set of coefficients used for computing the saturation vapor pressure:
    - 1: Coefficients inherited from the COSMO model (inaccurate at temperatures below {math}`-50^\circ`C)
    - 2: Coefficients used in the IFS (and the parameterizations taken over from the IFS)
  - inwp_satad > 0

* - (nwp_phy_nml-dt_conv)=
    **dt_conv**
  - R (max_ dom)
  - 600.0
  - s
  - time interval of convection call.
    by default, each subdomain has the same value
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-dt_ccov)=
    **dt_ccov**
  - R (max_ dom)
  - dt_conv
  - s
  - time interval of cloud-cover call.
    by default, dt_ccov equals dt_conv for each domain
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-dt_rad)=
    **dt_rad**
  - R (max_ dom)
  - 1800.0
  - s
  - time interval of radiation call
    by default, each subdomain has the same value
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-dt_sso)=
    **dt_sso**
  - R (max_ dom)
  - 1200.0
  - s
  - time interval of sso call
    by default, each subdomain has the same value
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-dt_gwd)=
    **dt_gwd**
  - R (max_ dom)
  - 1200.0
  - s
  - time interval of gwd call
    by default, each subdomain has the same value
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-lrtm_filename)=
    **lrtm_filename**
  - C(:)
  - "rrtmg_ lw.nc"
  -
  - NetCDF file containing longwave absorption coefficients and other data for RRTMG_LW k-distribution model.
  -

* - (nwp_phy_nml-cldopt_filename)=
    **cldopt_filename**
  - C(:)
  - "ECHAM 6_CldOpt Props.nc"
  -
  - NetCDF file with RRTM Cloud Optical Properties for ECHAM6.
  -

* - (nwp_phy_nml-icalc_reff)=
    **icalc_reff**
  - I (max_ dom)
  - 0
  -
  - Parameterization set for diagnostic calculations of effective radius:
    - 0: No calculation
    - 1,2,4,5,6,7: Consistent with microphysics given by icalc_reff (naming same convention as inwp_gscp)
    - 100: Consistent with current microphysics (it sets icalc_reff = inwp_gscp)
    - 101: Reff given by RRTM parameterization
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_phy_nml-icpl_rad_reff)=
    **icpl_rad_reff**
  - I (max_ dom)
  - 0
  -
  - Coupling of the effective radius with radiation:
    - 0: No coupling. The calculation of the effective radius happens at the radiation interface.
    - 1: Radiation uses the effective radius defined by icalc_reff. All hydrometeors are combined in a frozen and a liquid phase.
    - 2: Radiation uses the effective radius defined by icalc_reff for  all hydrometeors independently (given by ecrad_iset_genhyd).
  - [iforcing](run_nml-iforcing) = inwp
    inwp_radiation = 1 or 4
    icalc_reff > 0

* - (nwp_phy_nml-ithermo_water)=
    **ithermo_water**
  - I (max_ dom)
  - 0
  -
  - Latent Heat Function
    - 0: Temperature-dependent latent heat in saturation adjustment but constant in microphysics:
    - 1: Temperature-dependent latent heat in saturation adjustment and microphysics
  - [iforcing](run_nml-iforcing) = inwp
    inwp_gscp = 1,2,4,5,7

* - (nwp_phy_nml-lupatmo_phy)=
    **lupatmo_phy**
  - L (max_ dom)
  - .FALSE.
  -
  - Switch for upper-atmosphere physics.
    Examples of usage for multi-domain applications:
      - set lupatmo_phy = .TRUE. to switch on upatmo physics for all domains
      - set lupatmo_phy = .TRUE., .TRUE., .FALSE. to switch on upatmo physics for dom 1 and 2, but switch them off for dom 3
      - please note that "skipping" domains is currently not possible, i.e. lupatmo_phy = .TRUE., .FALSE., .TRUE. is transformed into lupatmo_phy = .TRUE., .FALSE., .FALSE.
    See [upatmo_nml](ref_buildrun_nml_upatmo_nml) for configuration of the upper-atmosphere physics parameterizations.
  - [iforcing](run_nml-iforcing) = inwp
    init_mode = 1, 2, 3, or 7
    inwp_turb > 0
    inwp_radiation > 0

* - (nwp_phy_nml-lcuda_graph_turb_tran)=
    **lcuda_graph_turb_tran**
  - L
  - .FALSE.
  -
  - Activate cuda graphs for turbulent transfers. Automatically set to .FALSE. if not compiled with the ICON_USE_CUDA_GRAPH cpp key.
  - ICON_USE_CUDA_GRAPH activated

* - (nwp_phy_nml-lstoch_pattern_generator)=
    **lstoch_pattern_generator**
  - L
  - .FALSE.
  -
  - Activate stochastic pattern generator using either spherical harmonics (global) or Fourier modes (LAM)
  -

* - (nwp_phy_nml-spg_fourier_modes)=
    **spg_fourier_modes**
  - L
  - .TRUE.
  -
  - Use Fourier modes for stochastic pattern generator in LAM configurations. For LAM simulations spherical harmonics can be used by setting this switch to .FALSE. but this is computationally inefficient and only useful for tests or consistency with a global configuration.
  - lstoch_pattern_generator=.TRUE.
    and l_limited_area=.TRUE.

* - (nwp_phy_nml-spg_use_asl)=
    **spg_use_asl**
  - L
  - .FALSE.
  -
  - Use Advanced Scientific Library (ASL) provided by NEC.
  - lstoch_pattern_generator=.TRUE.

* - (nwp_phy_nml-spg_num)=
    **spg_num**
  - I
  - 0
  -
  - Number of independent 2D stochastic patterns
  - lstoch_pattern_generator=.TRUE.
    Note that all patterns have the same length and time scales.

* - (nwp_phy_nml-spg_length_scale)=
    **spg_length_scale**
  - R(:)
  - 1000e3
  - m
  - Length scales of stochastic pattern
  - lstoch_pattern_generator=.TRUE.
    The maximum array size is 5. The number of length scales is independent from the number of patterns.

* - (nwp_phy_nml-spg_time_scale)=
    **spg_time_scale**
  - R(:)
  - 3600.0
  - s
  - Time scales of stochastic pattern
  - lstoch_pattern_generator=.TRUE.
    The maximum array size is 5. The number of time scales is independent from the number of patterns.

* - (nwp_phy_nml-spg_spec_modes)=
    **spg_spec_modes**
  - I(:)
  - 50
  -
  - Number of spectral modes used for stochastic pattern generator
  - lstoch_pattern_generator=.TRUE.
    The maximum array size is 5. The number of modes should adjusted to the length scales used.

* - (nwp_phy_nml-spg_variance)=
    **spg_variance**
  - R(:)
  - 1.0
  -
  - Variance of stochastic pattern in grid point space
  - lstoch_pattern_generator=.TRUE.
    The maximum array size is 5.

* - (nwp_phy_nml-itype_stoch_phys)=
    **itype_stoch_phys**
  - I
  - 0
  -
  - Type of stochastically perturbed physics based on spectral pattern generator:
    - 0: none
    - 1: CONV only
    - 2: GW only (SSO and non-oro.~GWD)
    - 3: CONV + GW using single pattern
    - 4: CONV + GW using independent patterns
  - lstoch_pattern_generator=.TRUE.
    Note that itype_stoch_phys=4 requires spg_num=2
```




Defined and used in: {{ '[src/namelists/mo_nwp_phy_nml.f90]({}/src/namelists/mo_nwp_phy_nml.f90)'.format(base_url) }}


### Notes on use of stochastic convection schemes

There are currently three stochastic convection schemes available, two versions for shallow convection and one for deep convection. Conceptually, these schemes attempt to represent that for grid box sizes smaller than the size of a typical cloud ensemble, the clouds actually populating the grid box will not be fully representative of that cloud ensemble.
The two stochastic shallow schemes (lstoch_expl, lstoch_sde) are therefore aimed at resolutions of a few kilometers (typically used for LAM simulations, where deep convection is resolved) and will in fact be automatically switched off for resolutions greater than 20km. The scheme converges to the standard Tiedtke-Bechtold mass flux scheme at resolutions sufficiently coarse, such that there is no additional gain from using the stochastic schemes. They should therefore be run with lshalloconv_only=T. A combination with the grayzone tuning (lgrayzone_deepconv) is technically possible, but not recommended as the grayzone tuning interferes with the intended behaviour of the stochastic scheme.

The stochastic deep convection scheme (lstoch_deep) is intended for resolutions where the deep convection parameterization is still active, but again, grid size is not large enough to contain a fully representative cloud ensemble (e.g. global runs with resolution on the order of 10s of kilometers). Thus the deep and shallow stochastic schemes are not intended to be used together, as the resolutions they are designed for are (mostly) mutually exclusive.

The shallow schemes should be run without resolution-dependent tuning of the convection parameters (lrestune_off=T) and with disabled mass flux limiters (lmflimiter_off=T). The mass flux limiters are in fact not fully disabled but set to values high enough to be rarely reached during shallow cloud simulations.
The deep stochastic scheme cannot be run without mass flux limiters or simulatons will become unstable.



(ref_buildrun_nml_nwp_tuning_nml)=
## nwp_tuning_nml

Please note: These tuning parameters are NOT domain specific.

**SSO (Lott and Miller)**

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nwp_tuning_nml-tune_gkwake)=
    **tune_gkwake**
  - R (max_ dom)
  - 1.5
  -
  - low level wake drag constant
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_gkdrag)=
    **tune_gkdrag**
  - R (max_ dom)
  - 0.075
  -
  - gravity wave drag constant
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_gkdrag_enh)=
    **tune_gkdrag_enh**
  - R (max_ dom)
  - 0.075
  -
  - enhanced value of gravity wave drag constant at low latitudes (needs to be actively set to a larger value than gkdrag to be effective)
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_gfrcrit)=
    **tune_gfrcrit**
  - R (max_ dom)
  - 0.4
  -
  - critical Froude number (controls depth of blocking layer)
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_grcrit)=
    **tune_grcrit**
  - R (max_ dom)
  - 0.25
  -
  - critical Richardson number (controls onset of wave breaking)
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_grcrit_enh)=
    **tune_grcrit_enh**
  - R (max_ dom)
  - 0.25
  -
  - enhanced value of critical Richardson number at low latitudes (needs to be actively set to a larger value than grcrit to be effective)
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_minsso)=
    **tune_minsso**
  - R (max_ dom)
  - 10.0
  - m
  - minimum SSO standard deviation for which SSO scheme is applied
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_minsso_gwd)=
    **tune_minsso_gwd**
  - R (max_ dom)
  - 0.0
  - m
  - minimum SSO standard deviation for which wave drag component if SSO scheme is applied (effective only if larger than minsso; the default of zero means that the parameter needs to be actively set)
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_blockred)=
    **tune_blockred**
  - R (max_ dom)
  - 100.0
  -
  - multiple of SSO standard deviation above which blocking tendency is reduced
  - [iforcing](run_nml-iforcing) = inwp
````

**GWD (Warner McIntyre)**

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nwp_tuning_nml-tune_gfluxlaun)=
    **tune_gfluxlaun**
  - R
  - 2.50e-3
  -
  - total launch momentum flux in each azimuth (rho_o x F_o)
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_gcstar)=
    **tune_gcstar**
  - R
  - 1.0
  -
  - constant in saturation wave spectrum
  - [iforcing](run_nml-iforcing) = inwp
```

**Grid scale microphysics**

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nwp_tuning_nml-tune_zceff_min)=
    **tune_zceff_min**
  - R
  - 0.01
  -
  - Minimum value for sticking efficiency
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_v0snow)=
    **tune_v0snow**
  - R
  - 25.0
  -
  - factor in the terminal velocity for snow
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_zvz0i)=
    **tune_zvz0i**
  - R
  - 1.25
  - m/s
  - Terminal fall velocity of ice
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_icesedi_exp)=
    **tune_icesedi_exp**
  - R
  - 0.33
  -
  - Exponent for density correction of cloud ice sedimentation
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_zcsg)=
    **tune_zcsg**
  - R
  - 0.5
  -
  - Efficiency for snow-to-graupel conversion by riming
  - inwp_gscp=2

* - (nwp_tuning_nml-tune_dice_conv)=
    **tune_dice_conv**
  - R
  - 100e-6
  - m
  - mean diameter of detrained cloud ice of parameterized convection
  - [iforcing](run_nml-iforcing) = inwp; inwp_convection = 1, inwp_gscp = 3

* - (nwp_tuning_nml-tune_sbmccn)=
    **tune_sbmccn**
  - R
  - 1.0
  -
  - Scaling factor (0,1] for initial aerosol concentration profile, used for comparison
    between simulations with two-moment and warm SBM microphysics
  - [iforcing](run_nml-iforcing) = inwp, inwp_gscp = 8
```

**Convection scheme (Tiedtke-Bechtold)**

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nwp_tuning_nml-tune_entrorg)=
    **tune_entrorg**
  - R
  - 1.95e-3
  - 1/m
  - Entrainment parameter valid for dx=20 km (depends on model resolution)
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_rprcon)=
    **tune_rprcon**
  - R
  - 1.4e-3
  -
  - Coefficient for conversion of cloud water into precipitation
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_rdepths)=
    **tune_rdepths**
  - R
  - 2.e4
  - Pa
  - Maximum allowed depth of shallow convection
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_capdcfac_et)=
    **tune_capdcfac_et**
  - R
  - 0.5
  -
  - Fraction of CAPE diurnal cycle correction applied in the extratropics
  - icapdcycl = 3

* - (nwp_tuning_nml-tune_capethresh)=
    **tune_capethresh**
  - R
  - 7000.0
  - J/kg
  - CAPE threshold above which the convective adjustment time scale and entrainment rate are reduced for numerical stability
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_rhebc_land)=
    **tune_rhebc_land**
  - R
  - 0.75
  -
  - RH threshold for onset of evaporation below cloud base over land
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_rhebc_land_trop)=
    **tune_rhebc_land_trop**
  - R
  - 0.75
  -
  - RH threshold for onset of evaporation below cloud base over land in the tropics
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_rhebc_ocean)=
    **tune_rhebc_ocean**
  - R
  - 0.85
  -
  - RH threshold for onset of evaporation below cloud base over sea
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_rhebc_ocean_trop)=
    **tune_rhebc_ocean_trop**
  - R
  - 0.80
  -
  - RH threshold for onset of evaporation below cloud base over sea in the tropics
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_rcucov)=
    **tune_rcucov**
  - R
  - 0.05
  -
  - Convective area fraction used for computing evaporation below cloud base
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_rcucov_trop)=
    **tune_rcucov_trop**
  - R
  - 0.05
  -
  - Convective area fraction used for computing evaporation below cloud base in the tropics
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_texc)=
    **tune_texc**
  - R
  - 0.125
  - K
  - Excess value for temperature used in test parcel ascent
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_qexc)=
    **tune_qexc**
  - R
  - 0.0125
  -
  - Excess fraction of grid-scale QV used in test parcel ascent
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_rmfdeps_land)=
    **tune_rmfdeps_land**
  - R
  - 0.25
  -
  - fractional mass flux for downdrafts over land
  - [iforcing](run_nml-iforcing) = inwp;  lshallowconv_only = .FALSE.; lgrayzone_deepconv = .FALSE;

* - (nwp_tuning_nml-tune_rmfdeps_ocean)=
    **tune_rmfdeps_ocean**
  - R
  - 0.15
  -
  - fractional mass flux for downdrafts over ocean
  - [iforcing](run_nml-iforcing) = inwp;  lshallowconv_only = .FALSE.; lgrayzone_deepconv = .FALSE;

* - (nwp_tuning_nml-tune_detrainment_profile)=
    **tune_detrainment_profile**
  - R(2)
  - 1.0, 0.0
  -
  - prefactor in RH-dependent detrainment profile, formulated as {math}`\mathrm{tunefac(1)} - \mathrm{tunefac(2)}*\mathrm{RH}`; difference between first and second value should be between 0.75 and 1
  - [iforcing](run_nml-iforcing) = inwp;

* - (nwp_tuning_nml-tune_entrainment_profile)=
    **tune_entrainment_profile**
  - R(2)
  - 1.3, 1.0
  -
  - prefactor in RH-dependent entrainment profile, formulated as {math}`\mathrm{tunefac(1)} - \mathrm{tunefac(2)}*\mathrm{RH}`; difference between first and second value should be about 0.3
  - [iforcing](run_nml-iforcing) = inwp;

* - (nwp_tuning_nml-tune_grzdc_offset)=
    **tune_grzdc_offset**
  - R
  - 0.0
  -
  - Scaling factor for offset in CAPE closure for grayzone deep convection. Positive values reduce the activity of the convection scheme and suppress convective drizzle (recommendation: 0.1--0.2)
  - [iforcing](run_nml-iforcing) = inwp; lgrayzone_deepconv = .TRUE.
```

**Cloud scheme**

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nwp_tuning_nml-tune_tau_shallow)=
    **tune_tau_shallow**
  - R
  - 1500.0
  - s
  - Decay time scale for shallow convective anvils in cloud cover scheme
  - [iforcing](run_nml-iforcing) = inwp; inwp_cldcover = 1

* - (nwp_tuning_nml-tune_tau_mid)=
    **tune_tau_mid**
  - R
  - 1500.0
  - s
  - Decay time scale for mid-level convective anvils in cloud cover scheme
  - [iforcing](run_nml-iforcing) = inwp; inwp_cldcover = 1

* - (nwp_tuning_nml-tune_tau_deep)=
    **tune_tau_deep**
  - R
  - 1500.0
  - s
  - Decay time scale for deep convective anvils in cloud cover scheme
  - [iforcing](run_nml-iforcing) = inwp; inwp_cldcover = 1

* - (nwp_tuning_nml-tune_box_liq)=
    **tune_box_liq**
  - R
  - 0.05
  -
  - Box width scale for liquid cloud diagnostic in cloud cover scheme
  - [iforcing](run_nml-iforcing) = inwp; inwp_cldcover = 1

* - (nwp_tuning_nml-tune_box_ice)=
    **tune_box_ice**
  - R
  - 0.05
  -
  - Box width scale for ice cloud diagnostic in cloud cover scheme
  - [iforcing](run_nml-iforcing) = inwp; inwp_cldcover = 1

* - (nwp_tuning_nml-tune_thicklayfac)=
    **tune_thicklayfac**
  - R
  - 0.005
  - 1/m
  - Factor for enhancing the box width for model layer thicknesses exceeding 150 m
  - [iforcing](run_nml-iforcing) = inwp; inwp_cldcover = 1

* - (nwp_tuning_nml-tune_box_liq_asy)=
    **tune_box_liq_asy**
  - R
  - 2.5
  -
  - Asymmetry factor for liquid cloud cover diagnostic
  - [iforcing](run_nml-iforcing) = inwp; inwp_cldcover = 1

* - (nwp_tuning_nml-tune_box_liq_sfc_fac)=
    **tune_box_liq_sfc_fac**
  - R (max_dom)
  - 1.0
  -
  - Tuning factor for box_liq reduction near the surface
  - [iforcing](run_nml-iforcing) = inwp; inwp_cldcover = 1

* - (nwp_tuning_nml-allow_overcast)=
    **allow_overcast**
  - R
  - 1.0
  -
  - Tuning factor for the dependence of liquid cloud cover on relative humidity. This is an unphysical ad-hoc parameter to improve the cloud cover in the Mediterranean
  - [iforcing](run_nml-iforcing) = inwp; inwp_cldcover = 1

* - (nwp_tuning_nml-tune_sgsclifac)=
    **tune_sgsclifac**
  - R
  - 0.0
  -
  - Scaling factor for parameterization of subgrid-scale (turbulence-induced) cloud ice (values > 0 not recommended for global configurations with RRTM radiation)
  - [iforcing](run_nml-iforcing) = inwp; inwp_cldcover = 1

* - (nwp_tuning_nml-icpl_turb_clc)=
    **icpl_turb_clc**
  - I
  - 1
  -
  - Mode of coupling between turbulence and cloud cover
    - 1: strong dependency of box width on rcld with upper and lower limit
    - 2: weak dependency of box width on rcld with additive term and upper limit
  - [iforcing](run_nml-iforcing) = inwp; inwp_cldcover = 1

* - (nwp_tuning_nml-lcalib_clcov)=
    **lcalib_clcov**
  - L
  - .TRUE.
  -
  - Apply calibration of layer-wise cloud cover diagnostics over land in order to improve scores against SYNOP reports
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-max_calibfac_clcl)=
    **max_calibfac_clcl**
  - R
  - 4.0
  -
  - Maximum allowed calibration factor for low clouds (CLCL)
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_eiscrit)=
    **tune_eiscrit**
  - R
  - 1000.0
  - K
  - Critical estimated inversion strength above which to switch off shallow convection (recommendation to activate: 7K)
  - [iforcing](run_nml-iforcing) = inwp; inwp_convection = 1

* - (nwp_tuning_nml-tune_sc_eis)=
    **tune_sc_eis**
  - R
  - 1000.0
  - K
  - Critical estimated inversion strength above which to use modified SGS cloud diagnostic for stratocumulus (recommendation to activate: 7K)
  - [iforcing](run_nml-iforcing) = inwp; inwp_clcover = 1

* - (nwp_tuning_nml-tune_sc_invmin)=
    **tune_sc_invmin**
  - R
  - 200.0
  - m
  - Minimum inversion height above which to apply modified SGS cloud diagnostic for stratocumulus
  - [iforcing](run_nml-iforcing) = inwp; inwp_clcover = 1

* - (nwp_tuning_nml-tune_sc_invmax)=
    **tune_sc_invmax**
  - R
  - 1500.0
  - m
  - Maximum inversion height below which to apply modified SGS cloud diagnostic for stratocumulus
  - [iforcing](run_nml-iforcing) = inwp; inwp_clcover = 1

* - (nwp_tuning_nml-tune_cu_alfa)=
    **tune_cu_alfa**
  - R
  - 0.0
  -
  - scaling parameter of modification of cloud fraction of shallow liquid clouds
  - [iforcing](run_nml-iforcing) = inwp; inwp_clcover = 1;  icpl_aero_gscp > 0

* - (nwp_tuning_nml-tune_cu_cdnc)=
    **tune_cu_cdnc**
  - R
  - 150e6
  - m{math}`^{-3}`
  - cloud droplet number threshold to identify regions for cloud cover modification
  - [iforcing](run_nml-iforcing) = inwp; inwp_clcover = 1
```

**Effective radius**

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nwp_tuning_nml-tune_reff_qi)=
    **tune_reff_qi**
  - R
  - 1.0
  -
  - Linear tuning factor for effective radius of clouds (affects both, grid-scale and sub-grid)
  - [iforcing](run_nml-iforcing) = inwp; icpl_rad_reff = 1; inwp_gscp=3
```

**Saturation adjustment**

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nwp_tuning_nml-tune_supsat_limfac)=
    **tune_supsat_limfac**
  - R
  - 0.0
  -
  - Limiting factor for parameterized supersaturation in updrafts during the saturation adjustment (0.0 for thermodynamic equilibrium)
  - [iforcing](run_nml-iforcing) = inwp; inwp_satad = 1
```

**Misc**

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nwp_tuning_nml-tune_gust_factor)=
    **tune_gust_factor**
  - R
  - 8.0
  -
  - Multiplicative factor for friction velocity in gust parameterization
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-tune_gustsso_lim)=
    **tune_gustsso_lim**
  - R
  - 100.0
  - m s{math}`^{-1}`
  - Basic gust speed at which the SSO correction starts to be reduced (recommendation to activate: 20 m s{math}`^{-1}`
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-itune_gust_diag)=
    **itune_gust_diag**
  - I
  - 1
  -
  - Method of SSO blocking correction used in the gust diagnostics
    - 1: Use level above "SSO envelope top" for gust enhancement over mountains
    - 2: Use "SSO envelope top" level for gust enhancement over mountains, combined with an adjusted nonlinearity factor (recommended for global configurations with MERIT/REMA orography)
    - 3: Variant of option 1, recommended for ICON-D2 with subgrid-scale condensation (do not use with ntiles=1)
    - 4: As 3, but using time-averaged 10-m wind speeds as input with additional limitation to resolved bouldary-layer wind speeds
  - [iforcing](run_nml-iforcing) = inwp; related switches are tune_gustlim_agl and tune_gustlim_fac

* - (nwp_tuning_nml-tune_gustlim_agl)=
    **tune_gustlim_agl**
  - R(max_dom)
  - 1500.0
  - m
  - Height range above ground, within which the maximum resolved wind speed is determined for gust limitation
  - itune_gust_diag = 4

* - (nwp_tuning_nml-tune_gustlim_fac)=
    **tune_gustlim_fac**
  - R(max_dom)
  - 0.0
  -
  - Tuning factor for gust limitation. The default value of zero means that the limitation is turned off. Otherwise, the difference between the 10-m wind speed and the maximum speed found below tune_gustlim_agl times tune_gustlim_fac is used to limit the excess gust speed
  - itune_gust_diag = 4

* - (nwp_tuning_nml-tune_ssolim_sfcfric)=
    **tune_ssolim_sfcfric**
  - R(max_dom)
  - 10000.0
  - m
  - Threshold value of SSO standard deviation above which the adaptive parameter tuning for surface friction is reduced (with a linear transition from 1 to 0 between the specified value and twice this value). Recommendation: use a value slightly above the threshold used in the data assimilation scheme for 10-m winds (70 m at DWD).
  - icpl_da_sfcfric {math}`\ge` 1; inwp_sso = 1

* - (nwp_tuning_nml-itune_vis_diag)=
    **itune_vis_diag**
  - I
  - 1
  -
  - Tuning variant of visibility diagnostics
    - 1: First operational implementation
    - 2: Optimized day-night conversion factor
    - 3: (2) plus optimized RH-dependency
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-itune_ceiling_diag)=
    **itune_ceiling_diag**
  - I
  - 1
  -
  - Method used for ceiling diagnostics
    - 1: Purely layer-wise diagnostic
    - 2: Vertical integral with overlap assumption like in cloud-cover diagnostic
  - [iforcing](run_nml-iforcing) = inwp

* - (nwp_tuning_nml-itune_albedo)=
    **itune_albedo**
  - I
  - 0
  -
  - MODIS albedo tuning
    - 0: None
    - 1: dimmed sahara
    - 2: dimmed sahara + brightened Antarctic (by {math}`5\%`)
  - [iforcing](run_nml-iforcing) = inwp
    albedo_type=2

* - (nwp_tuning_nml-tune_albedo_wso)=
    **tune_albedo_wso**
  - R(2)
  - 0.0, 0.0
  -
  - Add a correction to MODIS albedo over [dry,wet] soil for soil types 3-6. Valid range: [-0.03, 0.03]. Supposed to be negative for wet soil.
  - [iforcing](run_nml-iforcing) = inwp
    inwp_surface = 1
    itype_albedo = 2
    direct_albedo = 3,4

* - (nwp_tuning_nml-itune_slopecorr)=
    **itune_slopecorr**
  - I
  - 0
  -
  - Tuning measures for high-resolution configurations with mesh sizes around or below 1 km  
    - 0: None  
    - 1: Slope-dependent reduction of laminar transfer resistance and near-surface minimum vertical diffusion
  - [iforcing](run_nml-iforcing) = inwp  
    inwp_turb = 1

* - (nwp_tuning_nml-itune_o3)=
    **itune_o3**
  - I
  - 2
  -
  - Ozone tuning  
    - 0: no tuning  
    - 1: old tuning for RRTM radiation  
    - 2: standard tuning for EcRad with RRTM gas optics  
    - 3: updated variant of 2 for combination with icpl_o3_tp = 2  
    - 4: provisional tuning for EcRad with EcCKD gas optics
  - [iforcing](run_nml-iforcing) = inwp  
    irad_o3 = 7, 79 or 97

* - (nwp_tuning_nml-tune_difrad_3dcont)=
    **tune_difrad_3dcont**
  - R
  - 0.5
  -
  - Tuning factor for 3D contribution to diagnosed diffuse radiation (no impact on prognostic results!)
  - inwp_radiation = 1 or 4

* - (nwp_tuning_nml-tune_minsnowfrac)=
    **tune_minsnowfrac**
  - R
  - 0.2
  -
  - Minimum value to which the snow cover fraction is artificially reduced in case of melting snow
  - idiag_snowfrac = 20/30/40

* - (nwp_tuning_nml-tune_dursun_scaling)=
    **tune_dursun_scaling**
  - R
  - 1.0
  -
  - Tuning factor for direct solar irradiance in sunshine duration diagnostic to account for the delta-Eddington scaling in ecRad and other possible biases (e.g. liquid/ice water path)
  -

* - (nwp_tuning_nml-tune_urbahf)=
    **tune_urbahf**
  - R(4)
  - 0., 2., 2., 50.
  - W m{math}`^{-2}`
  - Tuning factors for specifying the anthropogenic heat flux (AHF) depending on climatological T2M.  
    First value: constant base value independent of temperature  
    Second value: gradient per K below T2M = 15{math}`^\circ`C for heating  
    Third value: gradient per K above T2M = 20{math}`^\circ`C for cooling  
    Fourth value: upper limit for AHF
  - lterra_urb = .TRUE.

* - (nwp_tuning_nml-tune_urbisa)=
    **tune_urbisa**
  - R(2)
  - 0.6, 1.0
  -
  - Lower and upper bound for variable ISA parameterization (fraction of impervious surface area on urban tiles) depending on smoothed urban fraction
  - lterra_urb = .TRUE.

* - (nwp_tuning_nml-tune_cdnc)=
    **tune_cdnc**
  - R(3)
  - 10e6, 300e6, 1.0
  - m{math}`^{-3}`, m{math}`^{-3}`, -
  - Lower and upper bound and scaling factor for external cloud droplet number concentration climatology
  - icpl_aero_gscp = 3.

* - (nwp_tuning_nml-tune_demax_hail_s)=
    **tune_demax_hail_s**
  - R
  - 1.6
  -
  - In case of two-moment microphysics, for optional diagnostic maximum estimated size of hail (output fields `demax_hail_s`, `demax_hail_tmax_s`): scaling factor for computing the estimated maximum hail diameter at the surface from the mean mass diameter of hail.
  - inwp_gscp = 4,5,6,7

* - (nwp_tuning_nml-prhthresh_demax_hail_s)=
    **prhthresh_demax_hail_s**
  - R
  - 5e-4
  - kg m{math}`^{-2}` s{math}`^{-1}`
  - In case of two-moment microphysics, for optional diagnostic maximum estimated size of hail (`demax_hail_s`, `demax_hail_tmax_s`): hail rate threshold above which to compute the estimated maximum hail diameter. For smaller hail rate, the diameter is set to 0.0.
  - inwp_gscp = 4,5,6,7

* - (nwp_tuning_nml-kefthresh_demax_hail_s)=
    **kefthresh_demax_hail_s**
  - R
  - 1.5e-1
  - W m{math}`^{-2}`
  - In case of two-moment microphysics, for optional diagnostic maximum estimated size of hail (`demax_hail_s`, `demax_hail_tmax_s`): hail kinetic energy flux threshold above which to compute the estimated maximum hail diameter. For smaller flux, the diameter is set to 0.0.
  - inwp_gscp = 4,5,6,7

* - (nwp_tuning_nml-qnhthresh_demax_hail_s)=
    **qnhthresh_demax_hail_s**
  - R
  - 1.5e-2
  - m{math}`^{-3}`
  - In case of two-moment microphysics, for optional diagnostic maximum estimated size of hail (`demax_hail_s`, `demax_hail_tmax_s`): hail number density threshold above which to compute the estimated maximum hail diameter. For smaller number, the diameter is set to 0.0.
  - inwp_gscp = 4,5,6,7

* - (nwp_tuning_nml-lwindeffect_kef_hail_s)=
    **lwindeffect_kef_hail_s**
  - L
  - .FALSE.
  -
  - In case of two-moment microphysics, for optional diagnostic hail kinetic energy (output fields `kef_hail_s`, `kef_hail_max_s`, `ke_hail_s`): switch to take into account near-surface wind speed effects.
  - inwp_gscp = 4,5,6,7
```

**IAU**

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nwp_tuning_nml-max_freshsnow_inc)=
    **max_freshsnow_inc**
  - R
  - 0.025
  -
  - Maximum allowed freshsnow increment per analysis cycle (positive or negative)
  - init_mode = 5 (MODE_IAU)
```

Defined and used in: {{ '[src/namelists/mo_nwp_tuning_nml.f90]({}/src/namelists/mo_nwp_tuning_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_twomom_mcrph_nml)=
## twomom_mcrph_nml

This namelist offers the possibility to adapt some configuration parameters of the two-moment cloud microphysical parameterisation by A. Seifert and K.D. Beheng. It is only effective if this scheme is used, i.e., if **inwp_gscp=4, 5, 6, or 7**.

The below set of parameters is a first reasonable choice to start with something. There might be coming more parameters in the future.

**Please note:** at the moment we do not support the option to have different configuration parameters on different domains.
We did not really test this up to now, but it cannot be ruled out for the future.
There are for sure parameters which could be optimized for different resolutions. Therefore, at the moment this possibility is not provided explicitly to the user via the below namelist parameters (they are scalars), but it is prepared internally in the ICON code.



```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (twomom_mcrph_nml-i2mom_solver)=
    **i2mom_solver**
  - I
  - 1
  -
  - Type of numerical time integration scheme for the two-moment scheme:
    - 0: explicit Euler-forward
    - 1: semi-implicit solver similar to that of the standard one-moment schemes
  - [iforcing](run_nml-iforcing)=3  inwp_gscp=4

* - (twomom_mcrph_nml-ccn_type)=
    **ccn_type**
  - I
  - Depends on inwp_gscp:
    for 4,7: 7
    for 5: 8
  -
  - Choice of the aerosol scenario for cloud nucleation (CCN):
    - 6: "low CCN" ("maritime")
    - 7: "intermediate CCN"
    - 8: "high CCN" ("continental")
    - 9: "very high CCN" ("polluted")
    If applied together with the ART aerosol physics **inwp_gscp=6**, this parameter has no effect.
  - [iforcing](run_nml-iforcing)=3  inwp_gscp=4,5,7

* - (twomom_mcrph_nml-ccn_Ncn0)=
    **ccn_Ncn0**
  - R
  - -999.99
  - m{math}`^{-3}`
  - CN concentration near ground. A value of < -900 indicates that the hardcoded value associated with the ccn_type will be used. If applied together with the ART aerosol physics **inwp_gscp=6**, this parameter has no effect.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,7

* - (twomom_mcrph_nml-ccn_wcb_min)=
    **ccn_wcb_min**
  - R
  - 0.1
  - m{math}`^{-3}`
  - Minimum updraft speed for Segal & Khain cloud nucleation parameterization. If applied together with the ART aerosol physics **inwp_gscp=6**, this parameter has no effect.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,7

* - (twomom_mcrph_nml-iicephase)=
    **iicephase**
  - I
  - 1
  -
  - Turning on/off mixed phase processes in the two-moment scheme:
    - 0: warm phase processes only
    - 1: mixed phase processes
  - [iforcing](run_nml-iforcing)=3,  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-alpha_spacefilling)=
    **alpha_spacefilling**
  - R
  - 0.01
  -
  - Parameter in conversion of snow or cloud ice to graupel by riming: degree of void filling by frozen supercooled droplets within the ice particle skeleton, above which the particle is converted to the graupel class. Smaller values lead to faster conversion. 0.01 means very fast conversion to graupel. A value of 0.68 is the theoretical limit for densest sphere packing and leads to rather slow conversion.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-D_conv_ii)=
    **D_conv_ii**
  - R
  - 75.0e-6
  - m
  - diameter threshold for the onset of conversion to snow by ice selfcollection
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-D_rainfrz_ig)=
    **D_rainfrz_ig**
  - R
  - 0.50e-3
  - m
  - Spectral size threshold below which freezing rain drops are converted to cloud ice. Larger drops are converted to graupel or hail, depending on parameter D_rainfrz_gh.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-D_rainfrz_gh)=
    **D_rainfrz_gh**
  - R
  - 1.25e-3
  - m
  - Spectral size threshold above which freezing rain drops are converted to hail. Smaller drops are converted to cloud ice or graupel, depending on parameter D_rainfrz_ig.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-luse_mu_Dm_rain)=
    **luse_mu_Dm_rain**
  - L
  - .FALSE.
  -
  - To switch on the usage of the dynamical {math}`\mu`-{math}`D_M`-relation for raindrops below cloud base.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-rain_cmu0)=
    **rain_cmu0**
  - R
  - 6.0
  -
  - Parameter of the {math}`\mu`-{math}`D`-relation in the rain size distribution for evaporation and sedimentation below cloud base:
    asymptotic {math}`\mu`-value for spectra with small mean diameter.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-rain_cmu1)=
    **rain_cmu1**
  - R
  - 30.0
  -
  - Parameter of the {math}`\mu`-{math}`D`-relation in the rain size distribution for evaporation and sedimentation below cloud base:
    asymptotic {math}`\mu`-value for spectra with large mean diameter.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-rain_cmu3)=
    **rain_cmu3**
  - R
  - 1.1e-3
  - m
  - Parameter of the {math}`\mu`-{math}`D`-relation in the rain size distribution for evaporation and sedimentation below cloud base:
    equilibrium mean spectral diameter for breakup and selfcollection.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-rain_cmu4)=
    **rain_cmu4**
  - R
  - 1.0
  -
  - Parameter of the {math}`\mu`-{math}`D`-relation in the rain size distribution for evaporation and sedimentation below cloud base:
    base value of {math}`\mu`.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-in_fact)=
    **in_fact**
  - R
  - 1.0
  -
  - Factor to tune the IN concentration for heterogeneous ice nucleation
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-nu_i)=
    **nu_i**
  - R
  - -999.99
  -
  - Shape parameter {math}`\nu` for cloud ice in the PSD
    {math}`f(x)=N_{0}x^{\nu}\exp(-\lambda x^{\mu})`
    A value of < -900  indicates that the background value {math}`\nu=0.0` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-mu_i)=
    **mu_i**
  - R
  - -999.99
  -
  - Shape parameter {math}`\mu` for cloud ice in the PSD
    {math}`f(x)=N_{0}x^{\nu}\exp(-\lambda x^{\mu})`
    A value of < -900  indicates that the background value {math}`\mu=1/3` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-ageo_i)=
    **ageo_i**
  - R
  - -999.99
  -
  - Prefactor of the assumed size-mass-relation for cloud ice {math}`D=a_{geo}x^{b_{geo}}` for {math}`x` in kg and {math}`D` in m.
    A value of < -900  indicates that the background value {math}`a_{geo}=0.835` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-bgeo_i)=
    **bgeo_i**
  - R
  - -999.99
  -
  - Exponent of the assumed size-mass-relation for cloud ice {math}`D=a_{geo}x^{b_{geo}}` for {math}`x` in kg and {math}`D` in m.
    A value of < -900  indicates that the background value {math}`b_{geo}=0.39` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-avel_i)=
    **avel_i**
  - R
  - -999.99
  -
  - Prefactor of the assumed fallspeed-mass-relation for cloud ice {math}`v=a_{vel}x^{b_{vel}}` for {math}`x` in kg and {math}`v` in m/s.
    A value of < -900  indicates that the background value {math}`a_{vel}=27.7` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-bvel_i)=
    **bvel_i**
  - R
  - -999.99
  -
  - Exponent of the assumed fallspeed-mass-relation for cloud ice {math}`v=a_{vel}x^{b_{vel}}` for {math}`x` in kg and {math}`v` in m/s.
    A value of < -900  indicates that the background value {math}`b_{vel}=0.21579` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-nu_s)=
    **nu_s**
  - R
  - -999.99
  -
  - Shape parameter {math}`\nu` for snow in the PSD
    {math}`f(x)=N_{0}x^{\nu}\exp(-\lambda x^{\mu})`
    A value of < -900  indicates that the background value {math}`\nu=0.0` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-mu_s)=
    **mu_s**
  - R
  - -999.99
  -
  - Shape parameter {math}`\mu` for snow in the PSD
    {math}`f(x)=N_{0}x^{\nu}\exp(-\lambda x^{\mu})`
    A value of < -900  indicates that the background value {math}`\mu=0.5` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-ageo_s)=
    **ageo_s**
  - R
  - -999.99
  -
  - Prefactor of the assumed size-mass-relation for snow {math}`D=a_{geo}x^{b_{geo}}` for {math}`x` in kg and {math}`D` in m.
    A value of < -900  indicates that the background value {math}`a_{geo}=5.13` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-bgeo_s)=
    **bgeo_s**
  - R
  - -999.99
  -
  - Exponent of the assumed size-mass-relation for snow {math}`D=a_{geo}x^{b_{geo}}` for {math}`x` in kg and {math}`D` in m.
    A value of < -900  indicates that the background value {math}`b_{geo}=1/2` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-avel_s)=
    **avel_s**
  - R
  - -999.99
  -
  - Prefactor of the assumed fallspeed-mass-relation for snow {math}`v=a_{vel}x^{b_{vel}}` for {math}`x` in kg and {math}`v` in m/s.
    A value of < -900  indicates that the background value {math}`a_{vel}=400.0` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-bvel_s)=
    **bvel_s**
  - R
  - -999.99
  -
  - Exponent of the assumed fallspeed-mass-relation for snow {math}`v=a_{vel}x^{b_{vel}}` for {math}`x` in kg and {math}`v` in m/s.
    A value of < -900  indicates that the background value {math}`b_{vel}=0.35` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-nu_r)=
    **nu_r**
  - R
  - -999.99
  -
  - Shape parameter {math}`\nu` of the rain mass distribution inside clouds. Refers to the generalized gamma distribution with respect to mass {math}`x`:
    {math}`f(x)=N_{0}x^{\nu}\exp(-\lambda x^{\mu})`
    A value < -900 indicates that the background value {math}`\nu=0.0` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-nu_g)=
    **nu_g**
  - R
  - -999.99
  -
  - Shape parameter {math}`\nu` for graupel in the PSD
    {math}`f(x)=N_{0}x^{\nu}\exp(-\lambda x^{\mu})`
    A value of < -900  indicates that the background value {math}`\nu=1.0` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-mu_g)=
    **mu_g**
  - R
  - -999.99
  -
  - Shape parameter {math}`\mu` for graupel in the PSD
    {math}`f(x)=N_{0}x^{\nu}\exp(-\lambda x^{\mu})`
    A value of < -900  indicates that the background value {math}`\mu=1/3` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3 inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-ageo_g)=
    **ageo_g**
  - R
  - -999.99
  -
  - Prefactor of the assumed size-mass-relation for graupel {math}`D=a_{geo}x^{b_{geo}}` for {math}`x` in kg and {math}`D` in m.
    A value of < -900  indicates that the background value {math}`a_{geo}=0.124` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-bgeo_g)=
    **bgeo_g**
  - R
  - -999.99
  -
  - Exponent of the assumed size-mass-relation for graupel {math}`D=a_{geo}x^{b_{geo}}` for {math}`x` in kg and {math}`D` in m.
    A value of < -900  indicates that the background value {math}`b_{geo}=0.314` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-avel_g)=
    **avel_g**
  - R
  - -999.99
  -
  - Prefactor of the assumed fallspeed-mass-relation for graupel {math}`v=a_{vel}x^{b_{vel}}` for {math}`x` in kg and {math}`v` in m/s.
    A value of < -900  indicates that the background value {math}`a_{vel}=100.0` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-bvel_g)=
    **bvel_g**
  - R
  - -999.99
  -
  - Exponent of the assumed fallspeed-mass-relation for graupel {math}`v=a_{vel}x^{b_{vel}}` for {math}`x` in kg and {math}`v` in m/s.
    A value of < -900  indicates that the background value {math}`b_{vel}=0.34` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-nu_h)=
    **nu_h**
  - R
  - -999.99
  -
  - Shape parameter {math}`\nu` for hail in the PSD
    {math}`f(x)=N_{0}x^{\nu}\exp(-\lambda x^{\mu})`
    A value of < -900  indicates that the background value {math}`\nu=1.0` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-mu_h)=
    **mu_h**
  - R
  - -999.99
  -
  - Shape parameter {math}`\mu` for hail in the PSD
    {math}`f(x)=N_{0}x^{\nu}\exp(-\lambda x^{\mu})`
    A value of < -900  indicates that the background value {math}`\mu=1/3` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-ageo_h)=
    **ageo_h**
  - R
  - -999.99
  -
  - Prefactor of the assumed size-mass-relation for hail {math}`D=a_{geo}x^{b_{geo}}` for {math}`x` in kg and {math}`D` in m.
    A value of < -900  indicates that the background value {math}`a_{geo}=0.1366` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-bgeo_h)=
    **bgeo_h**
  - R
  - -999.99
  -
  - Exponent of the assumed size-mass-relation for hail {math}`D=a_{geo}x^{b_{geo}}` for {math}`x` in kg and {math}`D` in m.
    A value of < -900  indicates that the background value {math}`b_{geo}=1/3` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-avel_h)=
    **avel_h**
  - R
  - -999.99
  -
  - Prefactor of the assumed fallspeed-mass-relation for hail {math}`v=a_{vel}x^{b_{vel}}` for {math}`x` in kg and {math}`v` in m/s.
    A value of < -900  indicates that the background value {math}`a_{vel}=39.3` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-bvel_h)=
    **bvel_h**
  - R
  - -999.99
  -
  - Exponent of the assumed fallspeed-mass-relation for hail {math}`v=a_{vel}x^{b_{vel}}` for {math}`x` in kg and {math}`v` in m/s.
    A value of < -900  indicates that the background value {math}`b_{vel}=1/6` from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-xmax_h)=
    **xmax_h**
  - R
  - -999.99
  - kg
  - Maximum allowed mean-mass diameter {math}`x_{max,h}` of the hail PSD. If the actual {math}`x_h` gets larger during time stepping, the hail number concentration is raised to keep {math}`x_h` below this limit. <br> A value of < -900  indicates that the background value {math}`x_{max,h}=0.005` kg from {{ '[src/atm_phy_schemes/mo_2mom_mcrph_main.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_main.f90)'.format(base_url) }} is used.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-melt_h_tune_fac)=
    **melt_h_tune_fac**
  - R
  - 1.0
  -
  - Tuning factor for the hail melting rate. Values larger than 1.0 enhance the hail melting, smaller values slow down the melting.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-melt_g_tune_fac)=
    **melt_g_tune_fac**
  - R
  - 1.0
  -
  - Tuning factor for the graupel melting rate. Values larger than 1.0 enhance the graupel melting, smaller values slow down the melting.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-Tmax_gr_rime)=
    **Tmax_gr_rime**
  - R
  - 270.16
  - K
  - Graupel formation by riming of snow and cloud ice is only active below this temperature threshold.
  - [iforcing](run_nml-iforcing)=2,3 inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-lturb_enhc)=
    **lturb_enhc**
  - L
  - .TRUE.
  -
  - To switch on the turbulent enhancement of collision processes involving water droplets (autoconversion, accretion, rain selfcollection).
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-lturb_len)=
    **lturb_len**
  - R
  - 300
  - m
  - Turbulent length scale used for lturb_enhc=.TRUE.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7 lturb_enhc=.TRUE.

* - (twomom_mcrph_nml-iice_stick)=
    **iice_stick**
  - I
  - 10
  -
  - Optional choice of the sticking efficiency parameterization for ice-ice-collisions as function of temperature:
    - 1: {math}`e = \min(\exp(0.09 (T-T_3)),1.0)` from Lin et al. (1983)
    - 2: {math}`e = \min(\exp(0.025 (T-T_3)),1.0)`, even larger as option 1
    - 3: same as 2, but {math}`e=0.01` for {math}`T<-40^{\circ}`C
    - 4: piecewise constant sticking efficiency with maxima at {math}`-7.5^{\circ}`C and {math}`-15^{\circ}`C and much smaller values at warm temperatures compared to 1, 2 and 3
    - 5: piecewise linear sticking efficiency with maximum at {math}`-15^{\circ}`C and smaller values elswhere compared to 5
    - 6: option 5 times factor 0.5
    - 7: constant {math}`e = 0.1`
    - 8: piecewise linear, a mix of 4 and 5
    - 9: option 5 times factor 0.75
    - 10: {math}`e = \min(10^{(0.035(T-T_3)-0.7)},0.2)` from Cotton et al. (1986)
    with {math}`T_3=273.16` K
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-isnow_stick)=
    **isnow_stick**
  - I
  - 5
  -
  - Optional choice of the sticking efficiency parameterization for snow-snow collisions as function of temperature.
    Same choices as for `iice_stick`.
  - [iforcing](run_nml-iforcing)=2,3 inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-iparti_stick)=
    **iparti_stick**
  - I
  - 5
  -
  - Optional choice of the sticking. efficiency parameterization for other frozen category collisions as function of temperature.
    Same choices as for `iice_stick`.
    Does not apply for graupel selfcollection. There, the below [ecoll_gg](twomom_mcrph_nml-ecoll_gg), `ecoll_gg_wet` and `Tcoll_gg_wet` apply.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-ecoll_gg)=
    **ecoll_gg**
  - R
  - 0.1
  -
  - Collision efficiency for graupel autoconversion (dry graupel). Value between 0 and 1.
  - [iforcing](run_nml-iforcing)=2,3 inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-ecoll_gg_wet)=
    **ecoll_gg_wet**
  - R
  - 0.4
  -
  - Collision efficiency for graupel autoconversion (wet graupel). Value between 0 and 1.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-Tcoll_gg_wet)=
    **Tcoll_gg_wet**
  - R
  - 270.16
  - K
  - Temperature limit above which graupel autoconversion is considered to be for wet surfaces.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-cap_ice)=
    **cap_ice**
  - R
  - -999.99
  -
  - Capacitance for clould ice depositional growth. A value < -900 indicates the usage of the code-internal backgroud value 3.0.
  - [iforcing](run_nml-iforcing)=2,3 inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-cap_snow)=
    **cap_snow**
  - R
  - -999.99
  -
  - Capacitance for snow depositional growth. A value < -900 indicates the usage of the code-internal backgroud value 3.0.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-vsedi_max_s)=
    **vsedi_max_s**
  - R
  - -999.99
  - m/s
  - Maximum allowed spectral mean sedimentation velocity of snow at sea level. A value < -900 indicates the usage of the code-internal backgroud value 1.2 m/s.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-itype_shedding_gh)=
    **itype_shedding_gh**
  - I
  - 0
  -
  - Choice of shedding parameterization for graupel and hail during collision processes with water droplets:
    - 0: no shedding,
    - 1: simpe,
    - 2: more physical. If applied together with the ART aerosol physics **inwp_gscp=6**, only options 0 and 1 are currently supported.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-D_shed_gh)=
    **D_shed_gh**
  - R
  - 0.009
  - m
  - Critical graupel/hail particle diameter for shedding during riming (wet growth) and melting. Shedding happens if:
    itype_shedding_gh = 1: {math}`D_{meanmass}>` D_shed_gh
    itype_shedding_gh = 2: in the spectral PSD-part where {math}`D>\max(D_{wetgr},\text{D\_shed\_gh})` - that is for wet growth but not below a stable diameter, e.g., 9 mm after Rasmussen and Heymfield
  - itype_shedding_gh=1,2 [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-llim_gr_prod_rain_riming)=
    **llim_gr_prod_rain_riming**
  - L
  - .FALSE.
  -
  - If .TRUE., limit the graupel production by rain riming of ice/snow by a bulk-density-based criterion on the mean-mass-particles of the collision partners.
  - [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-wgt_D_coll_limgrprod)=
    **wgt_D_coll_limgrprod**
  - R
  - 0.5
  -
  - Weight {math}`\in [0,1]` for the collided mean-mass-particle's diameter {math}`D_{coll}`: how much does the smaller collision partner contribute to the overall diameter?
  - llim_gr_prod_rain_riming=.TRUE. [iforcing](run_nml-iforcing)=2,3  inwp_gscp=4,5,6,7

* - (twomom_mcrph_nml-wgt_rho_coll_limgrprod)=
    **wgt_rho_coll_limgrprod**
  - R
  - 0.5
  -
  - Weight {math}`\in [0,1]` for the limit of the collided mean-mass-particle's bulk density: how near should it be to the bulk density of graupel in order to convert it to graupel?
  - llim_gr_prod_rain_riming=.TRUE.
    [iforcing](run_nml-iforcing)=2,3
    inwp_gscp=4,5,6,7
```



Defined and used in: {{ '[src/namelists/mo_2mom_mcrph_nml.f90]({}/src/namelists/mo_2mom_mcrph_nml.f90)'.format(base_url) }}

Internally these namelist parameters are stored in the container `atm_phy_nwp_config(jg)%cfg_2mom` of type `t_cfg_2mom`.

The defaults are defined in the container `cfg_2mom_default` in {{ '[src/atm_phy_schemes/mo_2mom_mcrph_config_default.f90]({}/src/atm_phy_schemes/mo_2mom_mcrph_config_default.f90)'.format(base_url) }}

Adding new parameters can easily be done along the lines of one of the above existing parameters.


(ref_buildrun_nml_octst_nml)=
## octst_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (octst_nml-h_val)=
    **h_val**
  -
  -
  -
  -
  -

* - (octst_nml-rlat_in)=
    **rlat_in**
  -
  -
  -
  -
  -

* - (octst_nml-rlon_in)=
    **rlon_in**
  -
  -
  -
  -
  -

* - (octst_nml-t_val)=
    **t_val**
  -
  -
  -
  -
  -
```


Defined and used in: {{ '[src/ocean/config/mo_ocean_nml.f90]({}/src/ocean/config/mo_ocean_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_oemctrl_nml)=
## oemctrl_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (oemctrl_nml-boundary_lambda_nc)=
    **boundary_lambda_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-boundary_regions_nc)=
    **boundary_regions_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-chem_init_nc)=
    **chem_init_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-chem_restart_nc)=
    **chem_restart_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-day_of_week_nc)=
    **day_of_week_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-ens_lambda_nc)=
    **ens_lambda_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-ens_reg_nc)=
    **ens_reg_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-gridded_emissions_nc)=
    **gridded_emissions_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-hour_of_day_nc)=
    **hour_of_day_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-hour_of_year_nc)=
    **hour_of_year_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-lat_cut_end)=
    **lat_cut_end**
  - REAL
  - 0.0
  -
  -
  -

* - (oemctrl_nml-lat_cut_start)=
    **lat_cut_start**
  - REAL
  - 0.0
  -
  -
  -

* - (oemctrl_nml-lcut_area)=
    **lcut_area**
  - LOGICAL
  - .FALSE.
  -
  -
  -

* - (oemctrl_nml-lon_cut_end)=
    **lon_cut_end**
  - REAL
  - 0.0
  -
  -
  -

* - (oemctrl_nml-lon_cut_start)=
    **lon_cut_start**
  - REAL
  - 0.0
  -
  -
  -

* - (oemctrl_nml-month_of_year_nc)=
    **month_of_year_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-restart_init_time)=
    **restart_init_time**
  - REAL
  - 0
  -
  -
  -

* - (oemctrl_nml-vegetation_indices_nc)=
    **vegetation_indices_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-vertical_profile_nc)=
    **vertical_profile_nc**
  - CHARACTER
  - ''
  -
  -
  -

* - (oemctrl_nml-vprm_alpha)=
    **vprm_alpha**
  -
  -
  -
  -
  -

* - (oemctrl_nml-vprm_beta)=
    **vprm_beta**
  -
  -
  -
  -
  -

* - (oemctrl_nml-vprm_lambda)=
    **vprm_lambda**
  -
  -
  -
  -
  -

* - (oemctrl_nml-vprm_par)=
    **vprm_par**
  -
  -
  -
  -
  -

* - (oemctrl_nml-vprm_tlow)=
    **vprm_tlow**
  -
  -
  -
  -
  -

* - (oemctrl_nml-vprm_tmax)=
    **vprm_tmax**
  -
  -
  -
  -
  -

* - (oemctrl_nml-vprm_tmin)=
    **vprm_tmin**
  -
  -
  -
  -
  -

* - (oemctrl_nml-vprm_topt)=
    **vprm_topt**
  -
  -
  -
  -
  -
```


Defined and used in: {{ '[src/namelists/mo_oem_nml.f90]({}/src/namelists/mo_oem_nml.f90)'.format(base_url) }}

(ref_buildrun_nml_output_nml)=
## output_nml

Relevant if [output](run_nml-output)='nml'

Please note: There may be several instances of
[output_nml](ref_buildrun_nml_output_nml) in the namelist file, every one defining a list of variables with
separate attributes for output.


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (output_nml-dom)=
    **dom**
  - I(:)
  - -1
  -
  - Array of domains for which this name-list is used. If not specified (or specified as -1 as the first array member), this name-list will be used for all domains. *Attention:* Depending on the setting of the parameter l_output_phys_patch these are either logical or physical domain numbers!
  -

* - (output_nml-file_interval)=
    **file_interval**
  - C
  - " "
  -
  - Defines the length of a file in terms of an ISO-8601 duration string. An example for this time stamp format is given below. This namelist parameter can be set instead of [steps_per_file](output_nml-steps_per_file).
  -

* - (output_nml-filename_format)=
    **filename_format**
  - C
  - see description.
  -
  - Output filename format. Includes keywords `path`, `output_filename`, `physdom`, etc. (see below). Default is `<output_filename>_DOM<physdom>_<levtype>_ <jfile>`
  -

* - (output_nml-filename_extn)=
    **filename_extn**
  - C
  - "default"
  -
  - User-specified filename extension (empty string also possible). If this namelist parameter is chosen as "default", then we have ".nc" for NetCDF output files, and ".grb" for GRIB2.
  -

* - (output_nml-filetype)=
    **filetype**
  - I
  - 4
  -
  - One of CDI's FILETYPE_XXX constants, or FILETYPE_NONE to prevent writing a file but adding the diagnostics (e.g. pressure level remap) to the variable list. Possible values:
    - 2: FILETYPE_GRB2
    - 4: FILETYPE_NC2
    - 5: FILETYPE_NC4
    - 999: FILETYPE_NONE
  -

* - (output_nml-filter_spec)=
    **filter_spec**
  - C
  - None
  -
  - Compression algorithm for all NetCDF-4 output fields. The required filter is defined via a unique identifier and parameters for controlling the action of the compression filter [(list of registered filters)](https://github.com/HDFGroup/hdf5_plugins/blob/master/docs/RegisteredFilterPlugins.md#list-of-filters-registered-with-the-hdf-group).  
  Here are two examples:  
  zstd with compression level 6: filter_spec='32015,6'  
  blosc with lz4 filter, compression level 5 and bitshuffle: filter_spec='32001,0,0,0,0,5,2,1'  
  Filters are available via HDF5 plugins. To be able to use them, the environment variable HDF5\_PLUGIN\_PATH must refer to a directory with an installation of these filters. Only filters installed in this directory are available. [More information about HDF5 filters](https://docs.unidata.ucar.edu/netcdf/NUG/md_filters.html).
  - filetype=5

* - (output_nml-chunk_size)=
    **chunk_size**
  - I
  - CDI default
  -
  - Chunk size for all NetCDF-4 output fields. chunk_size defines the number of data elements for a chunk. A chunk is a data block that is used for a compression unit. The size of a chunk depends on the size of the horizontal grid. For 32-bit floating point data, the optimum chunk size is between 32768 and 4194304. The default value in CDI is a maximum of 1048576. chunk_size is only used in combination with the filter_spec parameter.
  - filetype=5

* - (output_nml-number_of_bits)=
    **number_of_bits**
  - I
  - None
  -
  - Number of significant bits for all NetCDF-4 output fields. This leads to a loss of information and to higher compression rates. number_of_bits is only useful in combination with the filter_spec parameter.
  - filetype=5

* - (output_nml-m_levels)=
    **m_levels**
  - C
  - None
  -
  - Model level indices (optional). Allowed is a comma- (or semicolon-) separated list of integers, and of integer ranges like "10...20".  One may also use the keyword "nlev" to denote the maximum integer (or, equivalently, "n" or "N"). Furthermore, arithmetic expressions like "(nlev - 2)" are possible. Basic example: `m_levels = "1,3,5...10,20...(nlev-2)"`
  -

* - (output_nml-h_levels)=
    **h_levels**
  - R(:)
  - None
  - m
  - height levels (height above mean sea level, _not_ above ground)
  -

* - (output_nml-p_levels)=
    **p_levels**
  - R(:)
  - None
  - Pa
  - pressure levels
  -

* - (output_nml-i_levels)=
    **i_levels**
  - R(:)
  - None
  - K
  - isentropic levels
  -

* - (output_nml-ml_varlist)=
    **ml_varlist**
  - C(:)
  - None
  -
  - Name of model level fields to be output.
  -

* - (output_nml-hl_varlist)=
    **hl_varlist**
  - C(:)
  - None
  -
  - Name of height level fields to be output.
  -

* - (output_nml-pl_varlist)=
    **pl_varlist**
  - C(:)
  - None
  -
  - Name of pressure level fields to be output.
  -

* - (output_nml-il_varlist)=
    **il_varlist**
  - C(:)
  - None
  -
  - Name of isentropic level fields to be output.
  - itype_pres_msl < 3 (for technical reasons?)

* - (output_nml-include_last)=
    **include_last**
  - L
  - .TRUE.
  -
  - Flag whether to include the last time step
  -

* - (output_nml-mode)=
    **mode**
  - I
  - 2
  -
  - 1: forecast mode
    - 2: climate mode
    In climate mode the time axis of the output file is set to TAXIS_ABSOLUTE. In forecast mode it is set to TAXIS_RELATIVE. Till now the forecast mode only works if the output is at multiples of 1 hour
  -

* - (output_nml-taxis_tunit)=
    **taxis_tunit**
  - I
  - 2
  -
  - Time unit of the TAXIS_RELATIVE time axis.
    - 1: TUNIT_SECOND
    - 2: TUNIT_MINUTE
    - 5: TUNIT_HOUR
    - 9: TUNIT_DAY For a complete list of possible values see cdilib.c
  - mode=1

* - (output_nml-output_bounds)=
    **output_bounds**
  - R({math}`k \ast` 3)
  - None
  -
  - Post-processing times: start, end, increment. The increment (output interval) must be larger than the advection time step (`dtime`) and should be an integer multiple of it. Multiple triples are possible in order to define multiple starts/ends/intervals. See namelist parameters `output_start`, `output_end`, `output_interval` for an alternative specification of output events.
  -

* - (output_nml-output_time_unit)=
    **output_time_unit**
  - I
  - 1
  -
  - Units of output bounds specification:
    - 1: second
    - 2: minute
    - 3: hour
    - 4: day
    - 5: month
    - 6: year
  -

* - (output_nml-output_filename)=
    **output_filename**
  - C
  - None
  -
  - Output filename prefix (which may include path). Domain number, level type, file number and extension will be added, according to the format given in namelist parameter "filename_format".
  -

* - (output_nml-output_grid)=
    **output_grid**
  - L
  - .FALSE.
  -
  - Flag whether grid information is added to output.
  -

* - (output_nml-output_start)=
    **output_start**
  - C(:)
  - " "
  -
  - ISO8601 time stamp for begin of output. An example for this time stamp format is given below. More than one value is possible in order to define multiple start/end/interval triples. See namelist parameter `output_bounds` for an alternative specification of output events.
  -

* - (output_nml-output_end)=
    **output_end**
  - C(:)
  - " "
  -
  - ISO8601 time stamp for end of output. An example for this time stamp format is given below. More than one value is possible in order to define multiple start/end/interval triples. See namelist parameter `output_bounds` for an alternative specification of output events.
  -

* - (output_nml-output_interval)=
    **output_interval**
  - C(:)
  - " "
  -
  - ISO8601 time stamp for repeating output intervals. The output interval must be larger than the advection time step (`dtime`) and should be an integer multiple of it. An example for this time stamp format is given below. More than one value is possible in order to define multiple start/end/interval triples. See namelist parameter `output_bounds` for an alternative specification of output events.
  -

* - (output_nml-operation)=
    **operation**
  - C
  - None
  -
  - Use this variable for internal diagnostics applied on all given output variables or groups except time-constant ones: `mean` for generating time averaged, `square` for time averaged square values, `max` or `min` for maximum and minimum and `acc` for accumulated values within the corresponding interval, i.e.  `output_interval`. Supported are 2D, 3D and single values like global means on model levels of all components. All operations can be used on global and nested grids.
  -

* - (output_nml-pe_placement_il)=
    **pe_placement_il**
  - I(:)
  - -1
  -
  - Advanced output option: Explicit assignment of output MPI ranks to the isentropic level output file. At most [stream_partitions_il](output_nml-stream_partitions_il) different ranks can be specified. See namelist parameter `pe_placement_ml` for further details.
  -

* - (output_nml-pe_placement_hl)=
    **pe_placement_hl**
  - I(:)
  - -1
  -
  - Advanced output option: Explicit assignment of output MPI ranks to the height level output file. At most [stream_partitions_hl](output_nml-stream_partitions_hl) different ranks can be specified. See namelist parameter `pe_placement_ml` for further details.
  -

* - (output_nml-pe_placement_ml)=
    **pe_placement_ml**
  - I(:)
  - -1
  -
  - Advanced output option: Explicit assignment of output MPI ranks to the model level output file. At most [stream_partitions_ml](output_nml-stream_partitions_ml) different ranks can be specified, out of the following list: `0` {math}`\ldots` (`num_io_procs` - 1). If this namelist parameters is not provided, then the output ranks are chosen in a Round-Robin fashion among those ranks that are not occupied by explicitly placed output files.
  -

* - (output_nml-pe_placement_pl)=
    **pe_placement_pl**
  - I(:)
  - -1
  -
  - Advanced output option: Explicit assignment of output MPI ranks to the pressure level output file. At most [stream_partitions_pl](output_nml-stream_partitions_pl) different ranks can be specified. See namelist parameter `pe_placement_ml` for further details.
  -

* - (output_nml-ready_file)=
    **ready_file**
  - C
  - "default"
  -
  - A _ready file_ is a technique for handling dependencies between the NWP processes. The completion of the write process is signalled by creating a small file with name `ready_file`. Different [output_nml](ref_buildrun_nml_output_nml) may be joined together to form a single ready file event. The setting of `ready_file="default"` does not create a ready file. The ready file name may contain string tokens `<path>`, `<datetime>`, `<ddhhmmss>`, `<datetime2>` which are substituted as described for the namelist parameter `filename_format`.
  -

* - (output_nml-reg_def_mode)=
    **reg_def_mode**
  - I
  - 0
  -
  - Specify if the "delta" value prescribes an interval size or the total *number* of intervals:  
    - 0: switch automatically between increment and no. of grid points  
    - 1: `reg_lon/lat_def(2)` specifies increment  
    - 2: `reg_lon/lat_def(2)` specifies no. of grid points
  - [remap](output_nml-remap) = 1

* - (output_nml-remap)=
    **remap**
  - I
  - 0
  -
  - Interpolate horizontally:  
    - 0: none  
    - 1: to regular lat-lon grid
  -

* - (output_nml-north_pole)=
    **north_pole**
  - R(2)
  - 0,90
  -
  - Definition of north pole for rotated lon-lat grids (`[longitude, latitude]`).
  -

* - (output_nml-reg_lat_def)=
    **reg_lat_def**
  - R(3)
  - None
  -
  - Start, increment, end latitude in degrees. Alternatively, the user may set the number of grid points instead of an increment. Details for the setting of regular grids are given below together with an example.
  - [remap](output_nml-remap) = 1

* - (output_nml-reg_lon_def)=
    **reg_lon_def**
  - R(3)
  - None
  -
  - The regular grid points are specified by three values: start, increment, end given in degrees. Alternatively, the user may set the number of grid points instead of an increment. Details for the setting of regular grids are given below together with an example.
  - [remap](output_nml-remap) = 1

* - (output_nml-steps_per_file)=
    **steps_per_file**
  - I
  - -1
  -
  - Max number of output steps in one output file. If this number is reached, a new output file will be opened. Setting [steps_per_file](output_nml-steps_per_file) to 1 enforces a flush when writing is completed, so that the file is immediately accessible for reading.
  -

* - (output_nml-steps_per_file_inclfirst)=
    **steps_per_file_inclfirst**
  - L
  - see descr.
  -
  - Defines if first step is counted with respect to [steps_per_file](output_nml-steps_per_file) file count. The default is `.FALSE.` for GRIB2 output, and `.TRUE.` otherwise.
  -

* - (output_nml-stream_partitions_hl)=
    **stream_partitions_hl**
  - I
  - 1
  -
  - Splits height level output of this namelist into several concurrent alternating files. See namelist parameter [stream_partitions_ml](output_nml-stream_partitions_ml) for details.
  -

* - (output_nml-stream_partitions_il)=
    **stream_partitions_il**
  - I
  - 1
  -
  - Splits isentropic level output of this namelist into several concurrent alternating files. See namelist parameter [stream_partitions_ml](output_nml-stream_partitions_ml) for details.
  -

* - (output_nml-stream_partitions_ml)=
    **stream_partitions_ml**
  - I
  - 1
  -
  - Splits model level output of this namelist into several concurrent alternating files. The output is split into {math}`N` files, where the start date of part {math}`i` gets an offset of {math}`(i-1)*``output_interval`. The output interval is then replaced by {math}`N*``output_interval`, the `include_last` flag is set to `.FALSE.`, the [steps_per_file_inclfirst](output_nml-steps_per_file_inclfirst) flag is set to `.FALSE.`, and the [steps_per_file](output_nml-steps_per_file) counter is set to `1`.
  -

* - (output_nml-stream_partitions_pl)=
    **stream_partitions_pl**
  - I
  - 1
  -
  - Splits pressure level output of this namelist into several concurrent alternating files. See namelist parameter [stream_partitions_ml](output_nml-stream_partitions_ml) for details.
  -

* - (output_nml-rbf_scale)=
    **rbf_scale**
  - R
  - -1.
  -
  - Explicit setting of RBF shape parameter for interpolated lon-lat output. This namelist parameter is only active in combination with `rbf_scale_mode_ll = 3`.
  - rbf_scale_mode_ll = 3
```

Defined and used in: {{ '[src/io/shared/mo_name_list_output_init.f90]({}/src/io/shared/mo_name_list_output_init.f90)'.format(base_url) }}


### Interpolation onto regular grids


Horizontal interpolation onto regular grids is possible through the namelist setting  [remap](output_nml-remap)=1, where
the mesh is defined by the parameters

 - [reg_lon_def](output_nml-reg_lon_def): mesh latitudes in degrees,
 - [reg_lat_def](output_nml-reg_lat_def): mesh longitudes in degrees,
 - [north_pole](output_nml-north_pole): definition of north pole for rotated lon-lat grids.


The regular grid points in [reg_lon_def](output_nml-reg_lon_def), [reg_lat_def](output_nml-reg_lat_def) are each specified by three values, given in degrees:
_start_, _increment_, _end_.
The mesh then contains all grid points {math}`start + k * increment <= end`, where {math}`k` is an integer.
Instead of defining an increment it is also possible to prescribe the number of grid points.

 - Setting the namelist parameter [reg_def_mode](output_nml-reg_def_mode)=0:
   Switch automatically from increment specification to no. of grid points,
   when the  `reg_lon/lat_def(2)` value is larger than 5.0.
 - 1: `reg_lon/lat_def(2)` specifies increment
 - 2: `reg_lon/lat_def(2)` specifies no. of grid points

For longitude values the last grid point is omitted if the end point matches the start point, e.g. for 0 and 360 degrees.


**Examples**

```{list-table}
:header-rows: 1

* - Description
  - Parameter

* - local grid with 0.5 degree increment:
  - [reg_lon_def](output_nml-reg_lon_def) = -30.,0.5,30.

* -
  - [reg_lat_def](output_nml-reg_lat_def) = 90.,-0.5, -90.

* - global grid with 720x361 grid points:
  - [reg_lon_def](output_nml-reg_lon_def) = 0.,720,360.

* -
  - [reg_lat_def](output_nml-reg_lat_def) = -90.,360,90.
```



### Time stamp format


The namelist parameters `output_start`, `output_end`, `output_interval` allow
the specification of time stamps according to ISO 8601.
The general format for time stamps is `YYYY-MM-DDThh:mm:ss`
where `Y`: year, `M`: month, `D`: day for dates,
and   `hh`: hour, `mm`: minute, `ss`: second for time strings.
The general format for durations is `PnYnMnDTnHnMnS`.
See, for example, <https://en.wikipedia.org/wiki/ISO_8601> for details and further specifications.

**NOTE:** as the mtime library underlaying the output driver
currently has some restrictions concerning the specification of durations:

1. Any number `n` in `PnYnMnDTnHnMnS` must have two digits. For instance use `"PT06H"` instead of `"PT6H"`
2. In a duration string `PnyearYnmonMndayDTnhrHnminMnsecS` the numbers `nxyz` must not pass the carry over number to the next larger time unit: 0<=nmon<=12, 0<=nhr<=23, 0<=nmin<=59, 0<=nsec<=59.999. For instance use `"P01D"` instead of `"PT24H"`, or `"PT01M"` instead of `"PT60S"`.


Soon the formatting problem will be resolved and the valid number ranges will be enlarged.
(2013-12-16).

**Examples**

```{list-table}
:header-rows: 1

* - Description
  - Parameter

* - date and time representation (`output_start`, `output_end`)
  - **2013-10-27T13:41:00Z**

* - duration (`output_interval`)
  - **P00DT06H00M00S**
```



### Variable Groups


#### Keyword `group`

Using the `group:` keyword for the namelist parameters `ml_varlist`, `hl_varlist`, `pl_varlist`,
sets of common variables can be added to the output:

```{list-table}
:header-rows: 1

* - Parameter
  - Description

* - (group-all)=
    **all**
  - output of all variables (caution: do not combine with **mixed** vertical interpolation)

* - (group-atmo_ml_vars)=
    **atmo_ml_vars**
  - basic atmospheric variables on model levels

* - (group-atmo_pl_vars)=
    **atmo_pl_vars**
  - same set as atmo_ml_vars, but except pres

* - (group-atmo_zl_vars)=
    **atmo_zl_vars**
  - same set as atmo_ml_vars, but expect height

* - (group-nh_prog_vars)=
    **nh_prog_vars**
  - additional prognostic variables of the nonhydrostatic model

* - (group-atmo_derived_vars)=
    **atmo_derived_vars**
  - derived atmospheric variables

* - (group-rad_vars)=
    **rad_vars**
  -

* - (group-precip_vars)=
    **precip_vars**
  -

* - (group-cloud_diag)=
    **cloud_diag**
  -

* - (group-pbl_vars)=
    **pbl_vars**
  -

* - (group-phys_tendencies)=
    **phys_tendencies**
  -

* - (group-land_vars)=
    **land_vars**
  -

* - (group-snow_vars)=
    **snow_vars**
  - snow variables

* - (group-multisnow_vars)=
    **multisnow_vars**
  - multi-layer snow variables

* - (group-additional_precip_vars)=
    **additional_precip_vars**
  -

* - (group-dwd_fg_atm_vars)=
    **dwd_fg_atm_vars**
  - DWD first guess fields (atmosphere)

* - (group-dwd_fg_sfc_vars)=
    **dwd_fg_sfc_vars**
  - DWD first guess fields (surface/soil)

* - (group-ART_AERO_VOLC)=
    **ART_AERO_VOLC**
  - ART volcanic ash fields

* - (group-ART_AERO_RADIO)=
    **ART_AERO_RADIO**
  - ART radioactive tracer fields

* - (group-ART_AERO_DUST)=
    **ART_AERO_DUST**
  - ART mineral dust aerosol fields

* - (group-ART_AERO_SEAS)=
    **ART_AERO_SEAS**
  - ART sea salt aerosol fields

* - (group-prog_timemean)=
    **prog_timemean**
  - time mean output: temp, u, v, rho

* - (group-tracer_timemean)=
    **tracer_timemean**
  - time mean output: qv, qc, qi

* - (group-atmo_timemean)=
    **atmo_timemean**
  - time mean variables from `prog_timemean`,`tracer_timemean`
```




#### Keyword `tiles`

The `"tiles:"` keyword allows to add all tiles of a specific variable to the output, without the need to specify all
tile fields separately. E.g. `tiles:t_g` (read: "tiles of `t_g`") automatically adds all `t_g_t_X`
fields to the output. Here, `X` is a placeholder for the tile number. Make sure to specify the name of the aggregated
variable rather than the name of the corresponding tile container (i.e. in the given example it must be `t_g`, and not
`t_g_t`!).

_Note:_

There exists a special syntax which allows to remove variables from the output list, e.g. if
these undesired variables were contained in a previously selected group.<br>
Typing `-<varname>` (for example `-temp`) removes the
variable from the union set of group variables and other selected variables.
Note that typos are not detected but that the corresponding variable is
simply not removed!



#### Keyword substitution in output filename (`filename_format`)


```{list-table}
:header-rows: 1

* - Parameter
  - Description

* - (filename_format-path)=
    **path**
  - substituted by `model_base_dir`

* - (filename_format-output_filename)=
    **output_filename**
  - substituted by `output_filename`

* - (filename_format-physdom)=
    **physdom**
  - substituted by physical patch ID

* - (filename_format-levtype)=
    **levtype**
  - substituted by level type "ML", "PL", "HL", "IL"

* - (filename_format-levtype_l)=
    **levtype_l**
  - like `levtype`, but in lower case

* - (filename_format-jfile)=
    **jfile**
  - substituted by output file counter

* - (filename_format-datetime)=
    **datetime**
  - substituted by ISO-8601 date-time stamp in format `YYYY-MM-DDThh:mm:ss.sssZ`

* - (filename_format-datetime2)=
    **datetime2**
  - substituted by ISO-8601 date-time stamp in format `YYYYMMDDThhmmssZ`

* - (filename_format-datetime3)=
    **datetime3**
  - substituted by ISO-8601 date-time stamp in format `YYYYMMDDThhmmss.sssZ`

* - (filename_format-ddhhmmss)=
    **ddhhmmss**
  - substituted by _relative_ day-hour-minute-second string

* - (filename_format-dddhhmmss)=
    **dddhhmmss**
  - substituted by _relative_ three-digit day-hour-minute-second string

* - (filename_format-hhhmmss)=
    **hhhmmss**
  - substituted by _relative_ hour-minute-second string

* - (filename_format-npartitions)=
    **npartitions**
  - If namelist is split into concurrent files: number of stream partitions.

* - (filename_format-ifile_partition)=
    **ifile_partition**
  - If namelist is split into concurrent files: stream partition index of this file.

* - (filename_format-total_index)=
    **total_index**
  - If namelist is split into concurrent files: substituted by the file counter (like in `jfile`), which an "unsplit" namelist would have produced.

```




(ref_buildrun_nml_parallel_nml)=
## parallel_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (parallel_nml-nproma)=
    **nproma**
  - I
  - 1
  -
  - Loop chunk length. Only one of (_nproma_, _nblocks_c_, _nblocks_e_) may be specified in the namelist ({math}`> 0`) at any time.
  -

* - (parallel_nml-nblocks_c)=
    **nblocks_c**
  - I
  - 0
  -
  - Number of looping chunks used for cells. For values > 0, _nproma_ is recomputed according to the specified _nblocks_c_.
  -

* - (parallel_nml-nblocks_e)=
    **nblocks_e**
  - I
  - 0
  -
  - Number of looping chunks used for edges. For values > 0, _nproma_ is recomputed according to the specified _nblocks_e_.
  -

* - (parallel_nml-nproma_sub)=
    **nproma_sub**
  - I
  - nproma
  -
  - Chunk size of subblocks used for example by ecRad or rrtmgp, which is needed for the GPU port to reduce the memory footprint. May only specify one of (_nproma_sub_, _nblocks_sub_) in the namelist ({math}`> 0`) at any time.
  -

* - (parallel_nml-nblocks_sub)=
    **nblocks_sub**
  - I
  - 1
  -
  - Number of looping chunks used for subblocking. For values {math}`<=` 0 this is ignored. For bigger values, this overwrites _nproma_sub_. For reduced-grid radiation, we suggest explicitly specifying _nproma_sub_ instead of using _nblocks_sub_.
  -

* - (parallel_nml-n_ghost_rows)=
    **n_ghost_rows**
  - I
  - 1
  -
  - number of halo cell rows
  -

* - (parallel_nml-division_method)=
    **division_method**
  - I
  - 1
  -
  - Method of domain decomposition:
    - 0: read in from file
    - 1: use built-in geometric subdivision
  -

* - (parallel_nml-division_file_name)=
    **division_file_name**
  - C
  -
  -
  - Name of division file
  - division_method = 0

* - (parallel_nml-ldiv_phys_dom)=
    **ldiv_phys_dom**
  - L
  - .TRUE.
  -
  - .TRUE.: split into physical domains before computing domain decomposition (in case of merged domains)
    (This reduces load imbalance; turning off this option is not recommended except for very small processor numbers)
  - division_method = 1

* - (parallel_nml-p_test_run)=
    **p_test_run**
  - L
  - .FALSE.
  -
  - .TRUE. means verification run for MPI parallelization (PE 0 processes full domain)
  -

* - (parallel_nml-num_test_pe)=
    **num_test_pe**
  - I
  - -1
  -
  - If set to more than 1, use this many ranks for testing and switch to different consistency test. This enables tests for identity in setups which are too big to run on a single rank but is limited to comparing one MPI parallelization setup vs. another, obviously.
  - p_test_run = .TRUE.

* - (parallel_nml-l_test_openmp)=
    **l_test_openmp**
  - L
  - .FALSE.
  -
  - if .TRUE. is combined with p_test_run=.TRUE. and OpenMP parallelization, the test PE gets only 1 thread in order to verify the OpenMP parallelization
  - p_test_run = .TRUE.

* - (parallel_nml-l_log_checks)=
    **l_log_checks**
  - L
  - .FALSE.
  -
  - if .TRUE. messages are generated during each synchonization step (use for debugging only)
  -

* - (parallel_nml-l_fast_sum)=
    **l_fast_sum**
  - L
  - .FALSE.
  -
  - if .TRUE., use fast (not processor-configuration-invariant) global summation
  -

* - (parallel_nml-use_dycore_barrier)=
    **use_dycore_barrier**
  - L
  - .FALSE.
  -
  - if .TRUE., set an MPI barrier at the beginning of the nonhydrostatic solver (do not use for production runs!)
  -

* - (parallel_nml-itype_exch_barrier)=
    **itype_exch_barrier**
  - I
  - 0
  -
  - 1: set an MPI barrier at the beginning of each MPI exchange call
    - 2: set an MPI barrier after each MPI WAIT call
    - 3: 1+2 (do not use for production runs!)
  -

* - (parallel_nml-iorder_sendrecv)=
    **iorder_sendrecv**
  - I
  - 1
  -
  - Sequence of send/receive calls:
    - 1: irecv/send
    - 2: isend/recv
    - 3: isend/irecv
  -

* - **default_comm-
    _pattern_type**
  - I
  - 1
  -
  - Default implementation of mo_communication to be used:
    - 1: original
    - 2: YAXT
  -

* - (parallel_nml-num_io_procs)=
    **num_io_procs**
  - I
  - 0
  -
  - Number of I/O processors (running exclusively for doing I/O)
  -

* - (parallel_nml-num_io_procs_radar)=
    **num_io_procs_radar**
  - I
  - 0
  -
  - Number of dedicated I/O processors for the efficient radar forward operator EMVORADO. Choosing more I/O processors than the total number of simulated radar stations of all domains is not advisable, because one station is handled by one I/O processor. However, less I/O processors can be chosen, in which case one processor handles several stations. I/O tasks actually include much more than plain output for each station and can be very time consuming. More details can be found in the EMVORADO User's Guide, available from the COSMO web page (<https://www.cosmo-model.org> {math}`\rightarrow` Documentation {math}`\rightarrow` EMVORADO) or from the emvorado submodule `externals/emvorado/DOC/TEX/emvorado_userguide.pdf`. If `num_io_procs_radar`=0, a subset of the worker processors (=number of radar stations) are doing the I/O tasks, which may slow down the model considerably.
  - luse_radarfwo(\<idom>) =.TRUE., [iforcing](run_nml-iforcing)=3

* - (parallel_nml-num_restart_procs)=
    **num_restart_procs**
  - I
  - 0
  -
  - Number of restart processors (running exclusively for doing restart)
  -

* - (parallel_nml-num_prefetch_proc)=
    **num_prefetch_proc**
  - I
  - 1
  -
  - Number of processors for prefetching of boundary data asynchronously for a limited area run (running exclusively for reading Input boundary data. Maximum no of processors used for it is limited to 1).
  - Mandatory for itype_latbc = 1

* - (parallel_nml-proc0_shift)=
    **proc0_shift**
  - I
  - 0
  -
  - Number of processors at the beginning of the rank list that are excluded from the domain decomposition. Setting this parameter to 1 serves for offloading I/O to the vector hosts of the NEC Aurora, but it works technically on other platforms as well.
  -

* - (parallel_nml-use_omp_input)=
    **use_omp_input**
  - L
  - .FALSE.
  -
  - Setting this parameter to .TRUE. activates OpenMP sections in initicon that allow task parallelism for reading atmospheric input data, overlapping reading, sending, and statistics calculations.
  -

* - (parallel_nml-pio_type)=
    **pio_type**
  - I
  - 1
  -
  - Type of parallel I/O.
    - 1: Classical async I/O proccessors
    - 2: CDI-PIO (Experimental!)
  -

* - (parallel_nml-use_icon_comm)=
    **use_icon_comm**
  - L
  - .FALSE.
  -
  - Enable the use of MPI bulk communication through the icon_comm_lib
  -

* - (parallel_nml-icon_comm_debug)=
    **icon_comm_debug**
  - L
  - .FALSE.
  -
  - Enable debug mode for the icon_comm_lib
  -

* - **max_send_recv-
    _buffer_size**
  - I
  - 131072
  -
  - Size of the send/receive buffers for the icon_comm_lib.
  -

* - (parallel_nml-use_dp_mpi2io)=
    **use_dp_mpi2io**
  - L
  - .FALSE.
  -
  - Enable this flag if output fields shall be gathered by the output processes in DOUBLE PRECISION.
  -

* - (parallel_nml-restart_chunk_size)=
    **restart_chunk_size**
  - I
  - 1
  -
  - _(Advanced namelist parameter:)_ Number of levels to be buffered by the asynchronous restart process. The (asynchronous) restart is capable of writing and communicating more than one 2D slice at once.
  -

* - (parallel_nml-num_dist_array_replicas)=
    **num_dist_array_replicas**
  - I
  - 1
  -
  - _(Advanced namelist parameter:)_ Number of replicas of the distributed array used for the pre_patch.
  -

* - (parallel_nml-io_process_stride)=
    **io_process_stride**
  - I
  - -1
  -
  - _(Advanced namelist parameter:)_ Stride of processes taking part in reading of data. (Few reading processes, i.e. a large stride, often gives best performance.)
  -

* - (parallel_nml-process_stride_pgrib I)=
    **process_stride_pgrib I**
  - 1
  -
  -
  - _(Advanced namelist parameter:)_ Stride of processes taking part in parallel GRIB decoding (recommendation: select stride such as to have between 100 and 1000 decoder PEs), see parallel_grib_decoding
  -

* - (parallel_nml-io_process_rotate)=
    **io_process_rotate**
  - I
  - 0
  -
  - _(Advanced namelist parameter:)_ Rotate of processes taking part in reading of data. (Process taking part if p_pe_work % stride == rotate)
  -
```



Defined and used in: {{ '[src/namelists/mo_parallel_nml.f90]({}/src/namelists/mo_parallel_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_radiation_nml)=
## radiation_nml

Relevant if [iforcing](run_nml-iforcing)=3 (NWP)

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (radiation_nml-isolrad)=
    **isolrad**
  - I
  - 1
  -
  - Insolation scheme
    - 0: Use original insolation (from SRTM in case inwp_radiation=1 or from ecRad in case inwp_radiation=4)
    - 1: Use SSI values from Coddington et al. (2016) (inwp_radiation=1) or scale SSI values to Coddington et al. (2016) values (inwp_radiation=4)
    - 2: SSI from an external file containing monthly mean time series (inwp_radiation=4)
    - 4: SSI from AMIP (inwp_radiation=4)
  -

* - (radiation_nml-izenith)=
    **izenith**
  - I
  - 4
  -
  - Choice of zenith angle formula for the radiative transfer computation.
    - 0: Sun in zenith everywhere
    - 1: Zenith angle depends only on latitude
    - 2: Zenith angle depends only on latitude. Local time of day fixed at 07:14:15 for radiative transfer computation (sin(time of day) = 1/pi
    - 3: Zenith angle changing with latitude and time of day
    - 4: Zenith angle and irradiance changing with season, latitude, and time of day ([iforcing](run_nml-iforcing)=inwp only)
    - 5: Zenith angel for radiative convective equilibrium test: perpetual equinox with 340 W/m2
    - 6: Zenith angle with prescribed cosine of solar zenith angle (see parameter cos_zenith_fixed)
  -

* - (radiation_nml-cos_zenith_fixed)=
    **cos_zenith_fixed**
  - R
  - 0.5
  -
  - Cosine of zenith angle for test cases including SCM
  - izenith=6

* - (radiation_nml-islope_rad(max_dom))=
    **islope_rad(max_dom)**
  - I
  - 0
  -
  - Slope correction for surface radiation:
    - 0: None
    - 1: Slope correction for direct solar radiation without shading effects
    - 2: Removed
    - 3: Slope correction for direct solar radiation including shading, but no further consideration of sky-view factor effects.
  - 3 requires the additional field HORIZON to be present in the extpar data

* - (radiation_nml-albedo_type)=
    **albedo_type**
  - I
  - 1
  -
  - Type of surface albedo
    - 1: based on soil type specific tabulated values (dry soil)
    - 2: MODIS albedo
    - 3: fixed albedo for SCM and other testcases
  - [iforcing](run_nml-iforcing)=inwp

* - (radiation_nml-albedo_fixed)=
    **albedo_fixed**
  - R
  - 0.5
  -
  - Fixed albedo value for SCM and various testcases
  - [iforcing](run_nml-iforcing)=inwp
    albedo_type=3

* - (radiation_nml-direct_albedo)=
    **direct_albedo**
  - I
  - 4
  -
  - Direct beam surface albedo over land and sea-ice. Options mainly differ in terms of their solar zenith angle (SZA) dependency.
    - 1: Ritter-Geleyn (1992)
    - 2: Zaengl (pers. comm.): For 'rough surfaces' over land direct albedo is not allowed to exceed the corresponding broadband diffuse albedo. Ritter-Geleyn for ice.
    - 3: Yang et al (2008) for snow-free land points. Ritter-Geleyn for ice and Zaengl for snow.
    - 4: Briegleb and Ramanathan (1992) for snow-free land points. Ritter-Geleyn for ice and Zaengl for snow.
  - [iforcing](run_nml-iforcing)=inwp
    albedo_type=2

* - (radiation_nml-direct_albedo_water)=
    **direct_albedo_water**
  - I
  - 2
  -
  - Direct beam surface albedo over water (ocean or lake). Options mainly differ in terms of their solar zenith angle (SZA) dependency.
    - 1: Ritter-Geleyn (1992)
    - 2: Yang (2008), originally designed for land
    - 3: Taylor et al (1996) for direct and 0.06 for diffuse albedo as in the IFS.
  - [iforcing](run_nml-iforcing)=inwp
    albedo_type=2

* - (radiation_nml-albedo_whitecap)=
    **albedo_whitecap**
  - I
  - 0
  -
  - Ocean albedo increase by foam from breaking waves (whitecaps). Not applied over lakes.
    - 0: off
    - 1: whitecap describtion by Seferian et al 2018
  - [iforcing](run_nml-iforcing)=inwp
    albedo_type=2

* - (radiation_nml-icld_overlap)=
    **icld_overlap**
  - I
  - 2
  -
  - Method for cloud overlap calculation in shortwave part of RRTM
    - 1: maximum-random overlap
    - 2: generalized overlap (Hogan, Illingworth, 2000)
    - 3: maximum overlap
    - 4: random overlap
    - 5: exponential overlap
  - [iforcing](run_nml-iforcing)=inwp
    inwp_radiation=1 (1-4)
    inwp_radiation=4 (1,2,5)

* - (radiation_nml-irad_)=
    **irad_<gasname>**
  - I
  - 1
    2
    3
    3
    0
    2
    2
    2
  -
  - Switches for the concentration of radiative agents
    irad_xyz = -1: externally specified, e.g. by ComIn. This creates an additional mass mixing ratio field xyzrad_ext (e.g. co2rad_ext)
    irad_xyz = 0: set to zero
    irad_h2o = 1: vapor, cloud water and cloud ice from tracer variables
    irad_co2 = 1: CO{math}`_\text{2}` from tracer variable
    irad_co2/ch4/n2o/o2/cfc11/cfc12 = 2: concentration given by vmr_co2/ch4/n2o/o2/cfc11/cfc12
    irad_ch4/n2o = 3: tanh-profile with surface concentration given by vmr_ch4/n2o
    irad_co2/cfc11/cfc12 = 4: time dependent concentration from greenhouse gas file
    irad_ch4/n2o = 4: time dependent tanh-profile with surface concentration from greenhouse gas file
    irad_o3 = 2: ozone climatology from MPI
    irad_o3 = 4: ozone clim for Aqua Planet Exp
    irad_o3 = 5: 3-dim concentration, time dependent, monthly means from yearly files `bc_ozone_<year>.nc` or - with nesting - `bc_ozone_DOM<jg>_<year>.nc`
    irad_o3 = 6: ozone climatology with T5 geographical distribution and Fourier series for seasonal cycle for [iforcing](run_nml-iforcing) = 3 (NWP)
    irad_o3 = 7: GEMS ozone climatology (from IFS) for [iforcing](run_nml-iforcing) = 3 (NWP)
    irad_o3 = 9: MACC ozone climatology (from IFS) for [iforcing](run_nml-iforcing) = 3 (NWP)
    irad_o3 = 79: Blending between GEMS and MACC ozone climatologies (from IFS) for [iforcing](run_nml-iforcing) = 3 (NWP); MACC is used over Antarctica
    irad_o3 = 97: As 79, but MACC is also used above 1 hPa with transition zone between 5 hPa and 1 hPa
    irad_o3 = 10: Linearized ozone chemistry (ART extension necessary) for [iforcing](run_nml-iforcing) = 3 (NWP)
    irad_o3 = 11: Ozone from SCM input file
  -

* - (radiation_nml-vmr_)=
    **vmr_<agent>**
  - R
  - co2: 348.0e-6
    ch4: 1650.0e-9
    n2o: 306.0e-9
    o2: 0.20946
    cfc11: 214.5e-12
    cfc12: 371.1e-12
  -
  - Volume mixing ratio of the radiative agents
  -

* - (radiation_nml-irad_aero)=
    **irad_aero**
  - I
  - 2
  -
  - Aerosols
    - 0: none
    - 2: global constant
    - 3: externally specified. This creates additional 4D-fields for: optical depth long wave (od_lw), optical depth short wave (od_sw), single scattering albedo short wave (ssa_sw) and asymmetry parameter short wave (g_sw). These fields need to be filled externally (e.g. by a ComIn plugin).
    - 6: Tegen aerosol climatology for [iforcing](run_nml-iforcing) = 3 (NWP) .AND. itopo =1
    - 7: CAMS 3D aerosol climatology, the filename can be specified via cams_aero_filename
    - 8: CAMS 3D forecasted aerosol, the filename can be specified via cams_aero_filename
    - 9: ART online aerosol radiation interaction, uses Tegen for aerosols not chosen to be represented in ART for [iforcing](run_nml-iforcing) = 3 (NWP) .AND. itopo =1 .AND. lart=TRUE .AND. iart
    _ari=1
    - 12: tropospheric 'Kinne' aerosols, constant in time
    - 13: total tropospheric 'Kinne' aerosols, time dependent from file
    - 14: volcanic stratospheric aerosols, time dependent from file
    - 15: tropospheric 'Kinne' aerosols + volcanic stratospheric aerosols, time dependent, both from file. If the 1850--file of the 'Kinne' aerosols is given, only the natural background from Kinne aerosol is applied.
    - 18: tropospheric natural 'Kinne' aerosols from pre-industry (the 1850--file has to be linked for all simulations!) + time dep. volcanic stratospheric aerosols, both from file + param. time dep. anthropogenic 'simple plumes'
    - 19: tropospheric natural 'Kinne' aerosols from pre-industry (the 1850--file has to be linked for all simulations!) + param. time dep. anthropogenic 'simple plumes'
  -

* - (radiation_nml-lrad_aero_diag)=
    **lrad_aero_diag**
  - L
  - .FALSE.
  -
  - writes actual aerosol optical properties to output
  -

* - (radiation_nml-ecrad_data_path)=
    **ecrad_data_path**
  - C
  - "."
  -
  - Path to the folder containing ecRad optical properties files.
  - inwp_radiation=4 (ecRad)

* - (radiation_nml-cams_aero_filename)=
    **cams_aero_filename**
  - C
  - "CAMS_aero_
    R\<nroot0>B\<jlev>_
    DOM\<idom>.nc"
  -
  - Path to the file containing CAMS 3D data
    Climatology data can be prepared using the script scripts/preprocessing/
    make_camsclim_onICONgrid.sh
  - inwp_radiation=4 (ecRad)
    irad_aero=7 or 8 (CAMS 3D climatology or forecasted aerosols)

* - (radiation_nml-ecrad_isolver)=
    **ecrad_isolver**
  - I
  - 0
  -
  - Radiation solver
    - 0: McICA (Pincus et al. 2003)
    - 1: Tripleclouds (Shonk and Hogan 2008)
    - 2: McICA for OpenACC
    - 3: SPARTACUS (Hogan et al. 2016)
  - inwp_radiation=4 (ecRad)

* - (radiation_nml-ecrad_igas_model)=
    **ecrad_igas_model**
  - I
  - 0
  -
  - Gas model and spectral bands
    - 0: RRTMG (Iacono et al. 2008)
    - 1: ecckd (Hogan and Matricardi 2020)
  - inwp_radiation=4 (ecRad)

* - (radiation_nml-ecrad_llw_cloud_scat)=
    **ecrad_llw_cloud_scat**
  - L
  - .FALSE.
  -
  - Long-wave cloud scattering.
  - inwp_radiation=4 (ecRad)

* - (radiation_nml-ecrad_use_general_cloud_optics)=
    **ecrad_use_general_cloud_optics**
  - L
  - .FALSE.
  -
  - General cloud optics in ecrad.
    It allows for different optical properties of ice and liquid.
  - inwp_radiation=4 (ecRad)

* - (radiation_nml-ecrad_iliquid_scat)=
    **ecrad_iliquid_scat**
  - I
  - 0
  -
  - Optical properties for liquid cloud scattering.
    Different options depending on ecrad_use_general_cloud_optics (eugco)
    - 0: SOCRATES (eugco = .FALSE.) / Mie-droplet (eugco = .TRUE.)
    - 1: Slingo (1989) (eugco = .FALSE.)
  - inwp_radiation=4 (ecRad)

* - (radiation_nml-ecrad_iice_scat)=
    **ecrad_iice_scat**
  - I
  - 0
  -
  - Optical properties for ice cloud scattering.
    Different options depending on ecrad_use_general_cloud_optics (eugco)
    - 0: Fu et al. (1996) (eugco = .FALSE.) / Rough from Muskatel et al. (2021) (eugco = .TRUE.)
    - 1: Baran et al. (2016) (eugco=.FALSE.)
    - 2: Yi (2013) (eugco=.FALSE.)
    - 10: Smooth from Muskatel et al. (2021) (eugco=.TRUE.)
    - 11: Baum (eugco=.TRUE.)
  - inwp_radiation=4 (ecRad)

* - (radiation_nml-ecrad_isnow_scat)=
    **ecrad_isnow_scat**
  - I
  - -1
  -
  - Optical properties for snow scattering.
    - -1: No explicit snow in radiation.
    - 0: Rough from Muskatel et al. (2021)
    - 10: Smooth from Muskatel et al. (2021)
  - inwp_radiation=4 (ecRad)
    ecrad_use_general_cloud_optics = .TRUE.

* - (radiation_nml-ecrad_irain_scat)=
    **ecrad_irain_scat**
  - I
  - -1
  -
  - Optical properties for rain  scattering.
    - -1: No explicit rain in radiation.
    - 0: Mie-rain
  - inwp_radiation=4 (ecRad)
    ecrad_use_general_cloud_optics = .TRUE.

* - (radiation_nml-ecrad_igraupel_scat)=
    **ecrad_igraupel_scat**
  - I
  - -1
  -
  - Optical properties for graupel scattering.
    - -1: No explicit graupel in radiation.
    - 0: Rough from Muskatel et al. (2021)
    - 10: Smooth from Muskatel et al. (2021)
  - inwp_radiation=4 (ecRad)
    ecrad_use_general_cloud_optics = .TRUE.

* - (radiation_nml-decorr_pole)=
    **decorr_pole**
  - R
  - 2000
  - m
  - Decorrelation length scale at poles
  - inwp_radiation=4 (ecRad)

* - (radiation_nml-decorr_equator)=
    **decorr_equator**
  - R
  - 2000
  - m
  - Decorrelation length scale at equator
  - inwp_radiation=4 (ecRad)

* - (radiation_nml-ecrad_check_input)=
    **ecrad_check_input**
  - L
  - .FALSE.
  -
  - Debug mode for ecRad input:
      - Checks several input fields for physical consistency
      - Increases the verbosity of ecRad (setup and runtime)
  - inwp_radiation=4 (ecRad)

* - (radiation_nml-fsd_background)=
    **fsd_background**
  - R
  - 1
  -
  - Background value for fractional standard deviation (FSD)
  - inwp_radiation=4 (ecRad)

* - (radiation_nml-fsd_gridlen(max_dom))=
    **fsd_gridlen(max_dom)**
  - R
  - 80
  - km
  - Value for assumed horizontal grid spacing in FSD parameterization
  - inwp_radiation=4 (ecRad)

* - (radiation_nml-lcalculate_fsd)=
    **lcalculate_fsd**
  - L
  - .FALSE.
  -
  - Calculates regime-dependent FSD for radiation if .TRUE.
  - inwp_radiation=4 (ecRad)
```



Defined and used in: {{ '[src/namelists/mo_radiation_nml.f90]({}/src/namelists/mo_radiation_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_run_nml)=
## run_nml


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (run_nml-nsteps)=
    **nsteps**
  - I
  - -999
  -
  - Number of time steps of this run. Allowed range is {math}`\ge0`; setting a value of 0 allows writing initial output (including internal remapping) without calculating time steps.
  -

* - (run_nml-dtime)=
    **dtime**
  - R
  - 600.0
  - s
  - model time step
    For real case runs the maximum allowable time step can be estimated as
    {math}`1.8 \cdot ndyn_substeps \cdot \overline{\Delta x}\,\mathrm{s\, km^{-1}}`,
    where {math}`\overline{\Delta x}` is the average resolution in {math}`\mathrm{km}` and ndyn_substeps is the number of dynamics substeps.  ndyn_substeps should not be increased beyond the default value {math}`5`.
  -

* - (run_nml-modelTimeStep)=
    **modelTimeStep**
  - C
  - ''
  - ISO8601 formatted string
  - model time step (**should be preferred over the concurrent namelist parameter dtime**)
    For real case runs the maximum allowable time step can be estimated as
    {math}`1.8 \cdot` ndyn_substeps {math}`\cdot \overline{\Delta x}\,\mathrm{s\, km^{-1}}`,
    where {math}`\overline{\Delta x}` is the average resolution in {math}`\mathrm{km}` and ndyn_substeps is the number of dynamics substeps. ndyn_substeps should not be increased beyond the default value {math}`5`.
  -

* - (run_nml-ltestcase)=
    **ltestcase**
  - L
  - .TRUE.
  -
  - Idealized testcase runs
  -

* - (run_nml-ldynamics)=
    **ldynamics**
  - L
  - .TRUE.
  -
  - Compute adiabatic dynamic tendencies
  -

* - (run_nml-iforcing)=
    **iforcing**
  - I
  - 0
  -
  - Forcing of dynamics and transport by parameterized processes. Use positive indices for the atmosphere and negative indices for the ocean.
    - 0: no forcing
    - 1: Held-Suarez forcing
    - 2: AES forcing
    - 3: NWP forcing
    - 4: local diabatic forcing without physics
    - 5: local diabatic forcing with physics
    - -1: MPIOM forcing (to be done)
  -

* - (run_nml-ltransport)=
    **ltransport**
  - L
  - .FALSE.
  -
  - Compute large-scale tracer transport
  -

* - (run_nml-ntracer)=
    **ntracer**
  - I
  - 0
  -
  - Number of advected tracers handled by the large-scale transport scheme
  -

* - (run_nml-lvert_nest)=
    **lvert_nest**
  - L
  - .FALSE.
  -
  - If set to .TRUE. vertical nesting is switched on (i.e.
    variable number of vertical levels)
  -

* - (run_nml-num_lev)=
    **num_lev**
  - I(max_ dom)
  - 31
  -
  - Number of full levels (atm.) for each domain
  - lvert_nest=.TRUE.

* - (run_nml-nshift)=
    **nshift**
  - I(max_ dom)
  - 0
  -
  - vertical half level of parent domain which coincides with upper boundary of the current domain required for vertical refinement, which is **not yet implemented**
  - lvert_nest=.TRUE.

* - (run_nml-ltimer)=
    **ltimer**
  - L
  - .TRUE.
  -
  - TRUE: Timer for monitoring the runtime of specific routines is on (FALSE = off)
  -

* - (run_nml-timers_level)=
    **timers_level**
  - I
  - 1
  -
  -
  -

* - (run_nml-activate_sync_timers)=
    **activate_sync_timers**
  - L
  - F
  -
  - TRUE: Timer for monitoring runtime of communication routines (FALSE = off)
  -

* - (run_nml-msg_level)=
    **msg_level**
  - I
  - 10
  -
  - controls how much printout is written during runtime.
    For values less than 5, only the time step is written.
  -

* - (run_nml-msg_timestamp)=
    **msg_timestamp**
  - L
  - .FALSE.
  -
  - If .TRUE., precede output messages by time stamp.
  -

* - (run_nml-debug_check_level)=
    **debug_check_level**
  - I
  - 0
  -
  - Setting a value larger than 0 activates debug checks.
  -

* - (run_nml-output)=
    **output**
  - C(:)
  - "nml", "totint"
  -
  - Main switch for enabling/disabling components of the model output. One or more choices can be set (as an array of string constants). Possible choices are:
      - "none": switch off all output;
      - "nml": new output mode (cf.
    [output_nml](ref_buildrun_nml_output_nml));
      - "totint": computation of total integrals.
      - "maxwinds": write max. winds to separate ASCII file "maxwinds.log".
    If the `output` namelist parameter is not set explicitly, the default setting "nml","totint" is assumed.
  -

* - (run_nml-restart_filename)=
    **restart_filename**
  - C
  -
  -
  - File name for restart/checkpoint files (containing keyword substitution patterns `<gridfile>`, `<idom>`, `<rsttime>`, `<mtype>`). default: "`<gridfile>_restart_<mtype>_<rsttime>.nc`".
  -

* - (run_nml-profiling_output)=
    **profiling_output**
  - I
  - 1
  -
  - controls how profiling printout is written: TIMER_MODE_AGGREGATED=1, TIMER_MODE_DETAILED=2, TIMER_MODE_WRITE_FILES=3.
  -

* - (run_nml-lart)=
    **lart**
  - L
  - .FALSE.
  -
  - Main switch which enables the treatment of atmospheric aerosol and trace gases (The ART package of KIT is needed for this purpose)
  -

* - (run_nml-ldass_lhn)=
    **ldass_lhn**
  - L
  - .FALSE.
  -
  - Main switch which enables the assimilation of radar derived precipitation rate via Latent Heat Nudging
  -

* - (run_nml-check_uuid_gracefully)=
    **check_uuid_gracefully**
  - L
  - .FALSE.
  -
  - If this flag is set to `.TRUE.` we give only warnings for non-matching UUIDs.
  -

* - (run_nml-luse_radarfwo)=
    **luse_radarfwo**
  - L(max_ dom)
  - .FALSE.
  -
  - For each domain, switch to activate the efficient volume scan radar forward operator EMVORADO. The EMVORADO code is provided as a submodule named `emvorado`, which is part of the ICON distribution. ICON itself contains only some ICON specific interface modules.
    `./configure` (respectively the call to a configure wrapper script) needs the option `--enable-emvorado`.
    EMVORADO needs its own namelist(s) for each radar-active model domain in a separate namelist input file `RADARSIM_PARAMS`. More details can be found in the EMVORADO User's Guide available from the COSMO web page (<https://www.cosmo-model.org> {math}`\rightarrow` Documentation {math}`\rightarrow` EMVORADO) or from the submodule `externals/emvorado/DOC/TEX/emvorado_userguide.pdf`.
  - [iforcing](run_nml-iforcing)=3,
    ICON configure'd with `--enable-emvorado`

* - (run_nml-radarnmlfile)=
    **radarnmlfile**
  - C
  -
  -
  - The name of the file containing the EMVORADO namelist. If this is empty or not set, the Default from {{ '[externals/emvorado/src_emvorado/radar_data_namelist.f90]({}/externals/emvorado/src_emvorado/radar_data_namelist.f90)'.format(base_url) }} is used. Only used if luse_radarfwo is .TRUE. .
  -

* - (run_nml-lmsgwam)=
    **lmsgwam**
  - L
  - .FALSE.
  -
  - Main switch which enables the online 3D gravity waves parametrisation: Multi-Scale Gravity Wave Model MS-GWaM. The MS-GWaM submodule of GU is needed. See the detailed <https://gitlab.dkrz.de/atmodynamics-goethe-universitaet-frankfurt/msgwam-release/-/tree/3e6af3fa885aada7f960634d2334d81bb9fea87d/doc>{MSGWAM namelist documentation} at `./external/msgwam/doc/`.
  -
```



Defined and used in: {{ '[src/namelists/mo_run_nml.f90]({}/src/namelists/mo_run_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_scm_nml)=
## scm_nml

Relevant if l_scm_mode

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (scm_nml-i_scm_netcdf)=
    **i_scm_netcdf**
  - I
  - 1
  -
  - reading SCM input data from
    - 0: ASCII file
    - 1: normal ICON netcdf file
    - 2: DEPHY unified netcdf file
  -

* - (scm_nml-lscm_icon_ini)=
    **lscm_icon_ini**
  - L
  - .FALSE.
  -
  - read initial conditions produced by ICON on the native grid
  -

* - (scm_nml-lscm_random_noise)=
    **lscm_random_noise**
  - L
  - .FALSE.
  -
  - initialize with random noise - for LEM runs by ICON on the native grid
  -

* - (scm_nml-lscm_read_tke)=
    **lscm_read_tke**
  - L
  - .FALSE.
  -
  - read init. tke from netcdf
  -

* - (scm_nml-lscm_read_z0)=
    **lscm_read_z0**
  - L
  - .FALSE.
  -
  - read z0 from netcdf
  -

* - (scm_nml-scm_sfc_mom)=
    **scm_sfc_mom**
  - I
  - 0
  -
  - prescribed surface boundary condition for momentum using
    - 0: TERRA
    - 2: friction velocity
    - 4: drag coefficient
    - 5: Louis surface layer scheme
  -

* - (scm_nml-scm_sfc_qv)=
    **scm_sfc_qv**
  - I
  - 0
  -
  - prescribed surface boundary condition for moisture using
    - 0: TERRA
    - 1: surface moisture (qv_s)
    - 2: latent heat flux
    - 3: saturation
    - 4: draf coefficient
    - 5: Louis surface layer scheme
  -

* - (scm_nml-scm_sfc_temp)=
    **scm_sfc_temp**
  - I
  - 0
  -
  - prescribed surface boundary condition for temperature using
    - 0: TERRA
    - 1: surface temperature (t_g)
    - 2: sensible heat flux (shfl_s)
    - 4: drag coefficient
    - 5: Louis surface layer scheme
  -
```



Defined and used in: {{ '[src/namelists/mo_scm_nml.f90]({}/src/namelists/mo_scm_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_sleve_nml)=
## sleve_nml

Relevant if ivctype=2

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (sleve_nml-min_lay_thckn)=
    **min_lay_thckn**
  - R
  - 50
  - m
  - Layer thickness of lowermost layer; specifying zero or a negative value leads to constant layer thicknesses determined by top_height and nlev
  -

* - (sleve_nml-max_lay_thckn)=
    **max_lay_thckn**
  - R
  - 25000
  - m
  - Maximum layer thickness below the height given by htop_thcknlimit (NWP recommendation: 400 m)
    _Use with caution! Too ambitious settings may result in numerically unstable layer configurations._
  -

* - (sleve_nml-htop_thcknlimit)=
    **htop_thcknlimit**
  - R
  - 15000
  - m
  - Height below which the layer thickness does not exceed max_lay_thckn
  -

* - (sleve_nml-nshift_above_thcklay)=
    **nshift_above_thcklay**
  - I
  - 0
  -
  - Level shift above constant-thickness layer for further calculation of layer distribution. For strongly stretched grids with a deep constant-thickness layer, this parameter may be set to 1 in order to reduce the thickness jump right above the constant-thickness layer.
  -

* - (sleve_nml-itype_laydistr)=
    **itype_laydistr**
  - I
  - 1
  -
  - Type of analytical function used to specify the distribution of the vertical coordinate surfaces
    - 1: transformed cosine,
    - 2: third-order polynomial; in this case, stretch_fac should be less than 1, particularly for large numbers of model levels; the algorithm always works for stretch_fac=0.5
    - 3: second-order polynomial (see M. Baldauf COSMO-TR p. 33)
  -

* - (sleve_nml-top_height)=
    **top_height**
  - R
  - 23500.0
  - m
  - Height of model top
  -

* - (sleve_nml-stretch_fac)=
    **stretch_fac**
  - R
  - 1.0
  -
  - Stretching factor to vary distribution of model levels; values < 1 increase the layer thickness near the model top
  -

* - (sleve_nml-decay_scale_1)=
    **decay_scale_1**
  - R
  - 4000
  - m
  - Decay scale of large-scale topography component
  -

* - (sleve_nml-decay_scale_2)=
    **decay_scale_2**
  - R
  - 2500
  - m
  - Decay scale of small-scale topography component
  -

* - (sleve_nml-decay_exp)=
    **decay_exp**
  - R
  - 1.2
  -
  - Exponent of decay function
  -

* - (sleve_nml-flat_height)=
    **flat_height**
  - R
  - 16000
  - m
  - Height above which the coordinate surfaces are flat
  -

* - (sleve_nml-lread_smt)=
    **lread_smt**
  - L
  - .FALSE.
  -
  - read smoothed topography from file (TRUE) or compute internally (FALSE)
  -
```


Defined and used in: {{ '[src/namelists/mo_sleve_nml.f90]({}/src/namelists/mo_sleve_nml.f90)'.format(base_url) }}

(ref_buildrun_nml_sppt_nml)=
## sppt_nml

The Stochastic Perturbation of Physical Tendencies (SPPT) method is controlled by the following set of Namelist parameters. Note that SPPT is only available for the NWP physics package ([iforcing](run_nml-iforcing)=3).
In addition, SPPT is not supported on a global domain (hard exit) and is untested in limited area mode where the domain extends across the poles. Running the latter is currently not recommended.


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (sppt_nml-lsppt)=
    **lsppt**
  - L
  - .FALSE.
  -
  - TRUE: forecast with SPPT
  -

* - (sppt_nml-hinc_rn)=
    **hinc_rn**
  - R
  - 21600
  - second
  - time increment for drawing a new field of random numbers
  -

* - (sppt_nml-dlat_rn)=
    **dlat_rn**
  - R
  - 0.1
  - deg
  - random number coarse grid point distance in meridional direction
  -

* - (sppt_nml-dlon_rn)=
    **dlon_rn**
  - R
  - 0.1
  - deg
  - random number coarse grid point distance in zonal direction
  -

* - (sppt_nml-range_rn)=
    **range_rn**
  - R
  - 0.8
  -
  - max magnitude of random numbers
  -

* - (sppt_nml-stdv_rn)=
    **stdv_rn**
  - R
  - 1.0
  -
  - standard deviation of the gaussian distribution of random numbers
  -
```



Defined and used in: {{ '[src/namelists/mo_sppt_nml.f90]({}/src/namelists/mo_sppt_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_synsat_nml)=
## synsat_nml

(This feature is currently active for configuration `dwd+cray` only)

This namelist enables the RTTOV library incorporated into ICON for simulating satellite radiance and brightness temperatures.

RTTOV is a radiative transfer model for nadir-viewing passive visible, infrared and microwave satellite radiometers, spectrometers and interferometers.


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (synsat_nml-lsynsat)=
    **lsynsat**
  - L (max_dom)
  - .FALSE.
  -
  - Main switch: Enables/disables computation of synthetic satellite imagery for each model domain.
  -

* - (synsat_nml-nlev_rttov)=
    **nlev_rttov**
  - I
  - 51
  -
  - Number of RTTOV levels.
  -
```



Enabling the synsat module makes the following 32 two-dimensional output fields available:


```{list-table}
:header-rows: 1

* -
  -
  -
  -

* - `SYNMSG_RAD_CL_IR3.9`
  - `SYNMSG_BT_CL_IR3.9`
  - `SYNMSG_RAD_CL_WV6.2`
  - `SYNMSG_BT_CL_WV6.2`

* - `SYNMSG_RAD_CL_WV7.3`
  - `SYNMSG_BT_CL_WV7.3`
  - `SYNMSG_RAD_CL_IR8.7`
  - `SYNMSG_BT_CL_IR8.7`

* - `SYNMSG_RAD_CL_IR9.7`
  - `SYNMSG_BT_CL_IR9.7`
  - `SYNMSG_RAD_CL_IR10.8`
  - `SYNMSG_BT_CL_IR10.8`

* - `SYNMSG_RAD_CL_IR12.1`
  - `SYNMSG_BT_CL_IR12.1`
  - `SYNMSG_RAD_CL_IR13.4`
  - `SYNMSG_BT_CL_IR13.4`

* - `SYNMSG_RAD_CS_IR3.9`
  - `SYNMSG_BT_CS_IR3.9`
  - `SYNMSG_RAD_CS_WV6.2`
  - `SYNMSG_BT_CS_WV6.2`

* - `SYNMSG_RAD_CS_WV7.3`
  - `SYNMSG_BT_CS_WV7.3`
  - `SYNMSG_RAD_CS_IR8.7`
  - `SYNMSG_BT_CS_IR8.7`

* - `SYNMSG_RAD_CS_IR9.7`
  - `SYNMSG_BT_CS_IR9.7`
  - `SYNMSG_RAD_CS_IR10.8`
  - `SYNMSG_BT_CS_IR10.8`

* - `SYNMSG_RAD_CS_IR12.1`
  - `SYNMSG_BT_CS_IR12.1`
  - `SYNMSG_RAD_CS_IR13.4`
  - `SYNMSG_BT_CS_IR13.4`
```



Here,
`RAD` denotes radiance,
`BT` brightness temperature,
`CL` cloudy, and
`CS` clear sky,
supplemented by the channel name.

Defined and used in: {{ '[src/namelists/mo_synsat_nml.f90]({}/src/namelists/mo_synsat_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_synradar_nml)=
## synradar_nml

The list of diagnostic output variables in ICON incorporates some fields related to synthetic radar reflectivity on the model grid:

 - `'dbz', 'dbz_850', 'dbz_cmax', 'dbz_ctmax'`
 - `'echotop', 'echotopinm'`

By default, these are based on a simple analytic so-called Rayleigh-approximation for single-particle backscattering.

If ICON is configured with the flag `--enable-emvorado` and compiled with the pre-processor flag `-DHAVE_RADARFWO`, some alternative, more accurate Mie- or T-matrix methods from the radar forward operator EMVORADO can be used by namelist choice (see below), particularly for improving the simulation of the so-called "bright band", the enhanced reflectivity in the melting layer.


EMVORADO is the Efficient Modular VOlume RADar Operator for simulating radar volume scans for cloud- and weather radar wavelengths, see

 - EMVORADO User's Guide in ICON's EMVORADO submodule `externals/emvorado/DOC/TEX/emvorado_userguide.pdf` or from the COSMO webpage <https://www.cosmo-model.org/content/model/documentation/core/emvorado_userguide.pdf>
 - [COSMO Technical Report No. 28](https://www.cosmo-model.org/content/model/documentation/techReports/cosmo/docs/techReport28.pdf)

for detailed information.


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (synradar_nml-synradar_meta_)=
    **synradar_meta%...**
  - TYPE(dbzcalc_params)
    I
  - 4
  -
  - This type contains: synradar_meta%itype_refl and many other parameters which are only relevant if itype_refl is not the default (4). Instance of the derived type   **dbzcalc_params** from EMVORADO to specify details of the radar reflectivity calculation for related outputs (`'dbz'`, `'dbz_850'`, `'dbz_cmax'`, `'dbz_ctmax'`, `'echotop'`, `'echotopinm'`). The type is documented in detail in the EMVORADO User's Guide (link see above).
    The most important component is **itype_refl:**  
      - **1:** Mie-scattering from EMVORADO assuming spherical particles and including a detailed melting scheme for the radar "bright band".  
      - **3:** Rayleigh-Oguchi approximation from EMVORADO including a simple melting scheme, but not producing pronounced "bright bands".  
      - **4:** Traditional Rayleigh approximation from ICON, also without pronounced "bright bands". **This is the default.**  
      - **5:** T-matrix scattering from EMVORADO assuming oblate spheroids, otherwise similar to Mie-option 1.  
      - **6:** T-matrix scattering from EMVORADO assuming spherical particles, only for sanity-checks against Mie-option 1.  
    For options 1, 5, 6 there are many more relevant type components.
  - [iforcing](run_nml-iforcing)=3,
    ICON configure'd with `--enable-emvorado`

* - (synradar_nml-ydir_mielookup_write)=
    **ydir_mielookup_write**
  - C
  - ' '
  -
  - For reflectivity calculations: directory for storing new automatically created reflectivity lookup tables in case of EMVORADO-methods that employ reflectivity lookup tables to boost efficiency (synradar_meta%itype_refl=1, 5, 6 together with synradar_meta%llookup_mie=.TRUE.)
  - [iforcing](run_nml-iforcing)=3,
    ICON configure'd with `--enable-emvorado`,
    synradar_meta%itype_refl=1, 5, 6
    synradar_meta%llookup_mie=.TRUE.

* - (synradar_nml-ydir_mielookup_read)=
    **ydir_mielookup_read**
  - C
  - ' '
  -
  - For reflectivity calculations: directory for reading the reflectivity lookup tables in case of EMVORADO-methods that employ reflectivity lookup tables to boost efficiency (synradar_meta%itype_refl=1, 5, 6 together with synradar_meta%llookup_mie=.TRUE.)
  - [iforcing](run_nml-iforcing)=3,
    ICON configure'd with `--enable-emvorado`,
    synradar_meta%itype_refl=1, 5, 6
    synradar_meta%llookup_mie=.TRUE.

* - (synradar_nml-rain2mom_mu_incloud)=
    **rain2mom_mu_incloud**
  - R
  - -999.99
  -
  - For Mie- or T-matrix calculations together with the two-moment mircophysics: value to overwrite the shape parameter {math}`\mu` of rain drop size distribution in regions where {math}`q_c > 0` (the "incloud" value), but only in the EMVORADO, not in the microphysics itself! This is to provide a possibility for forward operator tuning with regards to polarimetric radar parameters.
    If rain2mom_mu_incloud{math}`\le`-1.0, then the original value of the two-moment microphysics is used.
  - [iforcing](run_nml-iforcing)=3,
    ICON configure'd with `--enable-emvorado`,
    synradar_meta%itype_refl=1, 5, 6
    inwp_gscp=4, 5, 6, 7

* - (synradar_nml-itype_Dlim_sgh)=
    **itype_Dlim_sgh**
  - I
  - 0
  -
  - Options (0,1, or 2) to limit excessive radar signals coming from the large tail of (some) hydrometeor PSDs. When integrating radar moments over the PSD, the mass of particles larger than a certain diameter {math}`D_c` (separate namelist parameters for hydrometeors, see below) is redistributed to
      - **1:** particles of exact size {math}`D_c`
      - **2:** to all particles with {math}`D < D_c` by scaling up the PSD there
    before integrating over the PSD. If set to 0 (default), no such redistribution is applied. More details in the EMVORADO User's Guide, see above.
  - [iforcing](run_nml-iforcing)=3,
    ICON configure'd with `--enable-emvorado`,
    synradar_meta%itype_refl=1, 5, 6

* - (synradar_nml-Dlim_rain)=
    **Dlim_rain**
  - R
  - 999.0
  -
  - For itype_Dlim_sgh>0: limiting diameter [m] for rain. Does not have any effect if set to a value {math}` > D_{max,r}` (upper integration limit in EMVORADO). The default is such a large value.
  - [iforcing](run_nml-iforcing)=3,
    ICON configure'd with `--enable-emvorado`,
    synradar_meta%itype_refl=1, 5, 6

* - (synradar_nml-Dlim_drysnow)=
    **Dlim_drysnow**
  - R
  - 999.0
  -
  - For itype_Dlim_sgh>0: limiting diameter [m] for dry snow. Does not have any effect if set to a value {math}` > D_{max,s}` (upper integration limit in EMVORADO). The default is such a large value.
  - [iforcing](run_nml-iforcing)=3,
    ICON configure'd with `--enable-emvorado`,
    synradar_meta%itype_refl=1, 5, 6

* - (synradar_nml-Dlim_meltsnow)=
    **Dlim_meltsnow**
  - R
  - 999.0
  -
  - For itype_Dlim_sgh>0: limiting diameter [m] for melting snow. Does not have any effect if set to a value {math}`> D_{max,s}` (upper integration limit in EMVORADO). The default is such a large value.
  - [iforcing](run_nml-iforcing)=3,
    ICON configure'd with `--enable-emvorado`,
    synradar_meta%itype_refl=1, 5, 6

* - (synradar_nml-Dlim_meltgraupel)=
    **Dlim_meltgraupel**
  - R
  - 999.0
  -
  - For itype_Dlim_sgh>0: limiting diameter [m] for melting graupel. Does not have any effect if set to a value {math}`> D_{max,g}` (upper integration limit in EMVORADO). The default is such a large value.
  - [iforcing](run_nml-iforcing)=3,
    ICON configure'd with `--enable-emvorado`,
    synradar_meta%itype_refl=1, 5, 6

* - (synradar_nml-Dlim_melthail)=
    **Dlim_melthail**
  - R
  - 999.0
  -
  - For itype_Dlim_sgh>0: limiting diameter [m] for melting hail. Does not have any effect if set to a value {math}`> D_{max,h}` (upper integration limit in EMVORADO). The default is such a large value.
  - [iforcing](run_nml-iforcing)=3,
    ICON configure'd with `--enable-emvorado`,
    synradar_meta%itype_refl=1, 5, 6
```



Defined and used in: {{ '[src/namelists/mo_synradar_nml.f90]({}/src/namelists/mo_synradar_nml.f90)'.format(base_url) }}




(ref_buildrun_nml_testbed_nml)=
## testbed_nml

Configuration for ICON's testbed mode. The testbed is independent of the other model components for atmosphere, ocean,
etc. It is primarily used for unit, function, and performance testing of infrastructure components.


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (testbed_nml-testbed_model)=
    **testbed_model**
  - I
  - 0
  -
  - Type of testbed model (non-zero to activate testbed):
    - 0: Null model
    - 3: Test jitter model
    - 4: Test mpi halo communication
    - 6: Test netcdf reading
    - 8: Test mpi gather communication
    - 9: Test mpi exchange communication
    - 10: Benchmark mpi multi-field exchange communication
    - 11: Test mpi sync_patch communication
  -

* - (testbed_nml-testbed_iterations)=
    **testbed_iterations**
  - I
  - 10
  -
  - Number of iterations to run testbed (if applicable)
  -

* - (testbed_nml-calculate_iterations)=
    **calculate_iterations**
  - I
  - 10
  -
  - Number of iterations to run calculate step (if applicable)
  -

* - (testbed_nml-no_of_blocks)=
    **no_of_blocks**
  - I
  - 16
  -
  - Number of blocks in field(if testbed_model==3)
  -

* - (testbed_nml-no_of_layers)=
    **no_of_layers**
  - I
  - 80
  -
  - Number of vertical layers in field(if testbed_model==3)
  -

* - (testbed_nml-testfile_3D_file)=
    **testfile_3D_file**
  - C
  - ' '
  -
  - Filename of 3D input file in first entry, and variable name in second entry (if testbed_model==6)
  -

* - (testbed_nml-testfile_2D_file)=
    **testfile_2D_file**
  - C
  - ' '
  -
  - Array: first entry is filename of 2D input file, and second entry is a variable name (if testbed_model==6)
  -
```



Defined and used in: {{ '[src/namelists/mo_icon_testbed_nml.f90]({}/src/namelists/mo_icon_testbed_nml.f90)'.format(base_url) }}. The testbed 'driver'
is in: {{ '[src/testbed/mo_icon_testbed.f90]({}/src/testbed/mo_icon_testbed.f90)'.format(base_url) }}.




(ref_buildrun_nml_time_nml)=
## time_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (time_nml-calendar)=
    **calendar**
  - I
  - 1
  -
  - Calendar type:
    - 0: Julian/Gregorian
    - 1: proleptic Gregorian
    - 2: 30day/month, 360day/year
  -

* - (time_nml-dt_restart)=
    **dt_restart**
  - R
  - 0.0
  - s
  - Length of restart cycle in seconds. This namelist parameter specifies how long the model runs until it saves its state to a file and stops. Later, the model run can be resumed, s.t.
    a simulation over a long period of time can be split into a chain of restarted model runs.  Note that the frequency of writing restart files is controlled by dt_checkpoint. Only if the value of `dt_checkpoint` resulting from model default or user's specification is longer than `dt_restart`, it will be reset (by the model) to `dt_restart` so that at least one restart file is generated during the restart cycle. If `dt_restart` is larger than but not a multiple of `dt_checkpoint`, restart file will _not_ be generated at the end of the restart cycle.
  -

* - (time_nml-ini_datetime_string)=
    **ini_datetime_string**
  - C
  - '2008- 09-01T 00:00:00Z'
  -
  - Initial date and time of the simulation
  -

* - (time_nml-end_datetime_string)=
    **end_datetime_string**
  - C
  - '2008- 09-01T 01:40:00Z'
  -
  - End date and time of the simulation
  -

* - (time_nml-is_relative_time)=
    **is_relative_time**
  - L
  - .FALSE.
  -
  - .TRUE., if time loop shall start with step 0 regardless whether we are in a standard run or in a restarted run (which means re-initialized run).
  -
```




### Length of the run


If nsteps is positive, then nsteps*dtime
is used to compute the end date and time of the run.
Else the initial date and time, the end date and time, dt_restart,
as well as the time step are used to compute nsteps.



(ref_buildrun_nml_transport_nml)=
## transport_nml

Used if ltransport=.TRUE.

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (transport_nml-lvadv_tracer)=
    **lvadv_tracer**
  - L
  - .TRUE.
  -
  - Main switch for vertical tracer transport.
    TRUE/FALSE : compute/do not compute vertical tracer advection.
    If vertical advection is switched off, the tracer mass fraction {math}`q` is kept constant.
  -

* - (transport_nml-ihadv_tracer)=
    **ihadv_tracer**
  - I(ntracer)
  - 2
  -
  - Tracer specific method to compute horizontal advection:
  -

* -
  -
  -
  -
  - 0: no horiz.
    transport. The tracer mass fraction {math}`q` is kept constant.
  -

* -
  -
  -
  -
  - 1: upwind (1st order)
  -

* -
  -
  -
  -
  - 2: Miura (2nd order, linear reconstr.)
  -

* -
  -
  -
  -
  - 3: Miura3 (quadr. or cubic reconstr.)
  - lsq_high_ord {math}`\in` [2,3]

* -
  -
  -
  -
  - 4: FFSL (quadr. or cubic reconstr.)
  - lsq_high_ord {math}`\in` [2,3]

* -
  -
  -
  -
  - 5: hybrid Miura3/FFSL (quadr. or cubic reconstr.)
  - lsq_high_ord {math}`\in` [2,3]

* -
  -
  -
  -
  - 20: miura (2nd order, lin. reconstr.) with subcycling
  -

* -
  -
  -
  -
  - 22: combination of miura and miura with subcycling
  -

* -
  -
  -
  -
  - 32: combination of miura3 and miura with subcycling
  -

* -
  -
  -
  -
  - 42: combination of FFSL and miura with subcycling
  -

* -
  -
  -
  -
  - 52: combination of hybrid FFSL/Miura3 with subcycling
    Subcycling means that the integration from time step n to n+1 is splitted into substeps to meet the stability requirements. For NWP runs, substepping is generally applied above {math}`z=22\,\mathrm{km}` (see hbot_qvsubstep).
  -

* - (transport_nml-ivadv_tracer)=
    **ivadv_tracer**
  - I(ntracer)
  - 3
  -
  - Tracer specific method to compute vertical advection:
  - lvadv_tracer=TRUE

* -
  -
  -
  -
  - 0: no vert.
    transport. The tracer mass fraction {math}`q` is kept constant.
  -

* -
  -
  -
  -
  - 1: upwind (1st order)
  -

* -
  -
  -
  -
  - 2: Parabolic Spline Method (PSM): allows for {math}`\mathrm{CFL}>1`
  -

* -
  -
  -
  -
  - 3: Piecewise parabolic method (PPM): allows for {math}`\mathrm{CFL}>1`
  -

* - (transport_nml-itype_hlimit)=
    **itype_hlimit**
  - I(ntracer)
  - 4
  -
  - Type of limiter for horizontal transport:
  -

* -
  -
  -
  -
  - 0: no limiter
  -

* -
  -
  -
  -
  - 3: monotonic flux limiter (FCT)
  -

* -
  -
  -
  -
  - 4: positive definite flux limiter
  -

* - (transport_nml-itype_vlimit)=
    **itype_vlimit**
  - I(ntracer)
  - 1
  -
  - Type of limiter for vertical transport:
  -

* -
  -
  -
  -
  - 0: no limiter
  -

* -
  -
  -
  -
  - 1: semi-monotonic reconstruction filter
  -

* -
  -
  -
  -
  - 2: monotonic reconstruction filter
  -

* -
  -
  -
  -
  - 3: positive definite flux limiter
  -

* - (transport_nml-ivlimit_selective)=
    **ivlimit_selective**
  - I(ntracer)
  - 0
  -
  - Reduce detrimental effect of vertical limiter by applying a method for identifying and avoiding spurious limiting of smooth extrema.
  -

* -
  -
  -
  -
  - 1: on
  - itype_vlimit=1, 2

* -
  -
  -
  -
  - 0: off
  -

* - (transport_nml-nadv_substeps)=
    **nadv_substeps**
  - I(max_dom)
  - 3
  -
  - Tracer substepping:
    Number of time integration substeps per fast physics/advective time step dtime.
    If only one value is specified, it is copied to all child domains, implying that the same value is used in all domains. If the number of values given in the namelist is larger than 1 but less than the number of model domains, then the settings from the highest domain ID are used for the remaining model domains.
  - only active for the schemes ihadv_tracer=20, 22, 32, 42, 52.
    Starts at minimum height height hbot_qv_substep for the schemes 22, 32, 42, 52, whereas it is applied throughout the entire domain for scheme 20.

* - (transport_nml-beta_fct)=
    **beta_fct**
  - R
  - 1.005
  -
  - global boost factor for range of permissible values {math}`\left[q_{max},q_{min}\right]` in the monotonic flux limiter. A value larger than 1 allows for (small) over and undershoots, while a value of {math}`1` gives strict monotonicity (at the price of increased diffusivity).
  - itype_hlimit = 3

* - (transport_nml-iadv_tke)=
    **iadv_tke**
  - I
  - 0
  -
  - Type of TKE advection
  - inwp_turb=1

* -
  -
  -
  -
  - 0: no TKE advection
  -

* -
  -
  -
  -
  - 1: vertical advection only
  -

* -
  -
  -
  -
  - 2: vertical and horizontal advection
  -

* - (transport_nml-tracer_names)=
    **tracer_names**
  - C(:)
  - 'Int2Str(i)'
  -
  - Tracer-specific name suffixes. When running idealized cases or the hydrostatic ICON, this variable is used to specify tracer names. If nothing is specified, the tracer name is given as `PREFIX+Int2String(i)`, where `i` is the tracer index. Note that this namelist variable has no effect for nonhydrostatic real-case runs, if the NWP- or AES physics packages are switched on.
  - [iforcing](run_nml-iforcing){math}`\ne` inwp, iaes'

* - (transport_nml-npassive_tracer)=
    **npassive_tracer**
  - I
  - 0
  -
  - number of additional passive tracers which have no sources and are transparent to any physical process (no effect).
    Passive tracers are named Qpassive_ID, where ID is a number between `ntracer` and `ntracer`{math}`+``npassive_tracer`.
    **NOTE:** By default, limiters are switched off for passive tracers and the scheme {math}`52` is selected for horizontal advection.
  -

* - (transport_nml-init_formula)=
    **init_formula**
  - C
  - ' '
  -
  - Comma-separated list of initialization formulas for additional passive tracers.
  - npassive_tracer > 0

* - (transport_nml-igrad_c_miura)=
    **igrad_c_miura**
  - I
  - 1
  -
  - Method for gradient reconstruction at cell center for 2nd order miura scheme
  -

* -
  -
  -
  -
  - 1: Least-squares (linear, non-consv)
  - ihadv_tracer=2, 20

* -
  -
  -
  -
  - 2: Green-Gauss
  -

* -
  -
  -
  -
  - 3: based on shape function derivatives for a three-node triangular element (Fish. J and T. Belytschko, 2007)
  -

* - (transport_nml-ivcfl_max)=
    **ivcfl_max**
  - I
  - 5
  -
  - determines stability range of vertical PPM/PSM-scheme in terms of the maximum allowable CFL-number
  - ivadv_tracer=3,4

* - (transport_nml-llsq_svd)=
    **llsq_svd**
  - L
  - .TRUE.
  -
  - use QR decomposition (FALSE) or SV decomposition (TRUE) for least squares design matrix A
  -

* - (transport_nml-lclip_tracer)=
    **lclip_tracer**
  - L
  - .FALSE.
  -
  - Clipping of negative values
  -
```



Defined and used in: {{ '[src/namelists/mo_advection_nml.f90]({}/src/namelists/mo_advection_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_turb_vdiff_nml)=
## turb_vdiff_nml

The parameterization of vertical diffusion (VDIFF) module is configured by a a set of parameters, each of which is a 1-dimensional array extending over all  domains.
The parameters provide control over some of the parametrized effects (only active when inwp_turb = 6):

### General

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (turb_vdiff_nml-lsfc_mom_flux)=
    **lsfc_mom_flux**
  - L
  - .TRUE.
  -
  - switch on/off surface momentum flux
  -

* - (turb_vdiff_nml-lsfc_heat_flux)=
    **lsfc_heat_flux**
  - L
  - .TRUE.
  -
  - switch on/off surface heat flux
  -

* - (turb_vdiff_nml-turb)=
    **turb**
  - S
  - 'tte'
  -
  - 'tte': TTE scheme
    '3dsmag': 3D Smagorinsky scheme
  -

* - (turb_vdiff_nml-z0m_min)=
    **z0m_min**
  - R
  - {math}`1.5\times 10^{-5}`
  - m
  - Minimum roughness length for momentum
  -

* - (turb_vdiff_nml-z0m_ice)=
    **z0m_ice**
  - R
  - 0.001
  - m
  - Roughness length for momentum over ice
  -

* - (turb_vdiff_nml-z0m_oce)=
    **z0m_oce**
  - R
  - 0.001
  - m
  - Roughness length for momentum over ocean
  -

* - (turb_vdiff_nml-fsl)=
    **fsl**
  - R
  - 0.4
  -
  - fraction of first-level height at which surface fluxes are nominally evaluated, tuning param for sfc stress
  -
```



### TTE Scheme

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (turb_vdiff_nml-pr0)=
    **pr0**
  - R
  - 1.0
  -
  - neutral limit Prandtl number, can be varied from about 0.6 to 1.0, fixes f_theta0
  -

* - (turb_vdiff_nml-f_tau0)=
    **f_tau0**
  - R
  - 0.17
  -
  - neutral non-dimensional stress factor (0.1 - 0.22)
  -

* - (turb_vdiff_nml-f_tau_limit_fraction)=
    **f_tau_limit_fraction**
  - R
  - 0.25
  -
  - Fraction of [f_tau0](turb_vdiff_nml-f_tau0) for large Ri numbers (0 - 0.6)
  -

* - (turb_vdiff_nml-f_theta_limit_fraction)=
    **f_theta_limit_fraction**
  - R
  - 0.0
  -
  - Fraction of f_theta0 for large Ri numbers (0 - 0.3)
  -

* - (turb_vdiff_nml-f_tau_decay)=
    **f_tau_decay**
  - R
  - 4.0
  -
  - Decay constant of [f_tau0](turb_vdiff_nml-f_tau0) for large Ri numbers (0.5 - 5)
  -

* - (turb_vdiff_nml-f_theta_decay)=
    **f_theta_decay**
  - R
  - 4.0
  -
  - Decay constant of f_theta0 for large Ri numbers (1 - 10)
  -

* - (turb_vdiff_nml-ek_ep_ratio_stable)=
    **ek_ep_ratio_stable**
  - R
  - 3.0
  -
  - Ratio of TKE to TPE for large positive Ri (Mauritsen: {math}`1/(0.3\pm 1) - 1`)
  -

* - (turb_vdiff_nml-ek_ep_ratio_unstable)=
    **ek_ep_ratio_unstable**
  - R
  - 2.0
  -
  - Ratio of TKE to TPE for large negative Ri (Mauritsen: 1)
  -

* - (turb_vdiff_nml-c_f)=
    **c_f**
  - R
  - 0.185
  -
  - mixing length: coriolis term tuning parameter
  -

* - (turb_vdiff_nml-c_n)=
    **c_n**
  - R
  - 2.0
  -
  - mixing length: stability term tuning parameter
  -

* - (turb_vdiff_nml-wmc)=
    **wmc**
  - R
  - 0.5
  -
  - ratio of typical horizontal velocity to wstar at free convection
  -

* - (turb_vdiff_nml-fbl)=
    **fbl**
  - R
  - 3.0
  -
  - 1/fbl: fraction of BL height at which lmix hat its max
  -

* - (turb_vdiff_nml-lmix_max)=
    **lmix_max**
  - R
  - 150
  - m
  - maximum mixing length
  -
```



### 3D Smagorinsky Scheme

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (turb_vdiff_nml-km_min)=
    **km_min**
  - R
  - 0.001
  - Pa s
  - minimum mass weighted turbulent viscosity
  -

* - (turb_vdiff_nml-turb_prandtl)=
    **turb_prandtl**
  - R
  - 1/3
  -
  - Turbulent Prandtl number
  -

* - (turb_vdiff_nml-min_sfc_wind)=
    **min_sfc_wind**
  - R
  - 1.0
  - m/s
  - minimum surface wind speed in free-convection limit
  -
```



The limit fractions {math}`L` and decay constants {math}`D` for {math}`f_\tau` and {math}`f_\theta` are defined with respect to the ansatz {math}`f_\tau(\mathrm{Ri}) = f_\tau(0) \left( L + \frac{1-L}{1 + D \,\mathrm{Ri}} \right)`

Defined and used in: {{ '[src/namelists/mo_turb_vdiff_nml.f90]({}/src/namelists/mo_turb_vdiff_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_turbdiff_nml)=
## turbdiff_nml


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (turbdiff_nml-imode_turb)=
    **imode_turb**
  - I
  - 1
  -
  - Mode of solving the TKE equation for atmosph. layers:
    - 0: diagnostic equation
    - 1: prognostic equation (current version)
    - 2: prognostic equation (intrinsically positive definite).
  -

* - (turbdiff_nml-imode_tran)=
    **imode_tran**
  - I
  - 0
  -
  - Same as _imode_turb_ but only for the transfer layer.
  -

* - (turbdiff_nml-icldm_turb)=
    **icldm_turb**
  - I
  - 2
  -
  - Mode of water cloud representation in turbulence for atmosph. layers:
    - -1: ignoring cloud water completely (pure dry scheme)
    - 0: no clouds considered (all cloud water is evaporated)
    - 1: only grid scale condensation possible
    - 2: also sub grid (turbulent) condensation considered.
  -

* - (turbdiff_nml-icldm_tran)=
    **icldm_tran**
  - I
  - 2
  -
  - Same as _icldm_turb_ but only for the transfer layer.
  -

* - (turbdiff_nml-itype_wcld)=
    **itype_wcld**
  - I
  - 2
  -
  - Type of water cloud diagnosis within the turbulence scheme:
    - 1: employing a scheme based on relative humitidy
    - 2: employing a statistical saturation adjustment.
  - icldm_turb=2 or icldm_tran=2

* - (turbdiff_nml-q_crit)=
    **q_crit**
  - R
  - 1.6
  -
  - Critical value for normalized super-saturation.
  - itype_wcld=2

* - (turbdiff_nml-itype_sher)=
    **itype_sher**
  - I
  - 0
  -
  - Type of shear forcing used in turbulence:
    - 0: only vertical shear of horizontal wind
    - 1: previous plus horizontal shear correction
    - 2: previous plus shear from vertical velocity.
  -

* - (turbdiff_nml-ltkeshs)=
    **ltkeshs**
  - L
  - .TRUE.
  -
  - Consider TKE-production by separated horizontal shear eddies.
  -

* - (turbdiff_nml-imode_shshear)=
    **imode_shshear**
  - I
  - 2
  -
  - Mode of calculat. the separated horizontal shear mode:
    - 0: with a constant length scale and based on 3D-shear and incompressibility
    - 1: with a constant length scale and considering the trace constraint for the 2D-strain tensor
    - 2: with a Ri-number depend. length-scale correct. and the trace constraint for the 2D-strain tensor.
  - ltkeshs=.TRUE. and a_hshr > 0

* - (turbdiff_nml-ltkesso)=
    **ltkesso**
  - L
  - .TRUE.
  -
  - Consider TKE-production by sub grid SSO wakes.
  - inwp_sso = 1

* - (turbdiff_nml-imode_tkesso)=
    **imode_tkesso**
  - I
  - 1
  -
  - Mode of calculat. the SSO source term for TKE production:
    - 1: original implementation
    - 2: with Ri-number dependent reduction factor for Ri > 1
    - 3: as "2", but with additional reduction for dx < 2 km.
  - ltkesso=.TRUE.

* - (turbdiff_nml-ltkecon)=
    **ltkecon**
  - L
  - .FALSE.
  -
  - Consider TKE-production by sub grid convective plumes.
  - inwp_conv = 1

* - (turbdiff_nml-ltmpcor)=
    **ltmpcor**
  - L
  - .FALSE.
  -
  - Consider thermal TKE sources in enthalpy equation.
  -

* - (turbdiff_nml-lcpfluc)=
    **lcpfluc**
  - L
  - .FALSE.
  -
  - Consideration of fluctuations of the heat capacity of air.
  -

* - (turbdiff_nml-tur_len)=
    **tur_len**
  - R
  - 500.0
  - m
  - Asymptotic maximal turbulent distance (where {math}`\kappa \cdot` _tur_len_ is the integral turbulent master length-scale).
  -

* - (turbdiff_nml-pat_len)=
    **pat_len**
  - R
  - 100.0
  - m
  - Effective length scale of thermal surface patterns controlling TKE-production by sub grid kata/ana-batic circulations. In case of _pat_len_=0, this production is switched off.
  -

* - (turbdiff_nml-c_diff)=
    **c_diff**
  - R
  - 0.2
  - 1
  - Length scale factor for vertical diffusion of TKE. In case of _c_diff_=0, TKE is not diffused vertically.
  -

* - (turbdiff_nml-a_stab)=
    **a_stab**
  - R
  - 0.0
  - 1
  - Factor for stability correction of turbulent master length-scale. In case of _a_stab_=0, this turbulent length scale is not reduced for stable stratification.
  -

* - (turbdiff_nml-a_hshr)=
    **a_hshr**
  - R
  - 0.20
  - 1
  - Length scale factor for the separated horizontal shear mode. In case of _a_hshr_=0, this shear mode has no effect.
  - ltkeshs=.TRUE.

* - (turbdiff_nml-tkhmin)=
    **tkhmin**
  - R
  - 0.75
  - {math}`\mathrm{m^2/s}`
  - Basic minimum vertical diffusion coefficient for scalar properties like heat and moisture (being corrected by an empirical factor proportional to {math}`{Ri}^{-2/3}`).
  -

* - (turbdiff_nml-tkhmin_strat)=
    **tkhmin_strat**
  - R
  - 0.75
  - {math}`\mathrm{m^2/s}`
  - Enhanced value of _tkhmin_ valid for the stratosphere above 17.5 km (tropics above 22.5 km) (being corrected by an empirical factor proportional to {math}`{Ri}^{-1/3}`).
  -

* - (turbdiff_nml-tkmmin)=
    **tkmmin**
  - R
  - 0.75
  - {math}`\mathrm{m^2/s}`
  - Basic minimum vertical diffusion coefficient for momentum (being corrected by an empirical factor proportional to {math}`{Ri}^{-2/3}`).
  -

* - (turbdiff_nml-tkmmin_strat)=
    **tkmmin_strat**
  - R
  - 4
  - {math}`\mathrm{m^2/s}`
  - Enhanced value of _tkmmin_ valid for the stratosphere above 17.5 km (tropics above 22.5 km) (being corrected by an empirical factor proportional to {math}`{Ri}^{-1/3}`).
  -

* - (turbdiff_nml-imode_tkemini)=
    **imode_tkemini**
  - I
  - 1
  -
  - Mode of adapting q=SQRT(2*TKE) and the entire Turbulence-Model (TMod) to Lower Limits for Diff. Coeffs. (LLDCs), which are the "minimal vertical diffusion coefficients" _tkhmin_ or _tkmmin_ (or their stratopheric extensions):
    - 1: LLDCs treated as corrections of stability lengths without any further adaptation.
    - 2: TKE adapted to that part of LLDCs representing so far missing shear forcing, while the assumed part of LLDCs, representing missing drag-forces, has no feedback to the TMod.
    - 3: Tuned variant of "2" that is suitable for operational forecasts, particularly excluding the effect of stratospheric extensions on the TKE adaptation.
  - any (extended) LLDC is active

* - (turbdiff_nml-alpha0)=
    **alpha0**
  - R
  - 0.0123
  - 1
  - Standard Charnock parameter.
  -

* - (turbdiff_nml-alpha0_max)=
    **alpha0_max**
  - R
  - 0.0335
  - 1
  - Upper bound of velocity-dependent Charnock parameter. Setting this parameter to 0.0335 or higher values, implies unconstrained velocity dependence.
  -

* - (turbdiff_nml-alpha1)=
    **alpha1**
  - R
  - 0.75
  - 1
  - Scaling parameter for molecular roughness of ocean waves.
  -

* - (turbdiff_nml-imode_charpar)=
    **imode_charpar**
  - I
  - 2
  -
  - Options for specifying the Charnock parameter:
    - 1: constant at _alpha0_
    - 2: wind-speed dependent with maximum at _alpha0_max_
    - 3: as "2", but decreasing again at speeds above about 25 m/s in order to improve pressure-speed relationship in tropical cyclones.
  -

* - (turbdiff_nml-lconst_z0)=
    **lconst_z0**
  - L
  - .FALSE.
  -
  - TRUE: horizontally homogeneous roughness length z0.
  -

* - (turbdiff_nml-const_z0)=
    **const_z0**
  - R
  - 0.001
  - m
  - Value for horizontally homogeneous roughness length z0.
  - lconst_z0=.TRUE.

* - (turbdiff_nml-itype_synd)=
    **itype_synd**
  - I
  - 2
  -
  - Type of diagnostics of synoptic near surface variables:
    - 1: Considering the mean surface roughness of a grid box
    - 2: Considering a fictive surface roughness of a SYNOP lawn.
  -

* - (turbdiff_nml-rlam_heat)=
    **rlam_heat**
  - R
  - 10.0
  - 1
  - Scaling factor of the laminar boundary layer for scalars (heat and vapor). The larger rlam_heat, the larger is the laminar resistance.
  -

* - (turbdiff_nml-rat_lam)=
    **rat_lam**
  - R
  - 1.0
  - 1
  - Vapour/Heat ratio of laminar scaling factors (over land). The larger _rat_lam_, the larger is the laminar resistance for evaporation compared to sensible heat.
  -

* - (turbdiff_nml-rat_sea)=
    **rat_sea**
  - R
  - 0.8
  - 1
  - Sea/Land ratio of laminar scaling factors for scalars (heat and vapor). The larger _rat_sea_, the larger is the laminar resistance for a sea surface compared to a land surface.
  -

* - (turbdiff_nml-rat_glac)=
    **rat_glac**
  - R
  - 3.0
  - 1
  - Glacier/Land ratio of laminar scaling factors for scalars (heat and vapor). The larger _rat_glac_, the larger is the laminar resistance over glaciers compared to other land surfaces.
  -

* - (turbdiff_nml-tkesmot)=
    **tkesmot**
  - R
  - 0.15
  - 1
  - Time smoothing factor within {math}`[0, 1]` for TKE. In case of _tkesmot_=0, no smoothing is active.
  -

* - (turbdiff_nml-frcsmot)=
    **frcsmot**
  - R
  - 0.0
  - 1
  - Vertical smoothing factor within {math}`[0, 1]` for TKE forcing terms. In case of _frcsmot_=0, no smoothing is active.
  -

* - (turbdiff_nml-imode_frcsmot)=
    **imode_frcsmot**
  - I
  - 1
  -
  - 1: Apply vertical smoothing uniformly over the globe
    - 2: Restrict vertical smoothing to the tropics (reduces the moist bias in the tropics while avoiding adverse effects on NWP skill scores in the extratropics).
  - frcsmot > 0

* - (turbdiff_nml-imode_snowsmot)=
    **imode_snowsmot**
  - I
  - 1
  -
  - Mode to treating the aerodynamic surface-smoothing by snow:
    - 0: no smoothing active at all
    - 1: no impact on SAI, but full smoothing of land-use R-length
    - 2: "1", but with full smoothing of SAI: full smoothing of R-length and SAI
    - 3: dynamical smoothing of R-length and SAI dependent on snow- and R-height.
  - itype_z0 {math}`\ge` 2

* - (turbdiff_nml-lsflcnd)=
    **lsflcnd**
  - L
  - .TRUE.
  -
  - Use lower flux condition for vertical diffusion calculation (TRUE) instead of a lower concentration condition (FALSE).
  -

* - (turbdiff_nml-lfreeslip)=
    **lfreeslip**
  - L
  - .FALSE.
  -
  - Use a free-slip lower boundary condition (TRUE), i.e. neither momentum nor heat/moisture fluxes (use for idealized runs only!).
  -

* - (turbdiff_nml-impl_s)=
    **impl_s**
  - R
  - 1.20
  - 1
  - Implicit weight near the surface (maximal value).
  -

* - (turbdiff_nml-impl_t)=
    **impl_t**
  - R
  - 0.75
  - 1
  - Implicit weight near top of the atmosphere (minimal value).
  -

* - (turbdiff_nml-ldiff_qi)=
    **ldiff_qi**
  - L
  - .FALSE.
  -
  - Turbulent diffusion of cloud ice (TRUE).
  -
```



Defined and used in: {{ '[src/namelists/mo_turbdiff_nml.f90]({}/src/namelists/mo_turbdiff_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_upatmo_nml)=
## upatmo_nml

### Extrapolation to determine the inital state of the upper atmosphere

Scope: itype_vert_expol = 2

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (upatmo_nml-expol_start_height)=
    **expol_start_height**
  - R
  - 70000
  - m
  - Height above which extrapolation of initial data starts.
  -

* - (upatmo_nml-expol_blending_scale)=
    **expol_blending_scale**
  - R
  - 10000
  - m
  - Vertical distance above expol_start_height within which blending of linearly extrapolated state and climatological state takes place.
  -

* - (upatmo_nml-expol_vn_decay_scale)=
    **expol_vn_decay_scale**
  - R
  - 10000
  - m
  - Scale height of vertically exponentially decaying factor multiplied to the extrapolated horizontal wind (to alleviate stability-endangering wind magnitudes).
  -

* - (upatmo_nml-expol_temp_infty)=
    **expol_temp_infty**
  - R
  - 400
  - K
  - Exospheric mean reference temperature of the climatology for the extrapolation blending.
  -

* - (upatmo_nml-lexpol_sanitycheck)=
    **lexpol_sanitycheck**
  - L
  - .FALSE.
  -
  - .TRUE.: Apply some rudimentary sanity check to the extrapolated atmospheric state in the region above expol_start_height (e.g., temperature values everywhere > 0). (Please, apply with care, since it is computationally relatively expensive.)
  -
```



### Upper-atmosphere physics

Scope: ([iforcing](run_nml-iforcing) = 2 (AES)) or ([iforcing](run_nml-iforcing) = 3 (NWP) & lupatmo_phy = .TRUE.)

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (upatmo_nml-orbit_type)=
    **orbit_type**
  - I
  - 1
  -
  - Orbit model for upper-atmosphere radiation (compare l_orbvsop87):
    - 1: vsop87 {math}`\rightarrow` standard and accurate model
    - 2: kepler {math}`\rightarrow` simple model appropriate for idealized work
  -

* - (upatmo_nml-solvar_type)=
    **solvar_type**
  - I
  - 1
  -
  - Solar activity:
    - 1: normal
    - 2: low
    - 3: high
  -

* - (upatmo_nml-solvar_data)=
    **solvar_data**
  - I
  - 2
  -
  - Data set for solar activity:
    - 1: G. Rottman data
    - 2: J. Lean data
  -

* - (upatmo_nml-solcyc_type)=
    **solcyc_type**
  - I
  - 2
  -
  - Solar cycle:
    - 1: standard cycle
    - 2: 27-day cycle
  -

* - (upatmo_nml-nwp_grp_)=
    **nwp_grp_\<groupname>%...**
  -
  -
  -
  - Configuration of the upper-atmosphere process groups under NWP-forcing (compare time control of processes in [aes_phy_nml](ref_buildrun_nml_aes_phy_nml)):
    \<groupname> = imf: ion drag, molecular diffusion and frictional heating
    \<groupname> = rad: radiation and chemical heating
  - [iforcing](run_nml-iforcing) = 3
    lupatmo_phy = .TRUE.

* - (upatmo_nml-nwp_grp_-imode)=
    **...imode**
  - I(max_dom)
  - 1
  -
  - Group mode:
    - 0: all processes clustered in the group \<groupname> are switched off
    - 1: all processes are switched on
    - 2: all processes run in offline-mode, i.e. tendencies are computed, but not coupled to the dynamics
    Example of usage for multi-domain applications:
      - set nwp_grp_imf%imode = 1 to switch on the IMF-group for all domains (default)
      - set nwp_grp_rad%imode = 1,1,0 to switch on the RAD-group for domain 1 and 2, but to switch it off for domain 3
    Please note: if imode = 1 or 2 for a domain, but lupatmo_phy = .FALSE. for this domain, imode is set to 0 and the group is switched off.
  -

* - (upatmo_nml-nwp_grp_-dt)=
    **...dt**
  - R(max_dom)
  - 300.0{math}`_{\text{imf}}`,
    600.0{math}`_{\text{rad}}`
  - s
  - Tendency update period. New tendencies from all processes of a group are computed every dt (temperature, wind and water vapor tendencies in case of IMF, and temperature tendencies in case of RAD). Please note: internal processing will round dt to the next multiple of the domain-adjusted value of dtime, which in turn might have been rescaled, if grid_rescale_factor {math}`\neq` 1. In case of a domain-wise assignment in a multi-domain application, dt(1) {math}`\geq` dt(2) {math}`\geq`
    ldots is required.
  -

* - (upatmo_nml-nwp_grp_-t_start)=
    **...t_start**
  - C
  - " "
  -
  - Tendencies from all processes of a group are computed within the time interval [t_start, t_end]. Outside this interval the tendencies are set to zero. Format as for ini_datetime_string, e.g. nwp_grp_imf%t_start = "2008-09-01T00:00:00Z". Empty strings will be replaced by the simulation start and/or end date and time of the domain. t_start and t_end apply to all domains, no domain-wise specification possible!
  -

* - (upatmo_nml-nwp_grp_-t_end)=
    **...t_end**
  - C
  - " "
  -
  - See t_start
  -

* - (upatmo_nml-nwp_grp_-start_height)=
    **...start_height**
  - R
  - -999.0
  - m
  - All processes of a group compute tendencies above start_height. Below start_height the processes are inactive and all tendencies are set to zero. A negative value means that the default start heights of each process, listed in {{ '[src/upper_atmosphere/mo_upatmo_impl_const.f90]({}/src/upper_atmosphere/mo_upatmo_impl_const.f90)'.format(base_url) }}: startHeightDef, are applied. Please note: start_height applies to all domains. If it is above the top of one domain, the group is switched off for that domain (imode(idom) is set to 0).
  -

* - (upatmo_nml-nwp_gas_)=
    **nwp_gas_\<gasname>%...**
  -
  -
  -
  - Configuration of the radiatively active gases in the upper atmosphere under NWP-forcing (compare [radiation_nml](ref_buildrun_nml_radiation_nml) and [aes_rad_nml](ref_buildrun_nml_aes_rad_nml)):
    \<gasname> = o3: ozone ({math}`\text{O}_3`)
    \<gasname> = o2: dioxygen ({math}`\text{O}_2`)
    \<gasname> = o: atomic oxygen ({math}`\text{O}`)
    \<gasname> = co2: carbon dioxide ({math}`\text{CO}_2`)
    \<gasname> = no: nitric oxide ({math}`\text{NO}`)
    (Dinitrogen ({math}`\text{N}_2`) is determined diagnostically.)
  - [iforcing](run_nml-iforcing) = 3
    lupatmo_phy = .TRUE.
    nwp_grp_rad%imode  > 0

* - (upatmo_nml-nwp_gas_-imode)=
    **...imode**
  - I
  - 2
  -
  - Gas mode (comparable, but generally not identical to the irad_\<gasname> in [radiation_nml](ref_buildrun_nml_radiation_nml) and [aes_rad_nml](ref_buildrun_nml_aes_rad_nml)).
    - 0: zero gas concentration
    - 1: constant gas concentration (independent of space and time), specified via nwp_gas_\<gasname>%vmr
    - 2: external data; meridionally, vertically and monthly varying gas concentrations are read from a file with name nwp_extdat_gases%filename
  -

* - (upatmo_nml-nwp_gas_-vmr)=
    **...vmr**
  - R
  - 0.0
  - {math}`\text{m}^3/\text{m}^3`
  - Constant volume mixing ratio for a radiatively active gas.
  - nwp_gas_\<gasname>%imode = 1

* - (upatmo_nml-nwp_gas_-fscale)=
    **...fscale**
  - R
  - 1.0
  -
  - Scaling factor the gas concentration in each grid cell is multiplied with.
  - nwp_gas_\<gasname>%imode > 0

* - (upatmo_nml-nwp_extdat_)=
    **nwp_extdat_\<extdatname>%...**
  -
  -
  -
  - Configuration of the external upper-atmosphere data:
    \<extdatname> = gases: concentrations of the radiatively active gases
    \<extdatname> = chemheat: temperature tendencies from chemical heating
    Please note: the standard NWP physics use other external gas data (e.g., for ozone)!
  - nwp_grp_rad%imode  > 0

* - (upatmo_nml-nwp_extdat_-dt)=
    **...dt**
  - R
  - 86400.0
  - s
  - Update period for the time interpolation of the external data. Currently, the external data provide monthly mean values. In order to avoid too strong jumps in the transition from one month to the next, the parameters are "smoothed" in time by a linear interpolation that is computed every dt. A value of the order of a day should be entirely sufficient for this purpose.
  -

* - (upatmo_nml-nwp_extdat_-filename)=
    **...filename**
  - C
  - "upatmo_gases_chemheat.nc"
  -
  - Name of the file containing the external data. The file of the default name can be found in the folder `data/`, to which a link has to be set in the run script, following the typical examples of lrtm_filename and cldopt_filename. May contain the keyword `<path>` which will be substituted by model_base_dir (e.g., `<path>upatmo_gases_chemheat.nc`). Please note: if you would like to use other external data files, their data structure has to follow _exactly_ the data structure of `data/upatmo_gases_chemheat.nc` (variable and dimension names and units, zonally averaged monthly mean gas concentrations on pressure levels, zonally averaged monthly mean temperature tendencies from chemical heating on geometric height levels etc.). Any other structure cannot be processed for the time being!
  -
```



Defined and used in: {{ '[src/namelists/mo_upatmo_nml.f90]({}/src/namelists/mo_upatmo_nml.f90)'.format(base_url) }}


### Some notes on the output of upper-atmosphere-specific variables (under NWP-forcing)

An output of upper-atmosphere fields is only possible,
if upper-atmosphere physics are switched on.

Upper-atmosphere fields cannot be output in the GRIB format
([filetype](output_nml-filetype) = 2).

Upper-atmosphere fields entered on [output_nml](ref_buildrun_nml_output_nml): m/h/pl_varlist
need the prefix "upatmo_".

The following fields can be output, if {math}`\ldots`

{math}`\ldots` lupatmo_phy = .TRUE.:

```{list-table}
:header-rows: 1

* - Parameter
  - Description

* - (upatmo_nml-upatmo_mdry)=
    **upatmo_mdry**
  - Mass of dry air

* - (upatmo_nml-upatmo_amd)=
    **upatmo_amd**
  - Molar mass of dry air

* - (upatmo_nml-upatmo_cpair)=
    **upatmo_cpair**
  - Heat capacity of (moist) air at constant pressure

* - (upatmo_nml-upatmo_grav)=
    **upatmo_grav**
  - Gravitational acceleration of Earth
```



{math}`\ldots` lupatmo_phy = .TRUE. & nwp_grp_rad%imode > 0:


```{list-table}
:header-rows: 1

* - Parameter
  - Description

* - (upatmo_nml-upatmo_sclrlw)=
    **upatmo_sclrlw**
  - Scaling factor for standard long-wave radiation heating rate from radiative processes

* -
  - out of local thermodynamic equilibrium

* - (upatmo_nml-upatmo_effrsw)=
    **upatmo_effrsw**
  - Efficiency factor for standard short-wave radiation heating rate from chemical heating

* - (upatmo_nml-upatmo_o3)=
    **upatmo_o3**
  - Mass mixing ratio of ozone (member of `group:upatmo_rad_gases`)

* - (upatmo_nml-upatmo_o2)=
    **upatmo_o2**
  - Mass mixing ratio of dioxygen (member of `group:upatmo_rad_gases`)

* - (upatmo_nml-upatmo_o)=
    **upatmo_o**
  - Mass mixing ratio of atomic oxygen (member of `group:upatmo_rad_gases`)

* - (upatmo_nml-upatmo_co2)=
    **upatmo_co2**
  - Mass mixing ratio of carbon dioxide (member of `group:upatmo_rad_gases`)

* - (upatmo_nml-upatmo_no)=
    **upatmo_no**
  - Mass mixing ratio of nitric oxide (member of `group:upatmo_rad_gases`)

* - (upatmo_nml-upatmo_n2)=
    **upatmo_n2**
  - Mass mixing ratio of dinitrogen (member of `group:upatmo_rad_gases`)

* - (upatmo_nml-upatmo_ddt_temp_srbc)=
    **upatmo_ddt_temp_srbc**
  - Temperature tendency due to absorbtion by O2 in Schumann-Runge band and continuum (member of `group:upatmo_tendencies`)

* - (upatmo_nml-upatmo_ddt_temp_nlte)=
    **upatmo_ddt_temp_nlte**
  - Temperature tendency due to radiative processes out of local thermodynamic equilibrium  (member of `group:upatmo_tendencies`)

* - (upatmo_nml-upatmo_ddt_temp_euv)=
    **upatmo_ddt_temp_euv**
  - Temperature tendency due to heating from extreme ultraviolet radiation  (member of `group:upatmo_tendencies`)

* - (upatmo_nml-upatmo_ddt_temp_no)=
    **upatmo_ddt_temp_no**
  - Temperature tendency due to NO heating at near infrared (member of `group:upatmo_tendencies`)

* - (upatmo_nml-upatmo_ddt_temp_chemheat)=
    **upatmo_ddt_temp_chemheat**
  - Temperature tendency due to chemical heating (member of `group:upatmo_tendencies`)
```



{math}`\ldots` lupatmo_phy = .TRUE. & nwp_grp_imf%imode > 0:


```{list-table}
:header-rows: 1

* - Parameter
  - Description

* - (upatmo_nml-upatmo_ddt_temp_vdfmol)=
    **upatmo_ddt_temp_vdfmol**
  - Temperature tendency due to molecular diffusion (member of `group:upatmo_tendencies`)

* - (upatmo_nml-upatmo_ddt_temp_fric)=
    **upatmo_ddt_temp_fric**
  - Temperature tendency due to frictional heating (member of `group:upatmo_tendencies`)

* - (upatmo_nml-upatmo_ddt_temp_joule)=
    **upatmo_ddt_temp_joule**
  - Temperature tendency due to Joule heating from ion drag (member of `group:upatmo_tendencies`)

* - (upatmo_nml-upatmo_ddt_u_vdfmol)=
    **upatmo_ddt_u_vdfmol**
  - Zonal component of wind tendency due to molecular diffusion (member of `group:upatmo_tendencies`)

* - (upatmo_nml-upatmo_ddt_v_vdfmol)=
    **upatmo_ddt_v_vdfmol**
  - Meridionl component of wind tendency due to molecular diffusion (member of `group:upatmo_tendencies`)

* - (upatmo_nml-upatmo_ddt_u_iondrag)=
    **upatmo_ddt_u_iondrag**
  - Zonal component of wind tendency due to ion drag (member of `group:upatmo_tendencies`)

* - (upatmo_nml-upatmo_ddt_v_iondrag)=
    **upatmo_ddt_v_iondrag**
  - Meridionl component of wind tendency due to ion drag (member of `group:upatmo_tendencies`)

* - (upatmo_nml-upatmo_ddt_qv_vdfmol)=
    **upatmo_ddt_qv_vdfmol**
  - Tendency of specific humidity due to molecular diffusion (member of `group:upatmo_tendencies`)
```



(ref_buildrun_nml_vertical_levels_nml)=
## vertical_levels_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (vertical_levels_nml-dz_full_level)=
    **dz_full_level**
  -
  -
  -
  -
  -

* - (vertical_levels_nml-zlevels)=
    **zlevels**
  -
  -
  -
  -
  -
```


Defined and used in: {{ '[src/io/icon_output_model/mo_read_icon_output_namelists.f90]({}/src/io/icon_output_model/mo_read_icon_output_namelists.f90)'.format(base_url) }}

# Ocean-specific namelist parameters


(ref_buildrun_nml_ocean_diagnostics_nml)=
## ocean_diagnostics_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ocean_diagnostics_nml-Green_tracer_width)=
    **Green_tracer_width**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-age_tracer_inv_relax_time)=
    **age_tracer_inv_relax_time**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-agulhas)=
    **agulhas**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-agulhas_long)=
    **agulhas_long**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-agulhas_longer)=
    **agulhas_longer**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-barentsOpening)=
    **barentsOpening**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-beringStrait)=
    **beringStrait**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-check_total_volume)=
    **check_total_volume**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-denmark_strait)=
    **denmark_strait**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-diagnose_age)=
    **diagnose_age**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-diagnose_for_heat_content)=
    **diagnose_for_heat_content**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-diagnose_for_horizontalVelocity)=
    **diagnose_for_horizontalVelocity**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-diagnose_for_tendencies)=
    **diagnose_for_tendencies**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-diagnose_green)=
    **diagnose_green**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-diagnostics_level)=
    **diagnostics_level**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-do_ts_budget)=
    **do_ts_budget**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-drake_passage)=
    **drake_passage**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-eddydiag)=
    **eddydiag**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-florida_strait)=
    **florida_strait**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-framStrait)=
    **framStrait**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-gibraltar)=
    **gibraltar**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-greenStartDate)=
    **greenStartDate**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-greenStopDate)=
    **greenStopDate**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-indonesian_throughflow)=
    **indonesian_throughflow**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-l_relaxage_ice)=
    **l_relaxage_ice**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-mode_layers)=
    **mode_layers**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-mozambique)=
    **mozambique**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-n_dlev)=
    **n_dlev**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-rho_lev_in)=
    **rho_lev_in**
  - REAL
  -
  -
  -
  -

* - (ocean_diagnostics_nml-scotland_iceland)=
    **scotland_iceland**
  -
  -
  -
  -
  -

* - (ocean_diagnostics_nml-use_layers)=
    **use_layers**
  -
  -
  -
  -
  -
```


Defined and used in: {{ '[src/ocean/config/mo_ocean_nml.f90]({}/src/ocean/config/mo_ocean_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_ocean_dynamics_nml)=
## ocean_dynamics_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ocean_dynamics_nml-HorizonatlVelocity_VerticalAdvection_form)=
    **HorizonatlVelocity_VerticalAdvection_form**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-KineticEnergy_type)=
    **KineticEnergy_type**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-MASS_MATRIX_INVERSION_TYPE)=
    **MASS_MATRIX_INVERSION_TYPE**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-MassMatrix_solver_tolerance)=
    **MassMatrix_solver_tolerance**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-NonlinearCoriolis_type)=
    **NonlinearCoriolis_type**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-ab_beta)=
    **ab_beta**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-ab_const)=
    **ab_const**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-ab_gam)=
    **ab_gam**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-basin_center_lat)=
    **basin_center_lat**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-basin_center_lon)=
    **basin_center_lon**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-basin_height_deg)=
    **basin_height_deg**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-basin_width_deg)=
    **basin_width_deg**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-cfl_check)=
    **cfl_check**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-cfl_stop_on_violation)=
    **cfl_stop_on_violation**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-cfl_threshold)=
    **cfl_threshold**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-cfl_write)=
    **cfl_write**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-coriolis_fplane_latitude)=
    **coriolis_fplane_latitude**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-coriolis_type)=
    **coriolis_type**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-createSolverMatrix)=
    **createSolverMatrix**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-dhdtw_abort)=
    **dhdtw_abort**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-discretization_scheme)=
    **discretization_scheme**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-dzlev_m)=
    **dzlev_m**
  - REAL
  -
  -
  -
  -

* - (ocean_dynamics_nml-fast_performance_level)=
    **fast_performance_level**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-i_bc_veloc_bot)=
    **i_bc_veloc_bot**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-i_bc_veloc_lateral)=
    **i_bc_veloc_lateral**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-i_bc_veloc_top)=
    **i_bc_veloc_top**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-iswm_oce)=
    **iswm_oce**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-l_RIGID_LID)=
    **l_RIGID_LID**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-l_edge_based)=
    **l_edge_based**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-l_inverse_flip_flop)=
    **l_inverse_flip_flop**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-l_lhs_direct)=
    **l_lhs_direct**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-l_max_bottom)=
    **l_max_bottom**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-l_partial_cells)=
    **l_partial_cells**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-l_solver_compare)=
    **l_solver_compare**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-lviscous)=
    **lviscous**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-minVerticalLevels)=
    **minVerticalLevels**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-n_zlev)=
    **n_zlev**
  - INTEGER
  - -1
  -
  -
  -

* - (ocean_dynamics_nml-press_grad_type)=
    **press_grad_type**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-select_lhs)=
    **select_lhs**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-select_solver)=
    **select_solver**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-select_transfer)=
    **select_transfer**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-solver_FirstGuess)=
    **solver_FirstGuess**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-solver_comp_nsteps)=
    **solver_comp_nsteps**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-solver_max_iter_per_restart)=
    **solver_max_iter_per_restart**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-solver_max_iter_per_restart_sp)=
    **solver_max_iter_per_restart_sp**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-solver_max_restart_iterations)=
    **solver_max_restart_iterations**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-solver_tolerance)=
    **solver_tolerance**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-solver_tolerance_comp)=
    **solver_tolerance_comp**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-solver_tolerance_sp)=
    **solver_tolerance_sp**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-surface_module)=
    **surface_module**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-threshold_vn)=
    **threshold_vn**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-use_absolute_solver_tolerance)=
    **use_absolute_solver_tolerance**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-use_continuity_correction)=
    **use_continuity_correction**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-use_smooth_ocean_boundary)=
    **use_smooth_ocean_boundary**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-use_ssh_in_momentum_eq)=
    **use_ssh_in_momentum_eq**
  -
  -
  -
  -
  -

* - (ocean_dynamics_nml-vert_cor_type)=
    **vert_cor_type**
  -
  -
  -
  -
  -
```


Defined and used in: {{ '[src/ocean/config/mo_ocean_nml.f90]({}/src/ocean/config/mo_ocean_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_ocean_forcing_nml)=
## ocean_forcing_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ocean_forcing_nml-forcing_center)=
    **forcing_center**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_enable_freshwater)=
    **forcing_enable_freshwater**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_fluxes_type)=
    **forcing_fluxes_type**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_frequency)=
    **forcing_frequency**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_set_runoff_to_zero)=
    **forcing_set_runoff_to_zero**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_timescale)=
    **forcing_timescale**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_windStress_u_amplitude)=
    **forcing_windStress_u_amplitude**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_windStress_v_amplitude)=
    **forcing_windStress_v_amplitude**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_windspeed_amplitude)=
    **forcing_windspeed_amplitude**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_windspeed_type)=
    **forcing_windspeed_type**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_windstress_merid_waveno)=
    **forcing_windstress_merid_waveno**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_windstress_u_type)=
    **forcing_windstress_u_type**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_windstress_v_type)=
    **forcing_windstress_v_type**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_windstress_zonalWavePhase)=
    **forcing_windstress_zonalWavePhase**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-forcing_windstress_zonal_waveno)=
    **forcing_windstress_zonal_waveno**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-heatflux_forcing_on_sst)=
    **heatflux_forcing_on_sst**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-lcheck_salt_content)=
    **lcheck_salt_content**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-lfix_salt_content)=
    **lfix_salt_content**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-lfwflux_enters_with_sst)=
    **lfwflux_enters_with_sst**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-sw_scaling_factor)=
    **sw_scaling_factor**
  -
  -
  -
  -
  -

* - (ocean_forcing_nml-zero_freshwater_flux)=
    **zero_freshwater_flux**
  -
  -
  -
  -
  -
```


Defined and used in: {{ '[src/ocean/config/mo_ocean_nml.f90]({}/src/ocean/config/mo_ocean_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_ocean_GentMcWilliamsRedi_nml)=
## ocean_GentMcWilliamsRedi_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ocean_GentMcWilliamsRedi_nml-BOLUS_VELOCITY_DIAGNOSTIC)=
    **BOLUS_VELOCITY_DIAGNOSTIC**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-GMREDI_COMBINED_DIAGNOSTIC)=
    **GMREDI_COMBINED_DIAGNOSTIC**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-GMRedi_configuration)=
    **GMRedi_configuration**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-GMRedi_usesRelativeMaxSlopes)=
    **GMRedi_usesRelativeMaxSlopes**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-INCLUDE_SLOPE_SQUARED_IMPLICIT)=
    **INCLUDE_SLOPE_SQUARED_IMPLICIT**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-Nmin)=
    **Nmin**
  - REAL
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-REVERT_VERTICAL_RECON_AND_TRANSPOSED)=
    **REVERT_VERTICAL_RECON_AND_TRANSPOSED**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-RossbyRadius_max)=
    **RossbyRadius_max**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-RossbyRadius_min)=
    **RossbyRadius_min**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-SLOPE_CALC_VIA_TEMPERTURE_SALINITY)=
    **SLOPE_CALC_VIA_TEMPERTURE_SALINITY**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-SWITCH_OFF_TAPERING)=
    **SWITCH_OFF_TAPERING**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-SWITCH_ON_REDI_BALANCE_DIAGONSTIC)=
    **SWITCH_ON_REDI_BALANCE_DIAGONSTIC**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-SWITCH_ON_TAPERING_HORIZONTAL_DIFFUSION)=
    **SWITCH_ON_TAPERING_HORIZONTAL_DIFFUSION**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-S_critical)=
    **S_critical**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-S_d)=
    **S_d**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-S_max)=
    **S_max**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-TEST_MODE_GM_ONLY)=
    **TEST_MODE_GM_ONLY**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-TEST_MODE_REDI_ONLY)=
    **TEST_MODE_REDI_ONLY**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-c_speed)=
    **c_speed**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-k_tracer_GM_kappa_parameter)=
    **k_tracer_GM_kappa_parameter**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-k_tracer_dianeutral_parameter)=
    **k_tracer_dianeutral_parameter**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-k_tracer_isoneutral_parameter)=
    **k_tracer_isoneutral_parameter**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-lvertical_GM)=
    **lvertical_GM**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-switch_off_diagonal_vert_expl)=
    **switch_off_diagonal_vert_expl**
  -
  -
  -
  -
  -

* - (ocean_GentMcWilliamsRedi_nml-tapering_scheme)=
    **tapering_scheme**
  -
  -
  -
  -
  -
```


Defined and used in: {{ '[src/ocean/config/mo_ocean_nml.f90]({}/src/ocean/config/mo_ocean_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_ocean_horizontal_diffusion_nml)=
## ocean_horizontal_diffusion_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ocean_horizontal_diffusion_nml-BiharmonicDiv_weight)=
    **BiharmonicDiv_weight**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-BiharmonicViscosity_background)=
    **BiharmonicViscosity_background**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-BiharmonicViscosity_reference)=
    **BiharmonicViscosity_reference**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-BiharmonicViscosity_scaling)=
    **BiharmonicViscosity_scaling**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-BiharmonicVort_weight)=
    **BiharmonicVort_weight**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-HarmonicDiv_weight)=
    **HarmonicDiv_weight**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-HarmonicViscosity_background)=
    **HarmonicViscosity_background**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-HarmonicViscosity_reference)=
    **HarmonicViscosity_reference**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-HarmonicViscosity_scaling)=
    **HarmonicViscosity_scaling**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-HarmonicVort_weight)=
    **HarmonicVort_weight**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-HorizontalViscosity_SmoothIterations)=
    **HorizontalViscosity_SmoothIterations**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-HorizontalViscosity_SpatialSmoothFactor)=
    **HorizontalViscosity_SpatialSmoothFactor**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-LeithBiharmonicViscosity_background)=
    **LeithBiharmonicViscosity_background**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-LeithBiharmonicViscosity_reference)=
    **LeithBiharmonicViscosity_reference**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-LeithBiharmonicViscosity_scaling)=
    **LeithBiharmonicViscosity_scaling**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-LeithClosure_form)=
    **LeithClosure_form**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-LeithClosure_order)=
    **LeithClosure_order**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-LeithHarmonicViscosity_background)=
    **LeithHarmonicViscosity_background**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-LeithHarmonicViscosity_reference)=
    **LeithHarmonicViscosity_reference**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-LeithHarmonicViscosity_scaling)=
    **LeithHarmonicViscosity_scaling**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-LeithViscosity_SmoothIterations)=
    **LeithViscosity_SmoothIterations**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-LeithViscosity_SpatialSmoothFactor)=
    **LeithViscosity_SpatialSmoothFactor**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-N_POINTS_IN_MUNK_LAYER)=
    **N_POINTS_IN_MUNK_LAYER**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-Salinity_HorizontalDiffusion_Background)=
    **Salinity_HorizontalDiffusion_Background**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-Salinity_HorizontalDiffusion_Reference)=
    **Salinity_HorizontalDiffusion_Reference**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-Temperature_HorizontalDiffusion_Background)=
    **Temperature_HorizontalDiffusion_Background**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-Temperature_HorizontalDiffusion_Reference)=
    **Temperature_HorizontalDiffusion_Reference**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-TracerDiffusion_LeithWeight)=
    **TracerDiffusion_LeithWeight**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-TracerHorizontalDiffusion_scaling)=
    **TracerHorizontalDiffusion_scaling**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-Tracer_HorizontalDiffusion_PTP_coeff)=
    **Tracer_HorizontalDiffusion_PTP_coeff**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-VelocityDiffusion_order)=
    **VelocityDiffusion_order**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-laplacian_form)=
    **laplacian_form**
  -
  -
  -
  -
  -

* - (ocean_horizontal_diffusion_nml-max_turbulenece_TracerDiffusion)=
    **max_turbulenece_TracerDiffusion**
  -
  -
  -
  -
  -
```


Defined and used in: {{ '[src/ocean/config/mo_ocean_nml.f90]({}/src/ocean/config/mo_ocean_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_ocean_initialConditions_nml)=
## ocean_initialConditions_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ocean_initialConditions_nml-InitialState_InputFileName)=
    **InitialState_InputFileName**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-ana_filename)=
    **ana_filename**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-ana_varnames_map_file_oce)=
    **ana_varnames_map_file_oce**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-check_ana_oce_)=
    **check_ana_oce(:)%list**
  - CHARACTER
  - ''
  -
  -
  -

* - (ocean_initialConditions_nml-check_fg_oce_)=
    **check_fg_oce(:)%list**
  - CHARACTER
  - ''
  -
  -
  -

* - (ocean_initialConditions_nml-dt_ana_oce)=
    **dt_ana_oce**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-dt_iau_oce)=
    **dt_iau_oce**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-dt_shift_oce)=
    **dt_shift_oce**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-fg_filename)=
    **fg_filename**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-init_mode_oce)=
    **init_mode_oce**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_age_type)=
    **initial_age_type**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_green_type)=
    **initial_green_type**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_perturbation_max_ratio)=
    **initial_perturbation_max_ratio**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_perturbation_waveNumber)=
    **initial_perturbation_waveNumber**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_salinity_bottom)=
    **initial_salinity_bottom**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_salinity_top)=
    **initial_salinity_top**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_salinity_type)=
    **initial_salinity_type**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_sst_type)=
    **initial_sst_type**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_temperature_VerticalGradient)=
    **initial_temperature_VerticalGradient**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_temperature_bottom)=
    **initial_temperature_bottom**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_temperature_north)=
    **initial_temperature_north**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_temperature_scale_depth)=
    **initial_temperature_scale_depth**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_temperature_shift)=
    **initial_temperature_shift**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_temperature_south)=
    **initial_temperature_south**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_temperature_top)=
    **initial_temperature_top**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_temperature_type)=
    **initial_temperature_type**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_velocity_amplitude)=
    **initial_velocity_amplitude**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initial_velocity_type)=
    **initial_velocity_type**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-initialize_fromRestart)=
    **initialize_fromRestart**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-lconsistency_checks_oce)=
    **lconsistency_checks_oce**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-lread_ana_oce)=
    **lread_ana_oce**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-sea_surface_height_type)=
    **sea_surface_height_type**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-smooth_initial_height_iterations)=
    **smooth_initial_height_iterations**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-smooth_initial_height_weights)=
    **smooth_initial_height_weights**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-smooth_initial_salinity_iterations)=
    **smooth_initial_salinity_iterations**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-smooth_initial_salinity_weights)=
    **smooth_initial_salinity_weights**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-smooth_initial_temperature_iterations)=
    **smooth_initial_temperature_iterations**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-smooth_initial_temperature_weights)=
    **smooth_initial_temperature_weights**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-smooth_initial_velocity_iterations)=
    **smooth_initial_velocity_iterations**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-smooth_initial_velocity_weights)=
    **smooth_initial_velocity_weights**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-topography_height_reference)=
    **topography_height_reference**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-topography_type)=
    **topography_type**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-type_iau_wgt_oce)=
    **type_iau_wgt_oce**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-use_file_initialConditions)=
    **use_file_initialConditions**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-use_fillValue)=
    **use_fillValue**
  -
  -
  -
  -
  -

* - (ocean_initialConditions_nml-use_initicono)=
    **use_initicono**
  -
  -
  -
  -
  -
```


Defined and used in: {{ '[src/ocean/config/mo_ocean_nml.f90]({}/src/ocean/config/mo_ocean_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_ocean_physics_nml)=
## ocean_physics_nml


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ocean_physics_nml-i_sea_ice)=
    **i_sea_ice**
  - I
  - 1
  -
  - 0: No sea ice, 1: Include sea ice
  -

* -
  -
  -
  -
  - .FALSE.: compute drag only
  -

* - (ocean_physics_nml-richardson_factor_tracer)=
    **richardson_factor_tracer**
  - I
  - 0.5e-5
  - m/s
  -
  -

* - (ocean_physics_nml-richardson_factor_veloc)=
    **richardson_factor_veloc**
  - I
  - 0.5e-5
  - m/s
  -
  -

* - (ocean_physics_nml-l_constant_mixing)=
    **l_constant_mixing**
  - L
  - .FALSE.
  -
  -
  -
```



(ref_buildrun_nml_ocean_tracer_transport_nml)=
## ocean_tracer_transport_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ocean_tracer_transport_nml-fct_high_order_flux)=
    **fct_high_order_flux**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-fct_limiter_horz)=
    **fct_limiter_horz**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-fct_low_order_flux)=
    **fct_low_order_flux**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-flux_calculation_horz)=
    **flux_calculation_horz**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-flux_calculation_vert)=
    **flux_calculation_vert**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-l_GRADIENT_RECONSTRUCTION)=
    **l_GRADIENT_RECONSTRUCTION**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-l_LAX_FRIEDRICHS)=
    **l_LAX_FRIEDRICHS**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-l_adpo_flowstrength)=
    **l_adpo_flowstrength**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-l_with_horz_tracer_advection)=
    **l_with_horz_tracer_advection**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-l_with_horz_tracer_diffusion)=
    **l_with_horz_tracer_diffusion**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-l_with_vert_tracer_advection)=
    **l_with_vert_tracer_advection**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-l_with_vert_tracer_diffusion)=
    **l_with_vert_tracer_diffusion**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-no_tracer)=
    **no_tracer**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-threshold_max_S)=
    **threshold_max_S**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-threshold_max_T)=
    **threshold_max_T**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-threshold_min_S)=
    **threshold_min_S**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-threshold_min_T)=
    **threshold_min_T**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-tracer_HorizontalAdvection_type)=
    **tracer_HorizontalAdvection_type**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-tracer_update_mode)=
    **tracer_update_mode**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-use_draftave_for_transport_h)=
    **use_draftave_for_transport_h**
  -
  -
  -
  -
  -

* - (ocean_tracer_transport_nml-use_tracer_x_height)=
    **use_tracer_x_height**
  -
  -
  -
  -
  -
```


Defined and used in: {{ '[src/ocean/config/mo_ocean_nml.f90]({}/src/ocean/config/mo_ocean_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_ocean_vertical_diffusion_nml)=
## ocean_vertical_diffusion_nml

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ocean_vertical_diffusion_nml-KappaH_min)=
    **KappaH_min**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-KappaM_max)=
    **KappaM_max**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-KappaM_min)=
    **KappaM_min**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-PPscheme_type)=
    **PPscheme_type**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-RichardsonDiffusion_threshold)=
    **RichardsonDiffusion_threshold**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-Salinity_VerticalDiffusion_background)=
    **Salinity_VerticalDiffusion_background**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-Temperature_VerticalDiffusion_background)=
    **Temperature_VerticalDiffusion_background**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-VerticalViscosity_TimeWeight)=
    **VerticalViscosity_TimeWeight**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-WindMixingDecayDepth)=
    **WindMixingDecayDepth**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-alpha_tke)=
    **alpha_tke**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-bottom_drag_coeff)=
    **bottom_drag_coeff**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-c_eps)=
    **c_eps**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-c_k)=
    **c_k**
  - REAL
  - 1
  -
  -
  -

* - (ocean_vertical_diffusion_nml-cd)=
    **cd**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-clc)=
    **clc**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-convection_InstabilityThreshold)=
    **convection_InstabilityThreshold**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-fpath_iwe_botforc)=
    **fpath_iwe_botforc**
  -
  - 'idemix_bottom_forcing.nc'
  -
  -
  -

* - (ocean_vertical_diffusion_nml-fpath_iwe_surforc)=
    **fpath_iwe_surforc**
  -
  - 'idemix_surface_forcing.nc'
  -
  -
  -

* - (ocean_vertical_diffusion_nml-gamma)=
    **gamma**
  - REAL
  - 3.e-4_wp
  -
  -
  -

* - (ocean_vertical_diffusion_nml-jstar)=
    **jstar**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-l_idemix_osborn_cox_kv)=
    **l_idemix_osborn_cox_kv**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-l_lc)=
    **l_lc**
  - LOGICAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-l_use_idemix_forcing)=
    **l_use_idemix_forcing**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-lambda_wind)=
    **lambda_wind**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-mu0)=
    **mu0**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-mxl_min)=
    **mxl_min**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-n_hor_iwe_prop_iter)=
    **n_hor_iwe_prop_iter**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-name_iwe_botforc)=
    **name_iwe_botforc**
  -
  - 'wave_dissipation'
  -
  -
  -

* - (ocean_vertical_diffusion_nml-name_iwe_surforc)=
    **name_iwe_surforc**
  -
  - 'niw_forc'
  -
  -
  -

* - (ocean_vertical_diffusion_nml-only_tke)=
    **only_tke**
  - LOGICAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-tau_h)=
    **tau_h**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-tau_v)=
    **tau_v**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-tke_min)=
    **tke_min**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-tke_mxl_choice)=
    **tke_mxl_choice**
  - INTEGER
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-tke_surf_min)=
    **tke_surf_min**
  - REAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-tracer_RichardsonCoeff)=
    **tracer_RichardsonCoeff**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-tracer_TopWindMixing)=
    **tracer_TopWindMixing**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-tracer_convection_MixingCoefficient)=
    **tracer_convection_MixingCoefficient**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-use_Kappa_min)=
    **use_Kappa_min**
  - LOGICAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-use_lbound_dirichlet)=
    **use_lbound_dirichlet**
  - LOGICAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-use_reduced_mixing_under_ice)=
    **use_reduced_mixing_under_ice**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-use_ubound_dirichlet)=
    **use_ubound_dirichlet**
  - LOGICAL
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-use_wind_mixing)=
    **use_wind_mixing**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-velocity_RichardsonCoeff)=
    **velocity_RichardsonCoeff**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-velocity_TopWindMixing)=
    **velocity_TopWindMixing**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-velocity_VerticalDiffusion_background)=
    **velocity_VerticalDiffusion_background**
  -
  -
  -
  -
  -

* - (ocean_vertical_diffusion_nml-vert_mix_type)=
    **vert_mix_type**
  -
  -
  -
  -
  -
```


Defined and used in: {{ '[src/ocean/config/mo_ocean_nml.f90]({}/src/ocean/config/mo_ocean_nml.f90)'.format(base_url) }}


(ref_buildrun_nml_sea_ice_nml)=
## sea_ice_nml

Relevant if [iforcing](run_nml-iforcing)=2 (ECHAM)

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (sea_ice_nml-i_ice_therm)=
    **i_ice_therm**
  - I
  - 2
  -
  - Switch for thermodynamic model:
    - 1: Zero-layer model
    - 2: Two layer Winton (2000) model
    - 3: Zero-layer model with analytical forcing (for diagnostics)
    - 4: Zero-layer model for atmosphere-only runs (for diagnostics)
  - In an ocean run i_sea_ice must be >=1. In an atmospheric run the ice surface type must be defined.

* - (sea_ice_nml-i_ice_dyn)=
    **i_ice_dyn**
  - I
  - 0
  -
  - Switch for sea-ice dynamics:
    - 0: No dynamics
    - 1: FEM dynamics (from AWI)
  -

* - (sea_ice_nml-i_ice_albedo)=
    **i_ice_albedo**
  - I
  - 1
  -
  - Switch for albedo model. Only one is implemented so far.
  -

* - (sea_ice_nml-i_Qio_type)=
    **i_Qio_type**
  - I
  - 2
  -
  - Switch for ice-ocean heat-flux calculation method:
    - 1: Proportional to ocean cell thickness (like MPI-OM)
    - 2: Proportional to speed difference between ice and ocean
  - Defaults to 1 when i_ice_dyn=0 and 2 otherwise.

* - (sea_ice_nml-kice)=
    **kice**
  - I
  - 1
  -
  - Number of ice classes (must be one for now)
  -

* - (sea_ice_nml-hnull)=
    **hnull**
  - R
  - 0.5
  - m
  - Hibler's {math}`h_0` parameter for new-ice growth.
  -

* - (sea_ice_nml-hmin)=
    **hmin**
  - R
  - 0.05
  - m
  - Minimum sea-ice thickness allowed.
  -

* - (sea_ice_nml-ramp_wind)=
    **ramp_wind**
  - R
  - 10
  - days
  - Number of days it takes the wind to reach correct strength. Only used at the start of an OMIP/NCEP simulation (not after restart).
  -
```




# Ocean waves specific namelist parameters

The following Namelists become active, if ICON is configured with `--enable_waves`,
and if model_type=98 is selected.




(ref_buildrun_nml_energy_propagation_nml)=
## energy_propagation_nml

Used if ltransport=.TRUE.

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (energy_propagation_nml-itype_limit)=
    **itype_limit**
  - I
  - 0
  -
  - Type of limiter for wave energy transport:
    - 0: no limiter
    - 3: monotonic flux limiter (FCT)
    - 4: positive definite flux limiter
  -

* - (energy_propagation_nml-beta_fct)=
    **beta_fct**
  - R
  - 1.005
  -
  - global boost factor for range of permissible values {math}`\left[q_{max},q_{min}\right]` in the monotonic flux limiter. A value larger than 1 allows for (small) over and undershoots, while a value of 1 gives strict monotonicity (at the price of increased diffusivity).
  - itype_limit = 3

* - (energy_propagation_nml-igrad_c_miura)=
    **igrad_c_miura**
  - I
  - 1
  -
  - Method for gradient reconstruction at cell center for 2nd order miura scheme
    - 1: Least-squares (linear, non-consv)
    - 2: Green-Gauss
    - 3: based on shape function derivatives for a three-node triangular element (Fish. J and T. Belytschko, 2007)
  -

* - (energy_propagation_nml-lgrid_refr)=
    **lgrid_refr**
  - L
  - .TRUE.
  -
  - .TRUE.: calculate grid refraction
  -
```



Defined and used in: {{ '[src/waves/config/mo_energy_propagation_nml.f90]({}/src/waves/config/mo_energy_propagation_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_initwave_nml)=
## initwave_nml


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (initwave_nml-dt_shift)=
    **dt_shift**
  - R
  - 0.0
  - s
  - time interval by which the actual model start date (tc_start_date) is shifted backwards in time.
  -

* - (initwave_nml-init_mode)=
    **init_mode**
  - I
  - 2
  -
  - 1 (MODE_ANA): read wave energy spectrum from analysis file
    2 (MODE_COLD): initialize wave energy spectrum by analytic parameterization (currently JONSWAP)
  -

* - (initwave_nml-init_spectrum_from_file)=
    **init_spectrum_from_file**
  - C
  - ''
  -
  - path and name of initial condition file
  - init_mode=1

* - (initwave_nml-lskip_inv_post_op)=
    **lskip_inv_post_op**
  - L
  - .FALSE.
  -
  - The initial condition file contains {math}`S(\omega)` (.FALSE.) or {math}`S(f)` (.TRUE.)
  - init_mode=1
```



Defined and used in: {{ '[src/waves/config/mo_initwave_nml.f90]({}/src/waves/config/mo_initwave_nml.f90)'.format(base_url) }}



(ref_buildrun_nml_wave_nml)=
## wave_nml


```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (wave_nml-ndirs)=
    **ndirs**
  - I
  - 24
  -
  - number of direction of wave spectrum
  -

* - (wave_nml-nfreqs)=
    **nfreqs**
  - I
  - 25
  -
  - number of frequencies of wave spectrum
  -

* - (wave_nml-fr1)=
    **fr1**
  - R
  - 0.04177248
  - Hz
  - first frequency of wave spectrum
  -

* - (wave_nml-co)=
    **co**
  - R
  - 1.1
  -
  - frequency ratio
  -

* - (wave_nml-iref)=
    **iref**
  - I
  - 1
  -
  - frequency bin number of reference frequency
  -

* - (wave_nml-tlength)=
    **Tlength**
  - I
  - 10800
  - s
  - Length of time series for max. wave height calculation
  -

* - (wave_nml-alpha)=
    **alpha**
  - R
  - 0.018
  -
  - Phillips parameter
  -

* - (wave_nml-fm)=
    **fm**
  - R
  - 0.2
  - Hz
  - peak frequency and/or maximum frequency
  -

* - (wave_nml-gamma_wave)=
    **gamma_wave**
  - R
  - 3.0
  -
  - overshoot factor
  -

* - (wave_nml-sigma_a)=
    **sigma_a**
  - R
  - 0.07
  -
  - left peak width of wave spectrum
  -

* - (wave_nml-sigma_b)=
    **sigma_b**
  - R
  - 0.09
  -
  - right peak width of wave spectrum
  -

* - (wave_nml-fetch)=
    **fetch**
  - R
  - 300000.0
  - m
  - fetch
  -

* - (wave_nml-fetch_min_energy)=
    **fetch_min_energy**
  - R
  - 25000.0
  - m
  - fetch used for calculation of minimum allowed energy level
  -

* - (wave_nml-roair)=
    **roair**
  - R
  - 1.225
  - kg/m3
  - air density
  -

* - (wave_nml-rnuair)=
    **rnuair**
  - R
  - 1.5e-5
  - m2/s
  - kinematic air viscosity
  -

* - (wave_nml-rnuairm)=
    **rnuairm**
  - R
  - 0.11*rnuair
  - m2/s
  - kinematic air viscosity for momentum transfer
  -

* - (wave_nml-rowater)=
    **rowater**
  - R
  - 1000.0
  - kg/m3
  - water density
  -

* - (wave_nml-xeps)=
    **xeps**
  - R
  - roair/rowater
  -
  - air water density ratio
  -

* - (wave_nml-xinveps)=
    **xinveps**
  - R
  - 1./xeps
  -
  - inverse air water density ratio
  -

* - (wave_nml-betamax)=
    **betamax**
  - R
  - 1.20
  -
  - parameter for wind input (ECMWF cy45r1)
  -

* - (wave_nml-zalp)=
    **zalp**
  - R
  - 0.0080
  -
  - shifts growth curve (ECMWF cy45r1)
  -

* - (wave_nml-jtot_tauhf)=
    **jtot_tauhf**
  - I
  - 19
  -
  - dimension of high freuency wave stress (**wtauhf**)
  - must be odd

* - (wave_nml-alpha_ch)=
    **alpha_ch**
  - R
  - 0.0075
  -
  - minimum Charnock constant (ECMWF cy45r1)
  -

* - (wave_nml-oce_vct_filename)=
    **oce_vct_filename**
  - C
  -
  -
  - Name of ocean vertical coordinates file. This file must contain the number of ocean vertical interfaces on the first line **oce_nifc** (type integer), followed by lines with the interface index (type integer) and interface depth in meters **oce_ifc** (type real), separated by a space.
    For example:
    `21`
    `1 0.0`
    `2 2.0`
    `3 4.0`
    ...
    `21 200.0`
  - the interface index must be in the range from 1 to **oce_nifc**

* - (wave_nml-depth)=
    **depth**
  - R
  - 0.0
  - m
  - ocean depth if not 0, then constant depth
  -

* - (wave_nml-depth_min)=
    **depth_min**
  - R
  - 0.2
  - m
  - allowed minimum of model depth
  -

* - (wave_nml-depth_max)=
    **depth_max**
  - R
  - 999.0
  - m
  - allowed maximum of model depth
  -

* - (wave_nml-niter_smooth)=
    **niter_smooth**
  - I
  - 1
  -
  - number of smoothing iterations for wave bathymetry
  -

* - (wave_nml-nsubs_refrac)=
    **nsubs_refrac**
  - I
  - 1
  -
  - number of substeps in wave refraction calculation (recommendation: use the same ratio as between atmosphere and wave model)
  -

* - (wave_nml-xkappa)=
    **xkappa**
  - R
  - 0.40
  -
  - von Karman constant
  -

* - (wave_nml-xnlev)=
    **xnlev**
  - R
  - 10.0
  - m
  - windspeed reference level
  -

* - (wave_nml-linput_sf1)=
    **linput_sf1**
  - L
  - .TRUE.
  -
  - .TRUE.: calculate wind input source function term
  -

* - (wave_nml-linput_sf2)=
    **linput_sf2**
  - L
  - .TRUE.
  -
  - .TRUE.: update wind input source function term
  -

* - (wave_nml-ldissip_sf)=
    **ldissip_sf**
  - L
  - .TRUE.
  -
  - .TRUE.: calculate dissipation source function term
  -

* - (wave_nml-lwave_brk_sf)=
    **lwave_brk_sf**
  - L
  - .TRUE.
  -
  - .TRUE.: calculate depth-induced wave breaking dissipation source function term
  -

* - (wave_nml-lnon_linear_sf)=
    **lnon_linear_sf**
  - L
  - .TRUE.
  -
  - .TRUE.: calculate non linear source function term
  -

* - (wave_nml-lbottom_fric_sf)=
    **lbottom_fric_sf**
  - L
  - .TRUE.
  -
  - .TRUE.: calculate bottom friction source function term
  -

* - (wave_nml-lwave_stress1)=
    **lwave_stress1**
  - L
  - .TRUE.
  -
  - .TRUE.: calculate wave stress
  -

* - (wave_nml-lwave_stress2)=
    **lwave_stress2**
  - L
  - .TRUE.
  -
  - .TRUE.: update wave stress
  -

* - (wave_nml-impl_fac)=
    **impl_fac**
  - R
  - 1.0
  -
  - Implicitness factor for time integration scheme of total source function
    Range of permissible values: {math}`\left[0.5,\dots,1\right]`
    0.5: second order Crank-Nicholson scheme
    1.0: first order Euler backward scheme
  -

* - (wave_nml-forc_file_prefix)=
    **forc_file_prefix**
  - C
  -
  -
  - Common prefix of forcing files. If not empty, the names of forcing files will be consctructed as:
    - forc_file_prefix+_wind.nc - for 10m wind
    - forc_file_prefix+_ice.nc - for sea ice concentration
    - forc_file_prefix+_slh.nc - for sea level height
    - forc_file_prefix+_osc.nc - for ocean surface currents
    Data for all time steps in the current simulation should be prepared in a single file. Variables should be named `u_10m`, `v_10m`, `fr_seaice`, `uosc`, `vosc`
  - - For forc_file_prefix+_wind.nc, coupled_mode= .FALSE.
    - For forc_file_prefix+_ice.nc, coupled_mode= .FALSE.

* - (wave_nml-peak_lat)=
    **peak_lat**
  - R
  - -60.0
  - degree
  - latitude of wind peak value
  - ltestcase=.TRUE.

* - (wave_nml-peak_lon)=
    **peak_lon**
  - R
  - -140.0
  - degree
  - longitude of wind peak value
  - ltestcase=.TRUE.

* - (wave_nml-peak_wsp10)=
    **peak_wsp10**
  - R
  - 25.0
  - m/s
  - peak value of 10m wind speed at (peak_lat, peak_lon)
  - ltestcase=.TRUE.

* - (wave_nml-dir_wsp10)=
    **dir_wsp10**
  - R
  - 45.0
  - degree
  - wind direction measured clockwise from true north
  - ltestcase=.TRUE.

* - (wave_nml-r_wsp10)=
    **r_wsp10**
  - R
  - 1000000.0
  - m
  - physical radius of wind field
  - ltestcase=.TRUE.


* - (wave_nml-stokes_method)=
    **stokes_method**
  - I
  - 1
  -
  - 1 - calculation of Stokes profile from the full spectrum
    2 - calculation based on Breivik (2016)
  -

* - (wave_nml-stokes_depth)=
    **stokes_depth**
  - R
  - 50.0
  - m
  - maximum Stokes layer depth, if set to zero, the depth of the last interface from **oce_vct_filename** will be used
  -
```



Defined and used in: {{ '[src/waves/config/mo_wave_nml.f90]({}/src/waves/config/mo_wave_nml.f90)'.format(base_url) }}



# Namelist parameters for testcases

(NAMELIST_ICON) The ICON model code includes several experiments, so-called test cases,
for the 2 and 3-dimensional atmosphere. Depending on the specified experiment,
initial conditions and boundary conditions are computed internally.




(ref_buildrun_nml_nh_testcase_nml)=
## nh_testcase_nml

Scope: ltestcase=.TRUE.

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (nh_testcase_nml-nh_test_name)=
    **nh_test_name**
  - C
  - jabw
  -
  - testcase selection
    **zero**: no orography
    **bell**: bell shaped mountain at 0E,0N
    **schaer**: hilly mountain at 0E,0N
    **jabw**: Initializes the full Jablonowski Williamson test
    **jabw_s**: Initializes the Jablonowski Williamson steady state test case.
    **jabw_m**: Initializes the Jablonowski Williamson test case with a mountain instead of the wind perturbation (specify mount_height).
    **mrw_nh**: Initializes the full Mountain-induced Rossby wave test case.
    **mrw2_nh**: Initializes the modified mountain-induced Rossby wave test case.
    **mwbr_const**: Initializes the mountain wave with two layers test case. The lower layer is isothermal and the upper layer has constant brunt vaisala frequency. The interface has constant pressure.
    **PA**: Initializes the pure advection test case.
    **HS_nh**: Initializes the Held-Suarez test case. At the moment with an isothermal atmosphere at rest (T=300K, ps=1000hPa, u=v=0, topography=0.0).
    **HS_jw**: Initializes the Held-Suarez test case with Jablonowski Williamson initial conditions and zero topography.
    **APE_nwp, APE_aes, APE_nh, APEc_nh**: Initializes the APE experiments. With the jabw test case, including moisture.
    **wk82**: Initializes the Weisman Klemp test case
    **g_lim_area**: Initializes a series of general limited area test cases: itype_atmos_ana determines the atmospheric profile, itype_anaprof_uv determines the wind profile and itype_topo_ana determines the topography
    **dcmip_bw_11**: Initializes (moist) baroclinic instability/wave (**DCMIP2016**)
    **dcmip_pa_12**: Initializes Hadley-like meridional circulation pure advection test case.
    **dcmip_rest_200**: atmosphere at rest test (Schaer-type mountain)
    **dcmip_mw_2x**: nonhydrostatic mountain waves triggered by Schaer-type mountain
    **dcmip_gw_31**: nonhydrostatic gravity waves triggered by a localized perturbation (nonlinear)
    **dcmip_gw_32**: nonhydrostatic gravity waves triggered by a localized perturbation (linear)
    **dcmip_tc_52**: tropical cyclone test case with with full physics in Aqua-planet mode
    **CBL**: convective boundary layer simulations for LES package on torus (doubly periodic) grid
    **bb13**: linear gravity- and sound-wave expansion in a channel (Baldauf, Brdar (2013) QJRMS)
    **SCM**: Single Column Mode
  - **schaer**: is_plane_torus=.TRUE.
    **wk82**: l_limited_area =.TRUE.
    **dcmip_rest_200**: [lcoriolis](dynamics_nml-lcoriolis) = .FALSE.
    **dcmip_mw_2x**: [lcoriolis](dynamics_nml-lcoriolis) = .FALSE.
    **dcmip_gw_32**: l_limited_area =.TRUE. and [lcoriolis](dynamics_nml-lcoriolis) = .FALSE.
    **dcmip_tc_52**: [lcoriolis](dynamics_nml-lcoriolis) = .TRUE.
    **CBL**: is_plane_torus= .TRUE.
    **bb13**: is_plane_torus= .TRUE.
    **SCM**: is_plane_torus= .TRUE.

* - (nh_testcase_nml-is_toy_chem)=
    **is_toy_chem**
  - L
  - .FALSE.
  -
  - Terminator toy chemistry activated when .TRUE.
  -

* - (nh_testcase_nml-tracer_inidist_list)=
    **tracer_inidist_list**
  - I(:)
  - 1
  -
  - For a subset of testcases pre-defined initial tracer distributions are available. This namelist parameter specifies the initial distribution for each tracer. In the following the testcases and the pre-defined numbers are given:
    `PA': 4,5,6,7,8
    'JABW':1,2,3,4
    'DF': 5,6,7,8,9
    For more details on the initial distributions, please have a look into the code.
  - nh_test_name='PA','JABW','DF'

* - (nh_testcase_nml-dcmip_bw_)=
    **dcmip_bw%...**
  -
  -
  -
  - DCMIP2016 baroclinic wave test
  - 'dcmip_bw_11'

* - (nh_testcase_nml-dcmip_bw_-deep)=
    **...deep**
  - I
  - 0
  -
  - deep atmosphere
    - 1: yes or
    - 0: no
  -

* - (nh_testcase_nml-dcmip_bw_-moist)=
    **...moist**
  - I
  - 0
  -
  - include moisture, i.e. {math}`qv\neq 0`
    - 1: yes or
    - 0: no
  -

* - (nh_testcase_nml-dcmip_bw_-pertt)=
    **...pertt**
  - I
  - 0
  -
  - type of initial perturbation
    - 0: exponential
    - 1: stream function
  -

* - (nh_testcase_nml-toy_chem_)=
    **toy_chem%...**
  -
  -
  -
  - terminator toy chemistry
  - is_toy_chem=.TRUE.

* - (nh_testcase_nml-toy_chem_-dt_chem)=
    **...dt_chem**
  - R
  - 300
  - s
  - chemistry tendency update interval
  -

* - (nh_testcase_nml-toy_chem_-dt_cpl)=
    **...dt_cpl**
  - R
  - 300
  - s
  - chemistry-transport coupling interval
  -

* - (nh_testcase_nml-toy_chem_-id_cl)=
    **...id_cl**
  - I
  - 1
  -
  - Tracer container slice index for species CL
  -

* - (nh_testcase_nml-toy_chem_-id_cl2)=
    **...id_cl2**
  - I
  - 2
  -
  - Tracer container slice index for species CL2
  -

* - (nh_testcase_nml-jw_up)=
    **jw_up**
  - R
  - 1.0
  - m/s
  - amplitude of the u-perturbation in jabw test case
  - nh_test_name='jabw'

* - (nh_testcase_nml-jw_u0)=
    **jw_u0**
  - R
  - 35.0
  - m/s
  - maximum zonal wind in jabw test case
  - nh_test_name='jabw'

* - (nh_testcase_nml-jw_temp0)=
    **jw_temp0**
  - R
  - 288.0
  - K
  - horizontal-mean temperature at surface in jabw test case
  - nh_test_name='jabw'

* - (nh_testcase_nml-u0_mrw)=
    **u0_mrw**
  - R
  - 20.0
  - m/s
  - wind speed for mrw(2) and mwbr_const cases
  - nh_test_name= 'mrw(2)_nh' and 'mwbr_const'

* - (nh_testcase_nml-mount_height_mrw)=
    **mount_height_mrw**
  - R
  - 2000.0
  - m
  - maximum mount height in mrw(2) and mwbr_const
  - nh_test_name= 'mrw(2)_nh' and 'mwbr_const'

* - (nh_testcase_nml-mount_half_width)=
    **mount_half_width**
  - R
  - 1500000.0
  - m
  - half width of mountain in mrw(2), mwbr_const and bell
  - nh_test_name= 'mrw(2)_nh', 'mwbr_const' and 'bell'

* - (nh_testcase_nml-mount_width)=
    **mount_width**
  - R
  - 1000.0
  - m
  - width of mountain
  -

* - (nh_testcase_nml-mount_width_2)=
    **mount_width_2**
  - R
  - 100.0
  - m
  - a 2nd width scale of mountain
  - nh_test_name='schaer'

* - (nh_testcase_nml-mount_lonctr_mrw_deg)=
    **mount_lonctr_mrw_deg**
  - R
  - 90.0
  - deg
  - lon of mountain center in mrw(2) and mwbr_const
  - nh_test_name= 'mrw(2)_nh' and 'mwbr_const'

* - (nh_testcase_nml-mount_latctr_mrw_deg)=
    **mount_latctr_mrw_deg**
  - R
  - 30.0
  - deg
  - lat of mountain center in mrw(2) and mwbr_const
  - nh_test_name= 'mrw(2)_nh' and 'mwbr_const'

* - (nh_testcase_nml-temp_i_mwbr_const)=
    **temp_i_mwbr_const**
  - R
  - 288.0
  - K
  - temp at isothermal lower layer for mwbr_const case
  - nh_test_name= 'mwbr_const'

* - (nh_testcase_nml-p_int_mwbr_const)=
    **p_int_mwbr_const**
  - R
  - 70000.0
  - Pa
  - pres at the interface of the two layers for mwbr_const case
  - nh_test_name= 'mwbr_const'

* - (nh_testcase_nml-bruntvais_u_mwbr_const)=
    **bruntvais_u_mwbr_const**
  - R
  - 0.025
  - 1/s
  - constant brunt vaissala frequency at upper layer for mwbr_const case
  - nh_test_name= 'mwbr_const'

* - (nh_testcase_nml-mount_height)=
    **mount_height**
  - R
  - 100.0
  - m
  - peak height of mountain
  - nh_test_name=  'bell'

* - (nh_testcase_nml-layer_thickness)=
    **layer_thickness**
  - R
  - -999.0
  - m
  - thickness of vertical layers
  - If layer_thickness < 0, the vertical level distribution is read in from externally given HYB_PARAMS_XX.

* - (nh_testcase_nml-n_flat_level)=
    **n_flat_level**
  - I
  - 2
  -
  - level number for which the layer is still flat and not terrain-following
  - layer_thickness > 0

* - (nh_testcase_nml-nh_u0)=
    **nh_u0**
  - R
  - 0.0
  - m/s
  - initial constant zonal wind speed
  - nh_test_name = 'bell'

* - (nh_testcase_nml-nh_t0)=
    **nh_t0**
  - R
  - 300.0
  - K
  - initial temperature at lowest level
  - nh_test_name = 'bell'

* - (nh_testcase_nml-nh_brunt_vais)=
    **nh_brunt_vais**
  - R
  - 0.01
  - 1/s
  - initial Brunt-Vaisala frequency
  - nh_test_name = 'bell'

* - (nh_testcase_nml-torus_domain_length)=
    **torus_domain_length**
  - R
  - 100000.0
  - m
  - length of slice domain
  - nh_test_name = 'bell', lplane=.TRUE.

* - (nh_testcase_nml-rotate_axis_deg)=
    **rotate_axis_deg**
  - R
  - 0.0
  - deg
  - Earth's rotation axis pitch angle
  - nh_test_name= 'PA'

* - (nh_testcase_nml-lhs_nh_vn_ptb)=
    **lhs_nh_vn_ptb**
  - L
  - .TRUE.
  -
  - Add random noise to the initial wind field in the Held-Suarez test.
  - nh_test_name= 'HS_nh'

* - (nh_testcase_nml-lhs_fric_heat)=
    **lhs_fric_heat**
  - L
  - .FALSE.
  -
  - add frictional heating from Rayleigh friction in the Held-Suarez test.
  - nh_test_name= 'HS_nh'

* - (nh_testcase_nml-hs_nh_vn_ptb_scale)=
    **hs_nh_vn_ptb_scale**
  - R
  - 1.0
  - m/s
  - Magnitude of the random noise added to the initial wind field in the Held-Suarez test.
  - nh_test_name= 'HS_nh'

* - (nh_testcase_nml-rh_at_1000hpa)=
    **rh_at_1000hpa**
  - R
  - 0.7
  - 1
  - relative humidity at 1000 hPa
  - nh_test_name= 'jabw', nh_test_name= 'mrw'

* - (nh_testcase_nml-qv_max)=
    **qv_max**
  - R
  - 20.e-3
  - kg/kg
  - specific humidity in the tropics
  - nh_test_name= 'jabw', nh_test_name= 'mrw'

* - (nh_testcase_nml-ape_sst_case)=
    **ape_sst_case**
  - C
  - 'sst1'
  -
  - SST distribution selection
    'sst1': Control experiment
    'sst2': Peaked experiment
    'sst3': Flat experiment
    'sst4': Control-5N experiment
    'sst_qobs': Qobs SST distribution exp.
    'sst_const': constant SST
  - nh_test_name='APE_nwp', 'APE_aes'

* - (nh_testcase_nml-ape_sst_val)=
    **ape_sst_val**
  - R
  - 29.0
  - degC
  - aqua planet SST  for ape_sst_case='sst_const'
  - nh_test_name= 'APE_nwp', 'APE_aes'

* - (nh_testcase_nml-linit_tracer_fv)=
    **linit_tracer_fv**
  - L
  - .TRUE.
  -
  - Finite volume initialization for tracer fields
  - pure advection tests, only

* - (nh_testcase_nml-lcoupled_rho)=
    **lcoupled_rho**
  - L
  - .FALSE.
  -
  - Integrate density equation 'offline'
  - pure advection tests, only

* - (nh_testcase_nml-qv_max_wk)=
    **qv_max_wk**
  - R
  - 0.014
  - Kg/kg
  - maximum specific humidity near the surface, range  0.012 - 0.016
    used to vary the buoyancy
  - nh_test_name='wk82'

* - (nh_testcase_nml-u_infty_wk)=
    **u_infty_wk**
  - R
  - 20.0
  - m/s
  - zonal wind at infinity height, range 0. - 45.0
    used to vary the wind shear
  - nh_test_name='wk82', 'bb13'

* - (nh_testcase_nml-bub_amp)=
    **bub_amp**
  - R
  - 2.0
  - K
  - maximum amplitud of the thermal perturbation
  - nh_test_name='wk82'

* - (nh_testcase_nml-bubctr_lat)=
    **bubctr_lat**
  - R
  - 0.0
  - deg
  - latitude of the center of the thermal perturbation
  - nh_test_name='wk82'

* - (nh_testcase_nml-bubctr_lon)=
    **bubctr_lon**
  - R
  - 90.0
  - deg
  - longitude of the center of the thermal perturbation
  - nh_test_name='wk82'

* - (nh_testcase_nml-bubctr_x)=
    **bubctr_x**
  - R
  - 0.0
  - m
  - x-position of the center of the thermal perturbation
  - is_plane_grid=.TRUE.

* - (nh_testcase_nml-bubctr_y)=
    **bubctr_y**
  - R
  - 0.0
  - m
  - y-position of the center of the thermal perturbation
  - is_plane_grid=.TRUE.

* - (nh_testcase_nml-bubctr_z)=
    **bubctr_z**
  - R
  - 1400.0
  - m
  - height of the center of the thermal perturbation
  - nh_test_name='wk82'

* - (nh_testcase_nml-bub_hor_width)=
    **bub_hor_width**
  - R
  - 10000.0
  - m
  - horizontal radius of the thermal perturbation
  - nh_test_name='wk82'

* - (nh_testcase_nml-bub_ver_width)=
    **bub_ver_width**
  - R
  - 1400.0
  - m
  - vertical radius of the thermal perturbation
  - nh_test_name='wk82'

* - (nh_testcase_nml-itype_atmo_ana)=
    **itype_atmo_ana**
  - I
  - 1
  -
  - kind of atmospheric profile:
    1 piecewise N constant layers
    2 piecewise polytropic layers
  - nh_test_name='g_lim_area'

* - (nh_testcase_nml-itype_anaprof_uv)=
    **itype_anaprof_uv**
  - I
  - 1
  -
  - kind of wind profile:
    - 1: piecewise linear wind layers
    - 2: constant wind
    - (3: constant meridional wind, deprecated)
  - nh_test_name='g_lim_area'

* - (nh_testcase_nml-itype_topo_ana)=
    **itype_topo_ana**
  - I
  - 1
  -
  - kind of orography:
    1 schaer test case mountain
    2 gaussian_2d mountain
    3 gaussian_3d mountain
    any other no orography
  - nh_test_name='g_lim_area'

* - (nh_testcase_nml-nlayers_nconst)=
    **nlayers_nconst**
  - I
  - 1
  -
  - Number of the desired layers with a constant Brunt-Vaisala-frequency
  - nh_test_name='g_lim_area' and itype_atmo_ana=1

* - (nh_testcase_nml-p_base_nconst)=
    **p_base_nconst**
  - R
  - 100000.0
  - Pa
  - pressure at the base of the first N constant layer
  - nh_test_name='g_lim_area' and itype_atmo_ana=1

* - (nh_testcase_nml-theta0_base_nconst)=
    **theta0_base_nconst**
  - R
  - 288.0
  - K
  - potential temperature at the base of the first N constant layer
  - nh_test_name='g_lim_area' and itype_atmo_ana=1

* - (nh_testcase_nml-h_nconst)=
    **h_nconst**
  - R(nlayers _nconst)
  - 0., 1500., 12000.
  - m
  - height of the base of each of the N constant layers
  - nh_test_name='g_lim_area' and itype_atmo_ana=1

* - (nh_testcase_nml-N_nconst)=
    **N_nconst**
  - R(nlayers _nconst)
  - 0.01
  - 1/s
  - Brunt-Vaisala-frequency at each of the N constant layers
  - nh_test_name='g_lim_area' and itype_atmo_ana=1

* - (nh_testcase_nml-rh_nconst)=
    **rh_nconst**
  - R(nlayers _nconst)
  - 0.5
  - %
  - relative humidity at the base of each N constant layers
  - nh_test_name='g_lim_area' and itype_atmo_ana=1

* - (nh_testcase_nml-rhgr_nconst)=
    **rhgr_nconst**
  - R(nlayers _nconst)
  - 0.0
  - %
  - relative humidity gradient at each of the N constant layers
  - nh_test_name='g_lim_area' and itype_atmo_ana=1

* - (nh_testcase_nml-nlayers_poly)=
    **nlayers_poly**
  - I
  - 2
  -
  - Number of the desired layers with constant gradient temperature
  - nh_test_name='g_lim_area' and itype_atmo_ana=2

* - (nh_testcase_nml-p_base_poly)=
    **p_base_poly**
  - R
  - 100000.0
  - Pa
  - pressure at the base of the first polytropic layer
  - nh_test_name='g_lim_area' and itype_atmo_ana=2

* - (nh_testcase_nml-h_poly)=
    **h_poly**
  - R(nlayers _poly)
  - 0., 12000.
  - m
  - height of the base of each of the polytropic layers
  - nh_test_name='g_lim_area' and itype_atmo_ana=2

* - (nh_testcase_nml-t_poly)=
    **t_poly**
  - R(nlayers _poly)
  - 288., 213.
  - K
  - temperature at the base of each of the polytropic layers
  - nh_test_name='g_lim_area' and itype_atmo_ana=2

* - (nh_testcase_nml-rh_poly)=
    **rh_poly**
  - R(nlayers _poly)
  - 0.8, 0.2
  - %
  - relative humidity at the base of each of the polytropic layers
  - nh_test_name='g_lim_area' and itype_atmo_ana=2

* - (nh_testcase_nml-rhgr_poly)=
    **rhgr_poly**
  - R(nlayers _poly)
  - 5.e-5, 0.
  - %
  - relative humidity gradient at each of the polytropic layers
  - nh_test_name='g_lim_area' and itype_atmo_ana=2

* - (nh_testcase_nml-nlayers_linwind)=
    **nlayers_linwind**
  - I
  - 2
  -
  - Number of the desired layers with constant U gradient
  - nh_test_name='g_lim_area' and itype_anaprof_uv=1

* - (nh_testcase_nml-h_linwind)=
    **h_linwind**
  - R(nlayers _linwind)
  - 0., 2500.
  - m
  - height of the base of each of the linear wind layers
  - nh_test_name='g_lim_area' and itype_anaprof_uv=1

* - (nh_testcase_nml-u_linwind)=
    **u_linwind**
  - R(nlayers _linwind)
  - 5,  10.
  - m/s
  - zonal wind at the base of each of the linear wind layers
  - nh_test_name='g_lim_area' and itype_anaprof_uv=1

* - (nh_testcase_nml-ugr_linwind)=
    **ugr_linwind**
  - R(nlayers _linwind)
  - 0., 0.
  - 1/s
  - zonal wind gradient at each of the linear wind layers
  - nh_test_name='g_lim_area' and itype_anaprof_uv=1

* - (nh_testcase_nml-vel_const)=
    **vel_const**
  - R
  - 20.0
  - m/s
  - constant wind speed (itype_anaprof_uv=2)
  - nh_test_name='g_lim_area' and itype_anaprof_uv=2

* - (nh_testcase_nml-dir_wind)=
    **dir_wind**
  - R
  - 270.0
  - deg
  - wind direction
    - 0: N,
    - 90: E
  - nh_test_name='g_lim_area'

* - (nh_testcase_nml-mount_lonc_deg)=
    **mount_lonc_deg**
  - R
  - 90.0
  - deg
  - longitud of the center of the  mountain
  - nh_test_name='g_lim_area'

* - (nh_testcase_nml-mount_latc_deg)=
    **mount_latc_deg**
  - R
  - 0.0
  - deg
  - latitud of the center of the  mountain
  - nh_test_name='g_lim_area'

* - (nh_testcase_nml-schaer_h0)=
    **schaer_h0**
  - R
  - 250.0
  - m
  - h0 parameter for the schaer mountain
  - nh_test_name='g_lim_area' and itype_topo_ana=1

* - (nh_testcase_nml-schaer_a)=
    **schaer_a**
  - R
  - 5000.0
  - m
  - -a- parameter for the schaer mountain, also half width in the north and south side of the finite ridge to round the sharp edges
  - nh_test_name='g_lim_area' and itype_topo_ana=1,2

* - (nh_testcase_nml-schaer_lambda)=
    **schaer_lambda**
  - R
  - 4000.0
  - m
  - lambda parameter for the schaer mountain
  - nh_test_name='g_lim_area' and itype_topo_ana=1

* - (nh_testcase_nml-lshear_dcmip)=
    **lshear_dcmip**
  - L
  - FALSE
  -
  - run dcmip_mw_2x with/without vertical wind shear
    FALSE: dcmip_mw_21: non-sheared
    TRUE : dcmip_mw_22: sheared
  - nh_test_name='dcmip_mw_2x'

* - (nh_testcase_nml-halfwidth_2d)=
    **halfwidth_2d**
  - R
  - 10000.0
  - m
  - half length of the finite ridge in the north-south direction
  - nh_test_name='g_lim_area' and itype_topo_ana=1,2

* - (nh_testcase_nml-m_height)=
    **m_height**
  - R
  - 1000.0
  - m
  - height of the mountain
  - nh_test_name='g_lim_area' and itype_topo_ana=2,3

* - (nh_testcase_nml-m_width_x)=
    **m_width_x**
  - R
  - 5000.0
  - m
  - half width of the gaussian mountain in the east-west direction
    half width in the north-south direction in the rounding of the finite ridge (gaussian_2d)
  - nh_test_name='g_lim_area' and itype_topo_ana=2,3

* - (nh_testcase_nml-m_width_y)=
    **m_width_y**
  - R
  - 5000.0
  - m
  - half width of the gaussian mountain in the north-south direction
  - nh_test_name='g_lim_area' and itype_topo_ana=2,3

* - (nh_testcase_nml-gw_u0)=
    **gw_u0**
  - R
  - 0.0
  - m/s
  - maximum amplitude of the zonal wind
  - nh_test_name='dcmip_gw_3X'

* - (nh_testcase_nml-gw_clat)=
    **gw_clat**
  - R
  - 90.0
  - deg
  - Lat of perturbation center
  - nh_test_name='dcmip_gw_3X'

* - (nh_testcase_nml-gw_delta_temp)=
    **gw_delta_temp**
  - R
  - 0.01
  - K
  - maximum temperature perturbation
  - nh_test_name='dcmip_gw_32'

* - (nh_testcase_nml-u_cbl(2))=
    **u_cbl(2)**
  - R
  - 0:0
  - m/s and 1/s
  - to prescribe initial zonal velocity profile for convective boundary layer simulations where u_cbl(1) sets the constant and u_cbl(2) sets the vertical gradient
  - nh_test_name=CBL

* - (nh_testcase_nml-v_cbl(2))=
    **v_cbl(2)**
  - R
  - 0:0
  - m/s and 1/s
  - to prescribe initial meridional velocity profile for convective boundary layer simulations where v_cbl(1) sets the constant and v_cbl(2) sets the vertical gradient
  - nh_test_name=CBL

* - (nh_testcase_nml-th_cbl(2))=
    **th_cbl(2)**
  - R
  - 290:0.006
  - K and K/m
  - to prescribe initial potential temperature profile for convective boundary layer simulations where th_cbl(1) sets the constant and th_cbl(2) sets the gradient
  - nh_test_name=CBL
```



Defined and used in: {{ '[src/testcases/mo_nh_testcases.f90]({}/src/testcases/mo_nh_testcases.f90)'.format(base_url) }}


# External data



(ref_buildrun_nml_extpar_nml)=
## extpar_nml

Scope: itopo=1 in [run_nml](ref_buildrun_nml_run_nml)

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (extpar_nml-itopo)=
    **itopo**
  - I
  - 0
  -
  - 0: analytical topography/ext. data
    - 1: topography/ext. data read from file
  -

* - (extpar_nml-itype_vegetation_cycle)=
    **itype_vegetation_cycle**
  - I
  - 1
  -
  - 1: annual cycle of LAI solely based on NDVI climatology
    - 2: additional use of monthly T2M climatology to get more realistic values in extratropics (requires external parameter data containing this field)
  -

* - (extpar_nml-n_iter_smooth_topo)=
    **n_iter_smooth_topo**
  - I(n_dom)
  - 0
  -
  - iterations of topography smoother
  - itopo = 1

* - (extpar_nml-fac_smooth_topo)=
    **fac_smooth_topo**
  - R
  - 0.015625
  -
  - pre-factor of topography smoother
  - n_iter_smooth_topo > 0

* - (extpar_nml-hgtdiff_max_smooth_topo)=
    **hgtdiff_max_smooth_topo**
  - R(n_dom)
  - 0.0
  - m
  - RMS height difference to neighbor grid points at which the smoothing pre-factor fac_smooth_topo reaches its maximum value (linear proportionality for weaker slopes)
  - n_iter_smooth_topo > 0

* - (extpar_nml-heightdiff_threshold)=
    **heightdiff_threshold**
  - R(n_dom)
  - 3000.0
  - m
  - height difference between neighboring grid points above which additional local nabla2 diffusion is applied
  -

* - (extpar_nml-pp_sso)=
    **pp_sso**
  - I
  - 1
  -
  - 1: Postprocess SSO standard deviation and slope over glaciers based on the ratio between grid-scale and subgrid-scale slope: both quantities are reduced if the subgrid-scale slope calculated in extpar largely reflects the grid-scale slope.
    - 2: Optimized tuning for MERIT/REMA orography data: the reduction is also applied at non-glacier points in the Arctic, and the adjustment of the SSO standard deviation to orography smoothing is turned off.
  - n_iter_smooth_topo > 0

* - (extpar_nml-lrevert_sea_height)=
    **lrevert_sea_height**
  - L
  - .FALSE.
  -
  - If .TRUE., sea point heights will be reverted to original (raw data) heights after topography smoothing was applied.
  - n_iter_smooth_topo > 0

* - (extpar_nml-itype_lwemiss)=
    **itype_lwemiss**
  - I
  - 1
  -
  - Type of data used for longwave surface emissivity:
    - 0: No data; use constant fallback value instead
    - 1: Read and use emissivities derived in extpar from landuse classes
    - 2: Read and use monthly climatologies derived from satellite measurements
  - itopo = 1

* - (extpar_nml-extpar_filename)=
    **extpar_filename**
  - C
  -
  -
  - Filename of external parameter input file, default: `<path>extpar_<gridfile>`. May contain the keyword `<path>` which will be substituted by `model_base_dir`.
  -

* - (extpar_nml-read_nc_via_cdi)=
    **read_nc_via_cdi**
  - L
  - .FALSE.
  -
  - .TRUE.: read NetCDF input data via cdi library
    .FALSE.: read NetCDF input data using parallel NetCDF library
    Note: GRIB2 input data is always read via cdi library / GRIB API. For NetCDF input, this switch allows optimizing the input performance, but there is no general rule which option is faster.
  -

* - (extpar_nml-extpar_varnames_map_ file)=
    **extpar_varnames_map_ file**
  - C
  - ' '
  -
  - Filename of external parameter dictionary, This is a text file with two columns separated by whitespace, where left column: NetCDF name, right column: GRIB2 short name. It is required, if external parameter are read from a file in GRIB2 format.
  -
```



Defined and used in: {{ '[src/namelists/mo_extpar_nml.f90]({}/src/namelists/mo_extpar_nml.f90)'.format(base_url) }}



# Serialization

(ref_buildrun_nml_ser_nml)=
## ser_nml

Some developments must not change model results. Serialbox allows reading and writing data at any point in ICON into savepoints. These savepoints can be used to restore model variables to some reference or compare different model versions. The simplest application of Serialbox is using {{ '[src/serialization/mo_ser_debug.f90]({}/src/serialization/mo_ser_debug.f90)'.format(base_url) }} (or writing a similar routine fitting ones needs). Following this method will allow reading and writing manually specified fields in ICON. This can be very useful for small subroutines where input and output are clearly specified (i.e. do not involve derived types) and can thus easily be translated to Serialbox read/write statements. For larger components (basically everything hanging from {{ '[src/atm_dyn_iconam/mo_nh_stepping.f90]({}/src/atm_dyn_iconam/mo_nh_stepping.f90)'.format(base_url) }}, e.g. nwp_physics) the interface is specified by the in and out types. The actual fields that are read or written to in these subroutines are not specified. For this purpose, serialize_all has been implemented. It provides a wrapper for Serialbox read and write statements by looping through variable lists. This approach does not require managing lists of fields to read or write by Serialbox. At the level of {{ '[src/atm_dyn_iconam/mo_nh_stepping.f90]({}/src/atm_dyn_iconam/mo_nh_stepping.f90)'.format(base_url) }} and {{ '[src/atm_phy_nwp.mo_nh_interface_nwp.f90]({}/src/atm_phy_nwp.mo_nh_interface_nwp.f90)'.format(base_url) }} many components are wrapped by such serialize_all calls that allow testing these components. Each of these hard-coded calls to serialize_all has a name and for each name there is a namelist switch specifying the following triplet (e.g. 0,12,12):

 - If 0 do not use this savepoint, else use this savepoint at every time step
 - the relative threshold for errors (given as {math}`N` for {math}`N` in {math}`10^{-N}`)
 - the absolute threshold for errors (given as {math}`N` for {math}`N` in {math}`10^{-N}`)



```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (ser_nml-ser_initialization)=
    **ser_initialization**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for initial data (Checked after regular initialization at model start as well as after initialization of nested domains during model run)
  -

* - (ser_nml-ser_output_diag_dyn)=
    **ser_output_diag_dyn**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for output diagnostics of dynamics fields
  -

* - (ser_nml-ser_output_diag)=
    **ser_output_diag**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for output diagnostics
  -

* - (ser_nml-ser_output_opt)=
    **ser_output_opt**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for optional output
  -

* - (ser_nml-ser_latbc_data)=
    **ser_latbc_data**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine recv_latbc_data
  -

* - (ser_nml-ser_nesting_save_progvars)=
    **ser_nesting_save_progvars**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine save_progvars which is related to nesting
  -

* - (ser_nml-ser_dynamics)=
    **ser_dynamics**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine perform_dyn_substepping
  -

* - (ser_nml-ser_diffusion)=
    **ser_diffusion**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine diffusion
  -

* - (ser_nml-ser_nesting_compute_tendencies)=
    **ser_nesting_compute_tendencies**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine compute_tendencies (related to nesting)
  -

* - (ser_nml-ser_nesting_boundary_interpolation)=
    **ser_nesting_boundary_interpolation**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine boundary_interpolation (related to nesting)
  -

* - (ser_nml-ser_nesting_relax_feedback)=
    **ser_nesting_relax_feedback**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine relax_feedback (related to nesting)
  -

* - (ser_nml-ser_step_advection)=
    **ser_step_advection**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine step_advection
  -

* - (ser_nml-ser_physics)=
    **ser_physics**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine nwp_nh_interface
  -

* - (ser_nml-ser_physics_init)=
    **ser_physics_init**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine nwp_nh_interface during initialization
  -

* - (ser_nml-ser_lhn)=
    **ser_lhn**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine organize_lhn
  -

* - (ser_nml-ser_nudging)=
    **ser_nudging**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the nudging computations
  -

* - (ser_nml-ser_surface)=
    **ser_surface**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine nwp_surface
  -

* - (ser_nml-ser_microphysics)=
    **ser_microphysics**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine nwp_microphysics
  -

* - (ser_nml-ser_turbtrans)=
    **ser_turbtrans**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine nwp_turbtrans
  -

* - (ser_nml-ser_turbdiff)=
    **ser_turbdiff**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine nwp_turbdiff
  -

* - (ser_nml-ser_convection)=
    **ser_convection**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine nwp_convection
  -

* - (ser_nml-ser_cover)=
    **ser_cover**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine cover_koe
  -

* - (ser_nml-ser_radiation)=
    **ser_radiation**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine nwp_radiation
  -

* - (ser_nml-ser_radheat)=
    **ser_radheat**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the computations involving radiative heating
  -

* - (ser_nml-ser_gwdrag)=
    **ser_gwdrag**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Serialization switch for the subroutine nwp_gwdrag
  -

* - (ser_nml-ser_time_loop_end)=
    **ser_time_loop_end**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Check the state at the end of the time loop (does not read in data)
  -

* - (ser_nml-ser_reset_to_initial_state)=
    **ser_reset_to_initial_state**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Check the reset to initial state after the first phase of IAU
  -

* - (ser_nml-ser_all_debug)=
    **ser_all_debug**
  - I (3)
  - 0,12,12
  - -, {math}`10^{-N}`, {math}`10^{-N}`
  - Additional calls to serialize_all (for debugging purposes) can be controlled using this switch.
  -

* - (ser_nml-ser_nfail)=
    **ser_nfail**
  - R
  - 1.0
  - %
  - Fields that fail more elements than the percentage specified by ser_nfail will be reported.
  -

* - (ser_nml-ser_nreport)=
    **ser_nreport**
  - I
  - 10
  -
  - The detailed serialization report will include the ser_nreport elements with largest relative differences to the reference
  -

* - (ser_nml-ser_debug)=
    **ser_debug**
  - L
  - .FALSE.
  -
  - Activates the debug serialization defined in {{ '[src/serialization/mo_ser_debug.f90]({}/src/serialization/mo_ser_debug.f90)'.format(base_url) }}
  -
```



Defined and used in: {{ '[src/namelists/mo_ser_nml.f90]({}/src/namelists/mo_ser_nml.f90)'.format(base_url) }}


# External packages

## Community Interface (ComIn)

### plugin_list

This namelist is only usable (and needed) if the [Community Interface (ComIn)](ref_tools_comin) has been enabled during ICON's configure process.
Several `plugin_list(pg)` definitions can be specified, each holding the following namelist settings. The numbering (pg) of the plugins has to be contiguous without any gaps.

```{list-table}
:header-rows: 1

* - Parameter
  - Type
  - Default
  - Unit
  - Description
  - Scope

* - (plugin_list-plugin_library)=
    **plugin_library**
  - C
  - ""
  -
  - Path to the plugin library file. If omitted the primary constructor specified by `primary_constructor` is loaded from the `icon` executable (static linkage)
  - `configure --enable-comin`

* - (plugin_list-name)=
    **name**
  - C
  - ""
  -
  - Name of the plugin -- currently only used for messages.
  - `configure --enable-comin`

* - (plugin_list-primary_constructor)=
    **primary_constructor**
  - C
  - "comin_main"
  -
  - Name of the symbol in the plugin library that holds the primary constructor.
  - `configure --enable-comin`

* - (plugin_list-options)=
    **options**
  - C
  - ""
  -
  - The options string passed to the plugin. This namelist parameter is necessary for certain plugins only, e.g. the Python adapter.
  - `configure --enable-comin`

* - (plugin_list-comm)=
    **comm**
  - C
  - ""
  -
  - Name of the MPI communicator used in the second MPI Handshake. (At the first MPI Handshake "comin" is used). This namelist parameter is necessary only when using ComIn plugins in combination with external processes.
  - `configure --enable-comin`
```



(ref_buildrun_nml_vlev_distr)=
# Information on vertical level distribution

The atmospheric model needs hybrid vertical level information (i.e. the so called vertical coordinate tables `vct_a`,
`vct_b` specifying the distribution of coordinate surfaces) to generate the terrain following height based coordinates. The 1D fields `vct_a`, `vct_b` are created within ICON
during the setup phase, given that no input file is provided (empty vct_filename).
For the SLEVE vertical coordinate (`ivctype=2`), the creation  of `vct_a`, `vct_b` is controlled by the
Namelist [sleve_nml](ref_buildrun_nml_sleve_nml) together with the parameter num_lev.
For the Gal-Chen vertical coordinate (`ivctype=1`), the user has only very limited control regarding its ICON internal creation.
It is e.g. possible to create an equidistant level distribution for idealized testcases, by specifying the parameters `layer_thickness`
and n_flat_level. For more general grids, it is recommended to read the vertical coordinate tables from file.
Example files and information on the required format can be found in `<icon home>/vertical_coord_tables`, as well as in the ICON tutorial.
Note that for the SLEVE coordinate, only `vct_a` must be provided in the input file. It is recommended to set `vct_b` to zero.



# Compile flag for mixed precision

To speed up code parts strongly limited by memory bandwidth (primarily the dynamical core and the tracer advection),
an option exists to use single precision for variables that are presumed to be insensitive to computational accuracy.
This affects most local arrays in the dynamical core routines (solve_nonhydro and velocity_advection), some local
arrays in the tracer transport routines, the metrics coefficients,
arrays used for storing tendencies or differenced fields (gradients, divergence etc.), reference atmosphere fields,
and interpolation coefficients. Prognostic variables and intermediate variables affecting the accuracy of mass conservation
are still treated in double precision. To activate the mixed-precision option, run the configure script with the
`--enable-mixed-precision` flag.




# Appendix A: Arithmetic expression evaluation

The {{ '[externals/fortran-support/src/mo_expression.F90]({}/externals/fortran-support/src/mo_expression.F90)'.format(base_url) }} module evaluates basic arithmetic
expressions specified by character-strings.

It is possible to include mathematical functions, operators, and
constants.
An application of this module is the evaluation of arithmetic
expressions povided as namelist parameters.

Besides, Fortran variables can be linked to the expression and used in
the evaluation.
The implementation supports scalar input variables as well as 2D and
3D fields.

From a users' point of view, the basic usage of this module is
described in subection below.
Technically, infix expressions are processed based on a Finite State
Machine (FSM) and Dijkstra's shunting yard algorithm.
A more detailed described of the Fortran interface is given in
Section [usage with fortran](ref_buildrun_nml_fortran).


(ref_buildrun_nml_arithmetic)=
## Examples for arithmetic expressions


Basic examples:

 - `sqrt(2.0)`
 - `sin(45*pi/180.) * 10 + 5`
 - `if(1. > 2, 99, -1.*pi)`
 - `min(1,2)`


Variables are used with a bracket notation:

 - `sqrt([u]^2 + [v]^2)`

Note that the use of variables requires that these are enabled
("linked") by the Fortran routine that calls the
{{ '[externals/fortran-support/src/mo_expression.F90]({}/externals/fortran-support/src/mo_expression.F90)'.format(base_url) }} module.



## Expression syntax


### List of functions


```{list-table}
:header-rows: 1

* - name
  - #args
  - description

* - `log()`, `exp()`
  - 1
  - natural logarithm and its inverse function.

* - `sin()`, `cos()`
  - 1
  - trigonometric functions

* - `sqrt()`
  - 1
  - square root

* - `erf()`
  - 1
  - Gauss error function

* - `min()`, `max()`
  - 2
  - minimum and maximum of two values

* - `if(value,then,else)`
  - 3
  - conditional expression (`value` > 0.`)
```




### List of operators


```{list-table}
:header-rows: 1

* - name
  - evaluates to

* - `a + b`, `a - b`, `a * b`, `a / b`
  - {math}`(a+b)`, {math}`(a-b)`, {math}`(a*b)`, {math}`(a/b)`

* - `a ^ b`
  - {math}`a^b`

* - `a > b`
  - `1` if `a > b`, `0` otherwise

* - `a < b`
  - `1` if `a < b`, `0` otherwise
```





### List of available constants


```{list-table}
:header-rows: 1

* - name of constant
  - assigned value
  - description

* - `pi`
  - {math}`4 \, \mathrm{atan}(1)`
  - mathematical constant equal to a circle's circumference divided by its diameter

* - `r`
  - {math}`6.371229\cdot 10^6`
  - Earth's radius. This number seems to be based on Hayford's 1910 estimate of the Earth. It is used in ICON as well as MPAS and was almost certainly taken from the Jablonowski and Williamson test case (QJRMS, 2006).
```



(ref_buildrun_nml_fortran)=
## Usage with Fortran


The minimal Fortran interface is as follows:

1. The `TYPE expression` which is initialized with the
   character-string that specifies the arithmetic expression.
2. The type-bound procedure `evaluate()`, which returns the
   result (scalar or array-shaped) as a `POINTER`.
3. The type-bound procedure `link()` connecting a variable
   to a name in the character-string expression.



### Fortran examples


The following examples illustrate the arithmetic expression parser.
The calls to `DEALLOCATE` the data structures have been
ommitted for the sake of brevity:

1. Scalar arithmetic expression:

```
formula = expression("sin(45*pi/180.) * 10 + 5")
CALL formula%evaluate(val)
  ... use "val" for some purpose ...
```

2. Masking of a 2D array as an example for the `link` procedure:

```
formula = expression("if([z_sfc] > 2., [z_sfc], 0. )")
CALL formula%link("z_sfc", z_sfc)
CALL formula%evaluate(val_2D)
  ... use "val_2D(:,:)" for some purpose ...
```



### Error handling


Invalid arithmetic expressions yield "empty" expression objects. When these
are evaluated, a `NULL()` pointer is returned.

A successful expression evaluation can be tested with the `err_no` variable:

```
IF (formula%err_no == ERR_NONE) THEN
  ...
END IF
```

In case of error, the `err_no` variable also provides the
reason for the aborted evaluation process.



## Remarks


 - Variable names are treated case-sensitive!
 - For 3D array input it is implicitly assumed that 2D fields are
   embedded in 3D fields as "3D(:,level,:) = 2D(:,:)".



# Appendix B: Changes incompatible with former versions of the model code

This page documents namelist changes that are incompatible with earlier
versions of the model code. Incompatible changes include renaming
namelist parameters or changing their type, removing parameters,
changing default settings, changing parameter scope, or introducing
new cross-check rules.

```{list-table}
:header-rows: 1

* - Update date
  - Parameter
  - Commit
  - Description

* - 2025-07-31
  - aes_cop_nml
  - `icon-mpim:master xxx`
  - Remove `aes_cop_config(:)%cinhoms` option: the cloud inhomogeneity factor for snow has been removed for the AES physics.

* - 2025-02-10
  - lnd_nml
  - `icon-nwp:master xxx`
  - Remove option `itype_interception=2`: the namelist switch `itype_interception` is still accepted, but the model
    now aborts when it is set to a value different from `1`.

* - 2023-12-11
  - dynamics_nml
  - `icon-nwp:master a1a19a69`
  - Remove Namelist switch `iequations`: this switch became obsolete, as there exists only one set of governing equations for the atmosphere and ocean each.

```
