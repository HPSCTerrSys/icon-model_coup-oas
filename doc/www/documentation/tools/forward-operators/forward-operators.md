```{eval-rst}
:orphan:
```

(ref_tools_fwo)=
# Forward operators (EMVORADO, RTTOV)

Forward operators are diagnostic tools to produce synthetic observations ("model equivalents") based on the simulated model state. They enable a direct comparison to the corresponding real observations and are commonly applied for data assimilation, for model verification and for forecast products in observation space (e.g., radar reflectivity maps, satellite brightness temperatures and reflectances).

In the most simple case of an in-situ observation of a prognostic state variable (e.g., a SYNOP temperature observation), the forward operator might be just a spatial interpolation operation. Much more involved are forward operators for ground- or space-based remote sensing instruments which do not directly observe state variables, but whose observables depend on the atmospheric state in a complicated and non-linear way. ICON offers some of those more complicated forward operators, which are directly ("online") coupled to the ICON model and may optionally be built and linked to the ICON executable and switched on.


## The **E**fficient **M**odular **VO**lume scan **RAD**ar **O**perator (EMVORADO)

EMVORADO computes synthetic radar volume scans
* radial winds
* horizontally polarized radar reflectivity factor

and dual-polarization parameters
* differential reflectivity
* horizontal attenuation coefficient
* specific differential phase shift
* total differential phase shift
* cross-correllation coefficient
* linear depolarization ratio
* differential attenuation

for several (many) ground-based radar stations distributed across the model domain in one go. It is modular in the sense that it offers different options of different computational complexity and accuracy for each involved physical process (scattering theory, atmospheric beam propagation and beam broadening), to enable a balancing of cost and accuracy for individual applications. It can be used for 3D radar data assimilation in DWD's ICON-KENDA LETKF system as well as for model verification, evaluation and visualization in radar observation space.

EMVORADO is COSMO-Software and has been developed by DWD and KIT. It is distributed with ICON as an external submodule under the ICON open-source license.

```{figure} emvorado_bspl.png
:align: center
:height: 350
:width: 450
<a name="#fig1-emvorado">Figure 1: example of simulated reflectivity volume scans for some German radar stations. For better visibility, we show only one elevation per station and only 5 stations, but EMVORADO is able to simulate much more elevations and stations in one go.
```

### How to enable building and linking to ICON:

To enable the building and linking of EMVORADO, the flag `--enable-emvorado` needs to be added in the [configure wrapper](ref_buildrun_configuration_wrappers) (see also [here](ref_buildrun_configuration_icondep)):

```sh
./configure ... --enable-emvorado ...
```


### How to switch it on:

EMVORADO may be applied in two different ways in ICON:

* **Full-fledged**: simulate volume scans of radar data for several stations / entire radar networks for realistic scan geometries, i.e., in polar coordinates range, azimuth and elevation centered around each radar station at the earth's surface and output them in netcdf, grib or ASCII format. This can be done for real case runs, but also for idealized runs with purely synthetic radar stations. It is also possible to produce so-called composites of simulated and observed PPI-scans in grib-format and feedback files for ICON KENDA radar data assimilation. For the latter, EMVORADO is able to ingest and process real observations.

  * See **{term}`EMVORADO User's Guide`** for instructions on how to set up the corresponding namelists for the different modes of operation.

* **Traditional grid point output**: enhance the traditional radar reflectivity output of the fields dbz, dbz_850, dbz_cmax, dbz_ctmax, echotop, echotop_in_m by advanced MIe- and T-matrix scattering methods provided by EMVORADO modules. This option may be configured by ICON's /synradar_nml/ namelist.

  * See the section on /synradar_nml/ in ICON's **Namelist_overview.pdf**.


### Documentation

* **{term}`EMVORADO User's Guide`**: [Download PDF from the COSMO webpage](https://www.cosmo-model.org/content/model/documentation/core/emvorado_userguide.pdf) or get the PDF directly from the submodule: **externals/emvorado/DOC/TEX/emvorado_userguide.pdf**. The User's Guide contains further references to scientific publications about EMVORADO.
* **Scientific documentation of the scattering computation part in {term}`COSMO Technical Report No. 28`**: [Download PDF from COSMO webpage](https://www.cosmo-model.org/content/model/cosmo/techReports/docs/techReport28.pdf)
* **Enhancement of the traditional grid point output**: section on /synradar_nml/ in ICON's **Namelist_overview.pdf**.
