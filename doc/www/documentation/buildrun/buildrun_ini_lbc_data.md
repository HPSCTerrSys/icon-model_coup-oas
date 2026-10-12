```{eval-rst}
:orphan:
```

(ref_buildrun_icbcifs)=
# Initial and Boundary Data from IFS

{material-regular}`warning;2em;pst-color-secondary` _Please note that the following description is targeted at [NWP mode](ref_atmosphere_nwp_physics) and might not fully apply to other modes of ICON._

ICON runs can be driven by analysis and / or forecast files from the ECMWF model IFS.
There are two sources for IFS data:

MARS:
: If you have a user account on ECMWF computational resources you can extract
  IFS fields from the Meteorological Archival and Retrieval System (MARS).
  Documentation for MARS can be found on the
  [ECMWF server](https://confluence.ecmwf.int/display/UDOC/MARS+user+documentation).
  Also interesting is the [Web MARS](https://apps.ecmwf.int/mars-catalogue) application,
  with which you can browse the contents of the archive.
  Please see the explanations for using the script [`mars4icon_smi`](ref_buildrun_icbcifs_mars4icon_smi) below.

ERA:
: Data from the ECMWF Re-Analysis projects are publicly available and can be accessed
  through the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/datasets). The latest
  re-analysis ERA5 dates back to the year 1940. Alas, access to these data is not trivial at all.
  Downloading data requires registration and installation of the CDS API. For more informations
  see the [official Copernicus Website](https://cds.climate.copernicus.eu)

ERA(5) for DKRZ Users:
: On Levante global ERA5 data is available in original ECMWF format (`/pool/data/ERA5`).
  See the `README_ERA5_POOL_DATA_v20240719.txt` for explanations how the data is organized.

The IFS data have to be converted to data that can be read by ICON.
There are basically two different tasks to do this conversion:

- Adjusting IFS field and NetCDF attribute names to ICON conventions.
- Horizontal remapping for all fields from the IFS grid to the target ICON grid.

{material-regular}`warning;2em;pst-color-secondary` _When processing data from IFS/ECMWF,
ICON only can read NetCDF format. The remapping process therefore has to write the data in NetCDF._

The necessary steps that have to be taken are described in the next sections.
See also Section 2.2.2 of the {term}`ICON Tutorial`, from which most of the
information here is taken.

## Necessary Analysis Fields

The tables below lists all fields necessary for the initial data set.
First the data necessary for the atmosphere, then data for the soil and surface.
The first three columns list the short names according to ECMWF, DWD, and ICON
conventions, resp. These short names are important for adjusting the field names
below.

:::{table} Necessary Atmospheric Fields for Initial Data Set. <span style="color:dodgerblue">Blue fields</span> are optional.
:width: 65
:widths: auto
:align: center

| ECMWF        | DWD           | ICON          | Unit          | Description  |
| :----------- | :------------ | :------------ | :------------ | :----------- |
| u, v         | U, V          | U, V          | m/s             | horizontal velocity components |
| omega        | W             | W             | Pa/s            | vertical velocity |
| t            | T             | T             | K               | temperature       |
| z            | FI            | GEOP_ML       | m2/s2           | geopotential      |
| qv           | QV            | QV            | kg/kg           | specific humidity |
| clwc         | QC            | QC            | kg/kg           | cloud liquid water content |
| ciwc         | QI            | QI            | kg/kg           | cloud ice content |
| <span style="color:dodgerblue">crwc</span> | QR            | QR            | kg/kg           | rain water content |
| <span style="color:dodgerblue">cswc</span> | QS            | QS            | kg/kg           | snow water content |
| lnsp         | LNPS          | LNPS          | -               | logarithm of surface pressure |
:::

:::{table} Necessary Soil and Surface Fields for Initial Data Set. <span style="color:dodgerblue">Blue fields</span> are optional.
:width: 65
:widths: auto
:align: center

| ECMWF        | DWD           | ICON          | Unit          | Description  |
| :----------- | :------------ | :------------ | :------------ | :----------- |
| <span style="color:dodgerblue">sst</span>     | SST           | SST           | K               | sea surface temperature |
| ci           | FR_ICE        | CI            | [0,1]           | sea ice cover |
| z            | FIS           | GEOP_SFC      | m2/s2           | surface geopotential |
| tsn          | T_SNOW        | T_SNOW        | K               | snow temperature |
| sd           | W_SNOW        | W_SNOW        | m of water eqv. | water content of snow |
| rsn          | RHO_SNOW      | RHO_SNOW      | kg / m3         | density of snow |
| <span style="color:dodgerblue">asn</span> | ALB_SNOW      | ALB_SNOW      | [0,1]           | snow albedo |
| skt          | SKT           | SKT           | K               | skin temperature |
| stl[1-4]     | T_SO[1-4]     | STL[1-4]      | K               | soil temperature level 1-4 |
| swvl[1-4]    | SMI[1-4]      | SMIL[1-4]     | m3 / m3         | soil moisture indes layer 1-4 |
| src          | W_I           | W_I           | m of water eqv. | water content of interception storage|
| lsm          | FR_LAND       | LSM           | [0,1]           | land-sea mask |
:::


## Necessary Boundary Fields

The next table lists all fields necessary for the lateral boundary data set.

:::{table} Necessary Fields for Lateral Boundary Data Set. <span style="color:dodgerblue">Blue fields</span> are optional.
:width: 65
:widths: auto
:align: center

| ECMWF        | DWD         | ICON       | Unit       | Description  |
| :----------- | :---------- | :--------- | :--------- | :----------  |
| u, v         | U, V        | U, V       | m/s        | horizontal velocity components |
| omega        | OMEGA       | W          | Pa/s       | vertical velocity |
| t            | T           | T          | K          | temperature       |
| qv           | QV          | QV         | kg/kg      | specific humidity |
| clwc         | QC          | QC         | kg/kg      | cloud liquid water content |
| ciwc         | QI          | QI         | kg/kg      | cloud ice content |
| <span style="color:dodgerblue">crwc</span> | QR          | QR         | kg/kg      | rain water content |
| <span style="color:dodgerblue">cswc</span> | QS          | QS         | kg/kg      | snow water content |
| z            | FIS         | GEOP_SFC   | m2/s2      | surface geopotential |
| sp           | PS          | PS         | Pa         | surface pressure |
:::

:::{admonition} The IFS Vertical Coordinate
:class: admonition-icontheme
IFS is a hydrostatic model with pressure as vertical coordinate. For vertical interpolation
to ICON levels the 3D height coordinate field for IFS must be computed. This is
done by integrating the hydrostatic equation for the geopotential {math}`{\Phi}`

```{math}
\Phi = -  \int (R_d \cdot T_v) \quad d(ln p)
```

Three fields are necessary for these calculations: the surface geopotential as starting point,
the logarithm of the surface pressure and the virtual temperature

As some IFS fields use orography in transformed space and some in grid point space, it is
important to choose consistent fields.
- For transformed space the 3D geopotential `z` (`FI` in DWD and `GEOP_ML` in ICON naming; including surface
  geopotential as lowest layer) and logarithm of surface pressure `lnsp` (`LNPS` in DWD/ICON naming)
  are archived only in analysis fields.
- In grid point space the surface pressure `sp` (`PS` in DWD/ICON naming) and the surface
  geopotential `z` (`FIS` in DWD and `GEOP_SFC` in ICON naming) are archived also in forecast fields every hour.

{material-outlined}`info;2em;pst-color-primary`
_When using analysis fields it is recommended to use `z` and `lnsp`. But analysis fields
are available only every 6 hours. For a higher temporal resolution for lateral boundaries
it is recommended to use forecast data, which is available every hour. Then `sp` and `z`
(surface geopotential) have to be used._

{material-regular}`warning;2em;pst-color-secondary` _Note that the 3D geopotential field and the surface geopotential field `z` have the same names. In GRIB format they can be distinguished by the height coordinate (type of first/second fixed surface). But the values of the lowest layer of the 3D field and the surface field differ slightly because they are represented in different spaces._
:::

(ref_buildrun_icbcifs_mars4icon_smi)=
## Script: mars4icon_smi

The ICON repository contains a script {{ '[`mars4icon_smi`]({}/scripts/preprocessing/mars4icon_smi)'.format(base_url) }}.
This script can be used to extract all fields listed above from MARS for a specific date.

The following command will retrieve all fields for the initial data listed above (even some more):

```shell
mars4icon_smi -a 70.0/-10.0/40.0/20.0 -r 1279 -l 1/to/137 -g 0.1/0.1 -d 2025100100 -O -L 1 -o ifs_2025100100.grb -p 5 -A
```

For a detailed list of options see `mars4icon_smi -h`. Options used here are:

```shell
-a North/West/South/East   area keyword of MARS to retrieve data on a regular lat/lon grid
-r resolution              spectral resolution (e.g. 511, 639, 799, 1023, 1279)
-l levellist               specifying the levels that have to be retrieved (with MARS syntax)
-g grd                     use dlon/dlat for regular lat/lon grids
-d date                    initial time in the format YYYYMMDDHH
-o grib_file               name of output grib file
-s step                    specify forecast steps (ranges start/to/stop/by/dh are also possible with option -E 0)
-A                         take geopotential / orography from MARS `type=an, step=0`
-E 0                       specify MARS parameter for soil data. The `0` indicates that NO soil data should be retrieved
-P                         retrieve surface pressure sp instead of lnsp on lowest model layer
```

It is also possible to retrieve only the necessary fields for the boundary data:

```shell
mars4icon_smi -a 70.0/-10.0/40.0/20.0 -r 1279 -l 1/to/137 -g 0.1/0.1 -d 2025100100 -s 0/to/6/by/3 -o ifs_[date]_[step]_lbc.grb -E 0 -P
```

{material-outlined}`info;2em;pst-color-primary`
_IFS still uses a mix of GRIB1 (soil and surface) and GRIB2 (atmosphere) data.
Because of different specifications for the horizontal grid this can lead to problems when
processing the grib output file from {{ '[`mars4icon_smi`]({}/scripts/preprocessing/mars4icon_smi)'.format(base_url) }}. The examples shown here do specify
the lat/lon grid in GRIB1 notations (where lon ranges from -180.0 to +180.0 degrees, so
lon west is specified with negative values). This usually works._

{material-outlined}`info;2em;pst-color-primary`
_ICON requires the soil moisture index `SMIL`, while IFS provides the volumetric soil moisture
content `SWVL`. A mathematical conversion is performed by `mars4icon_smi`, but this is not reflected
in the variable names. This has to be done when reformatting the GRIB data (see below)._

{material-outlined}`info;2em;pst-color-primary`
_Also ICON requires the geometric vertical velocity `W`, but IFS provides the pressure based
vertical velocity `OMEGA`. In the procedure described below we do change the name to `W`, but
the mathematical conversion is done in ICON._

{material-regular}`warning;2em;pst-color-secondary` _Note: Prior to 2013-06-25 12 UTC only 91 instead of 137 vertical levels were used by the operational
system at ECMWF. For more information see [changes in IFS](https://www.ecmwf.int/en/forecasts/documentation-and-support/changes-ecmwf-model)._

## Remapping IFS Data with CDO

### Reformatting GRIB data to NetCDF

When processing IFS data, ICON only reads NetCDF format. Therefore it is first necessary to reformat
the IFS GRIB data. This can be done by a `cdo` command. To map the GRIB metadata to NetCDF names,
it is important to choose proper `definition files` for the `eccodes` library. Note that there
are different definition files for ECMWF and DWD/ICON, which have different usage of `shortNames`
(see tables above).

Which definiton files are used can be set by the environment variable `ECCODES_DEFINITION_PATH`.

Empty `ECCODES_DEFINITION_PATH`
: If this environment variable is not set, `cdo` uses an intrinsic path to the ECMWF definitions
  and all NetCDF names will be `ECMWF shortNames`.

`ECCODES_DEFINITION_PATH=/path/to/definitions`
: If this environment variable is set to use only ECMWF definition files, the NetCDF names will
  also be only `ECMWF shortNames`.

`ECCODES_DEFINITION_PATH=/path/to/definitions.edzw:/path/to/definitions`
: If this environment variable is set to also use DWD definition files, the NetCDF names will
  be a mix of DWD and ECMWF shortNames.
  Note that DWD definitions _must_ precede the ECMWF definitions.

{material-regular}`warning;2em;pst-color-secondary`
_[`cdo`](ref_tools_cdo) is built with a special eccodes version. If the environment variable `ECCODES_DEFINITION_PATH`
is set, it is important to use eccodes definition files
from exactly the same version. This version can be identified with `cdo --version`._


In the following it is recommend to also use DWD definition files when doing the reformatting.
In this way at least some of the field names already have correct names.

Reformatting GRIB to NetCDF with `cdo`:

```shell
cdo -f nc copy <grb-file-name.grb> <nc-file-name.nc>
```

### Renaming Fields / Attributes in the NetCDF File

The names of fields in a NetCDF file can be changed by a `cdo` command:

```shell
cdo chname,name_old,name_new file_1.nc file_2.nc
```

It is possible to give a sequence of `name_old,name_new` lists. Note that a new file `file_2.nc`
has to be written by that command. The following block shows how the fields can be renamed
in two steps. We assume that also DWD definition files have been used when formatting the GRIB file to
NetCDF:

```shell
cdo chname,FI,GEOP_ML,z,GEOP_SFC,OMEGA,W,lnsp,LNPS,clwc,QC,ciwc,QI,crwc,QR,cswc,QS,tsn,T_SNOW,sd,W_SNOW,rsn,RHO_SNOW,asn,ALB_SNOW,skt,SKT,sst,SST <nc-file-name.nc> <nc-file-name_A.nc>
cdo chname,stl1,STL1,stl2,STL2,stl3,STL3,stl4,STL4,ci,CI,src,W_I,sr,Z0,lsm,LSM,swvl1,SMIL1,swvl2,SMIL2,swvl3,SMIL3,swvl4,SMIL4 <nc-file-name_A.nc> <nc-file-name_B.nc>
```

The resulting file `<nc-file-name_B.nc>` then should have all names according to the ICON conventions.

When processing boundary data the command for renaming the fields looks like:

```shell
cdo chname,OMEGA,W,clwc,QC,ciwc,QI,crwc,QR,cswc,QS,sp,PS,sst,SST,src,W_I,z,GEOP_SFC <lbc-file.nc> <lbc-file-new.nc>
```

For older ICON versions (before 2025.10) also a NetCDF attribute has to be renamed for initial data.
Reason is that the older versions expect the name `ncells`, but the grid files usually have the name `cell`.
To rename attributes the command `ncrename` from the NetCDF Operators `nco` has to be used:

```shell
ncrename -d cell,ncells <nc-file-name_B.nc> <nc-file-name-C.nc>
```


### Horizontal Remapping of the Data

The following `cdo` command interpolates the IFS data from the regular lat/lon grid to a
specified ICON target grid:

```shell
cdo remapcon,<path/to/target_grid.nc>:2 <nc-file-name_C.nc> <nc-file-ICON.nc>
```

If several IFS data files have to be remapped (e.g. also for boundary data), it is possible to
precalculate the interpolation weights from the source to the target grid:

```shell
cdo -P 4 gencon,<path/to/target_grid.nc>:2 <nc-file-name_C.nc> <target.weights>
```

This file with the interpolation weights can then be used for all other interpolations,
which can speed up the interpolation process.

```shell
cdo -s -r remap,<path/to/target_grid.nc>:2,<target.weights> <nc-file-name_X.nc> <target-file-final.nc>
```

{material-regular}`warning;2em;pst-color-secondary`
_Note: The workflow with precalculated weightings might fail when using the data in ICON due to missing values in the data._

ICON can now be run using the remapped data and setting the namelist variable `init_mode=2` in namelist
group `initicon_nml`.
