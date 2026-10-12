```{eval-rst}
:orphan:
```

(ref_buildrun_input_data)=
# Input Data

(ref_buildrun_grid_files)=
## Grid Files

The ICON model receives information about the horizontal grid from so-called **grid files** in the [NetCDF format](https://www.unidata.ucar.edu/software/netcdf/).
These files store coordinates and topological index relations between cells, edges and vertices of the chosen domain.
A detailed description of the content of these grid files is provided in the _Necessary Input Data_ section of the **{term}`ICON Tutorial`**.

The grid files for ICON usually follow the nomenclature `R<n>B<k>`, where `<n>` denotes the number of root divisions and `<k>` the number of subsequent bisections.
From `<n>` and `<k>` the resolution of the grid can be estimated by the formula:

```{math}
\Delta x \sim \frac{5050}{n \cdot 2^k} \quad km.
```

### Grid Point Search

There are many options to find the nearest grid point to given latitude/longitude coordinates.
Here, we provide an example which makes use of
[`KDTree`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html)
provided by [`SciPy`](https://docs.scipy.org/doc/scipy/index.html).

:::{dropdown} Grid Point Search using Python/SciPy
```python
import numpy as np
from scipy.spatial import KDTree
import xarray

def lonlat2xyz(lon, lat):
    clat = np.cos(lat)
    return clat * np.cos(lon), clat * np.sin(lon), np.sin(lat)

# Example grid and coordinates
gridFile = 'icon_grid_0026_R03B07_G.nc'
lon = np.deg2rad(8.7708333)
lat = np.deg2rad(50.099444)

grid = xarray.open_dataset(gridFile)
tree = KDTree(np.stack(lonlat2xyz(grid.clon.values, grid.clat.values), axis=-1))
NNdistance, NNindex = tree.query(np.stack(lonlat2xyz(lon,lat), axis=-1), k=1)

print(f'The closest cell to lon={np.rad2deg(lon):.4f} lat={np.rad2deg(lat):.4f} is:')
print(f'{NNindex=} lon={np.rad2deg(grid.clon[NNindex].values):.4f} lat={np.rad2deg(grid.clat[NNindex].values):.4f}')
```
:::

(ref_buildrun_external_param)=
## External Parameters (NWP)

{material-regular}`warning;2em;pst-color-secondary` _Please note that this description applies to the [](ref_atmosphere_nwp_physics)_.

External parameter datasets contain topological and climatological data that is assumed to be constant during a typical NWP integration.
These datasets are aggregated to a given ICON grid using the **[EXTPAR Software](http://www.cosmo-model.org/content/support/software/default.htm)**.
Further information is available in the **[EXTPAR Documentation](https://c2sm.github.io/extpar/)** and in the _Necessary Input Data_ section of the **{term}`ICON Tutorial`**.

### Obtaining Grid & External Parameter Files

Currently, there are two options to obtain [grid](ref_buildrun_grid_files) and [external parameter](ref_buildrun_external_param) data:

- A set of predefined grid and external parameter datasets is available at **[icon-downloads.mpimet.mpg.de/](http://icon-downloads.mpimet.mpg.de/)**.

- For users, who want to specify a custom domain, the **[Zonda Webinterface](ref_tools_gridextpargui)** provides all relevant options.

(ref_buildrun_ini_bound_data)=
## Initial & Boundary Data

Besides horizontal grid files and external parameters, ICON needs data describing
the initial state of the component to run. NWP runs require data for atmosphere,
land and sea. When running ICON in limited-area mode also lateral boundary data have
to be provided in regular time intervals.

ICON can take data from DWD's Data Assimilation Coding Environment (DACE), from
its own forecasts, and data interpolated from IFS forecasts or analysis.
Depending on which data is taken, several steps are necessary to process these
data in a way that they can be read by ICON.

* [Data Assimilation System]
* [ICON forecasts]
* [IFS analysis or forecasts](ref_buildrun_icbcifs)


(ref_lateral_boundary_update)=
## Lateral Boundary Update

{{author}}D. Reinert and G. Zängl{{endauthor}}

The lateral boundary zone update (zone 0 in Figure 6.1 of the {term}`ICON Tutorial`) follows the parent-child coupling approach for nesting, described in Section 3.9.1 of the {term}`ICON Tutorial`.

Let {math}`\psi_{bc}^{k} \in \{v_{n}, w,\rho,\theta_{v},q_{k}\}` denote the set of boundary data, which is provided externally at regular time intervals {math}`\Delta t_{bc}`. The index {math}`k` denotes discrete times for which a boundary data set is available, according to the definition {math}`t^{k}=t^{0} + k\Delta t_{bc}`. Further, let

```{math}
 \frac{\delta \psi_{bc}}{\delta t}= \frac{\psi_{bc}^{k+1} - \psi_{bc}^{k}}{\Delta t_{bc}}
```

define their discrete time tendency, which is assumed constant over the interval {math}`t^{k}\leq t < t^{k+1}`. In continuous form, the model state inside the lateral boundary zone  {math}`\psi_{lbz}(t)` is derived from

```{math}
:label: eq_lam_boundary_state
  \psi_{lbz}(t) = \psi(t^{0}) + \int\limits_{t^{0}}^{t} \frac{\partial \psi_{bc}}{\partial t} \,\mathrm{d}t
```

where {math}`\psi(t^{0})` is the initial condition originating from the first-guess or analysis. Equation {eq}`eq_lam_boundary_state` ensures that the model state {math}`\psi_{lbz}(t)` interpolates the external boundary data at times {math}`t=t^{k}` and transforms piecewise linearly between consecutive boundary data sets for intermediate times {math}`t^{k}< t < t^{k+1}`. If written in discrete form and applied repeatedly over {math}`n` dynamics time steps of size {math}`\Delta \tau` (or the fast physics time step {math}`\Delta t` in case of {math}`q_{k}`), the update formula {eq}`eq_lam_boundary_state` becomes

```{math}
:label: eq_lam_boundary_update
 \psi_{lbz}(t^{n+1}) = \psi_{lbz}(t^{n}) + \Delta \tau\, \frac{\delta \psi_{bc}}{\delta t}
```

Hence, the lateral boundary zone is updated _incrementally_ by adding boundary data
time tendencies, rather than applying the boundary data {math}`\psi_{bc}` directly.

### Known Pitfalls

Even though this approach is largely consistent with the boundary update approach for nesting, it can lead to unwanted side effects if used without care.

Upon inspection of Equation {eq}`eq_lam_boundary_state` it becomes clear that this approach relies on the implicit assumption

```{math}
\psi(t^{0}) = \psi_{bc}(t^{0})\,,
```

stating that the external boundary data at model start {math}`\psi_{bc}(t^{0})` match the model's initial conditions {math}`\psi(t^{0})`. Otherwise, Equation {eq}`eq_lam_boundary_state` and its discrete counterpart {eq}`eq_lam_boundary_update` will no longer interpolate external boundary data for {math}`t=t^{k}`, resulting in {math}`\psi_{lbz}(t^{k}) \neq \psi_{bc}(t^{k})`. Values in the boundary zone will deviate by a constant mismatch  {math}`\Delta \psi=\psi(t^{0}) - \psi_{bc}(t^{0})`, such that

```{math}
 \psi_{lbz}(t^{k}) = \psi_{bc}(t^{k}) + \Delta \psi\,, \quad \text{for } k\geq 1\,.
```

This error is persistent. It will not decay during the course of the simulation, and it will also be present at intermediate times {math}`t^{k}< t < t^{k+1}`. Several studies have shown that this error can introduce significant pressure biases throughout the entire model domain, thereby leading to a notable degradation of the overall forecast quality.

Significant mismatches will occur whenever initial conditions and boundary conditions originate
from different model runs or models. For example, mismatches occur when initial conditions are taken from ICON, while lateral boundary conditions are taken from IFS (or vice versa), or when ICON initial and boundary conditions from different resolutions or domains are mixed.

The recommended way to avoid these mismatches is to set the Namelist parameter `init_latbc_from_fg=.TRUE. (limarea_nml)` unconditionally. This ensures that the first boundary data set at {math}`t=t^{0}` is copied from the initial conditions file, thereby ensuring {math}`\psi(t^{0}) = \psi_{bc}(t^{0})`.

:::{admonition} Important note on lateral boundary conditions
:class: admonition-icontheme
Deviations of the boundary data at the model start date from the initial conditions produce significant surface pressure biases and degrade forecast quality throughout the simulation. To ensure consistent lateral boundary conditions, set the namelist parameter `init_latbc_from_fg=.TRUE. (limarea_nml)` for any type of limited area simulations. Lateral boundary conditions at the model start date are then copied from the analysis or first guess.

A lateral boundary data file must be provided for the start date of the experiment (see the Namelist parameter `experimentStartDate`), and must contain the full set of lateral boundary fields. The required fields include the 3D height field `HHL`, if the driving model uses a height-based vertical coordinate. The same requirement applies when `init_latbc_from_fg=.TRUE.`, regardless of whether the nominal model start date is shifted back in time (`dt_shift<0`, namelist `initicon_nml`) or not. For non-shifted runs (`dt_shift=0`) the boundary fields in this file are not actually used. Instead, they are replaced by the initial conditions during model initialization. Nevertheless, the file is still required for the following reasons:

* The 3D height information of the driving model (`HHL`) is read from this file by ICON. This information is required to remap the boundary data vertically onto the target grid (see Section 2.3 of the {term}`ICON Tutorial`).
* ICON inspects the file contents (e.g. variable names) to validate the completeness of the dataset. This process also identifies the necessary post-processing steps to convert the data into ICON's set of prognostic variables (see Section 6.4 of the {term}`ICON Tutorial`).
:::
