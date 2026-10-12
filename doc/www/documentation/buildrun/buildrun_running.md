```{eval-rst}
:orphan:
```

(ref_buildrun_running)=
# Running ICON

The ICON model is typically run through a *runscript*, which sets up the working directory, populates it with all required input files (grid files, namelists, etc.), sets environment variables, runs the model, and postprocesses its output.

You can generate a runscript with one of the two tools described in this document. The first tool [`mkexp`](ref_tools_mkexp) is well documented and continuously maintained. The only downside is that it does not have support for many environments and experiments that you can run with ICON yet. The second one [`make_runscript`](ref_buildrun_make_runscript) is a legacy, poorly documented set of shell scripts, which, however, supports a lot of (mainly HPC) environments and experiments.

(ref_buildrun_mkexp)=
## Using `mkexp` to prepare ICON experiments

[MakeExperiments! (`mkexp`)](ref_tools_mkexp) is a Python tool for preparing experiments with ICON. It helps users set up an experimental workflow, and generate the runscript needed to execute supported configurations. The tool presently supports setting up ICON configurations in the DKRZ environment (CPU and GPU) and is and can be adapted to other environments.

You can find the `mkexp` tool in the `utils/mkexp` directory of the ICON source code repository, which is managed as a git submodule.

### Requirements

Before you start using `mkexp`, make sure that your software environment meets the following requirements:

1. You need to make sure that [Jinja2](https://pypi.org/project/Jinja2/) and [six](https://pypi.org/project/six/) are available in your Python environment.
   Another option is installing `mkexp` in your local or virtual Python environment with

    ```sh
    python3 -m pip install utils/mkexp
    ```

1. If ICON is configured and built [out-of-source](ref_buildrun_configuration_oos), you need to set the `MKEXP_PATH` environment variable to include the absolute path of the "`run`" subdirectory in your root build directory, like

    ```sh
    export MKEXP_PATH=.:/path/to/icon/build/directory/run
    ```

### Steps to run an experiment

Running an experiment using `mkexp` generally consists of three steps briefly described below in this section. For more details, see the [mkexp docs](ref_tools_mkexp).

**Step 1: Create the configuration file**

The configuration file overrides the generic parameters specified in `run/mkexp/defaults/DEFAULT.config`. The experiment- and environment-specific parameters can be overridden or set via the `EXP_TYPE` and `ENVIRONMENT` variables of the configuration file. The values of those variables correspond to the `.config` files in the `run/mkexp/types` and `run/mkexp/environments` directories. The configuration file of an experiment can override or set additional parameters via the `EXP_OPTIONS` variable (the value corresponds to the files in the `run/mkexp/options` directory or directly. The resulting set of parameters is used by `mkexp` to create the experiment runscript, set the run environment and configure the required directory structure for the run.

The most simple way to create the configuration file for an experiment is to copy one of the `/run/examples` (e.g. `/run/examples/bubble.config`) and adjust it to your needs. First, you need to come up with the experiment identifier (e.g. `exp_id`) and copy the example configuration to the {{ '[`run`]({}/run)'.format(base_url) }} directory under the name that corresponds to that identifier:

```sh
cd ./run
cp ./examples/bubble.config ./exp_id.config
```

You may now review and edit the contents of `exp_id.config`. It should run on most personal computing devices without changes.

Running the experiment in an HPC environment may require adjustments. For example, if you intend to run the experiment on the DKRZ machine, make sure the `ACCOUNT` is set correctly (it should be set to the SLURM account that you normally submit jobs with):

```diff
@@ -15,2 +15,3 @@
 EXP_TYPE = torus
+ACCOUNT = xy1234

```

**Step 2: Generate the scripts and workflow environment**

Execute [`mkexp`](ref_tools_mkexp):

```sh
../utils/mkexp/mkexp exp_id.config
```

The command will create the required directory structure. For example:

```sh
Script directory: '/path/to/icon-srcdir/experiments/exp_id/scripts'
Data directory: '/work/your-account-number/your-user-number/master/experiments/exp_id/outdata'
Work directory: '/scratch/your-account-type/your-account-number/master/experiments/exp_id/work'
Log directory: '/work/your-account-number/your-user-number/master/experiments/exp_id/log'
```

The command will also create the `exp_id.run_start` runscript and place it to the `Script directory`. Review the script to make sure that the path to the grid file (`icon_grid_G.nc`) is set correctly for your environment.

**Step 3: Execute the runscript.**

For the last step, switch to the aforementioned `Script directory` and either execute or submit (if you are in the HPC environment) the runscript:

```sh
cd ../experiments/exp_id/scripts
sbatch exp_id.run_start
```

(ref_buildrun_make_runscript)=
## Using `make_runscript` to prepare ICON experiments

The `make_runscript` shell script is an ICON-specific tool for runscript generation. It takes the experiment template files from the {{ '[`run`]({}/run)'.format(base_url) }} directory, prepends the environment-specific shell snippets from the `run/create_target_header` file, adjusts the result based on how ICON is configured and built (see `run/collect.set-up.info.in`) and produces a shell script that is ready for the execution or submission.

### Requirements

The tool does not support the [out-of-source](ref_buildrun_configuration_oos) builds. To circumvent this, most of the existing [configure wrappers](ref_buildrun_configuration_wrappers) make the required files available in the build directory:

```sh
# Copy runscript-related files when building out-of-source:
if test $(pwd) != $(cd "${icon_dir}"; pwd); then
  rsync -uavz ${icon_dir}/run . --exclude='*.in' --exclude='.*' --exclude='standard_*' --exclude=mkexp
  ln -sf -t run/ ${icon_dir}/run/{standard_*,mkexp}
  for dir in \
    'externals/art/runctrl_examples' \
    'externals/ecrad/data' \
    'externals/jsbach/data'
  do
    src="${icon_dir}/${dir}"
    test -d "${src}" && mkdir -p "${dir}" && rsync -uavz "${src}/" "${dir}"
  done
  rsync -uavz ${icon_dir}/make_runscripts .
  ln -sf ${icon_dir}/data
  ln -sf ${icon_dir}/vertical_coord_tables
fi
```

### Steps to run an experiment

To generate a runscript based on a particular template (e.g. `run/exp.atm_tracer_Hadley`), switch to the root build directory of ICON and run `make_runscript` while providing the name of the experiment (without the prefix) as an argument:

```sh
./make_runscripts atm_tracer_Hadley
```

Alternatively, you can generate runscripts for all existing experiments:

```sh
./make_runscripts --all
```

The generated runscripts are saved to the {{ '[`run`]({}/run)'.format(base_url) }} subdirectory of the build directory. The headers of the runscripts containing arguments for the HPC workload manager, e.g. [SLURM](https://slurm.schedmd.com/), might require additional manual adjustments regarding CPU time accounting, node allocation, etc.

Once the runscript is created and adjusted, switch to the {{ '[`run`]({}/run)'.format(base_url) }} subdirectory of the root build directory of ICON and either execute or submit (if you are in the HPC environment) the runscript:

```sh
cd ./run
sbatch ./exp.atm_tracer_Hadley.run
```

:::{admonition} checksuite.nwp specifics
:class: admonition-icontheme
Experiments in the `run/checksuite.nwp` directory are prepared to be run in the hybrid mode, which uses two icon binaries. To generate those experiments the user can specify the `-r <run directory to process>` argument:
```sh
./make_runscripts --all -r run/checksuite.nwp
```
In the case of the NEC-Aurora, the second binary is the host (or x86 scalar) binary. It can be found automatically if the two build directories are either `something/vector` and `something/host` or `something_else/VH` and `something_else/VE`. In either case, `./make_runscripts` must be called in the `vector` (or `VE`) build directory. Alternatively, the x86-host binary can be specified by the `-s <secondary build dir>` option:
```sh
./make_runscripts run_ICON_01_R3B9_lam -r run/checksuite.nwp -s ../host_gcc-9.1.0/
```
:::

:::{admonition} make_target_runscript
:class: admonition-icontheme
Alternatively, the users can employ the low-level tool `/run/make_target_runscript` for runscript generation offering more fine-grained control over certain parameters. For example, the wall clock limit and the number of allocated nodes can be injected into the runscript as follows:
```sh
cd ./run && ln -sf ./checksuite.ocean_internal/omip/exp.ocean_omip_long exp.ocean_omip_long
./make_target_runscript in_script=exp.ocean_omip_long in_script=exec.iconrun \
  EXPNAME=ocean_omip_long cpu_time=08:00:00 no_of_nodes=20
```
:::


(ref_buildrun_gridextpar)=
## Grids & External Parameters

(ref_buildrun_grids)=
### Grid Files

The ICON model receives information about the horizontal grid from so-called **grid files** in the [NetCDF format](https://www.unidata.ucar.edu/software/netcdf/).
These files store coordinates and topological index relations between cells, edges and vertices of the chosen domain.
A detailed description of the content of these grid files is provided in the _Necessary Input Data_ section of the **{term}`ICON Tutorial`**.

The grid files for ICON usually follow the nomenclature `R<n>B<k>`, where `<n>` denotes the number of root divisions and `<k>` the number of subsequent bisections.
From `<n>` and `<k>` the resolution of the grid can be estimated by the formula:

```{math}
\Delta x \sim \frac{5050}{n \cdot 2^k} \quad km.
```

#### Grid Point Search

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


### External Parameters (NWP)

{material-regular}`warning;2em;pst-color-secondary` _Please note that this description applies to the [](ref_atmosphere_nwp_physics)_.

External parameter datasets contain topological and climatological data that is assumed to be constant during a typical NWP integration.
These datasets are aggregated to a given ICON grid using the **[EXTPAR Software](http://www.cosmo-model.org/content/support/software/default.htm)**.
Further information is available in the **[EXTPAR Documentation](https://c2sm.github.io/extpar/)** and in the _Necessary Input Data_ section of the **{term}`ICON Tutorial`**.

### Obtaining Grid & External Parameter Files

Currently, there are two options to obtain [grid](ref_buildrun_gridextpar) and [external parameter](ref_buildrun_external_param) data:

- A set of predefined grid and external parameter datasets is available at **[icon-downloads.mpimet.mpg.de/](http://icon-downloads.mpimet.mpg.de/)**.

- For users, who want to specify a custom domain, the **[Zonda Webinterface](ref_tools_gridextpargui)** provides all relevant options.

(ref_buildrun_icbc)=
## Initial & Boundary Data

## Using `hiopy` as an io component for ICON experiments
[Hiopy](https://dkrz-sw.gitlab-pages.dkrz.de/hiopy/index.html) can be used to write output from ICON simulations. Currently, the [Zarr](https://zarr.dev/) data format is supported and the output can be written in regular lat-lon grid, unstructured grid or the healpix grid.

Pre-requisites:
1. Mkexp (for configuration)
2. YAC built with a compiler which supports C++20 (eg: gcc 13)
3. Python 3.11 or greater
4. MPI 4.1.2 or greater

Supported environments via mkexp:
1. Levante-cpu
2. Dolpung-hybrid (cpu-gpu)

Current limitations for the experiment:
1. Only proleptic-gregorian calendar is supported
2. Only the ICON's native names (from add_var) is supported
3. Using hiopy parallel to the icon' async-output (pio-type 1) is not possible via mkexp

If you are interested in running hiopy without mkexp, refer to the [Hiopy](https://dkrz-sw.gitlab-pages.dkrz.de/hiopy/index.html) documentation or get in touch with its developers.

### Setting up hiopy

For HPC systems, the right modules will have to be known prior to building Hiopy alongside the ICON executable.
Issues could arrise if the HPC system limits the use of the required C++20 support in its compilers or uses older Python versions.
For further details on how to navigate potential limitations, visit [hiopy-installation](https://dkrz-sw.gitlab-pages.dkrz.de/hiopy/getting_started.html#installation).

#### Install YAC to get a $PATH_TO_YOUR_YAC_PKG_CONFIG
Build ICON with `--prefix ${PREFIX_PATH} --enable-bundled-python=mtime,yac` specified to an installation path of choice
The --enable-bundled-python=mtime,yac ensures that YAC and its dependencies are built with the
required -fPIC flag to allow Hiopy to use it.
```
cd $ICON_BUILD_DIR/externals/yac
make install # installs to ${PREFIX_PATH}
export PKG_CONFIG_PATH=${PREFIX_PATH}/lib/pkgconfig:${PKG_CONFIG_PATH}

# hiopy scripting supports only the use of python venv
python -m venv $MY_VENV # or use your own environment with >py3.11
source $MY_VENV/bin/activate

# ensure MPI is correctly set (or load appropriate modules on the HPC)
export MPI_ROOT=$MY_MPI_PATH
export CC="${MPI_ROOT}/bin/mpicc"
export CXX="${MPI_ROOT}/bin/mpicxx"

# clone and install hiopy
# with additional dependencies to configure hiopy
pip install git+https://gitlab.dkrz.de/dkrz-sw/hiopy isodate

# xarray and netcdf4 are needed only if you need output on model grid
pip install xarray netcdf4
```

#### Steps to configure an experiment with hiopy

Hiopy has a dedicated section which controls the following parameters:
1. Job: number of workers, pre-processing step for initialising the dataset
2. Data request: grid-type, desired hierarchy in case of healpix, paths to grid files in case of unstructured grids
3. User parameters: desired path to the dataset, virtual-env with hiopy installation

Example config files are under `run/examples/hiopy/*.config`

Running mkexp on it would create a pre-processing (.pre extension) python script alongside the regular run_start and run scripts. The pre-processing script can be used to initialise the Zarr store by running `./$EXP_NAME.pre`.

Once the Zarr store is initialised, the run scripts and their restart mechanism will ensure that they are filled up incrementally until the simulation ends. In case there are overlaps in time between restarts for any reason, the data will be over-written automatically.
