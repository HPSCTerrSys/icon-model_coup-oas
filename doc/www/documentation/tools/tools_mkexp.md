```{eval-rst}
:orphan:
```

(ref_tools_mkexp)=
# MakeExperiments! (mkexp)

[MakeExperiments!](https://gitlab.dkrz.de/esmenv/mkexp/-/blob/master/mkexp) is a tool for preparing experiments with MPI-M's earth system models. It provides a unified command line interface to perform experiments using ICON configurations maintained by the MPI-M.

Information on its use is provided in the _Running ICON_[^1] guide.

[^1]: `doc/www/buildrun/buildrun_running.md`, "Using `mkexp` to prepare ICON experiments"

## Basic structure

To set up an experiment, you create a config file eg.

```
EXP_TYPE = torus
[jobs]
  [[run_start]]
    time_limit = 00:05:00
```

and run _mkexp_ on that file. The configuration is then gathered from default settings, experiment type, possibly options, the computing environment settings, and finally the config file itself. Latter settings amend the former.

The final config is applied to the requested job templates, creating one script for each. Templates are sought in the same locations as configs except for the environment. The latest template found is chosen to create the script. In addition, the environment template is included at the beginning of each script.

For technical details, see {{ '[mkexp official documentation]({}/utils/mkexp/doc/mkexp.pdf)'.format(base_url) }}.

### Default settings and templates (`run/mkexp/defaults`)

* Global default configuration for all experiments (`DEFAULT.config`)
* Default job templates (`DEFAULT.*.tmpl`)
* Ancillary templates as required by formatting or job templates

### Experiment types (`run/mkexp/types`)

* Experiment specific configuration (`<EXP_TYPE>.config`)
* Specialized job templates if required (`<EXP_TYPE>.*.tmpl`)

### Packaged Options (`run/mkexp/options`)

* Collections of optional config settings per topic (`<EXP_OPTION>.config`)
* `link_output_files.config`: Create an inline `link` job which is added to each `run` job and links all output files in the work directory to the data directory. This is mutually exclusive to the standard post-processing job (`post`).

### Computing environment settings (`run/mkexp/environments`)

* Fallback settings for all environments (`DEFAULT.config`)
* Specific environment settings and job resource specifications (`<ENVIRONMENT>.config`)
* Environment template formatting information according to system needs (`<ENVIRONMENT>.tmpl`)

The run types can be set via the `ENVIRONMENT` variable. These are:

* `macOS`: run on macOS
* `levante`: run on the CPU partition of Levante (DKRZ). For coupled runs, the atmosphere and ocean processes can run on the same nodes (job:run:share_nodes = true), or they can be distributed as two blocks one after the other.
* `levante_gpu`: run on the GPU partition of Levante (DKRZ). IO processes are distributed equally among the CPUs of the GPU nodes.
* `dolpung_gpu`: run on the GPUs of the dolpung partition of Levante (DKR). IO processes are distributed equally among the CPUs of the nodes.
* `dolpung_hybrid`: run on the CPUs and the GPUs of the dolpung partition by using two binaries. IO processes run on the CPUs of the last node. For more than two nodes, the 1st node is dedicated to atmosphere processes only. Ocean processes are distributed among the remaining CPUs.

### Supported grids

* `mkexp` supports R2B4, R2B6, R2B8, R02B9, and R2B10 grids (check for `Change the resolution` in mkexp runscripts). See [here](ref_buildrun_grid_files) more information about grids.


## Frequently changed settings

### Input files

Input files required for an experiment are defined in `mkexp`'s `[files]` section as described in "Defining input files for an experiment" in `utils/mkexp/doc/mkexp.pdf`.
For ICON, we add several conventions to simplify common changes.

#### Model components use their own subsections

Files used by a specific ICON model component have their own `[[subsection]]` with a defined base directory; each base directory is controlled by a `GLOBAL_VARIABLE`:

* Atmosphere: `[[atmosphere]]`; `ATMO_INPUT_DIR`
* Land: `[[land]]`; `LAND_INPUT_DIR`
* Ocean: `[[ocean]]`; `OCEAN_INPUT_DIR`
* Ocean Biogeochemistry (HAMOCC): `[[obgc]]`; _set using `.base_dir` (see below)_

There may be subsections and subsubsections below this level to group files with common requirements.

#### Input files go into the working directory

Usually, variables defined in these sections define input files to be made available in the experiment's working directory.
Their value specifies where to find the corresponding data files and their original names.

The default is to create a symbolic link to the data file given as value, where the variable name is used as the link's new name.
Instead, the special variable `.method` can be set to define an arbitrary command or utility.
This is called with _value_ and _variable name_ as arguments, eg.
`.method = cp` and `aerosol.nc = aerosol_1950.nc` will run `cp aerosol_1950.nc aerosol.nc`.

`.method` is valid for all files of the current section. It's also used for its subsections, unless redefined there.

#### Sections may map to directory structure

The variable `.base_dir` specifies an absolute directory where to look for the _value_ files.
Like `.method`, this hold for all files of the current section and below.

Instead of redefining `.base_dir`, a subsection may define `.sub_dir` as a relative directory.
This is appended to the currently valid `.base_dir` from an ancestor section
such that files in this subsection are now expected in their respective subdirectory.

A section's settings will be passed to any subsection and in turn to its subsections,
taking into account amendments or redefinitions on the way.

#### Restart files from spin-up and own experiment

Restart file definitions use the subsection `[[[restart]]]` within the atmosphere and ocean section.
They define `restart_{atm,oce}_DOM01.nc` for single file restarts or `multifile_restart_{atm,oce}.mfr` for multi-file restarts.

For restart files from spin-up runs, their base dir is set to `$PARENT_DIR`.
During the experiment, its own restart files are based in the model's working directory.
