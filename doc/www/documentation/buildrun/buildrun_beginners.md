```{eval-rst}
:orphan:
```

(ref_buildrun_beginners)=
# ICON Tutorial for Beginners

This tutorial will guide you from cloning the ICON repository to successfully run an ICON simulation at the DKRZ Levante supercomputer. If you experience any issue or have some questions, check our [FAQ section](ref_buildrun_faq).

**Before starting:** make sure that git is available in your environment, and that you have an active python environment with `six` and `jinja2` installed.

## Cloning the Repository

Run the following command to clone the ICON repository, as well as to initialize and update the necessary submodules:

```
git clone --recurse-submodules git@gitlab.dkrz.de:icon/icon-mpim.git
```

The location where the repository is cloned becomes `<src-dir>`.


## Configuring and Building

In this tutorial, you will build ICON using an *in-source build*, i.e. building ICON in the main repository directory. This requires navigating to `<src-dir>`, and then calling a configuration script. To configure ICON to build with the `intel` compiler at `levante` supercomputer from DKRZ, do:

```
cd <src-dir>
./config/dkrz/levante.intel
```

Note the structure `config/<institution>/<machine>.<compiler>` of the config file system.

Alternatively, to use GPU, you can use the `nvhpc` compiler (note that `<compiler>` may include a `cpu.` or `gpu.` specification):

```
./config/dkrz/levante.gpu.nvhpc
```

The building step is done by running `make` command with an optional argument specifying the number of jobs to run simultaneously. Usually, 8 is a good choice:

```
make -j8
```

## Running

Here you will run an [AMIP experiment](ref_buildrun_AMIP), one of ICON [supported scientific configurations](ref_buildrun_supportedconf).

First, copy the example script of the desired configuration (`amip-aes.config` for CPUs or `amip-aes_gpu.config` for GPUs) to the run folder of the build directory. One good practice here is to add info about the compiler and commit id of icon version used for the build, as it will help a lot in revisiting the experiments in the future (run `git rev-parse --short HEAD` or `git log -1 --oneline` to get the commit id). So, at the time of writing this tutorial, the appropriate copying command would be:

```
cd run/
cp examples/amip-aes.config amip-aes-tutorial-intel-85bdcf4785.config
```

Open `amip-aes-tutorial.config` and adjust the experiment id:

```
EXP_ID = amip-aes-tutorial-intel-85bdcf4785
```

Now, choose the type of grid (`R<n>B<k>`) to use. To begin with, we use a small grid, `R2B4` (note that the type and size of grid used determines the number of nodes required for the run). Other grids are specified at the end of the file itself.

```
EXP_TYPE = amip-aes-R2B4
```

Add a new `ACCOUNT = ...` line for the account of the slurm project you are member of. `ACCOUNT` is a basic setting, so make sure this line is _above_ any `[...]` section. Accounts are composed of two letters and four numbers, e.g.:

```
ACCOUNT = mh0000
```
Adjust the number of nodes (R2B4 grids can run with a single node):

```
nodes = 1
```

Adjust the initial and last day of simulation as well as the interval for the restart file (here you will run three hours and print one restart file at the end of the simulation). Note that the timestamp format follows the [ISO 8601 for date- and time-related data](https://en.wikipedia.org/wiki/ISO_8601). Again, make sure that the following three lines are above any `[...]` section:

```
INITIAL_DATE = 2020-01-01T00:00:00
FINAL_DATE = 2020-01-01T03:00:00
INTERVAL = PT3H
```

Now, generate the run scripts using [mkexp](ref_tools_mkexp):
```
../utils/mkexp/mkexp amip-aes-tutorial-intel-85bdcf4785.config
```

This generates a runscript `amip-aes-tutorial-intel-85bdcf4785.run_start` at the `Script directory` (see path in output). Navigate there and finally submit the experiment to Levante:

```
cd ../experiments/amip-aes-tutorial-intel-85bdcf4785/scripts/
sbatch amip-aes-tutorial-intel-85bdcf4785.run_start
```

You can find runtime log of the simulation at the `amip-aes-tutorial-intel-85bdcf4785.run.log` file, and the output at the `Work directory` directory, in `run_<init-date>-<final-date>`.
