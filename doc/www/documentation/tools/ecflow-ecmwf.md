```{eval-rst}
:orphan:
```

(ref_tools_ec_ecflow)=
# ICON at ECMWF: ecflow suite

Below you find a quick and dirty itemized get-started page. This setup
does all the steps for you (initial data, ICON running, plotting,
sending plots to DWD) within the ecflow scheduler. Initial plots for the first 2 forecasts are available after 45min and all 2 months are
fully processed and deliverd after 3 hours.

For more details about ecflow look
[here](https://confluence.ecmwf.int/display/UDOC/HPC2020%3A+Using+ecFlow).

## Get ECMWF account for DWD ICON users

An application has to be sent to the DWD computer representative for ECMWF in the branch TI.

## Apply for ECFLOW server

Write a ticket on the ECMWF [confluence
support](https://confluence.ecmwf.int/site/support) web page applying
for a personal virtual ecflow server.

## Compile ICON

The default place to search for an executable is
`\$PERM/icon/icon-nwp`. See below for changing that.
Remember to configure with option `--disable-mixed-precision` before compiling.
The required tendency output on lat-lon grids is not support mixed-precision,
slowing the simulation by 20%.

## Branch for ecflow scripts

You have the option to have a separate branch or directory for all the ecflow scripts, for example, `icon-nwp-ecflow` (specify in icon_multi.py, see below). Open that branch.

    cd schedulers/ecmwf

There are 3 important directories:

-   **suites**: suite definition python files.
    -   icon_multi.py defines a month long set of 10d forecasts (called
        fcmon)
    -   play_suite.py is used once to set up ecflow suite (e.g. fcmon or
        fconly)
-   **cases**: namelist files, use links to the files below if you have many
    namelist files
    -   case_setup_icon: used for ICON initialisation
    -   case_setup_ifs: used for IFS initialisation
-   **tasks**: ksh scripts, each corresponds to one icon in ecflow
    -   model.ecf: uses namelist and runs ICON
    -   init_cp_binary: defines executable location to be used

## Setup ecflow suite

The following you can add to your .profile file.

    module load python3
    module load ecflow/new

## Namelists

The files `case_setup_icon` and `case_setup_ifs` in the directory
`cases` are the default namelist file names used for ICON and IFS
initialised setups. You have to copy or link an available case_setup
file to those locations.

Attention: if you later want to change the case_setup namelist files you
need to rerun the icon `once` in the ecflow_ui (see below).

## Suite definition creation

Edit the python script for suite definition **icon_multi.py**. Select
the experiment IDs **EXPNUM** and **EXPNUM2** (e.g. 001 and 001).
Optional you can use the input files, grids etc, from user dei2
(Martin):

    INPUTDIR  = '/perm/dei2/icon/icon-input'

Create the suite definition file **icon_multi.def** with

    python3 icon_multi.py

Possible errors will be detected.

## Load suite to ecflow server

Edit in directory suites the python scripts **play_suite.py** and
**play_exp.py**. Change user to your ECMWF user ID and experiment ID.

To create a **new suite** you need to run

    python3 play_suite.py

Once the suite (for example Damon) is already created and you want to
start a **new experiment** you use

    python3 play_exp.py

This can also be used to reload an experiment (if there are NO green
buttons).

## ecflow user interface

You can check if the ecflow suite has been correctly loaded in

    ecflow_ui

Here you can also follow the progress of the experiment, suspend, kill
and restart tasks.

You need to add a server:

-   select menu **Servers/Manage servers**
-   **Add server**
    -   Name: "dei2" or "Martin" or anything
    -   Host: ecflow-gen-????-001 (as given to you by ECMWF user
        support)
    -   Port: 3141
    -   Note: you can also add others servers

## Output location

In directory

    ${SCRATCH}/icon/fcmon

there are three subdirectories with init data, model output and plots:

    input
    output
    post

## Reload binary

After you recompile your code you need to requeue the ecflow icon
icon_cp_binary to make it active.

## Transfer plots to DWD

The transfer is done after plot generation in task `eom_archive`. You
need to setup entrans to allow for automatic data transfers between
ECMWF and DWD. Go to
[website](https://oflkd013.dwd.de/ecmwf/gateway/ECtrans) enter ECMWF UID
and token pw. Und "Ectrans" you can set up an association with the
following details:

-   host name: rcl.dwd.de
-   directory: /hpc/uwork/mkoehler/icon-plots (your UID, make that
    directory)
-   genericSftp
-   login: DWD UID (e.g. mkoehler)
-   password: password on rcl (attention change PW here when changing
    rcl PW)
-   association name: ECtoRCL_dei2 (your ECMWF UID)
-   sftp.execCmd = "/hpc/uhome/mkoehler/bin/ectrans_2oflxs464.s
    \$filename". (here you need to edit the file entrains_2oflxs464.s to
    change UIDs.

## Display plots at DWD in plot catalog

The plots from the ECMWF ecflow experiments are automatically displayed
in the DWD
[plot-catalog](http://oflxs464.dwd.de/~mkoehler/plot-catalog/index.php?cat=cat.koehler.icon.ec.metview_month.xml).

You can setup your own plot-catalog for viewing your experiments by
copying Martin's XML file:

    /fe1-daten/mkoehler/public_html/plot-catalog/cat.koehler.icon.ec.metview_month.xml

to your analogue directory, edit it and view it as follows: (edit only
USER!)

    http://oflxs464.dwd.de/~mkoehler/plot-catalog/index.php?cat=http://oflxs464.dwd.de/~USER/plot-catalog/cat.USER.icon.ec.metview_month.xml

You need to make two links on oflxs646 to make your xml files and data
files available for the web:

    public_html -> /fe1-daten/mkoehler/public_html
    /uwork1/USER/icon/plots -> public/html/plots_wrk

Also you need to copy the ectrans scripts:

workstation:

    cp /fe1-daten/mkoehler/bin/ectrans_operate.s /fe1-daten/$USER/bin/ectrans_operate.s
    mkdir /fe1-daten/$USER/temp

rcl:

    cp /hpc/uhome/mkoehler/bin/ectrans_2oflxs464.s /hpc/uhome/$USER/bin/ectrans_2oflxs464.s
    mkdir /hpc/uhome/$USER/temp

## Trouble shooting

-   The ecflow server requires password-free ssh between atos hosts!
