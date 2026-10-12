# ------------------------------------------
# Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

################################################################
# ICON ecflow suite
#   conversion of SMS icon_multi.def suite
# Martin Koehler, DWD, 2022-11-01
################################################################


import os

# from ecflow import Defs,Suite,Family,Task,Edit,Trigger,Event,Complete
from ecflow import *

# ---USER SETUP-------------------------------------------------

EXPNUM = "417"  # ---experiment number
EXPNUM2 = "416"  # ---experiment for comparison

FIRST_DATE = 20160801
LAST_DATE = 20160801
INIHOUR = "00"
MONTHS = "201608"  # ---months to be run e.g. 201201_201207
MON_IFS_ICON = "ifs"

RES = "R03B07"  # ---ICON resolution dynamics
RESRAD = "R03B06"  # ---ICON resolution radiation

# number of time steps (10x86400s/360s)
NDAYS = 45
DTIME = 120
NSTEPS = int(NDAYS * 86400 / DTIME)

# ntasks for model   R2B6: 88    R3B7: 620   408   [s/day]
NTASKS = 2048  # 1024        1024  2048
CPUSPERTASK = 16  # 4           16    16

# simulation modes: 1:oper, 2:monthly 10day forecasts, 3:climate (1year), 4:B-matrix
SIMMODE = 2

# no. of ensemble members (dummies here)
NENS = 0  # default for no ensemble run
NMEM = 0  # default for no ensemble run

# parallelisation of tasks: forecast, eom_prepare_2, eom_metview
npar_forecast = 1  # [1-10] number of parallel forecast jobs
npar_prepare2 = 1  # [1,8]  number of parallel data prepare jobs
npar_metview = 1  # [1-10] number of parallel metview jobs

# HPC environment:
perm = os.getenv("PERM")
hpcperm = os.getenv("HPCPERM")
user = os.getenv("USER")
home = os.getenv("HOME")

SCHOST = "hpc"
SCTEMP = "/ec/res4/scratch/" + user + "/icon"
SCPERM = perm

SCBASEDIR = perm + "/ecflow/icon-nwp/schedulers/ecmwf"
SCCODEDIR = perm + "/icon-nwp"
INPUTDIR = "/ec/res4/hpcperm/dei2/icon-input"


# ---SUITE DEFINITION-------------------------------------------

print("Creating suite definition")


def hpc_setup():
    return Edit(
        ECF_JOB_CMD="troika submit -o %ECF_JOBOUT% %SCHOST% %ECF_JOB%",
        ECF_KILL_CMD="troika kill                   %SCHOST% %ECF_JOB%",
        ECF_STATUS_CMD="troika monitor                %SCHOST% %ECF_JOB%",
        QUEUE="nf",
    )


def loop_forecast(i):
    fc_name = "fc{}".format(i)
    return Family(
        fc_name,
        RepeatDate("YMD", FIRST_DATE + i - 1, LAST_DATE, npar_forecast),
        Trigger("../init:YMD gt " + fc_name + ":YMD or ../init == complete"),
        Task(
            "testday",
            hpc_setup(),
            Event("time_eom"),
            Event("time_domonth"),
            Meter("time_ifs_icon", 1, 2),
        ),
        Task(
            "model",
            hpc_setup(),
            Trigger("testday==complete"),
            Complete("testday==complete and not testday:time_domonth"),
            Edit(QUEUE="np", NTASKS=NTASKS, CPUSPERTASK=CPUSPERTASK),
        ),
        Task(
            "check_progress",
            hpc_setup(),
            Trigger("model == active"),
            Complete(
                "model == complete or testday==complete and not testday:time_domonth"
            ),
            Edit(FC_NAME=fc_name),
            Meter("timesteps", 0, NSTEPS),
        ),
    )


def loop_post(i):
    fc_name = "fc{}".format(i)
    post_name = "post{}".format(i)
    return Family(
        post_name,
        RepeatDate("YMD", FIRST_DATE + i - 1, LAST_DATE, npar_forecast),
        Trigger(
            "../forecast_all/"
            + fc_name
            + ":YMD gt "
            + post_name
            + ":YMD or ../forecast_all == complete"
        ),
        Task(
            "testday",
            hpc_setup(),
            Event("time_eom"),
            Event("time_domonth"),
            Meter("time_ifs_icon", 1, 2),
        ),
        Task(
            "post_prepare",
            hpc_setup(),
            Trigger("testday==complete"),
            Complete("testday==complete and not testday:time_domonth"),
            Edit(FC_NAME=fc_name),
            Trigger("testday == complete"),
            Event("meteogram_data"),
        ),
    )


def trigger_eom(i):
    many_post_triggers = "post/post1:YMD gt endofmonth:YMD"
    for n in range(2, i + 1):
        post_name = "post{}".format(n)
        many_post_triggers = (
            many_post_triggers
            + " and post/"
            + post_name
            + ":YMD gt endofmonth:YMD"
        )
    return Trigger(many_post_triggers + " or post == complete")


def loop_prepare2(i):
    return Family(
        "prep2_{}".format(i),
        Task(
            "eom_prepare_2",
            hpc_setup(),
            Edit(prepall=npar_prepare2, prepnum=i),
            Meter("catvariables", 0, 100),
            Event("catdata"),
        ),
    )


def loop_metview(i):
    return Family(
        "met{}".format(i),
        Task("eom_metview", hpc_setup(), Edit(metnum=i, metall=npar_metview)),
    )


# --------------------------------------------------------------

suite0 = Suite(
    "fconly",
    Edit(
        ECF_INCLUDE=SCBASEDIR + "/include",
        ECF_FILES=SCBASEDIR + "/tasks",
        #       ECF_HOME   = home +'/ecflow_server',  # automatically set
        #       ECF_OUT    = home +'/ecflow_server',  # not set
        ECF_TRIES=1,
        SCHOST=SCHOST,
        SCTEMP=SCTEMP,
        SCPERM=SCPERM,
        SIMMODE=SIMMODE,
    ),
    # EXP --------------------------------------------------
    Family(
        user + "_" + EXPNUM,
        Edit(
            USER=user,
            SCCODEDIR=SCCODEDIR,
            SCBASEDIR=SCBASEDIR,
            INPUTDIR=INPUTDIR,
            EXPNUM=EXPNUM,
            EXPNUM2=EXPNUM2,
            INIHOUR=INIHOUR,
            MONTHS=MONTHS,
            MON_IFS_ICON=MON_IFS_ICON,
            RES=RES,
            RESRAD=RESRAD,
            NSTEPS=NSTEPS,
            NENS=NENS,
            NMEM=NMEM,
        ),
        # ONCE ---------------------------------------------
        Task("once", hpc_setup()),
        # INITIALIZATION -----------------------------------
        Family(
            "init",
            RepeatDate("YMD", FIRST_DATE, LAST_DATE, 1),
            Trigger("once == complete"),
            Task(
                "testday",
                hpc_setup(),
                Event("time_eom"),
                Event("time_domonth"),
                Meter("time_ifs_icon", 1, 2),
            ),
            Task(
                "setup",
                hpc_setup(),
                Trigger("testday:time_domonth"),
                Complete("testday==complete and not testday:time_domonth"),
                Event("enable_build"),
            ),
            Family(
                "build",
                Trigger("../once == complete and setup == complete"),
                Complete("testday==complete and not testday:time_domonth"),
                Task(
                    "init_cp_binary",
                    hpc_setup(),
                    Complete("../../init:YMD gt " + str(FIRST_DATE)),
                ),
            ),
            Task(
                "get_data",
                hpc_setup(),
                Trigger("../once == complete and setup == complete"),
                Complete("testday==complete and not testday:time_domonth"),
                Edit(QUEUE="nf", NTASKS=1, NTHREADS=8),  # MPI  # openMP),
            ),
            Task(
                "init_data",
                hpc_setup(),
                Trigger("get_data == complete and build == complete"),
                Complete(
                    "testday:time_ifs_icon == 2 or testday==complete and not testday:time_domonth"
                ),
                Edit(QUEUE="np", NTASKS=1, NTHREADS=8),  # MPI  # openMP
            ),
        ),
        # FORECAST -----------------------------------------
        Family(
            "forecast_all",
            Edit(NDAYS=NDAYS, DTIME=DTIME),
            [loop_forecast(i) for i in range(1, npar_forecast + 1)],
        ),
        # POST-PROCESS -------------------------------------
        Family("post", [loop_post(i) for i in range(1, npar_forecast + 1)]),
        # END OF MONTH -------------------------------------
        Family(
            "endofmonth",
            RepeatDate("YMD", FIRST_DATE, LAST_DATE, 1),
            trigger_eom(npar_forecast),
            Task(
                "testday",
                hpc_setup(),
                Event("time_eom"),
                Event("time_domonth"),
                Meter("time_ifs_icon", 1, 2),
            ),
            Family(
                "eom_prepare_all",
                Trigger("testday == complete"),
                Complete(
                    "testday==complete and (not testday:time_eom or not testday:time_domonth)"
                ),
                Task(
                    "eom_prepare_1",
                    hpc_setup(),
                    Meter("ecget", 0, 31),
                    Meter("splitdata", 0, 31),
                    Meter("deldata", 0, 31),
                    Event("get_pl_data"),
                    Event("get_hl_data"),
                    Event("get_ml_data"),
                ),
                Family(
                    "eom_prepare_2_all",
                    Trigger("eom_prepare_1 == complete"),
                    [loop_prepare2(i) for i in range(1, npar_prepare2 + 1)],
                ),
                Task(
                    "eom_prepare_3",
                    hpc_setup(),
                    Trigger("eom_prepare_2_all == complete"),
                    Event("cpdata"),
                ),
            ),
            Task(
                "eom_getmars",
                hpc_setup(),
                Trigger("testday == complete"),
                Complete(
                    "testday==complete and (not testday:time_eom or not testday:time_domonth)"
                ),
            ),
            Family(
                "eom_metview_all",
                Trigger(
                    "eom_prepare_all == complete and eom_getmars == complete"
                ),
                Complete(
                    "testday==complete and (not testday:time_eom or not testday:time_domonth)"
                ),
                [loop_metview(i) for i in range(1, npar_metview + 1)],
            ),
            Task(
                "eom_archive",
                hpc_setup(),
                Trigger("eom_metview_all == complete"),
                Complete(
                    "testday==complete and (not testday:time_eom or not testday:time_domonth)"
                ),
            ),
        ),
    ),
)


# ---Checking and DEF file generation---------------------------

defs = Defs(suite0)
print(defs)

print("Checking job creation: .ecf -> .job0")
print(defs.check_job_creation())

print("Saving definition to file 'icon_fc.def'")
defs.save_as_defs("icon_fc.def")
