# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

"""
Purpose
    Read and plot SCM and LES from NetCDF file with Python3
Details
    Variable-Vertical plot
Help
    module load unsupported
    module load python-extras/3.8.12
Author
    Martin Koehler, DWD
    Run as: python3 plot-scm+les-var-z.py RICO theta_v 15  (varname hour)
Revision history
    20200701 -- Initial version
"""

import datetime
import glob  # wildcard file completion
import sys

import matplotlib.pyplot as plt
import mpl_toolkits as mp
import numpy as np
from netCDF4 import Dataset, date2num, num2date

sys.path.append("/hpc/uhome/mkoehler/bin/python")
from autoscale import autoscale

# -------------------------------------------------------------------
# arguments

print("")
print("Number of arguments:", len(sys.argv), "arguments.")
print("Argument List:      ", str(sys.argv))
case = sys.argv[1]
varname = sys.argv[2]
if len(sys.argv) == 4:
    hour = np.array([float(sys.argv[3])])
else:
    hour = [float(sys.argv[3]), float(sys.argv[4])]
print("hour: ", hour)

# -------------------------------------------------------------------
# setup

mean = 0  # 0: single time level plots, 1: time mean between 2 hours

# hour     = np.array([15,15+24,15+48])
# varname  = 'temp'               # temp, u, tot_qv_dia
# npoints  = (0,                  # one point (0 is 1st point)
# npoints  = range(32)            # plot all 32 points
# npoints  = (-1,)                # mean of 32 points
modname = ""
# modname  = '_edmf'

# SCM output:
nc_files = [
    "/hpc/uwork/mkoehler/run-icon/scm/SCM_"
    + case
    + "_dephy"
    + modname
    + "/out_SCM_"
    + case
    + "_dephy"
    + modname
    + "_ML_*Z_mean.nc",
    "/hpc/uwork/mkoehler/run-icon/scm/LES_100x100_"
    + case
    + "_dephy/out_ml_LES_100x100_"
    + case
    + "_dephy_ML_*Z_mean.nc",
]
#           '/hpc/uwork/mkoehler/run-icon/scm/LES_50x50_'+case+'_dephy/out_ml_LES_50x50_'+case+'_dephy_ML_20041216T000000Z_mean.nc']
nc_zfiles = [
    "/hpc/uwork/mkoehler/run-icon/scm/SCM_"
    + case
    + "_dephy/out_SCM_"
    + case
    + "_dephy_ML_*Z_mean.nc",
    "/hpc/uwork/mkoehler/run-icon/scm/LES_100x100_"
    + case
    + "_dephy/out_z_LES_100x100_"
    + case
    + "_dephy_ML_*Z_mean.nc",
]


nfiles = len(nc_files)
labels = ["SCM", "LES 100x100"]  # , 'LES 50x50']


# -------------------------------------------------------------------
# read files


# -------------------------------------------------------------------
# plot loop over models

title = "ICON SCM and LES forecasts      " + case
subtitle = "forecast time: " + str(hour[0]) + "h"
fig, ax = plt.subplots(figsize=(5, 7))
fig.subplots_adjust(bottom=0.15, left=0.2)
ax.margins(y=0)

nn = 0
for file in nc_files:

    nc_file = glob.glob(file)  # wildcard completion
    nc_zfile = glob.glob(nc_zfiles[nn])
    print("reading file: ", nc_file[0])
    print("height file:  ", nc_zfile[0])
    nc_fid = Dataset(nc_file[0], "r")
    nc_fidz = Dataset(nc_zfile[0], "r")
    vardata = nc_fid.variables[varname][:]
    vardata = np.squeeze(vardata)
    dates = num2date(
        nc_fid.variables["time"][:], nc_fid.variables["time"].units
    )

    print("  shape data and date: ", vardata.shape, dates.shape)

    z_mc = nc_fidz.variables["z_mc"][:]

    # time conversion

    num_dates = [
        int(d.strftime("%Y%m%d%H")) for d in dates
    ]  # YYYYMMDDHH e.g. 2019081516
    hours = [
        (d.day - dates[0].day) * 24
        + (d.hour - dates[0].hour)
        + (d.minute / 60)
        + (d.second / 3600)
        for d in dates
    ]

    # find steps for hours_plot
    steps = np.empty((len(hour),), dtype=int)
    nnn = 0
    for h_plot in hour:
        steps[nnn] = h_plot / (hours[1] - hours[0])
        nnn = nnn + 1
    print("SCM hours and steps:       ", hours[steps[0]], steps)

    # plot

    if nn == 0:
        xtitle = "%s [%s]" % (
            nc_fid.variables[varname].standard_name,
            nc_fid.variables[varname].units,
        )
        ytitle = "height [m]"
        ytitle_right = "model levels"
        ax.set_title(title + "\n" + subtitle)
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)
        ax.set_ylim(0, 5000)

    # for ns in steps:
    #   for np in npoints:
    #     if np >= 0:
    #       datax = vardata[ns,:,np-1]
    #       datay = z_mc [ns,:,np-1]
    #     else:
    #       datax = vardata[ns,:,:].mean(axis=1)
    #       datay = z_mc   [ns,:,:].mean(axis=1)

    datay = np.squeeze(z_mc)
    if mean == 0:
        for ns in steps:
            datax = vardata[ns]
            print("plot data: ", datax.shape, datay.shape)
            ax.plot(
                datax, datay, label=labels[nn]
            )  #   str(hours[ns])+" h")   # (" + str(num_dates[ns]) + ")")
    else:
        print("hours and steps for mean", hour[0], hour[1], steps[0], steps[1])
        datax = vardata[steps[0] : steps[1]].mean(axis=0)
        print("plot data: ", datax.shape, datay.shape)
        ax.plot(datax, datay, label=labels[nn])

    # right axis with model level locations

    if nn == 0:
        ax2 = ax.twinx()  # second axes
        ax2.set_yticks(datay)
        ax2.tick_params(axis="y", labelright=False)
        ax2.set_ylabel(ytitle_right).set_fontsize(14)
        ax2.set_ylim(ax.get_ylim())  # use y-limits as in ax

    # label for each line

    ax.legend()

    autoscale(
        ax, "x", margin=0.1
    )  # automatic scaling of x-axis given fixed y_lim

    nn = nn + 1


dirplots = "plots/" + case + modname + "/"
if mean == 0:
    hourstxt = "+" + str(hour[0])
else:
    hourstxt = "+" + str(hour[0]) + "to" + str(hour[1])

fig.savefig(
    dirplots
    + "var-z_"
    + varname
    + "_"
    + case
    + modname
    + "_"
    + str(num_dates[0])
    + hourstxt
    + ".png"
)  # png
fig.savefig(
    dirplots
    + "var-z_"
    + varname
    + "_"
    + case
    + modname
    + "_"
    + str(num_dates[0])
    + hourstxt
    + ".pdf"
)  # pdf
# plt.show()                                # to screen
