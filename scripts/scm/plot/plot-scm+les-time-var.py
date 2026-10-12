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
    Read and plot SCM & LES from NetCDF file with Python3
Details
    Time-Variable plot
    Run as: python3 plot-scm+les-time-var.py RICO shfl_s
Author
    Martin Koehler, DWD
Revision history
    202201 -- Initial version
"""

import datetime
import getopt
import glob  # wildcard file completion
import math
import os
import sys  # file handling and operators
from datetime import datetime

import matplotlib.pyplot as plt
import mpl_toolkits as mp
import numpy as np
from netCDF4 import Dataset, date2num, num2date
from scipy.signal import savgol_filter

# -------------------------------------------------------------------
# arguments

print("")
print("Number of arguments:", len(sys.argv), "arguments.")
print("Argument List:      ", str(sys.argv))
case = sys.argv[1]
varname = sys.argv[2]

# -------------------------------------------------------------------
# setup

# varname      = 'shfl_s'      # temp, t_2m, qv_2m, u_10m, v_10m, shfl_s, lhfl_s, clct, tqc_dia
level = 1  #                   # 1000hPa, 90 (model level), 1 (e.g. t_2m)
# npoints      = (0,)          # one point (0 is 1st point)
# npoints      = range(32)     # plot all 32 points
# npoints      = (-1,)         # mean of 32 points
npoints = (-2,)  #             # mean input files
n2d3d = 2  #                   # 2/3: 2D or 3D variable
# case         = 'RICO'
# modname      = ''
# modname      = '_edmf'
modname = "_test"

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
    + "_dephy/out_sfc_LES_100x100_"
    + case
    + "_dephy_ML_*Z_mean.nc",
]
# '/hpc/uwork/mkoehler/run-icon/scm/LES_50x50_'+case+'_dephy/out_sfc_LES_50x50_'+case+'_dephy_ML_20041216T000000Z_mean.nc']
nfiles = len(nc_files)
labels = ["SCM", "LES 100x100"]  # , 'LES 50x50']


# -------------------------------------------------------------------
# read files

# find longest time-series
ntime = np.arange(nfiles)
nn = 0
for file in nc_files:
    file = glob.glob(file)
    print("length of file: ", file)
    nc_fid = Dataset(file[0], "r")
    data = nc_fid.variables[varname][:]
    ntime[nn] = data.shape[0]
    nn = nn + 1
maxtime = max(ntime)
print("ntimes: ", ntime)

# reading data and dates
vardata = np.reshape(
    np.arange(0, nfiles * maxtime, dtype=float), (nfiles, maxtime)
)
dates = np.reshape(
    np.arange(0, nfiles * maxtime, dtype=datetime), (nfiles, maxtime)
)
nn = 0
for file in nc_files:
    file = glob.glob(file)
    print("reading file: ", file)
    nc_fid = Dataset(file[0], "r")
    data = nc_fid.variables[varname][:]
    date = num2date(nc_fid.variables["time"][:], nc_fid.variables["time"].units)
    vardata[nn, 0 : ntime[nn]] = np.squeeze(data)
    dates[nn, 0 : ntime[nn]] = date
    nn = nn + 1

# -------------------------------------------------------------------
# variable specifics

# Variable names:
if varname == "t_2m":
    varname_axis = "T 2m"
elif varname == "shfl_s":
    varname_axis = "Sensible Heat Flux"
elif varname == "lhfl_s":
    varname_axis = "Latent Heat Flux"
elif varname == "sob_s":
    varname_axis = "Shortwave Surface Flux"
elif varname == "sob_t":
    varname_axis = "Shortwave Top Flux"
elif varname == "thb_s":
    varname_axis = "Thermal Surface Flux"
elif varname == "thb_t":
    varname_axis = "Thermal Top Flux"
elif varname == "clct":
    varname_axis = "Total Cloud Cover"
else:
    print("Variable name not available yet.")
    varname_axis = varname

if (
    varname == "lhfl_s"
    or varname == "shfl_s"
    or varname == "sob_t"
    or varname == "sob_s"
    or varname == "thb_t"
    or varname == "thb_s"
    or varname == "clct"
):
    vardata = vardata.reshape(vardata.shape[0], 1, vardata.shape[1])
if (
    varname == "lhfl_s"
    or varname == "shfl_s"
    or varname == "sob_t"
    or varname == "sob_s"
    or varname == "thb_t"
    or varname == "thb_s"
):
    vardata[0, 0, :] = vardata[
        1, 0, :
    ]  # first flux point is diagnostic and wrong


print("variable shapes", vardata.shape)

# flatten SCM model
if n2d3d == 3:
    vardata = vardata[:, level - 1, :]
# vardata = np.squeeze(vardata)


# -------------------------------------------------------------------
# time conversion

num_dates = [
    int(d.strftime("%Y%m%d%H")) for d in dates[0, 0 : ntime[0]]
]  # YYYYMMDDHH e.g. 2019081516
hours = [
    (d.day - dates[0, 0].day) * 24
    + (d.hour - dates[0, 0].hour)
    + (d.minute / 60)
    + (d.second / 3600)
    for d in dates[0, 0 : ntime[0]]
]

# deaccumulate
# if varname == 'lhfl_s' or varname == 'shfl_s' or varname == 'sob_t' or varname == 'sob_s' or \
#   varname == 'thb_t'  or varname == 'thb_s':
#  vardata3[1:,0,0] = ( vardata3[1:,0,0] - vardata3[0:vardata3.shape[0]-1,0,0] ) / hours3[1] / 3600.0
#  vardata3[0,0,0]  = vardata3[1,0,0]

# smoothing:
#   savitztky/golay filter: least squares to regress a small window of your data onto a polynomial
#   savgol_filter(x, window_length, polyorder, ...)

# tsmooth = 1.5         # smoothing in [h]
# print( 'smoothing with window width:', math.ceil(tsmooth/(hours3[1]-hours3[0])/2)*2+1 )
#
# vardata   = savgol_filter(vardata , math.ceil(tsmooth/(hours [1]-hours [0])/2)*2+1, 1, axis=0)
# vardata3  = savgol_filter(vardata3, math.ceil(tsmooth/(hours3[1]-hours3[0])/2)*2+1, 1, axis=0)

# -- unit formatting:  convert [W m-2] to [$W m^{-2}$] etc
unit = nc_fid.variables[varname].units
for r in (("-", "^{-"), ("1", "1}"), ("2", "2}"), ("**", "")):
    unit = unit.replace(*r)
if unit != "%":
    unit = "[$" + unit + "$]"

# -------------------------------------------------------------------
# plot

title = "ICON SCM and LES forecasts      " + case
xtitle = "Simulation time since %s [h]" % (num_dates[0])
# ytitle = "%s (%s)" % (nc_fid.variables[varname].standard_name,\
#                      nc_fid.variables[varname].units)
# ytitle = "%s [%s]" % (varname_axis,nc_fid.variables[varname].units)
ytitle = varname_axis + "  " + unit

dt_ticks = 6.0
ticks = np.arange(0.0, max(hours) + 0.1, dt_ticks)

fig, ax = plt.subplots(figsize=(7, 5))
fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(left=0.15)

ax.set_title(title)
# ax.set_title('(d)',fontsize=16, fontweight='bold',loc='right')
ax.set_xlabel(xtitle, fontsize=14)
ax.set_ylabel(ytitle, fontsize=14)
ax.set_xticks(ticks)
# ax.set_xlim(0,72)
# ax.set_ylim(281,297)       # t_2m
# ax.set_ylim(0.004,0.008)   # qv_2m
# ax.set_ylim(-4,1)          # v_10m
# ax.set_ylim(-300,10)       # lhfl_s
# ax.set_ylim(-400,50)       # shfl_s

plt.setp(ax.get_xticklabels(), fontsize=14)  # xtick labels fontsize
plt.setp(ax.get_yticklabels(), fontsize=14)  # ytick labels fontsize
plt.rc("legend", fontsize=13)  #             # legend fontsize

# plot SCM line

# for npt in npoints:
#   if npt >= 0:
#    #data = vardata[:,level-1,npt-1]
#     data = vardata[:,npt-1]
#   elif npt == -1:
#    #data = vardata[:,level-1,:].mean(axis=1)
#     data = vardata[:,:].mean(axis=1)
#   else:
#     print(vardata.shape)
#     data = np.squeeze(vardata)

vardata = np.squeeze(vardata)

# loop over models
nn = 0
for nt in ntime:
    hours = [
        (d.day - dates[0, 0].day) * 24
        + (d.hour - dates[0, 0].hour)
        + (d.minute / 60)
        + (d.second / 3600)
        for d in dates[nn, 0:nt]
    ]
    # optional smoothing
    plotdata = vardata[nn, 0:nt]
    tsmooth = 1.5  # smoothing in [h]
    plotdata = savgol_filter(
        plotdata,
        math.ceil(tsmooth / (hours[1] - hours[0]) / 2) * 2 + 1,
        1,
        axis=0,
    )

    ax.plot(hours[0:nt], plotdata, linewidth=2.0, label=labels[nn])
    nn = nn + 1

# label for each line

ax.legend()

dirplots = "plots/" + case + modname + "/"
fig.savefig(
    dirplots
    + "time-var_"
    + varname
    + "_"
    + case
    + modname
    + "_"
    + str(num_dates[0])
    + ".png"
)  # png
fig.savefig(
    dirplots
    + "time-var_"
    + varname
    + "_"
    + case
    + modname
    + "_"
    + str(num_dates[0])
    + ".pdf"
)  # pdf
# plt.show()                               # to screen
