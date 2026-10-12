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
    Read and plot SCM from NetCDF file with Python3
Details
    Time-Vertical plot
    Run as: python3 plot-scm-time-z.py RICO SCM tot_qc_dia
Author
    Martin Koehler, DWD
Revision history
    20200701 -- Initial version
"""

import glob  # wildcard file completion
import sys
from copy import copy  # color scale copy

import matplotlib.pyplot as plt
import mpl_toolkits as mp
import numpy as np
from netCDF4 import Dataset, date2num, num2date

# -------------------------------------------------------------------
# arguments

print("")
print("Number of arguments:", len(sys.argv), "arguments.")
print("Argument List:      ", str(sys.argv))
case = sys.argv[1]
model = sys.argv[2]
varname = sys.argv[3]

# -------------------------------------------------------------------
# setup

# varname = 'tot_qc_dia'      # temp, u, v, w, clc, rh, tot_qv_dia, tot_qc_dia, tot_qi_dia
# npoint  = 0                 # one point (0 is 1st point)
npoint = -2  # mean of 32 points
# case    = 'RICO'
# model   = 'LES' or 'SCM'
# modname = ''
# modname = '_edmf'
modname = "_test"


if model == "SCM":
    nc_file = (
        "/hpc/uwork/mkoehler/run-icon/scm/SCM_"
        + case
        + "_dephy"
        + modname
        + "/out_SCM_"
        + case
        + "_dephy"
        + modname
        + "_ML_*Z_mean.nc"
    )
    nc_zfile = nc_file
else:
    nc_file = (
        "/hpc/uwork/mkoehler/run-icon/scm/LES_100x100_"
        + case
        + "_dephy/out_ml_LES_100x100_"
        + case
        + "_dephy_ML_*Z_mean.nc"
    )
    nc_zfile = (
        "/hpc/uwork/mkoehler/run-icon/scm/LES_100x100_"
        + case
        + "_dephy/out_z_LES_100x100_"
        + case
        + "_dephy_ML_*Z_mean.nc"
    )

nc_file = glob.glob(nc_file)  # wildcard completion
nc_zfile = glob.glob(nc_zfile)  # wildcard completion
print("reading file: ", nc_file, "and", nc_zfile)
nc_fid = Dataset(nc_file[0], "r")
nc_fidz = Dataset(nc_zfile[0], "r")

# read file

dates = num2date(nc_fid.variables["time"][:], nc_fid.variables["time"].units)
vardata = nc_fid.variables[varname][:]
# geopot  = nc_fid.variables['geopot'][:]
z_mc = nc_fidz.variables["z_mc"][:]
print(
    "variable shapes: vardata, z_mc, dates ",
    vardata.shape,
    z_mc.shape,
    dates.shape,
)

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

# -------------------------------------------------------------------
# plot

title = (
    "ICON "
    + model
    + " forecast    %s" % (nc_fid.variables[varname].standard_name)
)
xtitle = "Simulation time since %s [h]" % (num_dates[0])
ytitle = "height above ground [m]"

fig, ax = plt.subplots(figsize=(8, 5))
fig.subplots_adjust(bottom=0.15, right=1.0)

top_height = 5000

ax.set_title(title).set_fontsize(16)
ax.set_xlabel(xtitle).set_fontsize(14)
ax.set_ylabel(ytitle).set_fontsize(14)
ax.set_ylim([0, top_height])

# grav  = 9.80665
# datay = geopot  [0,:,0] / grav           # height vertical axis
datay = np.squeeze(z_mc)  # height vertical axis
datay = datay - datay[len(datay) - 1]  # height above ground

if npoint >= 0:
    dataz = vardata[:, :, npoint]
elif npoint == -1:
    dataz = vardata[:, :, :].mean(axis=2)
else:
    dataz = np.squeeze(vardata)

[X, Y] = np.meshgrid(hours, datay)  # creating 2-D grid
Z = dataz
Z = np.swapaxes(Z, 0, 1)

print("shape of X, Y, Z: ", X.shape, Y.shape, Z.shape)

lev_vis = np.where(datay < top_height)  # levels that are visible in plot area
X = np.squeeze(X[lev_vis, :])
Y = np.squeeze(Y[lev_vis, :])
Z = np.squeeze(Z[lev_vis, :])

print("shape of X, Y, Z: ", X.shape, Y.shape, Z.shape)

# color map options                       # 'Spectral', 'jet', 'Blues_r'  (r for reversed, or cmap.reversed())

if varname == "tot_qc_dia" or varname == "clc":
    cmap = plt.cm.get_cmap("Blues")
    cmap.set_under("white", 1.0)  # 1.0 represents not transparent
    CS = ax.pcolormesh(
        X, Y, Z, cmap=cmap, vmin=0.00000001
    )  # set small values to white
else:
    cmap = plt.cm.get_cmap(
        "jet"
    )  # 'Spectral', 'jet', 'Blues_r'  (r for reversed, or cmap.reversed())
    CS = ax.pcolormesh(X, Y, Z, cmap=cmap)  # set small values to white

# CS = ax.contourf(X, Y, Z, 30, cmap=cmap, extend='both')            # plot filled contour plot
# CS = ax.pcolormesh(X,Y,Z, cmap=cmap, norm=colors.CenteredNorm())   # xsfor centered color range:
# CS = ax.pcolormesh(X,Y,Z, cmap=cmap, vmin=.00000001)               # set small values to white

# color bar

cbar = fig.colorbar(CS)
cbar.ax.set_ylabel("[" + nc_fid.variables[varname].units + "]")

# right axis with model level locations

ax2 = ax.twinx()  # second axes location
ax2.set_yticks(datay)
ax2.tick_params(axis="y", labelright=False)
ax2.set_ylim(ax.get_ylim())  # use y-limits as in ax

# screen, png, pdf

dirplots = "plots/" + case + modname + "/"
fig.savefig(
    dirplots
    + "time-z_"
    + varname
    + "_"
    + case
    + modname
    + "_"
    + model
    + "_"
    + str(num_dates[0])
    + ".png"
)  # png
fig.savefig(
    dirplots
    + "time-z_"
    + varname
    + "_"
    + case
    + modname
    + "_"
    + model
    + "_"
    + str(num_dates[0])
    + ".pdf"
)  # pdf
# plt.show()                                # to screen
