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
    Time-Variable plot
    Run as: python3 plot-scm+meteogram-time-var.py shfl_s
Author
    Martin Koehler, DWD
Revision history
    202112 -- Initial version
"""

import datetime
import getopt
import math
import os
import sys  # file handling and operators

import matplotlib.pyplot as plt
import mpl_toolkits as mp
import numpy as np
from netCDF4 import Dataset, date2num, num2date
from scipy.signal import savgol_filter

# -------------------------------------------------------------------
# functions: meteogram variables


def meteogram_station(nc2_fid, station):
    """Find meteogram station number
    input:  netcdf file ID, station name
    output: station number starting with 0"""
    station_char = nc2_fid.variables["station_name"][:, :]
    station_names = station_char.view(
        "S128"
    )  # convert 128 character array to string (byte)
    station_names = station_names.astype("U")  # convert byte to string
    # print(station_names)
    for n_station in range(len(station_names)):
        # print(n_station, station_names[n_station])
        if station_names[n_station] == station:
            break
    print("Station name:  ", station_names[n_station], n_station)
    return n_station


def meteogram_variable(nc2_fid, varname, n2d3d):
    """Find meteogram variable number
    input:  netcdf file ID, variable name, n2d3d=2 or 3
    output: variable number starting with 0"""
    if n2d3d == 2:
        variable_char = nc2_fid.variables["sfcvar_name"][:, :]
    else:
        variable_char = nc2_fid.variables["var_name"][:, :]
    variable_names = variable_char.view(
        "S128"
    )  # convert 128 character array to string (byte)
    variable_names = variable_names.astype("U")  # convert byte to string
    for n_variable in range(len(variable_names)):
        if variable_names[n_variable].item().lower() == varname:
            break
    print(
        "Variable name: ",
        variable_names[n_variable].item(),
        varname,
        n_variable,
    )
    return n_variable


def meteogram_longname(nc2_fid, n_variable, n2d3d):
    """Find meteogram variable long name
    input:  netcdf file ID, variable name, n2d3d=2 or 3
    output: variable long name"""
    if n2d3d == 2:
        variable_char = nc2_fid.variables["sfcvar_long_name"][n_variable, :]
    else:
        variable_char = nc2_fid.variables["var_long_name"][n_variable, :]
    varname_long = variable_char.view(
        "S128"
    )  # convert 128 character array to string (byte)
    varname_long = varname_long.astype("str")  # convert byte to string
    if varname == "tqc_dia":
        varname_long[0] = "tqc diagnostic"
    print("Variable long name: ", varname_long)
    return varname_long.item()


def meteogram_dates(nc2_fid):
    """Convert meteogram date to datetime format
    input:  netcdf file ID
    output: dates array in datetime format"""
    dates_mg = nc2_fid.variables["date"][:]
    dates_mg = dates_mg.view(
        "S128"
    )  # convert 128 character array to string (byte)
    dates_mg = dates_mg.astype("str")  # convert byte to string
    dates = np.empty((len(dates_mg),), dtype=datetime.datetime)
    for nn in range(len(dates)):
        dates[nn] = datetime.datetime.strptime(
            dates_mg[nn].item(), "%Y%m%dT%H%M%SZ"
        )
        # 20211119T000000Z   2021-11-19 00:00:00+00:00
    print("Meteogram last date: ", dates_mg[nn], dates[nn])
    return dates


# -------------------------------------------------------------------
# arguments

print("")
print("Number of arguments:", len(sys.argv), "arguments.")
print("Argument List:      ", str(sys.argv))
varname = sys.argv[1]

# -------------------------------------------------------------------
# setup

# varname      = 'shfl_s'      # temp, t_2m, qv_2m, u_10m, v_10m, shfl_s, lhfl_s, clct
level = 1  # 1000hPa, 90 (model level), 1 (e.g. t_2m)
# npoints      = (1,)          # one point (0 is 1st point)
# npoints      = range(32)     # plot all 32 points
npoints = (-1,)  # mean of 32 points

# ICON global meteogram output:
location = "PayerneNE"  # meteogram station_name
n2d3d = 2  # 2/3: 2D or 3D variable
nc2_file = "/hpc/uwork/mkoehler/run-icon/experiments/exp_003_2021111900/METEOGRAM_patch001.nc"

nc2_fid = Dataset(nc2_file, "r")
if n2d3d == 2:
    metgrm_2d3d = "sfcvalues"  # 2D variables
else:
    metgrm_2d3d = "values"  # 3D variables
mg_n_station = meteogram_station(nc2_fid, location)
mg_n_variable = meteogram_variable(nc2_fid, varname, n2d3d)
varname_axis = meteogram_longname(nc2_fid, mg_n_variable, n2d3d)

# SCM output:
nc_file = "/hpc/uwork/mkoehler/run-icon/scm/SCM_ICON_Payerne/scm_out_ML_20211119T000000Z.nc"
nc_fid = Dataset(nc_file, "r")


# -------------------------------------------------------------------
# read file

vardata = nc_fid.variables[varname][:]
vardata2 = nc2_fid.variables[metgrm_2d3d][:, mg_n_variable, mg_n_station]
vardata2 = vardata2.reshape(vardata2.shape[0], 1, 1)

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

if varname == "clct":
    vardata2 = vardata2 * 100.0

print("variable shapes", vardata.shape, vardata2.shape)

# flatten SCM model
if n2d3d == 3:
    vardata = vardata[:, level - 1, :]
vardata = np.squeeze(vardata)


# -------------------------------------------------------------------
# time conversion

dates = num2date(nc_fid.variables["time"][:], nc_fid.variables["time"].units)
dates2 = meteogram_dates(nc2_fid)

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
hours2 = [
    (d.day - dates2[0].day) * 24
    + (d.hour - dates2[0].hour)
    + (d.minute / 60)
    + (d.second / 3600)
    for d in dates2
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

title = "ICON SCM and global forecasts      " + location
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
ax.set_xlim(0, 72)
# ax.set_ylim(281,297)       # t_2m
# ax.set_ylim(0.004,0.008)   # qv_2m
# ax.set_ylim(-4,1)          # v_10m
# ax.set_ylim(-300,10)       # lhfl_s
# ax.set_ylim(-400,50)       # shfl_s

plt.setp(ax.get_xticklabels(), fontsize=14)  # xtick labels fontsize
plt.setp(ax.get_yticklabels(), fontsize=14)  # ytick labels fontsize
plt.rc("legend", fontsize=13)  # legend fontsize

# plot SCM line

for np in npoints:
    if np >= 0:
        # data = vardata[:,level-1,np-1]
        data = vardata[:, np - 1]
    else:
        # data = vardata[:,level-1,:].mean(axis=1)
        data = vardata[:, :].mean(axis=1)
        ax.plot(hours, data, linewidth=2.0, label="SCM")

# plot ICON global meteogram line

ax.plot(hours2, vardata2[:, level - 1], linewidth=2.0, label="global")

# label for each line

ax.legend()


fig.savefig(
    "time-var_" + varname + "_" + location + "_" + str(num_dates[0]) + ".png"
)  # png
# fig.savefig('time-var_'+varname+'.pdf')  # pdf
# plt.show()                               # to screen
