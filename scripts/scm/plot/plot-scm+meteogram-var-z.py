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
    Variable-Vertical plot
    Run as: python3 plot-scm+meteogram-var-z.py temp 15  (varname hour)
Author
    Martin Koehler, DWD
Revision history
    202112 -- Initial version
"""

import datetime
import sys

import matplotlib.pyplot as plt
import mpl_toolkits as mp
import numpy as np
from netCDF4 import Dataset, date2num, num2date

# sys.path.insert(1, '/path/to/application/app/folder')
sys.path.append("/hpc/uhome/mkoehler/bin/python")
from autoscale import autoscale

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
    if varname == "temp":
        varname = "t"  # different meteogram names
    if varname == "tot_qc_dia":
        varname = "qc_dia"
    if varname == "tot_qi_dia":
        varname = "qi_dia"
    if n2d3d == 2:
        variable_char = nc2_fid.variables["sfcvar_name"][:, :]
    else:
        variable_char = nc2_fid.variables["var_name"][:, :]
    variable_names = variable_char.view(
        "S128"
    )  # convert 128 character array to string (byte)
    variable_names = variable_names.astype("U")  # convert byte to string
    # print(variable_names)
    # print(varname)
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


def meteogram_longname(nc2_fid, varname, n_variable, n2d3d):
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
    varname_long = varname_long.item()
    if varname == "tot_qc_dia":
        varname_long = "qc diagnostic"
    if varname == "tot_qi_dia":
        varname_long = "qi diagnostic"
    print("Variable long name: ", varname_long)
    return varname_long


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
hour = float(sys.argv[2])

# -------------------------------------------------------------------
# setup

# varname      = 'tot_qc_dia'  # theta_v, temp, u, qv, tot_qc_dia, clc, rel_hum

# ICON global meteogram output:
location = "PayerneNE"  # meteogram station_name
nc2_file = "/hpc/uwork/mkoehler/run-icon/experiments/exp_003_2021111900/METEOGRAM_patch001.nc"
# npoints      = (0,)          # one point (0 is 1st point)
# npoints      = range(32)     # plot all 32 points
npoints = (-1,)  # mean of 32 points

nc2_fid = Dataset(nc2_file, "r")
metgrm_2d3d = "values"  # 3D variables
mg_n_station = meteogram_station(nc2_fid, location)
mg_n_variable = meteogram_variable(nc2_fid, varname, 3)
varname_axis = meteogram_longname(nc2_fid, varname, mg_n_variable, 3)

# SCM output:
nc_file = "/hpc/uwork/mkoehler/run-icon/scm/SCM_ICON_Payerne/scm_out_ML_20211119T000000Z.nc"
nc_fid = Dataset(nc_file, "r")

# output steps:
hours_plot = np.array([hour])
# 19UTC for first day
# hours_plot   = np.array([19])
# 15UTC for 3 days
# hours_plot  = np.array([15,15+24,15+48])
# steps        = (0,)                  # time step (0 is 1st step)
# steps        = (0,4,8,12)   # time step (0 is 1st step)
# steps_mtg    = (0,2,4,6)

# -------------------------------------------------------------------
# read file

vardata = nc_fid.variables[varname][:]
vardata2 = nc2_fid.variables[metgrm_2d3d][:, :90, mg_n_variable, mg_n_station]

# wind speed
# vardata__v  = nc_fid.variables['v'][:]
# vardata2_v = nc2_fid.variables[metgrm_2d3d][:,:90,7-1,mg_n_station]
# vardata  = ( vardata  ** 2 + vardata__v ** 2 ) ** 0.5
# vardata2 = ( vardata2 ** 2 + vardata2_v ** 2 ) ** 0.5

geopot = nc_fid.variables["geopot"][:]

print("variable shapes", vardata.shape, vardata2.shape)  # vardata3.shape,

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

# find steps for hours_plot
steps = np.empty((len(hours_plot),), dtype=int)
steps_mtg = np.empty((len(hours_plot),), dtype=int)
nn = 0
for h_plot in hours_plot:
    steps[nn] = h_plot / (hours[1] - hours[0])
    steps_mtg[nn] = h_plot / (hours2[1] - hours2[0])
    nn = nn + 1
print("SCM hours and steps:       ", hours[steps[0]], steps)
print("Meteogram hours and steps: ", hours2[steps_mtg[0]], steps_mtg)

# -- unit formatting:  convert [W m-2] to [$W m^{-2}$] etc
unit = nc_fid.variables[varname].units
for r in (("-", "^{-"), ("1", "1}"), ("2", "2}")):
    unit = unit.replace(*r)
if unit != "%":
    unit = "[$" + unit + "$]"

# -------------------------------------------------------------------
# plot

title = "ICON SCM and global forecasts"
# xtitle = "%s [%s]" % (nc_fid.variables[varname].standard_name,\
#                      nc_fid.variables[varname].units)
# xtitle = "%s [%s]" % (varname_axis,nc_fid.variables[varname].units)
xtitle = varname_axis + "  " + unit

ytitle = "height above ground [$m$]"
ytitle_right = "model levels"

fig, ax = plt.subplots(figsize=(6, 6))
fig.subplots_adjust(bottom=0.15, left=0.2)

ax.set_title(title).set_fontsize(14)
# ax.set_title('(c)',fontsize=16, fontweight='bold',loc='right')
ax.set_xlabel(xtitle).set_fontsize(14)
ax.set_ylabel(ytitle).set_fontsize(14)
# ax.set_xlim([285,300])    # theta_v
# ax.set_xlim([260,290])    # temp
# ax.set_xlim([0.0,0.007])  # qv
# ax.set_xlim([0.0,0.001])  # qc_dia
# ax.set_xlim([0,10])       # u
# ax.set_xlim([-10,2])      # v
# ax.set_xlim([-20,20])     # speed
ax.set_ylim([0, 3000])

plt.setp(ax.get_xticklabels(), fontsize=14)  # xtick labels fontsize
plt.setp(ax.get_yticklabels(), fontsize=14)  # ytick labels fontsize
plt.rc("legend", fontsize=13)  # legend fontsize

colors = ["mediumblue", "forestgreen", "red", "gold", "purple"]

grav = 9.80665
datay = geopot[0, :, 0] / grav  # height vertical axis
datay = datay - datay[len(datay) - 1]  # height above ground

# plot SCM lines

nc = 0
for ns in steps:
    for np in npoints:
        if np >= 0:
            datax = vardata[ns, :, np - 1]
        # datay = geopot [ns,:,np-1] / grav
        else:
            datax = vardata[ns, :, :].mean(axis=1)
        # datay = geopot [ns,:,:].mean(axis=1) / grav
        ax.plot(
            datax,
            datay,
            linewidth=1.5,
            # color=colors[nc], label=str(hours[ns])+" h")   # (" + str(num_dates[ns]) + ")")
            label="SCM",
        )  # (' + str(num_dates[ns]) + ')')
        ax.margins(y=0)
        nc = nc + 1

# right axis with model level locations

ax2 = ax.twinx()  # second axes
ax2.set_yticks(datay)
ax2.tick_params(axis="y", labelright=False)
# ax2.tick_params(axis='y')
ax2.set_ylabel(ytitle_right).set_fontsize(14)
ax2.set_ylim(ax.get_ylim())  # use y-limits as in ax


# plot ICON global meteogram lines

nc = 0
for ns in steps_mtg:
    print(vardata2.shape, geopot.shape)
    datax = vardata2[ns, :]
    # datay = geopot  [ns,:,0] / grav
    ax.plot(
        datax,
        datay,
        linewidth=2.0,
        # color=colors[nc]), linestyle="dotted"
        label="global",
    )  # (' + str(hours2[ns]) + ')')
    nc = nc + 1


# label for each line

ax.legend()

autoscale(ax, "x", margin=0.1)  # automatic scaling of x-axis given fixed y_lim

fig.savefig(
    "var-z_"
    + varname
    + "_"
    + location
    + "_"
    + str(num_dates[0])
    + "+"
    + str(hour)
    + "h.png"
)  # png
# fig.savefig('var-z_'+varname+'_'+hour+'h.pdf')      # pdf
# plt.show()                                 # to screen
