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

# Example Script of creating a animated gif in Python.
# I.Kroener 2022-01-26
#
# Requires: module load python/2021.3 (with resulting module switches)
#   module load gcc/9.1.0
#   module load hpcx/hpcx-ompi
#   module load netcdf4/4.7.3-x86-gnu
#   module load hdf5/1.10.5-x86-gnu
#
# Display as: animate file.gif
# ---------------------------------------------------------------

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib import animation

# -- data definition

datestring = "20041216"
filedir = "/hpc/uwork/mkoehler/run-icon/scm/"
EXP = "LES_100x100_RICO_dephy"
crmtemplate = "out_sfc_LES_100x100_RICO_dephy_ML_%sT000000Z.nc"
variable = "t_2m"
selheight = 0

vmin = 298
vmax = 302

# -- Create Inputfile with string replacement

inputfile = os.path.join(filedir, EXP, crmtemplate % (datestring[0:8]))
print("input file: ", inputfile)
data = xr.open_dataset(inputfile)

fig, ax = plt.subplots()

# -- Dictionary for integer selection in xarray

seldict = {"time": 0}
# If variable has a height dimension select only one level (selheight : int, indexing 0)
isheight = [dimname[0:5] in "height" for dimname in data[variable].dims]
if any(isheight):  # if yes
    dimname = np.array(data[variable].dims)[isheight][0]
    seldict[data[dimname].name] = selheight

Z = data[variable].isel(seldict)

cmap = plt.get_cmap("Spectral_r", 10000)
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap))
cnf = ax.tricontourf(
    data.clon.data, data.clat.data, Z, cmap=cmap, vmin=vmin, vmax=vmax
)
plt.title(data.time.isel(time=1).data)

# plt.savefig()

# -- Function to animate


def animate(i):
    seldict["time"] = i  # overwrite timestep in seldict
    z = data["t_2m"].isel(seldict).data
    cnf = ax.tricontourf(
        data.clon.data, data.clat.data, z, vmin=vmin, vmax=vmax, cmap=cmap
    )
    plt.title(data.time.isel(time=i).data)
    return cnf


from matplotlib import animation

Nsteps = len(data.time)  # length of animation
anim = animation.FuncAnimation(fig, animate, frames=Nsteps)
anim.save("LES_%s_movie.gif" % datestring[0:8], writer="imagemagick", fps=4)
