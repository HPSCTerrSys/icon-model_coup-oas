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
    Call python plotting program for multiple plots
Details
    Run as: python3 plot-all.py
Author
    Martin Koehler, DWD
Revision history
    202011 -- Initial version
"""

import subprocess

# cases = ['ARM', 'RICO', 'BOMEX', 'GABLS1', 'FIRE', 'MPACE', 'MAGIC']
# hours = ['14.5','24.0', '24.0' , '9.0'   , '37.0', '12.0' , '90.0' ]
cases = ["FIRE"]
hours = ["37.0"]

variables_2d = [
    "shfl_s",
    "lhfl_s",
    "u_10m",
    "v_10m",
    "qv_2m",
    "t_2m",
    "clct",
    "tqv_dia",
    "tqc_dia",
    "tqi_dia",
    "tot_prec",
    "sp_10m",
]
variables_2d = [
    "tqc_dia",
    "tqi_dia",
    "clct",
    "t_2m",
    "shfl_s",
    "lhfl_s",
    "sp_10m",
    "tot_prec",
]

variables_3d = [
    "u",
    "v",
    "theta",
    "tot_qv_dia",
    "tot_qc_dia",
    "tot_qi_dia",
    "clc",
]
variables_3d = ["tot_qc_dia", "clc", "theta_v", "u", "tot_qv_dia"]

nn = 0
for case in cases[:]:

    for var in variables_2d[:]:
        subprocess.call(["python3", "plot-scm+les-time-var.py", case, var])

    for var in variables_3d[:]:
        subprocess.call(
            ["python3", "plot-scm+les-var-z.py", case, var, hours[nn]]
        )

        for model in ["SCM", "LES"]:
            subprocess.call(
                ["python3", "plot-scm+les-time-z.py", case, model, var]
            )

    nn = nn + 1


# works only on workstation:
# pdfjam --outfile RICO-plots.pdf --suffix nup --nup 2x4 time-z_clc*pdf time-z_tot*pdf time-z_theta_v*pdf time-z_u*pdf time-var*pdf var-z*pdf


# convert -append time-var_tqc_dia_RICO_2004121600.png  time-z_tot_qc_dia_RICO_SCM_2004121600.png var-z_theta_v_RICO_2004121600.png out_left.png

# convert -append time-var_clct_RICO_2004121600.png time-z_tot_qc_dia_RICO_LES_2004121600.png var-z_tot_qc_dia_RICO_2004121600.png out_right.png

# convert +append out_left.png out_right.png out.png


# pdfjam --suffix nup --nup 2x2 input.pdf

# montage -mode concatenate -tile NxM in-*.pdf out.pdf

# works only on workstation:
# pdfjam --outfile RICO-plots.pdf --suffix nup --nup 2x4 time-z_clc*pdf time-z_tot*pdf time-z_theta_v*pdf time-var*pdf var-z*pdf
