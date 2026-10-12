#!/bin/bash

# ICON
#
# ----------------------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
#
# Merge plots from SCM and LES runs
#
# runs on workstation
#
# Martin Koehler, Jan. 2022
# ----------------------------------------------------------------------------

set -ex

cd /home/mkoehler/rcl/home/mkoehler/icon/icon-nwp-test2/scripts/scm/plot/plots

modname=''
#modname='_edmf'

#for case in RICO ARM BOMEX GABLS1 FIRE MPACE MAGIC ; do
for case in MAGIC ; do

  cd ${case}${modname}

  pdfjam -q --outfile ${case}-plots-1.pdf --nup 2x4 time-z_clc*pdf time-z_tot*pdf time-z_theta_v*pdf time-z_u*pdf time-var*pdf
  pdfjam -q --outfile ${case}-plots-2.pdf --nup 2x2 var-z*pdf

  pdfunite ${case}-plots-1.pdf ${case}-plots-2.pdf ${case}${modname}-plots.pdf

  cp ${case}${modname}-plots.pdf ..
  cd ..

done




# works only on workstation:
#pdfjam --outfile RICO-plots.pdf --suffix nup --nup 2x4 time-z_clc*pdf time-z_tot*pdf time-z_theta_v*pdf time-z_u*pdf time-var*pdf var-z*pdf


# convert -append time-var_tqc_dia_RICO_2004121600.png  time-z_tot_qc_dia_RICO_SCM_2004121600.png var-z_theta_v_RICO_2004121600.png out_left.png

# convert -append time-var_clct_RICO_2004121600.png time-z_tot_qc_dia_RICO_LES_2004121600.png var-z_tot_qc_dia_RICO_2004121600.png out_right.png

# convert +append out_left.png out_right.png out.png


#pdfjam --suffix nup --nup 2x2 input.pdf

#montage -mode concatenate -tile NxM in-*.pdf out.pdf
