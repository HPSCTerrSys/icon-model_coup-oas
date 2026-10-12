#!/bin/ksh

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

# Interpolate tropospheric aerosol data (Kinne, CMIP6) for ICON runs
# Applicable for irad_aero = 12,13,15,18,19 (for 18 + 19, only the year 1850 is necessary)
##########################################################################################

## User settings

# Begin and end year for interpolation:
BYEAR=1850
EYEAR=1850

RES=<YOUR_RESOLUTION>
OUTDIR=<YOUR_OUTPUT_PATH>/${RES}/aeropt_kinne_${RES}
TARGETGRID=<PATH_TO_YOUR_TARGETGRID>

##########################################################################################

DATADIR='/pool/data/ICON/grids/public/mpim/independent/aerosol_kinne/'
SOURCEGRID=$DATADIR/aeropt_kinne_lw_b16_coa_rast.nc

source /sw/etc/profile.levante
module unload cdo
module load cdo/2.0.6-gcc-11.2.0
module load nco/5.0.6-gcc-11.2.0

mkdir -p $OUTDIR && cd ${OUTDIR}
echo ${OUTDIR}

## gencon or genlaf ???
REMAP_WEIGHTS=weights.nc
cdo selvar,cell_area $TARGETGRID cell_area.nc
cdo genlaf,cell_area.nc $SOURCEGRID $REMAP_WEIGHTS

params="sw_b14_coa sw_b14_fin lw_b16_coa"

for param in $params ; do

   if [ $param = sw_b14_fin ]; then
      lastyear=$EYEAR
   else
      lastyear=1850
   fi

   ## only for sw_b14_fin, loop over all years
   for ((y=$BYEAR; y<=$lastyear; y++)) ; do

     if [ $param = sw_b14_fin ]; then
       IFILE=aeropt_kinne_sw_b14_fin_${y}_rast.nc
     else
       IFILE=aeropt_kinne_${param}_rast.nc
     fi

     OFILE=${OUTDIR}/${RES}_${IFILE}
     cdo -f nc5 -P 8 remap,cell_area.nc,$REMAP_WEIGHTS ${DATADIR}/${IFILE} ${OFILE}

     if [ $param = lw_b16_coa ]; then
     	ncks -A -v lnwl,wl_lo,wl_up,zbot_abs,ztop_abs,delta_z ${DATADIR}/${IFILE} ${OFILE}
     elif [ $param = sw_b14_fin ]; then
           ncks -A -v lnwl,wl_lo,wl_up,wn_lo,wn_up,zbot_abs,ztop_abs,delta_z ${DATADIR}/aeropt_kinne_sw_b14_fin_1850_rast.nc ${OFILE}
     else
     	ncks -A -v lnwl,wl_lo,wl_up,wn_lo,wn_up,zbot_abs,ztop_abs,delta_z ${DATADIR}/${IFILE} ${OFILE}
     fi

     ls ${OFILE}
   done
done

## cleanup
rm cell_area.nc $REMAP_WEIGHTS
