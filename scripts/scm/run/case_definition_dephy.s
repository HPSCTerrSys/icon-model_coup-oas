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
# SCM definition for DEPHY cases:
#
# Martin Koehler, Jan. 2022
# ----------------------------------------------------------------------------

# argument:

dephy_case=$1
echo "--- Selected case: "${dephy_case}

set -x

# ----------------------------------------------------------------------------
# default parameters

echo "--- Setting default namelist values"

lscm_read_tke=.FALSE.            # read initial tke from netcdf

inwp_gscp=1
inwp_convection=1                # 1:Tiedtke/Bechtold 0:apdf
inwp_radiation=4                 # 1:RRTM radiation
inwp_cldcover=1                  # 3: clouds from COSMO SGS cloud scheme 0:apdf
inwp_turb=1                      # 1: COSMO diffusion and transfer
inwp_satad=1
inwp_sso=1
inwp_gwd=1

nhours=24
start_date="2000-01-01T00:00:00Z"
end_date="2000-01-02T00:00:00Z"
corio_lat=0.0

lscm_read_tke=.FALSE.            # read initial tke from netcdf
lscm_read_z0=.FALSE.             # read initial z0 from netcdf
scm_sfc_mom=0                    # 0: TURBTRANS: no prescribed u*
scm_sfc_temp=0                   # 0: TERRA:     no prescribed sensible heat flux at surface
scm_sfc_qv=0                     # 0: TERRA:     no prescribed latent heat flux at surface

is_subsidence_moment=.FALSE.
is_subsidence_heat=.FALSE.
is_advection=.FALSE.
is_advection_uv=.FALSE.
is_advection_tq=.FALSE.
is_geowind=.FALSE.
is_rad_forcing=.FALSE.
is_nudging=.FALSE.
is_nudging_uv=.FALSE.
is_nudging_t=.FALSE.
is_nudging_q=.FALSE.
nudge_start_height_t=0.0
nudge_full_height_t=0.0
nudge_start_height_q=0.0
nudge_full_height_q=0.0
nudge_start_height_uv=0.0
nudge_full_height_uv=0.0
dt_relax_uv=43200.0
dt_relax_t=43200.0
dt_relax_q=43200.0
is_sim_rad=.FALSE.

# ----------------------------------------------------------------------------
# case dependent parameters

echo "--- Setting case dependent namelist values"

case $dephy_case in
   # -------------------------------------------------------------------------
   ARM )
     nhours=14.5
     start_date="1997-06-21T00:00:00Z"
     end_date="1997-06-21T14:30:00Z"
     file_dephy=ARMCU_REF_SCM_driver.nc
     corio_lat=35.00
     inwp_radiation=0
     inwp_sso=0
     inwp_gwd=0
     lscm_read_z0=.TRUE.              # read initial z0 from netcdf
     scm_sfc_mom=0                    # 0: TURBTRANS: no prescribed u*
     scm_sfc_temp=1                   # 0: TERRA:     no prescribed sensible heat flux at surface
     scm_sfc_qv=3                     # 0: TERRA:     no prescribed latent heat flux at surface
     is_advection=.TRUE.
     is_advection_tq=.TRUE.
     is_geowind=.TRUE.
     ;;

   # -------------------------------------------------------------------------
   RICO )   # Question: which surface forcing?
     nhours=24
     start_date="2004-12-16T00:00:00Z"
     end_date="2004-12-17T00:00:00Z"
     file_dephy=RICO_SHORT_SCM_driver.nc
     corio_lat=18.00
     lscm_read_z0=.FALSE.             # read initial z0 from netcdf
     scm_sfc_mom=0                    # 0: no prescribed surface flux
     scm_sfc_temp=1                   # 1: prescribed Ts
     scm_sfc_qv=3                     # 3: qv_s based on saturation
     is_subsidence_moment=.TRUE.
     is_subsidence_heat=.TRUE.
     is_advection=.TRUE.
     is_advection_tq=.TRUE.
     is_geowind=.TRUE.
     ;;

   # -------------------------------------------------------------------------
   BOMEX )   # Question: which surface forcing?
     nhours=24
     start_date="1969-06-26T00:00:00Z"
     end_date="1969-06-27T00:00:00Z"
     file_dephy=BOMEX_REF_SCM_driver.nc
     corio_lat=15.0
     lscm_read_tke=.TRUE.             # read initial tke from netcdf
     lscm_read_z0=.FALSE.             # read initial z0 from netcdf
     scm_sfc_mom=2                    # prescribed u*
     scm_sfc_temp=2                   # prescribed sensible heat flux at surface
     scm_sfc_qv=2                     # prescribed latent heat flux at surface
     is_subsidence_moment=.TRUE.
     is_subsidence_heat=.TRUE.
     is_advection=.TRUE.
     is_advection_tq=.TRUE.
     is_geowind=.TRUE.
     ;;

   # -------------------------------------------------------------------------
   GABLS1 )
     nhours=9
     start_date="1969-06-26T00:00:00Z"
     end_date="1969-06-26T09:00:00Z"
     file_dephy=GABLS1_REF_SCM_driver.nc
     corio_lat=73.00
     lscm_read_z0=.TRUE.              # read initial z0 from netcdf
     scm_sfc_mom=5                    # prescribed u at surface - 0 :MO used
     scm_sfc_temp=5                   # prescribed temperature at surface: MO used
     scm_sfc_qv=5                     # prescribed moisture - dry = 0: MO used
     is_geowind=.TRUE.
     inwp_gscp=0
     inwp_convection=0                # 1:Tiedtke/Bechtold 0:apdf
     inwp_radiation=0                 # 1:RRTM radiation
     inwp_cldcover=0                  # 3: clouds from COSMO SGS cloud scheme 0:apdf
     inwp_turb=1                      # 1: COSMO diffusion and transfer
     inwp_satad=0
     inwp_sso=0
     inwp_gwd=0
     icapdcycl=0                      # 3: apply CAPE modification to improve diurnalcycle over tropical land (optimizes NWP scores)
     ;;

   # -------------------------------------------------------------------------
   FIRE )     # QUESTION: how to specify z0=0.0002m (Dynkerke, et al, 2004)
     nhours=37
     start_date="1987-07-14T08:00:00Z"
     end_date="1987-07-15T21:00:00Z"
     file_dephy=FIRE_REF_SCM_driver.nc
     corio_lat=33.3
     lscm_read_z0=.FALSE.             # read initial z0 from netcdf
     scm_sfc_mom=0                    # 0: no prescribed surface flux
     scm_sfc_temp=1                   # 1: prescribed Ts
     scm_sfc_qv=3                     # 3: qv_s based on saturation
     is_subsidence_moment=.TRUE.
     is_subsidence_heat=.TRUE.
     is_advection=.TRUE.
     is_advection_tq=.TRUE.
     is_geowind=.TRUE.
     ;;

   # -------------------------------------------------------------------------
   MPACE )    # not tested
     nhours=13
     start_date="2004-10-09T17:00:00Z"
     end_date="2004-10-10T05:00:00Z"
     file_dephy=MPACE_REF_SCM_driver.nc
     corio_lat=71.75
     lscm_read_z0=.FALSE.             # read initial z0 from netcdf
     scm_sfc_mom=0                    # 0: no prescribed surface flux
     scm_sfc_temp=1                   # 1: prescribed Ts
     scm_sfc_qv=3                     # 3: qv_s based on saturation
     is_subsidence_moment=.TRUE.
     is_subsidence_heat=.TRUE.
     is_advection=.TRUE.
     is_advection_tq=.TRUE.
     is_geowind=.TRUE.
     ;;

   # -------------------------------------------------------------------------
   MAGIC )   # Question: which surface forcing?
     nhours=106
     start_date="2013-07-20T17:30:00Z"
     end_date="2013-07-25T03:30:00Z"
##   end_date="2013-07-20T18:00:00Z"
     file_dephy=MAGIC_LEG15A_SCM_driver.nc
     corio_lat=18.00
     scm_sfc_mom=0                    # 0: no prescribed surface flux
     scm_sfc_temp=1                   # 1: prescribed Ts
     scm_sfc_qv=3                     # 3: qv_s based on saturation
#    is_subsidence_moment=.TRUE.
     is_subsidence_heat=.TRUE.
     is_advection=.TRUE.
     is_advection_tq=.TRUE.
     is_geowind=.TRUE.
     is_nudging=.TRUE.
     is_nudging_uv=.TRUE.
     is_nudging_t=.TRUE.
     is_nudging_q=.TRUE.
     nudge_start_height_uv=0.0
     nudge_full_height_uv=0.0
     nudge_start_height_t=3000.0
     nudge_full_height_t=3000.0
     nudge_start_height_q=3000.0
     nudge_full_height_q=3000.0
     dt_relax_uv=43200.0
     dt_relax_t=1800.0
     dt_relax_q=1800.0
     ;;

   # -------------------------------------------------------------------------
   COMBLE )     # QUESTION: how to specify z0=0.0002m (Dynkerke, et al, 2004)
     nhours=37
     start_date="2020-03-12 22:00:00"
     end_date="2020-03-13 18:00:00"
     file_dephy=COMBLE_INTERCOMPARISON_FORCING_V2.5.nc
     inwp_gscp=4
     corio_lat=74.5
     lscm_read_z0=.FALSE.             # read initial z0 from netcdf
     scm_sfc_mom=0                    # 0: no prescribed surface flux
     scm_sfc_temp=1                   # 1: prescribed Ts
     scm_sfc_qv=3                     # 3: qv_s based on saturation
     is_subsidence_moment=.FALSE.  #.TRUE.
     is_subsidence_heat=.FALSE.    #.TRUE.
     is_advection=.FALSE.          #.TRUE.
     is_advection_tq=.FALSE.       #.TRUE.
     is_geowind=.TRUE.
     #is_nudging=.TRUE.
     #is_nudging_t=.TRUE.
     #dt_relax_t=3600.0
     #nudge_full_height_t=3000.0
     ;;

   # -------------------------------------------------------------------------
   CP-MIP )
     nhours=34
     start_date="2020-02-01T16:00:00Z"
     end_date="2020-02-01T17:00:00Z"  # short test
#     end_date="2020-02-03T02:00:00Z" # full experiment
     file_dephy=FLOWER_SCM_medium_driver_u_uadv_ug.nc
     corio_lat=13.2
     scm_sfc_mom=5                    # 0: no prescribed surface flux
     scm_sfc_temp=5                   # 1: prescribed Ts
     scm_sfc_qv=6                     # 3: qv_s based on saturation
     is_subsidence_moment=.TRUE.
     is_subsidence_heat=.TRUE.
     is_advection=.TRUE.
     is_advection_tq=.TRUE.
     is_advection_uv=.TRUE.
     is_geowind=.TRUE.
     is_nudging=.TRUE.
     is_nudging_uv=.TRUE.
     is_nudging_t=.TRUE.
     is_nudging_q=.TRUE.
     nudge_start_height_uv=5250.0
     nudge_full_height_uv=5250.0
     nudge_start_height_t=5250.0
     nudge_full_height_t=5250.0
     nudge_start_height_q=5250.0
     nudge_full_height_q=5250.0
     dt_relax_uv=1800.0
     dt_relax_t=1800.0
     dt_relax_q=1800.0
     ;;

   # -------------------------------------------------------------------------
   # AMMA, AYOTTE, IHOP, GABLS4, DYNAMO, SCMS, SANDU

   * )
     echo "---------- Case not yet defined! ----------"
     ;;

esac
