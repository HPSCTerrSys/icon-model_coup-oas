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

: "${ICON_DIR:?ICON_DIR must be defined before sourcing this script.}"

# Join all arguments using the separator given by the first argument.
# Usage: join_arr <sep> <item0> <item1>...
function join_arr {
  local IFS="$1"; shift
  printf '%s' "$*"
}

# ICON
INLINE_LIST_ICON=(
  src/advection/mo_advection_utils.f90
  src/atm_phy_aes/mo_aes_thermo.f90
  src/atm_phy_nwp/mo_util_phys.f90
  src/atm_phy_schemes/cloud_random_numbers.f90
  src/atm_phy_schemes/mo_2mom_mcrph_driver.f90
  src/atm_phy_schemes/mo_2mom_mcrph_processes.f90
  src/atm_phy_schemes/mo_2mom_mcrph_setup.f90
  src/atm_phy_schemes/mo_2mom_mcrph_util.f90
  src/atm_phy_schemes/mo_aerosol_sources.f90
  src/atm_phy_schemes/mo_cpl_aerosol_microphys.f90
  src/atm_phy_schemes/mo_albedo.f90
  src/atm_phy_schemes/mo_cufunctions.f90
  src/atm_phy_schemes/mo_thdyn_functions.f90
  src/atm_phy_schemes/mo_turb_vdiff.f90
  src/atm_phy_schemes/random_rewrite.f90
  src/atm_phy_schemes/turb_utilities.f90
  src/configure_model/mo_parallel_config.f90
  src/lnd_phy_nwp/mo_nwp_sfc_interp.f90
  src/lnd_phy_schemes/sfc_flake.f90
  src/lnd_phy_schemes/sfc_seaice.f90
  src/parallel_infrastructure/mo_extents.f90
  src/shared/mo_statistics.f90
  src/shared/mo_loopindices.f90
)
INLINE_LIST_ICON=("${INLINE_LIST_ICON[@]/#/${ICON_DIR}/}")

# ICONMATH
INLINE_LIST_ICONMATH=(
  externals/iconmath/src/support/mo_math_utilities.F90
  externals/iconmath/src/support/mo_lib_loopindices.f90
)
INLINE_LIST_ICONMATH=("${INLINE_LIST_ICONMATH[@]/#/${ICON_DIR}/}")

# ART
INLINE_LIST_ART=(
  externals/art/aerosol_dynamics/mo_art_aerosol_utilities.f90
  externals/art/shared/mo_art_modes.f90
  externals/art/tools/mo_art_clipping.f90
)
INLINE_LIST_ART=("${INLINE_LIST_ART[@]/#/${ICON_DIR}/}")

# DACE
INLINE_LIST_DACE=(
  externals/dace/fetch/src/src_for_icon/mo_physics.f90
)
INLINE_LIST_DACE=("${INLINE_LIST_DACE[@]/#/${PWD}/}")

# ECRAD
INLINE_LIST_ECRAD=(
  externals/ecrad/radiation/radiation_two_stream.F90
  externals/ecrad/radiation/radiation_liquid_optics_socrates.F90
  externals/ecrad/radiation/radiation_ice_optics_fu.F90
)
INLINE_LIST_ECRAD=("${INLINE_LIST_ECRAD[@]/#/${ICON_DIR}/}")

# EMVORADO
INLINE_LIST_EMVORADO=(
  externals/emvorado/src_emvorado/radar_gamma_functions_vec.f90
  externals/emvorado/src_emvorado/radar_mie_meltdegree.f90
  externals/emvorado/src_emvorado/radar_mie_specint.f90
  externals/emvorado/src_emvorado/radar_mie_utils.f90
  externals/emvorado/src_emvorado/radar_mielib_vec.f90
  externals/emvorado/src_emvorado/radar_model2rays.f90
  externals/emvorado/src_emvorado/radar_utilities.f90
  externals/emvorado/src_emvorado/radar_dmin_wetgrowth.f90
  externals/emvorado/src_iface_icon/radar_interface.f90
)
INLINE_LIST_EMVORADO=("${INLINE_LIST_EMVORADO[@]/#/${ICON_DIR}/}")

# OCEAN
INLINE_LIST_OCEAN=(
  externals/iconmath/src/support/mo_math_utilities.F90
  src/advection/mo_advection_utils.f90
  src/ocean/physics/mo_ocean_thermodyn.f90
)
INLINE_LIST_OCEAN=("${INLINE_LIST_OCEAN[@]/#/${ICON_DIR}/}")
