#! /bin/bash

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

# The two arrays expected to be filled with the same amount of entries (variations) to run.
# Each entry in directories will be used as the name of the subdirectory for the respective
# test run and each line in patchlines is part of a Python command to modify the respective
# namelist. Although multiple commands in a patchlines entry are separated by a semicolon,
# no semicolon is allowed at the end

directories=()
patchlines=()

# Assumed to be the default (strictly speaking no changes should be necessary)
directories+=("tke")
patchlines+=("nml['ocean_vertical_diffusion_nml']['vert_mix_type'] = 2")


directories+=("pp0")
patchlines+=("nml['ocean_vertical_diffusion_nml']['vert_mix_type'] = 1; nml['ocean_vertical_diffusion_nml']['ppscheme_type'] = 0")


directories+=("pp4")
patchlines+=("nml['ocean_vertical_diffusion_nml']['vert_mix_type'] = 1; nml['ocean_vertical_diffusion_nml']['ppscheme_type'] = 4")


directories+=("agetracer")
patchlines+=("nml['ocean_tracer_transport_nml']['no_tracer'] = 4; nml['ocean_diagnostics_nml']['diagnose_age'] = True; nml['ocean_diagnostics_nml']['diagnose_green'] = True")


directories+=("layers1")
patchlines+=("nml['ocean_diagnostics_nml']['use_layers'] = True; nml['ocean_diagnostics_nml']['mode_layers'] = 1")


directories+=("layers2")
patchlines+=("nml['ocean_diagnostics_nml']['use_layers'] = True; nml['ocean_diagnostics_nml']['mode_layers'] = 2")
