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

if(NOT RGK_ENABLE_STANDALONE)
  add_subdirectory(icon_bridge)
endif()

add_subdirectory(timer)
add_subdirectory(testing)
