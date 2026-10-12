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

target_sources(ragnarok PRIVATE thermo.hpp thermo_constants.hpp)

if(BUILD_TESTING)
  set(kokkos_dependents thermo_test.cpp)

  set_source_files_properties(
    ${kokkos_dependents} TARGET_DIRECTORY ragnarok_test
    PROPERTIES LANGUAGE ${Kokkos_COMPILE_LANGUAGE}
  )

  target_sources(ragnarok_test PRIVATE ${kokkos_dependents})
endif()
