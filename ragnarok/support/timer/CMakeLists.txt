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

if(RGK_ENABLE_STANDALONE)
  target_sources(ragnarok PRIVATE rgk_timer.hpp)
else()
  target_sources(
    ragnarok PRIVATE icon_timer.hpp icon_timer.cpp mo_ragnarok_timer.f90
  )
endif()

if(BUILD_TESTING AND RGK_ENABLE_STANDALONE)
  set(kokkos_dependents timer_test.cpp)

  set_source_files_properties(
    ${kokkos_dependents} TARGET_DIRECTORY ragnarok_test
    PROPERTIES LANGUAGE ${Kokkos_COMPILE_LANGUAGE}
  )

  target_sources(ragnarok_test PRIVATE ${kokkos_dependents})
endif()
