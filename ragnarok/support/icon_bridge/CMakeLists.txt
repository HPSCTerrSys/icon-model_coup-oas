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

set(kokkos_dependents
    # cmake-format: sort
    icon_bridge.cpp
    icon_domain.cpp
    icon_sync.cpp
    ragnarok.cpp
)

set_source_files_properties(
  ${kokkos_dependents} TARGET_DIRECTORY ragnarok
  PROPERTIES LANGUAGE ${Kokkos_COMPILE_LANGUAGE}
)

target_sources(
  ragnarok
  PRIVATE # cmake-format: sort
          ${kokkos_dependents}
          icon_exception.cpp
          icon_exception.hpp
          icon_f2c.cpp
          icon_f2c.hpp
          icon_sync.hpp
          icon_sync.ipp
          mo_ragnarok.f90
          mo_ragnarok_f2c.f90
          ragnarok.hpp
)
