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
add_executable(aes_graupel aes_graupel.cpp)

set_source_files_properties(
  aes_graupel.cpp PROPERTIES LANGUAGE ${Kokkos_COMPILE_LANGUAGE}
)

target_include_directories(
  aes_graupel
  PUBLIC
    $<BUILD_INTERFACE:$<$<NOT:$<COMPILE_LANGUAGE:Fortran>>:${PROJECT_SOURCE_DIR}>>
)

target_link_libraries(
  aes_graupel PRIVATE Kokkos::kokkos NetCDF::NetCDF_CXX NetCDF::NetCDF_C
)
