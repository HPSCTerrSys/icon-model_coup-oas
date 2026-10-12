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

set(kokkos_dependents bindings/fortran/aes_muphys.cpp)
set_source_files_properties(
  ${kokkos_dependents} TARGET_DIRECTORY ragnarok
  PROPERTIES LANGUAGE ${Kokkos_COMPILE_LANGUAGE}
)

target_sources(
  ragnarok
  PRIVATE graupel.hpp
          graupel.ipp
          common/utils.hpp
          constants.hpp
          physics.hpp
          transition.hpp
          # always enable Fortran bindings
          bindings/fortran/aes_muphys.hpp
          bindings/fortran/mo_ragnarok_microphysics.f90
          ${kokkos_dependents}
)

if(RGK_ENABLE_PYTHON_BINDINGS)
  pybind11_add_module(aes_muphys_py bindings/python/aes_muphys_py.cpp)
  target_link_libraries(aes_muphys_py PRIVATE ragnarok)
  set_target_properties(aes_muphys_py PROPERTIES OUTPUT_NAME "aes_muphys_py")
  set_target_properties(aes_muphys_py PROPERTIES SUFFIX ".so")
  set_target_properties(ragnarok PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

if(RGK_ENABLE_STANDALONE)
  add_subdirectory(standalone)
endif()

if(BUILD_TESTING)
  set(kokkos_dependents_test
      physics_test.cpp transition_test.cpp graupel_test.cpp
  )

  set_source_files_properties(
    ${kokkos_dependents_test} TARGET_DIRECTORY ragnarok_test
    PROPERTIES LANGUAGE ${Kokkos_COMPILE_LANGUAGE}
  )

  target_sources(ragnarok_test PRIVATE ${kokkos_dependents_test})
endif()
