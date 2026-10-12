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

cmake_minimum_required(VERSION 3.18)

if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

project(
  ragnarok
  LANGUAGES CXX Fortran
  VERSION 0.0.1
)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_TESTING "Build tests" ON)
option(RGK_ENABLE_STANDALONE "Enable standalone" ON)
option(RGK_ENABLE_PYTHON_BINDINGS "Enable python bindings" OFF)
option(RGK_ENABLE_SINGLE_PRECISION "Use single precision" OFF)

include(CTest)

if(BUILD_TESTING)
  find_package(GTest) # note that the version constraint won't work without the
                      # CONFIG keyword
  if(NOT GTest_FOUND)
    message(CHECK_START "Fetching external GTest")
    include(FetchContent)
    FetchContent_Declare(
      GTest
      URL https://github.com/google/googletest/releases/download/v1.16.0/googletest-1.16.0.tar.gz
      URL_HASH MD5=9a75eb2ac97300cdb8b65b1a5833f411
    )
    # Function to avoid namespace pollution:
    function(configure_gtest)
      set(BUILD_GMOCK OFF)
      set(INSTALL_GTEST OFF)
      set(gtest_disable_pthreads ON)
      FetchContent_MakeAvailable(GTest)
    endfunction()
    configure_gtest()
    message(CHECK_PASS "done")
  endif()
endif()

# KokkosConfig.cmake overrides the CMAKE_CXX_EXTENSIONS cache variable. We do
# not want our targets (some of which might not even depend on Kokkos) to be
# affected by that. Below, we save the current value of the variable (in case
# the user specified it) and restore it later.
if(DEFINED CMAKE_CXX_EXTENSIONS)
  set(rgk_save_CMAKE_CXX_EXTENSIONS "${CMAKE_CXX_EXTENSIONS}")
else()
  unset(rgk_save_CMAKE_CXX_EXTENSIONS)
endif()

# Suppress the warning from KokkosConfig.cmake:
set(CMAKE_CXX_EXTENSIONS
    OFF
    CACHE BOOL ""
)

# If Kokkos that we find below is built with the CUDA device backend, it will
# try to find the CUDAToolkit package. The correct way to ensure the consistency
# between the CUDA compiler and the CUDAToolkit package is to enable the CUDA
# language first. However, we do not know yet if we need the CUDA language. To
# resolve this circular issue, we enable the language but do not fail the
# configuration if that is not possible. Hopefully, the users will not be
# confused too much by the respective messages when they build against Kokkos
# that is known to have no CUDA support.
include(CheckLanguage)
# The check_language macro ignores CMAKE_CUDA_FLAGS, which might affect the
# result. A way to work around the problem is to temporarily expose the flags
# via the CUDAFLAGS environment variable.
set(rgk_unset_CUDAFLAGS OFF)
if(DEFINED CMAKE_CUDA_FLAGS)
  # If ENV{CUDAFLAGS} is set, it already equals to CMAKE_CUDA_FLAGS and we do
  # not need to do anything. Otherwise, we modify the environment:
  if(NOT DEFINED ENV{CUDAFLAGS})
    set(ENV{CUDAFLAGS} "${CMAKE_CUDA_FLAGS}")
    set(rgk_unset_CUDAFLAGS ON)
  endif()
endif()
check_language(CUDA)
if(rgk_unset_CUDAFLAGS)
  unset(ENV{CUDAFLAGS})
endif()
unset(rgk_unset_CUDAFLAGS)

if(CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
endif()

find_package(Kokkos 4.4.1 CONFIG)

if(NOT Kokkos_FOUND)

  message(CHECK_START "Fetching and configuring Kokkos")

  include(FetchContent)

  FetchContent_Declare(
    kokkos
    URL https://github.com/kokkos/kokkos/releases/download/4.4.01/kokkos-4.4.01.tar.gz
    URL_HASH MD5=eafd0d42c9831858aa84fde78576644c
  )

  set(Kokkos_ENABLE_IMPL_MDSPAN
      OFF
      CACHE BOOL "Experimental mdspan support"
  )

  set(Kokkos_ENABLE_LIBDL
      OFF
      CACHE BOOL "Avoid dependency on libdl"
  )

  set(Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE
      ON
      CACHE BOOL "Whether to use native cmake language support"
  )

  FetchContent_MakeAvailable(kokkos)

  message(CHECK_PASS "done")

  # Kokkos does not expose the variable when configured as a subproject:
  if(NOT DEFINED Kokkos_COMPILE_LANGUAGE)
    set(Kokkos_COMPILE_LANGUAGE "CXX")
    if(Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE)
      if(Kokkos_ENABLE_CUDA)
        set(Kokkos_COMPILE_LANGUAGE "CUDA")
      elseif(Kokkos_ENABLE_HIP)
        set(Kokkos_COMPILE_LANGUAGE "HIP")
      endif()
    endif()
  endif()

endif()

if(DEFINED rgk_save_CMAKE_CXX_EXTENSIONS)
  set(CMAKE_CXX_EXTENSIONS
      "${rgk_save_CMAKE_CXX_EXTENSIONS}"
      CACHE BOOL "" FORCE
  )
  unset(rgk_save_CMAKE_CXX_EXTENSIONS)
else()
  unset(CMAKE_CXX_EXTENSIONS CACHE)
endif()

if(Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE)
  if(Kokkos_ENABLE_CUDA)
    enable_language(CUDA)
  elseif(Kokkos_ENABLE_HIP)
    if(CMAKE_VERSION VERSION_LESS "3.21")
      message(FATAL_ERROR "HIP language requires CMake 3.21 or later")
    endif()
    enable_language(HIP)
  endif()
endif()

# Create a single library comprising all components that ICON can link to:
add_library(ragnarok)

set(Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(
  ragnarok
  PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}"
             Fortran_PREPROCESS ON
)

target_include_directories(
  ragnarok
  PUBLIC
    $<BUILD_INTERFACE:$<$<NOT:$<COMPILE_LANGUAGE:Fortran>>:${CMAKE_CURRENT_SOURCE_DIR}>>
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${Fortran_MODULE_DIRECTORY}>>
)

target_link_libraries(ragnarok PUBLIC Kokkos::kokkos)

if(RGK_ENABLE_SINGLE_PRECISION)
  target_compile_definitions(ragnarok PRIVATE __SINGLE_PRECISION)
endif()

if(RGK_ENABLE_STANDALONE)
  target_compile_definitions(ragnarok PUBLIC __STANDALONE)

  list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
  find_package(
    NetCDF
    COMPONENTS C CXX
    REQUIRED
  )
  include_directories(${NetCDF_CXX_INCLUDE_DIR})
endif()

if(RGK_ENABLE_PYTHON_BINDINGS)
  # Function to avoid namespace pollution:
  function(find_pybind11)
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 REQUIRED)
  endfunction()
  find_pybind11()
endif()

if(BUILD_TESTING)

  add_executable(ragnarok_test)

  target_link_libraries(ragnarok_test PRIVATE ragnarok GTest::gtest)

  # Delay test discovery until runtime (needed for cross-compilation):
  gtest_discover_tests(ragnarok_test DISCOVERY_MODE PRE_TEST)
else()
  # Allow for 'make test' even if the tests are disabled:
  enable_testing()
endif()

add_subdirectory(common)
add_subdirectory(support)
add_subdirectory(aes_thermodynamics)
add_subdirectory(aes_microphysics)
