// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------

///
/// @file
/// @brief This file contains the implementation for
/// accessing some of ICON's Fortran functionality
///
//----------------------------

#include "ragnarok.hpp"

#include <string.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <cassert>
#include <iostream>

namespace ragnarok {

enum class InitState { pre_init = 0, initialized, finalized, problem };

InitState init_state = InitState::pre_init;

bool is_initialized() { return init_state == InitState::initialized; }

namespace {

std::string construct_kokkos_version() {
  return std::to_string(KOKKOS_VERSION_MAJOR) + "." + std::to_string(KOKKOS_VERSION_MINOR) + "." +
         std::to_string(KOKKOS_VERSION_PATCH);
}

}  // namespace

void init() {
  assert(init_state == InitState::pre_init);
  Kokkos::initialize();
  if (ragnarok::isCPU()) {
    initialize_serial_backend();
  }
  init_state = InitState::initialized;
}

void finalize() {
  assert(init_state == InitState::initialized);
  serial_spaces.clear();
  Kokkos::finalize();
  // enter final state - we are not even allowed to call Kokkos::initialize again:
  init_state = InitState::finalized;
}

const std::string& retrieve_kokkos_version() {
  const static std::string version = construct_kokkos_version();
  return version;
}

}  // namespace ragnarok

extern "C" {

void init_ragnarok() { ragnarok::init(); }
void finalize_ragnarok() { ragnarok::finalize(); }

const char* retrieve_kokkos_version_c(int& length) {
  const std::string& version = ragnarok::retrieve_kokkos_version();
  length                     = static_cast<int>(version.length());
  return version.c_str();
}

}  // extern "C"
