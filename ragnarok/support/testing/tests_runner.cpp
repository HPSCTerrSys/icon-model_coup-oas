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
#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#ifndef __STANDALONE
#include "support/icon_bridge/ragnarok.hpp"
#endif

int main(int argc, char** argv) {
  // Initialize Kokkos before any tests run.
  Kokkos::initialize(argc, argv);

#ifndef __STANDALONE
  ragnarok::initialize_serial_backend();
#endif

  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();

  // Finalize Kokkos after all tests have completed.
  Kokkos::finalize();
  return result;
}
