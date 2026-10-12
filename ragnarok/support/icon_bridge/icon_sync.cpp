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
/// @brief This file provides boundary exchange support
///
//----------------------------

#include "icon_sync.hpp"

#include <iostream>

#include "icon_domain.hpp"
#include "icon_exception.hpp"

namespace icon {

// for testing only
extern "C" void ragnarok_sync_patch_array_r3_sync_c(int pid, double arr[], int arr_shape[3]) {
  auto& patch = get_patch(pid);
  sync_patch_array_r3(sync_c, patch, arr, arr_shape, false);
}

}  // namespace icon
