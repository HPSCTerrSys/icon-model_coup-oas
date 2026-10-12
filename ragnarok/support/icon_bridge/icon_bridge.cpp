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

#include "icon_bridge.hpp"

#include <cassert>

#include "icon_domain.hpp"
#include "icon_f2c.hpp"
#include "ragnarok.hpp"

namespace icon {

extern "C" void init_ragnarok_f2c(f2c::FunTable* fun_tab) {
  assert(ragnarok::is_initialized());
  f2c::init(fun_tab);
  init_domain();
}

}  // namespace icon
