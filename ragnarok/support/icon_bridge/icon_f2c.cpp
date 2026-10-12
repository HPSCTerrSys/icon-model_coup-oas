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
/// @brief implements support for Fortran <-> C access
///
//----------------------------

#include "icon_f2c.hpp"

#include <cassert>
#include <sstream>
#include <string>

namespace icon {
namespace f2c {

// init static data of Internal class:
bool Internal::init_state;
FunTable Internal::fun_table;
DomainInfo Internal::dom_info;

void init(const FunTable* funtab) {
  auto& init_state = Internal::init_state;
  if (init_state) return;
  init_state      = true;
  // copy funtab:
  auto& fun_table = Internal::fun_table;
  fun_table       = *funtab;
  // provide dom_info:
  auto& dom_info  = Internal::dom_info;
  fun_table.get_domain_info(&dom_info);
  if (dom_info.id_max < dom_info.id_min) {
    const std::string context = "f2c::init";
    const std::string text    = "no domain available";
    fun_table.finish(context.c_str(), context.size(), text.c_str(), text.size());
  }
}

}  // namespace f2c
}  // namespace icon
