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
/// @brief implements top level t_patch initialization
///
//----------------------------

#include <cassert>

#include "icon_domain.hpp"

namespace icon {

namespace {

std::vector<Patch> p_patch;
int p_patch_id_offset = 0;

}  // namespace

void Patch::init(const f2c::PatchDescr f2c_patch_descr) {
  this->f2c_descr            = f2c_patch_descr;
  f2c::PatchInfo& patch_info = *this;
  auto& fun_table            = f2c::get_fun_table();
  fun_table.get_patch_info(f2c_patch_descr, &patch_info);
  fun_table.get_comm_patch(f2c_patch_descr, &this->comm_patch);
}

const Patch& get_patch(int dom_id) {
  int pos = p_patch_id_offset + dom_id;
  assert(pos >= 0 and pos < p_patch.size());
  return p_patch[pos];
}

void init_domain() {
  assert(f2c::is_initialized());
  assert(p_patch.size() == 0);
  auto& fun_table = f2c::get_fun_table();
  auto& dom_info  = f2c::get_domain_info();

  int nids        = dom_info.id_max - dom_info.id_min + 1;
  assert(nids > 0);
  p_patch.resize(nids);
  p_patch_id_offset = -dom_info.id_min;
  for (int id = dom_info.id_min; id <= dom_info.id_max; id++) {
    auto f2c_patch_descr = fun_table.get_mo_model_domain_p_patch_descr(id);
    auto& patch          = p_patch[p_patch_id_offset + id];
    patch.init(f2c_patch_descr);
    assert(patch.id == id);
  }
}

}  // namespace icon
