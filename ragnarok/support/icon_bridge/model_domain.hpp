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

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_DOMAIN_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_DOMAIN_H_

#include <Kokkos_Core.hpp>

#include "icon_bridge.hpp"
#include "icon_f2c.hpp"

namespace icon {

struct Patch : public f2c::PatchInfo {
  f2c::CommPatch comm_patch;
  void init(const f2c::PatchDescr f2c_patch_descr);

 private:
  f2c::PatchDescr f2c_descr;
};

void init_domain();

const Patch& get_patch(int dom_id);

}  // namespace icon

#endif  // RAGNAROK_SUPPORT_ICON_BRIDGE_MODEL_ICON_DOMAIN_H_
