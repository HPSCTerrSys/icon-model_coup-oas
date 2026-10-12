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

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_SYNC_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_SYNC_H_

#include "icon_f2c.hpp"

namespace icon {

void sync_patch_array(int typ, f2c::PatchDescr patch, double arr[], int shape[3]);

struct SyncC {
  static constexpr int value = 1;
};

struct SyncE {
  static constexpr int value = 2;
};

struct SyncV {
  static constexpr int value = 3;
};

struct SyncC1 {
  static constexpr int value = 4;
};

constexpr SyncC sync_c;
constexpr SyncE sync_e;
constexpr SyncV sync_v;
constexpr SyncC1 sync_c1;

inline const f2c::ProcessInfo& get_process_info() {
  static const f2c::ProcessInfo process_info = f2c::get_fun_table().get_process_info();
  return process_info;
}

}  // namespace icon

#include "icon_sync.ipp"

#endif  // RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_SYNC_H_
