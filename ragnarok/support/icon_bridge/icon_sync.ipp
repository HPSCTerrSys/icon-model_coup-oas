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

#include "icon_domain.hpp"
#include "icon_exception.hpp"

namespace icon {

template <typename S>
inline const f2c::CommPattern& comm_pat_of_type(const f2c::CommPatch& comm_patch, S typ) {
  if constexpr (std::is_same_v<S, SyncC>) {
    return comm_patch.comm_pat_c;
  } else if constexpr (std::is_same_v<S, SyncE>) {
    return comm_patch.comm_pat_e;
  } else if constexpr (std::is_same_v<S, SyncV>) {
    return comm_patch.comm_pat_v;
  } else if constexpr (std::is_same_v<S, SyncC1>) {
    return comm_patch.comm_pat_c1;
  }
  finish("comm_pat_of_type", "unsupported typ variable");
  return f2c::comm_pat_null;  // never reached
}

template <typename S>
void sync_patch_array_r3(S typ, const Patch& patch, double* arr, int arr_shape[], bool lacc) {
  auto& fun_table    = f2c::get_fun_table();
  auto pat           = comm_pat_of_type(patch.comm_patch, typ);
  auto& my_proc_info = get_process_info();
  if (my_proc_info.is_mpi_parallel) {
    fun_table.exchange_data_r3d(pat, lacc, arr, arr_shape);
  }
}

}  // namespace icon
