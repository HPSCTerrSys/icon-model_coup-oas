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
/// @brief This header file contains the definitions required
/// for supporting Fortran <=> C access.
///
//----------------------------

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_F2C_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_F2C_H_

#include <cassert>
#include <cstddef>
namespace icon {
namespace f2c {

// opaque pointer to a Fortran t_patch instance
struct PatchDescr {
  void* cptr = NULL;
};

// portable part of top level t_patch components
struct PatchInfo {
  int id      = -1;
  int nlev    = 0;
  int nblks_c = 0;
  int nblks_e = 0;
  int nblks_v = 0;
  // to be extended
};

// domain info:
struct DomainInfo {
  int nproma = 0;
  int id_min = -1;
  int id_max = -2;
};

// opaque pointer to a Fortran ICON CommPattern instance
struct CommPatternDescr {
  void* cptr = NULL;
};

constexpr CommPatternDescr comm_pat_null = {NULL};

using CommPattern                        = CommPatternDescr;

// communication patch
struct CommPatch {
  CommPattern comm_pat_c;
  CommPattern comm_pat_e;
  CommPattern comm_pat_v;
  CommPattern comm_pat_c1;
};

struct ProcessInfo {
  bool is_mpi_parallel;
};

struct FunTable {
  void (*get_domain_info)(DomainInfo* dom_info)                                               = NULL;
  void (*get_patch_info)(const PatchDescr patch_descr, PatchInfo* patch_info)                 = NULL;
  PatchDescr (*get_mo_model_domain_p_patch_descr)(int id)                                     = NULL;
  void (*get_comm_patch)(const PatchDescr patch_descr, CommPatch* comm_patch)                 = NULL;
  void (*message)(const char* name, int name_len, const char* text, int text_len)             = NULL;
  void (*finish)(const char* name, int name_len, const char* text, int text_len)              = NULL;
  void (*exchange_data_r3d)(CommPatternDescr descr, bool lacc, double* recv, int* recv_shape) = NULL;
  ProcessInfo (*get_process_info)()                                                           = NULL;
  int (*new_timer)(const char* name, int name_len)                                            = NULL;
  void (*timer_start)(int timer)                                                              = NULL;
  void (*timer_stop)(int timer)                                                               = NULL;
  double (*timer_value)(int timer)                                                            = NULL;
  void (*get_timer_config)(bool& ltimer, int& timers_level)                                   = NULL;
};

class Internal {
  friend void init(const FunTable* fun_tab);
  friend bool is_initialized();
  friend const FunTable& get_fun_table();
  friend const DomainInfo& get_domain_info();
  static bool init_state;
  static FunTable fun_table;
  static DomainInfo dom_info;
};

inline bool is_initialized() { return Internal::init_state; }

inline const FunTable& get_fun_table() {
  assert(Internal::init_state);
  return Internal::fun_table;
}

inline const DomainInfo& get_domain_info() {
  assert(Internal::init_state);
  return Internal::dom_info;
}

void init(const FunTable* fun_tab);

}  // namespace f2c
}  // namespace icon

#endif  // RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_F2C_H_
