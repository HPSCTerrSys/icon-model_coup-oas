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
#ifndef RAGNAROK_AES_MICROPHYSICS_COMMON_UTILS_H_
#define RAGNAROK_AES_MICROPHYSICS_COMMON_UTILS_H_

#include <vector>

template <typename T>
using array_2d_t = std::vector<std::vector<T, std::allocator<T>>>;

namespace utils {

template <typename T>
static void calc_dz(const std::vector<T>& z, std::vector<T>& dz, int ncells, int nlev) {
  dz.resize(ncells * nlev);
  array_2d_t<T> zh(nlev + 1, std::vector<T>(ncells));

  for (int i = 0; i < ncells; i++) {
    zh[nlev][i] =
        (static_cast<T>(3.0) * z[i + (nlev - 1) * (ncells)] - z[i + (nlev - 2) * (ncells)]) * static_cast<T>(0.5);
  }

  for (int i = nlev - 1; i >= 0; --i) {
    for (int j = 0; j < ncells; j++) {
      zh[i][j]           = static_cast<T>(2.0) * z[j + (i * ncells)] - zh[i + 1][j];
      dz[i * ncells + j] = -zh[i + 1][j] + zh[i][j];
    }
  }
}
}  // namespace utils

#endif  // RAGNAROK_AES_MICROPHYSICS_COMMON_UTILS_H_
