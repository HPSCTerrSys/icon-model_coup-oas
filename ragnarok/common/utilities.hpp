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
#ifndef RAGNAROK_COMMON_UTILITIES_H_
#define RAGNAROK_COMMON_UTILITIES_H_

#include <Kokkos_Core.hpp>

namespace ragnarok {

template <typename T>
KOKKOS_INLINE_FUNCTION T fmax(const T a, const T b) {
  return Kokkos::fmax(a, b);
}
template <typename T, typename... Args>
KOKKOS_INLINE_FUNCTION T fmax(const T a, const T b, const Args... args) {
  return fmax(a, fmax(b, args...));
}

}  // namespace ragnarok

#endif  // RAGNAROK_COMMON_UTILITIES_H_
