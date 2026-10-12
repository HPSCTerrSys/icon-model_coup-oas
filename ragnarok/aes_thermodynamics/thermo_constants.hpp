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
#ifndef RAGNAROK_AES_THERMODYNAMICS_CONSTANTS_H_
#define RAGNAROK_AES_THERMODYNAMICS_CONSTANTS_H_

#include "common/physical_constants.hpp"

namespace thermo {
/// constants for saturation vapor pressure
template <typename T>
constexpr T c1es = T{610.78};
template <typename T>
constexpr T c2es = c1es<T> * rd<T> / rv<T>;
template <typename T>
constexpr T c3les = T{17.269};
template <typename T>
constexpr T c3ies = T{21.875};
template <typename T>
constexpr T c4les = T{35.86};
template <typename T>
constexpr T c4ies = T{7.66};
template <typename T>
constexpr T c5les = c3les<T> * (tmelt<T> - c4les<T>);
template <typename T>
constexpr T c5ies = c3ies<T> * (tmelt<T> - c4ies<T>);

/// Specific heat capacity for ice
/// IAPW-06 standard for a reference temperature of the triple point of water (273.16K)
template <typename T>
constexpr T ci = T{2108.0};  // T{2096.8};  //[J/kg/K]
/// invariant part of vaporization enthalpy
template <typename T>
constexpr T lvc = alv<T> - (cpv<T> - clw<T>)*tmelt<T>;
/// invariant part of sublimation enthalpy
template <typename T>
constexpr T lsc = als<T> - (cpv<T> - ci<T>)*tmelt<T>;
}  // namespace thermo

#endif  // RAGNAROK_AES_THERMODYNAMICS_CONSTANTS_H_
