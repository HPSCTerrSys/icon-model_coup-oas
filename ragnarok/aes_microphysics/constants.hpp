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
#ifndef RAGNAROK_AES_MICROPHYSICS_CONSTANTS_H_
#define RAGNAROK_AES_MICROPHYSICS_CONSTANTS_H_

#include "aes_thermodynamics/thermo_constants.hpp"

namespace graupel {

/// reference air density
template <typename T>
constexpr T rho_00 = T{1.225};
template <typename T>
constexpr T q1 = T{8.e-6};
template <typename T>
constexpr T qmin = T{1.0e-15};  // threshold for computation

/// Formfactor in the mass-size relation of snow particles
template <typename T>
constexpr T ams = T{0.069};

/// Exponent in the mass-size relation of snow particles
template <typename T>
constexpr T bms = T{2.0};

/// prefactor in snow fall speed
template <typename T>
constexpr T v0s = T{25.0};

/// Exponent in the terminal velocity for snow
template <typename T>
constexpr T v1s = T{0.5};

/// Snitial crystal mass for cloud ice nucleation
template <typename T>
constexpr T m0_ice = T{1.0e-12};

template <typename T>
constexpr T tx = T{3339.5};

/// temperature for het. freezing of cloud water with supersat

template <typename T>
constexpr T tfrz_het1 = tmelt<T> - static_cast<T>(6.0);

/// temperature for het. freezing of cloud water
template <typename T>
constexpr T tfrz_het2 = tmelt<T> - static_cast<T>(25.0);

/// temperature for hom. freezing of cloud water
template <typename T>
constexpr T tfrz_hom = tmelt<T> - static_cast<T>(37.0);

}  // namespace graupel

namespace idx {
constexpr short nx  = 6;  // number of water species
constexpr short np  = 4;  // number of precipitating water species
constexpr short lqr = 0;  // index for rain
constexpr short lqi = 1;  // index for ice
constexpr short lqs = 2;  // index for snow
constexpr short lqg = 3;  // index for graupel
constexpr short lqc = 4;  // index for cloud
constexpr short lqv = 5;  // index for vapor
}  // namespace idx

template <typename T>
constexpr T ZERO = T{0.0};

#endif  // RAGNAROK_AES_MICROPHYSICS_CONSTANTS_H_
