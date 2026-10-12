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
#ifndef RAGNAROK_AES_THERMODYNAMICS_THERMO_H_
#define RAGNAROK_AES_THERMODYNAMICS_THERMO_H_

#include <Kokkos_Core.hpp>

#include "common/types.hpp"
#include "thermo_constants.hpp"

namespace thermo {

template <typename T>
KOKKOS_INLINE_FUNCTION T internal_energy(const T tk, const T qv, const T qliq, const T qice, const T rho, const T dz);

template <typename T>
KOKKOS_INLINE_FUNCTION T T_from_internal_energy(const T U, const T qv, const T qliq, const T qice, const T rho,
                                                const T dz);

template <typename T>
KOKKOS_INLINE_FUNCTION T specific_humidity(const T pvapor, const T ptotal);

template <typename T>
KOKKOS_INLINE_FUNCTION T sat_pres_water(const T tk);

template <typename T>
KOKKOS_INLINE_FUNCTION T sat_pres_ice(const T tk);

template <typename T>
KOKKOS_INLINE_FUNCTION T qsat_rho(const T tk, const T rho);

template <typename T>
KOKKOS_INLINE_FUNCTION T qsat_ice_rho(const T tk, const T rho);

template <typename T>
KOKKOS_INLINE_FUNCTION T dqsatdT_rho(const T qs, const T tk);

template <typename T>
KOKKOS_INLINE_FUNCTION T dqsatdT(const T qs, const T tk);

template <typename T>
KOKKOS_INLINE_FUNCTION T dqsatdT_ice(const T qs, const T tk);

template <typename T>
KOKKOS_INLINE_FUNCTION T vaporization_energy(const T tk);

template <typename T>
KOKKOS_INLINE_FUNCTION T sublimation_energy(const T tk);

template <typename T>
KOKKOS_INLINE_FUNCTION T potential_temperature(const T tk, const T pres);

template <typename T>
KOKKOS_INLINE_FUNCTION T dewpoint_temperature(const T tk, const T qv, const T pres);

template <typename T>
KOKKOS_INLINE_FUNCTION void saturation_adjustment(T& te, T& qve, T& qce, const T qre, const T qti, const T rho);

template <typename T>
void saturation_adjustment_2d(const int nproma, const int nlev, const int kstart, const int kend, const int jcs,
                              const int jce, T* ta, T* qv, T* qc, T* qr, T* total_ice, T* rho);

template <typename T>
void saturation_adjustment_2d(const int kstart, const int kend, const int jcs, const int jce, View2D<T> ta_view,
                              View2D<T> qv_view, View2D<T> qc_view, View2D<T> qr_view, View2D<T> total_ice_view,
                              View2D<T> rho_view);

}  // namespace thermo

// include the implementation
#include "thermo.ipp"

#endif  // RAGNAROK_AES_THERMODYNAMICS_THERMO_H_
