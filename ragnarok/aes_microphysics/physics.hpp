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
#ifndef RAGNAROK_AES_MICROPHYSICS_PHYSICS_H_
#define RAGNAROK_AES_MICROPHYSICS_PHYSICS_H_

#include "Kokkos_Core.hpp"
#include "aes_microphysics/constants.hpp"

namespace physics {

///
/// @param [in] qi Ice specific mass
/// @param [in] m_ice Ice crystal mass
/// @param [in] ice_dep Rate of ice deposition (some to snow)
///
template <typename T>
KOKKOS_INLINE_FUNCTION T deposition_auto_conversion(const T qi, const T m_ice, const T ice_dep) {
  constexpr T m0_s  = T{3.0e-9};  // initial mass of snow crystals
  constexpr T b_dep = static_cast<T>(2.0) / static_cast<T>(3.0);
  constexpr T xcrit = T{1.0};  // threshold parameter

  T result          = ZERO<T>;

  if (qi > graupel::qmin<T>) {
    T tau_inv = b_dep / (Kokkos::pow((m0_s / m_ice), b_dep) - xcrit);
    result    = Kokkos::fmax(ZERO<T>, ice_dep) * tau_inv;
  }

  return result;
}

///
/// @param [in] t Temperature
/// @param [in] qvsi Saturation (ice) specific vapor mass
/// @return Deposition factor
///
template <typename T>
KOKKOS_INLINE_FUNCTION T deposition_factor(const T t, const T qvsi) {
  constexpr T kappa = T{2.40e-2};  // thermal conductivity of dry air
  constexpr T b     = T{1.94};
  constexpr T a     = als<T> * als<T> / (kappa * rv<T>);

  T cx              = static_cast<T>(2.22e-5) * Kokkos::pow(tmelt<T>, (-b)) * static_cast<T>(101325.0);
  T x               = cx / rd<T> * Kokkos::pow(t, b - static_cast<T>(1.0));
  return x / (static_cast<T>(1.0) + a * x * qvsi / (t * t));
}

///
/// @brief Ice deposition nucleation
/// @param [in] t Temperature
/// @param [in] qc Specific humidity of ice
/// @param [in] qi Specific humidity of ice
/// @param [in] ni Ice crystal number
/// @param [in] dvsi Vapor excess with respect to ice sat
/// @param [in] dt Time step
/// @return Rate of vapor deposition for new ice
///
template <typename T>
KOKKOS_INLINE_FUNCTION T ice_deposition_nucleation(const T t, const T qc, const T qi, const T ni, const T dvsi,
                                                   const T dt) {
  return (qi <= graupel::qmin<T> &&
          ((t < graupel::tfrz_het2<T> && dvsi > ZERO<T>) || (t <= graupel::tfrz_het1<T> && qc > graupel::qmin<T>)))
             ? Kokkos::fmin(graupel::m0_ice<T> * ni, Kokkos::fmax(ZERO<T>, dvsi)) / dt
             : ZERO<T>;
}

///
/// @param [in] qi Ice specific mass
/// @param [in] ni Ice crystal number
/// @return ice mass
///
template <typename T>
KOKKOS_INLINE_FUNCTION T ice_mass(const T qi, const T ni) {
  constexpr T mi_max = T{1.0e-09};  // maximum mass of cloud ice crystals
  return Kokkos::fmax(graupel::m0_ice<T>, Kokkos::fmin(qi / ni, mi_max));
}

///
/// @brief Ice number following cooper
/// @param [in] t Ambient temperature (kelvin)
/// @param [in] rho Ambient density
/// @return Ice number
///
template <typename T>
KOKKOS_INLINE_FUNCTION T ice_number(const T t, const T rho) {
  constexpr T a     = T{5.000};    // parameter in cooper fit
  constexpr T b     = T{0.304};    // parameter in cooper fit
  constexpr T nimax = T{250.e+3};  // maximal number of ice crystals
  return Kokkos::fmin(nimax, a * Kokkos::exp(b * (tmelt<T> - t))) / rho;
}

///
/// @brief sticking efficiency of ice
/// @param [in] t Temperature
/// @return Ice sticking
//
template <typename T>
KOKKOS_INLINE_FUNCTION T ice_sticking(const T t, const T cia) {
  constexpr T a       = T{0.09};                         // scale factor for freezing depression
  constexpr T b       = T{1.00};                         // maximum for exponential temperature factor
  constexpr T eff_min = T{0.075};                        // minimum sticking efficiency
  constexpr T eff_fac = T{3.5E-3};                       // Scaling factor [1/K] for cloud ice sticking efficiency
  constexpr T tcrit   = tmelt<T> - static_cast<T>(85.);  //   Temperature at which cloud ice autoconversion starts

  // per original code seems like aggregation is allowed even with no snow
  // present
  return Kokkos::fmax(Kokkos::fmax(Kokkos::fmin(Kokkos::exp(a * (t - tmelt<T>)), b), eff_min), eff_fac * (t - tcrit)) *
         cia;
}

///
/// @param [in] rho_s Snow specific density
/// @param [in] ns Snow number
/// @return riming snow rate
///
template <typename T>
KOKKOS_INLINE_FUNCTION T snow_lambda(const T rho_s, const T ns) {
  constexpr T lmd_0 = T{1.0e+10};  // no snow value of lambda

  return (rho_s > graupel::qmin<T>) ? Kokkos::pow((static_cast<T>(2.) * graupel::ams<T> * ns / rho_s),
                                                  (static_cast<T>(1.) / (graupel::bms<T> + static_cast<T>(1.))))
                                    : lmd_0;
}

/// Snow number concentration as a function of temperature and snow specific density
/// @param [in] t Temperature
/// @param [in] rho_s Snow specific density
/// @return Snow number
///
template <typename T>
KOKKOS_INLINE_FUNCTION T snow_number(const T t, const T rho_s) {
  constexpr T tmin     = tmelt<T> - static_cast<T>(40.);
  constexpr T tmax     = tmelt<T>;
  constexpr T rho_s_mn = T{2.0e-7};
  constexpr T xa1      = T{-1.65e+0};
  constexpr T xa2      = T{5.45e-2};
  constexpr T xa3      = T{3.27e-4};
  constexpr T xb1      = T{1.42e+0};
  constexpr T xb2      = T{1.19e-2};
  constexpr T xb3      = T{9.60e-5};
  constexpr T n0s0     = T{8.00e+5};
  constexpr T n0s1     = static_cast<T>(13.5) * static_cast<T>(5.65e+05);
  constexpr T n0s2     = T{-0.107};
  constexpr T n0s3     = T{13.5};
  constexpr T n0s4     = static_cast<T>(0.5) * n0s1;
  constexpr T n0s5     = T{1.e6};
  constexpr T n0s6     = static_cast<T>(1.e2) * n0s1;
  constexpr T n0s7     = T{1.e9};

  if (rho_s > graupel::qmin<T>) {
    T tc  = Kokkos::fmax(Kokkos::fmin(t, tmax), tmin) - tmelt<T>;
    T alf = Kokkos::pow(static_cast<T>(10.), (xa1 + tc * (xa2 + tc * xa3)));
    T bet = xb1 + tc * (xb2 + tc * xb3);
    T n0s = n0s3 *
            Kokkos::pow((Kokkos::fmax(rho_s, rho_s_mn) / graupel::ams<T>),
                        (static_cast<T>(4.0) - static_cast<T>(3) * bet)) /
            (alf * alf * alf);
    T y     = Kokkos::exp(n0s2 * tc);
    T n0smn = Kokkos::fmax(n0s4 * y, n0s5);
    T n0smx = Kokkos::fmin(n0s6 * y, n0s7);
    return Kokkos::fmin(n0smx, Kokkos::fmax(n0smn, n0s));
  } else {
    return n0s0;
  }
}

///
/// @param [in] iqx Index of hydrometeor
/// @param [in] rho_x Hydrometeor density
/// @param [in] rho Air density
/// @param [in] t Temperature
///
template <typename T>
KOKKOS_INLINE_FUNCTION T vm(int iqx, const T rho_x, const T rho, const T t) {
  const T rho_mx = T{6.97604e-03};
  const T rho_mn = T{3.26216e-08};

  T x            = Kokkos::fmin(rho_mx, Kokkos::fmax(rho_mn, rho_x));

  switch (iqx) {
    case idx::lqr: {
      const T a_r[5] = {T{-5.91051e-01}, T{-5.37440e+00}, T{-1.00459e+00}, T{-6.44895e-02}, T{-1.40361e-03}};
      x              = Kokkos::log(x);
      return a_r[0] + x * (a_r[1] + x * (a_r[2] + x * (a_r[3] + x * a_r[4]))) * Kokkos::sqrt(graupel::rho_00<T> / rho);
    }
    case idx::lqi: {
      const T a_i[2]  = {T{0.80}, T{0.160}};
      constexpr T b_i = static_cast<T>(1.0) / static_cast<T>(3.0);
      return a_i[0] * Kokkos::pow(x, a_i[1]) * Kokkos::pow(graupel::rho_00<T> / rho, b_i);
    }
    case idx::lqs: {
      const T a_s[2]  = {T{2.0} * T{57.80}, T{1.0} / T{6.0}};
      constexpr T b_s = -static_cast<T>(1.0) / static_cast<T>(6.0);
      return a_s[0] * Kokkos::pow(x, a_s[1]) * Kokkos::sqrt(graupel::rho_00<T> / rho) *
             Kokkos::pow(snow_number(t, x), b_s);
    }
    case idx::lqg: {
      const T a_g[2] = {T{12.24}, T{0.217}};
      return a_g[0] * Kokkos::pow(x, a_g[1]) * Kokkos::sqrt(graupel::rho_00<T> / rho);
    }
    default:
      return ZERO<T>;
  }
}

}  // namespace physics

#endif  // RAGNAROK_AES_MICROPHYSICS_PHYSICS_H_
