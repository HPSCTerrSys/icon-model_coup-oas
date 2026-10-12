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
#ifndef RAGNAROK_AES_MICROPHYSICS_TRANSITION_H_
#define RAGNAROK_AES_MICROPHYSICS_TRANSITION_H_

#include "Kokkos_Core.hpp"
#include "aes_microphysics/constants.hpp"

namespace transition {

///
/// @brief Accretion
/// @param [in] t Temperature
/// @param [in] rho Ambient density
/// @param [in] qc Snow specific mass
/// @param [in] qg Graupel specific mass
/// @return Graupel riming rate
///
template <typename T>
KOKKOS_INLINE_FUNCTION T cloud_to_graupel(const T t, const T rho, const T qc, const T qg) {
  constexpr T a_rim = T{4.43};
  constexpr T b_rim = T{0.94878};

  return (Kokkos::fmin(qc, qg) > graupel::qmin<T> && t > graupel::tfrz_hom<T>)
             ? a_rim * qc * Kokkos::pow(qg * rho, b_rim)
             : ZERO<T>;
}

/// @brief Conversion from cloud water to rain water
/// @details
/// Kessler (1969) autoconversion rate
///    scau = zccau * MAX( qc_ik - qc0, 0.0 )
///    scac = zcac  * qc_ik * zeln7o8qrk
/// Seifert and Beheng (2001) autoconversion rate
/// with constant cloud droplet number concentration qnc
///
/// @param [in] t Temperature
/// @param [in] rho Cloud water specific density
/// @param [in] qc Cloud water specific mass
/// @param [in] qr Rain water specific mass
/// @param [in] nc Cloud water number concentration
/// @return conversion rate
///
template <typename T>
KOKKOS_INLINE_FUNCTION T cloud_to_rain(const T t, const T rho, const T qc, const T qr, const T nc) {
  constexpr T qmin_ac = T{1.00e-06};  // threshold for auto conversion
  constexpr T tau_max = T{0.90e+00};  // maximum allowed value of tau
  constexpr T tau_min = T{1.00e-30};  // maximum allowed value of tau
  constexpr T a_phi   = T{6.00e+02};  // constant in phi-function for autoconversion
  constexpr T b_phi   = T{0.68e+00};  // exponent in phi-function for autoconversion
  constexpr T c_phi   = T{5.00e-05};  // exponent in phi-function for accretion
  constexpr T x3      = T{2.00e+00};  // gamma exponent for cloud distribution
  constexpr T x2      = T{2.60e-10};  // separating mass between cloud and rain
  constexpr T x1      = T{9.44e+09};  // kernel coeff for SB2001 autoconversion
  constexpr T rho_mx  = T{6.97604e-03};
  constexpr T rho_mn  = T{3.26216e-08};

  const T au_kernel   = x1 / (static_cast<T>(20.0) * x2) * (x3 + static_cast<T>(2.0)) * (x3 + static_cast<T>(4.0)) /
                      Kokkos::pow((x3 + static_cast<T>(1.0)), static_cast<T>(2.0));
  const T a_ac[5] = {T{-2.155543e+00}, T{-1.148491e+00}, T{-1.882563e-02}, T{2.941391e-03}, T{5.575598e-05}};

  T result        = ZERO<T>;

  if (qc > qmin_ac && t > graupel::tfrz_hom<T>) {
    T x         = Kokkos::log(Kokkos::fmin(rho_mx, Kokkos::fmax(rho_mn, rho * qr)));
    T ac_kernel = a_ac[0] + x * (a_ac[1] + x * (a_ac[2] + x * (a_ac[3] + x * a_ac[4])));
    T tau       = Kokkos::fmax(tau_min, Kokkos::fmin(static_cast<T>(1.0) - qc / (qc + qr),
                                                     tau_max));  // time-scale
    T phi       = Kokkos::pow(tau, b_phi);                       // similarity function for autoconversion
    phi         = a_phi * phi * Kokkos::pow((static_cast<T>(1.0) - phi), static_cast<T>(3.0));
    T xau       = au_kernel * Kokkos::pow(qc * qc / nc, static_cast<T>(2.)) *
            (static_cast<T>(1.0) + phi / Kokkos::pow(static_cast<T>(1.0) - tau,
                                                     static_cast<T>(2.0)));                  // autoconversion rate
    T xac  = ac_kernel * qc * qr * Kokkos::pow((tau / (tau + c_phi)), static_cast<T>(4.0));  // accretion rate
    result = xau + xac;
  }

  return result;
}

///
/// @param [in] t Temperature
/// @param [in] qc Cloud specific mass
/// @param [in] qs Snow specific mass
/// @param [in] ns Snow number
/// @param [in] lambda   Snow slope parameter (lambda)
/// @return Riming snow rate
///
template <typename T>
KOKKOS_INLINE_FUNCTION T cloud_to_snow(const T t, const T qc, const T qs, const T ns, const T lambda) {
  /// Collection efficiency for snow collecting cloud water
  constexpr T ecs   = T{0.9};
  constexpr T b_rim = -(graupel::v1s<T> + static_cast<T>(3.0));
  // (with pi*gam(v1s+3)/4 = 2.610) and tunning factor 3
  constexpr T c_rim = static_cast<T>(2.61) * ecs * graupel::v0s<T> * static_cast<T>(3.);

  return (Kokkos::fmin(qc, qs) > graupel::qmin<T> && t > graupel::tfrz_hom<T>)
             ? (c_rim * ns) * qc * Kokkos::pow(lambda, b_rim)
             : ZERO<T>;
}

///
/// @param [in] t Temperature
/// @param [in] qc Cloud specific mass
/// @param [in] qi Ice specific mass
/// @param [in] dt Time step
/// @return Homogeneous freezing rate
///
template <typename T>
KOKKOS_INLINE_FUNCTION T cloud_x_ice(const T t, const T qc, const T qi, const T dt) {
  T result = ZERO<T>;

  if (qc > graupel::qmin<T> && t < graupel::tfrz_hom<T>) result = qc / dt;

  if (qi > graupel::qmin<T> && t > tmelt<T>) result = -qi / dt;

  return result;
}

///
/// @brief Melting of graupel to form rain
/// @param [in] t Ambient temperature
/// @param [in] p Ambient pressure
/// @param [in] rho Ambient density
/// @param [in] dvsw0 qv-qsat_water(T0)
/// @param [in] qg graupel specific mass
///
template <typename T>
KOKKOS_INLINE_FUNCTION T graupel_to_rain(const T t, const T p, const T rho, const T dvsw0, const T qg) {
  constexpr T c1_melt = T{12.31698};
  constexpr T c2_melt = T{7.39441e-05};
  constexpr T a_melt  = graupel::tx<T> - static_cast<T>(389.5);     // melting prefactor
  constexpr T b_melt  = static_cast<T>(3.0) / static_cast<T>(5.0);  // melting exponent

  return (t > Kokkos::fmax(tmelt<T>, tmelt<T> - graupel::tx<T> * dvsw0) && qg > graupel::qmin<T>)
             ? (c1_melt / p + c2_melt) * (t - tmelt<T> + a_melt * dvsw0) * Kokkos::pow(qg * rho, b_melt)
             : ZERO<T>;
}

///
/// @param [in] rho Density
/// @param [in] qr Rain specific mass
/// @param [in] qg Graupel specific mass
/// @param [in] qi Ice specific mass
/// @param [in] sticking_eff Sticking effiency
/// @return aggregation of ice by graupel
///
template <typename T>
KOKKOS_INLINE_FUNCTION T ice_to_graupel(const T rho, const T qr, const T qg, const T qi, const T sticking_eff) {
  constexpr T a     = T{1.72};  //  (15/32)*(PI**0.5)*(EIR/RHOW)*V0R*AR**(1/8)
  constexpr T b     = static_cast<T>(7.0) / static_cast<T>(8.0);
  constexpr T c_agg = T{2.46};
  constexpr T b_agg = T{0.94878};

  T result          = ZERO<T>;
  if (qi > graupel::qmin<T>) {
    if (qg > graupel::qmin<T>) {
      result = sticking_eff * qi * c_agg * Kokkos::pow(rho * qg, b_agg);
    }
    if (qr > graupel::qmin<T>) {
      result = result + a * qi * Kokkos::pow(rho * qr, b);
    }
  }

  return result;
}

///
/// @brief Conversion rate of ice to snow
/// @param [in] qi Ice specific mass
/// @param [in] ns Snow number
/// @param [in] lambda Snow intercept parameter, lambda
/// @param [in] sticking_eff Ice sticking effiency
/// @return conversion rate of ice to snow
///
template <typename T>
KOKKOS_INLINE_FUNCTION T ice_to_snow(const T qi, const T ns, const T lambda, const T sticking_eff) {
  constexpr T qi0   = ZERO<T>;                                   // critical ice required for autoconversion
  constexpr T c_iau = T{1.0E-3};                                 // coefficient of auto conversion
  constexpr T c_agg = static_cast<T>(2.61) * graupel::v0s<T>;    // coeff of aggregation (2.610 = pi*gam(v1s+3)/4)
  constexpr T b_agg = -(graupel::v1s<T> + static_cast<T>(3.0));  // aggregation exponent

  return (qi > graupel::qmin<T>) ? sticking_eff * (c_iau * Kokkos::fmax(ZERO<T>, (qi - qi0)) +
                                                   qi * (c_agg * ns) * Kokkos::pow(lambda, b_agg))
                                 : ZERO<T>;
}

///
/// @brief Freezing rain
/// @param [in] t Temperature
/// @param [in] rho Ambient density
/// @param [in] qr Specific humidity of rain
/// @param [in] qc Cloud liquid specific mass
/// @param [in] qi Cloud ice specific mass
/// @param [in] qs Snow specific mass
/// @param [in] mi Ice crystal mass
/// @param [in] dvsw qv-qsat_water(T)
/// @param [in] dt Time step
/// @return convertion rate from graupel to rain
///
template <typename T>
KOKKOS_INLINE_FUNCTION T rain_to_graupel(const T t, const T rho, const T qc, const T qr, const T qi, const T qs,
                                         const T mi, const T dvsw, const T dt) {
  constexpr T tfrz_rain = tmelt<T> - static_cast<T>(2.0);
  constexpr T a1        = T{9.95e-5};  // FR: 1. coefficient for immersion raindrop freezing: alpha_if
  constexpr T b1 =
      static_cast<T>(7.0) / static_cast<T>(4.0);  // FR: 2. coefficient for immersion raindrop freezing: a_if
  constexpr T c2      = T{0.66};                  // FR: 2. coefficient for immersion raindrop freezing: a_if
  constexpr T c3      = T{1.0};                   // FR: 2. coefficient for immersion raindrop freezing: a_if
  constexpr T c4      = T{0.1};                   // FR: 2. coefficient for immersion raindrop freezing: a_if
  constexpr T a2      = T{1.24E-3};               //  (PI/24)*EIR*V0R*Gamma(6.5)*AR**(-5/8)
  constexpr T b2      = static_cast<T>(13.0) / static_cast<T>(8.0);
  constexpr T qs_crit = T{1.e-7};

  T result            = ZERO<T>;

  if (qr > graupel::qmin<T> && t < tfrz_rain) {
    if (t > graupel::tfrz_hom<T>) {
      if (dvsw + qc <= ZERO<T> || qr > c4 * qc) {
        result = (Kokkos::exp(c2 * (tfrz_rain - t)) - c3) * (a1 * Kokkos::pow((qr * rho), b1));
      }
    } else {
      result = qr / dt;
    }
  }

  if (fmin(qi, qr) > graupel::qmin<T> && qs > qs_crit) {
    //  rain + ice creating graupel
    result = result + a2 * (qi / mi) * Kokkos::pow((rho * qr), b2);
  }

  return result;
}

///
/// @param [in] t Temperature
/// @param [in] rho Ambient density
/// @param [in] qc Specific humidity of cloud
/// @param [in] qr Specific humidity of rain
/// @param [in] dvsw qv - qsat_water(T)
/// @param [in] dt Time step
/// @return Mass from qc to qr
///
template <typename T>
KOKKOS_INLINE_FUNCTION T rain_to_vapor(const T t, const T rho, const T qc, const T qr, const T dvsw, const T dt) {
  constexpr T c1     = T{0.61};
  constexpr T c2     = T{-0.0163};
  constexpr T c3     = T{1.111e-4};
  constexpr T rho_mx = T{6.97604e-03};
  constexpr T rho_mn = T{3.26216e-08};
  const T a_ev[5]    = {T{-5.532194e+00}, T{2.432848e-01}, T{-4.145391e-02}, T{-1.798439e-03}, T{-1.405764e-05}};

  if (qr > graupel::qmin<T> && (dvsw + qc <= ZERO<T>)) {
    T tc       = t - tmelt<T>;
    T evap_max = (c1 + tc * (c2 + c3 * tc)) * (-dvsw) / dt;
    T x        = Kokkos::log(Kokkos::fmin(rho_mx, Kokkos::fmax(rho_mn, qr * rho)));
    T evap     = -Kokkos::exp(a_ev[0] + x * (a_ev[1] + x * (a_ev[2] + x * (a_ev[3] + x * a_ev[4])))) * dvsw;
    return Kokkos::fmin(evap, evap_max);
  }

  return ZERO<T>;
}

///
/// @param [in] t ambient temperature
/// @param [in] rho ambient density
/// @param [in] qc cloud specific mass
/// @param [in] qs snow specific mass
/// @returns convertion rate
///
template <typename T>
KOKKOS_INLINE_FUNCTION T snow_to_graupel(const T t, const T rho, const T qc, const T qs) {
  /// Constants in riming formula
  constexpr T a_rim = T{0.5};
  constexpr T b_rim = static_cast<T>(3.0) / static_cast<T>(4.0);

  return (Kokkos::fmin(qc, qs) > graupel::qmin<T> && t > graupel::tfrz_hom<T>)
             ? a_rim * qc * Kokkos::pow(qs * rho, b_rim)
             : ZERO<T>;
}

///
/// @brief Melting of snow to form rain
/// @param [in] t Temperature
/// @param [in] p Ambient pressure
/// @param [in] rho Ambient density
/// @param [in] dvsw0  qv-qsat_water(T0)
/// @param [in] qs Snow specific mass
/// @return conversion rate from snow to rain
///
template <typename T>
KOKKOS_INLINE_FUNCTION T snow_to_rain(const T t, const T p, const T rho, const T dvsw0, const T qs) {
  constexpr T c1 = T{79.6863};                                 // Constants in melting formula
  constexpr T c2 = T{0.612654E-3};                             // Constants in melting formula
  constexpr T a  = graupel::tx<T> - static_cast<T>(389.5);     // melting prefactor
  constexpr T b  = static_cast<T>(4.0) / static_cast<T>(5.0);  // melting exponent

  return (t > Kokkos::fmax(tmelt<T>, tmelt<T> - graupel::tx<T> * dvsw0) && qs > graupel::qmin<T>)
             ? (c1 / p + c2) * (t - tmelt<T> + a * dvsw0) * Kokkos::pow(qs * rho, b)
             : ZERO<T>;
}

///
/// @brief Graupel-vapor exchange rate
/// @param [in] t Temperature
/// @param [in] p Ambient pressure
/// @param [in] rho Ambient density
/// @param [in] qg Graupel specific mass
/// @param [in] dvsw qv-qsat_water(T)
/// @param [in] dvsi qv-qsat_ice(T)
/// @param [in] dvsw0 qv-qsat_water(T0)
/// @param [in] dt Time step
///
template <typename T>
KOKKOS_INLINE_FUNCTION T vapor_x_graupel(const T t, const T p, const T rho, const T qg, const T dvsw, const T dvsi,
                                         const T dvsw0, const T dt) {
  constexpr T a1 = T{0.398561};
  constexpr T a2 = T{-0.00152398};
  constexpr T a3 = T{2554.99};
  constexpr T a4 = T{2.6531E-7};
  constexpr T a5 = T{0.153907};
  constexpr T a6 = T{-7.86703e-07};
  constexpr T a7 = T{0.0418521};
  constexpr T a8 = T{-4.7524E-8};
  constexpr T b  = T{0.6};

  T result       = ZERO<T>;

  if (qg > graupel::qmin<T>) {
    if (t < tmelt<T>) {
      result = (a1 + a2 * t + a3 / p + a4 * p) * dvsi * Kokkos::pow(qg * rho, b);
    } else {
      if (t > (tmelt<T> - graupel::tx<T> * dvsw0)) {
        result = (a5 + a6 * p) * Kokkos::fmin(ZERO<T>, dvsw0) * Kokkos::pow(qg * rho, b);
      } else {
        result = (a7 + a8 * p) * dvsw * Kokkos::pow(qg * rho, b);
      }
    }
    result = Kokkos::fmax(result, -qg / dt);
  }

  return result;
}

///
/// @param [in] qi Specific humidity of ice
/// @param [in] mi Ice crystal mass
/// @param [in] eta Deposition factor
/// @param [in] dvsi Vapor excess with respect to ice sat
/// @param [in] dt Time step
/// @return Rate of vapor deposition to ice
///
template <typename T>
KOKKOS_INLINE_FUNCTION T vapor_x_ice(const T qi, const T mi, const T eta, const T dvsi, const T rho, const T dt) {
  constexpr T ami   = T{130.0};  // Formfactor for mass-size relation of cld ice
  constexpr T b_exp = T{-0.67};  // exp. for conv. (-1 + 0.33) of ice mass to sfc area
  const T a_fact    = static_cast<T>(4.0) * Kokkos::pow(ami, static_cast<T>(-1.0) / static_cast<T>(3.0));
  T result          = ZERO<T>;

  if (qi > graupel::qmin<T>) {
    result = (a_fact * eta) * rho * qi * Kokkos::pow(mi, b_exp) * dvsi;

    if (result > ZERO<T>) {
      result = Kokkos::fmin(result, dvsi / dt);
    } else {
      result = Kokkos::fmax(result, dvsi / dt);
      result = Kokkos::fmax(result, -qi / dt);
    }
  }

  return result;
}

///
/// @param [in] t Temperature
/// @param [in] p Ambient pressure
/// @param [in] rho Ambient density
/// @param [in] qs Snow specific mass
/// @param [in] ns Snow number
/// @param [in] lambda  Slope parameter (lambda) snow
/// @param [in] eta Deposition factor
/// @param [in] ice_dep Limiter for vapor dep on snow
/// @param [in] dvsw qv-qsat_water(T)
/// @param [in] dvsi  qv-qsat_ice(T)
/// @param [in] dvsw0 qv-qsat_water(T0)
/// @param [in] dt Time step
/// @return Rate of vapor deposition to snow
///
template <typename T>
KOKKOS_INLINE_FUNCTION T vapor_x_snow(const T t, const T p, const T rho, const T qs, const T ns, const T lambda,
                                      const T eta, const T ice_dep, const T dvsw, const T dvsi, const T dvsw0,
                                      const T dt) {
  constexpr T nu     = T{1.75e-5};  // kinematic viscosity of air
  constexpr T a0     = T{1.0};
  constexpr T a2     = -(graupel::v1s<T> + static_cast<T>(1.0)) / static_cast<T>(2.0);
  constexpr T eps    = T{1.e-15};
  constexpr T qs_lim = T{1.e-7};
  constexpr T cnx    = T{4.0};
  constexpr T b      = T{0.8};
  constexpr T c1     = T{31282.3};
  constexpr T c2     = T{0.241897};
  constexpr T c3     = T{0.28003};
  constexpr T c4     = T{-0.146293E-6};
  const T a1         = static_cast<T>(0.4182) * Kokkos::sqrt(graupel::v0s<T> / nu);

  T result           = ZERO<T>;

  if (qs > graupel::qmin<T>) {
    if (t < tmelt<T>) {
      result = (cnx * ns * eta / rho) * (a0 + a1 * Kokkos::pow(lambda, a2)) * dvsi / (lambda * lambda + eps);

      // GZ: This limitation, which was missing in the original graupel scheme,
      // is crucial for numerical stability in the tropics!
      // a meaningful distiction between cloud ice and snow

      if (result > ZERO<T>) result = Kokkos::fmin(result, dvsi / dt - ice_dep);
      if (qs <= qs_lim) result = Kokkos::fmin(result, ZERO<T>);
    } else {
      if (t > (tmelt<T> - graupel::tx<T> * dvsw0)) {
        result = (c1 / p + c2) * Kokkos::fmin(ZERO<T>, dvsw0) * Kokkos::pow(qs * rho, b);
      } else {
        result = (c3 + c4 * p) * dvsw * Kokkos::pow(qs * rho, b);
      }
    }
    result = Kokkos::fmax(result, -qs / dt);
  }

  return result;
}

}  // namespace transition

#endif  // RAGNAROK_AES_MICROPHYSICS_TRANSITION_H_
