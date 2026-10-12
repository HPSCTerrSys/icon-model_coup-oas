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
namespace thermo {

/// @brief Calculates internal energy
///
/// @param [in] tk Temperature (Kelvin)
/// @param [in] qliq Specific mass of liquid phases
/// @param [in] qice Specific mass of solid phases
/// @param [in] rho Density
/// @param [in] dz Extend of the grid cell
/// @return Internal energy
template <typename T>
KOKKOS_INLINE_FUNCTION T internal_energy(const T tk, const T qv, const T qliq, const T qice, const T rho, const T dz) {
  // total water specific mass
  T qtot = qliq + qice + qv;

  // moist isometric specific heat
  T cv   = cvd<T> * (static_cast<T>(1.0) - qtot) + cvv<T> * qv + clw<T> * qliq + ci<T> * qice;

  return rho * dz * (cv * tk - qliq * lvc<T> - qice * lsc<T>);
}

/// @brief Calculates temperature from internal energy
///
/// @param [in] U  Internal energy (extensive)
/// @param [in] qv Water vapor specific humidity
/// @param [in] qliq Specific mass of liquid phases
/// @param [in] qice Specific mass of solid phases
/// @param [in] rho Density
/// @param [in] dz Density
/// @return Temperature
template <typename T>
KOKKOS_INLINE_FUNCTION T T_from_internal_energy(const T U, const T qv, const T qliq, const T qice, const T rho,
                                                const T dz) {
  // total water specific mass
  T qtot = qliq + qice + qv;

  // moist isometric specific heat
  T cv   = (cvd<T> * (static_cast<T>(1.0) - qtot) + cvv<T> * qv + clw<T> * qliq + ci<T> * qice) * rho * dz;

  return (U + rho * dz * (qliq * lvc<T> + qice * lsc<T>)) / cv;
}

/// @brief Calculates specific humidity from vapor and total pressure
///
/// @param [in] pvapor Vapor pressure
/// @param [in] ptotal Total pressure
/// @return Humidity
template <typename T>
KOKKOS_INLINE_FUNCTION T specific_humidity(const T pvapor, const T ptotal) {
  constexpr T rdv     = rd<T> / rv<T>;
  constexpr T o_m_rdv = static_cast<T>(1.0) - rdv;

  return rdv * pvapor / (ptotal - o_m_rdv * pvapor);
}

/// @brief Calculates saturation pressure over water
///
/// @param [in] tk Temperature (Kelvin)
/// @return Saturation pressure
template <typename T>
KOKKOS_INLINE_FUNCTION T sat_pres_water(const T tk) {
  return c1es<T> * Kokkos::exp(c3les<T> * (tk - tmelt<T>) / (tk - c4les<T>));
}

/// @brief Calculates saturation pressure over ice
///
/// @param [in] tk Temperature (Kelvin)
/// @return Saturation pressure
template <typename T>
KOKKOS_INLINE_FUNCTION T sat_pres_ice(const T tk) {
  return c1es<T> * Kokkos::exp(c3ies<T> * (tk - tmelt<T>) / (tk - c4ies<T>));
}

/// @brief Calculates saturation vapor pressure (over liquid) at constant density
///
/// @param [in] tk Temperature (Kelvin)
/// @param [in] rho Density
/// @return saturation pressure
template <typename T>
KOKKOS_INLINE_FUNCTION T qsat_rho(const T tk, const T rho) {
  return sat_pres_water(tk) / (rho * rv<T> * tk);
}

/// @brief Saturation vapor pressure (over ice) at constant density
///
/// @param [in] tk Temperature (Kelvin)
/// @param [in] rho Density
/// @return saturation pressure
template <typename T>
KOKKOS_INLINE_FUNCTION T qsat_ice_rho(const T tk, const T rho) {
  return sat_pres_ice(tk) / (rho * rv<T> * tk);
}

/// @brief Computes the derivative d(qsat_rho)/dT
///
/// @param [in] qs Saturation vapor pressure (over liquid)
/// @param [in] tk Temperature (Kelvin)
/// @return derivative d(qsat_rho)/dT
template <typename T>
KOKKOS_INLINE_FUNCTION T dqsatdT_rho(const T qs, const T tk) {
  return qs * (c5les<T> / Kokkos::pow(tk - c4les<T>, static_cast<T>(2.0)) - static_cast<T>(1.0) / tk);
}

/// @brief:  Computes the derivative of the saturation over water vapour by temperature
///
/// @param [in] qs Saturation vapor pressure
/// @param [in] tk Temperature (Kelvin)
/// @return derivative d(qsat)/dT
template <typename T>
KOKKOS_INLINE_FUNCTION T dqsatdT(const T qs, const T tk) {
  return c5les<T> * (static_cast<T>(1.0) + vtmpc1<T> * qs) * qs / Kokkos::pow((tk - c4les<T>), static_cast<T>(2));
}

/// @brief : Computes the derivative of the saturation over ice by temperature
/// @param [in] qs Saturation vapor pressure
/// @param [in] tk Temperature (Kelvin)
/// @return derivative d(qsat)/dT
template <typename T>
KOKKOS_INLINE_FUNCTION T dqsatdT_ice(const T qs, const T tk) {
  return c5ies<T> * (static_cast<T>(1.0) + vtmpc1<T> * qs) * qs / Kokkos::pow((tk - c4ies<T>), static_cast<T>(2.0));
}

/// @brief Computes internal energy of vaporization
///
/// @param [in] tk Temperature (Kelvin)
/// @return Energy of vaporization
template <typename T>
KOKKOS_INLINE_FUNCTION T vaporization_energy(const T tk) {
  constexpr T tmp = cvv<T> - clw<T>;
  return lvc<T> + tmp * tk;
}

/// @brief Computes internal energy of sublimation
///
/// @param [in] t Temperature (Kelvin)
/// @return Energy of sublimation
template <typename T>
KOKKOS_INLINE_FUNCTION T sublimation_energy(const T tk) {
  return als<T> + (cpv<T> - ci<T>)*(tk - tmelt<T>)-rv<T> * tk;
}

/// @brief Calculate potential temperature
///
/// @param [in] tk Temperature (Kelvin)
/// @param [in] pres Pressure
/// @return Potential temperature
template <typename T>
KOKKOS_INLINE_FUNCTION T potential_temperature(const T tk, const T pres) {
  return tk * Kokkos::exp(rd_o_cpd<T> * Kokkos::log(p0ref<T> / pres));
}

/// @brief Calculate dewpoint temperature
///
/// @param [in] tk Temperature (Kelvin)
/// @param [in] qv Specific humidity
/// @param [in] pres Pressure
/// @return Potential temperature
template <typename T>
KOKKOS_INLINE_FUNCTION T dewpoint_temperature(const T tk, const T qv, const T pres) {
  T zfrac, zcvm3, zcvm4;

  if (tk > tmelt<T>) {
    zcvm3 = c3les<T>;
    zcvm4 = c4les<T>;
  } else {
    zcvm3 = c3ies<T>;
    zcvm4 = c4ies<T>;
  }
  zfrac = Kokkos::log(pres * qv / (c2es<T> * (static_cast<T>(1.) + vtmpc1<T> * qv))) / zcvm3;

  return Kokkos::fmin(tk, (tmelt<T> - zfrac * zcvm4) / (static_cast<T>(1.) - zfrac));
}

/// @brief Partitions water mass to maintain saturation
///
/// @param [inout] te Temperature (Kelvin)
/// @param [inout] qve specific humidity
/// @param [inout] qce specific cloud water content
/// @param [in] qre specific rain water
/// @param [in] qti specific mass of all ice species (total-ice)
/// @param [in] rho density containing dry air and water constituents
///
/// Description:
///   This routine performs the saturation adjustment to find the combination
///   of cloud water and temperature that is in equilibirum at the same internal
///   energy as the initial fields.
///
/// Method:
///   Saturation adjustment in the presence of non-zero cloud water requires
///   solving a non-linear equation, which is done using a Newton-Raphson
///   method.  The procedure first checks for the special case of sub-saturation
///   in which case the solution can be Calculated directly.   If not then the
///   the solver looks for the zero of the function f(T) denoted fT, whereby T
///   is temperature and f is the difference between the internal energy at T
///   and the internal energy at the initial T and qc, denoted ue.
template <typename T>
KOKKOS_INLINE_FUNCTION void saturation_adjustment(T& te, T& qve, T& qce, const T qre, const T qti, const T rho) {
  T qt  = qve + qce + qre + qti;
  T cvc = cvd<T> * (static_cast<T>(1.0) - qt) + clw<T> * qre + ci<T> * qti;
  T cv  = cvc + cvv<T> * qve + clw<T> * qce;
  T ue  = cv * te - qce * lvc<T>;
  T Tx  = ue / (cv + qce * (cvv<T> - clw<T>));
  T qx  = qsat_rho(Tx, rho);

  //  If subsaturated upon evaporating all cloud water, T can be diagnosed
  //  explicitly, ! so test for this.  If not then T needs to be solved for
  //  iteratively.
  if (qve + qce <= qx) {
    qve = qve + qce;
    qce = static_cast<T>(0.);
  } else {
    Tx = te;
    for (size_t i = 0; i < 6; ++i) {
      qx    = qsat_rho(Tx, rho);
      T dqx = dqsatdT_rho(qx, Tx);
      T qcx = qve + qce - qx;
      cv    = cvc + cvv<T> * qx + clw<T> * qcx;
      T ux  = cv * Tx - qcx * lvc<T>;
      T dux = cv + dqx * (lvc<T> + (cvv<T> - clw<T>)*Tx);
      Tx    = Tx - (ux - ue) / dux;
    }
    qx  = qsat_rho(Tx, rho);
    qce = Kokkos::fmax(qve + qce - qx, static_cast<T>(0.));
    qve = qx;
  }
  te = Tx;
}

template <typename T>
void saturation_adjustment_2d(const int kstart, const int kend, const int jcs, const int jce, View2D<T> ta_view,
                              View2D<T> qv_view, View2D<T> qc_view, View2D<T> qr_view, View2D<T> total_ice,
                              View2D<T> rho_view) {
  Kokkos::parallel_for(
      "saturation_adjustment_2d", RangePolicy2D({kstart, jcs}, {kend, jce}), KOKKOS_LAMBDA(const int jk, const int jc) {
        saturation_adjustment(ta_view(jk, jc), qv_view(jk, jc), qc_view(jk, jc), qr_view(jk, jc), total_ice(jk, jc),
                              rho_view(jk, jc));
      });
  Kokkos::fence();
}

template <typename T>
void saturation_adjustment_2d(const int nproma, const int nlev, const int kstart, const int kend, const int jcs,
                              const int jce, T* zta, T* qv, T* qc, T* qr, T* total_ice, T* rho) {
  View2D<T> ta_view(zta, nlev, nproma);
  View2D<T> qv_view(qv, nlev, nproma);
  View2D<T> qc_view(qc, nlev, nproma);
  View2D<T> qr_view(qr, nlev, nproma);
  View2D<T> total_ice_view(total_ice, nlev, nproma);
  View2D<T> rho_view(rho, nlev, nproma);

  saturation_adjustment_2d<T>(kstart, kend, jcs, jce, ta_view, qv_view, qc_view, qr_view, total_ice_view, rho_view);
}

}  // namespace thermo
