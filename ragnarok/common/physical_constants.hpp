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
/// Thermodynamic constants for the dry and moist atmosphere

#ifndef RAGNAROK_COMMON_PHYSICAL_CONSTANTS_H_
#define RAGNAROK_COMMON_PHYSICAL_CONSTANTS_H_

// dry air

/// [J/K/kg] gas constant
template <typename T>
constexpr T rd = T{287.04};
/// [J/K/kg] specific heat at constant pressure
template <typename T>
constexpr T cpd = T{1004.64};
/// [J/K/kg] specific heat at constant volume
template <typename T>
constexpr T cvd = cpd<T> - rd<T>;
/// [m^2/s]  kinematic viscosity of dry air
template <typename T>
constexpr T con_m = T{1.50E-5};
/// [m^2/s] scalar conductivity of dry air
template <typename T>
constexpr T con_h = T{2.20E-5};
/// [J/m/s/K] thermal conductivity of dry air
template <typename T>
constexpr T con0_h = T{2.40e-2};
/// [N*s/m2] dyn viscosity of dry air at tmelt
template <typename T>
constexpr T eta0d = T{1.717e-5};

// H2O
// gas
/// [J/K/kg] gas constant for water vapor
template <typename T>
constexpr T rv = T{461.51};
/// [J/K/kg] specific heat at constant pressure
template <typename T>
constexpr T cpv = T{1869.46};
/// [J/K/kg] specific heat at constant volume
template <typename T>
constexpr T cvv = cpv<T> - rv<T>;
/// [m^2/s] diff coeff of H2O vapor in dry air at tmelt
template <typename T>
constexpr T dv0 = T{2.22e-5};

// liquid / water
/// [kg/m3] density of liquid water
template <typename T>
constexpr T rhoh2o = T{1000.};

// solid / ice
/// [kg/m3] density of pure ice
template <typename T>
constexpr T rhoice = T{916.7};
template <typename T>
constexpr T cv_i = T{2000.0};

// phase changes
/// [J/kg] latent heat for vaporisation
template <typename T>
constexpr T alv = T{2.5008e6};
/// [J/kg] latent heat for sublimation
template <typename T>
constexpr T als = T{2.8345e6};
/// [J/kg] latent heat for fusion
template <typename T>
constexpr T alf = als<T> - alv<T>;
/// [K] melting temperature of ice/snow
template <typename T>
constexpr T tmelt = T{273.15};
/// [K] Triple point of water at 611hPa
template <typename T>
constexpr T t3 = T{273.16};

// Auxiliary constants
template <typename T>
constexpr T rdv = rd<T> / rv<T>;
template <typename T>
constexpr T vtmpc1 = rv<T> / rd<T> - static_cast<T>(1.);
template <typename T>
constexpr T vtmpc2 = cpv<T> / cpd<T> - static_cast<T>(1.);
template <typename T>
constexpr T rcpv = cpd<T> / cpv<T> - static_cast<T>(1.);
template <typename T>
constexpr T alvdcp = alv<T> / cpd<T>;  // [K]
template <typename T>
constexpr T alsdcp = als<T> / cpd<T>;  // [K]
template <typename T>
constexpr T rcpd = static_cast<T>(1.) / cpd<T>;  // [K*kg/J]
template <typename T>
constexpr T rcvd = static_cast<T>(1.) / cvd<T>;  // [K*kg/J]
template <typename T>
constexpr T rcpl = T{3.1733};

/// specific heat capacity of liquid water
template <typename T>
constexpr T clw = (rcpl<T> + static_cast<T>(1.0)) * cpd<T>;
template <typename T>
constexpr T cv_v = (rcpv<T> + static_cast<T>(1.0)) * cpd<T> - rv<T>;
template <typename T>
constexpr T rd_o_cpd = rd<T> / cpd<T>;

/// reference pressure for Exner function
template <typename T>
constexpr T p0ref = T{100000.0};  // [Pa]

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

#endif  // RAGNAROK_COMMON_PHYSICAL_CONSTANTS_H_
