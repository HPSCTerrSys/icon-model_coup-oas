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
#include "thermo.hpp"

#include <gtest/gtest.h>

namespace {
using testing::Types;

template <class T>
class ThermoTest : public testing::Test {
 public:
  void validate(T actual, T expected) {
    if constexpr (std::is_same_v<T, float>) {
      EXPECT_FLOAT_EQ(expected, actual);
    } else {
      EXPECT_DOUBLE_EQ(expected, actual);
    }
  }

  void validate(T actual, T expected_float, T expected_double) {
    if constexpr (std::is_same_v<T, float>) {
      EXPECT_FLOAT_EQ(expected_float, actual);
    } else {
      EXPECT_DOUBLE_EQ(expected_double, actual);
    }
  }
};  // class ThermoTest

TYPED_TEST_SUITE_P(ThermoTest);

TYPED_TEST_P(ThermoTest, CheckInternalEnergy) {
  TypeParam tk        = TypeParam{255.756};
  TypeParam qv        = TypeParam{0.00122576};
  TypeParam qliq      = TypeParam{1.63837e-20};
  TypeParam qice      = TypeParam{1.09462e-08};
  TypeParam rho       = TypeParam{0.83444};
  TypeParam dz        = TypeParam{249.569};
  TypeParam reference = TypeParam{38265357.270336017};

  TypeParam result    = thermo::internal_energy(tk, qv, qliq, qice, rho, dz);
  this->validate(result, reference);
}

TYPED_TEST_P(ThermoTest, CheckTInternalEnergy) {
  TypeParam u         = TypeParam{38265357.270336017};
  TypeParam qv        = TypeParam{0.00122576};
  TypeParam qliq      = TypeParam{1.63837e-20};
  TypeParam qice      = TypeParam{1.09462e-08};
  TypeParam rho       = TypeParam{0.83444};
  TypeParam dz        = TypeParam{249.569};
  TypeParam reference = TypeParam{255.75599999999997};

  TypeParam result    = thermo::T_from_internal_energy(u, qv, qliq, qice, rho, dz);
  this->validate(result, reference);
}

TYPED_TEST_P(ThermoTest, CheckSpecificHumidity) {
  TypeParam pvapor    = TypeParam{84.52};
  TypeParam ptotal    = TypeParam{94.52};
  TypeParam reference = TypeParam{0.84017368667728631};

  TypeParam result    = thermo::specific_humidity(pvapor, ptotal);
  this->validate(result, reference);
}

TYPED_TEST_P(ThermoTest, CheckSatPresWater) {
  TypeParam tk        = TypeParam{281.787};
  TypeParam reference = TypeParam{1120.1604149806028};

  TypeParam result    = thermo::sat_pres_water(tk);
  this->validate(result, reference);
}

TYPED_TEST_P(ThermoTest, CheckSatPresIce) {
  TypeParam tk               = TypeParam{281.787};
  TypeParam reference_float  = TypeParam{1216.774};
  TypeParam reference_double = TypeParam{1216.7746246067475};

  TypeParam result           = thermo::sat_pres_ice(tk);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(ThermoTest, CheckQSatRho) {
  TypeParam tk               = TypeParam{281.787};
  TypeParam rho              = TypeParam{1.24783};
  TypeParam reference_float  = TypeParam{0.0069027566};
  TypeParam reference_double = TypeParam{0.0069027592942577506};

  TypeParam result           = thermo::qsat_rho(tk, rho);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(ThermoTest, CheckQSatIceRho) {
  TypeParam tk               = TypeParam{281.787};
  TypeParam rho              = TypeParam{1.24783};
  TypeParam reference_float  = TypeParam{0.0074981214};
  TypeParam reference_double = TypeParam{0.0074981245870634101};

  TypeParam result           = thermo::qsat_ice_rho(tk, rho);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(ThermoTest, CheckdQSatdTRho) {
  TypeParam tk               = TypeParam{273.909};
  TypeParam qx               = TypeParam{0.00448941};
  TypeParam reference_float  = TypeParam{0.00030825072};
  TypeParam reference_double = TypeParam{0.00030825070286492049};

  TypeParam result           = thermo::dqsatdT_rho(qx, tk);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(ThermoTest, CheckdQSatdT) {
  TypeParam tk               = TypeParam{273.909};
  TypeParam qx               = TypeParam{0.00448941};
  TypeParam reference_float  = TypeParam{0.00032552675};
  TypeParam reference_double = TypeParam{0.00032552672594464161};

  TypeParam result           = thermo::dqsatdT(qx, tk);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(ThermoTest, CheckdQSatdTIce) {
  TypeParam tk               = TypeParam{273.909};
  TypeParam qx               = TypeParam{0.00448941};
  TypeParam reference_float  = TypeParam{0.00036880185};
  TypeParam reference_double = TypeParam{0.00036880177774940582};

  TypeParam result           = thermo::dqsatdT_ice(qx, tk);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(ThermoTest, CheckVaporizationEnergy) {
  TypeParam tk               = TypeParam{273.909};
  TypeParam reference_float  = TypeParam{2372625};
  TypeParam reference_double = TypeParam{2372624.9454889921};

  TypeParam result           = thermo::vaporization_energy(tk);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(ThermoTest, CheckSublimationEnergy) {
  TypeParam tk               = TypeParam{273.909};
  TypeParam reference_float  = TypeParam{2707907.2};
  TypeParam reference_double = TypeParam{2707907.2055500001};

  TypeParam result           = thermo::sublimation_energy(tk);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(ThermoTest, CheckPotentialTemperature) {
  TypeParam tk               = TypeParam{302.496350495765};
  TypeParam pres             = TypeParam{100529.557094430798};
  TypeParam reference_float  = TypeParam{302.0402192};
  TypeParam reference_double = TypeParam{302.04021921529369};

  TypeParam result           = thermo::potential_temperature(tk, pres);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(ThermoTest, CheckDewpointTemperature) {
  TypeParam tk               = TypeParam{271.4664669994};  // below 273.15K
  TypeParam qv               = TypeParam{0.0031244326};
  TypeParam pres             = TypeParam{98156.4651914430};
  TypeParam reference_float  = TypeParam{270.55487};
  TypeParam reference_double = TypeParam{270.55486409672693};

  TypeParam result           = thermo::dewpoint_temperature(tk, qv, pres);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(ThermoTest, CheckDewpointTemperatureMelt) {
  TypeParam tk               = TypeParam{295.8042511680};  // above 273.15K
  TypeParam qv               = TypeParam{0.0116970540};
  TypeParam pres             = TypeParam{101562.2742229487};
  TypeParam reference_float  = TypeParam{289.81256};
  TypeParam reference_double = TypeParam{289.81257531027944};

  TypeParam result           = thermo::dewpoint_temperature(tk, qv, pres);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(ThermoTest, CheckSaturationAdjustment) {
  TypeParam te                   = TypeParam{273.90911754406039};
  TypeParam qve                  = TypeParam{4.4913424511676030E-003};
  TypeParam qce                  = TypeParam{6.0066941654987605E-013};
  TypeParam qre                  = TypeParam{2.5939378002267028E-004};
  TypeParam qti                  = TypeParam{1.0746937601645517E-005};
  TypeParam rho                  = TypeParam{1.1371657035251757};
  TypeParam qce_reference_double = TypeParam{9.5724552280369163e-07};
  TypeParam qce_reference_single = TypeParam{9.5600262e-07};

  thermo::saturation_adjustment(te, qve, qce, qre, qti, rho);

  this->validate(te, static_cast<TypeParam>(273.91226488486984));
  this->validate(qve, static_cast<TypeParam>(0.004490385206245469));
  this->validate(qce, qce_reference_single, qce_reference_double);
}

TYPED_TEST_P(ThermoTest, CheckSaturationAdjustment2D) {
  const int nproma             = 2;
  const int nlev               = 1;

  TypeParam zta[nproma * nlev] = {TypeParam{273.90911754406039}, TypeParam{271.90911754406039}};
  auto h_ta                    = HostView2D<TypeParam>(zta, nlev, nproma);
  auto d_ta                    = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_ta);

  TypeParam zqv[nproma * nlev] = {TypeParam{4.4913424511676030E-003}, TypeParam{3.4913424511676030E-003}};
  auto h_qv                    = HostView2D<TypeParam>(zqv, nlev, nproma);
  auto d_qv                    = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qv);

  TypeParam zqc[nproma * nlev] = {TypeParam{6.0066941654987605E-013}, TypeParam{5.0066941654987605E-013}};
  auto h_qc                    = HostView2D<TypeParam>(zqc, nlev, nproma);
  auto d_qc                    = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qc);

  TypeParam zqr[nproma * nlev] = {TypeParam{2.5939378002267028E-004}, TypeParam{2.8939378002267028E-004}};
  auto h_qr                    = HostView2D<TypeParam>(zqr, nlev, nproma);
  auto d_qr                    = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qr);

  TypeParam zti[nproma * nlev] = {
      TypeParam{4.4913424511676030E-005 + 2.5939378002267028E-005 + 6.0066941654987605E-005},
      TypeParam{3.4913424511676030E-005 + 2.8939378002267028E-005 + 2.8939378002267028E-005}};
  auto h_ti                    = HostView2D<TypeParam>(zti, nlev, nproma);
  auto d_ti                    = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_ti);

  TypeParam rho[nproma * nlev] = {TypeParam{1.1371657035251757}, TypeParam{1.1871657035251757}};
  auto h_rho                   = HostView2D<TypeParam>(rho, nlev, nproma);
  auto d_rho                   = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_rho);

  thermo::saturation_adjustment_2d<TypeParam>(0, nlev, 0, nproma, d_ta, d_qv, d_qc, d_qr, d_ti, d_rho);

  Kokkos::deep_copy(h_ta, d_ta);
  Kokkos::deep_copy(h_qv, d_qv);
  Kokkos::deep_copy(h_qc, d_qc);

  this->validate(h_ta(0, 0), static_cast<TypeParam>(273.91226), static_cast<TypeParam>(273.91226452301134));
  this->validate(h_ta(0, 1), static_cast<TypeParam>(271.90912), static_cast<TypeParam>(271.90911754240938));
  this->validate(h_qv(0, 0), static_cast<TypeParam>(0.0044903862), static_cast<TypeParam>(0.0044903850946812571));
  this->validate(h_qv(0, 1), static_cast<TypeParam>(0.0034913425), static_cast<TypeParam>(0.0034913424516682724));
  this->validate(h_qc(0, 0), static_cast<TypeParam>(9.5600262e-07), static_cast<TypeParam>(9.5735708701555344e-07));
  this->validate(h_qc(0, 1), static_cast<TypeParam>(0), static_cast<TypeParam>(0.));
}

REGISTER_TYPED_TEST_SUITE_P(ThermoTest, CheckInternalEnergy, CheckTInternalEnergy, CheckSpecificHumidity,
                            CheckSatPresWater, CheckSatPresIce, CheckQSatRho, CheckQSatIceRho, CheckdQSatdT,
                            CheckdQSatdTRho, CheckdQSatdTIce, CheckVaporizationEnergy, CheckSublimationEnergy,
                            CheckPotentialTemperature, CheckDewpointTemperature, CheckDewpointTemperatureMelt,
                            CheckSaturationAdjustment, CheckSaturationAdjustment2D);
using MyTypes = ::testing::Types<float, double>;
INSTANTIATE_TYPED_TEST_SUITE_P(RagnarokTests, ThermoTest, MyTypes);

};  // namespace
