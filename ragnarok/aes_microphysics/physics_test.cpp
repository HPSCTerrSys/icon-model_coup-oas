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
#include "aes_microphysics/physics.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <type_traits>

namespace {
using testing::Types;

template <typename T>
class GraupelPhysicsTest : public testing::Test {
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
};  // class GraupelPhysicsTest

TYPED_TEST_SUITE_P(GraupelPhysicsTest);

TYPED_TEST_P(GraupelPhysicsTest, DepAutoConversionDefault) {
  TypeParam qi      = TypeParam{2.02422e-23};
  TypeParam m_ice   = TypeParam{1e-12};
  TypeParam ice_dep = TypeParam{-2.06276e-05};

  TypeParam result  = physics::deposition_auto_conversion(qi, m_ice, ice_dep);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(GraupelPhysicsTest, DepAutoConversion) {
  TypeParam qi        = TypeParam{2.02422e-2};
  TypeParam m_ice     = TypeParam{1e-12};
  TypeParam ice_dep   = TypeParam{-2.06276e-05};
  TypeParam reference = TypeParam{6.6430804299795412e-08};

  TypeParam result    = physics::deposition_auto_conversion(qi + graupel::qmin<TypeParam>, m_ice, -ice_dep);
  this->validate(result, reference);
}

TYPED_TEST_P(GraupelPhysicsTest, IceDepositionNucleationDefault) {
  TypeParam t      = TypeParam{272.731};
  TypeParam qc     = TypeParam{0.0};
  TypeParam qi     = TypeParam{2.02422e-23};
  TypeParam ni     = TypeParam{5.05089};
  TypeParam dvsi   = TypeParam{-0.000618828};
  TypeParam dt     = TypeParam{30.0};

  TypeParam result = physics::ice_deposition_nucleation(t, qc, qi, ni, dvsi, dt);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(GraupelPhysicsTest, IceDepositionNucleation) {
  TypeParam t         = TypeParam{160.9};
  TypeParam qc        = TypeParam{1.0e-2};
  TypeParam qi        = TypeParam{2.02422e-23};
  TypeParam ni        = TypeParam{5.05089};
  TypeParam dvsi      = TypeParam{0.0001};
  TypeParam dt        = TypeParam{30.0};
  TypeParam reference = TypeParam{1.6836299999999999e-13};

  TypeParam result    = physics::ice_deposition_nucleation(t, qc, qi, ni, dvsi, dt);
  this->validate(result, reference);
}

TYPED_TEST_P(GraupelPhysicsTest, DepFactor) {
  TypeParam t         = TypeParam{272.731};
  TypeParam qvsi      = TypeParam{0.00416891};
  TypeParam reference = TypeParam{1.3234329478493952e-05};

  TypeParam result    = physics::deposition_factor(t, qvsi);
  this->validate(result, reference);
}

TYPED_TEST_P(GraupelPhysicsTest, IceSticking) {
  TypeParam reference = TypeParam{1.0};
  TypeParam cia       = TypeParam{1.0};

  TypeParam result    = physics::ice_sticking(tmelt<TypeParam>, cia);
  this->validate(result, reference);
}

TYPED_TEST_P(GraupelPhysicsTest, IceNumber) {
  TypeParam t                = TypeParam{272.731};
  TypeParam rho              = TypeParam{1.12442};
  TypeParam reference_float  = TypeParam{5.0508094};
  TypeParam reference_double = TypeParam{5.0507995893464388};

  TypeParam result           = physics::ice_number(t, rho);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(GraupelPhysicsTest, IceMass) {
  TypeParam qi        = TypeParam{2.02422e-23};
  TypeParam ni        = TypeParam{5.05089};
  TypeParam reference = TypeParam{1e-12};

  TypeParam result    = physics::ice_mass(qi, ni);
  this->validate(result, reference);
}

TYPED_TEST_P(GraupelPhysicsTest, SnowLambdaDefault) {
  TypeParam rho_s  = TypeParam{0.12204} * graupel::qmin<TypeParam>;
  TypeParam ns     = TypeParam{1.76669e+07};
  TypeParam lmd_0  = TypeParam{1.0e+10};

  TypeParam result = physics::snow_lambda(rho_s, ns);
  this->validate(result, lmd_0);
}

TYPED_TEST_P(GraupelPhysicsTest, SnowLambda) {
  TypeParam rho_s            = TypeParam{1.12204} * graupel::qmin<TypeParam>;
  TypeParam ns               = TypeParam{1.76669e+07};
  TypeParam reference_float  = TypeParam{12952211};
  TypeParam reference_double = TypeParam{12952204.691943912};

  TypeParam result           = physics::snow_lambda(rho_s, ns);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(GraupelPhysicsTest, SnowNumberDefault) {
  TypeParam t      = TypeParam{276.302};
  TypeParam rho_s  = TypeParam{1.17797} * TypeParam{8.28451e-24};
  TypeParam n0s0   = TypeParam{8.00e+5};

  TypeParam result = physics::snow_number(t, rho_s);
  this->validate(result, n0s0);
}

TYPED_TEST_P(GraupelPhysicsTest, SnowNumber) {
  TypeParam t         = TypeParam{276.302};
  TypeParam rho_s     = TypeParam{1.17797} * TypeParam{8.28451e-4};
  TypeParam reference = TypeParam{3813750};

  TypeParam result    = physics::snow_number(t, rho_s);
  this->validate(result, reference);
}

TYPED_TEST_P(GraupelPhysicsTest, VmRain) {
  TypeParam xrho      = TypeParam{1.16163};
  TypeParam rho       = TypeParam{0.907829};
  TypeParam t         = TypeParam{257.501};
  TypeParam reference = TypeParam{4.8407268925025155};

  TypeParam result    = physics::vm(idx::lqg, xrho, rho, t);
  this->validate(result, reference);
}

TYPED_TEST_P(GraupelPhysicsTest, VmIce) {
  TypeParam xrho      = TypeParam{1.17873};
  TypeParam rho       = TypeParam{0.881673};
  TypeParam t         = TypeParam{257.458};
  TypeParam reference = TypeParam{0.40334524060776378};

  TypeParam result    = physics::vm(idx::lqi, xrho, rho, t);
  this->validate(result, reference);
}

TYPED_TEST_P(GraupelPhysicsTest, VmSnow) {
  TypeParam xrho             = TypeParam{1.17787};
  TypeParam rho              = TypeParam{0.882961};
  TypeParam t                = TypeParam{257.101};
  TypeParam reference_float  = TypeParam{2.5447843};
  TypeParam reference_double = TypeParam{2.5447819734203958};

  TypeParam result           = physics::vm(idx::lqs, xrho, rho, t);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(GraupelPhysicsTest, VmGraupel) {
  TypeParam xrho      = TypeParam{1.40017};
  TypeParam rho       = TypeParam{0.907829};
  TypeParam t         = TypeParam{257.501};
  TypeParam reference = TypeParam{4.8407268925025155};

  TypeParam result    = physics::vm(idx::lqg, xrho, rho, t);
  this->validate(result, reference);
}

REGISTER_TYPED_TEST_SUITE_P(GraupelPhysicsTest, DepAutoConversionDefault, DepAutoConversion,
                            IceDepositionNucleationDefault, IceDepositionNucleation, DepFactor, IceSticking, IceNumber,
                            IceMass, SnowLambdaDefault, SnowLambda, SnowNumberDefault, SnowNumber, VmRain, VmIce,
                            VmSnow, VmGraupel);
using MyTypes = ::testing::Types<float, double>;
INSTANTIATE_TYPED_TEST_SUITE_P(RagnarokTests, GraupelPhysicsTest, MyTypes);

}  // namespace
