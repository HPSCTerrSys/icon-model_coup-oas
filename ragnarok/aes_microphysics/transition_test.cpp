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
#include "aes_microphysics/transition.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <type_traits>

namespace {
using testing::Types;

template <typename T>
class HydrometeorTransitionTest : public testing::Test {
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
};  // class HydrometeorTransitionTest

TYPED_TEST_SUITE_P(HydrometeorTransitionTest);

TYPED_TEST_P(HydrometeorTransitionTest, CloudToGraupelDefault) {
  TypeParam t      = TypeParam{281.787};
  TypeParam rho    = TypeParam{1.24783};
  TypeParam qc     = ZERO<TypeParam>;
  TypeParam qg     = TypeParam{1.03636e-25};

  TypeParam result = transition::cloud_to_graupel(t, rho, qc, qg);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, CloudToGraupel) {
  TypeParam t         = TypeParam{256.983};
  TypeParam rho       = TypeParam{0.909677};
  TypeParam qc        = TypeParam{8.60101e-06};
  TypeParam qg        = TypeParam{4.11575e-06};
  TypeParam reference = TypeParam{2.7054723496793982e-10};

  TypeParam result    = transition::cloud_to_graupel(t, rho, qc, qg);
  this->validate(result, reference);
}

TYPED_TEST_P(HydrometeorTransitionTest, CloudToRainDefault) {
  TypeParam t      = TypeParam{281.787};
  TypeParam rho    = TypeParam{0.909677};
  TypeParam qc     = TypeParam{0.0};
  TypeParam qr     = TypeParam{5.2312e-07};
  TypeParam nc     = TypeParam{100};

  TypeParam result = transition::cloud_to_rain(t, rho, qc, qr, nc);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, CloudToRain) {
  TypeParam t                = TypeParam{267.25};
  TypeParam rho              = TypeParam{0.909677};
  TypeParam qc               = TypeParam{5.52921e-05};
  TypeParam qr               = TypeParam{2.01511e-12};
  TypeParam nc               = TypeParam{100};
  TypeParam reference_float  = TypeParam{0.0045578829};
  TypeParam reference_double = TypeParam{0.0045484481075162512};

  TypeParam result           = transition::cloud_to_rain(t, rho, qc, qr, nc);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(HydrometeorTransitionTest, CloudToSnowDefault) {
  TypeParam t      = TypeParam{281.787};
  TypeParam qc     = TypeParam{0.0};
  TypeParam qs     = TypeParam{3.63983e-40};
  TypeParam ns     = TypeParam{800000};
  TypeParam lambda = TypeParam{1e+10};

  TypeParam result = transition::cloud_to_snow(t, qc, qs, ns, lambda);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, CloudToSnow) {
  TypeParam t         = TypeParam{256.571};
  TypeParam qc        = TypeParam{3.31476e-05};
  TypeParam qs        = TypeParam{7.47365e-06};
  TypeParam ns        = TypeParam{3.37707e+07};
  TypeParam lambda    = TypeParam{8989.78};
  TypeParam reference = TypeParam{2.86295623693317e-09};

  TypeParam result    = transition::cloud_to_snow(t, qc, qs, ns, lambda);
  this->validate(result, reference);
}

TYPED_TEST_P(HydrometeorTransitionTest, CloudXIceDefault) {
  TypeParam t      = TypeParam{256.835};
  TypeParam qc     = TypeParam{0.0};
  TypeParam qi     = TypeParam{4.50245e-07};
  TypeParam dt     = TypeParam{30};

  TypeParam result = transition::cloud_x_ice(t, qc, qi, dt);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, CloudXIceI) {
  TypeParam t         = tmelt<TypeParam> + static_cast<TypeParam>(1.0);
  TypeParam qc        = TypeParam{0.0};
  TypeParam qi        = TypeParam{4.50245e-07};
  TypeParam dt        = TypeParam{30};
  TypeParam reference = TypeParam{-1.5008166666666666e-08};

  TypeParam result    = transition::cloud_x_ice(t, qc, qi, dt);
  this->validate(result, reference);
}

TYPED_TEST_P(HydrometeorTransitionTest, CloudXIceC) {
  TypeParam t         = TypeParam{236};
  TypeParam qc        = TypeParam{0.000198441};
  TypeParam qi        = TypeParam{1.55543e-19};
  TypeParam dt        = TypeParam{30};
  TypeParam reference = TypeParam{6.6146999999999995e-06};

  TypeParam result    = transition::cloud_x_ice(t, qc, qi, dt);
  this->validate(result, reference);
}

TYPED_TEST_P(HydrometeorTransitionTest, GraupelToRainDefault) {
  TypeParam t      = TypeParam{280.156};
  TypeParam p      = TypeParam{98889.4};
  TypeParam rho    = TypeParam{1.22804};
  TypeParam dvsw0  = TypeParam{-0.00167867};
  TypeParam qg     = TypeParam{1.53968e-17};

  TypeParam result = transition::graupel_to_rain(t, p, rho, dvsw0, qg);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, GraupelToRain) {
  TypeParam t                = TypeParam{280.156};
  TypeParam p                = TypeParam{98889.4};
  TypeParam rho              = TypeParam{1.22804};
  TypeParam dvsw0            = TypeParam{-0.00167867};
  TypeParam qg               = TypeParam{1.53968e-15};
  TypeParam reference_float  = TypeParam{5.9748441e-13};
  TypeParam reference_double = TypeParam{5.9748142538569357e-13};

  TypeParam result           = transition::graupel_to_rain(t, p, rho, dvsw0, qg);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(HydrometeorTransitionTest, IceToGraupelDefault) {
  TypeParam rho          = TypeParam{1.12442};
  TypeParam qr           = TypeParam{1.34006e-17};
  TypeParam qg           = TypeParam{1.22571e-13};
  TypeParam qi           = TypeParam{2.02422e-23};
  TypeParam sticking_eff = TypeParam{0.962987};

  TypeParam result       = transition::ice_to_graupel(rho, qr, qg, qi, sticking_eff);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, IceToGraupelR) {
  TypeParam rho          = TypeParam{1.04848};
  TypeParam qr           = TypeParam{6.00408e-13};
  TypeParam qg           = TypeParam{1.19022e-18};
  TypeParam qi           = TypeParam{1.9584e-08};
  TypeParam sticking_eff = TypeParam{0.518393};
  TypeParam reference    = TypeParam{7.1049436957697864e-19};

  TypeParam result       = transition::ice_to_graupel(rho, qr, qg, qi, sticking_eff);
  this->validate(result, reference);
}

TYPED_TEST_P(HydrometeorTransitionTest, IceToGraupelG) {
  TypeParam rho          = TypeParam{1.04848};
  TypeParam qr           = TypeParam{6.00408e-16};
  TypeParam qg           = TypeParam{1.19022e-05};
  TypeParam qi           = TypeParam{1.9584e-08};
  TypeParam sticking_eff = TypeParam{0.518393};
  TypeParam reference    = TypeParam{5.557203647060206e-13};

  TypeParam result       = transition::ice_to_graupel(rho, qr, qg, qi, sticking_eff);
  this->validate(result, reference);
}

TYPED_TEST_P(HydrometeorTransitionTest, IceToSnowDefault) {
  TypeParam qi           = TypeParam{7.95122e-25};
  TypeParam ns           = TypeParam{2.23336e+07};
  TypeParam lambda       = TypeParam{61911.1};
  TypeParam sticking_eff = TypeParam{0.241568};

  TypeParam result       = transition::ice_to_snow(qi, ns, lambda, sticking_eff);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, IceToSnow) {
  TypeParam qi           = TypeParam{6.43223e-08};
  TypeParam ns           = TypeParam{1.93157e+07};
  TypeParam lambda       = TypeParam{10576.8};
  TypeParam sticking_eff = TypeParam{0.511825};
  TypeParam reference    = TypeParam{3.3262745200740486e-11};

  TypeParam result       = transition::ice_to_snow(qi, ns, lambda, sticking_eff);
  this->validate(result, reference);
}

TYPED_TEST_P(HydrometeorTransitionTest, RainToGraupelDefault) {
  TypeParam t      = TypeParam{272.731};
  TypeParam rho    = TypeParam{1.12442};
  TypeParam qc     = TypeParam{0.0};
  TypeParam qr     = TypeParam{1.34006e-17};
  TypeParam qi     = TypeParam{2.02422e-23};
  TypeParam qs     = TypeParam{1.02627e-19};
  TypeParam mi     = TypeParam{1e-12};
  TypeParam dvsw   = TypeParam{-0.000635669};
  TypeParam dt     = TypeParam{30};

  TypeParam result = transition::rain_to_graupel(t, rho, qc, qr, qi, qs, mi, dvsw, dt);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, RainToGraupel1) {
  TypeParam t                = TypeParam{258.542};
  TypeParam rho              = TypeParam{0.956089};
  TypeParam qc               = TypeParam{0.0};
  TypeParam qr               = TypeParam{3.01332e-11};
  TypeParam qi               = TypeParam{5.57166e-06};
  TypeParam qs               = TypeParam{3.55432e-05};
  TypeParam mi               = TypeParam{1e-09};
  TypeParam dvsw             = TypeParam{0.0};
  TypeParam dt               = TypeParam{30};
  TypeParam reference_float  = TypeParam{5.1570357e-17};
  TypeParam reference_double = TypeParam{5.1570340525841922e-17};

  TypeParam result           = transition::rain_to_graupel(t, rho, qc, qr, qi, qs, mi, dvsw, dt);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(HydrometeorTransitionTest, RainToGraupel2) {
  TypeParam t         = TypeParam{230.542};
  TypeParam rho       = TypeParam{0.956089};
  TypeParam qc        = TypeParam{8.6157e-05};
  TypeParam qr        = TypeParam{3.01332e-11};
  TypeParam qi        = TypeParam{5.57166e-06};
  TypeParam qs        = TypeParam{3.55432e-05};
  TypeParam mi        = TypeParam{1e-09};
  TypeParam dvsw      = TypeParam{0.0};
  TypeParam dt        = TypeParam{30};
  TypeParam reference = TypeParam{1.0044914238516472e-12};

  TypeParam result    = transition::rain_to_graupel(t, rho, qc, qr, qi, qs, mi, dvsw, dt);
  this->validate(result, reference);
}

TYPED_TEST_P(HydrometeorTransitionTest, RainToGraupel3) {
  TypeParam t         = TypeParam{258.542};
  TypeParam rho       = TypeParam{0.956089};
  TypeParam qc        = TypeParam{8.6157e-05};
  TypeParam qr        = TypeParam{3.01332e-11};
  TypeParam qi        = TypeParam{5.57166e-06};
  TypeParam qs        = TypeParam{3.55432e-05};
  TypeParam mi        = TypeParam{1e-09};
  TypeParam dvsw      = TypeParam{0.0};
  TypeParam dt        = TypeParam{30};
  TypeParam reference = TypeParam{5.1423851647153399e-17};

  TypeParam result    = transition::rain_to_graupel(t, rho, qc, qr, qi, qs, mi, dvsw, dt);
  this->validate(result, reference);
}

TYPED_TEST_P(HydrometeorTransitionTest, RainToVaporDefault) {
  TypeParam t      = TypeParam{258.542};
  TypeParam rho    = TypeParam{0.956089};
  TypeParam qc     = TypeParam{8.6157e-05};
  TypeParam qr     = TypeParam{3.01332e-11};
  TypeParam dvsw   = TypeParam{0.0};
  TypeParam dt     = TypeParam{30};

  TypeParam result = transition::rain_to_vapor(t, rho, qc, qr, dvsw, dt);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, RainToVapor) {
  TypeParam t                = TypeParam{258.542};
  TypeParam rho              = TypeParam{0.956089};
  TypeParam qc               = TypeParam{0.0};
  TypeParam qr               = TypeParam{3.01332e-11};
  TypeParam dvsw             = TypeParam{-1e-10};
  TypeParam dt               = TypeParam{30};
  TypeParam reference_float  = TypeParam{7.7282399e-17};
  TypeParam reference_double = TypeParam{7.7282338679917529e-17};

  TypeParam result           = transition::rain_to_vapor(t, rho, qc, qr, dvsw, dt);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(HydrometeorTransitionTest, SnowToGraupelDefault) {
  TypeParam t      = TypeParam{258.157};
  TypeParam rho    = TypeParam{0.93171};
  TypeParam qc     = TypeParam{0.0};
  TypeParam qs     = TypeParam{4.34854e-05};

  TypeParam result = transition::snow_to_graupel(t, rho, qc, qs);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, SnowToGraupel) {
  TypeParam t         = TypeParam{265.85};
  TypeParam rho       = TypeParam{1.04848};
  TypeParam qc        = TypeParam{7.02792e-05};
  TypeParam qs        = TypeParam{4.44664e-07};
  TypeParam reference = TypeParam{6.2696154545048011e-10};

  TypeParam result    = transition::snow_to_graupel(t, rho, qc, qs);
  this->validate(result, reference);
}

TYPED_TEST_P(HydrometeorTransitionTest, SnowToRainDefault) {
  TypeParam t      = TypeParam{265.83};
  TypeParam p      = TypeParam{80134.5};
  TypeParam rho    = TypeParam{1.04892};
  TypeParam dvsw0  = TypeParam{-0.00258631};
  TypeParam qs     = TypeParam{1.47687e-06};

  TypeParam result = transition::snow_to_rain(t, p, rho, dvsw0, qs);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, SnowToRain) {
  TypeParam t                = TypeParam{275.83};
  TypeParam p                = TypeParam{80134.5};
  TypeParam rho              = TypeParam{1.04892};
  TypeParam dvsw0            = TypeParam{0.00258631};
  TypeParam qs               = TypeParam{1.47687e-06};
  TypeParam reference_float  = TypeParam{3.7268515e-07};
  TypeParam reference_double = TypeParam{3.7268547760462804e-07};

  TypeParam result           = transition::snow_to_rain(t, p, rho, dvsw0, qs);
  this->validate(result, reference_float, reference_double);
}

TYPED_TEST_P(HydrometeorTransitionTest, VaporXGraupelDefault) {
  TypeParam t      = TypeParam{278.026};
  TypeParam p      = TypeParam{95987.1};
  TypeParam rho    = TypeParam{1.20041};
  TypeParam qg     = TypeParam{2.05496e-16};
  TypeParam dvsw   = TypeParam{-0.00234674};
  TypeParam dvsi   = TypeParam{-0.00261576};
  TypeParam dvsw0  = TypeParam{-0.00076851};
  TypeParam dt     = TypeParam{30};

  TypeParam result = transition::vapor_x_graupel(t, p, rho, qg, dvsw, dvsi, dvsw0, dt);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, VaporXGraupel) {
  TypeParam t         = TypeParam{278.026};
  TypeParam p         = TypeParam{95987.1};
  TypeParam rho       = TypeParam{1.20041};
  TypeParam qg        = TypeParam{2.05496e-11};
  TypeParam dvsw      = TypeParam{-0.00234674};
  TypeParam dvsi      = TypeParam{-0.00261576};
  TypeParam dvsw0     = TypeParam{-0.00076851};
  TypeParam dt        = TypeParam{30};
  TypeParam reference = TypeParam{-6.8498666666666675e-13};

  TypeParam result    = transition::vapor_x_graupel(t, p, rho, qg, dvsw, dvsi, dvsw0, dt);
  this->validate(result, reference);
}

TYPED_TEST_P(HydrometeorTransitionTest, VaporXSnowDefault) {
  TypeParam t       = TypeParam{278.748};
  TypeParam p       = TypeParam{95995.5};
  TypeParam rho     = TypeParam{1.19691};
  TypeParam qs      = TypeParam{1.25653e-20};
  TypeParam ns      = TypeParam{800000};
  TypeParam lambda  = TypeParam{1e+10};
  TypeParam eta     = TypeParam{0.0};
  TypeParam ice_dep = TypeParam{0.0};
  TypeParam dvsw    = TypeParam{-0.00196781};
  TypeParam dvsi    = TypeParam{-0.00229367};
  TypeParam dvsw0   = TypeParam{-0.000110022};
  TypeParam dt      = TypeParam{30};

  TypeParam result  = transition::vapor_x_snow(t, p, rho, qs, ns, lambda, eta, ice_dep, dvsw, dvsi, dvsw0, dt);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, VaporXSnow) {
  TypeParam t         = TypeParam{278.748};
  TypeParam p         = TypeParam{95995.5};
  TypeParam rho       = TypeParam{1.19691};
  TypeParam qs        = TypeParam{1.25653e-10};
  TypeParam ns        = TypeParam{800000};
  TypeParam lambda    = TypeParam{1e+10};
  TypeParam eta       = TypeParam{0.0};
  TypeParam ice_dep   = TypeParam{0.0};
  TypeParam dvsw      = TypeParam{-0.00196781};
  TypeParam dvsi      = TypeParam{-0.00229367};
  TypeParam dvsw0     = TypeParam{-0.000110022};
  TypeParam dt        = TypeParam{30};
  TypeParam reference = TypeParam{-8.6584296264775935e-13};

  TypeParam result    = transition::vapor_x_snow(t, p, rho, qs, ns, lambda, eta, ice_dep, dvsw, dvsi, dvsw0, dt);
  this->validate(result, reference);
}

TYPED_TEST_P(HydrometeorTransitionTest, VaporXIceDefault) {
  TypeParam qi     = TypeParam{2.02422e-23};
  TypeParam mi     = TypeParam{1e-12};
  TypeParam eta    = TypeParam{1.32343e-05};
  TypeParam dvsi   = TypeParam{-0.000618828};
  TypeParam rho    = TypeParam{1.19691};
  TypeParam dt     = TypeParam{30};

  TypeParam result = transition::vapor_x_ice(qi, mi, eta, dvsi, rho, dt);
  this->validate(result, ZERO<TypeParam>);
}

TYPED_TEST_P(HydrometeorTransitionTest, VaporXIce) {
  TypeParam qi               = TypeParam{9.53048e-07};
  TypeParam mi               = TypeParam{1e-09};
  TypeParam eta              = TypeParam{1.90278e-05};
  TypeParam dvsi             = TypeParam{0.000120375};
  TypeParam rho              = TypeParam{1.19691};
  TypeParam dt               = TypeParam{30};
  TypeParam reference_float  = TypeParam{2.210617e-09};
  TypeParam reference_double = TypeParam{2.2106162342610385e-09};

  TypeParam result           = transition::vapor_x_ice(qi, mi, eta, dvsi, rho, dt);
  this->validate(result, reference_float, reference_double);
}

REGISTER_TYPED_TEST_SUITE_P(HydrometeorTransitionTest, CloudToGraupelDefault, CloudToGraupel, CloudToRainDefault,
                            CloudToRain, CloudToSnowDefault, CloudToSnow, CloudXIceDefault, CloudXIceI, CloudXIceC,
                            GraupelToRainDefault, GraupelToRain, IceToGraupelDefault, IceToGraupelR, IceToGraupelG,
                            IceToSnowDefault, IceToSnow, RainToGraupelDefault, RainToGraupel1, RainToGraupel2,
                            RainToGraupel3, RainToVaporDefault, RainToVapor, SnowToGraupelDefault, SnowToGraupel,
                            SnowToRainDefault, SnowToRain, VaporXGraupelDefault, VaporXGraupel, VaporXSnowDefault,
                            VaporXSnow, VaporXIceDefault, VaporXIce);
using MyTypes = ::testing::Types<float, double>;
INSTANTIATE_TYPED_TEST_SUITE_P(RagnarokTests, HydrometeorTransitionTest, MyTypes);

}  // namespace
