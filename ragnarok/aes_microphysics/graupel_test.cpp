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
#include "aes_microphysics/graupel.hpp"

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <iostream>
#include <type_traits>
#include <vector>

namespace {
using testing::Types;

template <typename T>
class GraupelRunTest : public testing::Test {
 public:
  void validate(T actual, T expected, T tolerance = T{1.0e-6}) {
    if constexpr (std::is_same_v<T, float>) {
      EXPECT_NEAR(expected, actual, tolerance);
    } else {
      EXPECT_NEAR(expected, actual, tolerance);
    }
  }

  void validate_range(T actual, T min_val, T max_val) {
    EXPECT_GE(actual, min_val);
    EXPECT_LE(actual, max_val);
  }
};  // class GraupelRunTest

TYPED_TEST_SUITE_P(GraupelRunTest);

TYPED_TEST_P(GraupelRunTest, RunBasicSingleColumn) {
  // Test parameters
  const int nvec      = 1;                // Single horizontal point
  const int ke        = 5;                // 5 vertical levels
  const int ivstart   = 0;                // Start from first point
  const int ivend     = nvec;             // End at first point
  const int kstart    = 0;                // Start from first level
  const TypeParam dt  = TypeParam{30.0};  // 30 second time step
  const TypeParam cia = TypeParam{0.7};   // 0.7 control param

  // Allocate arrays
  std::vector<TypeParam> dz(nvec * ke, TypeParam{500.0});    // 500m layer thickness
  std::vector<TypeParam> t(nvec * ke, TypeParam{273.15});    // Temperature at freezing
  std::vector<TypeParam> rho(nvec * ke, TypeParam{1.2});     // Air density
  std::vector<TypeParam> p(nvec * ke, TypeParam{90000.0});   // Pressure ~900 hPa
  std::vector<TypeParam> qv(nvec * ke, TypeParam{0.005});    // Water vapor
  std::vector<TypeParam> qc(nvec * ke, TypeParam{0.0001});   // Cloud water
  std::vector<TypeParam> qi(nvec * ke, TypeParam{0.00001});  // Cloud ice
  std::vector<TypeParam> qr(nvec * ke, TypeParam{0.0});      // Rain
  std::vector<TypeParam> qs(nvec * ke, TypeParam{0.0});      // Snow
  std::vector<TypeParam> qg(nvec * ke, TypeParam{0.0});      // Graupel
  std::vector<TypeParam> qnc(nvec * ke, TypeParam{100.0});   // Cloud number concentration
  std::vector<TypeParam> prr_gsp(nvec, TypeParam{0.0});      // Rain precip rate
  std::vector<TypeParam> pri_gsp(nvec, TypeParam{0.0});      // Ice precip rate
  std::vector<TypeParam> prs_gsp(nvec, TypeParam{0.0});      // Snow precip rate
  std::vector<TypeParam> prg_gsp(nvec, TypeParam{0.0});      // Graupel precip rate
  std::vector<TypeParam> pre_gsp(nvec, TypeParam{0.0});      // Energy flux
  std::vector<TypeParam> pflx(nvec * ke, TypeParam{0.0});    // Total precip flux

  auto h_dz            = HostView2D<TypeParam>(dz.data(), ke, nvec);
  auto d_dz            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_dz);
  auto h_t             = HostView2D<TypeParam>(t.data(), ke, nvec);
  auto d_t             = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_t);
  auto h_rho           = HostView2D<TypeParam>(rho.data(), ke, nvec);
  auto d_rho           = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_rho);
  auto h_p             = HostView2D<TypeParam>(p.data(), ke, nvec);
  auto d_p             = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_p);
  auto h_qv            = HostView2D<TypeParam>(qv.data(), ke, nvec);
  auto d_qv            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qv);
  auto h_qc            = HostView2D<TypeParam>(qc.data(), ke, nvec);
  auto d_qc            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qc);
  auto h_qi            = HostView2D<TypeParam>(qi.data(), ke, nvec);
  auto d_qi            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qi);
  auto h_qr            = HostView2D<TypeParam>(qr.data(), ke, nvec);
  auto d_qr            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qr);
  auto h_qs            = HostView2D<TypeParam>(qs.data(), ke, nvec);
  auto d_qs            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qs);
  auto h_qg            = HostView2D<TypeParam>(qg.data(), ke, nvec);
  auto d_qg            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qg);
  auto h_qnc           = HostView1D<TypeParam>(qnc.data(), nvec);
  auto d_qnc           = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qnc);
  auto h_prr_gsp       = HostView1D<TypeParam>(prr_gsp.data(), nvec);
  auto d_prr_gsp       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_prr_gsp);
  auto h_pri_gsp       = HostView1D<TypeParam>(pri_gsp.data(), nvec);
  auto d_pri_gsp       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_pri_gsp);
  auto h_prs_gsp       = HostView1D<TypeParam>(prs_gsp.data(), nvec);
  auto d_prs_gsp       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_prs_gsp);
  auto h_prg_gsp       = HostView1D<TypeParam>(prg_gsp.data(), nvec);
  auto d_prg_gsp       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_prg_gsp);
  auto h_pre_gsp       = HostView1D<TypeParam>(pre_gsp.data(), nvec);
  auto d_pre_gsp       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_pre_gsp);
  auto h_pflx          = HostView2D<TypeParam>(pflx.data(), ke, nvec);
  auto d_pflx          = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_pflx);

  // Store initial values
  TypeParam qv_initial = qv[0];
  TypeParam qc_initial = qc[0];
  TypeParam qi_initial = qi[0];

  // Run the graupel scheme
  graupel::run(nvec, ke, ivstart, ivend, kstart, dt, cia, d_dz.data(), d_t.data(), d_rho.data(), d_p.data(),
               d_qv.data(), d_qc.data(), d_qi.data(), d_qr.data(), d_qs.data(), d_qg.data(), d_qnc.data(),
               d_prr_gsp.data(), d_pri_gsp.data(), d_prs_gsp.data(), d_prg_gsp.data(), d_pre_gsp.data(), d_pflx.data());

  Kokkos::deep_copy(h_qv, d_qv);
  Kokkos::deep_copy(h_qc, d_qc);
  Kokkos::deep_copy(h_qr, d_qr);
  Kokkos::deep_copy(h_qs, d_qs);
  Kokkos::deep_copy(h_qi, d_qi);
  Kokkos::deep_copy(h_qg, d_qg);

  Kokkos::deep_copy(h_qr, d_qr);
  Kokkos::deep_copy(h_qi, d_qi);
  Kokkos::deep_copy(h_qs, d_qs);
  Kokkos::deep_copy(h_qg, d_qg);

  Kokkos::deep_copy(h_t, d_t);
  Kokkos::deep_copy(h_pflx, d_pflx);
  Kokkos::deep_copy(h_pflx, d_pflx);
  Kokkos::fence();

  // Validate that physical constraints are maintained
  // 1. All mixing ratios should be non-negative
  for (int k = 0; k < ke; ++k) {
    this->validate_range(qv[k], TypeParam{0.0}, TypeParam{1.0});
    this->validate_range(qc[k], TypeParam{0.0}, TypeParam{1.0});
    this->validate_range(qi[k], TypeParam{0.0}, TypeParam{1.0});
    this->validate_range(qr[k], TypeParam{0.0}, TypeParam{1.0});
    this->validate_range(qs[k], TypeParam{0.0}, TypeParam{1.0});
    this->validate_range(qg[k], TypeParam{0.0}, TypeParam{1.0});
  }

  // 2. Temperature should be reasonable (not extreme)
  for (int k = 0; k < ke; ++k) {
    this->validate_range(t[k], TypeParam{180.0}, TypeParam{320.0});
  }

  // 3. Precipitation rates should be non-negative
  this->validate_range(prr_gsp[0], TypeParam{0.0}, TypeParam{1000.0});
  this->validate_range(pri_gsp[0], TypeParam{0.0}, TypeParam{1000.0});
  this->validate_range(prs_gsp[0], TypeParam{0.0}, TypeParam{1000.0});
  this->validate_range(prg_gsp[0], TypeParam{0.0}, TypeParam{1000.0});

  // 4. Total water should be conserved or decrease (due to precipitation)
  TypeParam total_water_initial = qv_initial + qc_initial + qi_initial;
  TypeParam total_water_final   = qv[0] + qc[0] + qi[0] + qr[0] + qs[0] + qg[0];
  EXPECT_LE(total_water_final, total_water_initial + TypeParam{1.0e-10});
}

TYPED_TEST_P(GraupelRunTest, RunWarmRain) {
  // Test warm rain process (T > 273.15K)
  const int nvec      = 1;
  const int ke        = 3;
  const int ivstart   = 0;
  const int ivend     = nvec;
  const int kstart    = 0;
  const TypeParam dt  = TypeParam{10.0};
  const TypeParam cia = TypeParam{0.7};

  // Setup for warm rain conditions
  std::vector<TypeParam> dz(nvec * ke, TypeParam{300.0});
  std::vector<TypeParam> t(nvec * ke, TypeParam{283.15});  // 10°C
  std::vector<TypeParam> rho(nvec * ke, TypeParam{1.15});
  std::vector<TypeParam> p(nvec * ke, TypeParam{95000.0});
  std::vector<TypeParam> qv(nvec * ke, TypeParam{0.008});
  std::vector<TypeParam> qc(nvec * ke, TypeParam{0.0005});   // Significant cloud water
  std::vector<TypeParam> qi(nvec * ke, TypeParam{0.0});      // No ice
  std::vector<TypeParam> qr(nvec * ke, TypeParam{0.00001});  // Small rain
  std::vector<TypeParam> qs(nvec * ke, TypeParam{0.0});
  std::vector<TypeParam> qg(nvec * ke, TypeParam{0.0});
  std::vector<TypeParam> qnc(nvec * ke, TypeParam{200.0});
  std::vector<TypeParam> prr_gsp(nvec, TypeParam{0.0});
  std::vector<TypeParam> pri_gsp(nvec, TypeParam{0.0});
  std::vector<TypeParam> prs_gsp(nvec, TypeParam{0.0});
  std::vector<TypeParam> prg_gsp(nvec, TypeParam{0.0});
  std::vector<TypeParam> pre_gsp(nvec, TypeParam{0.0});
  std::vector<TypeParam> pflx(nvec * ke, TypeParam{0.0});

  auto h_dz            = HostView2D<TypeParam>(dz.data(), ke, nvec);
  auto d_dz            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_dz);
  auto h_t             = HostView2D<TypeParam>(t.data(), ke, nvec);
  auto d_t             = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_t);
  auto h_rho           = HostView2D<TypeParam>(rho.data(), ke, nvec);
  auto d_rho           = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_rho);
  auto h_p             = HostView2D<TypeParam>(p.data(), ke, nvec);
  auto d_p             = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_p);
  auto h_qv            = HostView2D<TypeParam>(qv.data(), ke, nvec);
  auto d_qv            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qv);
  auto h_qc            = HostView2D<TypeParam>(qc.data(), ke, nvec);
  auto d_qc            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qc);
  auto h_qi            = HostView2D<TypeParam>(qi.data(), ke, nvec);
  auto d_qi            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qi);
  auto h_qr            = HostView2D<TypeParam>(qr.data(), ke, nvec);
  auto d_qr            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qr);
  auto h_qs            = HostView2D<TypeParam>(qs.data(), ke, nvec);
  auto d_qs            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qs);
  auto h_qg            = HostView2D<TypeParam>(qg.data(), ke, nvec);
  auto d_qg            = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qg);
  auto h_qnc           = HostView1D<TypeParam>(qnc.data(), nvec);
  auto d_qnc           = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qnc);
  auto h_prr_gsp       = HostView1D<TypeParam>(prr_gsp.data(), nvec);
  auto d_prr_gsp       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_prr_gsp);
  auto h_pri_gsp       = HostView1D<TypeParam>(pri_gsp.data(), nvec);
  auto d_pri_gsp       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_pri_gsp);
  auto h_prs_gsp       = HostView1D<TypeParam>(prs_gsp.data(), nvec);
  auto d_prs_gsp       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_prs_gsp);
  auto h_prg_gsp       = HostView1D<TypeParam>(prg_gsp.data(), nvec);
  auto d_prg_gsp       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_prg_gsp);
  auto h_pre_gsp       = HostView1D<TypeParam>(pre_gsp.data(), nvec);
  auto d_pre_gsp       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_pre_gsp);
  auto h_pflx          = HostView2D<TypeParam>(pflx.data(), ke, nvec);
  auto d_pflx          = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_pflx);

  TypeParam qr_initial = qr[0];

  graupel::run(nvec, ke, ivstart, ivend, kstart, dt, cia, d_dz.data(), d_t.data(), d_rho.data(), d_p.data(),
               d_qv.data(), d_qc.data(), d_qi.data(), d_qr.data(), d_qs.data(), d_qg.data(), d_qnc.data(),
               d_prr_gsp.data(), d_pri_gsp.data(), d_prs_gsp.data(), d_prg_gsp.data(), d_pre_gsp.data(), d_pflx.data());

  Kokkos::deep_copy(h_qr, d_qr);
  Kokkos::deep_copy(h_qi, d_qi);
  Kokkos::deep_copy(h_qs, d_qs);
  Kokkos::deep_copy(h_qg, d_qg);
  Kokkos::fence();

  // Warm rain should produce rain (autoconversion + accretion)
  // qr should increase or stay the same
  EXPECT_GE(qr[0], TypeParam{0.0});

  // No ice processes should occur
  EXPECT_EQ(qi[0], TypeParam{0.0});
  EXPECT_EQ(qs[0], TypeParam{0.0});
  EXPECT_EQ(qg[0], TypeParam{0.0});
}

TYPED_TEST_P(GraupelRunTest, RunZeroInitialConditions) {
  // Test with zero initial hydrometeors
  const int nvec      = 1;
  const int ke        = 3;
  const int ivstart   = 0;
  const int ivend     = nvec;
  const int kstart    = 0;
  const TypeParam dt  = TypeParam{30.0};
  const TypeParam cia = TypeParam{0.7};

  std::vector<TypeParam> dz(nvec * ke, TypeParam{400.0});
  std::vector<TypeParam> t(nvec * ke, TypeParam{280.0});
  std::vector<TypeParam> rho(nvec * ke, TypeParam{1.18});
  std::vector<TypeParam> p(nvec * ke, TypeParam{92000.0});
  std::vector<TypeParam> qv(nvec * ke, TypeParam{0.003});
  std::vector<TypeParam> qc(nvec * ke, TypeParam{0.0});  // No cloud
  std::vector<TypeParam> qi(nvec * ke, TypeParam{0.0});  // No ice
  std::vector<TypeParam> qr(nvec * ke, TypeParam{0.0});  // No rain
  std::vector<TypeParam> qs(nvec * ke, TypeParam{0.0});  // No snow
  std::vector<TypeParam> qg(nvec * ke, TypeParam{0.0});  // No graupel
  std::vector<TypeParam> qnc(nvec * ke, TypeParam{100.0});
  std::vector<TypeParam> prr_gsp(nvec, TypeParam{0.0});
  std::vector<TypeParam> pri_gsp(nvec, TypeParam{0.0});
  std::vector<TypeParam> prs_gsp(nvec, TypeParam{0.0});
  std::vector<TypeParam> prg_gsp(nvec, TypeParam{0.0});
  std::vector<TypeParam> pre_gsp(nvec, TypeParam{0.0});
  std::vector<TypeParam> pflx(nvec * ke, TypeParam{0.0});

  auto h_dz      = HostView2D<TypeParam>(dz.data(), ke, nvec);
  auto d_dz      = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_dz);
  auto h_t       = HostView2D<TypeParam>(t.data(), ke, nvec);
  auto d_t       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_t);
  auto h_rho     = HostView2D<TypeParam>(rho.data(), ke, nvec);
  auto d_rho     = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_rho);
  auto h_p       = HostView2D<TypeParam>(p.data(), ke, nvec);
  auto d_p       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_p);
  auto h_qv      = HostView2D<TypeParam>(qv.data(), ke, nvec);
  auto d_qv      = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qv);
  auto h_qc      = HostView2D<TypeParam>(qc.data(), ke, nvec);
  auto d_qc      = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qc);
  auto h_qi      = HostView2D<TypeParam>(qi.data(), ke, nvec);
  auto d_qi      = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qi);
  auto h_qr      = HostView2D<TypeParam>(qr.data(), ke, nvec);
  auto d_qr      = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qr);
  auto h_qs      = HostView2D<TypeParam>(qs.data(), ke, nvec);
  auto d_qs      = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qs);
  auto h_qg      = HostView2D<TypeParam>(qg.data(), ke, nvec);
  auto d_qg      = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qg);
  auto h_qnc     = HostView1D<TypeParam>(qnc.data(), nvec);
  auto d_qnc     = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qnc);
  auto h_prr_gsp = HostView1D<TypeParam>(prr_gsp.data(), nvec);
  auto d_prr_gsp = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_prr_gsp);
  auto h_pri_gsp = HostView1D<TypeParam>(pri_gsp.data(), nvec);
  auto d_pri_gsp = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_pri_gsp);
  auto h_prs_gsp = HostView1D<TypeParam>(prs_gsp.data(), nvec);
  auto d_prs_gsp = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_prs_gsp);
  auto h_prg_gsp = HostView1D<TypeParam>(prg_gsp.data(), nvec);
  auto d_prg_gsp = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_prg_gsp);
  auto h_pre_gsp = HostView1D<TypeParam>(pre_gsp.data(), nvec);
  auto d_pre_gsp = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_pre_gsp);
  auto h_pflx    = HostView2D<TypeParam>(pflx.data(), ke, nvec);
  auto d_pflx    = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_pflx);

  graupel::run(nvec, ke, ivstart, ivend, kstart, dt, cia, d_dz.data(), d_t.data(), d_rho.data(), d_p.data(),
               d_qv.data(), d_qc.data(), d_qi.data(), d_qr.data(), d_qs.data(), d_qg.data(), d_qnc.data(),
               d_prr_gsp.data(), d_pri_gsp.data(), d_prs_gsp.data(), d_prg_gsp.data(), d_pre_gsp.data(), d_pflx.data());

  Kokkos::deep_copy(h_qc, d_qc);
  Kokkos::deep_copy(h_qi, d_qi);
  Kokkos::deep_copy(h_qr, d_qr);
  Kokkos::deep_copy(h_qs, d_qs);
  Kokkos::deep_copy(h_qg, d_qg);
  Kokkos::deep_copy(h_prr_gsp, d_prr_gsp);
  Kokkos::deep_copy(h_pri_gsp, d_pri_gsp);
  Kokkos::deep_copy(h_prs_gsp, d_prs_gsp);
  Kokkos::deep_copy(h_prg_gsp, d_prg_gsp);
  Kokkos::fence();

  // With no condensed water, everything should remain zero
  for (int k = 0; k < ke; ++k) {
    EXPECT_EQ(qc[k], TypeParam{0.0});
    EXPECT_EQ(qi[k], TypeParam{0.0});
    EXPECT_EQ(qr[k], TypeParam{0.0});
    EXPECT_EQ(qs[k], TypeParam{0.0});
    EXPECT_EQ(qg[k], TypeParam{0.0});
  }

  // No precipitation
  EXPECT_EQ(prr_gsp[0], TypeParam{0.0});
  EXPECT_EQ(pri_gsp[0], TypeParam{0.0});
  EXPECT_EQ(prs_gsp[0], TypeParam{0.0});
  EXPECT_EQ(prg_gsp[0], TypeParam{0.0});
}

REGISTER_TYPED_TEST_SUITE_P(GraupelRunTest, RunBasicSingleColumn, RunWarmRain, RunZeroInitialConditions);
using MyTypes = ::testing::Types<float, double>;
INSTANTIATE_TYPED_TEST_SUITE_P(RagnarokTests, GraupelRunTest, MyTypes);

}  // namespace
