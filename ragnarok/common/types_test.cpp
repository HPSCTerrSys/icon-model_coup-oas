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
#include "types.hpp"

#include <gtest/gtest.h>

namespace {
using testing::Types;

template <typename T>
class TypesTest : public testing::Test {
 public:
  // workaround for nvcc + kokkos + gtest bug: cannot dispatch a KOKKOS lambda from a
  // non-public parent function when compiled with nvcc
  static void update_view_2d(ManagedView2D<T> view2d, size_t dim1, size_t dim2) {
    using RangePolicy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;

    Kokkos::parallel_for(
        "init2d", RangePolicy2D({0, 0}, {dim1, dim2}),
        KOKKOS_LAMBDA(const size_t i, const size_t j) { view2d(i, j) = static_cast<T>(i * dim2 + j); });
    Kokkos::fence();
  }

  // workaround for nvcc + kokkos + gtest bug: cannot dispatch a KOKKOS lambda from a
  // non-public parent function when compiled with nvcc
  static void update_view_3d(ManagedView3D<T> view3d, size_t dim1, size_t dim2, size_t dim3) {
    using RangePolicy3D = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

    Kokkos::parallel_for(
        "init3d", RangePolicy3D({0, 0, 0}, {dim1, dim2, dim3}),
        KOKKOS_LAMBDA(const size_t i, const size_t j, const size_t k) {
          view3d(i, j, k) = static_cast<T>(i * dim2 + j * dim3 + k);
        });
    Kokkos::fence();
  }

};  // end class

TYPED_TEST_SUITE_P(TypesTest);

TYPED_TEST_P(TypesTest, CheckManagedView2DLayout) {
  const size_t dim1 = 2;
  const size_t dim2 = 3;

  ManagedView2D<TypeParam> view2d("2d", dim1, dim2);
  Kokkos::View<TypeParam**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryManaged> hostView("host", dim1, dim2);
  TypesTest<TypeParam>::update_view_2d(view2d, dim1, dim2);
  Kokkos::deep_copy(hostView, view2d);

  // check the layout setups
  EXPECT_EQ(hostView(0, 2), TypeParam{2});
}

TYPED_TEST_P(TypesTest, CheckManagedView3DLayout) {
  const size_t dim1 = 2;
  const size_t dim2 = 3;
  const size_t dim3 = 2;

  ManagedView3D<TypeParam> view3d("3d", dim1, dim2, dim3);
  Kokkos::View<TypeParam***, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryManaged> hostView("host", dim1, dim2,
                                                                                                     dim3);
  TypesTest<TypeParam>::update_view_3d(view3d, dim1, dim2, dim3);
  Kokkos::deep_copy(hostView, view3d);

  EXPECT_EQ(hostView(0, 1, 1), TypeParam{3});
}

TYPED_TEST_P(TypesTest, CheckUnmanagedView2DLayout) {
  const size_t dim1              = 2;
  const size_t dim2              = 3;

  // emulate ICON's action to allocate the pointer in device memory
  TypeParam array2D[dim1 * dim2] = {TypeParam{0}, TypeParam{1}, TypeParam{2}, TypeParam{3}, TypeParam{4}, TypeParam{5}};
  ManagedView2D<TypeParam> icon2d("2d", dim1, dim2);
  TypesTest<TypeParam>::update_view_2d(icon2d, dim1, dim2);

  // create the view based on the existing pointer in the execution space
  auto view2d = View2D<TypeParam>(icon2d.data(), dim1, dim2);

  // expect reuse of the pointer in the same memory space
  EXPECT_EQ(icon2d.data(), view2d.data());
}

TYPED_TEST_P(TypesTest, CheckUnmanagedView3DLayout) {
  const size_t dim1                     = 2;
  const size_t dim2                     = 3;
  const size_t dim3                     = 2;

  // emulate ICON's action to allocate the pointer in device memory
  TypeParam array3D[dim1 * dim2 * dim3] = {TypeParam{0}, TypeParam{1}, TypeParam{2},  TypeParam{3},
                                           TypeParam{4}, TypeParam{5}, TypeParam{6},  TypeParam{7},
                                           TypeParam{8}, TypeParam{9}, TypeParam{10}, TypeParam{11}};
  ManagedView3D<TypeParam> icon3d("3d", dim1, dim2, dim3);
  TypesTest<TypeParam>::update_view_3d(icon3d, dim1, dim2, dim3);

  // create the view based on the existing pointer in the execution space
  auto view3d = View3D<TypeParam>(icon3d.data(), dim1, dim2, dim3);

  // expect reuse of the pointer in the same memory space
  EXPECT_EQ(icon3d.data(), view3d.data());
}

TYPED_TEST_P(TypesTest, CheckHostView) {
  const size_t dim1       = 2;
  TypeParam array1D[dim1] = {TypeParam{0}, TypeParam{1}};
  auto view1d             = HostView1D<TypeParam>(array1D, dim1);

  // expect always the same pointer, regardless of the Kokkos backend
  EXPECT_TRUE(array1D == view1d.data());
}

REGISTER_TYPED_TEST_SUITE_P(TypesTest, CheckManagedView2DLayout, CheckManagedView3DLayout, CheckUnmanagedView2DLayout,
                            CheckUnmanagedView3DLayout, CheckHostView);

using MyTypes = ::testing::Types<float, double>;
INSTANTIATE_TYPED_TEST_SUITE_P(RagnarokTests, TypesTest, MyTypes);

}  // namespace
