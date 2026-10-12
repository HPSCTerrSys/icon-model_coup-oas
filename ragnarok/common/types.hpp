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
#ifndef RAGNAROK_COMMON_TYPES_H_
#define RAGNAROK_COMMON_TYPES_H_

#include <Kokkos_Core.hpp>

using UnmanagedRestrictedTrait = Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Restrict>;
using ManagedRestrictedTrait   = Kokkos::MemoryTraits<Kokkos::Restrict>;

// unmanaged mutable views
template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
using View1D = Kokkos::View<T*, ExecutionSpace, UnmanagedRestrictedTrait>;

template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
using View2D = Kokkos::View<T**, Kokkos::LayoutRight, ExecutionSpace, UnmanagedRestrictedTrait>;

template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
using View3D = Kokkos::View<T***, Kokkos::LayoutRight, ExecutionSpace, UnmanagedRestrictedTrait>;

//  unmanaged immutable views
template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
using ConstView1D = Kokkos::View<const T*, ExecutionSpace, UnmanagedRestrictedTrait>;

template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
using ConstView2D = Kokkos::View<const T**, Kokkos::LayoutRight, ExecutionSpace, UnmanagedRestrictedTrait>;

template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
using ConstView3D = Kokkos::View<const T***, Kokkos::LayoutRight, ExecutionSpace, UnmanagedRestrictedTrait>;

//  unmanaged mutable host views
template <typename T>
using HostView1D = Kokkos::View<T*, Kokkos::HostSpace, UnmanagedRestrictedTrait>;

template <typename T>
using HostView2D = Kokkos::View<T**, Kokkos::LayoutRight, Kokkos::HostSpace, UnmanagedRestrictedTrait>;

template <typename T>
using HostView3D = Kokkos::View<T***, Kokkos::LayoutRight, Kokkos::HostSpace, UnmanagedRestrictedTrait>;

//  managed mutable views
template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
using ManagedView1D = Kokkos::View<T*, ExecutionSpace, ManagedRestrictedTrait>;

template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
using ManagedView2D = Kokkos::View<T**, Kokkos::LayoutRight, ExecutionSpace, ManagedRestrictedTrait>;

template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
using ManagedView3D = Kokkos::View<T***, Kokkos::LayoutRight, ExecutionSpace, ManagedRestrictedTrait>;

// Kokkos default memory space
using MemorySpace   = Kokkos::DefaultExecutionSpace::memory_space;

using RangePolicy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>;

#endif  // RAGNAROK_COMMON_TYPES_H_
