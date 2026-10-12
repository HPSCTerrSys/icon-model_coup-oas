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
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <Kokkos_Core.hpp>
#include <cstdio>

#include "aes_microphysics/common/utils.hpp"
#include "aes_microphysics/graupel.hpp"
#include "aes_thermodynamics/thermo.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

typedef int Integer_t;

bool is_initialized = false;

template <typename T>
static auto get_ptr(py::array_t<T> a) {
  py::buffer_info xinfo = a.request();
  return static_cast<T*>(xinfo.ptr);
}

template <typename T>
static auto full_size(py::array_t<T> a) {
  py::buffer_info xinfo = a.request();
  int n                 = 1;
  for (auto r : xinfo.shape) {
    n *= r;
  }
  return n;
}

template <typename T>
void calc_dz(py::array_t<T> z, py::array_t<T> dz, const Integer_t ncells, const Integer_t nlev) {
  // Convert numpy array to std::vector for input
  auto z_ptr = get_ptr(z);
  const std::vector<T> z_vec(z_ptr, z_ptr + ncells * nlev);

  // Create output vector
  std::vector<T> dz_vec;

  // Call the actual calc_dz function
  utils::calc_dz<T>(z_vec, dz_vec, ncells, nlev);

  // Copy result back to numpy array
  auto dz_ptr = get_ptr(dz);
  assert(full_size(dz) == ncells * nlev);
  assert(dz_vec.size() == ncells * nlev);
  std::memcpy(dz_ptr, dz_vec.data(), ncells * nlev * sizeof(T));
}

template <typename T>
void saturation_adjustment(const Integer_t ncells, const Integer_t nlev, py::array_t<T> ta, py::array_t<T> qv,
                           py::array_t<T> qc, py::array_t<T> qr, py::array_t<T> total_ice, py::array_t<T> rho) {
  Integer_t nproma = ncells;
  Integer_t kstart = 0;
  Integer_t kend   = nlev;
  Integer_t jcs    = 0;
  Integer_t jce    = ncells;

  assert(full_size(ta) == ncells * nlev);
  thermo::saturation_adjustment_2d<T>(nproma, nlev, kstart, kend, jcs, jce, get_ptr(ta), get_ptr(qv), get_ptr(qc),
                                      get_ptr(qr), get_ptr(total_ice), get_ptr(rho));
}

void finalize() {}

void initialize() {
  if (is_initialized) return;
  Kokkos::initialize();
  std::atexit(Kokkos::finalize);
  is_initialized = true;
}

template <typename T>
void run(const Integer_t ncells, const Integer_t nlev, const T dt, py::array_t<T> dz, py::array_t<T> t,
         py::array_t<T> rho, py::array_t<T> p, py::array_t<T> qv, py::array_t<T> qc, py::array_t<T> qi,
         py::array_t<T> qr, py::array_t<T> qs, py::array_t<T> qg, const T qnc, py::array_t<T> prr_gsp,
         py::array_t<T> pri_gsp, py::array_t<T> prs_gsp, py::array_t<T> prg_gsp, py::array_t<T> pre_gsp,
         py::array_t<T> pflx, const bool lrain) {
  assert(is_initialized);
  Integer_t nvec    = ncells;
  Integer_t kend    = nlev;
  Integer_t ivstart = 0;
  Integer_t ivend   = ncells;
  Integer_t kstart  = 0;
  T qnc_vec[ncells];
  qnc_vec[0]     = T{qnc};
  const T cia    = T{1.0};

  auto execSpace = Kokkos::DefaultExecutionSpace();
  graupel::run<decltype(execSpace), T>(execSpace, nvec, kend, ivstart, ivend, kstart, dt, cia, get_ptr(dz), get_ptr(t),
                                       get_ptr(rho), get_ptr(p), get_ptr(qv), get_ptr(qc), get_ptr(qi), get_ptr(qr),
                                       get_ptr(qs), get_ptr(qg), qnc_vec, get_ptr(prr_gsp), get_ptr(pri_gsp),
                                       get_ptr(prs_gsp), get_ptr(prg_gsp), get_ptr(pre_gsp), get_ptr(pflx), lrain);
}

PYBIND11_MODULE(aes_muphys_py, m) {
  m.doc() = "pybind11 AES microphysics plugin";
  m.def("calc_dz", &calc_dz<float>, "utils::calc_dz", "z"_a, "dz"_a, "ncells"_a, "nlev"_a);
  m.def("calc_dz", &calc_dz<double>, "utils::calc_dz", "z"_a, "dz"_a, "ncells"_a, "nlev"_a);
  m.def("saturation_adjustment", &saturation_adjustment<float>, "saturation_adjustment_2d", "ncells"_a, "nlev"_a,
        "ta"_a, "qv"_a, "qc"_a, "qr"_a, "total_ice"_a, "rho"_a);
  m.def("saturation_adjustment", &saturation_adjustment<double>, "saturation_adjustment_2d", "ncells"_a, "nlev"_a,
        "ta"_a, "qv"_a, "qc"_a, "qr"_a, "total_ice"_a, "rho"_a);
  m.def("initialize", &initialize);
  m.def("finalize", &finalize);
  m.def("run", &run<float>, "graupel::run", "ncells"_a, "nlev"_a, "dt"_a, "dz"_a, "t"_a, "rho"_a, "p"_a, "qv"_a, "qc"_a,
        "qi"_a, "qr"_a, "qs"_a, "qg"_a, "qnc"_a, "prr_gsp"_a, "pri_gsp"_a, "prs_gsp"_a, "prg_gsp"_a, "pre_gsp"_a,
        "pflx"_a, "lrain"_a);
  m.def("run", &run<double>, "graupel::run", "ncells"_a, "nlev"_a, "dt"_a, "dz"_a, "t"_a, "rho"_a, "p"_a, "qv"_a,
        "qc"_a, "qi"_a, "qr"_a, "qs"_a, "qg"_a, "qnc"_a, "prr_gsp"_a, "pri_gsp"_a, "prs_gsp"_a, "prg_gsp"_a,
        "pre_gsp"_a, "pflx"_a, "lrain"_a);
}
