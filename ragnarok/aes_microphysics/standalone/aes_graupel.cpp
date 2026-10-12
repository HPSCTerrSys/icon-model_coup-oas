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
#include <Kokkos_Core.hpp>
#include <chrono>
#include <iostream>
#include <string>

#include "aes_microphysics/common/utils.hpp"
#include "aes_microphysics/graupel.hpp"
#include "aes_microphysics/standalone/io/netcdf_client.hpp"
#include "aes_thermodynamics/thermo.hpp"

void parse_args(std::string& file, std::string& precision, int argc, char** argv) {
  if (argc != 3) {
    std::cout << "Usage: ./aes_graupel <input-file> <precision>" << std::endl;
    std::cout << "E.g.:  ./aes_graupel input.nc single" << std::endl;
    exit(1);
  }

  file = argv[1];
  std::cout << "Input file: " << file << std::endl;

  precision = argv[2];
  if (precision.compare("single") != 0 && precision.compare("double") != 0) {
    std::cout << "Valid values for precision are: single and double" << std::endl;
    exit(1);
  }

  std::cout << "running in " << precision << " precision" << std::endl;
}

template <typename T, typename NCT>
void run_standalone(const std::string input_file) {
  // timers
  std::chrono::time_point<std::chrono::steady_clock> start_graupel, end_graupel;

  const std::string output_file = "output.nc";

  const T dt                    = T{30};
  const T cia                   = T{0.7};
  const int itime               = 0;

  std::cout << "Using default setup: dt=" << dt << ", cia=" << cia << ", itime=" << itime << std::endl;

  // Parameters from the input file
  int ncells, nlev;

  // Pre-calculated parameters and field data
  std::vector<T> dz;
  std::vector<T> z, t, p, rho, qv, qc, qi, qr, qs, qg;

  io::read_fields<T>(input_file, itime, ncells, nlev, z, t, p, rho, qv, qc, qi, qr, qs, qg);
  utils::calc_dz<T>(z, dz, ncells, nlev);

  const auto grid_size = nlev * ncells;
  std::vector<T> prr_gsp(ncells);
  std::vector<T> pri_gsp(ncells);
  std::vector<T> prs_gsp(ncells);
  std::vector<T> prg_gsp(ncells);
  std::vector<T> pflx(grid_size);
  std::vector<T> pre_gsp(ncells);
  std::vector<T> qnc(ncells, static_cast<T>(100));

  const int kbeg  = 0;
  const int kend  = nlev;
  const int ivbeg = 0;
  const int ivend = ncells;
  const int nvec  = ncells;

  // assume ICON scenario -> data is already allocated on the GPU when graupel is called
  auto h_dz       = HostView2D<T>(dz.data(), kend, nvec);
  auto d_dz       = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_dz);

  auto h_t        = HostView2D<T>(t.data(), kend, nvec);
  auto d_t        = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_t);

  auto h_rho      = HostView2D<T>(rho.data(), kend, nvec);
  auto d_rho      = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_rho);

  auto h_p        = HostView2D<T>(p.data(), kend, nvec);
  auto d_p        = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_p);

  auto h_pflx     = HostView2D<T>(pflx.data(), kend, nvec);
  auto d_pflx     = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_pflx);

  auto h_qx_lqc   = HostView2D<T>(qc.data(), kend, nvec);
  auto d_qx_lqc   = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qx_lqc);

  auto h_qx_lqi   = HostView2D<T>(qi.data(), kend, nvec);
  auto d_qx_lqi   = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qx_lqi);

  auto h_qx_lqr   = HostView2D<T>(qr.data(), kend, nvec);
  auto d_qx_lqr   = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qx_lqr);

  auto h_qx_lqs   = HostView2D<T>(qs.data(), kend, nvec);
  auto d_qx_lqs   = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qx_lqs);

  auto h_qx_lqg   = HostView2D<T>(qg.data(), kend, nvec);
  auto d_qx_lqg   = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qx_lqg);

  auto h_qx_lqv   = HostView2D<T>(qv.data(), kend, nvec);
  auto d_qx_lqv   = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qx_lqv);

  auto h_qp_lqi   = HostView1D<T>(pri_gsp.data(), nvec);
  auto d_qp_lqi   = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qp_lqi);

  auto h_qp_lqr   = HostView1D<T>(prr_gsp.data(), nvec);
  auto d_qp_lqr   = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qp_lqr);

  auto h_qp_lqs   = HostView1D<T>(prs_gsp.data(), nvec);
  auto d_qp_lqs   = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qp_lqs);

  auto h_qp_lqg   = HostView1D<T>(prg_gsp.data(), nvec);
  auto d_qp_lqg   = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qp_lqg);

  auto h_qp_flx   = HostView1D<T>(pre_gsp.data(), nvec);
  auto d_qp_flx   = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qp_flx);

  auto h_qnc      = HostView1D<T>(qnc.data(), nvec);
  auto d_qnc      = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_qnc);

  start_graupel   = std::chrono::steady_clock::now();

  // run the computations
  Kokkos::DefaultExecutionSpace execSpace;
  graupel::run<Kokkos::DefaultExecutionSpace, T>(execSpace, nvec, kend, ivbeg, ivend, kbeg, dt, cia, d_dz, d_t, d_rho,
                                                 d_p, d_qx_lqv, d_qx_lqc, d_qx_lqi, d_qx_lqr, d_qx_lqs, d_qx_lqg, d_qnc,
                                                 d_qp_lqr, d_qp_lqi, d_qp_lqs, d_qp_lqg, d_qp_flx, d_pflx, true);

  end_graupel = std::chrono::steady_clock::now();

  Kokkos::deep_copy(h_qx_lqv, d_qx_lqv);
  Kokkos::deep_copy(h_qx_lqc, d_qx_lqc);
  Kokkos::deep_copy(h_qx_lqr, d_qx_lqr);
  Kokkos::deep_copy(h_qx_lqs, d_qx_lqs);
  Kokkos::deep_copy(h_qx_lqi, d_qx_lqi);
  Kokkos::deep_copy(h_qx_lqg, d_qx_lqg);

  Kokkos::deep_copy(h_qp_lqr, d_qp_lqr);
  Kokkos::deep_copy(h_qp_lqi, d_qp_lqi);
  Kokkos::deep_copy(h_qp_lqs, d_qp_lqs);
  Kokkos::deep_copy(h_qp_lqg, d_qp_lqg);

  Kokkos::deep_copy(h_t, d_t);
  Kokkos::deep_copy(h_pflx, d_pflx);
  Kokkos::deep_copy(h_qp_flx, d_qp_flx);
  Kokkos::fence();

  const std::chrono::duration<double> graupel_time{end_graupel - start_graupel};

  std::cout << "Total graupel [ms]: " << graupel_time.count() * 1000 << std::endl;

  if (std::getenv("OUTPUT_CF"))
    io::write_fields<T, NCT>(output_file, input_file, ncells, nlev, t, qv, qc, qi, qr, qs, qg, prr_gsp, pri_gsp,
                             prs_gsp, prg_gsp, pflx, pre_gsp);
  else
    io::write_fields<T, NCT>(output_file, ncells, nlev, t, qv, qc, qi, qr, qs, qg, prr_gsp, pri_gsp, prs_gsp, prg_gsp,
                             pflx, pre_gsp);

  std::cout << "Results saved to " << output_file << std::endl;
}

int main(int argc, char* argv[]) {
  // Arguments from the command line
  std::string file;
  std::string precision;

  parse_args(file, precision, argc, argv);

  Kokkos::initialize();

  if (precision == "single")
    run_standalone<float, netCDF::NcFloat>(file);
  else if (precision == "double")
    run_standalone<double, netCDF::NcDouble>(file);
  else
    return 1;

  Kokkos::finalize();

  return 0;
}
