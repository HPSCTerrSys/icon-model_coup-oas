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
#ifndef RAGNAROK_STANDALONE_AES_MICROPHYSICS_IO_NETCDF_CLIENT_H_
#define RAGNAROK_STANDALONE_AES_MICROPHYSICS_IO_NETCDF_CLIENT_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <netcdf>
#include <numeric>
#include <string>
#include <vector>

namespace io {

constexpr int NC_ERR           = 2;
constexpr const char* BASE_VAR = "zg";

template <typename T>
void input_vector(netCDF::NcFile& datafile, std::vector<T>& v, const std::string& input, int ncells, int nlev,
                  int itime);

template <typename T>
void input_vector(netCDF::NcFile& datafile, std::vector<T>& v, const std::string& input, int ncells, int nlev);

template <typename T, typename NCT>
void output_vector(netCDF::NcFile& datafile, std::vector<netCDF::NcDim>& dims, const std::string& output,
                   const std::vector<T>& v, int ncells, int nlev, int deflate_level);

template <typename T, typename NCT>
void output_vector(netCDF::NcFile& datafile, std::vector<netCDF::NcDim>& dims, const std::string& output,
                   std::map<std::string, netCDF::NcVarAtt>, const std::vector<T>& v, int ncells, int nlev,
                   int deflate_level);

template <typename T>
void read_fields(const std::string& input_file, int itime, int& ncells, int& nlev, std::vector<T>& z, std::vector<T>& t,
                 std::vector<T>& p, std::vector<T>& rho, std::vector<T>& qv, std::vector<T>& qc, std::vector<T>& qi,
                 std::vector<T>& qr, std::vector<T>& qs, std::vector<T>& qg);

template <typename T, typename NCT>
void write_fields(const std::string& output_file, int ncells, int nlev, const std::vector<T>& t,
                  const std::vector<T>& qv, const std::vector<T>& qc, const std::vector<T>& qi,
                  const std::vector<T>& qr, const std::vector<T>& qs, const std::vector<T>& qg,
                  const std::vector<T>& prr_gsp, const std::vector<T>& pri_gsp, const std::vector<T>& prs_gsp,
                  const std::vector<T>& prg_gsp, const std::vector<T>& pflx, const std::vector<T>& pre_gsp);

template <typename T, typename NCT>
void write_fields(const std::string& output_file, const std::string& input_file, int ncells, int nlev,
                  const std::vector<T>& t, const std::vector<T>& qv, const std::vector<T>& qc, const std::vector<T>& qi,
                  const std::vector<T>& qr, const std::vector<T>& qs, const std::vector<T>& qg,
                  const std::vector<T>& prr_gsp, const std::vector<T>& pri_gsp, const std::vector<T>& prs_gsp,
                  const std::vector<T>& prg_gsp, const std::vector<T>& pflx, const std::vector<T>& pre_gsp);

}  // namespace io

// include implementation
#include "netcdf_client.ipp"

#endif  // RAGNAROK_STANDALONE_AES_MICROPHYSICS_IO_NETCDF_CLIENT_H_
