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
#ifndef RAGNAROK_AES_MICROPHYSICS_BINDINGS_FORTRAN_H_
#define RAGNAROK_AES_MICROPHYSICS_BINDINGS_FORTRAN_H_

#include "aes_microphysics/graupel.hpp"

#ifdef __SINGLE_PRECISION
typedef float real_t;
#else
typedef double real_t;
#endif

extern "C" {

void run(const int nvec, const int ke, const int ivstart, const int ivend, const int kstart, const real_t dt,
         const real_t cia, real_t* dz, real_t* t, real_t* rho, real_t* p, real_t* qv, real_t* qc, real_t* qi,
         real_t* qr, real_t* qs, real_t* qg, const real_t* qnc, real_t* prr_gsp, real_t* pri_gsp, real_t* prs_gsp,
         real_t* prg_gsp, real_t* pre_gsp, real_t* pflx);
}

#endif  // RAGNAROK_AES_MICROPHYSICS_BINDINGS_FORTRAN_H_
