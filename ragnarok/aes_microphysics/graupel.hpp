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
#ifndef RAGNAROK_AES_MICROPHYSICS_GRAUPEL_H_
#define RAGNAROK_AES_MICROPHYSICS_GRAUPEL_H_

#include <Kokkos_Core.hpp>

#include "aes_microphysics/physics.hpp"
#include "aes_microphysics/transition.hpp"
#include "aes_thermodynamics/thermo.hpp"
#include "aes_thermodynamics/thermo_constants.hpp"
#include "common/types.hpp"
#include "common/utilities.hpp"
#include "support/icon_bridge/ragnarok.hpp"

namespace graupel {

/// Below are coefficients to parameterize the radar reflectivity and the reflectivity
/// weighted fall speed of the hydrometeors.  They are not presently used in the Muphys
/// code, but are documented here for possible diagnostic use, or to allow the extension
/// of the code to faciliate the output of these variables in the future.
/// author Bjorn Stevens (MPIM)
// template <typename T>
// constexpr T rain_re[5] = {T{2.17470e+07}, T{8.47148e+06}, T{1.19189e+06}, T{7.18530e+04}, T{1.56936e+03}};
// template <typename T>
// constexpr T rain_vd[5] = {T{-6.81479e+00}, T{-7.53619e+00}, T{-1.20634e+00}, T{-7.05880e-02}, T{-1.42711e-03}};

///
/// @brief
/// @param [in] nvec Number of horizontal points
/// @param [in] ke Number of grid points in vertical direction
/// @param [in] ivstart Start index for horizontal direction
/// @param [in] ivend End index for horizontal direction
/// @param [in] kstart Start index for vertical direction
/// @param [in] dt Time step for integration of microphysics (s)
/// @param [in] cia Parameter to change the ice sticking efficiency
/// @param [in] dz Layer thickness of full levels (m)
/// @param [inout] t Temperature in Kelvin
/// @param [in] rho Density of moist air (kg/m3)
/// @param [in] p Pressure (Pa)
/// @param [inout] qv Specific water vapor content (kg/kg)
/// @param [inout] qc Specific cloud water content (kg/kg)
/// @param [inout] qi Specific cloud ice content (kg/kg)
/// @param [inout] qr Specific rain content (kg/kg)
/// @param [inout] qs Specific snow content  kg/kg)
/// @param [inout] qg Specific graupel content (kg/kg)
/// @param [in] qnc Cloud number concentration
/// @param [out] prr_gsp Precipitation rate of rain, grid-scale (kg/(m2*s))
/// @param [out] pri_gsp Precipitation rate of ice, grid-scale (kg/(m2*s))
/// @param [out] prs_gsp Precipitation rate of snow, grid-scale (kg/(m2*s))
/// @param [out] prg_gsp Precipitation rate of graupel, grid-scale (kg/(m2*s))
/// @param [out] pre_gsp Energy flux at sfc from precipitation (W/m2)
/// @param [out] pflx Total precipitation flux
///
template <typename T>
void run(const int nvec, const int ke, const int ivstart, const int ivend, const int kstart, const T dt, const T cia,
         T* dz, T* t, T* rho, T* p, T* qv, T* qc, T* qi, T* qr, T* qs, T* qg, const T* qnc, T* prr_gsp, T* pri_gsp,
         T* prs_gsp, T* prg_gsp, T* pre_gsp, T* pflx);

template <class ExecutionSpace, typename T>
void run(ExecutionSpace execSpace, const int nvec, const int ke, const int ivstart, const int ivend, const int kstart,
         const T dt, const T cia, T* dz, T* t, T* rho, T* p, T* qv, T* qc, T* qi, T* qr, T* qs, T* qg, const T* qnc,
         T* prr_gsp, T* pri_gsp, T* prs_gsp, T* prg_gsp, T* pre_gsp, T* pflx, const bool lrain);

template <class ExecutionSpace, typename T>
void run(ExecutionSpace execSpace, const int nvec, const short ke, const int ivstart, const int ivend,
         const short kstart, const T dt, const T cia, View2D<T> d_dz, View2D<T> d_t, View2D<T> d_rho, View2D<T> d_p,
         View2D<T> d_qx_lqv, View2D<T> d_qx_lqc, View2D<T> d_qx_lqi, View2D<T> d_qx_lqr, View2D<T> d_qx_lqs,
         View2D<T> d_qx_lqg, ConstView1D<T> qnc, View1D<T> d_qp_lqr, View1D<T> d_qp_lqi, View1D<T> d_qp_lqs,
         View1D<T> d_qp_lqg, View1D<T> d_qp_flx, View2D<T> d_plfx, const bool lrain);
}  // namespace graupel

// include implementation
#include "graupel.ipp"

#endif  // RAGNAROK_AES_MICROPHYSICS_GRAUPEL_H_
