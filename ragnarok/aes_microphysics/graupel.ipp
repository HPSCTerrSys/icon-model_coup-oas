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
template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
using MatrixView = Kokkos::View<T**, Kokkos::LayoutRight, ExecutionSpace>;

template <class ExecutionSpace = Kokkos::DefaultExecutionSpace>
using MinIndexView = Kokkos::View<short**, Kokkos::LayoutRight, ExecutionSpace>;

template <typename T, int size>
KOKKOS_INLINE_FUNCTION void precip(T (&p)[size], const short iqx, const T zeta, const T t, const T flx, const T vt,
                                   const T q, const T q_kp1, const T rho) {
  T rho_x       = q * rho;
  T flx_eff     = rho_x / zeta + static_cast<T>(2.0) * flx;
  T flx_partial = Kokkos::fmin(rho_x * physics::vm(iqx, rho_x, rho, t), flx_eff);

  p[0]          = zeta * (flx_eff - flx_partial) / ((static_cast<T>(1.0) + zeta * vt) * rho);  // q update
  p[1]          = (p[0] * rho * vt + flx_partial) * static_cast<T>(0.5);                       // flx
  rho_x         = (p[0] + q_kp1) * static_cast<T>(0.5) * rho;
  p[2]          = physics::vm(iqx, rho_x, rho, t);
}

template <typename T, short rows, short cols>
KOKKOS_INLINE_FUNCTION T sum_on_line(T (&matrix)[rows][cols], short index) {
  T result = ZERO<T>;
  for (short i = 0; i < cols; ++i) {
    result = result + matrix[index][i];
  }
  return result;
}

template <typename T, short rows, short cols>
KOKKOS_INLINE_FUNCTION T sum_on_column(T (&matrix)[rows][cols], short index) {
  T result = ZERO<T>;
  for (short i = 0; i < rows; ++i) {
    result = result + matrix[i][index];
  }
  return result;
}

template <typename T, short rows, short cols>
KOKKOS_INLINE_FUNCTION void multiply_on_extent(T (&matrix)[rows][cols], short index, T value1, T value2) {
  for (short i = 0; i < rows; ++i) {
    matrix[index][i] = matrix[index][i] * value1 / value2;
  }
}

template <class ExecutionSpace = Kokkos::DefaultExecutionSpace>
KOKKOS_INLINE_FUNCTION short min_on_extent_0(MinIndexView<ExecutionSpace> matrix, int index) {
  short result = matrix(0, index);
  for (short i = 1; i < matrix.extent(0); ++i) {
    short tmp = matrix(i, index);
    result    = (tmp < result) ? tmp : result;
  }
  return result;
}

template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
KOKKOS_INLINE_FUNCTION void compute_sink(T (&d_sink)[idx::nx], T (&d_sx2x)[idx::nx][idx::nx], short iqx, T value, T dt,
                                         bool is_sig) {
  d_sink[iqx] = ZERO<T>;
  if (is_sig || iqx == idx::lqc || iqx == idx::lqv || iqx == idx::lqr) {
    d_sink[iqx] = sum_on_line<T, idx::nx, idx::nx>(d_sx2x, iqx);
    T stot      = value / dt;
    if (d_sink[iqx] > stot && value > graupel::qmin<T>) {
      multiply_on_extent<T, idx::nx, idx::nx>(d_sx2x, iqx, stot, d_sink[iqx]);
      d_sink[iqx] = sum_on_line<T, idx::nx, idx::nx>(d_sx2x, iqx);
    }
  }
}

template <typename T>
KOKKOS_INLINE_FUNCTION T update_qx(T (&d_dqdt)[idx::nx], T (&d_sink)[idx::nx], T (&d_sx2x)[idx::nx][idx::nx], short iqx,
                                   T value, T dt) {
  d_dqdt[iqx] = sum_on_column<T>(d_sx2x, iqx) - d_sink[iqx];
  value       = Kokkos::fmax(ZERO<T>, value + d_dqdt[iqx] * dt);
  return value;
}

template <typename T>
KOKKOS_INLINE_FUNCTION void update_qx(T (&d_dqdt)[idx::nx], T (&d_sink)[idx::nx], T (&d_sx2x)[idx::nx][idx::nx],
                                      T (&d_qx)[idx::nx], T dt, bool is_sig) {
  for (short i = 0; i < idx::nx; ++i) {
    compute_sink(d_sink, d_sx2x, i, d_qx[i], dt, is_sig);
  }

  for (short i = 0; i < idx::nx; ++i) {
    d_qx[i] = update_qx(d_dqdt, d_sink, d_sx2x, i, d_qx[i], dt);
  }
}

template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
KOKKOS_INLINE_FUNCTION void update_qx_qp(MinIndexView<ExecutionSpace> d_kmin, short iqx, int iv, short k,
                                         View2D<T> d_rho, View2D<T> d_t, View2D<T> d_qx, View1D<T> d_qp,
                                         MatrixView<T, ExecutionSpace> d_vt, T zeta, short kp1) {
  if (k >= d_kmin(iqx, iv)) {
    T d_update[3] = {ZERO<T>, ZERO<T>, ZERO<T>};
    precip<T, 3>(d_update, iqx, zeta, d_t(k, iv), d_qp(iv), d_vt(iqx, iv), d_qx(k, iv), d_qx(kp1, iv), d_rho(k, iv));
    d_qx(k, iv)   = d_update[0];
    d_qp(iv)      = d_update[1];
    d_vt(iqx, iv) = d_update[2];
  }
}

template <typename T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
KOKKOS_INLINE_FUNCTION void update_qx_qp(int iv, short k, MinIndexView<ExecutionSpace> d_kmin,
                                         const View2D<T> d_qx[idx::nx], const View1D<T> d_qp[idx::np], View2D<T> d_rho,
                                         View2D<T> d_t, MatrixView<T, ExecutionSpace> d_vt, T zeta, short kp1) {
  for (short i = 0; i < idx::np; ++i) {
    update_qx_qp(d_kmin, i, iv, k, d_rho, d_t, d_qx[i], d_qp[i], d_vt, zeta, kp1);
  }
}

template <typename T>
KOKKOS_INLINE_FUNCTION void compute_sx2x(const T& dt, const T& cia, const T& qnc, T (&d_sx2x)[idx::nx][idx::nx],
                                         const T (&d_qx)[idx::nx], const T& d_rho, const T& d_t, const T& d_p,
                                         const bool is_sig_present) {
  T eta                      = ZERO<T>;
  T ice_dep                  = ZERO<T>;
  T dvsw                     = d_qx[idx::lqv] - thermo::qsat_rho(d_t, d_rho);
  T qvsi                     = thermo::qsat_ice_rho(d_t, d_rho);
  T dvsi                     = d_qx[idx::lqv] - qvsi;
  T rho_x                    = d_rho * d_qx[idx::lqs];
  T n_snow                   = physics::snow_number(d_t, rho_x);
  T l_snow                   = physics::snow_lambda(rho_x, n_snow);

  d_sx2x[idx::lqc][idx::lqr] = transition::cloud_to_rain(d_t, d_rho, d_qx[idx::lqc], d_qx[idx::lqr], qnc);
  d_sx2x[idx::lqr][idx::lqv] = transition::rain_to_vapor(d_t, d_rho, d_qx[idx::lqc], d_qx[idx::lqr], dvsw, dt);
  d_sx2x[idx::lqc][idx::lqi] = transition::cloud_x_ice(d_t, d_qx[idx::lqc], d_qx[idx::lqi], dt);
  d_sx2x[idx::lqi][idx::lqc] = -Kokkos::fmin(d_sx2x[idx::lqc][idx::lqi], ZERO<T>);
  d_sx2x[idx::lqc][idx::lqi] = Kokkos::fmax(d_sx2x[idx::lqc][idx::lqi], ZERO<T>);
  d_sx2x[idx::lqc][idx::lqs] = transition::cloud_to_snow(d_t, d_qx[idx::lqc], d_qx[idx::lqs], n_snow, l_snow);
  d_sx2x[idx::lqc][idx::lqg] = transition::cloud_to_graupel(d_t, d_rho, d_qx[idx::lqc], d_qx[idx::lqg]);

  if (d_t < tmelt<T>) {
    T n_ice = physics::ice_number(d_t, d_rho);
    T m_ice = physics::ice_mass(d_qx[idx::lqi], n_ice);
    T x_ice = physics::ice_sticking(d_t, cia);
    if (is_sig_present) {
      eta                        = physics::deposition_factor(d_t, qvsi);
      d_sx2x[idx::lqv][idx::lqi] = transition::vapor_x_ice(d_qx[idx::lqi], m_ice, eta, dvsi, d_rho, dt);
      d_sx2x[idx::lqi][idx::lqv] = -Kokkos::fmin(d_sx2x[idx::lqv][idx::lqi], ZERO<T>);
      d_sx2x[idx::lqv][idx::lqi] = Kokkos::fmax(d_sx2x[idx::lqv][idx::lqi], ZERO<T>);
      ice_dep                    = Kokkos::fmin(d_sx2x[idx::lqv][idx::lqi], dvsi / dt);

      d_sx2x[idx::lqi][idx::lqs] = physics::deposition_auto_conversion(d_qx[idx::lqi], m_ice, ice_dep);
      d_sx2x[idx::lqi][idx::lqs] =
          d_sx2x[idx::lqi][idx::lqs] + transition::ice_to_snow(d_qx[idx::lqi], n_snow, l_snow, x_ice);
      d_sx2x[idx::lqi][idx::lqg] =
          transition::ice_to_graupel(d_rho, d_qx[idx::lqr], d_qx[idx::lqg], d_qx[idx::lqi], x_ice);
      d_sx2x[idx::lqs][idx::lqg] = transition::snow_to_graupel(d_t, d_rho, d_qx[idx::lqc], d_qx[idx::lqs]);
      d_sx2x[idx::lqr][idx::lqg] = transition::rain_to_graupel(d_t, d_rho, d_qx[idx::lqc], d_qx[idx::lqr],
                                                               d_qx[idx::lqi], d_qx[idx::lqs], m_ice, dvsw, dt);
    }
    d_sx2x[idx::lqv][idx::lqi] = d_sx2x[idx::lqv][idx::lqi] + physics::ice_deposition_nucleation(
                                                                  d_t, d_qx[idx::lqc], d_qx[idx::lqi], n_ice, dvsi, dt);
  } else {
    d_sx2x[idx::lqc][idx::lqr] = d_sx2x[idx::lqc][idx::lqr] + d_sx2x[idx::lqc][idx::lqs] + d_sx2x[idx::lqc][idx::lqg];
    d_sx2x[idx::lqc][idx::lqs] = ZERO<T>;
    d_sx2x[idx::lqc][idx::lqg] = ZERO<T>;
    ice_dep                    = ZERO<T>;
    eta                        = ZERO<T>;
  }

  if (is_sig_present) {
    T dvsw0 = d_qx[idx::lqv] - thermo::qsat_rho(tmelt<T>, d_rho);
    d_sx2x[idx::lqv][idx::lqs] =
        transition::vapor_x_snow(d_t, d_p, d_rho, d_qx[idx::lqs], n_snow, l_snow, eta, ice_dep, dvsw, dvsi, dvsw0, dt);
    d_sx2x[idx::lqs][idx::lqv] = -Kokkos::fmin(d_sx2x[idx::lqv][idx::lqs], ZERO<T>);
    d_sx2x[idx::lqv][idx::lqs] = Kokkos::fmax(d_sx2x[idx::lqv][idx::lqs], ZERO<T>);
    d_sx2x[idx::lqv][idx::lqg] = transition::vapor_x_graupel(d_t, d_p, d_rho, d_qx[idx::lqg], dvsw, dvsi, dvsw0, dt);
    d_sx2x[idx::lqg][idx::lqv] = -Kokkos::fmin(d_sx2x[idx::lqv][idx::lqg], ZERO<T>);
    d_sx2x[idx::lqv][idx::lqg] = Kokkos::fmax(d_sx2x[idx::lqv][idx::lqg], ZERO<T>);
    d_sx2x[idx::lqs][idx::lqr] = transition::snow_to_rain(d_t, d_p, d_rho, dvsw0, d_qx[idx::lqs]);
    d_sx2x[idx::lqg][idx::lqr] = transition::graupel_to_rain(d_t, d_p, d_rho, dvsw0, d_qx[idx::lqg]);
  }
}

template <typename T>
KOKKOS_INLINE_FUNCTION void update_temp(const T& dt, T (&d_dqdt)[idx::nx], const T (&d_qx)[idx::nx], T& d_t) {
  T qice = d_qx[idx::lqs] + d_qx[idx::lqi] + d_qx[idx::lqg];
  T qliq = d_qx[idx::lqc] + d_qx[idx::lqr];
  T qtot = d_qx[idx::lqv] + qice + qliq;
  T cv   = cvd<T> + (cvv<T> - cvd<T>)*qtot + (clw<T> - cvv<T>)*qliq + (ci<T> - cvv<T>)*qice;
  d_t    = d_t + dt *
                  ((d_dqdt[idx::lqc] + d_dqdt[idx::lqr]) * (lvc<T> - (clw<T> - cvv<T>)*d_t) +
                   (d_dqdt[idx::lqi] + d_dqdt[idx::lqs] + d_dqdt[idx::lqg]) * (lsc<T> - (ci<T> - cvv<T>)*d_t)) /
                  cv;
}

template <class ExecutionSpace, typename T>
void graupel::run(ExecutionSpace execSpace, const int nvec, const short ke, const int ivstart, const int ivend,
                  const short kstart, const T dt, const T cia, View2D<T> d_dz, View2D<T> d_t, View2D<T> d_rho,
                  View2D<T> d_p, View2D<T> d_qx_lqv, View2D<T> d_qx_lqc, View2D<T> d_qx_lqi, View2D<T> d_qx_lqr,
                  View2D<T> d_qx_lqs, View2D<T> d_qx_lqg, ConstView1D<T> d_qnc, View1D<T> d_qp_lqr, View1D<T> d_qp_lqi,
                  View1D<T> d_qp_lqs, View1D<T> d_qp_lqg, View1D<T> d_qp_flx, View2D<T> d_pflx, const bool lrain) {
  // create working data structures
  View2D<T> d_qx[idx::nx] = {d_qx_lqr, d_qx_lqi, d_qx_lqs, d_qx_lqg, d_qx_lqc, d_qx_lqv};
  View1D<T> d_qp[idx::np] = {d_qp_lqr, d_qp_lqi, d_qp_lqs, d_qp_lqg};

  // kokkos managed views (host & device)
  Kokkos::View<int, ExecutionSpace> d_jmx_("jmx_");
  Kokkos::View<int, Kokkos::Serial>::HostMirror h_jmx_ = Kokkos::create_mirror_view(d_jmx_);

  // global kokkos::views in DefaultExecutionSpace ONLY, no mirror
  // no need for initialization since the Kokkos::View<T> constructor calls
  // the T default constructor by default
  MatrixView<T, ExecutionSpace> d_vt("vt", idx::np, nvec);
  MinIndexView<ExecutionSpace> d_kmin("kmin", idx::np, nvec);

  ManagedView1D<int, ExecutionSpace> index_iv("index_iv", nvec * ke);
  ManagedView1D<short, ExecutionSpace> index_k("index_k", nvec * ke);
  ManagedView1D<bool, ExecutionSpace> is_sig("is_sig", nvec * ke);

  Kokkos::parallel_for(
      "init_prognostics", Kokkos::RangePolicy(execSpace, ivstart, ivend), KOKKOS_LAMBDA(const int iv)->void {
        d_qp[idx::lqr](iv) = ZERO<T>;
        d_qp[idx::lqi](iv) = ZERO<T>;
        d_qp[idx::lqs](iv) = ZERO<T>;
        d_qp[idx::lqg](iv) = ZERO<T>;
        d_qp_flx(iv)       = ZERO<T>;
      });

  Kokkos::parallel_for(
      "compute_indexes", Kokkos::RangePolicy(execSpace, ivstart, ivend), KOKKOS_LAMBDA(const int iv)->void {
        for (short k = ke - 1; k >= kstart; --k) {
          T freeze_max = ragnarok::fmax(d_qx[idx::lqs](k, iv), d_qx[idx::lqi](k, iv), d_qx[idx::lqg](k, iv));
          if ((ragnarok::fmax(freeze_max, d_qx[idx::lqc](k, iv), d_qx[idx::lqr](k, iv)) > graupel::qmin<T>) ||
              (d_t(k, iv) < graupel::tfrz_het2<T> &&
               d_qx[idx::lqv](k, iv) > thermo::qsat_ice_rho(d_t(k, iv), d_rho(k, iv)))) {
            int index       = Kokkos::atomic_fetch_add(&d_jmx_(), 1);
            index_k(index)  = k;
            index_iv(index) = iv;
            is_sig(index)   = freeze_max > graupel::qmin<T>;
          }
          for (short ix = 0; ix < idx::np; ix++) {
            if (k == ke - 1) {
              d_kmin(ix, iv) = ke + static_cast<short>(1);
            }
            if (d_qx[ix](k, iv) > graupel::qmin<T>) {
              d_kmin(ix, iv) = k;
            }
          }
        }
      });

  // copy jmx_ to the host (if computed on device)
  Kokkos::deep_copy(execSpace, h_jmx_, d_jmx_);
  execSpace.fence();

  Kokkos::parallel_for(
      "compute_qx", Kokkos::RangePolicy(execSpace, 0, *h_jmx_.data()), KOKKOS_LAMBDA(const int j) {
        const int iv               = index_iv(j);
        const short k              = index_k(j);
        const bool val             = is_sig(j);

        T d_sx2x[idx::nx][idx::nx] = {ZERO<T>};
        T d_sink[idx::nx]          = {ZERO<T>};
        T d_dqdt[idx::nx]          = {ZERO<T>};

        T d_qx_tmp[idx::nx]        = {d_qx_lqr(k, iv), d_qx_lqi(k, iv), d_qx_lqs(k, iv),
                                      d_qx_lqg(k, iv), d_qx_lqc(k, iv), d_qx_lqv(k, iv)};

        compute_sx2x(dt, cia, d_qnc(ivstart), d_sx2x, d_qx_tmp, d_rho(k, iv), d_t(k, iv), d_p(k, iv), val);
        update_qx(d_dqdt, d_sink, d_sx2x, d_qx_tmp, dt, val);
        update_temp(dt, d_dqdt, d_qx_tmp, d_t(k, iv));

        for (short i = 0; i < idx::nx; ++i) {
          d_qx[i](k, iv) = d_qx_tmp[i];
        }
      });  // end parallel loop

  if (lrain) {
    Kokkos::parallel_for(
        "update_qx_qp", Kokkos::RangePolicy(execSpace, ivstart, ivend), KOKKOS_LAMBDA(const int iv) {
          const short min_k = min_on_extent_0(d_kmin, iv);
          T d_eflx          = ZERO<T>;
          for (short k = min_k; k < ke; ++k) {
            short kp1 = Kokkos::fmin(ke - 1, k + 1);
            T qliq    = d_qx[idx::lqc](k, iv) + d_qx[idx::lqr](k, iv);
            T qice    = d_qx[idx::lqs](k, iv) + d_qx[idx::lqi](k, iv) + d_qx[idx::lqg](k, iv);
            T e_int =
                thermo::internal_energy(d_t(k, iv), d_qx[idx::lqv](k, iv), qliq, qice, d_rho(k, iv), d_dz(k, iv)) +
                d_eflx;
            T zeta = dt / (static_cast<T>(2.) * d_dz(k, iv));
            update_qx_qp<T, ExecutionSpace>(iv, k, d_kmin, d_qx, d_qp, d_rho, d_t, d_vt, zeta, kp1);
            d_pflx(k, iv) = d_qp[idx::lqs](iv) + d_qp[idx::lqi](iv) + d_qp[idx::lqg](iv);
            d_eflx = dt * (d_qp[idx::lqr](iv) * (clw<T> * d_t(k, iv) - cvd<T> * d_t(kp1, iv) - lvc<T>)+d_pflx(k, iv) *
                           (ci<T> * d_t(k, iv) - cvd<T> * d_t(kp1, iv) - lsc<T>));
            d_pflx(k, iv) = d_pflx(k, iv) + d_qp[idx::lqr](iv);
            qliq          = d_qx[idx::lqc](k, iv) + d_qx[idx::lqr](k, iv);
            qice          = d_qx[idx::lqs](k, iv) + d_qx[idx::lqi](k, iv) + d_qx[idx::lqg](k, iv);
            e_int         = e_int - d_eflx;

            d_t(k, iv) =
                thermo::T_from_internal_energy(e_int, d_qx[idx::lqv](k, iv), qliq, qice, d_rho(k, iv), d_dz(k, iv));
          }
          d_qp_flx(iv) = d_eflx / dt;
        });  // end loop
  }

  execSpace.fence();
}

template <class ExecutionSpace, typename T>
void graupel::run(ExecutionSpace execSpace, const int nvec, const int ke, const int ivstart, const int ivend,
                  const int kstart, const T dt, const T cia, T* dz, T* t, T* rho, T* p, T* qv, T* qc, T* qi, T* qr,
                  T* qs, T* qg, const T* qnc, T* prr_gsp, T* pri_gsp, T* prs_gsp, T* prg_gsp, T* pre_gsp, T* pflx,
                  const bool lrain) {
  auto d_dz     = View2D<T, ExecutionSpace>(dz, ke, nvec);
  auto d_t      = View2D<T, ExecutionSpace>(t, ke, nvec);
  auto d_rho    = View2D<T, ExecutionSpace>(rho, ke, nvec);
  auto d_p      = View2D<T, ExecutionSpace>(p, ke, nvec);
  auto d_pflx   = View2D<T, ExecutionSpace>(pflx, ke, nvec);
  auto d_qx_lqv = View2D<T, ExecutionSpace>(qv, ke, nvec);
  auto d_qx_lqc = View2D<T, ExecutionSpace>(qc, ke, nvec);
  auto d_qx_lqi = View2D<T, ExecutionSpace>(qi, ke, nvec);
  auto d_qx_lqr = View2D<T, ExecutionSpace>(qr, ke, nvec);
  auto d_qx_lqs = View2D<T, ExecutionSpace>(qs, ke, nvec);
  auto d_qx_lqg = View2D<T, ExecutionSpace>(qg, ke, nvec);

  auto d_qp_lqr = View1D<T, ExecutionSpace>(prr_gsp, nvec);
  auto d_qp_lqi = View1D<T, ExecutionSpace>(pri_gsp, nvec);
  auto d_qp_lqs = View1D<T, ExecutionSpace>(prs_gsp, nvec);
  auto d_qp_lqg = View1D<T, ExecutionSpace>(prg_gsp, nvec);
  auto d_qp_flx = View1D<T, ExecutionSpace>(pre_gsp, nvec);

  auto d_qnc    = ConstView1D<T, ExecutionSpace>(qnc, nvec);

  graupel::run<ExecutionSpace, T>(execSpace, nvec, static_cast<short>(ke), ivstart, ivend, static_cast<short>(kstart),
                                  dt, cia, d_dz, d_t, d_rho, d_p, d_qx_lqv, d_qx_lqc, d_qx_lqi, d_qx_lqr, d_qx_lqs,
                                  d_qx_lqg, d_qnc, d_qp_lqr, d_qp_lqi, d_qp_lqs, d_qp_lqg, d_qp_flx, d_pflx, lrain);
}

template <typename T>
void graupel::run(const int nvec, const int ke, const int ivstart, const int ivend, const int kstart, const T dt,
                  const T cia, T* dz, T* t, T* rho, T* p, T* qv, T* qc, T* qi, T* qr, T* qs, T* qg, const T* qnc,
                  T* prr_gsp, T* pri_gsp, T* prs_gsp, T* prg_gsp, T* pre_gsp, T* pflx) {
  // switch to enable precipitation
  const bool lrain = true;

  // force the Kokkos kernels to run Serial backend even if OpenMP is enabled
  // this is a workaround since the graupel is called from the omp_block_loop of ICON
  // and would induce a high oversubscription of threads which impacts performance;
  // if/when the omp_block_loop can be bypassed, this workaround can be deleted
  if constexpr (ragnarok::isCPU()) {
    Kokkos::Serial space = ragnarok::get_serial_exec_space();
    graupel::run<decltype(space), T>(space, nvec, static_cast<short>(ke), ivstart, ivend, static_cast<short>(kstart),
                                     dt, cia, dz, t, rho, p, qv, qc, qi, qr, qs, qg, qnc, prr_gsp, pri_gsp, prs_gsp,
                                     prg_gsp, pre_gsp, pflx, lrain);
  } else {
    auto space = Kokkos::DefaultExecutionSpace();
    graupel::run<decltype(space), T>(space, nvec, static_cast<short>(ke), ivstart, ivend, static_cast<short>(kstart),
                                     dt, cia, dz, t, rho, p, qv, qc, qi, qr, qs, qg, qnc, prr_gsp, pri_gsp, prs_gsp,
                                     prg_gsp, pre_gsp, pflx, lrain);
  }
}
