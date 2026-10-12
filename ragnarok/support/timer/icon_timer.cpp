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
///
/// @file
/// @brief This file provides access to timer functionality from ICON
///
//----------------------------

#include "icon_timer.hpp"

namespace timer {

// define static Timer data:
int Timer::timers_level = 0;
TimerBackend Timer::backend{};

namespace {
void do_nothing(int timer) {};
}

void init_timer(const icon::f2c::FunTable& ftab) {
  static bool is_initialized = false;
  if (is_initialized) return;

  // set timer backend:
  bool ltimer;  //  == (time measurements are active), set below
  ftab.get_timer_config(ltimer, Timer::timers_level);
  if (ltimer) {
    Timer::set_backend(TimerBackend(ftab.new_timer, ftab.timer_start, ftab.timer_stop, ftab.timer_value));
  } else {
    Timer::set_backend(TimerBackend(ftab.new_timer, do_nothing, do_nothing, ftab.timer_value));
  }
  is_initialized = true;
}

}  // namespace timer

extern "C" void init_ragnarok_timer() {
  auto& ftab = icon::f2c::get_fun_table();
  timer::init_timer(ftab);
}
