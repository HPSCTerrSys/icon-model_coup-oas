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

/// *The public Timer interfaces*:
///  These are the common interfaces provided by
///
///  - the standalone implementation (defined RGK_ENABLE_STANDALONE), and
///  - the icon-integrated implementation (undefined RGK_ENABLE_STANDALONE)
///
///  The main focus lies on the icon-integrated aspect where timing
///  information is relevant for production. Here our implementation
///  delegates all functionality to the existing ICON Fortran
///  implementation (see mo_timer.f90 and mo_real_timer.f90). The
///  standalone timer implementation is a means to keep the instrumented
///  C++ source code valid for, e.g., unit tests. Multithreaded
///  measurements are currently not supported in the standalone version.
///
///  *timer definition*:
///  Timer mytimer("myname");
///
///  Deleted timers will not be reported (see below). Therefore timers
///  should be declared with static storage duration.
///
///  *start of a new measurement interval*:
///  mytimer.start()
///
///  For the icon-integrated case we need to have the icon-backend
///  available. Therefore the ICON subroutine init_ragnarok_bridge() must
///  have been processed.
///
///  *end of the measurement interval*:
///  mytimer.stop()
///
///  Multiple measurement intervals are accumulated. There are assert
///  guards in the source code which help to detect ill-formed start/stop
///  sequences (not active in Release build mode).
///
///  *get current total measurement time*:
///  double current_total = mytimer.value()
///
///  This includes the ongoing measurement without stopping it.
///
///  *timer report*:
///  void report_all_timers();
///
///  The timer-report prints statistics about all existing timers, e.g.,
///  the number of calls and the total measurement time for each timer.
///  Timers which have never been called will be excluded
///  from the timer report. This is only effective for the standalone build.
///
///  In the icon-integrated case the report is triggered by the main
///  Fortran code. Therefore the report_all_timers() function should not
///  be used (it only prints a warning that it is not applicable). The
///  timer report will appear in the end of the logfile of the experiment
///  together with other ICON timers.
#ifndef RAGNAROK_SUPPORT_TIMER_ICON_TIMER_H_
#define RAGNAROK_SUPPORT_TIMER_ICON_TIMER_H_

#include <cassert>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "support/icon_bridge/icon_f2c.hpp"

namespace timer {

// Represents timer functionality used from, e.g., the existing ICON implementation.
// This is a helper class for the Timer class below.
class TimerBackend {
  friend class Timer;

 public:
  using gen_id_t       = int (*)(const char* name, int name_len);
  using timer_switch_t = void (*)(int timer);
  using timer_value_t  = double (*)(int timer);

  // default backend, incomplete support, only valid function is gen_id():
  TimerBackend()
      : gen_id(gen_undef_id), start(invalid_start), stop(invalid_stop), value(invalid_value), is_complete(false) {}

  // full functional backend:
  TimerBackend(gen_id_t gen_id, timer_switch_t start, timer_switch_t stop, timer_value_t value)
      : gen_id(gen_id), start(start), stop(stop), value(value), is_complete(true) {}

 private:
  gen_id_t gen_id;
  timer_switch_t start;
  timer_switch_t stop;
  timer_value_t value;
  bool is_complete = false;
  bool is_valid() const { return is_complete; }

  static constexpr int undef_id = -1;  // represents undefined timer
  // id generation always succeeds; id must be checked by caller
  static int gen_undef_id(const char* name, int name_len) { return undef_id; }
  // these invalid... functions need to abort in debug mode to make a correct code instrumentation easier
  static void invalid_start(int timer) { assert(false); }
  static void invalid_stop(int timer) { assert(false); }
  static double invalid_value(int timer) {
    assert(false);
    return 0.;
  }
};

/*
  class Timer represents a single timer.
  Usage:
    static Timer mytimer("myname");
    ...
    mytimer.start();
    ...
    mytimer.stop();
  Note:
    - The measurement will be included in the ICON timer report.
    - start()/stop() may only be used after a backend has been attached via set_backend().
    - A timer definition is allowed without backend.
*/
class Timer {
  // provide access to set_backend:
  friend void init_timer(const icon::f2c::FunTable& fun_tab);

 public:
  Timer(const std::string& name) : name(name) {
    if (backend.is_valid()) {
      id = backend.gen_id(name.c_str(), name.length());
    } else {
      id = TimerBackend::undef_id;
    }
  }

  void start() {
    if (id == TimerBackend::undef_id) {
      // [[unlikely]] // C++20
      id = backend.gen_id(name.c_str(), name.length());
      assert(id != TimerBackend::undef_id);
    }
    backend.start(id);
  }

  void start(int level) {
    if (level >= timers_level) start();
  }

  void stop() { backend.stop(id); }

  void stop(int level) {
    if (level >= timers_level) stop();
  }

  double value() { return backend.value(id); }

 private:
  int id = TimerBackend::undef_id;
  std::string name;
  static int timers_level;
  static TimerBackend backend;

  // The backend can only be set once and it must be valid.
  static void set_backend(TimerBackend tbe) {
    assert(tbe.is_valid());
    assert(!backend.is_valid());
    backend = tbe;
  }
};

inline void report_all_timers() {
  std::cout << "icon_timer.report_all() not applicable for ragnarok\n"
            << "\n";
}

void init_timer(icon::f2c::FunTable* fun_tab);

}  // namespace timer

#endif  // RAGNAROK_SUPPORT_TIMER_ICON_TIMER_H_
