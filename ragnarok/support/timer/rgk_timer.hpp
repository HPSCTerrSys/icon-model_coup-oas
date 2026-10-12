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

// standalone and reference timer for testing

#ifndef RAGNAROK_SUPPORT_TIMER_RGK_TIMER_HPP_
#define RAGNAROK_SUPPORT_TIMER_RGK_TIMER_HPP_

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "Kokkos_Timer.hpp"

namespace timer {

// Ragnarok standalone Timer
// serial version: multithreaded time measurements are not supported

class Timer;

// Helper class for Timer.
class TimerTable {
  friend class Timer;
  friend void report_all_timers();

  static inline std::vector<Timer*> table;

  static int add_timer(Timer* timer) {
    // enter new timer into timer table:
    table.push_back(timer);
    int id = table.size() - 1;
    return id;
  }

  static void delete_timer(int id) {
    // delete timer from timer table:
    assert(id >= 0 and id < table.size());
    table[id] = nullptr;
    trim_table();
  }

  // const access to table:
  static const auto& get_all_timers() { return table; }

  static void trim_table() {
    while (!table.empty() && table.back() == nullptr) {
      table.pop_back();
    }
  }
};

class Timer {
  friend void report_all_timers();

 public:
  // start measurement interval:
  void start() {
    if (s_enable_measurements) {
      assert(not is_on);
      is_on = true;
      start_timestamp.reset();
    }
  }

  // start conditional measurement:
  void start(int level) {
    if (level >= s_minlevel) start();
  }

  // end measurement interval:
  void stop() {
    if (s_enable_measurements) {
      time_sum += start_timestamp.seconds();
      assert(is_on);
      is_on = false;
      call_num++;
    }
  }

  // end conditional measurement:
  void stop(int level) {
    if (level >= s_minlevel) stop();
  }

  // return current time sum plus time spent in ongoing measurement:
  double value() const {
    double val = time_sum;
    if (is_on) val += start_timestamp.seconds();
    return val;
  }

  Timer(const std::string& name) : name(name) { table_pos = TimerTable::add_timer(this); }

  ~Timer() { TimerTable::delete_timer(table_pos); }

 private:
  Kokkos::Timer start_timestamp;  // represents start time of measurement
  double time_sum = 0.0;          // sum over all measurements
  bool is_on      = false;        // true if measurement is ongoing
  int call_num    = 0;            // number of measurements
  std::string name;               // timer description in report
  int table_pos = -1;             // position within timer table

  bool is_used() const { return call_num > 0; }

  std::string get_report() const {
    std::stringstream ss;
    ss << name << ": ncalls=" << call_num << ", total=" << value();
    return ss.str();
  }

  // write timer report to stdout:
  void report() const { std::cout << get_report() << std::endl; }

  static const bool s_enable_measurements = true;  // switch for all measurements
  static const int s_minlevel             = 1;     // extra threshold for conditional measurements
};

inline void report_all_timers() {
  std::cout << "Timer report:\n";
  auto& all_timers = TimerTable::get_all_timers();
  for (const Timer* t : all_timers) {
    if (t and t->is_used()) t->report();
  }
}

}  // namespace timer

#endif  // RAGNAROK_SUPPORT_TIMER_RGK_TIMER_HPP_
