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
#include <gtest/gtest.h>

#include <chrono>
#include <thread>

#include "rgk_timer.hpp"

namespace {
using namespace timer;

TEST(TimerTests, CheckTimerMeasurement) {
  // test rough precision of timer measurement:

  Timer timer1("timer1");
  timer1.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(42));
  timer1.stop();

  auto value1 = timer1.value();
  EXPECT_TRUE(value1 >= 0.042 && value1 < 1.0);
}

};  // namespace
