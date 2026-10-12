# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

run_cmd=
test -z "$MPI_LAUNCH" || run_cmd="$MPI_LAUNCH -n 1"

test_prog_status=0
$run_cmd $TEST_PROG || test_prog_status=$?

# The test program is expected to fail:
test 0 -ne "$test_prog_status"
