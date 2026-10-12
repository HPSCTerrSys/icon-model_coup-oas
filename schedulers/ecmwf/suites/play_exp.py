# ------------------------------------------
# Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

import ecflow

suite = "fcmon"
script = "icon_multi.def"
# suite = 'fconly'
# script= 'icon_fc.def'
exp = "dei2_493"

try:
    ci = ecflow.Client()
    ci.sync_local()  # get the defs from the server, and place on ci
    defs = ci.get_defs()  # retrieve the defs from ci
    if len(defs) == 0:
        print("No suites in server, loading defs from disk")
        ci.load(script)

        print("Restarting the server. This starts job scheduling")
        ci.restart_server()
    else:
        print("read definition from disk and replace on the server")
        ci.replace("/" + suite + "/" + exp, script)

#   print("Begin the suite named "+suite)
#   ci.begin_suite(suite)

except RuntimeError as e:
    print("Failed:", e)
