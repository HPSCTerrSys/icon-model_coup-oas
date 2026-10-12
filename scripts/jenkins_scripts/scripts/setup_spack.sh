#!/bin/bash

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

set -e -x

GIT_REMOTE='https://github.com/C2SM/spack-c2sm.git'

git clone --depth 1 --recurse-submodules --shallow-submodules -b ${SPACK_VERSION} ${GIT_REMOTE}
. spack-c2sm/setup-env.sh
