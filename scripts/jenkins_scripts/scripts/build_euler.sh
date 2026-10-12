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

. spack-c2sm/setup-env.sh

spack env activate config/ethz/spack/$SPACK_VERSION/${ENV_NAME}

NUM_CORES=12

srun_cmd="srun -n ${NUM_CORES} --mem-per-cpu=1G"

${srun_cmd} spack install -v --show-log-on-error -j ${NUM_CORES} icon
