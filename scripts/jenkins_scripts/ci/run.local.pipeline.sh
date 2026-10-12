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

set -euo pipefail

[[ $(git rev-parse --show-toplevel 2>/dev/null) = $(pwd) ]] || error "$0 not launched from toplevel of repository"

echo "==> Load module eth_proxy"
module load eth_proxy

echo "==> Setting environment"
export ENV_NAME="euler_cpu_gcc"

echo "==> Reading SPACK_VERSION from config"
SPACK_VERSION=$(cat config/ethz/SPACK_TAG_EULER | tr -d ' \t\n')
export SPACK_VERSION
echo "SPACK_VERSION=$SPACK_VERSION"

echo "==> Running setup_spack.sh"
./scripts/jenkins_scripts/scripts/setup_spack.sh

echo "==> Running build_euler.sh"
./scripts/jenkins_scripts/scripts/build_euler.sh

echo "==> Running jenkins_euler.sh"
./scripts/jenkins_scripts/scripts/jenkins_euler.sh

echo "==> Cleaning up workspace"
# Simulate Jenkins' deleteDir(): clean all files except the script itself (optional)
# Uncomment the following if you want to really clean up:
# find . -mindepth 1 ! -name "$(basename "$0")" -exec rm -rf {} +

echo "==> Done"
