#!/usr/local/bin/bash -l

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

LOCAL_RANK=$SLURM_LOCALID
GLOBAL_RANK=$SLURM_PROCID

N_SOCKETS=$(nvidia-smi --list-gpus | wc -l)

NUMA_NODE=$((LOCAL_RANK % N_SOCKETS))

export CUDA_VISIBLE_DEVICES=$NUMA_NODE

ulimit -s unlimited
numactl --cpunodebind=$NUMA_NODE --membind=$NUMA_NODE bash -c "$@"
