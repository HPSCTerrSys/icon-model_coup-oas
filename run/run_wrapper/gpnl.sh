#! /bin/bash

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

set +eu

nvsmi_logger_PID=0
function kill_nvsmi()
{
    set +x
    if (( nvsmi_logger_PID != 0 ))
    then
        kill $nvsmi_logger_PID
    fi
}
trap kill_nvsmi ERR
trap kill_nvsmi EXIT

lrank=$OMPI_COMM_WORLD_RANK
compute_tasks=$OMPI_COMM_WORLD_SIZE

echo "compute_tasks: $compute_tasks, lrank: $lrank, CUDA_VISIBLE_DEVICES: '$CUDA_VISIBLE_DEVICES'"

# Local Task 0 runs always the nvidia-smi logger
# To enable logging of nvidia-smi by setting the following in your run script
# (either directly or via create_target_header)
#     export ENABLE_NVIDIA_SMI_LOGGER=yes
if [[ "$lrank" == 0 ]] && [[ ${ENABLE_NVIDIA_SMI_LOGGER:-"no"} == "yes" ]]
then
    set +x
    # Start logger in background. It will be killed by the ERR trap or kill_nvsmi.
    loop_repetition_time=500000000 # in nano seconds
    while sleep 0.$(( ( 1999999999 - 1$(date +%N) ) % loop_repetition_time ))
    do
        LC_TIME=en_US date -Ins
        nvidia-smi --format=csv --query-gpu=index,power.draw,utilization.gpu,temperature.gpu,memory.used
    done > "nvsmi.log.${lrank}" &
    nvsmi_logger_PID=$!
    set -x
fi

# Use CUDA_VISIBLE_DEVICES, if set, to respect the GPUs assigned by NQSV
IFS=',' read -r -a gpus <<< "${CUDA_VISIBLE_DEVICES:-0,1,2,3,4,5,6,7}"
# gpus is an array of the available GPU IDs

if [[ 0 == ${#gpus[@]} ]]
then
    echo "No GPU visible. CUDA_VISIBLE_DEVICES: '$CUDA_VISIBLE_DEVICES'"
    echo "Run Process $OMPI_COMM_WORLD_RANK on $(hostname)"
else
    # split the blocks for later use of gpnl nodes as I/O server
    if [[ $lrank < $compute_tasks ]]
    then
        export CUDA_VISIBLE_DEVICES=${gpus[lrank % ${#gpus[@]}]}
        echo "Compute process ${OMPI_COMM_WORLD_RANK} on $(hostname) with GPU ${CUDA_VISIBLE_DEVICES}"
    else
        echo "IO process $OMPI_COMM_WORLD_RANK on $(hostname)"
    fi
fi

export KMP_AFFINITY=scatter

$@
return=$?

echo "nvsmi_logger_PID $nvsmi_logger_PID"
kill_nvsmi
exit $return
