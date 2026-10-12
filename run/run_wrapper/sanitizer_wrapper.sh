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

# Run a MPI task with the compute-sanitizer attached to certain tasks
# Execution:  srun -srun_options ./sanitizer-wrapper.sh path/to/icon
#
set -o pipefail

if [[ -n "${SLURM_PROCID}" ]]; then
    # Use rank information from SLURM
    MPI_RANK=${SLURM_PROCID}
elif [[ -n "${OMPI_COMM_WORLD_RANK}" ]]; then
    # Use rank information from open MPI
    MPI_RANK=${OMPI_COMM_WORLD_RANK}
else
    echo "sanitizer_wrapper.sh: Can not detect current tasks MPI rank." >&2
    exit 1
fi

# chose the MPI rank here
# if [[ "${MPI_RANK}" == 0 ]]; then
# Frequently used options of compute-sanitizer
# ============================================
#
# --tool arg    # Set the tool to use.
#                 memcheck  : Memory access checking
#                 initcheck : Global memory initialization checking
#                 racecheck : Shared memory hazard checking
#                 synccheck : Synchronization checking
#
# --print-limit 10  # Limit the printing to 10 issues.
#
# Filters to control which kernels are checked:
#    --kernel-name kns=mo_nwp_sfc_interface               <-- kns = kernel name substring, kne = full kernel name
#    --kernel-name kns=terra                              <-- check only kernels inside TERRA
#    --kernel-name-exclude kns=mo_communication_orig      <-- to exclude kernels
# Notes:
#   - Generally, the function names in which the kernels are located, are substrings in the full kernel name.
#   - More than one filter can be specified (e.g., --kernel-name kns=terra,kns=sfc)
#
#    --check-api-memory-access no \
#    --check-device-heap no \
#
# For more options check the printout of `compute-sanitizer --help`

if [[ "${MPI_RANK}" < 4 ]]; then # change me to run the sanitizer on more than one node
    compute-sanitizer \
    --tool initcheck \
    --print-limit 10 \
    --error-exitcode 1 \
    --show-backtrace no \
    --launch-timeout 180 \
    "$@" 2>&1 | tee "LOG.sanitizer.${MPI_RANK}.${SLURM_JOBID}.${MPI_RANK}"
else
    "$@"
fi
