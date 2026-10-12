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

# run command passed as argument and return exit status
function run_command {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
        echo "*** Error with $1" >&2
        # mark command as failed
        (( error_count++ ))
    fi
    return $status
}

[[ ! -f "run/runscript_list" ]] ||  rm run/runscript_list
# loop over all test defined in testlist
echo "Make runscripts for test..."
for test in $(cat $JENKINS_TESTLIST/$TESTLIST); do

    echo $test
    # make runscript for test
    run_command ./run/make_target_runscript in_script=checksuite.icon-dev/check.$test in_script=exec.iconrun out_script=check.$test.run EXPNAME=$test
    echo check.$test.run >> run/runscript_list

    # remove set -x for nicer output
    sed -i '/set -x/d' run/check.$test.run
done
echo "...done"


export PATH=/cluster/software/stacks/2024-06/spack/opt/spack/linux-ubuntu22.04-x86_64_v3/gcc-12.2.0/cdo-2.2.2-542ffdodwastvi2nzovpqxvzjwbr4mpp/bin:$PATH
export HOST='euler'
source .venv/bin/activate
# launch tests
run_command ./scripts/jenkins_scripts/scripts/runexp_euler.bash

echo "Summary of all test:"
if [[ -f "LOOP_STATUS_EXP_FILE" ]];then
    cat LOOP_STATUS_EXP_FILE
else
    echo Summary file not found!
fi
if [[ $error_count -ne "0" ]]; then
    exit 1
else
    exit 0
fi
