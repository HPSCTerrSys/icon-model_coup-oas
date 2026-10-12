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

####### utility functions for Jenkins ######

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

function error {
    echo "*** Error: $@" >&2
    exit 1
}

function warning {
    echo "*** WARNING: $@" >&2
}

display_help() {
    echo
    echo
    echo "DISCLAIMER:"
    echo " JENKINS 4 ICON is more a wrapper for the ICON buildbot testing framework than a stand-alone testing tool!"
    echo " Changes in the runscripts of buildbot can affect JENKINS 4 ICON in unpredicted way!"
    echo
    echo
    echo
    echo "Usage: $0 [option...] " >&2
    echo
    echo "   -t, --testlist <name of list>                use custom testlist default: test.<host>.cpu"
    echo "  -nc, --no-of-cores <integer>                  number of cores"
    echo "  -nt, --node-type <integer>                    type of Euler-nodes: [6:7]"
    echo "   -h, --help                                   print what you are currently looking at"
    echo
    exit 1
}

###### start main script ######

# simple argparser
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -t|--testlist)      TESTLIST=$2; shift ;;
        -nc|--no-of-cores)  JENKINS_NO_OF_CORES=$2; shift ;;
        -nt|--node-type)    JENKINS_NODE_TYPE=$2; shift ;;
        -h|--help)          display_help ;;
        *) echo "Unknown parameter passed: $1! Use -h or --help for more information"; exit 1 ;;
    esac
    shift
done

# get name of HPC-Centre
host=$(hostname)
if [[ $host == *eu* ]]; then
    centre="ethz"
    host="euler"
else
    error "Unknown HPC-centre with hostname $host!"
fi

# define defaults for testing
[[ -z $TESTLIST ]] && TESTLIST="test.$host.cpu"
[[ -z $JENKINS_NO_OF_CORES ]] && JENKINS_NO_OF_CORES=12
[[ -z $JENKINS_NODE_TYPE ]] && JENKINS_NODE_TYPE=7

# export important directories as environment variables
export JENKINS_ROOTDIR=$(pwd)
export JENKINS_SCRIPTS=${JENKINS_ROOTDIR%%/}/scripts/jenkins_scripts/scripts
export JENKINS_CONFIGURE=${JENKINS_ROOTDIR%%/}/config/${centre}
export JENKINS_TESTLIST=${JENKINS_SCRIPTS}
export TESTLIST
export COMPILER
export BB_NAME=$host.gcc.cpu

# needed for mpi on Euler
export JENKINS_NO_OF_CORES
export JENKINS_NODE_TYPE

# flag for error tracking
error_count=0

echo ========== Start Jenkins ==========

# currently Jenkins only supports normal build -> out-of-source is aborted
[[ $(git rev-parse --show-toplevel 2>/dev/null) = $(pwd) ]] || error "$0 not launched from toplevel of repository"

echo ========== Check user input =======

# Euler nodes type
if [[ $host == "euler" ]]; then
    if [[ "$JENKINS_NODE_TYPE" -lt 8 && "$JENKINS_NODE_TYPE" -gt 5 ]]; then
        echo "Run on Euler $JENKINS_NODE_TYPE"
    else
        error "Invalid Euler node: $JENKINS_NODE_TYPE"
    fi
fi


# existence of testlist and configure
[[ -f ${JENKINS_TESTLIST}/${TESTLIST} ]] || error "Cannot find testlist $TESTLIST"

echo "========== *** VALID *** =========="


echo ========== Start tests ============

echo Testlist $TESTLIST contains these tests:
cat ${JENKINS_TESTLIST}/${TESTLIST}

# activate python venv
run_command ${JENKINS_SCRIPTS}/setup_python.sh

# launch tests
run_command ${JENKINS_SCRIPTS}/test_euler.sh

if [[ $error_count -ne "0" ]]; then
    error "====== *** TEST FAIL *** ========"
else
    echo "====== *** TEST SUCCESS *** ======="
fi

echo ========== End Jenkins ============
