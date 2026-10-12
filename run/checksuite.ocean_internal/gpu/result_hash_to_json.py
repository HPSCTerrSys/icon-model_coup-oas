#! /usr/bin/env python

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


import sys


def get_dirs(args):
    import os

    failreturn = [False, None, None]
    if len(sys.argv) < 2:
        print("No directory given")
        return failreturn

    target_dir = sys.argv[1]
    if not os.path.isdir(target_dir):
        print(f"Directory {target_dir} does not exist/is not accessible")
        return failreturn

    skip_dir = ""
    if len(sys.argv) > 2:
        skip_dir = sys.argv[2]

    return [True, target_dir, skip_dir]


def check_cdo():
    import subprocess

    retval = subprocess.run(["which", "cdo"], capture_output=True)
    if retval.returncode != 0:
        print("cdo check failed")
        return False

    return True


def hash_for_cdo_output(filename):
    import hashlib
    import subprocess

    retval = subprocess.run(
        ["cdo", "-s", "infon", filename], capture_output=True
    )
    if retval.returncode != 0:
        print(f"cdo failed on file {filename}")
        return [False, False]

    result_hash = hashlib.sha256(retval.stdout).hexdigest()
    return [True, result_hash]


def process_dir(target_dir, skip_dir="", level=0, max_levels=1):
    import os

    if level > max_levels:
        return [True, {}]

    curdir = os.path.abspath(os.path.curdir)
    os.chdir(target_dir)
    hash_dict = dict()
    status = True
    for elem in os.listdir("."):
        if elem in [".", ".."]:
            continue

        if os.path.isdir(elem):
            status, subdict = process_dir(
                target_dir=elem,
                skip_dir=skip_dir,
                level=level + 1,
                max_levels=max_levels,
            )
            if subdict != {}:
                hash_dict.update({elem: subdict})

            if not status:
                break

        elif os.path.isfile(elem) and not os.path.islink(elem):
            if elem.endswith(".nc"):
                status, file_hash = hash_for_cdo_output(filename=elem)
                mod_filename = elem
                if elem.startswith(skip_dir):
                    mod_filename = elem[len(skip_dir) :]

                hash_dict.update({mod_filename: file_hash})
                if not status:
                    break

    os.chdir(curdir)
    return [status, hash_dict]


def main(args):
    import json

    status, target_dir, skip_dir = get_dirs(args=args)
    if not status:
        sys.exit(1)

    if not check_cdo():
        sys.exit(1)

    status, hash_dict = process_dir(target_dir=target_dir, skip_dir=skip_dir)
    print(json.dumps(hash_dict, indent=2, sort_keys=True))
    if not status:
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main(args=sys.argv)
