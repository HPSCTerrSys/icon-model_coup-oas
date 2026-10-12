#!/usr/bin/env python3

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

import argparse
import collections
import fnmatch
import multiprocessing
import os
import re
import sys
import textwrap

# REUSE-IgnoreStart
ICON_LICENSE = """\
ICON

---------------------------------------------------------------
Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
Contact information: icon-model.org

See AUTHORS.TXT for a list of authors
See LICENSES/ for license information
SPDX-License-Identifier: BSD-3-Clause
---------------------------------------------------------------\
"""
# REUSE-IgnoreEnd

# ICON directories:
ICON_DIRECTORIES = [
    ".gitlab",
    "config",
    "data",
    "doc",
    "etc",
    "ragnarok",
    "run",
    "scripts",
    "src",
    "support",
    "test",
    "utils",
]

# ICON ignored patterns:
ICON_IGNORED_PATTERNS = [
    "ragnarok/cmake/FindNetCDF.cmake",
    # Atmospheric and Environmental Research and Regents of the University of Colorado
    # (BSD-3-Clause):
    "src/atm_phy_rte_rrtmgp/mo_cloud_optics.f90",
    "src/atm_phy_rte_rrtmgp/mo_load_cloud_coefficients.f90",
    "src/atm_phy_rte_rrtmgp/mo_load_coefficients.f90",
    "src/atm_phy_rte_rrtmgp/mo_simple_netcdf.f90",
    # Atmospheric & Environmental Research, Inc. (BSD-3-Clause):
    "src/atm_phy_schemes/mo_lrtm_*.f90",
    "src/atm_phy_schemes/mo_srtm*.f90",
    # ECMWF (Apache-2.0):
    "scripts/preprocessing/compute_full_level_pressure_on_ml.py",
    "scripts/preprocessing/compute_half_level_pressure_on_ml.py",
    "src/atm_phy_schemes/cloud_random_numbers.f90",
    "src/atm_phy_schemes/data_gwd.f90",
    "src/atm_phy_schemes/mo_adjust.f90",
    "src/atm_phy_schemes/mo_cuascn.f90",
    "src/atm_phy_schemes/mo_cudescn.f90",
    "src/atm_phy_schemes/mo_cuflxtends.f90",
    "src/atm_phy_schemes/mo_cufunctions.f90",
    "src/atm_phy_schemes/mo_cuinit.f90",
    "src/atm_phy_schemes/mo_cumaster.f90",
    "src/atm_phy_schemes/mo_cuparameters.f90",
    "src/atm_phy_schemes/mo_gwd_wms.f90",
    "src/atm_phy_schemes/mo_o3_*.f90",
    "src/atm_phy_schemes/mo_sso_cosmo.f90",
    "src/atm_phy_schemes/mo_sso_ifs.f90",
    "src/atm_phy_schemes/mo_vdftofdc.f90",
    "src/lnd_phy_schemes/mo_voskin.f90",
    "test/microphysics_1mom_schemes/cmake/Modules/FindNetCDF.cmake",
    # Alan Miller (ACM):
    "src/atm_phy_schemes/random_rewrite.f90",
    # AER, Rebecca Adams-Selin (BSD-3-Clause):
    "src/diagnostics/atmosphere/mo_diag_hailcast.f90",
    # AWI (BSD-3-Clause):
    "src/sea_ice/dynamics_fem/mo_ice_fem_advection.f90",
    "src/sea_ice/dynamics_fem/mo_ice_fem_evp.f90",
    "src/sea_ice/dynamics_fem/mo_ice_fem_init.f90",
    "src/sea_ice/dynamics_fem/mo_ice_fem_mesh.f90",
    # SHTOOLS (BSD-3-Clause):
    "src/shr_horizontal/mo_math_legendre.f90",
    # "Unimportant" technical files:
    ".gitlab/ci/data/system-checksum-table.yml",
    ".gitlab/iconbot/notifier.yml",
    # External projects:
    "utils/mkhelper/*",
    "utils/fpp-wrappers/*",
]

FileType = collections.namedtuple(
    "FileType",
    [
        "name",
        "glob_patterns",
        "line_comment_starts",
        "re_license_prefix",
        "license_format_message",
    ],
)

FILE_TYPES = [
    FileType(
        name="Fortran",
        glob_patterns=["*.F90", "*.f90", "*.inc", "*.incf"],
        line_comment_starts=["!"],
        re_license_prefix=None,
        license_format_message="must start on the first line of the file",
    ),
    FileType(
        name="C/C++/CUDA/HIP",
        glob_patterns=[
            "*.c",
            "*.cpp",
            "*.cu",
            "*.h",
            "*.hip.cc",
            "*.hpp",
            "*.ipp",
        ],
        line_comment_starts=["//"],
        re_license_prefix=None,
        license_format_message="must start on the first line of the file",
    ),
    FileType(
        name="Shell/Python",
        glob_patterns=["*.sh", "*.py"],
        line_comment_starts=["#"],
        re_license_prefix=r"(?:#![^\n]+\n\n?)?",
        license_format_message="should start on the first line of the file but can be "
        "prefixes with a shebang and an empty line",
    ),
    FileType(
        name="CMake/YAML",
        glob_patterns=["*/CMakeLists.txt", "*.cmake", "*.yml"],
        line_comment_starts=["#"],
        re_license_prefix=None,
        license_format_message="must start on the first line of the file",
    ),
    FileType(
        name="mkexp-jinja",
        glob_patterns=["*.tmpl"],
        line_comment_starts=["#%#-"],
        re_license_prefix=r"^(?:(?:#!|#%# -\*-).*?\n)?",
        license_format_message="should start on the first line of the file but can be "
        "prefixed with a mode line",
    ),
    FileType(
        name="mkexp-configs",
        glob_patterns=["*.config"],
        line_comment_starts=["#", "#-"],
        re_license_prefix=None,
        license_format_message="must start on the first line of the file",
    ),
]


def get_file_type(filepath):
    for file_type in FILE_TYPES:
        if any(
            fnmatch.fnmatch(filepath, pattern)
            for pattern in file_type.glob_patterns
        ):
            return file_type
    return None


def get_license_header(file_type, line_comment_start):
    return "\n".join(
        (f"{line_comment_start} {line}" if line else line_comment_start)
        for line in ICON_LICENSE.split("\n")
    )


def get_license_regexp(file_type, line_comment_start):
    return "^{0}{1}".format(
        file_type.re_license_prefix or "",
        re.escape(get_license_header(file_type, line_comment_start)),
    )


def get_forbidden_doxygen_regexp(file_type, line_comment_start):
    return r"(?im)^\s*{0}.*@(?:par\s+revision\s+history|author).*$".format(
        line_comment_start
    )


def list_files(dirs_or_files, ignored_patterns):
    for dir_or_file in dirs_or_files:
        if not os.path.exists(dir_or_file):
            raise Exception(f"ERROR: '{dir_or_file}' does not exist")
        elif os.path.isfile(dir_or_file):
            yield dir_or_file
        elif os.path.isdir(dir_or_file):
            for subdir, _, filenames in os.walk(dir_or_file):
                for filename in filenames:
                    filepath = os.path.join(subdir, filename)
                    if get_file_type(filepath) and (
                        ignored_patterns is None
                        or not any(
                            fnmatch.fnmatch(filepath, pattern)
                            for pattern in ignored_patterns
                        )
                    ):
                        yield filepath
        else:
            print(
                f"WARNING: input argument '{dir_or_file}' "
                f"is neither a directory nor a file",
                file=sys.stderr,
            )


def check_file(filepath):
    file_type = get_file_type(filepath)

    if file_type is None:
        print(
            f"WARNING: file '{filepath}' has unknown type",
            file=sys.stderr,
        )
        return filepath, True, []

    with open(filepath, "rb") as f:
        raw = f.read(-1)
        txt = raw.decode("utf-8", errors="replace")

        license_match = False
        for line_comment_start in file_type.line_comment_starts:
            if bool(
                re.match(get_license_regexp(file_type, line_comment_start), txt)
            ):
                license_match = True

        return (
            filepath,
            file_type.name,
            license_match,
            re.findall(
                get_forbidden_doxygen_regexp(
                    file_type, file_type.line_comment_starts
                ),
                txt,
            ),
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Checks that the input files contain the license header (prefixed "
        "with the respective language-specific comment character)"
    )

    parser.add_argument(
        "files_or_dirs",
        nargs="*",
        metavar="FILE_OR_DIRECTORY",
        help="path to a file or directory to check",
    )

    args = parser.parse_args()

    if not args.files_or_dirs:
        icon_dir = os.path.normpath(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        )
        args.files_or_dirs.extend(
            os.path.join(icon_dir, icon_subdir)
            for icon_subdir in ICON_DIRECTORIES
        )
        args.ignored_patterns = [
            os.path.join(icon_dir, pattern) for pattern in ICON_IGNORED_PATTERNS
        ]
    else:
        args.ignored_patterns = None

    return args


def main():
    args = parse_args()

    wrong_license_files = collections.defaultdict(list)
    wrong_doxygen_files = []
    with multiprocessing.Pool() as pool:
        for (
            filepath,
            file_type_name,
            license_correct,
            invalid_doxygen,
        ) in pool.imap_unordered(
            check_file, list_files(args.files_or_dirs, args.ignored_patterns)
        ):
            if not license_correct:
                wrong_license_files[file_type_name].append(filepath)

            if invalid_doxygen:
                wrong_doxygen_files.append((filepath, invalid_doxygen))

    exit_code = 0

    if wrong_license_files:
        exit_code = 1

        for file_type in FILE_TYPES:
            files = wrong_license_files.get(file_type.name, None)
            if not files:
                continue

            print(
                "{0}\n\t{1}\n{2}\n{3}".format(
                    "\n".join(
                        textwrap.wrap(
                            f"ERROR: the following {file_type.name} files do not "
                            f"contain the expected license header:",
                            width=80,
                        )
                    ),
                    "\n\t".join(sorted(files)),
                    "\n".join(
                        textwrap.wrap(
                            f"the expected license header for {file_type.name} files "
                            f"is one of ({file_type.license_format_message}):",
                            width=80,
                        )
                    ),
                    "\n\n".join(
                        get_license_header(file_type, line_comment_start)
                        for line_comment_start in file_type.line_comment_starts
                    ),
                ),
                file=sys.stderr,
            )

    if wrong_doxygen_files:
        exit_code = 1
        print(
            "ERROR: the following files contain forbidden Doxygen directives:\n\t"
            "{0}".format(
                "\n\t".join(
                    "{0}\n{1}".format(x[0], "\n".join(x[1]))
                    for x in sorted(wrong_doxygen_files)
                )
            ),
            file=sys.stderr,
        )

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
