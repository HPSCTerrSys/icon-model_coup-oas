# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: CC0-1.0
# ---------------------------------------------------------------

with section("parse"):
    additional_commands = {
        "gtest_discover_tests": {
            "pargs": {
                "nargs": "+",
                "flags": ["NO_PRETTY_TYPES", "NO_PRETTY_VALUES"],
            },
            "kwargs": {
                "EXTRA_ARGS": "+",
                "WORKING_DIRECTORY": 1,
                "TEST_PREFIX": 1,
                "TEST_SUFFIX": 1,
                "PROPERTIES": "+",
                "TEST_LIST": 1,
                "DISCOVERY_TIMEOUT": 1,
                "XML_OUTPUT_DIR": 1,
                "DISCOVERY_MODE": {
                    "pargs": {
                        "nargs": 1,
                        "flags": ["POST_BUILD", "PRE_TEST"],
                    }
                },
            },
        },
    }


with section("format"):
    dangle_parens = True
    max_lines_hwrap = 0
    keyword_case = "upper"
    autosort = True

with section("markup"):
    first_comment_is_literal = True
