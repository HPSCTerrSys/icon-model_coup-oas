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

import os
import re
import sys

sys.path.insert(0, "scripts/buildbot_scripts")
sys.path.insert(0, "scripts/experiments")
sys.path.insert(0, "doc/www")

import conf as sphinx_config
from buildbot_config import BuildbotConfig
from yaml_experiment_test_processor import cscs_ci_to_data

# Get base_url directly from Sphinx conf.py
base_url = sphinx_config.myst_substitutions.get("base_url")


def to_html_table(data, base_url):
    """
    Convert a dictionary with machines, builders, experiments, and status
    into an HTML table with vertical headers. Builders with only negative
    tests are removed.
    """
    out = ""
    for midx in data["machines"]:
        experiments = data[midx]["experiments"]
        all_builders = data[midx]["builders"]
        status_all = data[midx]["status"]

        # Keep only builders that have at least one True in any experiment column
        indices_keep = [
            j
            for j in range(len(all_builders))
            if any(row[j] for row in status_all)
        ]
        builders = [all_builders[j] for j in indices_keep]
        status = [[row[j] for j in indices_keep] for row in status_all]

        # Skip if no experiments or builders
        if not experiments or not builders:
            continue

        out += '<div class="pst-scrollable-table-container style="text-align: left;">\n'
        out += '<table class="table" style="table-layout: fixed; width: auto; margin-left: 0; margin-right: auto;">\n'

        # Header
        out += "  <thead>\n"
        out += '    <tr class="row-odd">\n'
        out += f'      <th class="head" style="width: 200px;"><p><h1>Machine: {midx}</h1></p></th>\n'
        for bidx in builders:
            out += (
                '      <th class="head" style="width: 5px;">'
                f'<div style="writing-mode: vertical-rl; transform: rotate(180deg); white-space: nowrap;">{bidx}</div>'
                "</th>\n"
            )
        out += "    </tr>\n"
        out += "  </thead>\n"

        # Body
        out += "  <tbody>\n"
        row_class = "row-even"
        for exp, row in zip(experiments, status):
            exp_name = os.path.basename(exp)
            out += f'    <tr class="{row_class}">\n'
            out += f'      <td><p><a class="reference external" href="{base_url}/run/{exp}"><code class="docutils literal notranslate"><span class="pre">{exp_name}</span></code></a></p></td>\n'
            for val in row:
                out += (
                    "      <td><p>✔︎</p></td>\n"
                    if val
                    else "      <td><p> </p></td>\n"
                )
            out += "    </tr>\n"
            row_class = "row-odd" if row_class == "row-even" else "row-even"

        out += "  </tbody>\n</table>\n</div>\n\n"

    return out


# Load experiments list
exp_list = BuildbotConfig.from_pickle(
    "scripts/buildbot_scripts/experiment_lists/merge2rc"
)

data_bb = exp_list.buildbot_to_data()
data_ci = cscs_ci_to_data("merge2rc")

# Merge data_bb and data_ci
data = data_bb
for midx in data_ci["machines"]:
    data["machines"].append(midx)
    data["machines"].sort()
    data[midx] = {
        "builders": data_ci[midx]["builders"],
        "experiments": data_ci[midx]["experiments"],
        "status": data_ci[midx]["status"],
    }

# Construct HTML table of test list
ci_test_list_table = to_html_table(data, base_url)

# Read documentation file
doc_path = "doc/www/documentation/infrastructure/testing/system_tests.md"
if not os.path.exists(doc_path):
    print(f"Documentation file not found: {doc_path}")
    sys.exit(1)

with open(doc_path, "r", encoding="utf-8") as f:
    content = f.read()

# Replace the placeholder with the HTML table
new_content = re.sub(
    r"<!-- EXTERNAL-CI-SYSTEM-TESTS -->",
    "<details>\n"
    "<summary>CI System Tests on External Machines by Builder (Click to Expand)</summary>\n\n"
    + ci_test_list_table
    + "</details>\n",
    content,
    count=1,
)

with open(doc_path, "w", encoding="utf-8") as f:
    f.write(new_content)
