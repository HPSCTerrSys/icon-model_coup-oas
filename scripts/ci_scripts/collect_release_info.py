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
import sys

# Read documentation release file
doc_path = "doc/www/release/release.md"
if not os.path.exists(doc_path):
    print(f"Documentation file not found: {doc_path}")
    sys.exit(1)

release_notes = "RELEASE_NOTES.md"
if not os.path.exists(release_notes):
    print(f"Release notes not found: {release_notes}")
    sys.exit(1)

with open(release_notes, "r", encoding="utf-8") as f:
    release_notes_content = f.read()
    release_notes_content = release_notes_content.replace(
        """# Release notes""", """## Release notes"""
    )

with open(doc_path, "a", encoding="utf-8") as f:
    f.write("""\n""")
    f.write(release_notes_content)
