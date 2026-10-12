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

# ~~~
# This script fetches files from a private Git repository. The files are a
# subset of sources from the DACE (Data Assimilation Coding Environment) code,
# required for the data assimilation online forward operators in ICON. The
# repository is fetched:
#   - at the configure time of ICON when the DACE modules for data assimilation
#     are enabled (--enable-dace);
#   - via the SSH protocol by default, or via the URL provided in the
#     ICON_DACE_GIT_REPOSITORY environment variable.
# The specific reference to fetch is controlled via the ICON_DACE_GIT_TAG
# variable defined below.
# ~~~

# Reference (tag/branch/commit) to fetch:
set(ICON_DACE_GIT_TAG "a56af966bfc607d64b96580d1f35d70996800547")

# All required functionality works as expected to work with CMake 3.14, except
# for the status messages, which requires CMake 3.17 to make them look as
# expected:
cmake_minimum_required(VERSION 3.14)

project(IconDaceFetchOnly NONE)

find_package(Git REQUIRED)

# Initialize from environment and fallback to default if empty:
set(git_repository "$ENV{ICON_DACE_GIT_REPOSITORY}")
if(NOT git_repository)
  set(git_repository "git@gitlab.dkrz.de:dwd-sw/dace-icon-interface.git")
endif()

# Make it possible to override the values on the command line:
set(ICON_DACE_SOURCE_DIR
    "${PROJECT_BINARY_DIR}/src"
    CACHE PATH "Path to the directory to the repository to"
)
set(ICON_DACE_FETCH_STAMP_FILE
    "${PROJECT_BINARY_DIR}/.fetchstamp"
    CACHE FILEPATH
          "Path to the file signaling the update of the repository update"
)

message(
  CHECK_START
  "Fetching '${ICON_DACE_GIT_TAG}' from '${git_repository}' to '${ICON_DACE_SOURCE_DIR}'"
)

# We want to fetch as little as possible while still allowing switching to
# arbitrary commits even after the initial clone. The latter is impossible with
# older Git versions and the way shallow cloning is implemented in CMake.
# Therefore, we do the shallow cloning and minimize the amount of transferred
# data only for newer Git versions:
set(git_config "")
set(git_shallow OFF)
if(GIT_VERSION_STRING VERSION_GREATER_EQUAL "2.19")
  list(
    APPEND
    git_config
    "remote.origin.promisor=true"
    "remote.origin.partialclonefilter=blob:none"
  )
  set(git_shallow ON)
endif()

include(FetchContent)
FetchContent_Declare(
  dace-icon-interface
  SOURCE_DIR "${ICON_DACE_SOURCE_DIR}"
  GIT_REPOSITORY "${git_repository}"
  GIT_TAG "${ICON_DACE_GIT_TAG}"
  GIT_SUBMODULES ""
  GIT_SHALLOW "${git_shallow}"
  GIT_CONFIG "${git_config}"
)
FetchContent_MakeAvailable(dace-icon-interface)

message(CHECK_PASS "done")

# Git-ignore all files in the build directory:
set(gitignore_file "${PROJECT_BINARY_DIR}/.gitignore")
if(NOT EXISTS "${gitignore_file}")
  file(WRITE "${gitignore_file}" "*")
endif()

# Create a dummy makefile include file, which is touched when ICON_DACE_GIT_TAG
# is changed:
file(
  GENERATE
  OUTPUT "${ICON_DACE_FETCH_STAMP_FILE}"
  CONTENT "# ${ICON_DACE_GIT_TAG}\n"
)

# Allow for 'make test' run by the ICON build system for all CMake-based
# subprojects:
enable_testing()
