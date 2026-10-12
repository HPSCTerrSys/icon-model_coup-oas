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
#
# Reading ERA5 grib files via gribscan
#
# https://gribscan.readthedocs.io/
#
# ------------------------------------------

import argparse
import calendar
import concurrent.futures
import datetime as dt
import logging
import os
import time
from contextlib import contextmanager

import isodate
import numpy as np
import xarray as xr


@contextmanager
def timer(name, logger):
    """
    Context manager for timing the execution of a code block.
    """
    start = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        logger.debug("Elapsed time for %s: %.3f seconds", name, elapsed)


# Argument parser for configuration
def parse_args():
    parser = argparse.ArgumentParser(description="ERA5 OMIP Runoff Provider")
    parser.add_argument(
        "--loglevel",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level (default: INFO)",
    )
    parser.add_argument(
        "--dryrun",
        action="store_true",
        help="Enable DRYRUN mode (default: False)",
    )
    parser.add_argument(
        "--dataPath",
        default="/pool/data/ICON/oes/",
        help="Base data path  (default: /pool/data/ICON/oes/)",
    )
    parser.add_argument(
        "--era5Path",
        default=None,
        help="ERA5 forcing path (default: <dataPath>/ERA5_forcing/)",
    )
    parser.add_argument(
        "--componentName",
        default="era5_provider",
        help="Component name (default: era5_provider)",
    )
    parser.add_argument(
        "--era5GridName",
        default="era5_grid",
        help="ERA5 grid name (default: era5_grid)",
    )
    parser.add_argument(
        "--omipGridName",
        default="omip_grid",
        help="OMIP grid name (default: omip_grid)",
    )
    return parser.parse_args()


def get_dataset(era5Path, year):
    """
    Return an xarray.Dataset for the requested year.
    """

    incat = f"reference::{era5Path}/ERA5_INDEX_{year}.json/atm2d.json"
    logger.info(f"Opening dataset: {incat}")
    dataset = xr.open_dataset(incat, engine="zarr", decode_timedelta=False).sel(
        time=str(year)
    )

    return dataset


# Centralized ERA5 variable metadata
ERA5_VARIABLES = [
    # (var_key, yac_field_name, era5_short_name, scale, offset)
    ("u10", "u10", "10u", 1.0, 0.0),
    ("v10", "v10", "10v", 1.0, 0.0),
    ("ustress", "ustress", "ewss", 1.0 / 3600.0, 0.0),
    ("vstress", "vstress", "nsss", 1.0 / 3600.0, 0.0),
    ("ldown", "ldown", "strd", 1.0 / 3600.0, 0.0),
    ("swdown", "swdown", "ssrd", 1.0 / 3600.0, 0.0),
    ("precip", "precip", "tp", 1.0 / 3600.0, 0.0),
    ("slp", "sea_level_pressure", "msl", 1.0, 0.0),
    ("t2m", "t2m", "2t", 1.0, -273.15),
    ("tcc", "tcc", "tcc", 1.0, 0.0),
    ("tdew", "tdew", "2d", 1.0, 0.0),
]


def read_var(var, idx, shape, scale=None, offset=None):
    """
    Read variable data for a given time index and reshape it.
    Apply scaling and offset to the data.
    """
    arr = var[{"time": idx}].values.reshape(shape)
    if scale is None:
        if offset is None:
            return arr
        else:
            return arr + offset
    else:
        if offset is None:
            return arr * scale
        else:
            return arr * scale + offset


class DummyExecutor(concurrent.futures.Executor):
    """
    A dummy executor that runs tasks synchronously.
    Useful when only one thread is desired.
    """

    def submit(self, fn, *args, **kwargs):
        future = concurrent.futures.Future()
        try:
            result = fn(*args, **kwargs)
            future.set_result(result)
        except Exception as exc:
            future.set_exception(exc)
        return future


# Parse command-line arguments
args = parse_args()

# Set DRYRUN from argument
DRYRUN = args.dryrun
if not DRYRUN:
    from yac import *

# Setup logger
logger = logging.getLogger("ERA5g_provider")
handler = logging.StreamHandler()
formatter = logging.Formatter("%(name)s %(asctime)s %(levelname)s %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(getattr(logging, args.loglevel))

dataPath = args.dataPath
era5Path = (
    args.era5Path if args.era5Path is not None else f"{dataPath}/ERA5_forcing/"
)

# Get number of threads from environment variable, default to 4 if not set
num_threads = int(os.environ.get("OMP_NUM_THREADS", 4))

# Initialize global thread pool
GLOBAL_THREAD_POOL = (
    DummyExecutor()
    if num_threads == 1
    else concurrent.futures.ThreadPoolExecutor(max_workers=num_threads)
)

if not DRYRUN:
    """
    Initialize YAC and define the ERA5 provider component.
    """
    yac = YAC()

    def_calendar(Calendar.PROLEPTIC_GREGORIAN)
    comp = yac.def_comp(args.componentName)

    # Get the number of provider ranks and the rank of the local process
    num_provider_rank = comp.size
    provider_rank = comp.rank

    if provider_rank == 0:
        logger.info(
            f"""
            Initialized YAC component '{args.componentName}' with
            {num_provider_rank} ranks using {num_threads} threads each
            """
        )
else:
    num_provider_rank = 1
    provider_rank = 0

# get sea-land mask and ERA5 data grid information
filename = era5Path + "ERA5_lsm.nc"
logger.info(f"Opening dataset: {filename}")
ds_era5_mask = xr.open_dataset(filename)
is_valid = ds_era5_mask["var172"][{"time": 0}].values
lon = np.deg2rad(ds_era5_mask["lon"])
lat = np.deg2rad(ds_era5_mask["lat"])
ds_era5_mask.close()

# get OMIP runoff grid information
filename = dataPath + "OMIP_runoff/runoff.nc"
logger.info(f"Opening dataset: {filename}")
ds_runoff = xr.open_dataset(filename, decode_timedelta=False)
clon_runoff = np.deg2rad(ds_runoff["lon_edge"])
clat_runoff = np.deg2rad(ds_runoff["lat_edge"])
lon_runoff = np.deg2rad(ds_runoff["lonpt"])
lat_runoff = np.deg2rad(ds_runoff["latpt"])
coast_mask = ds_runoff["coast_mask"]

if not DRYRUN:
    """
    Define the grids for ERA5 and OMIP runoff data.

    - The ERA5 grid is a regular latitude-longitude grid.
    - The OMIP runoff grid is similar but the last longitude edge is a duplicate
      and skipped.
    - Global indices are set for corners and cells to ensure correct mapping.
    - Masks are applied to the OMIP grid to identify coastal cells.

    This setup ensures that the spatial structure of both datasets is correctly
    represented for coupling and data exchange.
    """

    # Generate global indices for ERA5 grid corners
    # (corners at the poles share the same coordinate -> same global index)
    era5_gid_corner = np.arange(lat.size * lon.size, dtype=np.int32).reshape(
        lat.size, lon.size
    )
    era5_gid_corner[0, :] = era5_gid_corner[0, 0]
    era5_gid_corner[-1, :] = era5_gid_corner[-1, -1]

    # Define ERA5 grid and points at corners
    era5_grid = Reg2dGrid(args.era5GridName, lon, lat, cyclic=[True, False])
    era5_grid.set_global_index(era5_gid_corner.ravel(), Location.CORNER)
    era5_points = era5_grid.def_points(Location.CORNER, lon, lat)

    # Generate global indices for OMIP runoff grid corners and cells
    # (corners at the poles share the same coordinate -> same global index)
    # (In contrast to ERA5, the OMIP grid longitude edges include a duplicate
    #  last edge)
    omip_gid_corner = np.arange(
        (clon_runoff.size - 1) * clat_runoff.size, dtype=np.int32
    ).reshape(clat_runoff.size, clon_runoff.size - 1)
    omip_gid_corner[0, :] = omip_gid_corner[0, 0]
    omip_gid_corner[-1, :] = omip_gid_corner[-1, -1]
    omip_gid_cell = np.arange(
        lon_runoff.size * lat_runoff.size, dtype=np.int32
    ).reshape(lat_runoff.size, lon_runoff.size)

    # Define OMIP grid, points at cells centers, and set coast mask
    # (In the omip grid, the last corner is a duplication of the first one)
    omip_grid = Reg2dGrid(
        args.omipGridName, clon_runoff[:-1], clat_runoff, cyclic=[True, False]
    )
    omip_grid.set_global_index(omip_gid_corner.ravel(), Location.CORNER)
    omip_grid.set_global_index(omip_gid_cell.ravel(), Location.CELL)
    runoff_points = omip_grid.def_points(Location.CELL, lon_runoff, lat_runoff)
    runoff_points.set_mask(np.ravel(coast_mask))

if not DRYRUN:
    """
    Create YAC fields for all required variables.

    Each entry in `field_defs` maps a variable name (used as a global)
    to a tuple of (field_name, grid_points). The field is created with
    the specified name, component, grid points, and time interval.

    This approach avoids repetitive code and makes it easy to add or
    modify fields in one place.
    """
    era5_fields = {}
    for i, (var_key, yac_field_name, _, _, _) in enumerate(ERA5_VARIABLES):
        # Distribute variables evenly among provider ranks
        if i % num_provider_rank == provider_rank:
            era5_fields[var_key] = Field.create(
                yac_field_name,
                comp,
                era5_points,
                1,
                "PT1H",
                TimeUnit.ISO_FORMAT,
            )

    # OMIP runoff field is handled by provider_rank 0
    if provider_rank == 0:
        runoff_field = Field.create(
            "river_runoff", comp, runoff_points, 1, "PT1H", TimeUnit.ISO_FORMAT
        )

if not DRYRUN:
    """
    Finalize the YAC component and field definitions.

    This call marks the end of the definition phase for the YAC coupling setup.
    After this, no further components or fields should be defined.

    Collective call that includes weight file computation.
    """
    # synchronize definitions across all processes
    # (helps to improve wall-clock time measurements of enddef by synchronizing
    #  all ranks)
    with timer("synchronizing definitions", logger):
        yac.sync_def()

    with timer("weight computation", logger):
        yac.enddef()

if not DRYRUN:
    start_date = isodate.parse_datetime(yac.start_datetime)
    end_date = isodate.parse_datetime(yac.end_datetime)
else:
    start_date = isodate.parse_datetime("1970-12-30T00:00:00.000")
    end_date = isodate.parse_datetime("1971-01-02T00:00:00.000")

# Initialize ERA5 dataset
ds_era5 = None

# Initialize variables to None
era5_vars = {var_key: None for var_key, *_ in ERA5_VARIABLES}

# OMIP runoff has only one year
var_runoff = ds_runoff["runoff"]

logger.debug("going into time loop from %s to %s", start_date, end_date)

# Get the upper limit for the year loop (exclusive)
end_year = end_date.year
if end_date.month > 1 or end_date.day > 1 or end_date.hour > 0:
    end_year = end_date.year + 1

# Loop over years from start_date to end_date
for year in range(start_date.year, end_year):

    logger.info(f"Processing year {year}")

    # Determine current year's start and end date
    curr_start_date = (
        start_date if year == start_date.year else dt.datetime(year, 1, 1)
    )
    curr_end_date = (
        end_date if year == end_date.year else dt.datetime(year + 1, 1, 1)
    )

    # Determine start and end hour indices for the current year
    start_hour_idx = int(
        (curr_start_date - dt.datetime(year, 1, 1)) / dt.timedelta(hours=1)
    )
    end_hour_idx = start_hour_idx + int(
        (curr_end_date - curr_start_date) / dt.timedelta(hours=1)
    )

    # Load ERA5 dataset for the current year
    with timer("getting ERA5 dataset", logger):
        ds_era5 = get_dataset(era5Path, year)

    # Get variables in a loop
    with timer("getting variables", logger):
        for var_key, _, era5_short_name, _, _ in ERA5_VARIABLES:
            era5_vars[var_key] = ds_era5[era5_short_name]

    # check that all ERA5 variables have the same start time and ensure that
    # it matches curr_year_data
    with timer("getting file start times", logger):
        file_start_times = [
            isodate.parse_datetime(str(era5_vars[var_key].time[0].values))
            for var_key, *_ in ERA5_VARIABLES
        ]
        if len(set(file_start_times)) != 1:
            raise RuntimeError(
                f"Time axis start times mismatch: {file_start_times}"
            )
        if file_start_times[0] != dt.datetime(year, 1, 1):
            raise RuntimeError(
                f"""
                File start time {file_start_times[0]} does not match
                expected {dt.datetime(year, 1, 1)}
                """
            )

    # Loop over hours in the current year
    for hour_idx in range(start_hour_idx, end_hour_idx):

        # Read ERA5 variables in parallel
        with timer("reading ERA5 field data", logger):
            # Prepare tasks for each ERA5 variable handled by this provider rank
            tasks = [
                (
                    var_key,
                    era5_vars[var_key],
                    hour_idx,
                    (lat.size, lon.size),
                    scale,
                    offset,
                )
                for i, (var_key, _, _, scale, offset) in enumerate(
                    ERA5_VARIABLES
                )
                if i % num_provider_rank == provider_rank
            ]
            era5_data = {}
            # Submit tasks to the global thread pool
            future_to_name = {
                GLOBAL_THREAD_POOL.submit(
                    read_var, var, idx, shape, scale, offset
                ): var_key
                for var_key, var, idx, shape, scale, offset in tasks
            }
            # Collect results as they complete
            for future in concurrent.futures.as_completed(future_to_name):
                var_key = future_to_name[future]
                era5_data[var_key] = future.result()

        # Read OMIP runoff data (only provider_rank 0 handles this)
        if provider_rank == 0:
            with timer("reading OMIP runoff field data", logger):

                # OMIP runoff is daily, ERA5 is hourly; idx_runoff is the day-of-year index.
                day_idx_runoff = hour_idx // 24  # 0-based index

                # Handle leap years: if day_idx_runoff is 59 (Feb 28), reuse for 60 (Feb 29)
                if calendar.isleap(year) and day_idx_runoff > 59:
                    day_idx_runoff = day_idx_runoff - 1

                runoff = var_runoff[{"Time": day_idx_runoff}].values.reshape(
                    lat_runoff.size, lon_runoff.size
                )

        logger.debug(
            "sending data for year=%s hour_idx=%s time=%s",
            year,
            hour_idx,
            isodate.parse_datetime(str(era5_vars["u10"].time[hour_idx].values)),
        )

        # Send data via YAC
        if not DRYRUN:
            with timer("field data put", logger):
                for var_key in era5_data:
                    era5_fields[var_key].put(era5_data[var_key])
                if provider_rank == 0:
                    runoff_field.put(runoff)

logger.info("End of time loop")
