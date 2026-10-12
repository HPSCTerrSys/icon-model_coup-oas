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
# some early definition changable by user
#
VERBOSE = 3
DRYRUN = False
UPDATE_CASE = 1
#
# ---------------------------------------------------------------

import socket
import sys
from datetime import date, datetime
from glob import glob

import f90nml
import isodate
import numpy as np
import pandas as pd
import xarray as xr

if not DRYRUN:
    from yac import *


def get_hostname():

    fqdn = socket.getfqdn().split(".", 3)
    while len(fqdn) < 4:
        fqdn.append("")

    if "nid" == fqdn[0][:3]:
        hostname = "Lumi"
    elif "lvt.dkrz.de" == fqdn[1] + "." + fqdn[2] + "." + fqdn[3]:
        hostname = "Levante"
    elif "jupiter" == parts[1]:
        hostname = "Jupiter"
    else:
        hostname = "unknown"
        raise ValueError("aero_provider: Host cannot be detected")

    return hostname


def filename_dflt(dataPath, fileRoot, band):
    filename = dataPath + fileRoot + "_" + band + "_rast.nc"
    return filename


def filename_year(dataPath, fileRoot, band, year):
    if year >= 1845 and year <= 2100:
        filename = (
            dataPath + fileRoot + "_" + band + "_" + str(year) + "_rast.nc"
        )
        return filename
    else:
        raise ValueError("aero_provider: ", filename, " not available!")


if get_hostname() == "Levante":
    dataPath = "/pool/data/ICON/grids/public/mpim/independent/aerosol_kinne/"

if get_hostname() == "Lumi":
    dataPath = "/appl/local/climatedt/pool/data/ICON/grids/public/mpim/common/aerosol_kinne/"

if get_hostname() == "Jupiter":
    dataPath = "/p/data1/slmet/model_data/ICON/common/aerosol_kinne/"

fileRoot = "aeropt_kinne"

iso_data_interval = "P1M"

NAMELIST = sys.argv[1]
nml_fname = glob(NAMELIST)[0]
nml = f90nml.read(nml_fname)


iforcing = nml["run_nml"]["iforcing"]

if iforcing < 2 or iforcing > 3:
    print("aero_provider: iforcing ", iforcing, "is not supported", nml_fname)
    sys.exit(1)

if iforcing == 2:

    try:
        irad_aero = nml["aes_rad_nml"]["aes_rad_config"][0]["irad_aero"]
        if VERBOSE > 0:
            print(
                "aero_provider: found irad_aero =",
                irad_aero,
                "in",
                nml_fname,
                flush=True,
            )
    except:
        print("aero_provider: irad_aero for aes not set in", nml_fname)
        sys.exit(1)

    try:
        iso_coupling_interval = nml["aes_phy_nml"]["aes_phy_config"][0][
            "dt_rad"
        ]
        if VERBOSE > 0:
            print(
                "aero_provider: found dt_rad =",
                iso_coupling_interval,
                "in",
                nml_fname,
                flush=True,
            )
    except:
        print("aero_provider: dt_rad for aes not set in", nml_fname)
        sys.exit(1)

    try:
        lyr_perp = nml["aes_rad_nml"]["aes_rad_config"][0]["lyr_perp"]
        if VERBOSE > 0:
            print(
                "aero_provider: found lyr_perp =",
                lyr_perp,
                "in",
                nml_fname,
                flush=True,
            )
    except:
        print("aero_provider: lyr_perp not set in", nml_fname)
        lyr_perp = False

    if lyr_perp:

        try:
            yr_perp = nml["aes_rad_nml"]["aes_rad_config"][0]["yr_perp"]
            if VERBOSE > 0:
                print(
                    "aero_provider: found yr_perp =",
                    yr_perp,
                    "in",
                    nml_fname,
                    flush=True,
                )
        except:
            print("aero_provider: yr_perp for aes not set in", nml_fname)
            sys.exit(1)

if iforcing == 3:

    try:
        irad_aero = nml["radiation_nml"]["irad_aero"]
        if VERBOSE > 0:
            print(
                "aero_provider: found irad_aero =",
                irad_aero,
                "in",
                nml_fname,
                flush=True,
            )
    except:
        print("aero_provider: irad_aero for nwp not set in", nml_fname)
        sys.exit(1)

    dt_rad = 86400.0
    dr = pd.Timedelta(dt_rad, "s")
    iso_coupling_interval = dr.isoformat()
    if VERBOSE > 0:
        print(
            "aero_provider: found dt_rad =", dt_rad, "in", nml_fname, flush=True
        )
        print(
            "aero_provider: dt_rad converted to",
            iso_coupling_interval,
            flush=True,
        )

    try:
        nwp_dtime = nml["run_nml"]["dtime"]
        if VERBOSE > 0:
            print(
                "aero_provider: found dtime =",
                nwp_dtime,
                "in",
                nml_fname,
                flush=True,
            )

    except:
        print("aero_provider: dtime for nwp not set in", nml_fname)
        sys.exit(1)

    lyr_perp = False

if (
    irad_aero != 12
    and irad_aero != 13
    and irad_aero != 15
    and irad_aero != 18
    and irad_aero != 19
):
    raise ValueError(
        "aero_provider: irad_aero =", irad_aero, "is not supported"
    )

try:
    coupled_to_aero = nml["coupling_mode_nml"]["coupled_to_aero"]
    if VERBOSE > 0:
        print(
            "aero_provider: found coupled_to_aero =",
            coupled_to_aero,
            "in",
            nml_fname,
            flush=True,
        )
except:
    print("aero_provider: coupled_to_aero not set in", nml_fname)
    coupled_to_aero = False

if not coupled_to_aero:
    raise ValueError(
        "aero_provider: coupled_to_aero = .FALSE. cannot be used when running aero_provider "
    )

if irad_aero == 12 or irad_aero == 18 or irad_aero == 19:
    scenario = "picontrol"  # Kinne background from 1850

if irad_aero == 13 or irad_aero == 15:
    scenario = "historical"
    if lyr_perp:
        raise ValueError(
            "aero_provider: perpetual year conflicts with irad_aero"
        )

if lyr_perp:
    scenario = "perpetual"
    perpetual_year = yr_perp

switch_year = False

watch_this = 0
watch_other = 0

data_interval = isodate.parse_duration(iso_data_interval)
coupling_interval = isodate.parse_duration(iso_coupling_interval)

if not DRYRUN:
    yac = YAC()
    def_calendar(Calendar.PROLEPTIC_GREGORIAN)
    comp = yac.def_comp(f"aero_provider")

filename = filename_dflt(dataPath, fileRoot, "lw_b16_coa")

if VERBOSE > 0:
    print("aero_provider: Reading \n", filename, flush=True)

lw_b16_coa = xr.open_dataset(filename, decode_times=False)

filename = filename_dflt(dataPath, fileRoot, "sw_b14_coa")

if VERBOSE > 0:
    print("aero_provider: Reading \n", filename, flush=True)

sw_b14_coa = xr.open_dataset(filename, decode_times=False)

deg2rad = np.pi / 180
lon = deg2rad * sw_b14_coa["lon"]
lat = deg2rad * sw_b14_coa["lat"]

if VERBOSE > 0:
    print("aero_provider: size of lw  is", lw_b16_coa.sizes["lnwl"], flush=True)
    print("aero_provider: size of sw  is", sw_b14_coa.sizes["lnwl"], flush=True)
    print("aero_provider: size of lev is", lw_b16_coa.sizes["lev"], flush=True)

if not DRYRUN:

    grid = Reg2dGrid("aero_grid", lon, lat, cyclic=[True, False])
    grid.set_global_index(np.arange(grid.nbr_corners), Location.CORNER)

    points = grid.def_points(Location.CORNER, lon, lat)

    aod_lw_b16_coa_field = Field.create(
        "aod_lw_b16_coa",
        comp,
        points,
        lw_b16_coa.sizes["lnwl"],
        iso_coupling_interval,
        TimeUnit.ISO_FORMAT,
    )
    ssa_lw_b16_coa_field = Field.create(
        "ssa_lw_b16_coa",
        comp,
        points,
        lw_b16_coa.sizes["lnwl"],
        iso_coupling_interval,
        TimeUnit.ISO_FORMAT,
    )
    aer_lw_b16_coa_field = Field.create(
        "aer_lw_b16_coa",
        comp,
        points,
        lw_b16_coa.sizes["lev"],
        iso_coupling_interval,
        TimeUnit.ISO_FORMAT,
    )

    aod_sw_b14_coa_field = Field.create(
        "aod_sw_b14_coa",
        comp,
        points,
        sw_b14_coa.sizes["lnwl"],
        iso_coupling_interval,
        TimeUnit.ISO_FORMAT,
    )
    ssa_sw_b14_coa_field = Field.create(
        "ssa_sw_b14_coa",
        comp,
        points,
        sw_b14_coa.sizes["lnwl"],
        iso_coupling_interval,
        TimeUnit.ISO_FORMAT,
    )
    asy_sw_b14_coa_field = Field.create(
        "asy_sw_b14_coa",
        comp,
        points,
        sw_b14_coa.sizes["lnwl"],
        iso_coupling_interval,
        TimeUnit.ISO_FORMAT,
    )

    aod_sw_b14_fin_field = Field.create(
        "aod_sw_b14_fin",
        comp,
        points,
        sw_b14_coa.sizes["lnwl"],
        iso_coupling_interval,
        TimeUnit.ISO_FORMAT,
    )
    ssa_sw_b14_fin_field = Field.create(
        "ssa_sw_b14_fin",
        comp,
        points,
        sw_b14_coa.sizes["lnwl"],
        iso_coupling_interval,
        TimeUnit.ISO_FORMAT,
    )
    asy_sw_b14_fin_field = Field.create(
        "asy_sw_b14_fin",
        comp,
        points,
        sw_b14_coa.sizes["lnwl"],
        iso_coupling_interval,
        TimeUnit.ISO_FORMAT,
    )
    aer_sw_b14_fin_field = Field.create(
        "aer_sw_b14_fin",
        comp,
        points,
        sw_b14_coa.sizes["lev"],
        iso_coupling_interval,
        TimeUnit.ISO_FORMAT,
    )

if not DRYRUN:
    yac.enddef()

if not DRYRUN:
    start_date = isodate.parse_datetime(yac.start_datetime)
    end_date = isodate.parse_datetime(yac.end_datetime)
else:
    start_date = isodate.parse_datetime("1970-12-28T00:00:00.000")
    end_date = isodate.parse_datetime("1971-01-02T00:00:00.000")

model_date = start_date

if scenario == "picontrol":
    filename = filename_year(dataPath, fileRoot, "sw_b14_fin", 1850)
if scenario == "perpetual":
    filename = filename_year(dataPath, fileRoot, "sw_b14_fin", perpetual_year)
if scenario == "historical":
    filename = filename_year(dataPath, fileRoot, "sw_b14_fin", start_date.year)

if VERBOSE > 0:
    print("aero_provider: Reading \n", filename, "for", scenario, flush=True)

this_sw_b14_fin = xr.open_dataset(filename, decode_times=False)
other_sw_b14_fin = this_sw_b14_fin

if VERBOSE > 0:
    print("aero_provider: Reading", filename, flush=True)

while model_date < end_date:

    this_date = model_date
    this_month = this_date.month
    this_year = this_date.year

    if this_month < 12:
        mid_of_this_month = (
            date(this_year, this_month + 1, 1) - date(this_year, this_month, 1)
        ).days * 43200
    else:
        mid_of_this_month = 31 * 43200

    sec_of_this_month = (
        (this_date.day - 1) * 86400
        + this_date.hour * 3600
        + this_date.minute * 60
        + this_date.second
    )

    if sec_of_this_month <= mid_of_this_month:
        switch = -1
    else:
        switch = 1

    other_date = this_date + switch * data_interval
    other_month = other_date.month
    other_year = other_date.year

    if scenario == "historical":

        if other_year != this_year and not switch_year:
            filename = filename_year(
                dataPath, fileRoot, "sw_b14_fin", other_year
            )
            other_sw_b14_fin = xr.open_dataset(filename, decode_times=False)
            switch_year = True
            if VERBOSE > 0:
                print(
                    "aero_provider: Opened other dataset for year",
                    other_year,
                    flush=True,
                )

        if other_year == this_year and switch_year:
            switch_year = False
            other_sw_b14_fin = this_sw_b14_fin
            if VERBOSE > 0:
                print(
                    "aero_provider: Complete switch to year ",
                    other_year,
                    flush=True,
                )

    if other_month < 12:
        mid_of_other_month = (
            date(this_year, other_month + 1, 1)
            - date(this_year, other_month, 1)
        ).days * 43200
    else:
        mid_of_other_month = 31 * 43200

    other_wght = (
        switch
        * (sec_of_this_month - mid_of_this_month)
        / (mid_of_this_month + mid_of_other_month)
    )
    this_wght = 1 - other_wght

    if this_month != watch_this:
        if watch_this == 0:

            if VERBOSE > 1:
                print(
                    "aero_provider: Reading this month", this_month, flush=True
                )

            this_aod_lw_b16_coa = lw_b16_coa["aod"][
                {"time": this_month - 1}
            ].values.reshape((1, -1))
            this_ssa_lw_b16_coa = lw_b16_coa["ssa"][
                {"time": this_month - 1}
            ].values.reshape((1, -1))
            this_aer_lw_b16_coa = lw_b16_coa["z_aer_coarse_mo"][
                {"time": this_month - 1}
            ].values.reshape((1, -1))

            this_aod_sw_b14_coa = sw_b14_coa["aod"][
                {"time": this_month - 1}
            ].values.reshape((1, -1))
            this_ssa_sw_b14_coa = sw_b14_coa["ssa"][
                {"time": this_month - 1}
            ].values.reshape((1, -1))
            this_asy_sw_b14_coa = sw_b14_coa["asy"][
                {"time": this_month - 1}
            ].values.reshape((1, -1))

            this_aod_sw_b14_fin = this_sw_b14_fin["aod"][
                {"time": this_month - 1}
            ].values.reshape((1, -1))
            this_ssa_sw_b14_fin = this_sw_b14_fin["ssa"][
                {"time": this_month - 1}
            ].values.reshape((1, -1))
            this_asy_sw_b14_fin = this_sw_b14_fin["asy"][
                {"time": this_month - 1}
            ].values.reshape((1, -1))
            this_aer_sw_b14_fin = this_sw_b14_fin["z_aer_fine_mo"][
                {"time": this_month - 1}
            ].values.reshape((1, -1))

        else:

            if VERBOSE > 1:
                print(
                    "aero_provider: Flipp this month to", this_month, flush=True
                )

            this_aod_lw_b16_coa = other_aod_lw_b16_coa
            this_ssa_lw_b16_coa = other_ssa_lw_b16_coa
            this_aer_lw_b16_coa = other_aer_lw_b16_coa

            this_aod_sw_b14_coa = other_aod_sw_b14_coa
            this_ssa_sw_b14_coa = other_ssa_sw_b14_coa
            this_asy_sw_b14_coa = other_asy_sw_b14_coa

            this_aod_sw_b14_fin = other_aod_sw_b14_fin
            this_ssa_sw_b14_fin = other_ssa_sw_b14_fin
            this_asy_sw_b14_fin = other_asy_sw_b14_fin
            this_aer_sw_b14_fin = other_aer_sw_b14_fin

        watch_this = this_month

    if other_month != watch_other:

        if VERBOSE > 1:
            print("aero_provider: Reading other month", other_month, flush=True)

        other_aod_lw_b16_coa = lw_b16_coa["aod"][
            {"time": other_month - 1}
        ].values.reshape((1, -1))
        other_ssa_lw_b16_coa = lw_b16_coa["ssa"][
            {"time": other_month - 1}
        ].values.reshape((1, -1))
        other_aer_lw_b16_coa = lw_b16_coa["z_aer_coarse_mo"][
            {"time": other_month - 1}
        ].values.reshape((1, -1))

        other_aod_sw_b14_coa = sw_b14_coa["aod"][
            {"time": other_month - 1}
        ].values.reshape((1, -1))
        other_ssa_sw_b14_coa = sw_b14_coa["ssa"][
            {"time": other_month - 1}
        ].values.reshape((1, -1))
        other_asy_sw_b14_coa = sw_b14_coa["asy"][
            {"time": other_month - 1}
        ].values.reshape((1, -1))

        other_aod_sw_b14_fin = other_sw_b14_fin["aod"][
            {"time": other_month - 1}
        ].values.reshape((1, -1))
        other_ssa_sw_b14_fin = other_sw_b14_fin["ssa"][
            {"time": other_month - 1}
        ].values.reshape((1, -1))
        other_asy_sw_b14_fin = other_sw_b14_fin["asy"][
            {"time": other_month - 1}
        ].values.reshape((1, -1))
        other_aer_sw_b14_fin = other_sw_b14_fin["z_aer_fine_mo"][
            {"time": other_month - 1}
        ].values.reshape((1, -1))

        watch_other = other_month

    if VERBOSE > 2:
        print(
            "aero_provider: ",
            this_date.isoformat(),
            " wght ",
            "%8.6f" % this_wght,
            "*",
            "%2i" % this_month,
            "+",
            "%8.6f" % other_wght,
            "*",
            "%2i" % other_month,
            flush=True,
        )

    aod_lw_b16_coa_array = (
        this_wght * this_aod_lw_b16_coa + other_wght * other_aod_lw_b16_coa
    )
    ssa_lw_b16_coa_array = (
        this_wght * this_ssa_lw_b16_coa + other_wght * other_ssa_lw_b16_coa
    )
    aer_lw_b16_coa_array = (
        this_wght * this_aer_lw_b16_coa + other_wght * other_aer_lw_b16_coa
    )

    aod_sw_b14_coa_array = (
        this_wght * this_aod_sw_b14_coa + other_wght * other_aod_sw_b14_coa
    )
    ssa_sw_b14_coa_array = (
        this_wght * this_ssa_sw_b14_coa + other_wght * other_ssa_sw_b14_coa
    )
    asy_sw_b14_coa_array = (
        this_wght * this_asy_sw_b14_coa + other_wght * other_asy_sw_b14_coa
    )

    aod_sw_b14_fin_array = (
        this_wght * this_aod_sw_b14_fin + other_wght * other_aod_sw_b14_fin
    )
    ssa_sw_b14_fin_array = (
        this_wght * this_ssa_sw_b14_fin + other_wght * other_ssa_sw_b14_fin
    )
    asy_sw_b14_fin_array = (
        this_wght * this_asy_sw_b14_fin + other_wght * other_asy_sw_b14_fin
    )
    aer_sw_b14_fin_array = (
        this_wght * this_aer_sw_b14_fin + other_wght * other_aer_sw_b14_fin
    )

    if not DRYRUN:
        aod_lw_b16_coa_field.put(aod_lw_b16_coa_array)
        ssa_lw_b16_coa_field.put(ssa_lw_b16_coa_array)
        aer_lw_b16_coa_field.put(aer_lw_b16_coa_array)

        aod_sw_b14_coa_field.put(aod_sw_b14_coa_array)
        ssa_sw_b14_coa_field.put(ssa_sw_b14_coa_array)
        asy_sw_b14_coa_field.put(asy_sw_b14_coa_array)

        aod_sw_b14_fin_field.put(aod_sw_b14_fin_array)
        ssa_sw_b14_fin_field.put(ssa_sw_b14_fin_array)
        asy_sw_b14_fin_field.put(asy_sw_b14_fin_array)
        aer_sw_b14_fin_field.put(aer_sw_b14_fin_array)

    if iforcing == 3:
        # mimicks strange intervals in nwp case
        if UPDATE_CASE == 1:
            model_date = model_date
            UPDATE_CASE = 2
        elif UPDATE_CASE == 2:
            UPDATE_CASE = 3
            model_date = (
                np.datetime64(model_date) + np.timedelta64(nwp_dtime, "s")
            ).astype(datetime)
        else:
            model_date = model_date + coupling_interval

    else:
        model_date = model_date + coupling_interval

if VERBOSE > 0:
    print("aero_provider: Done \n", flush=True)
