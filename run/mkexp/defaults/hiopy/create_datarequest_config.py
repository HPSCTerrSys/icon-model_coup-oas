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
import isodate


def params_3d(
    defaults_2d={},
    name_of_zaxis="",
    positive_direction="up",
    is_pressure_level=False,
    levels=[],
    source_yac_collection_size=0,
    lev_indices_to_save=[],
):
    defaults_3d = {
        **defaults_2d,
        "levels": levels,
        "name_of_level": name_of_zaxis,
        "positive_direction": positive_direction,
        "is_pressure_level": is_pressure_level,
    }
    if source_yac_collection_size > len(levels):
        defaults_3d = {
            **defaults_3d,
            "yac_collection_size": source_yac_collection_size,
            "indices_to_save": lev_indices_to_save,
        }
    return defaults_3d


def params_2d(
    icon_source,
    icon_grid,
    time_res,
    time_method,
    frac_mask=None,
):
    return {
        "yac_source_comp": icon_source,
        "yac_source_grid": icon_grid,
        "time_method": time_method,
        "time_res_needed": time_res,
        "frac_mask": frac_mask,
    }


def atmo_2d_defaults(time_res, time_method):
    return params_2d("atm_output", "icon_atmos_grid", time_res, time_method)


def ocean_2d_defaults(time_res, time_method, frac_mask="ocean_frac_mask_sfc"):
    return params_2d(
        "oce_output",
        "icon_ocean_grid",
        time_res,
        time_method,
        frac_mask=frac_mask,
    )


def parse_data_request(user_data_request_config, time_aggs_needed, grid_info):
    data_request = {}
    for data_group_name in user_data_request_config:
        data_group = user_data_request_config[data_group_name]
        if (
            ".hide" in data_group and data_group[".hide"].lower() == "true"
        ) or "variables" not in data_group:
            continue

        assert type(data_group["variables"]) == list

        grid_details = {}
        if grid_info["type"] == "healpix":
            grid_details = {
                "type": "healpix",
                "zoom_levels": [int(z) for z in grid_info["zoom_levels"]],
            }
        elif "icon" in grid_info["type"]:
            if "atm" in data_group_name or "land" in data_group_name:
                grid_details["type"] = grid_info["type"] + "_atmo"
                grid_details["path_to_grid_file"] = grid_info["atmo"][
                    "path_to_grid_file"
                ]
            elif "ocean" in data_group_name or "hamocc" in data_group_name:
                grid_details["type"] = grid_info["type"] + "_ocean"
                grid_details["path_to_grid_file"] = grid_info["ocean"][
                    "path_to_grid_file"
                ]

        settings_2d = {}
        if "atm" in data_group_name or "land" in data_group_name:
            settings_2d = atmo_2d_defaults(
                data_group["time_res"], data_group["time_method"]
            )
        elif "ocean" in data_group_name or "hamocc" in data_group_name:
            settings_2d = ocean_2d_defaults(
                data_group["time_res"], data_group["time_method"]
            )

        if "2d" in data_group_name:
            data_request[data_group_name] = {
                "variables": data_group["variables"],
                "grid_details": grid_details,
                "settings": settings_2d,
            }

        if "3d" in data_group_name:
            is_pressure_level = False
            if (
                "is_pressure_level" in data_group
                and data_group["is_pressure_level"].lower() == "true"
            ):
                is_pressure_level = True

            # land sea mask requires special handling for ocean 3d fields
            if "ocean" in data_group_name or "hamocc" in data_group_name:
                nr_levels_needed = len(
                    list(data_group.get("levels_needed", []))
                )
                if nr_levels_needed == 0:
                    nr_levels_needed = data_group.get(
                        "total_levels_available", 0
                    )
                if nr_levels_needed == 0:
                    raise ValueError(
                        f"Number of levels needed for ocean variable group {data_group_name} cannot be determined from the config. Please specify either 'levels_needed' or 'total_levels_available'."
                    )

                settings_2d = ocean_2d_defaults(
                    data_group["time_res"],
                    data_group["time_method"],
                    frac_mask=f"ocean_frac_mask_{nr_levels_needed}_lev",
                )

            data_request[data_group_name] = {
                "variables": data_group["variables"],
                "grid_details": grid_details,
                "settings": params_3d(
                    settings_2d,
                    data_group["name_of_zaxis"],
                    data_group["positive_direction"],
                    is_pressure_level,
                    levels=list(data_group.get("levels_needed", [])),
                    source_yac_collection_size=int(
                        data_group.get("total_levels_available", 0)
                    ),
                    lev_indices_to_save=data_group.get(
                        "level_indices_to_save", []
                    ),
                ),
            }

    for v in data_request.values():
        v["settings"]["time_res_needed"] = list(
            {v["settings"]["time_res_needed"]} | set(time_aggs_needed)
        )

    return data_request
