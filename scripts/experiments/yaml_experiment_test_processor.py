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
import shutil
import subprocess
import sys
import warnings

import numpy as np
import yaml

sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "buildbot_scripts")
)
from exp_utils import adddep, addexp, rmexp


class ExperimentTestCollection:
    def __init__(self):
        test_yml = os.path.join(
            os.path.dirname(__file__), "../experiments/all_tests.yml"
        )
        exp_yml = os.path.join(
            os.path.dirname(__file__), "../experiments/all_experiments.yml"
        )
        self.items = self._expand_tests_with_experiment_config(
            self._load_yaml_with_key(test_yml, "tests"),
            self._load_yaml_with_key(exp_yml, "experiments"),
        )
        self.defaults = self._load_defaults()

        # supported machines and builders, any other machine or builder will raise an error
        # BuildBot
        self.supported_machines_BB = {"balfrin", "horeka", "levante"}
        self.supported_builders_BB = {
            "balfrin_cpu_nvidia",
            "balfrin_gpu_nvidia",
            "balfrin_cpu_nvidia_mixed",
            "balfrin_gpu_nvidia_mixed",
            "alps_mch_test_cpu",
            "alps_mch_test_gpu",
            "horeka_cpu_nvhpc",
            "horeka_gpu_nvhpc",
        }
        # CSCS-CI
        self.supported_machines_CSCS_CI = {"santis"}
        self.supported_builders_CSCS_CI = {
            "santis_cpu_nvhpc",
            "santis_gpu_nvhpc",
        }

        self.supported_builders = self.supported_builders_BB.union(
            self.supported_builders_CSCS_CI
        )
        self.supported_machines = self.supported_machines_BB.union(
            self.supported_machines_CSCS_CI
        )

    def get_items_by_tag(self, tag_name, mode):
        items_by_tag = []
        # select-members and tolerance use a common tag and run under the probtest mode
        if mode == "probtest":
            tag_name = mode
        for item in self.items["tests"]:
            if tag_name in item["tags"]:
                items_by_tag.append(item)
        if not items_by_tag:
            raise Exception(f"Tag {tag_name} not found in items")
        return {"tests": items_by_tag}

    def get_items_by_name(self, name):
        return {"tests": [self.get_item_by_name(name)]}

    def get_items_by_names(self, name):
        return {"tests": [self.get_item_by_name(n.strip()) for n in name]}

    def get_item_by_name(self, name):
        item_by_name = next(
            (i for i in self.items["tests"] if i["name"] == name), None
        )
        if not item_by_name:
            warnings.warn(
                f"Entry with name {name} not found in items. Using default values.",
                UserWarning,
            )
        return item_by_name

    def print_file_ids_for_exp(self, name):
        # print to stdout for usage in bash scripts
        print(self.get_file_ids_for_exp_as_string(name))

    def print_param_for_exp_by_machine(self, exp, param, bb_name):
        # print to stdout for usage in bash scripts
        print(self.get_param_for_exp_by_machine_as_string(exp, param, bb_name))

    def print_checksuite_param_for_exp(self, param, exp):
        # print to stdout for usage in bash scripts
        print(self._get_checksuite_param_for_exp_as_string(param, exp))

    def check(self):
        basepath = os.path.join(os.path.dirname(__file__), "../../run")

        # Check that tests has entries
        if not self.items["tests"]:
            raise Exception("No entries in tests")

        checks = set()
        for experiment in self.items["tests"]:

            # Check that 'check' is unique
            if experiment["check"] in checks:
                raise Exception(f"Duplicate check: {experiment['check']}")
            else:
                checks.add(experiment["check"])

            # Check that 'check' is a file
            if not os.path.isfile(os.path.join(basepath, experiment["check"])):
                raise Exception(f"Check {experiment['check']} is not a file")

            for machine in experiment.get("machines", []):
                # Check that machine is in supported machines
                if machine["name"] not in self.supported_machines:
                    raise Exception(
                        f"Machine {machine['name']} not in supported machines"
                    )

                # Check that 'include_only' and 'exclude' are not both present
                if "include_only" in machine and "exclude" in machine:
                    raise Exception(
                        f"Machine {machine['name']} in experiment {experiment['name']} "
                        "has both 'include_only' and 'exclude'"
                    )
                # Check that 'include_only' has only supported builders
                if "include_only" in machine:
                    for builder in machine["include_only"]:
                        if builder not in self.supported_builders:
                            raise Exception(
                                f"Builder {builder} in machine {machine['name']} "
                                "is not in supported builders"
                            )
                # Check that 'exclude' has only supported builders
                if "exclude" in machine:
                    for builder in machine["exclude"]:
                        if builder not in self.supported_builders:
                            raise Exception(
                                f"Builder {builder} in machine {machine['name']} "
                                "is not in supported builders"
                            )

            # Check that 'refgen' is present when 'probtest' is in 'tags'
            if "probtest" in experiment["tags"] and "refgen" not in experiment:
                raise Exception(
                    f"Experiment {experiment['name']} has 'probtest' in 'tags' "
                    "but no 'refgen'"
                )
            # Check that builders for which there is an ensemble_num give are in refgen
            if "refgen" in experiment and "ensemble_num" in experiment:
                for builder in experiment["ensemble_num"]:
                    builder_name = list(builder.keys())[0]
                    if builder_name not in experiment["refgen"]:
                        raise Exception(
                            f"Builder {builder_name} in ensemble_num but not in refgen "
                            f"for experiment {experiment['name']}"
                        )

            # Check that builders for which there is a checksuite_mode "t" are in refgen
            checksuite_modes = experiment.get("checksuite_modes", [])
            for mode in checksuite_modes:
                for machine, cm in mode.items():
                    if "t" in cm:
                        refgen = experiment.get("refgen", [])
                        if not any(
                            machine.lower() in ref.lower() for ref in refgen
                        ):
                            raise Exception(
                                f"Builder {machine} has checksuite mode {cm}, but no corresponding refgen "
                                f"for experiment {experiment['name']}"
                            )

    def _load_defaults(self):
        exp_defaults = self._load_yaml(
            os.path.join(
                os.path.dirname(__file__),
                "../experiments/defaults_experiments.yml",
            )
        )
        test_defaults = self._load_yaml(
            os.path.join(
                os.path.dirname(__file__), "../experiments/defaults_tests.yml"
            )
        )
        return {**exp_defaults, **test_defaults}

    def _load_yaml(self, file_path):
        with open(file_path, "r") as file:
            data = yaml.safe_load(file)
        return data

    def _load_yaml_with_key(self, file_path, valid_key):
        data = self._load_yaml(file_path)

        if "include" in data:
            for include_file in data["include"]:
                include_file_path = os.path.join(
                    os.path.dirname(file_path), include_file
                )
                include_data = self._load_yaml_with_key(
                    include_file_path, valid_key
                )  # recursive call
                if valid_key not in include_data:
                    raise Exception(
                        f"The included file {include_file} does not contain "
                        f"{valid_key} key"
                    )
                if valid_key not in data:
                    data[valid_key] = []

                data[valid_key].extend(include_data[valid_key])

        return data

    def _get_perturb_amplitude_as_string(self, name, member_type):
        return str(self._get_perturb_amplitude(name, member_type))

    def _get_perturb_amplitude(self, name, member_type):
        default = next(
            (
                item.get(member_type)
                for item in self.defaults.get("perturb_amplitude")
                if member_type in item
            )
        )

        exp_data = self.get_item_by_name(name)
        tolexp = exp_data.get("tolerance") if exp_data else None
        # use default value if there is no 'tolerance' section in yml
        if tolexp is None:
            return default
        else:
            actual = next(
                (
                    item.get(member_type)
                    for item in tolexp.get("perturb_amplitude", [])
                    if member_type in item
                ),
                None,
            )
            return actual if actual is not None else default

    def get_file_ids_for_exp_as_string(self, name):
        return " ".join(self._get_file_ids_for_exp_as_list(name))

    def _get_file_ids_for_exp_as_list(self, name):
        file_ids = self.defaults["file_id"]
        exp_data = self.get_item_by_name(name)
        tolexp = exp_data.get("tolerance") if exp_data else None
        if tolexp:
            file_id_list = tolexp.get("file_id")
            if file_id_list is not None:
                file_ids = []
                for id_dict in file_id_list:
                    for file_group, glob_patterns in id_dict.items():
                        for glob_pattern in glob_patterns:
                            file_ids.extend(
                                ["--file-id", file_group, glob_pattern]
                            )
        return file_ids

    def get_ensemble_num_for_exp_as_string(self, name, bb_name):
        return ",".join(map(str, self.get_ensemble_num_for_exp(name, bb_name)))

    def get_ensemble_num_for_exp(self, name, bb_name):
        num = self.defaults["ensemble_num"]
        test_item = self.get_item_by_name(name)
        if test_item is not None:
            if "ensemble_num" in test_item:
                ensemble_num = test_item["ensemble_num"]
                for item in ensemble_num:
                    if bb_name in item:
                        num = item[bb_name]
                        break
        return num

    def get_param_for_exp_by_machine_as_string(self, exp, param, bb_name):
        return str(self._get_param_for_exp_by_machine(exp, param, bb_name))

    def _get_checksuite_param_for_exp_as_string(self, param, exp):
        return " ".join(self._get_checksuite_param_for_exp(param, exp))

    def _get_param_for_exp_by_machine(self, exp, param, bb_name):
        default = self.defaults[param]
        tests_item = self.get_item_by_name(exp)
        if tests_item is not None:
            if param in tests_item:
                cm = tests_item[param]
                for item in cm:
                    # check if the key is in the bb_name
                    # i.e. balfrin in balfrin_cpu_nvidia
                    for key, value in item.items():
                        # lowercase, i.e BALFRIN_CPU_nvidia -> balfrin_cpu_nvidia
                        if key in bb_name.lower():
                            return value
        return default

    def _get_checksuite_param_for_exp(self, param, exp):
        return self.get_item_by_name(exp).get(param)

    def _expand_tests_with_experiment_config(self, tests, experiments):
        expanded = []
        for test in tests["tests"]:
            name = test["name"]
            exp_config = next(
                (e for e in experiments["experiments"] if e["name"] == name),
                None,
            )
            if not exp_config:
                raise Exception(f"Experiment {name} not found in yml-config")

            expanded.append({**exp_config, **test})

        return {"tests": expanded}


class BuildBotInterface(ExperimentTestCollection):
    def __init__(self, tag_name):
        super().__init__()
        self.tag_name = tag_name
        self.bb_name = os.getenv("BB_NAME")
        if not self.bb_name:
            raise Exception("Environment variable BB_NAME is not set")

    def items_to_bb(self):
        self.items["tests"] = self._get_experiments_with_supported_machines()
        if self.tag_name == "tolerance" or self.tag_name == "tolerance-update":
            self._register_tolerance_list()
        elif self.tag_name == "select-members":
            self._register_select_members_list()
        else:
            self._register_default_list()

    def _get_experiments_with_supported_machines(self):
        valid_experiments = []
        # Add experiments on supported machines only
        for exp in self.items["tests"]:
            # Check intersection of supported machines and machines in the experiment
            valid_machines = self.supported_machines_BB.intersection(
                {machine["name"] for machine in exp.get("machines", [])}
            )

            if valid_machines:
                exp.update(
                    {
                        "machines": [
                            machine
                            for machine in exp.get("machines", [])
                            if machine["name"] in valid_machines
                        ]
                    }
                )
                valid_experiments.append(exp)
        return valid_experiments

    def _register_tolerance_list(self):
        all_tolerance_exps = []
        all_pp_gentol = []
        for exp in self.items["tests"]:
            for builder in exp["refgen"]:
                if builder != self.bb_name:
                    continue  # we are not interested in other builders
                pp_gentol = self._register_tolerance_for_current_builder(
                    exp, builder
                )
                # safety check
                if pp_gentol is not None:
                    all_pp_gentol.append(pp_gentol)
                    all_tolerance_exps.append(exp["name"])

        # now we need to collect the hashes
        pp = ["tolerance/pp.collect_tolerance_hashes"]
        self._add_to_bb_list(
            pp,
            builders=[self.bb_name],
            runflags=f"tolerance_experiments={','.join(all_tolerance_exps)}",
        )
        self._add_dep_to_bb_list(
            builders=[self.bb_name],
            from_experiment=pp,
            to_experiment=all_pp_gentol,
        )

    def _register_select_members_list(self):
        all_select_members_exps = []
        all_pp_selmem = []
        for exp in self.items["tests"]:
            for builder in exp["refgen"]:
                if builder != self.bb_name:
                    continue  # we are not interested in other builders
                pp_selmem = self._register_select_members_for_current_builder(
                    exp, builder
                )
                # safety check
                if pp_selmem is not None:
                    all_pp_selmem.append(pp_selmem)
                    all_select_members_exps.append(exp["name"])

        # now we need to collect the hashes and the selected members
        pp = ["tolerance/pp.collect_tolerance_hashes"]
        self._add_to_bb_list(
            pp,
            builders=[self.bb_name],
            runflags=f"tolerance_experiments={','.join(all_select_members_exps)}",
        )
        self._add_dep_to_bb_list(
            builders=[self.bb_name],
            from_experiment=pp,
            to_experiment=all_pp_selmem,
        )

        pp = ["tolerance/pp.collect_selected_members"]
        self._add_to_bb_list(
            pp,
            builders=[self.bb_name],
            runflags=f"tolerance_experiments={','.join(all_select_members_exps)}",
        )
        self._add_dep_to_bb_list(
            builders=[self.bb_name],
            from_experiment=pp,
            to_experiment=all_pp_selmem,
        )

    def _register_default_list(self):
        for exp in self.items["tests"]:
            for machine in exp.get("machines", []):
                if "include_only" in machine:
                    self._add_to_bb_list(
                        [exp["check"]],
                        builders=machine["include_only"],
                        runflags=machine.get("runflags"),
                    )
                elif "exclude" in machine:
                    self._add_to_bb_list(
                        [exp["check"]],
                        machines=[machine["name"]],
                        runflags=machine.get("runflags"),
                    )
                    self._remove_from_bb_list(
                        [exp["check"]], builders=machine.get("exclude")
                    )
                else:
                    self._add_to_bb_list(
                        [exp["check"]],
                        machines=[machine["name"]],
                        runflags=machine.get("runflags"),
                    )

    def _add_to_bb_list(
        self, experiment_name, builders=None, machines=None, runflags=None
    ):
        runflags = self._convert_types_for_bb(runflags)
        addexp(
            experiment_name,
            builders,
            None,
            None,
            machines,
            runflags,
            self.tag_name,
        )

    def _add_dep_to_bb_list(
        self,
        from_experiment=None,
        to_experiment=None,
        from_builder=None,
        to_builder=None,
        builders=None,
        machines=None,
    ):
        adddep(
            from_builder,
            from_experiment,
            to_builder,
            to_experiment,
            builders,
            None,
            None,
            machines,
            self.tag_name,
        )

    def _remove_from_bb_list(
        self, experiment_name, builders=None, machines=None
    ):
        rmexp(experiment_name, builders, None, None, machines, self.tag_name)

    def _convert_types_for_bb(self, runflags=None):
        if runflags:
            return dict(item.split("=") for item in "".join(runflags).split())

    def _add_probtest_ensemble_for_member_ids(
        self, basedir, exp, builder, member_ids
    ):
        # set path to probtest entry script
        PROBTEST = os.path.join(basedir, "externals/probtest/probtest.py")
        # define member_type (must be in sync with add_refgen_routines)
        if "mixed" in builder:
            member_type = "mixed"
        else:
            member_type = "double"

        probtest_config = os.path.join(basedir, exp["name"] + "-config.json")
        # Set environment variable for probtest to identify json file for each experiment by probtest
        os.environ["PROBTEST_CONFIG"] = probtest_config
        # initialize probtest namelist (most of it is unused, but makes life easier)
        cmd = [
            "python3",
            PROBTEST,
            "init",
            "--config",
            probtest_config,
            "--codebase-install",
            basedir,
            "--experiment-name",
            exp["name"],
            "--member-type",
            member_type,
            "--perturb-amplitude",
            self._get_perturb_amplitude_as_string(exp["name"], member_type),
            "--member-ids",
            member_ids,
        ]
        cmd.extend(self._get_file_ids_for_exp_as_list(exp["name"]))
        subprocess.run(cmd, check=True)

        # create the runscripts for the ensemble (exp.<EXP>_seed_N)
        # needs to be overwritten from namelist because here we deal with the templates
        subprocess.run(
            [
                "python3",
                PROBTEST,
                "run-ensemble",
                "--dry",
                "--run-script-name",
                f"exp.{exp['name']}",
                "--perturbed-run-script-name",
                f"exp.{exp['name']}_member_id_{{member_id}}",
            ],
            check=True,
        )

        self._add_to_bb_list(
            [f"exp.{exp['name']}"],
            builders=[builder],
            runflags="tolerance=true",
        )

        perturbed_experiments = []
        member_ids = map(int, member_ids.split(","))
        for member_id in member_ids:
            perturbed_experiments.append(
                f"exp.{exp['name']}_member_id_{member_type}_{member_id}"
            )
        return perturbed_experiments

    def _register_tolerance_for_current_builder(self, exp, builder):
        basedir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "../..")
        )
        # Add probtest ensemble runs
        member_ids = self.get_ensemble_num_for_exp_as_string(
            exp["name"], self.bb_name
        )
        perturbed_experiments = self._add_probtest_ensemble_for_member_ids(
            basedir, exp, builder, member_ids
        )

        # Make a copy of the generic pp.generate_tolerance for each exp.
        # This copy is needed as unique identifier in adddep.
        pp_gentol = f"tolerance/pp.generate_tolerance_{exp['name']}"
        shutil.copy(
            os.path.join(basedir, "run/tolerance/pp.generate_tolerance"),
            os.path.join(basedir, f"run/{pp_gentol}"),
        )
        self._add_to_bb_list([pp_gentol], builders=[builder])

        self._add_to_bb_list(
            perturbed_experiments,
            builders=[builder],
            runflags="tolerance_run=true",
        )
        self._add_dep_to_bb_list(
            from_experiment=[pp_gentol],
            to_experiment=perturbed_experiments,
            builders=[builder],
        )
        return pp_gentol

    def _register_select_members_for_current_builder(self, exp, builder):
        basedir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "../..")
        )
        # Add probtest ensemble runs
        member_ids = ",".join(map(str, range(1, 50)))
        perturbed_experiments = self._add_probtest_ensemble_for_member_ids(
            basedir, exp, builder, member_ids
        )

        # Make a copy of the generic pp.select_members for each exp.
        # This copy is needed as unique identifier in adddep.
        pp_selmem = f"tolerance/pp.select_members_{exp['name']}"
        shutil.copy(
            os.path.join(basedir, "run/tolerance/pp.select_members"),
            os.path.join(basedir, f"run/{pp_selmem}"),
        )
        self._add_to_bb_list(
            [pp_selmem], builders=[builder], runflags="cpu_time=01:00:00"
        )

        self._add_to_bb_list(
            perturbed_experiments,
            builders=[builder],
            runflags="tolerance_run=true",
        )
        self._add_dep_to_bb_list(
            from_experiment=[pp_selmem],
            to_experiment=perturbed_experiments,
            builders=[builder],
        )
        return pp_selmem


class CscsCiInterface(ExperimentTestCollection):
    def __init__(self):
        super().__init__()

    def get_supported_builders(self):
        supported_builders = {
            "santis": ["santis_cpu_nvhpc", "santis_gpu_nvhpc"]
        }
        return supported_builders

    def get_items_for_builder(self, builder, mode):
        supported_builders = self.get_supported_builders()
        items = []
        for exp in self.items["tests"]:
            # in case of no machines section, we assume not needed to run on this builder
            for machine in exp.get("machines", []):
                name = machine["name"]

                # only consider the builder we are currently interested in
                if name not in builder:
                    continue
                elif (
                    "include_only" in machine
                    and builder not in machine["include_only"]
                ):
                    continue
                elif "exclude" in machine and builder in machine["exclude"]:
                    continue
                elif builder not in supported_builders[name]:
                    continue
                elif ("gpu" in builder) and (mode == "probtest"):
                    continue
                elif (mode == "probtest") and builder not in exp["refgen"]:
                    continue
                items.append(exp)

        return {"tests": items}

    def get_items_for_machine(self, machine):
        supported_builders = self.get_supported_builders()
        items = []
        for exp in self.items["tests"]:
            # in case of no machines section, we assume not needed to run on this machine
            for m in exp.get("machines", []):
                name = m["name"]

                # only consider the machine we are currently interested in
                if name != machine:
                    continue
                else:
                    items.append(exp)

        return {"tests": items}

    def items_to_cscs_ci(self, builder, tag_name, mode):
        pipeline = self._gen_header()
        pipeline["Build ICON"] = self._gen_step_build(builder)
        basedir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "../..")
        )

        for test in self.items["tests"]:
            pipeline[f'{test["name"]}'] = self._gen_step_build_for_test(
                test, builder
            )
            if mode == "probtest":
                pipeline.update(
                    self._gen_step_build_for_probtest(test, builder, tag_name)
                )

        if mode == "probtest":
            tests = self.items["tests"]
            pipeline.update(
                {
                    f"collect_tolerance_hashes": {
                        "extends": ".collect_tolerance_hashes_santis",
                        "variables": {
                            "BB_NAME": builder,
                            "tolerance_experiments": " ".join(
                                test["name"] for test in tests
                            ),
                        },
                    }
                }
            )
            if tag_name == "select-members":
                tests = self.items["tests"]
                pipeline.update(
                    {
                        f"collect_selected_members": {
                            "extends": ".collect_selected_members_santis",
                            "variables": {
                                "BB_NAME": builder,
                                "BB_SYSTEM": "santis",
                                "tolerance_experiments": " ".join(
                                    test["name"] for test in tests
                                ),
                            },
                        }
                    }
                )

        with open(os.path.join(basedir, "pipeline.yml"), "w") as outfile:
            yaml.dump(
                pipeline,
                outfile,
                default_flow_style=False,
                sort_keys=False,
                indent=2,
            )

    def _gen_step_build(self, builder):
        images = {
            "santis_cpu_nvhpc": ".build_santis_cpu_nvhpc",
            "santis_gpu_nvhpc": ".build_santis_gpu_nvhpc",
        }
        return {
            "extends": images[builder],
        }

    def _gen_header(self):
        return {
            "include": [
                ".gitlab/ci/cscs/remote.yml",
                ".gitlab/ci/cscs/recipes.yml",
            ],
            "variables": {
                "GIT_DEPTH": 1,
            },
            "stages": ["build", "run", "post"],
        }

    def _gen_step_build_for_probtest(self, test, builder, tag_name):
        member_type = "mixed" if "mixed" in builder else "double"

        if tag_name.startswith("tolerance"):
            ensemble_num = self.get_ensemble_num_for_exp(test["name"], builder)
            pipeline_name = "create_tolerance"
        elif tag_name == "select-members":
            ensemble_num = list(range(1, 50))
            pipeline_name = "select_members"

        output = {
            test["name"]: {
                "extends": ".run_ensemble_member_santis",
                "variables": {"EXPERIMENT": test["name"]},
            }
        }

        output.update(
            {
                f'{test["name"]}_member_id_{num}': {
                    "extends": ".run_ensemble_member_santis",
                    "variables": {
                        "EXPERIMENT": test["name"],
                        "MEMBER_NAME": f'{test["name"]}_member_id_{num}',
                        "MEMBER_ID": num,
                        "PERTURB_AMPLITUDE": self._get_perturb_amplitude(
                            test["name"], member_type
                        ),
                    },
                }
                for num in ensemble_num
            }
        )

        output.update(
            {
                f'{test["name"]}_{pipeline_name}': {
                    "extends": f".{pipeline_name}_santis",
                    "variables": {
                        "EXPERIMENT": test["name"],
                        "BB_NAME": builder,
                    },
                    "needs": [
                        f'{test["name"]}',
                        *[
                            f'{test["name"]}_member_id_{num}'
                            for num in ensemble_num
                        ],
                    ],
                }
            }
        )

        return output

    def _gen_step_build_for_test(self, test, builder):
        return {
            "extends": ".run_common_santis",
            "variables": {
                "EXPERIMENT": test["name"],
                "CHECK": test["check"],
                "BB_NAME": builder,
            },
        }


def set_mode(tag_name):
    # The tags "tolerance", "select-members" and "tolerance-update" have common functionality
    # and use the same tag "probtest" in the yaml lists.
    if (
        tag_name == "tolerance"
        or tag_name == "select-members"
        or tag_name == "tolerance-update"
    ):
        mode = "probtest"
    else:
        mode = "run_test"
    return mode


# main entrypoint for CSCS CI
def register_experiments_for_cscs_ci(tag_name, builder, exp=None):
    mode = set_mode(tag_name)
    cci = CscsCiInterface()
    if exp:
        if isinstance(exp, list):
            cci.items = cci.get_items_by_names(exp)
        else:
            cci.items = cci.get_items_by_name(exp)

    cci.items = cci.get_items_by_tag(tag_name, mode)

    cci.items = cci.get_items_for_builder(builder, mode)
    if not cci.items["tests"]:
        raise Exception(
            "No tests left after filtering, check GitLab CI variables CSCS_CI_EXPERIMENTS and CSCS_CI_TAG"
        )

    cci.items_to_cscs_ci(builder, tag_name, mode)


# main entrypoint for BuildBot
def register_experiments_for_bb(tag_name, exp=None):

    bbi = BuildBotInterface(tag_name)
    # only keep entry with name of single_exp
    if exp:
        bbi.items = bbi.get_items_by_name(exp)

    mode = set_mode(tag_name)
    # only keep the relevent entries for the tag
    bbi.items = bbi.get_items_by_tag(tag_name, mode)

    bbi.items_to_bb()


def cscs_ci_to_data(tag_name):
    """
    Returns a dictionary of machines with builders and experiments and their status (True/False).
    """
    cci = CscsCiInterface()
    mode = set_mode(tag_name)
    cci.items = cci.get_items_by_tag(tag_name, mode)
    supported_builders = cci.get_supported_builders()

    data = {"machines": list(supported_builders.keys())}

    for midx, builders in supported_builders.items():
        # Get experiments on this machine
        items_on_machine = cci.get_items_for_machine(midx)["tests"]
        experiments = [exp["check"] for exp in items_on_machine]

        # Get experiments per builder
        items_on_builder = {
            builder: cci.get_items_for_builder(builder, mode)["tests"]
            for builder in builders
        }

        # Build status matrix
        status = [
            [exp in items_on_builder[builder] for builder in builders]
            for exp in items_on_machine
        ]

        # Order experiments alphabetically
        experiments, status = (
            zip(*sorted(zip(experiments, status), key=lambda x: x[0]))
            if experiments
            else ([], [])
        )

        data[midx] = {
            "builders": builders,
            "experiments": experiments,
            "status": status,
        }

    return data


if __name__ == "__main__":
    tests = ExperimentTestCollection()
    tests.check()
