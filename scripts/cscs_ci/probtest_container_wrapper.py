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
import copy
import os
import shutil
import subprocess
import sys
from pathlib import Path

import click
import toml

# Add yaml_experiment_test_processor to path
script_dir = Path(__file__).resolve().parent
icon_dir = script_dir.parents[1]
sys.path.insert(0, str(icon_dir / "scripts/experiments"))
from yaml_experiment_test_processor import ExperimentTestCollection

# Set probtest container executable
probtest_container = "srun --container-writable --environment=probtest"
probtest = "/probtest/probtest.py"
# Full bash command that runs inside the container

# Set EDF_PATH if not set to know the path to the toml file
if "EDF_PATH" not in os.environ:
    home = os.environ["HOME"]
    os.environ["EDF_PATH"] = f"{home}/.edf"


def get_env_var(var_name):
    if var_name not in os.environ:
        raise EnvironmentError(
            f"Missing required environment variable: {var_name}"
        )
    else:
        return os.environ[var_name]


# Two different names are used for the reference file.
# Ensure that the required one is available and copy in case it is not.
def check_ref_file(experiment, build_dir):
    reference_file_path = build_dir / f"stats_{experiment}_ref.csv"
    if not os.path.exists(reference_file_path):
        alt_reference_file_path = build_dir / f"{experiment}_reference.csv"
        if os.path.exists(alt_reference_file_path):
            shutil.copy(alt_reference_file_path, reference_file_path)
        else:
            raise FileNotFoundError(
                f"No reference file at {reference_file_path} or {alt_reference_file_path}"
            )


def prepare_probtest_mounts(base_toml: Path, mounts: list[str]) -> Path:
    """Return a probtest_mounts.toml file containing the given mount entries,
    only adding entries that are not already present.
    """
    probtest_toml_mounts = base_toml.with_name("probtest_mounts.toml")
    shutil.copy(base_toml, probtest_toml_mounts)

    with base_toml.open("r") as f:
        config = toml.load(f)

    config_mounts = copy.deepcopy(config)
    existing = set(config_mounts.get("mounts", []))

    new_mounts = [m for m in mounts if m not in existing]
    if existing or new_mounts:
        config_mounts["mounts"] = list(existing | set(new_mounts))

    with probtest_toml_mounts.open("w") as f:
        toml.dump(config_mounts, f)

    return probtest_toml_mounts


def generate_stats_file(
    etc,
    parent_experiment,
    build_dir,
    member_name=None,
    stats_file_name=None,
    file_id=None,
):
    """Generates the stats file for the given experiment."""
    if member_name:
        experiment = member_name
    else:
        experiment = parent_experiment

    model_output_dir = build_dir / "experiments" / experiment

    # Get file IDs from YAML files if not give
    if not file_id:
        file_id = etc.get_file_ids_for_exp_as_string(parent_experiment)

    # Determine stats file name
    if not stats_file_name:
        stats_file_name = build_dir / f"stats_{experiment}.csv"

    # Commands to run inside container
    probtest_init = f"{probtest} init {file_id}"
    probtest_stats = (
        f"{probtest} stats --no-ensemble "
        f"--stats-file-name {stats_file_name} "
        f"--model-output-dir {model_output_dir}"
    )

    # If the model_output_dir is not a subdirectory of EDF_PATH, we need to add a mount
    EDF_PATH = Path(get_env_var("EDF_PATH"))
    if EDF_PATH not in model_output_dir.parents:
        mounts = [f"{model_output_dir}:{model_output_dir}"]
        probtest_toml = EDF_PATH / "probtest.toml"
        probtest_toml_mounts = prepare_probtest_mounts(probtest_toml, mounts)
        container_cmd = (
            "srun --container-writable --environment=probtest_mounts"
        )
    else:
        container_cmd = probtest_container

    subprocess.run(
        f'{container_cmd} bash -c "{probtest_init} && {probtest_stats}"',
        check=True,
        shell=True,
    )

    # Delete probtest_mounts.toml after use
    if "probtest_toml_mounts" in locals():
        os.remove(probtest_toml_mounts)


def generate_tolerance_file(
    etc, experiment, build_dir, file_id=None, member_ids=None
):
    """Generates the tolerance file from the stats files of an ensemble run using members specified in YAML files."""

    # Get file IDs from YAML files if not give
    if not file_id:
        file_id = etc.get_file_ids_for_exp_as_string(experiment)
    else:
        file_id = "--file-id " + file_id[0][0] + " " + file_id[0][1]

    # Get ensemble member numbers
    if not member_ids:
        bb_name = get_env_var("BB_NAME")
        reference_builder = bb_name.replace("_gpu", "_cpu")
        member_ids = etc.get_ensemble_num_for_exp_as_string(
            experiment, reference_builder
        )

    # Initialize probtest
    probtest_init = (
        probtest
        + " init "
        + file_id
        + " --experiment-name "
        + experiment
        + " --member-ids "
        + member_ids
    )

    # Generate tolerance file
    ensemble_files = build_dir / f"stats_{experiment}_{{member_id}}.csv"
    tolerance_files = build_dir / f"{experiment}_tolerance.csv"

    # Ensure correct reference file is available
    check_ref_file(experiment, build_dir)

    probtest_tolerance = (
        probtest
        + " tolerance"
        + " --ensemble-files "
        + str(ensemble_files)
        + " --tolerance-files "
        + str(tolerance_files)
    )

    subprocess.run(
        f'{probtest_container} bash -c "{probtest_init} && {probtest_tolerance}"',
        check=True,
        shell=True,
    )

    # Rename file to name used in pp.collect_tolerance_hashes
    shutil.move(
        build_dir / f"stats_{experiment}_ref.csv",
        build_dir / f"{experiment}_reference.csv",
    )


def run_tolerance_check(
    etc,
    experiment,
    file_id,
    current_files,
    reference_files,
    tolerance_files,
    factor,
):
    """Runs the tolerance check using the given input file, reference file, and tolerance file."""
    bb_name = get_env_var("BB_NAME")
    EDF_PATH = Path(get_env_var("EDF_PATH"))

    current_file = Path(current_files).resolve()
    reference_file = Path(reference_files).resolve()
    tolerance_file = Path(tolerance_files).resolve()

    # Mount input, references and tolerances files in probtest container
    probtest_toml = EDF_PATH / f"probtest.toml"
    probtest_toml = EDF_PATH / "probtest.toml"
    mounts = [
        f"{current_file}:{current_file}",
        f"{reference_file}:{reference_file}",
        f"{tolerance_file}:{tolerance_file}",
    ]

    probtest_toml_mounts = prepare_probtest_mounts(probtest_toml, mounts)
    probtest_container_mounts = (
        "srun --container-writable --environment=probtest_mounts"
    )

    # Get file IDs from YAML files if not given
    if not file_id:
        file_id = etc.get_file_ids_for_exp_as_string(experiment)
    else:
        file_id = "--file-id " + file_id[0][0] + " " + file_id[0][1]

    # Initialize probtest
    probtest_init = probtest + " init " + file_id

    if not factor:
        factor = etc.get_param_for_exp_by_machine_as_string(
            experiment, "tolerance_factor", bb_name
        )
    # Validate generated stats against reference
    probtest_check = (
        probtest
        + " check"
        + " --current-files "
        + str(current_file)
        + " --reference-files "
        + str(reference_file)
        + " --tolerance-files "
        + str(tolerance_file)
        + " --factor "
        + str(factor)
    )

    # Use toml file with mounts
    probtest_container_mounts = (
        "srun --container-writable --environment=probtest_mounts"
    )
    subprocess.run(
        f'{probtest_container_mounts} bash -c "{probtest_init} && {probtest_check}"',
        check=True,
        shell=True,
    )

    # Delete probtest_mounts.toml after use
    os.remove(probtest_toml_mounts)


def select_members(etc, experiment, build_dir, file_id):
    """Selects the members for tolerance generation from 49 members"""
    # Get file IDs from YAML files if not given
    if not file_id:
        file_id = etc.get_file_ids_for_exp_as_string(experiment)
    else:
        file_id = "--file-id " + file_id[0][0] + " " + file_id[0][1]

    # Initialize probtest
    probtest_init = (
        probtest + " init " + file_id + " --experiment-name " + experiment
    )

    ensemble_files = build_dir / f"stats_{experiment}_{{member_id}}.csv"
    selected_members_file_name = (
        build_dir / f"{experiment}_selected_members.csv"
    )
    tolerance_file = build_dir / f"{experiment}_tolerance.csv"
    reference_file_path = build_dir / f"stats_{experiment}_ref.csv"

    # Ensure correct reference file is available
    check_ref_file(experiment, build_dir)

    # Select members
    probtest_select = (
        probtest
        + " select-members"
        + " --ensemble-files "
        + str(ensemble_files)
        + " --selected-members-file-name "
        + str(selected_members_file_name)
        + " --max-member-count 20"
        # GitLab-CI does not allow more than 50 dependents for a single job
        # use 49 members instead
        + " --total-member-count 49"
        + " --min-factor 5"
        + " --max-factor 50"
        + " --tolerance-files "
        + str(tolerance_file)
    )

    # Generate tolerance file from all members
    probtest_tolerance = (
        probtest
        + " tolerance"
        + " --ensemble-files "
        + str(ensemble_files)
        + " --tolerance-files "
        + str(tolerance_file)
        + " --member-ids '"
        + ",".join(str(i) for i in range(1, 50))
        + "'"
    )

    subprocess.run(
        f'{probtest_container} bash -c "{probtest_init}  && {probtest_select} && {probtest_tolerance}"',
        check=True,
        shell=True,
    )

    # Rename file to name used in pp.collect_tolerance_hashes
    shutil.move(
        build_dir / f"stats_{experiment}_ref.csv",
        build_dir / f"{experiment}_reference.csv",
    )


def ensemble_member(
    etc,
    parent_experiment,
    build_dir,
    member_name=None,
    member_id=None,
    perturb_amplitude=None,
    file_id=None,
):
    """
    Creates a runscript from the parent experiment with the respective amplitude and seed
    for perturbation and generates the stats file after running the ensemble member.
    """
    if member_name:
        experiment = member_name
    else:
        experiment = os.getenv("MEMBER_NAME")
        if not experiment:
            experiment = parent_experiment

    run_dir = build_dir / "run"
    perturbed_run_script_name = f"exp.{experiment}.run"
    run_script_name = f"exp.{parent_experiment}.run"

    # Create member runscript
    if parent_experiment != experiment:
        if not member_id:
            member_id = int(get_env_var("MEMBER_ID"), 0)
        if not perturb_amplitude:
            perturb_amplitude = get_env_var("PERTURB_AMPLITUDE")
        rhs_new = f"{{seed}},{perturb_amplitude}"

        # Generate probtest.json
        probtest_init = f"{probtest} init"
        # Generate runscript for member
        probtest_member = (
            probtest
            + " run-ensemble"
            + " --dry"
            + " --run-dir"
            + f" {run_dir}"
            + " --run-script-name"
            + f" {run_script_name}"
            + " --perturbed-run-script-name"
            + f" {perturbed_run_script_name}"
            + " --member-ids"
            + f" {member_id}"
            + " --experiment-name"
            + f" {parent_experiment}"
            + " --perturbed-experiment-name"
            + f" {experiment}"
            + " --rhs-new"
            + f" {rhs_new}"
        )
        subprocess.run(
            f'{probtest_container} bash -c "{probtest_init} && {probtest_member}"',
            check=True,
            shell=True,
        )
        run_script_path = build_dir / "run" / perturbed_run_script_name
        stats_file_name = (
            build_dir / f"stats_{parent_experiment}_{member_id}.csv"
        )
    else:
        run_script_path = build_dir / "run" / run_script_name
        stats_file_name = build_dir / f"stats_{experiment}_ref.csv"

    # Create experiments folder
    os.makedirs(build_dir / "experiments", exist_ok=True)

    # Run member
    os.chmod(run_script_path, 0o755)
    subprocess.run([run_script_path], check=True, stdout=subprocess.PIPE)

    generate_stats_file(
        etc,
        parent_experiment,
        build_dir=build_dir,
        member_name=experiment,
        stats_file_name=stats_file_name,
        file_id=file_id,
    )


def run_ensemble(
    etc,
    experiment,
    build_dir,
    file_id=None,
    member_ids=None,
    perturb_amplitude=1e-14,
):
    """Runs the ensemble with perturbated members and creates the stats files."""
    if not member_ids:
        # Get ensemble member numbers
        bb_name = get_env_var("BB_NAME")
        reference_builder = bb_name.replace("_gpu", "_cpu")
        member_ids = etc.get_ensemble_num_for_exp(experiment, reference_builder)
    else:
        member_ids = [int(x) for x in member_ids.split(",")]

    # Run reference
    ensemble_member(
        etc,
        parent_experiment=experiment,
        build_dir=build_dir,
        file_id=file_id,
    )

    # Run ensemble members
    for i in member_ids:
        ensemble_member(
            etc,
            parent_experiment=experiment,
            build_dir=build_dir,
            member_name=f"{experiment}_member_id_{i}",
            member_id=i,
            perturb_amplitude=perturb_amplitude,
            file_id=file_id,
        )


@click.command()
@click.argument(
    "task",
    type=click.Choice(
        [
            "stats",
            "tolerance",
            "check",
            "select-members",
            "ensemble-member",
            "ensemble",
        ]
    ),
)
@click.argument("experiment", required=True)
@click.option(
    "--build-dir",
    help="ICON build directory.",
)
@click.option(
    "--member-name", help="Name of member in case of an ensemble run."
)
@click.option("--stats-file-name", help="Name of stats file to be created.")
@click.option(
    "--file-id",
    nargs=2,
    type=str,
    multiple=True,
    metavar="FILE_TYPE FILE_PATTERN",
    help="Unique identifier and file pattern.",
)
@click.option(
    "--member-ids",
    type=str,
    help="Comma separated list of members (e.g. '1,3,14').",
)
@click.option(
    "--current-files", help="Path to stats file to run the tolerance check."
)
@click.option(
    "--reference-files",
    help="Path to the reference file for the tolerance check.",
)
@click.option(
    "--tolerance-files",
    help="Path to the tolerance file for the tolerance check.",
)
@click.option(
    "--factor",
    type=float,
    help="Tolerance factor.",
)
@click.option(
    "--perturb-amplitude",
    type=float,
    default=1e-14,
    help="Perturbation amplitude.",
)
def main(
    task,
    experiment,
    build_dir,
    member_name,
    stats_file_name,
    file_id,
    member_ids,
    current_files,
    reference_files,
    tolerance_files,
    factor,
    perturb_amplitude,
):
    if not build_dir:
        build_dir = icon_dir
    else:
        build_dir = Path(build_dir)

    etc = ExperimentTestCollection()
    if task == "stats":
        generate_stats_file(
            etc, experiment, build_dir, member_name, stats_file_name, file_id
        )
    elif task == "tolerance":
        generate_tolerance_file(etc, experiment, build_dir, file_id, member_ids)
    elif task == "check":
        run_tolerance_check(
            etc,
            experiment,
            file_id,
            current_files,
            reference_files,
            tolerance_files,
            factor,
        )
    elif task == "select-members":
        select_members(etc, experiment, build_dir, file_id)
    elif task == "ensemble-member":
        ensemble_member(etc, experiment, build_dir, member_name, file_id)
    elif task == "ensemble":
        run_ensemble(
            etc,
            experiment,
            build_dir,
            file_id,
            member_ids,
            perturb_amplitude=1e-14,
        )


if __name__ == "__main__":
    main()
