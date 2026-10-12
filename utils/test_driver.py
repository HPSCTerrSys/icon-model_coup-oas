#!/bin/sh

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

""":"
for cmd in python3 python; do
  if command -v > /dev/null "${cmd}"; then
    exec "${cmd}" "$0" "$@"
  fi
done
echo "Error: could not find a python interpreter!" >&2
exit 1
":"""

__doc__ = "Unit Test Driver"

import argparse
import collections
import json
import os
import subprocess
import sys

_term_style_codes = {
    "blue": "\033[1;34m",
    "bold": "\033[1m",
    "green": "\033[0;32m",
    "red": "\033[0;31m",
}

_test_status_term_styles = {
    "FAIL": "red",
    "PASS": "green",
    "SKIP": "blue",
}


def get_test_config(src_file, head_size=20, json_tag="! 4EXEC-JSON"):
    """
    Parses the source file of a test program and returns a dictionary of its
    runtime configuration parameters. The dictionary is expected to be specified
    in the source file as a JSON string. Each line of the JSON must appear in
    the first `head_size` lines of the source file and prefixed with the
    `json_tag` string. If no such strings are found in the file, the result is
    an empty Python dictionary.

    :param src_file: path to the source file of the test program
    :param head_size: number of lines in the source file to process
    :param json_tag:
        marker indicating that the rest of the line must be interpreted as part
        of the multi-line JSON with the runtime configuration parameters
    :return:
        Python dictionary with runtime configuration parameters of the test
        program
    """
    with open(src_file) as f:
        test_config_json = "".join(
            line[len(json_tag) :].strip()
            for line in f.readlines()[:head_size]
            if line.startswith(json_tag)
        )
        return json.loads(test_config_json) if test_config_json else {}


def get_test_cmd_and_env(
    name, src_dir, test_subdir, test_exec_ext, mpi_launch, shell
):
    """
    Generates a shell command and a dictionary of extra environment variables
    that are required to run a test.

    :param name: name of the test
    :param src_dir: path to the root source directory of the project
    :param test_subdir:
        relative (to `src_dir`) path to the root directory containing the source
        files of the test programs
    :param test_exec_ext:
        filename extensions (with the leading dot) of the test binary
        executables
    :param mpi_launch: string with the interactive MPI launcher command
    :param shell: string with the shell interpreter command
    :return:
        tuple of a string containing the shell command of the test and a
        dictionary of extra environment variables to be exported before running
        the command
    """
    executable = os.path.abspath(
        os.path.join(test_subdir, name + test_exec_ext)
    )

    script = os.path.join(src_dir, test_subdir, name + ".sh")
    if os.path.isfile(script):
        return "{0} -ex {1}".format(shell, script), {
            "TEST_PROG": os.path.abspath(executable),
            "MPI_LAUNCH": mpi_launch,
        }

    src_file = os.path.join(src_dir, test_subdir, name + ".f90")
    if os.path.isfile(src_file):
        config = get_test_config(src_file)

        mpi_task_count = config.get("MPI_NTASKS", None)

        if not (mpi_task_count is None or mpi_launch):
            return None, None

        env = {k: str(v) for k, v in config.get("ENV", {}).items()} or None

        if mpi_task_count is None:
            return executable, env
        else:
            return (
                "{0} -n {1} {2}".format(mpi_launch, mpi_task_count, executable),
                env,
            )

    return executable, None


def run_cmd(cmd, extra_env):
    """
    Runs a shell command in the modified environment.

    :param cmd: string containing the shell command
    :param extra_env:
        dictionary with environment variables to be exported before running the
        command
    :return:
        tuple of an integer that equals to the exit code of the command and a
        string to be dumped to the log file
    """
    env = None
    if extra_env:
        env = dict(os.environ)
        env.update(extra_env)

    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=env,
        shell=True,
    )

    out, _ = p.communicate()

    return p.returncode, "{0}+ {1}\n{2}".format(
        (
            "\n".join(
                "+ export {0}='{1}'".format(k, v.replace("'", "'\"'\"'"))
                for k, v in extra_env.items()
            )
            + "\n"
            if extra_env
            else ""
        ),
        cmd,
        out.decode("utf-8"),
    )


def colorize_string(string, color):
    """
    Wraps the provided string with the ANSI color escape sequences.

    :param string:
        arbitrary string
    :param color:
        name of the color or style to apply to the string (see keys of the
        `_term_style_codes` dictionary for the possible values)
    :return:
        copy of the input string surrounded with ANSI color escape sequences
    """
    return "{0}{1}\033[m".format(_term_style_codes[color], string)


def is_test_suite_successful(summary):
    """
    Checks whether the test suite run is successful based on the summary.

    :param summary:
        dictionary that maps possible test final statuses to the number of tests
        that have the respective status
    :return:
        `True` if the test suite run is considered successful and `False`
        otherwise
    """
    return summary["FAIL"] == 0


def generate_test_suite_summary_table(summary, colorize):
    """
    Generates a table with the test suite summary to be dumped to the logs and
    the terminal output.

    :param summary:
        dictionary that maps possible test final statuses to the number of tests
        that have the respective status
    :param colorize:
        boolean denoting whether the output string must be colorized
    :return: string containing the test suite summary table
    """
    table = [("TOTAL", sum(summary.values()), "bold")]
    table.extend(
        (status, summary[status], _test_status_term_styles[status])
        for status in ["PASS", "SKIP", "FAIL"]
    )

    longest_key_len = max(len(item[0]) for item in table)

    lines = []
    for key, value, style in table:
        line = "# {0}:{1} {2}".format(
            key, " " * (longest_key_len - len(key)), value
        )
        if colorize and value > 0:
            line = colorize_string(line, style)
        lines.append(line)

    return "\n".join(lines)


def generate_test_suite_log(title, summary, skipped_or_failed_test_logs):
    """
    Generates the contents of the test suite log file.

    :param title: string containing the title of the log file
    :param summary:
        dictionary that maps possible test final statuses to the number of tests
        that have the respective status
    :param skipped_or_failed_test_logs:
        enumerable of strings with the logs of unsuccessful tests (including the
        skipped ones)
    :return: string with the contents of the test suite log file
    """
    log = """\
{ruler}
{title}
{ruler}

{table}

{logs}
""".format(
        ruler="=" * len(title),
        title=title,
        table=generate_test_suite_summary_table(summary, colorize=False),
        logs="\n".join(skipped_or_failed_test_logs),
    )

    return log


def generate_test_suite_report(title, summary, failure_epilogue, colorize):
    """
    Generates the contents of the test suite report to be dumped to the terminal
    output.

    :param title: string containing the title of the report
    :param summary:
        dictionary that maps possible test final statuses to the number of tests
        that have the respective status
    :param failure_epilogue:
        string to be dumped to the terminal output in case any of the tests
        failed
    :param colorize:
        boolean denoting whether the output string must be colorized
    :return: the contents of the test suite report
    """
    if is_test_suite_successful(summary):
        color = "green"
        epilogue = None
    else:
        color = "red"
        epilogue = failure_epilogue

    ruler = "=" * 76
    if colorize:
        title = colorize_string(title, color)
        ruler = colorize_string(ruler, color)
        epilogue = colorize_string(epilogue, color) if epilogue else epilogue

    report = """\
{ruler}
{title}
{ruler}
{table}
{epilogue}
""".format(
        ruler=ruler,
        title=title,
        table=generate_test_suite_summary_table(summary, colorize),
        epilogue="{0}\n{1}\n{0}".format(ruler, epilogue) if epilogue else ruler,
    )

    return report


def parse_args():
    """
    Parses the command line arguments of the test driver.

    :return: instance of the argparse.Namespace object
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
    )
    parser.add_argument(
        "--tests",
        nargs="*",
        default=[],
        help="space-separated list of tests names to run",
    )
    parser.add_argument(
        "--skip-tests",
        nargs="*",
        default=[],
        help="space-separated list of tests names to skip",
    )
    parser.add_argument(
        "--src-dir",
        metavar="SRC_DIR",
        help="root source code directory",
    )
    parser.add_argument(
        "--test-subdir",
        help="root test subdirectory relative to SRC_DIR",
    )
    parser.add_argument(
        "--test-exec-ext",
        default="",
        help="filename extension of the test binary executables",
    )
    parser.add_argument(
        "--shell",
        help="shell interpreter to run the test execution scripts",
    )
    parser.add_argument(
        "--mpi-launch",
        help="interactive MPI launcher command to run the MPI tests",
    )
    parser.add_argument(
        "--colorize",
        choices=["always", "no", "auto"],
        default="auto",
        help="colorize output (default: `%(default)s`)",
    )
    parser.add_argument(
        "--output-on-failure",
        action="store_true",
        help="print the test suite log if any of the tests fails",
    )
    parser.add_argument(
        "--package-name-version",
        help="name and version of the tested package to be used in logs and "
        "reports",
    )

    args = parser.parse_args()

    args.skip_tests = set(args.skip_tests)

    args.colorize = (
        sys.stdout.isatty()
        if args.colorize == "auto"
        else args.colorize == "always"
    )

    return args


def main():
    """
    The main function of the test driver.

    :return: `True` if the test suite run is successful and `False` otherwise
    """
    args = parse_args()

    test_suite_summary = collections.defaultdict(lambda: 0)

    skipped_or_failed_test_logs = []

    for test_name in args.tests:
        test_exitcode = test_output = None

        if test_name not in args.skip_tests:
            test_cmd, test_env = get_test_cmd_and_env(
                test_name,
                args.src_dir,
                args.test_subdir,
                args.test_exec_ext,
                args.mpi_launch,
                args.shell,
            )
            if test_cmd:
                test_exitcode, test_output = run_cmd(test_cmd, test_env)

        if test_exitcode == 0:
            test_status = "PASS"
        elif test_exitcode in (77, None):
            test_status = "SKIP"
        else:
            test_status = "FAIL"

        test_log = "{0}{1} {2} {3}\n".format(
            test_output or "",
            test_status,
            test_name,
            (
                "(without running)"
                if test_exitcode is None
                else "(exit status: {0})".format(test_exitcode)
            ),
        )

        with open(os.path.join(args.test_subdir, test_name + ".log"), "w") as f:
            f.write(test_log)

        if test_status != "PASS":
            test_log_title = "{0}: {1}".format(test_status, test_name)
            skipped_or_failed_test_logs.extend(
                [
                    test_log_title,
                    "=" * len(test_log_title),
                    test_log,
                ]
            )

        test_suite_summary[test_status] += 1
        sys.stdout.write(
            "{0}: {1}\n".format(
                (
                    colorize_string(
                        test_status, _test_status_term_styles[test_status]
                    )
                    if args.colorize
                    else test_status
                ),
                test_name,
            )
        )

    test_suite_log_file = os.path.join(args.test_subdir, "test-suite.log")
    test_suite_log = generate_test_suite_log(
        "   {0}{1}   ".format(
            (
                args.package_name_version + " "
                if args.package_name_version
                else ""
            ),
            test_suite_log_file,
        ),
        test_suite_summary,
        skipped_or_failed_test_logs,
    )

    with open(test_suite_log_file, "w") as f:
        f.write(test_suite_log)

    test_suite_success = is_test_suite_successful(test_suite_summary)

    if not test_suite_success and args.output_on_failure:
        sys.stdout.write(test_suite_log)

    sys.stdout.write(
        generate_test_suite_report(
            "Testsuite summary{0}".format(
                " for " + args.package_name_version
                if args.package_name_version
                else ""
            ),
            test_suite_summary,
            "See {0}".format(test_suite_log_file),
            args.colorize,
        )
    )

    return test_suite_success


if __name__ == "__main__":
    exit(int(not main()))
