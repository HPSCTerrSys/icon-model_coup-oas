```{eval-rst}
:orphan:
```

(ref_infrastructure_testing_unit_testing)=
# Unit Testing in ICON

The ICON build system supports compiling and running unit tests.

## What is Unit Testing

:::topic
_Module testing (or unit testing) is a process of testing the individual subprograms, subroutines, classes, or procedures in a program. More specifically, rather than initially testing the program as a whole, testing is first focused on the smaller building blocks of the program. &mdash; The Art of Software Testing - G.J. Myers, T. Badgett, C. Sandler_
:::

## Running Unit Tests

To build and run all tests, use the following command:
```sh
make -j8 check-icon
```

What it does:
- builds the bundled packages, compiles the test programs, and the relevant subset of the ICON source files in parallel;
- runs the default set of unit tests via the {{ '[test driver]({}/utils/test_driver.py)'.format(base_url) }} sequentially.

Sample output:
```console
$ make check-icon
...
PASS: test_divide_cell
SKIP: test_divide_cell.mpi
PASS: test_index_list
PASS: test_insert_dimension
PASS: test_kind
PASS: test_slice_array
PASS: test_start_mpi.mpi
SKIP: test_start_mpi.nompi
PASS: test_sync_sums.mpi
============================================================================
Testsuite summary for icon 2025.04
============================================================================
# TOTAL: 9
# PASS:  7
# SKIP:  2
# FAIL:  0
============================================================================
...
```

Each unit-test is classified as either `PASS`, `SKIP`, or `FAIL`. The command exits with a non-zero code if any of the tests fail.

### Special Makefile Variables

This subsection describes several makefile variables that control the unit test execution. They can be either exported to the environment before running `make` or specified as arguments of the command.

By default, the test driver produces colorized output if its standard output stream is connected to a terminal. The colorization can be enforced or suppressed by setting the `ICON_COLOR_TESTS` variable to the `always` and `no` values, respectively. Any other value of the variable results in the default behavior.

By default, the test driver does not dump log files of the individual tests to the standard output. The default behavior can be overridden by setting the `ICON_VERBOSE_TESTS` variable to a non-empty value. In that case, if any of the tests fail, the output of all unsuccessful tests is emitted to the standard output.

By default, the test driver runs all tests. It is possible to run only a specific subset of tests by setting the `TESTS` variable to a space-separated list of test names.

## Adding Unit Tests

To add a test, you need to implement a Fortran program that performs the required checks. The name of the source file of the program must have the `test_` prefix and the `.f90` suffix (extension). The file itself must be saved in the {{ '[`test/unit-tests`]({}/test/unit-tests)'.format(base_url) }} directory (subdirectories are allowed).

The status of a unit test depends on the exit code of its test program: the zero code indicates success, code 77 classifies the test as skipped, which does not affect the final status of the whole test suite, and any other code is considered as a failure of the test.

To ensure the portability across different Fortran compilers, it is recommended to terminate test programs by calling the `test_pass`, `test_skip` and `test_fail` subroutines of the {{ '[`mo_test_common`]({}/test/unit-tests/common/mo_test_common.f90)'.format(base_url) }} module.

For example:

```fortran
PROGRAM test_example
    USE mo_test_common, ONLY: test_fail, test_pass, test_skip

    IMPLICIT NONE

    LOGICAL :: lsuccess

#ifdef NOMPI
    ! Skip the test if ICON is configured without the MPI support:
    CALL test_skip()
#endif

    ! Set the variable to the result of the test:
    lsuccess = ...

    IF (lsuccess) THEN
      CALL test_pass()
    ELSE
      CALL test_fail()
    END IF
END PROGRAM test_example
```

If a test program cannot be built in a certain configuration, you can indicate that to the build system of ICON by adding a special infix to the name of the source file:
- `.mpi.`: the test program is built only when the MPI support is enabled;
- `.nompi.`: the test program is built only when the MPI support is disabled.

If the building of a test program is omitted, the respective unit test is classified as skipped with no attempt to run any executable. For example, if the source file of a test has name `test_example.mpi.f90`, it will not be compiled when ICON is configured with the `--disable-mpi` option and the `test_example` test will get the `SKIP` mark when running the test suite.

You can control the execution of a test program by adding the runtime configuration to its source file. All lines among the first 20 that start with the `! 4EXEC-JSON` prefix are interpreted as a single mult-line JSON string. The test driver converts the string into a dictionary and the contents of the dictionary affect the way the test program is run.

Currently, the test driver honors the following entries:
- `MPI_NTASKS`: the integer number of MPI ranks to run the test program with; if the key is specified (regardless of the value), the test program will be run using the `MPI_LAUNCH` command detected by the {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON; if the configure script fails to detect a working interactive MPI launcher command, the test will get the `SKIP` mark when running the test suite;
- `ENV`: the dictionary of environment variables that must be set when running the test program.

For example, if a test program must be run with two MPI ranks and the `OMP_NUM_THREADS` variable set to five, you can add the following comment to the source file of the program (must be specified in the first 20 lines of the file):

```fortran
! 4EXEC-JSON {
! 4EXEC-JSON   "MPI_NTASKS": 2,
! 4EXEC-JSON   "ENV": {
! 4EXEC-JSON     "OMP_NUM_THREADS" : 5
! 4EXEC-JSON   }
! 4EXEC-JSON }
```

If you need more control over how a test program is run, you can implement the required logic as a portable shell script. If the directory containing the `<name>.f90` source file of the test program also contains the `<name>.sh` file, the test driver will ignore the test runtime configuration specified in the source file and run the script using the shell interpreter detected by the configure script of ICON. Before running the script, the driver will export the following environment variables:
- `TEST_PROG`: the absolute path to the test program;
- `MPI_LAUNCH`: the interactive MPI launcher command.

Since the script is executed by passing its absolute path to the shell interpreter, it does not need to have the [`shebang`](https://en.wikipedia.org/wiki/Shebang_(Unix)) and the execute permissions. The exit code convention for the test script is the same as for the test programs: the zero code indicates success, code 77 classifies the test as skipped, and any other code is considered as a failure. Below, you can find an example of a test script:

```sh
# Skip the test if the MPI launcher command is not available:
test -n "${MPI_LAUNCH}" || exit 77

for n in 1 2 3 4; do
  ${MPI_LAUNCH} -n ${n} "${TEST_PROG}"
done
```
