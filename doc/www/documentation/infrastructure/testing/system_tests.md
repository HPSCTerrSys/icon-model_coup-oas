```{eval-rst}
:orphan:
```

(ref_infrastructure_system_tests)=
# System Tests

<!-- EXTERNAL-CI-SYSTEM-TESTS -->

## ICON Development Checksuite

The {{ '[ICON development checksuite]({}/run/checksuite.icon-dev/icon-dev.checksuite)'.format(base_url) }} (`icon-dev.checksuite`) defines a set of generic system tests that can be applied to any ICON configuration/experiment. One or several of the checksuite flags (see definitions below) defined in the `check.<exp-name>` file indicates which of the system tests are applied to the experiment. `check.<exp-name>` may also be used to add specific input to the experiments. See a list of checksuite experiments {{ '[here]({}/run/checksuite.icon-dev)'.format(base_url) }}.

**Note:** Checksuite is only compatible with `make_runscripts`, but not with `mkexp` at the moment.

The central element of the checksuite is the *base* (`b`) test, which runs a simulation of the experiment and only fails if the simulation crashes (i.e., a **smoke test**). Most of the other tests use this first simulation (referred hereon as the *base simulation*) as a reference for comparison of results, see e.g. restart test below. Note that:
- If `check.<exp-name>` does not define `'b'` among its flags, the base simulation still runs in order to have a reference for other tests.
- If there are several tests defined in `check.<exp-name>`, the base simulation runs once and the same base simulation is used as reference for all tests.


### Performance Tests

- Performance (`p`) test: Measures the total runtime of the simulation and compares it to a stored reference value. The test fails if the current run is more than 10% slower than the reference.


### Regression Tests

- Tolerance (`t`) test: Uses [probtest](https://github.com/MeteoSwiss/probtest) to generate statistics on the output of the base simulation. These statistics are then compared to a set of stored reference values, and the test checks that all deviations fall within predefined tolerance intervals. Both the reference values and their associated tolerances are stored in advance.

- Update (`u`) test: Checks bit-identity of results of the base simulation with stored reference values.


### Runtime Configuration Tests

- CUDA Graph (`g`) test: Runs the simulation with CUDA Graphs enabled and compares the output to the base simulation to verify correctness. This tests compatibility and stability of the code when GPU graph execution is active.


### Sanitizer Tests

- Compute Sanitizer (`c`) test: Runs the simulation under NVIDIA's [Compute Sanitizer](https://docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html) to detect memory access errors, uninitialized memory usage, synchronization issues, and potential data races on the GPU. This test is useful for debugging low-level GPU issues but can significantly increase runtime and produce very large log files.


### Technical Feature Tests

- Restart (`r`) test: Tests the "restart" feature, by writing a checkpoint file at the middle of the base simulation runtime, and then running a second simulation starting from that checkpoint. Then, the test checks bit-identity of results between the second half of the base simulation and the second simulation.


### Technical Decomposition Tests

- mpi (`m`) test: Runs a simulation with modified MPI settings. Specifically, it reduces the number of MPI processes per node by 1 if the number is greater than 1; otherwise, it reduces the number of nodes by 1. After running the simulation, it checks bit-identity of results between this and the base simulation.

- nproma (`n`) test: Runs a simulation with modified nproma settings (_nproma modified = nproma base + 1_), then checks bit-identity of results between this and the base simulation.

- omp (`o`) test: Runs a simulation with modified  omp settings (one thread less per MPI processor), then checks bit-identity of results between this and the base simulation.
