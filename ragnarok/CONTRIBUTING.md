# ICON Ragnarok Contribution Guide

## Design Decisions

### Goals

C++/Kokkos development in ICON should aim to meet the following goals, listed in order of importance:

1. **Readability** - Favor simple, clear implementations and avoid unnecessary complexity or obscurity. Use explicit include paths whenever this improves clarity. A useful set of practices can be found in the [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines).
2. **Correctness** - Unit tests are mandatory, and including functional tests in a merge request is strongly encouraged. For unit-testing, ragnarok uses the [GoogleTest](https://google.github.io/googletest/) framework, which is either provided by the system (e.g., on Levante) or fetched automatically during the build.
3. **Performance portability** - New code must perform well across different target architectures.

ICON’s GitLab CI infrastructure supports both correctness and performance testing. If you are unsure how to validate your changes, feel free to reach out.

### Building Ragnarok

Ragnarok uses the [CMake](https://cmake.org) build system and is automatically configured by ICON through [`configure.ac`](../configure.ac). To build ICON with Ragnarok, pass the `--enable-ragnarok` option to the appropriate configure wrapper. Ragnarok can also be built as a standalone project to allow isolated testing of physics components. Details on dependencies, project structure, compilation, and testing are provided in the project [README](README.md).

When ICON is configured with `--enable-ragnarok`, the [Kokkos](https://github.com/kokkos/kokkos) submodule is built automatically using the selected backend. For information about using or developing with Kokkos, refer to the official [Kokkos tutorials](https://github.com/kokkos/kokkos-tutorials).

### Adding a New Physics Component to Ragnarok

To add a new physics component, create a directory named (`<component_name>`) in `ragnarok/` and place the component’s source and test files there.
Add the component to the build system by calling `add_subdirectory(<component_name>)` in the project's [CMakeLists.txt](CMakeLists.txt).
For reference, the `aes_thermodynamics` structure looks like this:
```
└── 📁ragnarok
    └── 📁aes_thermodynamics
        ├── CMakeLists.txt
        ├── thermo_constants.hpp
        ├── thermo_test.cpp
        ├── thermo.hpp
        ├── thermo.ipp
```

## Coding Style

The ICON C++/Kokkos codebase follows the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html), with minor ICON-specific adaptations. We use [`pre-commit`](https://pre-commit.com) hooks to ensure that:

* code formatting follows the Google style via [`clang-format` tool](https://clang.llvm.org/docs/ClangFormat.html) and
* common technical issues are identified early through static analysis using [`clang-tidy` tool](https://clang.llvm.org/extra/clang-tidy/).

These tools automatically modify staged files. The same checks run in GitLab CI and are required for any merge request. Instructions for enabling `pre-commit` in ICON are available in the [ICON Contribution Guidelines](../CONTRIBUTING.md#code-contributions).

### General Suggestions

* Ragnarok supports the C++20 standard.
* Code should be as self-explanatory as possible. Public functions and interfaces must include inline documentation using Doxygen syntax `///`. Documentation for private functions is encouraged but optional.
* Place documentation close to the implementation where appropriate— for example, in `.cpp` rather than the `.hpp` file.

> **Note:** Ragnarok aims for a rigorous review process via GitLab merge requests. Ensure that your code follows the principles above, includes appropriate documentation, and is covered by automated tests before marking a merge request as 'Ready'.
