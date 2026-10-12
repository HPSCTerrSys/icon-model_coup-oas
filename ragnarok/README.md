# ICON RAGNAROK

Ragnarok is the C++ implementation of the AES physics based on the [Kokkos](https://github.com/kokkos/kokkos) programming model.


## Project Structure

```
└── 📁ragnarok
    └── 📁cmake
        └── 📁caches
            ├── levante_cpu.cmake
            ├── levante_gpu.cmake
    └── 📁common
        ├── CMakeLists.txt
        ├── types_test.cpp
        ├── types.hpp
    └── 📁support
        └── 📁icon_bridge
        └── 📁testing
        └── 📁timer
        ├── CMakeLists.txt
    └── 📁aes_component
        ├── CMakeLists.txt
    ├── CMakeLists.txt
    ├── CONTRIBUTING.md
    └── README.md
```

Shared functionality that is used across physics components belongs in the `common/` directory.

The `support/` directory provides infrastructure used throughout Ragnarok:

* `icon_bridge` provides access to the data and functionality available in the main ICON code. Additionally, it also provides initialize and finalize routines for Ragnarok. Its support will be extended as needed.
* `timer` mirrors the ICON timer infrastructure. When Ragnarok is built as part of ICON, calls to `start` and `stop` are forwarded to `mo_timer.f90` so that Ragnarok timers appear in the experiment logfile. When built standalone, Ragnarok provides an ICON-compatible timer implementation based on `Kokkos::Timer`.
* `testing` provides the entry point for standalone unit testing and is also used in the GitLab CI pipelines.

## Usage

Ragnarok is integrated into ICON and is built automatically when ICON is configured with `--enable-ragnarok`. The ICON build system also ensures that the correct compilation flags and the appropriate Kokkos backend (from the `externals/kokkos` submodule) are used.

## Available Compile Options

> **Note:** default options are marked in __bold__

* _Unit-test_ - compile unit-tests
  * BUILD_TESTING=__ON__/OFF
* Standalone - compile standalone components
  * RGK_ENABLE_STANDALONE=__ON__/OFF
* Bindings - interoperability with other programming languages; Fortran is supported by default.
  * RGK_ENABLE_PYTHON_BINDINGS=ON/**OFF**
* Precision - either single or double precision
  * RGK_ENABLE_SINGLE_PRECISION=ON/**OFF**

### Build as part of ICON

The ICON configuration scripts automatically configure both Ragnarok and Kokkos.

### Ragnarok standalone

The standalone has a dependency of [NetCDF for CXX](https://github.com/Unidata/netcdf-cxx4). For Levante: `spack load netcdf-cxx4@4.3.1`.

To build ragnarok standalone, one can use the cache files provided for levante under `ragnarok/cmake/caches/`, e.g.:
  ```
  cmake -B build \
        -S ragnarok \
        -C ragnarok/cmake/caches/levante_cpu.cmake \
        && cmake --build build --parallel
  ```

The aes_microphysics executable can be run in single or double precision, e.g.:
  ```
    ./build/aes_microphysics/standalone/aes_graupel <input.nc> double  
  ```
The results are saved to NetCDF file called `output.nc`.

### Interoperability with python

The aes_microphysics provides python bindings for connecting the library to python environments. To ensure that, python3 (&pybind11) and netcdf-cxx4 to be accessible. For Levante:
  ```
  module load python3
  python3 -m venv venv
  source venv/bin/activate
  pip install -r ragnarok/aes_microphysics/bindings/python/requirements.txt
  cmake -B build_py \
        -S ragnarok \
        -C ragnarok/cmake/caches/levante_cpu.cmake \
        -DRGK_ENABLE_PYTHON_BINDINGS=ON \
        -Dpybind11_DIR=<path-to-venv>/venv/lib/python3.10/site-packages/pybind11/share/cmake/pybind11/
  cmake --build build_py --parallel
  ```

The following shared library will be generated `build_py/aes_microphysics/aes_muphys_py.so`, which can be later imported by a python source code.

## Testing

Ragnarok uses GoogleTest for unit testing. To run the tests:

  ```
  cd build && ctest
  ```
