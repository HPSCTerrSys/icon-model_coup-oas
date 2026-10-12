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

# system paths on levante
set(GTest_ROOT
    "/sw/spack-levante/googletest-1.10.0-opzgcq"
    CACHE PATH "Path to Googletest on Levante" FORCE
)
set(NetCDF_ROOT
    "/sw/spack-levante/netcdf-cxx4-4.3.1-42ju4n;/sw/spack-levante/netcdf-c-main-k4lh4v"
    CACHE PATH "Path to NetCDF C++ on Levante" FORCE
)

set(GCC_ROOT "/sw/spack-levante/gcc-11.2.0-bcn7mb")
set(NVHPC_ROOT "/sw/spack-levante/nvhpc-24.7-py26uc/Linux_x86_64/24.7")

# cmake options
set(CMAKE_CXX_COMPILER
    "${NVHPC_ROOT}/compilers/bin/nvc++"
    CACHE STRING "C++ compiler" FORCE
)
set(CMAKE_C_COMPILER
    "${NVHPC_ROOT}/compilers/bin/nvc"
    CACHE STRING "C compiler" FORCE
)
set(CMAKE_Fortran_COMPILER
    "${NVHPC_ROOT}/compilers/bin/nvfortran"
    CACHE STRING "Fortran compiler" FORCE
)
set(CMAKE_CUDA_COMPILER
    "${NVHPC_ROOT}/compilers/bin/nvcc"
    CACHE STRING "CUDA compiler" FORCE
)
set(CMAKE_CUDA_HOST_COMPILER
    "${GCC_ROOT}/bin/g++"
    CACHE STRING "CUDA host compiler" FORCE
)
set(CMAKE_CUDA_ARCHITECTURES
    80
    CACHE STRING "Cuda Arch" FORCE
)
set(CMAKE_BUILD_TYPE
    "RelWithDebInfo"
    CACHE STRING "Build type" FORCE
)

# A workaround for the compiler bug because NVHPC is supposed to add RPATHs to
# the gcc libraries, but doesn't do it
set(CMAKE_BUILD_RPATH
    "${CMAKE_BUILD_RPATH};${GCC_ROOT}/lib64"
    CACHE INTERNAL "Extending -rpath"
)

# ragnarok options
set(BUILD_TESTING
    ON
    CACHE BOOL "Enable testing" FORCE
)
set(RGK_ENABLE_STANDALONE
    ON
    CACHE BOOL "Enable standalone" FORCE
)
set(RGK_ENABLE_PYTHON_BINDINGS
    OFF
    CACHE BOOL "Enable python bindings" FORCE
)

# kokkos options
set(Kokkos_ENABLE_SERIAL
    ON
    CACHE BOOL "Build Kokkos Serial backend" FORCE
)
set(Kokkos_ENABLE_CUDA
    ON
    CACHE BOOL "Build Kokkos CUDA backend" FORCE
)
set(Kokkos_ARCH_AMPERE80
    ON
    CACHE BOOL "Build Kokkos for AMPERE arch" FORCE
)
