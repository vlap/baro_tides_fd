#!/bin/bash
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MKLROOT="/opt/intel/oneapi/mkl/latest"
export IOMPROOT="/opt/intel/oneapi/compiler/latest/lib"

# Default to first available compiler
if command -v gfortran &> /dev/null; then FC_ACTIVE="gfortran"; else FC_ACTIVE="ifx"; fi

export NETCDFF_PREFIX="${PROJECT_DIR}/vendor/netcdf-${FC_ACTIVE}"
[ -d "${PROJECT_DIR}/vendor/cmake/bin" ] && export PATH="${PROJECT_DIR}/vendor/cmake/bin:$PATH"
export PATH="${NETCDFF_PREFIX}/bin:$PATH"
export LD_LIBRARY_PATH="${NETCDFF_PREFIX}/lib:${IOMPROOT}:$LD_LIBRARY_PATH"
export LIBRARY_PATH="${IOMPROOT}:$LIBRARY_PATH"
export CMAKE_PREFIX_PATH="${NETCDFF_PREFIX}:${CMAKE_PREFIX_PATH}"

echo ">>> Baro Tides FD Environment Activated (${FC_ACTIVE})"
echo ">>> To use Intel, run: export FC=ifx && export NETCDFF_PREFIX=${PROJECT_DIR}/vendor/netcdf-ifx"
echo ">>> Build with: cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build"
