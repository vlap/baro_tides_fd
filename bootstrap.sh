#!/bin/bash
# ==============================================================================
# Baro Tides FD - Bootstrap Script
# ==============================================================================
# Goal: One-click installation of all dependencies for HPC and Laptops.
# This script uses Conda to create a portable, reproducible environment.
# ==============================================================================

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${PROJECT_DIR}/.env_baro"

echo ">>> Starting Baro Tides FD Bootstrap..."

# 0. HPC Detection
if command -v module &> /dev/null; then
    echo ">>> HPC environment detected (module command available)."
    echo ">>> Note: You can also use 'module load intel-oneapi-mkl netcdf-fortran' if preferred."
fi

# 1. Toolchain Detection
if command -v cmake &> /dev/null; then
    echo ">>> CMake detected: $(cmake --version | head -n 1)"
else
    echo ">>> CMake not found. Installing portable version to ${PROJECT_DIR}/vendor/cmake..."
    CMAKE_VER="3.26.4"
    mkdir -p "${PROJECT_DIR}/vendor/cmake"
    cd "${PROJECT_DIR}/vendor"
    if [ ! -f "cmake-${CMAKE_VER}-linux-x86_64.tar.gz" ]; then
        wget "https://github.com/Kitware/CMake/releases/download/v${CMAKE_VER}/cmake-${CMAKE_VER}-linux-x86_64.tar.gz"
    fi
    tar -xzf "cmake-${CMAKE_VER}-linux-x86_64.tar.gz" --strip-components=1 -C "${PROJECT_DIR}/vendor/cmake"
    export PATH="${PROJECT_DIR}/vendor/cmake/bin:$PATH"
    cd "${PROJECT_DIR}"
fi

if command -v ifx &> /dev/null || command -v gfortran &> /dev/null; then
    echo ">>> System compilers detected. Skipping Conda installation."
else
    echo ">>> No compilers found. Attempting to install Miniconda locally..."
    # ... (Conda installation code)
fi

# 2. Build NetCDF locally (Compiler Specific)
BUILD_DIR="${PROJECT_DIR}/vendor"
mkdir -p "${BUILD_DIR}/src"

# Detect all available compilers
COMPILERS=""
command -v gfortran &> /dev/null && COMPILERS="$COMPILERS gfortran"
command -v ifx &> /dev/null && COMPILERS="$COMPILERS ifx"
command -v ifort &> /dev/null && COMPILERS="$COMPILERS ifort"

for FC_EXEC in $COMPILERS; do
    INSTALL_PREFIX="${BUILD_DIR}/netcdf-${FC_EXEC}"
    
    if [ ! -f "${INSTALL_PREFIX}/lib/libnetcdff.a" ]; then
        echo ">>> Building NetCDF for ${FC_EXEC}..."
        
        # Set C compiler to match Fortran
        if [[ "$FC_EXEC" == "ifx" ]]; then CC_EXEC="icx"
        elif [[ "$FC_EXEC" == "ifort" ]]; then CC_EXEC="icc"
        else CC_EXEC="gcc"; fi

        # Build C library first
        cd "${BUILD_DIR}/src"
        NC_VER="4.8.1"
        [ ! -f "netcdf-c-${NC_VER}.tar.gz" ] && wget --no-check-certificate "https://github.com/Unidata/netcdf-c/archive/refs/tags/v${NC_VER}.tar.gz" -O "netcdf-c-${NC_VER}.tar.gz"
        tar -xzf "netcdf-c-${NC_VER}.tar.gz"
        cd "netcdf-c-${NC_VER}"
        export CC=$CC_EXEC
        ./configure --prefix="${INSTALL_PREFIX}" --disable-shared --disable-dap --disable-hdf5 --disable-libxml2
        make -j$(nproc) && make install
        
        # Build Fortran library
        cd "${BUILD_DIR}/src"
        NF_VER="4.5.3"
        [ ! -f "netcdf-fortran-${NF_VER}.tar.gz" ] && wget --no-check-certificate "https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${NF_VER}.tar.gz" -O "netcdf-fortran-${NF_VER}.tar.gz"
        tar -xzf "netcdf-fortran-${NF_VER}.tar.gz"
        cd "netcdf-fortran-${NF_VER}"
        export FC=$FC_EXEC
        export CC=$CC_EXEC
        export CPPFLAGS="-I${INSTALL_PREFIX}/include"
        export LDFLAGS="-L${INSTALL_PREFIX}/lib"
        ./configure --prefix="${INSTALL_PREFIX}" --disable-shared
        make -j1 && make install
        
        echo ">>> NetCDF stack installed to ${INSTALL_PREFIX}"
    else
        echo ">>> NetCDF for ${FC_EXEC} already exists."
    fi
done
cd "${PROJECT_DIR}"

# 3. Finalize
echo ">>> Environment ready!"

# 4. Generate a convenience activation script
cat > "${PROJECT_DIR}/activate.sh" <<EOF
#!/bin/bash
PROJECT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
export MKLROOT="/opt/intel/oneapi/mkl/latest"
export IOMPROOT="/opt/intel/oneapi/compiler/latest/lib"

# Default to first available compiler
if command -v gfortran &> /dev/null; then FC_ACTIVE="gfortran"; else FC_ACTIVE="ifx"; fi

export NETCDFF_PREFIX="\${PROJECT_DIR}/vendor/netcdf-\${FC_ACTIVE}"
[ -d "\${PROJECT_DIR}/vendor/cmake/bin" ] && export PATH="\${PROJECT_DIR}/vendor/cmake/bin:\$PATH"
export PATH="\${NETCDFF_PREFIX}/bin:\$PATH"
export LD_LIBRARY_PATH="\${NETCDFF_PREFIX}/lib:\${IOMPROOT}:\$LD_LIBRARY_PATH"
export LIBRARY_PATH="\${IOMPROOT}:\$LIBRARY_PATH"
export CMAKE_PREFIX_PATH="\${NETCDFF_PREFIX}:\${CMAKE_PREFIX_PATH}"

echo ">>> Baro Tides FD Environment Activated (\${FC_ACTIVE})"
echo ">>> To use Intel, run: export FC=ifx && export NETCDFF_PREFIX=\${PROJECT_DIR}/vendor/netcdf-ifx"
echo ">>> Build with: cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build"
EOF

chmod +x "${PROJECT_DIR}/activate.sh"
echo ">>> Done. Use './activate.sh' to start."
