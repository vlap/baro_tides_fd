cmake_minimum_required(VERSION 3.18)
project(baro_tides_fd LANGUAGES Fortran C)

# --- Options ---
option(USE_OPENMP "Enable OpenMP" ON)

# --- Compiler Flags ---
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ffast-math -funroll-loops")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbounds-check -fbacktrace -Wall")
    add_compile_options("-mcmodel=medium")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -ipo -fp-model fast=2")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback -check all -fpe0")
    add_compile_options("-mcmodel=medium" "-shared-intel")
endif()

# --- Dependencies ---

# OpenMP
if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

# NetCDF
# NetCDF Fortran usually depends on NetCDF C. 
# We'll use pkg-config or standard search.
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(NETCDF_F netcdf-fortran)
endif()

if(NOT NETCDF_F_FOUND)
    # Fallback to manual search
    find_path(NETCDF_F_INCLUDE_DIRS netcdf.mod HINTS ${NETCDFF_PREFIX}/include)
    find_library(NETCDF_F_LIBRARIES netcdff HINTS ${NETCDFF_PREFIX}/lib)
    if(NETCDF_F_INCLUDE_DIRS AND NETCDF_F_LIBRARIES)
        set(NETCDF_F_FOUND TRUE)
    endif()
endif()

if(NOT NETCDF_F_FOUND)
    message(FATAL_ERROR "NetCDF-Fortran not found. Please set NETCDFF_PREFIX or ensure it's in your environment.")
endif()

# MKL
if(DEFINED ENV{MKLROOT})
    set(MKL_ROOT $ENV{MKLROOT})
elseif(DEFINED MKLROOT)
    set(MKL_ROOT ${MKLROOT})
endif()

if(MKL_ROOT)
    set(MKL_INCLUDE_DIRS "${MKL_ROOT}/include")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        set(MKL_IFACE "gf")
        set(MKL_THREAD "gnu_thread")
    else()
        set(MKL_IFACE "intel")
        set(MKL_THREAD "intel_thread")
    endif()

    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(MKL_LIB_DIR "${MKL_ROOT}/lib/intel64")
        set(MKL_LIBRARIES 
            -Wl,--start-group
            "${MKL_LIB_DIR}/libmkl_${MKL_IFACE}_lp64.a"
            "${MKL_LIB_DIR}/libmkl_core.a"
            "${MKL_LIB_DIR}/libmkl_${MKL_THREAD}.a"
            -Wl,--end-group
        )
    else()
        set(MKL_LIB_DIR "${MKL_ROOT}/lib/ia32")
        set(MKL_LIBRARIES 
            -Wl,--start-group
            "${MKL_LIB_DIR}/libmkl_${MKL_IFACE}.a"
            "${MKL_LIB_DIR}/libmkl_core.a"
            "${MKL_LIB_DIR}/libmkl_${MKL_THREAD}.a"
            -Wl,--end-group
        )
    endif()
    set(MKL_FOUND TRUE)
else()
    find_package(MKL QUIET)
endif()

message(STATUS "MKL_ROOT: ${MKL_ROOT}")
message(STATUS "MKL_INCLUDE_DIRS: ${MKL_INCLUDE_DIRS}")
message(STATUS "MKL_LIBRARIES: ${MKL_LIBRARIES}")

if(NOT MKL_FOUND)
    message(FATAL_ERROR "Intel MKL not found. Please set MKLROOT.")
endif()

# --- Executable ---
add_executable(baro_v1)

target_sources(baro_v1 PRIVATE
    # Low level / Misc
    src/defs/precisions.f90
    src/defs/f90_kind.f90
    src/misc/config_mod.f90
    src/misc/iso_varying_string.f90
    src/misc/dispmodule.f90
    src/misc/my_trigs.f90
    src/misc/err_manager.f90
    src/misc/control_file.f90
    src/misc/save_load.f90
    src/misc/interpolate.f90
    src/misc/spline.f90
    
    # Sparse / Solvers
    src/sparse/blassm.f
    src/sparse/formats.f
    src/sparse/unary.f
    src/sparse/ilut.f
    src/sparse/mkl_spblas.f90
    src/sparse/mkl_legacy_wrapper.f90
    src/sparse/my_sparse.f90
    src/sparse/my_sparse_aggregate.f90
    src/sparse/iter_solvers.f90
    src/sparse/unsym_solvers.f90
    
    # Spherepack
    src/spherepack/sphcom.f
    src/spherepack/hrfft.f
    src/spherepack/shaes.f
    src/spherepack/shaec.f
    src/spherepack/shses.f
    src/spherepack/shsec.f
    src/spherepack/shags.f
    src/spherepack/shsgs.f
    src/spherepack/gaqd.f
    src/shtools/SHRotate.f90
    src/spherepack/spherepack_iface.f90
    src/spherepack/testshrotate.f90
    
    # Generate
    src/generate/read_etopo.f90
    src/generate/generate_global_grid.f90
    src/generate/generate_matrices.f90
    
    # Baro Utils / Main
    src/baro_utils/itd.f90
    src/baro_utils/sal.f90
    src/baro_utils/projection.f90
    src/baro_utils/baro_domain_integrals.f90
    src/baro_utils/baro_solver.f90
    src/baro_utils/baro_v1.f90
)

# --- Includes and Libraries ---
target_include_directories(baro_v1 PRIVATE 
    src/defs
    src/misc
    src/sparse
    src/spherepack
    src/generate
    src/baro_utils
    ${NETCDF_F_INCLUDE_DIRS}
    ${MKL_INCLUDE_DIRS}
)

target_link_libraries(baro_v1 PRIVATE
    ${NETCDF_F_LIBRARIES}
    ${MKL_LIBRARIES}
    $<$<BOOL:${USE_OPENMP}>:OpenMP::OpenMP_Fortran>
    m
    dl
    pthread
)

# Output directory
set_target_properties(baro_v1 PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod"
)

# Special handling for spherepack (needs -fdefault-real-8 in original makefile)
set_source_files_properties(
    src/shtools/SHRotate.f90
    src/spherepack/spherepack_iface.f90
    src/spherepack/testshrotate.f90
    src/spherepack/sphcom.f
    src/spherepack/hrfft.f
    src/spherepack/shaes.f
    src/spherepack/shaec.f
    src/spherepack/shses.f
    src/spherepack/shsec.f
    src/spherepack/shags.f
    src/spherepack/shsgs.f
    src/spherepack/gaqd.f
    PROPERTIES COMPILE_FLAGS "-fdefault-real-8"
)
# Note: Intel's equivalent is -r8. Let's make it compiler-aware if needed.
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set_source_files_properties(
        src/shtools/SHRotate.f90
        src/spherepack/spherepack_iface.f90
        src/spherepack/testshrotate.f90
        src/spherepack/sphcom.f
        src/spherepack/hrfft.f
        src/spherepack/shaes.f
        src/spherepack/shaec.f
        src/spherepack/shses.f
        src/spherepack/shsec.f
        src/spherepack/shags.f
        src/spherepack/shsgs.f
        src/spherepack/gaqd.f
        PROPERTIES COMPILE_FLAGS "-r8"
    )
endif()
