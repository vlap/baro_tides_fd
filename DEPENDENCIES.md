# Baro Tides FD - Dependency Tracking

## Build Tools
- **GNU Make** (v3.81 or later)
- **Compilers:**
  - GNU Fortran (`gfortran`) v9.4+
  - Intel Fortran (`ifx` or `ifort`) v2021.1+

## Core Libraries
- **Intel MKL (Math Kernel Library):**
  - Used for PARDISO direct solver.
  - Used for BLAS/LAPACK and Sparse BLAS operations.
- **NetCDF (Network Common Data Form):**
  - `netcdf-c` (v4.7+)
  - `netcdf-fortran` (v4.5+) - *Must be built with the same compiler used for the project.*

## Optional Libraries
- **SuiteSparse (UMFPACK):** Used as an alternative sparse solver.

## Data Requirements
- **ETOPO:** Topography data in NetCDF format.
- **WOA05:** World Ocean Atlas 2005 stratification data in NetCDF format.

## System Requirements
- **OS:** Linux (Ubuntu 20.04+ recommended, standard HPC distributions like RHEL/CentOS supported).
- **Environment Management:** Conda/Miniconda (Recommended for portability).
