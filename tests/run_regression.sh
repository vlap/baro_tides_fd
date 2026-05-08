#!/bin/bash
# Barotropic Tides Solver - Regression Test Script
# Goal: Ensure numerical consistency during refactoring.

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"

source ./activate.sh

echo ">>> Building..."
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

echo ">>> Running regression test..."
./build/bin/baro_v1 tests/regression/control_file_test.txt > tests/regression/current_run.log 2>&1

echo ">>> Verifying results..."
# Compare domain integrals with tolerance for small numerical drift
# di.txt structure: component, KE, PE, D, D_BL, D_IT, D_SAL, D_f, etc.
# We'll use a simple diff for now, but in the future we could use a Python script for epsilon comparison.
if diff tests/regression/master_di.txt data/LAG/baro_fd/out/0000_00_00__00_00/global/sols/di.txt; then
    echo ">>> REGRESSION TEST PASSED: Numerical results match Golden Master."
else
    echo ">>> REGRESSION TEST FAILED: Numerical drift detected!"
    exit 1
fi
