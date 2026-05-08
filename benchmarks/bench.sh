#!/bin/bash
# Barotropic Tides Solver - Benchmarking Script
# Goal: Record performance metrics to track optimization progress.

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"

source ./activate.sh

echo ">>> Building in Release mode..."
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

echo ">>> Running benchmark (Medium Case: 180x90)..."
START_TIME=$(date +%s.%N)
./build/bin/baro_v1 tests/regression/control_file_test.txt > benchmarks/last_run.log 2>&1
END_TIME=$(date +%s.%N)

# Calculate elapsed time
ELAPSED=$(python3 -c "print(round($END_TIME - $START_TIME, 3))")

echo ">>> BENCHMARK RESULT: ${ELAPSED} seconds"

# Append to log
echo "$(date '+%Y-%m-%d %H:%M:%S') | $(gfortran --version | head -n 1) | ${ELAPSED}s" >> benchmarks/performance_history.log
