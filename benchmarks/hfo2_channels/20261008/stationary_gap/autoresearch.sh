#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
output=${1:?provide a new output JSON path}
python -m scripts.check_stationary_gap --output "$output"
