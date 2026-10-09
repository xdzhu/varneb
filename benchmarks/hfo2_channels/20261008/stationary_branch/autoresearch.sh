#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then
    printf 'Usage: bash autoresearch.sh /new/analytic-check.json\n' >&2
    exit 2
fi
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -m scripts.check_stationary_branch --output "$1"
