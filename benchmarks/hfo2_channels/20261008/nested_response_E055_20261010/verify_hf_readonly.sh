#!/bin/bash
set -euo pipefail
analysis_root=/public/home/iai806/abacus/agent-runs/20261010-varneb-nested-controls-E055-r1
archive=/public/home/iai806/abacus/agent-runs/varneb-E055-clean-86d2637.tar
[[ ! -e "$analysis_root" ]] || exit 2
[[ $(sha256sum "$archive" | cut -d' ' -f1) == cbe6552a3fa4237fec2098c08c8056735a9a95df4e18691d09fed6b4f24a613c ]] || exit 2
mkdir -p "$analysis_root/source"
tar -xf "$archive" -C "$analysis_root/source"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$analysis_root/source"
cd "$analysis_root/source"
py=/public/home/iai806/.conda/envs/icu/bin/python
"$py" -m scripts.check_nested_response --output "$analysis_root/analytic_hf.json"
"$py" -m pytest -q tests/test_nested_response.py tests/test_response_cost.py \
 tests/test_nested_response_benchmark.py tests/test_stationary_gap.py \
 tests/test_stationary_branch.py tests/test_prediction_controls.py \
 --junitxml="$analysis_root/hf_readonly_junit.xml"
"$py" -m scripts.export_hfo2_scf_costs \
 --run-root /public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E053-r1 \
 --output "$analysis_root/prior_pilot_scf_costs.json"
# No sbatch, mpirun, ABACUS input writes, live-code replacement or held-out labels.
