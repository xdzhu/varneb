#!/bin/bash
# Local analytic/parser checks only; no implicit DFT submission.
set -euo pipefail
python -m pytest -q tests/test_relaxed_curvature.py tests/test_joint_curvature.py tests/test_hfo2_static_replica_audit.py \
  tests/test_hfo2_atomic_reduction.py tests/test_hfo2_reduction_holdout.py \
  tests/test_hfo2_fixed_input_factory.py tests/test_hfo2_M_endpoint_audit.py tests/test_hfo2_PO_M_chain.py
