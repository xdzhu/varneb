#!/bin/bash
# Local analytic/parser checks only; no implicit DFT submission.
set -euo pipefail
python -m pytest -q tests/test_relaxed_curvature.py tests/test_joint_curvature.py tests/test_hfo2_static_replica_audit.py
