"""Regression guards for the QE static-audit production gate."""

from __future__ import annotations

import runpy
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "validate_qe_static_convergence_gate.py"


def _audit() -> dict:
    return {
        "status": "ok",
        "n_static_points": 3,
        "fixed_conditions": {
            "git_revision": "3016099",
            "initial_endpoint_sha256": "qe-initial",
            "calculator_parameters_except_cutoffs": {
                "kpts": [4, 4, 4],
                "input_data": {"electrons": {"conv_thr": 1e-8}},
            },
        },
        "accepted_final_pair": {"higher_cutoffs_Ry": [100.0, 600.0]},
    }


def _endpoint_gate() -> dict:
    return {"matches": True, "endpoints": {"initial": {"matches": True, "candidate_sha256": "qe-initial"}}}


def test_qe_static_gate_accepts_matching_static_basis() -> None:
    validate = runpy.run_path(str(SCRIPT))["validate"]
    validate(_audit(), _endpoint_gate(), ecutwfc=100, ecutrho=600, kpts=(4, 4, 4), scf_thr=1e-8, git_revision="3016099")


def test_qe_static_gate_rejects_unknown_or_changed_basis() -> None:
    validate = runpy.run_path(str(SCRIPT))["validate"]
    with pytest.raises(ValueError, match="explicit VCNEB_GIT_REVISION"):
        validate(_audit(), _endpoint_gate(), ecutwfc=100, ecutrho=600, kpts=(4, 4, 4), scf_thr=1e-8, git_revision="remote-sync-unknown")
    with pytest.raises(ValueError, match="highest cutoff"):
        validate(_audit(), _endpoint_gate(), ecutwfc=120, ecutrho=720, kpts=(4, 4, 4), scf_thr=1e-8, git_revision="3016099")
