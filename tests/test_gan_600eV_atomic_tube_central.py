"""Guard the predeclared central GaN atomic-tube coordinate choice."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.compare_gan_600eV_atomic_tube_transports import analyze
from scripts.prepare_gan_600eV_atomic_tube_central import AMPLITUDE_A, ANCHORS


ROOT = Path(__file__).resolve().parents[1]
TRAJECTORY = ROOT / "paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927/gan_vasp_45p7_final_chain.traj"
ATOMIC_FEASIBILITY = ROOT / "benchmarks/numerical_integrity/gan_600eV_atomic_tube_feasibility_20260928.json"
HESSIAN = ROOT / "benchmarks/numerical_integrity/gan_600eV_ts_hessian_0p02_20260928/joint_hessian.npz"
TRANSPORT = ROOT / "benchmarks/numerical_integrity/gan_600eV_atomic_tube_transport_comparison_20260928.json"


def test_central_chart_is_smooth_but_global_mode_identity_is_not() -> None:
    calculated = analyze(TRAJECTORY, ATOMIC_FEASIBILITY, HESSIAN)
    pinned = json.loads(TRANSPORT.read_text(encoding="utf-8"))
    assert calculated["status"] == pinned["status"]
    assert calculated["fixed_seed_projection"]["minimum_central_5_to_22_overlap"] == pytest.approx(
        pinned["fixed_seed_projection"]["minimum_central_5_to_22_overlap"]
    )
    assert calculated["fixed_seed_projection"]["minimum_central_5_to_22_overlap"] > 0.95
    assert calculated["minimum_fixed_vs_recursive_normal_overlap"] < 0
    assert ANCHORS == (5, 8, 11, 14, 17, 18, 20, 22)
    assert AMPLITUDE_A == 0.05
