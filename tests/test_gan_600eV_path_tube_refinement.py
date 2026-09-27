"""Geometry-only guard for the same-600-eV GaN tube refinement."""

from __future__ import annotations

import json
from pathlib import Path

from scripts.prepare_gan_600eV_path_tube_refinement import REFINEMENT_ANCHORS
from scripts.screen_gan_600eV_path_tube_refinement import screen
from scripts.analyze_gan_600eV_tube_residual import analyze


ROOT = Path(__file__).resolve().parents[1]
TRAJECTORY = ROOT / "paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927/gan_vasp_45p7_final_chain.traj"
HESSIAN = ROOT / "benchmarks/numerical_integrity/gan_600eV_ts_hessian_0p02_20260928/joint_hessian.npz"
FEASIBILITY = ROOT / "benchmarks/numerical_integrity/gan_600eV_path_tube_feasibility_20260928.json"
SCREEN = ROOT / "benchmarks/numerical_integrity/gan_600eV_path_tube_refinement_screen_20260928.json"
FIRST = ROOT / "benchmarks/numerical_integrity/gan_600eV_path_tube_20260928.json"
REFINED = ROOT / "benchmarks/numerical_integrity/gan_600eV_path_tube_refinement_20260928.json"
TRANSPORT = ROOT / "benchmarks/numerical_integrity/gan_600eV_path_tube_transport_comparison_20260928.json"


def test_all_refinement_anchors_pass_both_signed_geometry_checks() -> None:
    report = screen(TRAJECTORY, HESSIAN, FEASIBILITY)
    pinned = json.loads(SCREEN.read_text(encoding="utf-8"))
    assert report["status"] == pinned["status"]
    assert report["rows"] == pinned["rows"]
    safe = {row["image_index"] for row in report["rows"]
            if row["amplitude_A"] == 0.015 and row["both_sides_pass_empirical_screen"]}
    assert set(REFINEMENT_ANCHORS) <= safe
    assert 3 not in safe


def test_refined_tube_fails_contour_gate_for_localized_physical_response() -> None:
    report = analyze(FIRST, REFINED, TRANSPORT)
    assert 0.10 < report["refined_max_LOO_meV_per_GaN"] < report["first_max_LOO_meV_per_GaN"]
    assert (5, 6) in {(turn["left_image"], turn["right_image"])
                      for turn in report["analytic_sign_turns"]}
    assert (19, 20) in {(turn["left_image"], turn["right_image"])
                       for turn in report["analytic_sign_turns"]}
    near_peak = next(row for row in report["anchor_diagnostics"] if row["image_index"] == 18)
    assert near_peak["quadratic_minimum_outside_measured_width"]
