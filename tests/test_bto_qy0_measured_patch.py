"""Offline consistency checks for the preregistered restricted BTO patch."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from examples.bto_q1q2_reference import sha256


ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "benchmarks" / "numerical_integrity"


def _load(name: str) -> dict:
    return json.loads((BENCH / name).read_text(encoding="utf-8"))


def test_archived_qy0_reuse_source_and_auditor_are_frozen() -> None:
    report = _load("bto_qy0_reuse_with_holdout_audit_20260928.json")
    assert report["n_reused_DFT_points"] == 5
    assert report["n_new_DFT_evaluations"] == 0
    assert report["source_manifest_sha256"] == sha256(
        BENCH / "bto_qy0_reuse_sources_with_holdout_2026-09-28.json"
    )
    assert report["auditor_sha256"] == sha256(ROOT / "scripts" / "audit_bto_qy0_reuse.py")
    assert [row["role"] for row in report["points"]].count(
        "independent_cell_center_holdout"
    ) == 1


def test_ten_measured_nodes_and_holdout_predictions_are_internally_consistent() -> None:
    measured_path = BENCH / "bto_qy0_ten_measured_nodes_20260928.json"
    measured = json.loads(measured_path.read_text(encoding="utf-8"))
    plan = _load("bto_qy0_three_holdout_plan_20260928.json")
    assert measured["n_measured_nodes"] == 10
    assert measured["n_reused_archived_nodes"] == 5
    assert measured["n_new_pilot_nodes"] == 5
    assert measured["auditor_sha256"] == sha256(
        ROOT / "scripts" / "assemble_bto_qy0_measured_patch.py"
    )
    assert plan["source_measured_patch_sha256"] == sha256(measured_path)
    assert plan["n_holdouts"] == 3
    assert plan["energy_absolute_error_gate_meV_per_BTO"] == 2.0
    assert plan["full_atom_plus_strain_metric_coordinate_error_gate_sqrt_amu_A"] == 0.10
    nodes = {tuple(row["q"]): row for row in measured["points"]}
    assert len(nodes) == 10
    assert all(abs(row["third_soft_y_amplitude_sqrt_amu_A"]) < 1e-8
               and row["restricted_gradient_eV_per_sqrt_amu_A"] <= 0.003
               and row["maximum_absolute_stress_kbar"] <= 2.0
               and len(row["coordinates_u_A_eta_voigt"]) == 21
               for row in nodes.values())
    for holdout in plan["holdouts"]:
        q = tuple(holdout["q_parallel_q_transverse_sqrt_amu_A"])
        corners = [nodes[tuple(corner)] for corner in holdout["source_corners_q"]]
        weights = np.asarray(holdout["source_corner_weights"], dtype=float)
        assert q not in nodes
        np.testing.assert_allclose(
            weights @ np.asarray(holdout["source_corners_q"]), q,
            atol=1e-12, rtol=0,
        )
        np.testing.assert_allclose(
            weights @ np.asarray([row["energy_minus_c_eV_per_BTO"] for row in corners]),
            holdout["predicted_energy_minus_c_eV_per_BTO"], atol=1e-12, rtol=0,
        )
        np.testing.assert_allclose(
            weights @ np.asarray([row["coordinates_u_A_eta_voigt"] for row in corners]),
            holdout["predicted_coordinates_u_A_eta_voigt"], atol=1e-12, rtol=0,
        )
    center = nodes[(0.75, 0.15)]["energy_minus_c_eV_per_BTO"]
    corners = [nodes[q]["energy_minus_c_eV_per_BTO"]
               for q in ((0.6, 0.0), (0.6, 0.3), (0.9, 0.0), (0.9, 0.3))]
    assert abs(1000.0 * (center - sum(corners) / 4.0) - 0.8112019398822667) < 1e-7
