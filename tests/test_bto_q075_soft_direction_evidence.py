"""Regression gate for the audited BTO holdout soft-direction cross-check.

The five-point line does not certify the full conditional surface.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest


EVIDENCE = (Path(__file__).resolve().parents[1] / "paper" / "VARNEB_CPC"
            / "evidence" / "bto_q075_q015_holdout_2026-09-27")
REPO = Path(__file__).resolve().parents[1]


def _read(name: str) -> dict:
    return json.loads((EVIDENCE / name).read_text(encoding="utf-8"))


def _sha(name: str) -> str:
    return hashlib.sha256((EVIDENCE / name).read_bytes()).hexdigest()


def test_selected_soft_direction_has_raw_provenance_and_two_positive_steps() -> None:
    preflight = _read("soft_direction_preflight.json")
    summary = _read("soft_direction_result.json")
    all_points = _read("all_points_145_raw_audit.json")
    audit = _read("soft_direction_raw_audit.json")
    work = _read("soft_direction_work_integral.json")

    assert preflight["status"] == "four_soft_direction_geometries_passed_no_dft"
    assert summary["status"] == "four_static_probes_complete_pending_independent_SCF_audit"
    assert audit["status"] == "four_original_static_DFT_points_verified_soft_direction_curvature_screen_only"
    assert all_points["n_individually_audited_DFT_points"] == 145
    assert len({point["directory"] for point in all_points["evaluations"]}) == 145
    assert audit["n_individually_audited_DFT_points"] == 4
    assert summary["slurm_job_id"] == audit["slurm_job_id"] == "27791225"
    assert audit["source_sha256"]["summary"] == _sha("soft_direction_result.json")
    assert audit["source_sha256"]["soft_preflight"] == _sha("soft_direction_preflight.json")
    assert audit["source_sha256"]["all_points_audit"] == _sha("all_points_145_raw_audit.json")
    assert work["status"] == "work_residual_quantified_cause_unresolved"
    for key, filename in {
        "soft_preflight": "soft_direction_preflight.json",
        "soft_summary": "soft_direction_result.json",
        "soft_audit": "soft_direction_raw_audit.json",
        "branch_result": "conditional_q075_q015_result.json",
        "center_result": "soft_direction_center_result.json",
    }.items():
        assert work["source_sha256"][key] == _sha(filename)

    assert sorted((point["step_sqrt_amu_A"], point["sign"])
                  for point in audit["probes"]) == [
                      (0.05, -1), (0.05, 1), (0.1, -1), (0.1, 1)]
    assert [point["step_sqrt_amu_A"] for point in audit["curvatures"]] == [0.05, 0.1]
    for point in audit["curvatures"]:
        assert point["energy_curvature_eV_per_amu_A2"] > 0
        assert point["gradient_curvature_eV_per_amu_A2"] > 0
        assert point["energy_gradient_disagreement_eV_per_amu_A2"] < 0.0004
    assert work["maximum_absolute_work_residual_meV"] == pytest.approx(
        max(abs(item["E_minus_integrated_gradient_meV"]) for item in work["intervals"])
    )
    assert work["maximum_absolute_work_residual_meV"] < 0.01


def test_line_integral_rebuilds_from_frozen_audited_records(tmp_path: Path) -> None:
    output = tmp_path / "line_integral.json"
    subprocess.run([
        sys.executable, str(REPO / "scripts" / "audit_bto_q1q2_soft_direction_work.py"),
        "--soft-preflight", str(EVIDENCE / "soft_direction_preflight.json"),
        "--soft-summary", str(EVIDENCE / "soft_direction_result.json"),
        "--soft-audit", str(EVIDENCE / "soft_direction_raw_audit.json"),
        "--branch-result", str(EVIDENCE / "conditional_q075_q015_result.json"),
        "--center-result", str(EVIDENCE / "soft_direction_center_result.json"),
        "--output", str(output),
    ], check=True, cwd=REPO, capture_output=True, text=True)
    assert json.loads(output.read_text(encoding="utf-8")) == _read(
        "soft_direction_work_integral.json")
