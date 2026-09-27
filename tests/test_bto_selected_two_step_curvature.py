"""Regression gate for the raw-audited BTO local holdout Hessians."""

from __future__ import annotations

import json
from pathlib import Path
from shutil import copy2

import pytest

from scripts.audit_bto_selected_two_step_curvature import join_two_steps


ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "paper/VARNEB_CPC/evidence/bto_q075_q015_holdout_2026-09-27"


def test_two_step_soft_screen_does_not_overclaim_a_minimum() -> None:
    report = join_two_steps(FOLDER)
    frozen = json.loads((FOLDER / "two_step_curvature_screen.json").read_text(encoding="utf-8"))
    assert report == frozen
    assert report["status"] == "positive_force_hessian_sign_reproduced_energy_force_certification_open"
    assert report["absolute_lowest_eigenvector_overlap"] > 0.9999
    assert report["minimum_eigenvalue_step_difference_eV_per_amu_A2"] < 3e-5
    assert report["energy_hessian_diagonal_max_abs_differences_eV_per_amu_A2"][0] > min(
        report["minimum_eigenvalues_eV_per_amu_A2"]
    )
    assert report["n_raw_audited_DFT_points_by_step"] == [109, 141]


def test_two_step_screen_rejects_swapped_source(tmp_path: Path) -> None:
    for label in ("0p05", "0p10"):
        for kind in ("result", "raw_audit", "reconstruction"):
            name = f"curvature_{label}_{kind}.json"
            copy2(FOLDER / name, tmp_path / name)
    small_result = tmp_path / "curvature_0p05_result.json"
    content = json.loads(small_result.read_text(encoding="utf-8"))
    content["minimum_eigenvalue_eV_per_amu_A2"] = -1.0
    small_result.write_text(json.dumps(content), encoding="utf-8")
    with pytest.raises(ValueError, match="inconsistent raw-audited"):
        join_two_steps(tmp_path)
