"""Offline replay of the audited BTO conditional five-point snapshot."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.audit_bto_conditional_five_point_patch import audit_patch


SOURCE = (Path(__file__).resolve().parents[1] / "benchmarks" / "numerical_integrity"
          / "bto_conditional_five_point_patch_2026-09-28.json")


def _modified_copy(tmp_path: Path, edit) -> Path:
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    edit(data)
    destination = tmp_path / "changed.json"
    destination.write_text(json.dumps(data), encoding="utf-8")
    return destination


def test_measured_patch_metrics_reproduce_without_any_dft():
    report = audit_patch(SOURCE)
    assert report["status"] == "five_point_metrics_recomputed_not_conditional_PES_certificate"
    assert report["measured_minus_predicted_energy_meV_per_BTO"] == pytest.approx(
        -0.47612973094, abs=1e-8,
    )
    assert report["predicted_positive_branch_Qy_sqrt_amu_A"] == pytest.approx(
        1.01631812960, abs=1e-9,
    )
    assert report["full_coordinate_error_sqrt_amu_A"] == pytest.approx(
        0.12542030938, abs=1e-9,
    )
    assert report["strain_coordinate_error_sqrt_amu_A"] > 5 * report[
        "atomic_coordinate_error_sqrt_amu_A"
    ]
    assert report["raw_selected_holdout_coordinate_error_sqrt_amu_A"] > 1.0
    assert max(report["aligned_neighbor_off_plane_jumps_sqrt_amu_A"][:2]) < 1.5
    assert min(report["raw_selected_neighbor_off_plane_jumps_sqrt_amu_A"][:2]) > 2.6
    assert len(report["mirror_checks"]) == 2
    assert all(check["energy_difference_eV_per_BTO"] < 1e-10
               for check in report["mirror_checks"].values())


def test_alternative_branch_requires_actual_mirror_equivalence(tmp_path):
    def break_mirror(data):
        data["points"][1]["aligned_positive_qy_coordinates_u_A_eta_voigt"][4] += 0.01

    with pytest.raises(ValueError, match="mirror-equivalent"):
        audit_patch(_modified_copy(tmp_path, break_mirror))


def test_alternative_branch_requires_matching_energy(tmp_path):
    def break_energy(data):
        data["points"][3]["aligned_positive_qy_energy_eV_per_BTO"] += 0.001

    with pytest.raises(ValueError, match="mirror-equivalent"):
        audit_patch(_modified_copy(tmp_path, break_energy))


def test_production_calculator_contract_is_pinned(tmp_path):
    def change_cutoff(data):
        data["calculator_contract"]["ecutwfc_Ry"] = 120

    with pytest.raises(ValueError, match="fixed calculator"):
        audit_patch(_modified_copy(tmp_path, change_cutoff))
