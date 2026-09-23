from __future__ import annotations

import json

from ase import Atoms
import pytest

from examples.run_vcneb_ase import _validate_static_endpoint_identity
from vcneb import endpoint_structure_record
from scripts.validate_ase_static_gate import validate_summary


def _summary(force: float, stress: float) -> dict:
    return {
        "status": "static_completed",
        "endpoint_structures": {
            "initial": {"n_atoms": 4, "composition": {"Ba": 1, "Ti": 1, "O": 2}},
            "final": {"n_atoms": 4, "composition": {"Ba": 1, "Ti": 1, "O": 2}},
        },
        "static_endpoints": [
            {
                "label": "initial",
                "max_force_eV_per_A": force,
                "stress_eV_per_A3": [[stress, 0.0, 0.0], [0.0, stress, 0.0], [0.0, 0.0, stress]],
                "n_atoms": 4,
            },
            {
                "label": "final",
                "max_force_eV_per_A": force,
                "stress_eV_per_A3": [[stress, 0.0, 0.0], [0.0, stress, 0.0], [0.0, 0.0, stress]],
                "n_atoms": 4,
            },
        ],
    }


def test_static_gate_passes_only_when_force_and_stress_are_below_targets() -> None:
    report = validate_summary(_summary(0.05, 0.00001))
    assert report["status"] == "passed"
    assert report["issues"] == []


def test_default_static_gate_accepts_sub_kbar_residual_stress() -> None:
    report = validate_summary(_summary(0.05, 0.00057))
    assert report["status"] == "passed"
    assert report["endpoints"][0]["max_abs_stress_residual_kbar"] == pytest.approx(
        0.91324068138
    )


@pytest.mark.parametrize(
    "summary",
    [_summary(0.10, 0.00001), _summary(0.05, 0.001)],
)
def test_static_gate_rejects_boundary_or_high_stress(summary: dict) -> None:
    assert validate_summary(summary)["status"] == "failed"


def test_static_gate_rejects_endpoint_identity_mismatch() -> None:
    summary = _summary(0.05, 0.00001)
    summary["endpoint_structures"]["final"]["n_atoms"] = 5
    report = validate_summary(summary)
    assert report["status"] == "failed"
    assert "atom counts differ" in " ".join(report["issues"])


def test_static_gate_compares_stress_with_nonzero_target_pressure() -> None:
    pressure_gpa = 45.7
    target_stress = -pressure_gpa / 160.21766208
    summary = _summary(0.02, target_stress)
    summary["external_pressure_gpa"] = pressure_gpa
    report = validate_summary(summary, pressure_gpa=pressure_gpa)
    assert report["status"] == "passed"
    assert report["endpoints"][0]["max_abs_stress_kbar"] > 450.0
    assert report["endpoints"][0]["max_abs_stress_residual_kbar"] == pytest.approx(0.0)
    assert validate_summary(summary, pressure_gpa=0.0)["status"] == "failed"
    summary["external_pressure_gpa"] = 45.0
    assert validate_summary(summary, pressure_gpa=pressure_gpa)["status"] == "failed"


def test_static_gate_identity_must_match_effective_aligned_endpoints(tmp_path) -> None:
    initial = Atoms("H", positions=[[0.0, 0.0, 0.0]], cell=[4.0, 4.0, 4.0], pbc=True)
    final = Atoms("H", positions=[[1.0, 0.0, 0.0]], cell=[4.0, 4.0, 4.0], pbc=True)
    summary = {
        "endpoint_structures": {
            "initial": endpoint_structure_record(initial),
            "final": endpoint_structure_record(final),
        }
    }
    path = tmp_path / "static.json"
    path.write_text(json.dumps(summary), encoding="utf-8")
    assert _validate_static_endpoint_identity(path, [initial, final])["status"] == "passed"

    translated = final.copy()
    translated.positions += [0.25, 0.0, 0.0]
    with pytest.raises(ValueError, match="mapped/aligned"):
        _validate_static_endpoint_identity(path, [initial, translated])
