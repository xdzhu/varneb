from __future__ import annotations

import pytest

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
