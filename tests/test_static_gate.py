from __future__ import annotations

import pytest

from scripts.validate_ase_static_gate import validate_summary


def _summary(force: float, stress: float) -> dict:
    return {
        "status": "static_completed",
        "static_endpoints": [
            {
                "label": "initial",
                "max_force_eV_per_A": force,
                "stress_eV_per_A3": [[stress, 0.0, 0.0], [0.0, stress, 0.0], [0.0, 0.0, stress]],
            },
            {
                "label": "final",
                "max_force_eV_per_A": force,
                "stress_eV_per_A3": [[stress, 0.0, 0.0], [0.0, stress, 0.0], [0.0, 0.0, stress]],
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
