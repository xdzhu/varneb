from __future__ import annotations

import json

import pytest

from scripts.build_ase_static_summary import build_summary


def _write_summary(tmp_path, name: str, *, converged: bool = True, n_atoms: int = 4):
    path = tmp_path / name
    path.write_text(
        json.dumps(
            {
                "status": "completed",
                "converged": converged,
                "endpoint": {
                    "n_atoms": n_atoms,
                    "composition": {"Ba": 1, "Ti": 1, "O": 2},
                },
                "potential_energy_eV": -10.0,
                "max_atomic_force_eV_per_A": 0.04,
                "max_generalized_force_eV_per_A": 0.05,
                "stress_eV_per_A3": [[0.0, 0.0, 0.0]] * 3,
            }
        ),
        encoding="utf-8",
    )
    return path


def test_build_summary_creates_static_gate_contract(tmp_path) -> None:
    initial = _write_summary(tmp_path, "initial.json")
    final = _write_summary(tmp_path, "final.json")

    summary = build_summary(initial, final)

    assert summary["status"] == "static_completed"
    assert summary["endpoint_evaluation_policy"] == "independent_bfgs_relaxations"
    assert [entry["label"] for entry in summary["static_endpoints"]] == [
        "initial",
        "final",
    ]
    assert summary["endpoint_structures"]["initial"]["n_atoms"] == 4


def test_build_summary_rejects_unconverged_endpoint(tmp_path) -> None:
    initial = _write_summary(tmp_path, "initial.json", converged=False)
    final = _write_summary(tmp_path, "final.json")

    with pytest.raises(ValueError, match="not converged"):
        build_summary(initial, final)


def test_build_summary_rejects_mismatched_atom_count(tmp_path) -> None:
    initial = _write_summary(tmp_path, "initial.json")
    final = _write_summary(tmp_path, "final.json", n_atoms=5)

    with pytest.raises(ValueError, match="different atom counts"):
        build_summary(initial, final)
