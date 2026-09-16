"""Regression checks for the non-DFT QE static-series audit."""

from __future__ import annotations

import runpy
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_qe_static_convergence.py"


def _summary(ecutwfc: float, *, energy: float, force: float, stress: float, endpoint: str = "endpoint-a") -> dict:
    return {
        "_source": f"cutoff-{ecutwfc}.json",
        "status": "completed",
        "execution_mode": "fixed_initial_endpoint_static_scf",
        "evaluated_image_index": 0,
        "n_images": 7,
        "n_interior_images": 5,
        "git_revision": "4eae070",
        "endpoint_structures": {"initial": {"sha256": endpoint, "n_atoms": 5}, "final": {"sha256": "endpoint-final"}},
        "calculator_parameters": {
            "kpts": [4, 4, 4],
            "input_data": {"system": {"ecutwfc": ecutwfc, "ecutrho": ecutwfc * 6}, "electrons": {"conv_thr": 1e-8}, "control": {"calculation": "scf"}},
        },
        "pseudopotential_reports": [{"species": "Ba", "filename": "Ba.UPF", "md5": "a", "expected_md5": "a"}],
        "potential_energy_eV": energy,
        "forces_eV_per_A": [[force, 0.0, 0.0]] * 5,
        "stress_eV_per_A3_voigt": [stress, 0.0, 0.0, 0.0, 0.0, 0.0],
    }


def test_static_qe_audit_accepts_only_declared_final_cutoff_pair() -> None:
    audit = runpy.run_path(str(SCRIPT))["audit"]
    report = audit(
        [_summary(80, energy=-10.0000, force=0.10, stress=0.010), _summary(100, energy=-10.0050, force=0.102, stress=0.011), _summary(120, energy=-10.0055, force=0.1025, stress=0.0115)],
        energy_tol_mev_per_atom=1.0,
        force_tol_ev_per_a=0.001,
        stress_tol_gpa=0.1,
    )
    assert report["status"] == "ok"
    assert report["accepted_final_pair"]["energy_difference_meV_per_atom"] == pytest.approx(0.1)
    assert report["cutoff_points_Ry"] == [(80.0, 480.0), (100.0, 600.0), (120.0, 720.0)]


def test_static_qe_audit_rejects_mixed_endpoints() -> None:
    audit = runpy.run_path(str(SCRIPT))["audit"]
    with pytest.raises(ValueError, match="identical endpoints"):
        audit(
            [_summary(80, energy=-10, force=0.1, stress=0.01), _summary(100, energy=-10, force=0.1, stress=0.01), _summary(120, energy=-10, force=0.1, stress=0.01, endpoint="different")],
            energy_tol_mev_per_atom=1.0,
            force_tol_ev_per_a=0.01,
            stress_tol_gpa=1.0,
        )
