"""Offline contract tests for the prospective even-mode holdout reader."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from scripts.audit_bto_qy0_even_mode_holdouts import evaluate


ROOT = Path(__file__).resolve().parents[1]


def test_one_prospective_point_passes_only_when_model_and_raw_audit_match() -> None:
    frozen = json.loads((ROOT / "benchmarks" / "numerical_integrity" /
                         "bto_qy0_even_mode_two_holdout_plan_20260928.json").read_text())
    point = frozen["new_holdouts"][0]
    q = point["q_sqrt_amu_A"]
    energy = point["model_prediction_energy_minus_C_eV_per_BTO"] + 0.001
    summary = {
        "status": "orthogonal_gradient_and_stress_converged_curvature_unchecked",
        "kind": "bto_symmetry_restricted_qy_zero_variable_cell_local_candidate_not_PES_or_barrier",
        "preflight_sha256": "pre", "phonon_supercell": [1, 1, 1],
        "electronic_kpoints": [4, 4, 4], "stress_target_kbar": 2.0,
        "stress_target_passed": True, "q_parallel_q_transverse_sqrt_amu_A": q,
        "energy_minus_c_eV_per_BTO": energy,
        "orthogonal_gradient_norm_eV_per_sqrt_amu_A": 0.001,
        "maximum_absolute_stress_kbar": 1.0,
        "coordinates_u_A_eta_voigt": point["predicted_structure_u_A_eta_voigt"],
    }
    audit = {
        "status": "verified_gradient_stationary_candidate_curvature_and_branches_unchecked",
        "third_soft_mode_restricted_at_zero": True,
        "final_third_soft_y_amplitude_sqrt_amu_A": 0.0,
        "n_individually_audited_DFT_points": 1, "evaluations": [{}],
        "source_sha256": {"summary": "sum", "preflight": "pre"},
        "q_parallel_q_transverse_sqrt_amu_A": q,
        "final_energy_minus_c_eV_per_BTO": energy,
        "orthogonal_gradient_norm_eV_per_sqrt_amu_A": 0.001,
        "final_maximum_absolute_stress_kbar": 1.0,
    }
    preflight = {"no_dft_launched": True,
                 "q_parallel_q_transverse_sqrt_amu_A": q,
                 "strain_metric_weights_amu_A2": [100.0] * 6}
    masses = np.array([137.327, 47.867, 15.999, 15.999, 15.999])
    row = evaluate(point, summary, audit, preflight, "sum", "pre", masses, frozen)
    assert row["all_gates_pass"]
    assert row["audited_minus_predicted_meV_per_BTO"] == pytest.approx(1.0)
    audit["source_sha256"]["summary"] = "wrong"
    with pytest.raises(ValueError, match="raw-audit contract mismatch"):
        evaluate(point, summary, audit, preflight, "sum", "pre", masses, frozen)
