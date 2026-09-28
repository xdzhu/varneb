"""Offline checks for the BTO Qy=0 holdout gate (no ssh or DFT)."""

from __future__ import annotations

import pytest
import numpy as np

from scripts.audit_bto_qy0_preregistered_holdouts import gate_one


def fixture() -> tuple[dict, dict, dict, dict, dict]:
    q = [0.15, 0.15]
    plan = {"q_parallel_q_transverse_sqrt_amu_A": q,
            "interpolation": "bilinear_center",
            "predicted_energy_minus_c_eV_per_BTO": -0.013,
            "predicted_coordinates_u_A_eta_voigt": [0.0] * 21}
    summary = {
        "kind": "bto_symmetry_restricted_qy_zero_variable_cell_local_candidate_not_PES_or_barrier",
        "status": "orthogonal_gradient_and_stress_converged_curvature_unchecked",
        "preflight_sha256": "pre", "phonon_supercell": [1, 1, 1],
        "electronic_kpoints": [4, 4, 4], "stress_target_kbar": 2.0,
        "stress_target_passed": True, "q_parallel_q_transverse_sqrt_amu_A": q,
        "energy_minus_c_eV_per_BTO": -0.007,
        "orthogonal_gradient_norm_eV_per_sqrt_amu_A": 0.001,
        "maximum_absolute_stress_kbar": 1.0,
        "coordinates_u_A_eta_voigt": [0.0] * 21,
    }
    audit = {
        "status": "verified_gradient_stationary_candidate_curvature_and_branches_unchecked",
        "source_sha256": {"summary": "sum", "preflight": "pre"},
        "third_soft_mode_restricted_at_zero": True,
        "n_individually_audited_DFT_points": 1, "evaluations": [{}],
        "final_third_soft_y_amplitude_sqrt_amu_A": 0.0,
        "q_parallel_q_transverse_sqrt_amu_A": q,
        "final_energy_minus_c_eV_per_BTO": -0.007,
        "orthogonal_gradient_norm_eV_per_sqrt_amu_A": 0.001,
        "final_maximum_absolute_stress_kbar": 1.0,
    }
    preflight = {
        "no_dft_launched": True, "q_parallel_q_transverse_sqrt_amu_A": q,
        "strain_metric_weights_amu_A2": [100.0] * 6,
    }
    gates = {"energy_absolute_error_gate_meV_per_BTO": 2.0,
             "full_atom_plus_strain_metric_coordinate_error_gate_sqrt_amu_A": 0.1,
             "gradient_gate_eV_per_sqrt_amu_A": 0.003,
             "maximum_absolute_stress_gate_kbar": 2.0}
    return plan, summary, audit, preflight, gates


def test_energy_gate_fails_without_rejecting_valid_raw_calculation() -> None:
    plan, summary, audit, preflight, gates = fixture()
    row = gate_one(plan, summary, audit, preflight, "sum", "pre",
                   np.array([137.327, 47.867, 15.999, 15.999, 15.999]), gates)
    assert row["audited_minus_predicted_energy_meV_per_BTO"] == pytest.approx(6.0)
    assert row["checks"] == {"energy": False, "full_coordinate": True,
                             "gradient": True, "stress": True}
    assert not row["all_gates_pass"]


def test_gate_rejects_unmatched_raw_audit() -> None:
    plan, summary, audit, preflight, gates = fixture()
    audit["source_sha256"]["summary"] = "wrong"
    with pytest.raises(ValueError, match="raw-audited holdout contract"):
        gate_one(plan, summary, audit, preflight, "sum", "pre",
                 np.array([137.327, 47.867, 15.999, 15.999, 15.999]), gates)
