"""Keep the first audited BTO holdout Hessian at its actual claim level."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = (ROOT / "paper" / "VARNEB_CPC" / "evidence"
            / "bto_q075_q015_holdout_2026-09-27")


def _load(name: str) -> tuple[dict, str]:
    raw = (EVIDENCE / name).read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def test_first_holdout_curvature_is_raw_audited_but_not_stability_certified() -> None:
    parent, parent_sha = _load("conditional_q075_q015_result.json")
    replay, replay_sha = _load("audit-q075_q015-branch-replay-general-27787487.json")
    preflight, _ = _load("curvature_0p05_preflight.json")
    result, result_sha = _load("curvature_0p05_result.json")
    raw_audit, raw_audit_sha = _load("curvature_0p05_raw_audit.json")
    reconstruction, _ = _load("curvature_0p05_reconstruction.json")

    assert parent["selected_start"] == replay["selected_start"]
    assert preflight["status"] == "all_signed_probes_geometry_safe_not_a_curvature_result"
    assert preflight["n_signed_probes"] == 32
    assert preflight["source_sha256"]["summary"] == parent_sha
    assert preflight["source_sha256"]["branch_replay"] == replay_sha
    assert result["status"] == "orthogonal_curvature_screen_complete_not_branch_certification"
    assert result["slurm_job_id"] == "27787714"
    assert result["start_from_result_sha256"] == parent_sha
    assert result["n_new_evaluations"] == result["n_gradient_evaluations"] == 32

    assert raw_audit["n_individually_audited_DFT_points"] == 109
    assert len(raw_audit["evaluations"]) == 109
    assert len({row["directory"] for row in raw_audit["evaluations"]}) == 109
    assert sum(row["slurm_job_id"] == "27787714"
               for row in raw_audit["evaluations"]) == 32
    assert raw_audit["source_sha256"]["summary"] == parent_sha
    assert reconstruction["source_sha256"]["refined_result"] == parent_sha
    assert reconstruction["source_sha256"]["branch_replay"] == replay_sha
    assert reconstruction["source_sha256"]["curvature_result"] == result_sha
    assert reconstruction["source_sha256"]["all_points_audit"] == raw_audit_sha
    assert reconstruction["n_audited_DFT_points"] == 109
    assert reconstruction["n_matched_signed_probes"] == 32
    assert len(set(reconstruction["probe_evaluation_directories"])) == 32
    np.testing.assert_allclose(
        reconstruction["eigenvalues_eV_per_amu_A2"],
        result["orthogonal_hessian_eigenvalues_eV_per_amu_A2"],
        atol=1e-10, rtol=0.0,
    )
    assert reconstruction["status"] == "no_negative_curvature_at_one_difference_step_only"
    minimum = reconstruction["eigenvalues_eV_per_amu_A2"][0]
    assert 0.0 < minimum < 0.01
    assert reconstruction["energy_hessian_diagonal_max_abs_difference_eV_per_amu_A2"] > minimum


def test_second_step_geometry_is_preflight_only() -> None:
    parent, parent_sha = _load("conditional_q075_q015_result.json")
    preflight, _ = _load("curvature_0p10_preflight.json")
    assert parent["status"] == "orthogonal_gradient_and_stress_converged_curvature_unchecked"
    assert preflight["status"] == "all_signed_probes_geometry_safe_not_a_curvature_result"
    assert preflight["source_sha256"]["summary"] == parent_sha
    assert preflight["n_signed_probes"] == 32
    assert preflight["curvature_step_sqrt_amu_A"] == 0.10
    assert preflight["no_dft_launched"] is True
