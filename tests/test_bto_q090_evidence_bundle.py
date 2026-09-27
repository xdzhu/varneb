"""The archived BTO Qz=0.9 single-step screen is not a PES certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest


BUNDLE = (Path(__file__).resolve().parents[1] / "benchmarks" /
          "numerical_integrity" / "bto_q090_hf_2026-09-27")


def _load(name: str) -> tuple[dict, str]:
    path = BUNDLE / name
    payload = path.read_bytes()
    return json.loads(payload), hashlib.sha256(payload).hexdigest()


@pytest.mark.parametrize("point,q,n_eval", [
    ("q090_q000", [0.9, 0.0], 77),
    ("q090_q030", [0.9, 0.3], 78),
])
def test_raw_audit_replay_and_curvature_preflight_chain(point, q, n_eval):
    summary, summary_sha = _load(f"conditional_{point}_result.json")
    preflight, preflight_sha = _load(
        f"bto_transverse_soft_conditional_{point}_preflight_2026-09-27.json")
    audit, audit_sha = _load(f"audit-{point}-production-raw-final.json")
    replay, replay_sha = _load(f"replay-{point}-final.json")
    probe, _ = _load(
        f"bto_transverse_soft_{point}_selected_curvature_0p05_preflight_2026-09-27.json")

    assert summary["status"] == "orthogonal_gradient_and_stress_converged_curvature_unchecked"
    assert summary["q_parallel_q_transverse_sqrt_amu_A"] == q
    assert summary["phonon_supercell"] == [1, 1, 1]
    assert summary["electronic_kpoints"] == [4, 4, 4]
    assert summary["stress_target_passed"] is True
    assert summary["stress_target_kbar"] == 1.0
    assert summary["preflight_sha256"] == preflight_sha
    assert preflight["calculator_id"].startswith("abacus-bto-pbe100-dzp10au-4x4x4-")

    assert audit["status"] == "verified_gradient_stationary_candidate_curvature_and_branches_unchecked"
    assert audit["source_sha256"]["summary"] == summary_sha
    assert audit["source_sha256"]["preflight"] == preflight_sha
    assert audit["n_individually_audited_DFT_points"] == n_eval
    assert len(audit["evaluations"]) == n_eval
    assert all(item["raw_log_sha256"] and len(item["raw_input_sha256"]) == 3
               and item["mpi_dsize"] == 32 for item in audit["evaluations"])

    assert replay["status"] == "all_three_branch_outcomes_reproduced_from_raw_audited_cache"
    assert replay["source_sha256"]["summary"] == summary_sha
    assert replay["source_sha256"]["raw_audit"] == audit_sha
    assert replay["selected_start"] == summary["selected_start"] == 2
    outcomes = replay["branch_outcomes"]
    assert [row["label"] for row in outcomes] == ["frozen", "+Q_y", "-Q_y"]
    assert all(row["converged"] for row in outcomes)
    assert outcomes[0]["energy_minus_c_eV_per_BTO"] > (
        outcomes[2]["energy_minus_c_eV_per_BTO"] + 0.02)
    assert outcomes[1]["energy_minus_c_eV_per_BTO"] == pytest.approx(
        outcomes[2]["energy_minus_c_eV_per_BTO"], abs=1e-7)
    assert outcomes[2]["energy_minus_c_eV_per_BTO"] == pytest.approx(
        summary["energy_minus_c_eV_per_BTO"], abs=1e-10)

    assert probe["status"] == "all_signed_probes_geometry_safe_not_a_curvature_result"
    assert probe["source_sha256"]["summary"] == summary_sha
    assert probe["source_sha256"]["audit"] == audit_sha
    assert probe["source_sha256"]["branch_replay"] == replay_sha
    assert probe["n_signed_probes"] == 32
    assert probe["curvature_step_sqrt_amu_A"] == 0.05
    assert probe["minimum_probe_atomic_distance_A"] > 1.6


@pytest.mark.parametrize("point,q,n_eval", [
    ("q090_q000", [0.9, 0.0], 109),
    ("q090_q030", [0.9, 0.3], 110),
])
def test_single_step_curvature_is_raw_audited_but_not_certified(point, q, n_eval):
    parent, parent_sha = _load(f"conditional_{point}_result.json")
    replay, replay_sha = _load(f"replay-{point}-final.json")
    all_points, all_points_sha = _load(
        f"audit-{point}-all-points-through-curvature-0p05-27787086.json")
    result, result_sha = _load(f"conditional_{point}_selected_curvature_0p05.json")
    audit, _ = _load(f"audit-{point}-curvature-0p05-27787086.json")

    assert all_points["status"] == "verified_gradient_stationary_candidate_curvature_and_branches_unchecked"
    assert all_points["n_individually_audited_DFT_points"] == n_eval
    assert len(all_points["evaluations"]) == n_eval
    assert all(item["raw_log_sha256"] and len(item["raw_input_sha256"]) == 3
               and item["mpi_dsize"] == 32 for item in all_points["evaluations"])
    assert result["status"] == "orthogonal_curvature_screen_complete_not_branch_certification"
    assert result["q1_q2_sqrt_amu_A"] == q
    assert result["curvature_step_sqrt_amu_A"] == 0.05
    assert result["n_gradient_evaluations"] == 32
    assert result["start_from_result_sha256"] == parent_sha
    assert result["evaluator_contract_sha256"] == parent["evaluator_contract_sha256"]
    assert audit["status"] == "no_negative_curvature_at_one_difference_step_only"
    assert audit["n_audited_DFT_points"] == n_eval
    assert audit["n_matched_signed_probes"] == 32
    assert audit["source_sha256"]["refined_result"] == parent_sha
    assert audit["source_sha256"]["curvature_result"] == result_sha
    assert audit["source_sha256"]["all_points_audit"] == all_points_sha
    assert audit["source_sha256"]["branch_replay"] == replay_sha
    assert audit["negative_directions"] == []
    assert min(audit["eigenvalues_eV_per_amu_A2"]) > 0
    assert audit["energy_hessian_diagonal_max_abs_difference_eV_per_amu_A2"] > min(
        audit["eigenvalues_eV_per_amu_A2"])


@pytest.mark.parametrize("point,q", [
    ("q090_q000", [0.9, 0.0]),
    ("q090_q030", [0.9, 0.3]),
])
def test_second_step_preflight_is_geometry_only(point, q):
    parent, parent_sha = _load(f"conditional_{point}_result.json")
    initial_audit, initial_audit_sha = _load(f"audit-{point}-production-raw-final.json")
    replay, replay_sha = _load(f"replay-{point}-final.json")
    preflight, _ = _load(
        f"bto_transverse_soft_{point}_selected_curvature_0p10_preflight_2026-09-27.json")
    assert preflight["status"] == "all_signed_probes_geometry_safe_not_a_curvature_result"
    assert preflight["q1_q2_sqrt_amu_A"] == q
    assert preflight["n_signed_probes"] == 32
    assert preflight["curvature_step_sqrt_amu_A"] == 0.10
    assert preflight["minimum_probe_atomic_distance_A"] > 1.6
    assert preflight["minimum_probe_volume_A3"] > 0
    assert preflight["no_dft_launched"] is True
    assert preflight["source_sha256"]["summary"] == parent_sha
    assert preflight["source_sha256"]["audit"] == initial_audit_sha
    assert preflight["source_sha256"]["branch_replay"] == replay_sha
    assert parent["stress_target_passed"] and initial_audit["source_sha256"]["summary"] == parent_sha
    assert replay["selected_start"] == parent["selected_start"]
