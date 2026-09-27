"""Score the pre-declared BTO (0.75, 0.15) conditional holdout, without DFT.

This is a *single-point interpolation* check, not a certification of a
continuous conditional potential-energy surface or of local curvature.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "paper" / "VARNEB_CPC" / "evidence" / "bto_q075_q015_holdout_2026-09-27"
PROTOCOL = ROOT / "docs" / "VARNEB_BTO_CONDITIONAL_HOLDOUT_PROTOCOL_2026-09-27.md"
FOUR_CORNERS = ROOT / "paper" / "VARNEB_CPC" / "figures" / "bto_conditional_four_point_stage_2026-09-27_source_data.csv"
PROTOCOL_SHA256 = "5e85bf7c93812f32a93942a1f54ed569d7f1de0215e7f341e1d9d6355f1d5b7a"
FOUR_CORNERS_SHA256 = "c1d2c39f0e8fae56660f0b2dc7aca385b916967b0c5fd73063cb0fcb16abbd63"
PREDICTED_LOW_EV = -0.106411735308
PREDICTED_ZERO_EV = -0.063858619978
LOW_ERROR_LIMIT_MEV = 2.0
SIGNED_BRANCH_LIMIT_MEV = 0.5
GRADIENT_LIMIT = 0.003
STRESS_LIMIT_KBAR = 1.0


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def analyze(evidence: Path = EVIDENCE, protocol: Path = PROTOCOL,
            corners: Path = FOUR_CORNERS) -> dict:
    if _sha256(protocol) != PROTOCOL_SHA256 or _sha256(corners) != FOUR_CORNERS_SHA256:
        raise ValueError("pre-result protocol or four-corner prediction source has changed")
    result_path = evidence / "conditional_q075_q015_result.json"
    audit_path = evidence / "audit-q075_q015-relax-27787487.json"
    canary_path = evidence / "audit-q075_q015-canary-27787471.json"
    replay_path = evidence / "audit-q075_q015-branch-replay-general-27787487.json"
    result, audit, canary, replay = map(_read, (result_path, audit_path, canary_path, replay_path))
    if (result.get("slurm_job_id") != "27787487"
            or result.get("q1_q2_sqrt_amu_A") != [0.75, 0.15]
            or result.get("branch_start_labels") != ["frozen", "+Q_y", "-Q_y"]
            or result.get("n_starts") != 3
            or result.get("status") != "orthogonal_gradient_and_stress_converged_curvature_unchecked"
            or result.get("stress_target_passed") is not True
            or result.get("n_new_evaluations") != 74
            or audit.get("status") != "verified_gradient_stationary_candidate_curvature_and_branches_unchecked"
            or audit.get("n_individually_audited_DFT_points") != 77
            or audit.get("source_sha256", {}).get("summary") != _sha256(result_path)
            or audit.get("source_sha256", {}).get("canary_audit") != _sha256(canary_path)
            or canary.get("status") != "all_three_raw_DFT_points_validated"
            or replay.get("status") != "all_three_branch_outcomes_reproduced_from_raw_audited_cache"
            or replay.get("n_cached_DFT_points") != 77
            or replay.get("n_replay_cache_hits") != 78
            or replay.get("source_sha256", {}).get("summary") != _sha256(result_path)
            or replay.get("source_sha256", {}).get("raw_audit") != _sha256(audit_path)):
        raise ValueError("holdout result, source hashes, or independent raw audit do not match")
    evaluations = audit["evaluations"]
    if (len(evaluations) != 77
            or sum(row["slurm_job_id"] == "27787471" for row in evaluations) != 3
            or sum(row["slurm_job_id"] == "27787487" for row in evaluations) != 74
            or any(row["mpi_dsize"] != 32 for row in evaluations)):
        raise ValueError("raw DFT point count, Slurm identity, or MPI contract differs")
    eligible = [row for row in evaluations
                if row["orthogonal_gradient_norm_eV_per_sqrt_amu_A"] <= GRADIENT_LIMIT
                and row["maximum_absolute_stress_kbar"] <= STRESS_LIMIT_KBAR]
    if len(eligible) != 3:
        raise ValueError("the three converged cache candidates are not unique")
    branches = {}
    for row in eligible:
        q_y = row["third_soft_y_amplitude_sqrt_amu_A"]
        label = "Q_y=0" if abs(q_y) < 1e-8 else ("+Q_y" if q_y > 0.5 else "-Q_y" if q_y < -0.5 else None)
        if label is None or label in branches:
            raise ValueError("converged cache point has ambiguous soft-mode basin")
        branches[label] = {key: row[key] for key in (
            "directory", "energy_minus_c_eV_per_BTO",
            "orthogonal_gradient_norm_eV_per_sqrt_amu_A",
            "maximum_atomic_force_eV_per_A", "maximum_absolute_stress_kbar",
            "minimum_distance_A", "third_soft_y_amplitude_sqrt_amu_A",
            "raw_log_sha256", "raw_input_sha256",
        )}
    if set(branches) != {"Q_y=0", "+Q_y", "-Q_y"}:
        raise ValueError("one or more expected basins are missing")
    replayed = {row["label"]: row for row in replay["branch_outcomes"]}
    if (set(replayed) != {"frozen", "+Q_y", "-Q_y"}
            or replayed["frozen"]["final_evaluation_directory"] != branches["Q_y=0"]["directory"]
            or any(replayed[label]["final_evaluation_directory"] != branches[label]["directory"]
                   for label in ("+Q_y", "-Q_y"))
            or any(row["converged"] is not True for row in replayed.values())):
        raise ValueError("deterministic branch replay differs from the unique audited basin points")
    plus, minus, zero = (branches[label] for label in ("+Q_y", "-Q_y", "Q_y=0"))
    low = min((plus, minus), key=lambda row: row["energy_minus_c_eV_per_BTO"])
    energy = low["energy_minus_c_eV_per_BTO"]
    zero_energy = zero["energy_minus_c_eV_per_BTO"]
    low_error_mev = 1000 * (energy - PREDICTED_LOW_EV)
    zero_error_mev = 1000 * (zero_energy - PREDICTED_ZERO_EV)
    signed_difference_mev = 1000 * abs(
        plus["energy_minus_c_eV_per_BTO"] - minus["energy_minus_c_eV_per_BTO"])
    single_holdout_pass = (abs(low_error_mev) <= LOW_ERROR_LIMIT_MEV
                           and signed_difference_mev <= SIGNED_BRANCH_LIMIT_MEV
                           and energy < zero_energy
                           and result["final_evaluation_directory"] == low["directory"])
    return {
        "kind": "BTO_single_internal_conditional_holdout_analysis_not_continuous_PES_or_barrier",
        "status": ("single_holdout_energy_prediction_and_three_branch_replay_passed_curvature_pending"
                   if single_holdout_pass else "single_holdout_requires_review"),
        "single_holdout_energy_prediction_passed": single_holdout_pass,
        "continuous_conditioned_PES_certified": False,
        "source_branch_terminals_explicit": bool(result.get("branch_outcomes")),
        "branch_terminal_mapping_replayed_from_audited_cache": True,
        "curvature_checked_at_holdout": False,
        "q_parallel_q_transverse_sqrt_amu_A": [0.75, 0.15],
        "prediction_low_eV_per_BTO": PREDICTED_LOW_EV,
        "observed_low_eV_per_BTO": energy,
        "low_prediction_error_meV_per_BTO": low_error_mev,
        "low_prediction_absolute_error_limit_meV_per_BTO": LOW_ERROR_LIMIT_MEV,
        "prediction_zero_eV_per_BTO": PREDICTED_ZERO_EV,
        "observed_zero_eV_per_BTO": zero_energy,
        "zero_prediction_error_meV_per_BTO": zero_error_mev,
        "observed_low_minus_zero_gap_meV_per_BTO": 1000 * (zero_energy - energy),
        "signed_branch_energy_difference_meV_per_BTO": signed_difference_mev,
        "branches_from_unique_audited_cache_points": branches,
        "branch_replay": replayed,
        "n_audited_DFT_points": len(evaluations),
        "n_new_DFT_points": result["n_new_evaluations"],
        "source_sha256": {
            "protocol": _sha256(protocol), "four_corners": _sha256(corners),
            "result": _sha256(result_path), "raw_audit": _sha256(audit_path),
            "canary_audit": _sha256(canary_path), "branch_replay": _sha256(replay_path),
        },
        "limitations": ("The frozen job preserved only its selected branch in the top-level result. "
                        "Three branch terminals were reconstructed by an exact deterministic replay "
                        "using only independently audited cached DFT points; no new DFT was run. "
                        "No holdout Hessian, global branch map, or continuous PES certificate exists."),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = analyze()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "status", "low_prediction_error_meV_per_BTO",
        "zero_prediction_error_meV_per_BTO", "signed_branch_energy_difference_meV_per_BTO",
    )}, indent=2))


if __name__ == "__main__":
    main()
