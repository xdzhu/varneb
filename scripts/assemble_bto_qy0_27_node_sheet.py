"""Assemble the restricted BTO 9x3 sheet only from individually raw-audited points."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def assemble(points_root: Path, evidence_root: Path, measured_csv: Path) -> dict:
    plan_path = points_root / "plan.json"
    plan = load(plan_path)
    if (plan.get("status") != "15_missing_Qy0_inputs_preflighted_no_DFT"
            or plan.get("n_grid_nodes") != 27
            or plan.get("n_reused_measured_nodes") != 12
            or plan.get("n_new_conditional_points") != 15
            or plan.get("source_sha256", {}).get("measured_csv") != sha256(measured_csv)):
        raise ValueError("restricted-sheet source plan changed")
    grid: dict[tuple[int, int], dict] = {}
    for row in plan["reused"]:
        key = tuple(row["grid_index"])
        if key in grid:
            raise ValueError("duplicate previously measured grid coordinate")
        grid[key] = {
            "grid_index": list(key), "q_sqrt_amu_A": row["q_sqrt_amu_A"],
            "energy_minus_C_meV_per_BTO": row["measured_minus_C_meV_per_BTO"],
            "source": "previously_raw_audited_DFT", "class": row["class"],
        }
    new = []
    contracts = set()
    for row in plan["new"]:
        point = row["name"]
        preflight_path = points_root / "preflights" / f"{point}.json"
        preflight = load(preflight_path)
        stress_refined = point in {"qz01_qx02", "qz02_qx01"}
        initial_path = evidence_root / f"run-{point}" / f"{point}_result.json"
        summary_name = (f"{point}_stress_refined_result.json" if stress_refined
                        else f"{point}_result.json")
        summary_path = evidence_root / f"run-{point}" / summary_name
        audit_path = evidence_root / f"audit-{point}.json"
        summary, audit = load(summary_path), load(audit_path)
        q = row["q_sqrt_amu_A"]
        if (sha256(preflight_path) != row["preflight_sha256"]
                or preflight.get("calculator_id") != plan["calculator_id"]
                or preflight.get("q_parallel_q_transverse_sqrt_amu_A") != q
                or preflight.get("kind")
                   != "bto_symmetry_restricted_soft_qy_zero_preflight_no_dft"
                or preflight.get("phonon_supercell") != [1, 1, 1]
                or preflight.get("electronic_kpoints") != [4, 4, 4]
                or summary.get("status")
                   != "orthogonal_gradient_and_stress_converged_curvature_unchecked"
                or summary.get("preflight_sha256") != sha256(preflight_path)
                or summary.get("q_parallel_q_transverse_sqrt_amu_A") != q
                or summary.get("stress_target_kbar") != 2.0
                or summary.get("stress_target_passed") is not True
                or audit.get("status")
                   != "verified_gradient_stationary_candidate_curvature_and_branches_unchecked"
                or audit.get("source_sha256", {}).get("summary") != sha256(summary_path)
                or audit.get("source_sha256", {}).get("preflight") != sha256(preflight_path)
                or audit.get("q_parallel_q_transverse_sqrt_amu_A") != q
                or audit.get("third_soft_mode_restricted_at_zero") is not True
                or abs(audit.get("final_third_soft_y_amplitude_sqrt_amu_A", np.inf)) > 1e-8
                or audit.get("n_individually_audited_DFT_points", 0) < 1):
            raise ValueError(f"new restricted-sheet point not raw-audited: {point}")
        if stress_refined:
            initial = load(initial_path)
            if (initial.get("status")
                    != "orthogonal_gradient_converged_stress_target_failed_curvature_unchecked"
                    or summary.get("start_from_result_sha256") != sha256(initial_path)
                    or summary.get("gradient_tolerance_eV_per_sqrt_amu_A") != 0.0015):
                raise ValueError(f"failed-stress continuation provenance invalid: {point}")
        elif summary.get("gradient_tolerance_eV_per_sqrt_amu_A") != 0.003:
            raise ValueError(f"unplanned optimizer tolerance at {point}")
        gradient = float(audit["orthogonal_gradient_norm_eV_per_sqrt_amu_A"])
        stress = float(audit["final_maximum_absolute_stress_kbar"])
        energy = float(audit["final_energy_minus_c_eV_per_BTO"]) * 1000
        if (gradient > (0.0015 if stress_refined else 0.003) + 1e-12
                or stress > 2.0 + 1e-9
                or abs(energy - 1000 * float(summary["energy_minus_c_eV_per_BTO"])) > 1e-5
                or not np.isfinite([gradient, stress, energy]).all()):
            raise ValueError(f"new point fails gradient/stress/energy audit: {point}")
        contracts.add(summary["evaluator_contract_sha256"])
        predicted = float(row["frozen_model_prediction_minus_C_meV_per_BTO"])
        key = tuple(row["grid_index"])
        if key in grid:
            raise ValueError(f"duplicate new grid coordinate: {point}")
        result = {
            "name": point, "grid_index": list(key), "q_sqrt_amu_A": q,
            "energy_minus_C_meV_per_BTO": energy,
            "frozen_nine_node_model_prediction_meV_per_BTO": predicted,
            "DFT_minus_model_meV_per_BTO": energy - predicted,
            "orthogonal_gradient_eV_per_sqrt_amu_A": gradient,
            "maximum_absolute_stress_kbar": stress,
            "n_individually_raw_audited_DFT_evaluations": audit[
                "n_individually_audited_DFT_points"],
            "source": "new_raw_audited_DFT", "stress_refined": stress_refined,
            "summary_sha256": sha256(summary_path), "audit_sha256": sha256(audit_path),
        }
        new.append(result)
        grid[key] = result
    if len(contracts) != 1 or len(grid) != 27 or len(new) != 15:
        raise ValueError("the 27-point sheet is incomplete or uses mixed evaluator contracts")
    energy = [[grid[(i, j)]["energy_minus_C_meV_per_BTO"] for j in range(3)]
              for i in range(9)]
    errors = np.array([row["DFT_minus_model_meV_per_BTO"] for row in new])
    return {
        "status": "BTO_Qy0_restricted_27_node_sheet_raw_audited",
        "claim_limit": "one symmetry-restricted Qy=0 local branch; not the unrestricted conditional PES",
        "calculator_id": plan["calculator_id"],
        "evaluator_contract_sha256": next(iter(contracts)),
        "grid_qz_sqrt_amu_A": plan["grid_qz_sqrt_amu_A"],
        "grid_qx_sqrt_amu_A": plan["grid_qx_sqrt_amu_A"],
        "energy_minus_C_meV_per_BTO": energy,
        "n_reused_measured_nodes": 12, "n_new_raw_audited_nodes": 15,
        "n_total_measured_nodes": 27,
        "n_new_raw_audited_DFT_evaluations": sum(
            row["n_individually_raw_audited_DFT_evaluations"] for row in new),
        "predeclared_model_error_gate_meV_per_BTO": 2.0,
        "new_node_model_max_abs_error_meV_per_BTO": float(np.max(np.abs(errors))),
        "new_node_model_rms_error_meV_per_BTO": float(np.sqrt(np.mean(errors**2))),
        "new_node_model_gate_pass": bool(np.max(np.abs(errors)) <= 2.0),
        "reused": plan["reused"], "new": new,
        "source_sha256": {
            "plan": sha256(plan_path), "measured_csv": sha256(measured_csv),
            "assembler": sha256(Path(__file__)),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("points-root", "evidence-root", "measured-csv", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = assemble(args.points_root, args.evidence_root, args.measured_csv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "status", "n_new_raw_audited_nodes", "n_new_raw_audited_DFT_evaluations",
        "new_node_model_max_abs_error_meV_per_BTO", "new_node_model_gate_pass",
    )}))


if __name__ == "__main__":
    main()
