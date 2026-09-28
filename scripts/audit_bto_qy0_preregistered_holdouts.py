"""Gate the three pre-registered BTO Qy=0 interpolation holdouts.

Reads only the frozen plan and previously completed, independently audited
ABACUS output. Prints a compact JSON report; never submits DFT or edits source
calculations. The source audits themselves check each raw SCF/force/stress.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from io import StringIO
import json
from pathlib import Path
import shlex
import subprocess

import numpy as np
from ase.io import read


CASES = (
    ("q015_q015", "run-holdout-q015_q015/q015_q015_holdout_tight_result.json",
     "audit-holdout-q015_q015-27795015.json",
     "inputs/bto_qy0_q015_q015_holdout_preflight_2026-09-28.json"),
    ("q045_q015", "run-holdout-q045_q015/q045_q015_holdout_result.json",
     "audit-holdout-q045_q015-27794991_1.json",
     "inputs/bto_qy0_q045_q015_holdout_preflight_2026-09-28.json"),
    ("q105_q020", "run-holdout-q105_q020/q105_q020_holdout_result.json",
     "audit-holdout-q105_q020-27794991_2.json",
     "inputs/bto_qy0_q105_q020_holdout_preflight_2026-09-28.json"),
)


def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def fetch(host: str, root: str, relative: str) -> bytes:
    remote = f"{root.rstrip('/')}/{relative}"
    return subprocess.run(
        ["ssh", host, f"cat {shlex.quote(remote)}"],
        check=True, capture_output=True,
    ).stdout


def gate_one(plan_point: dict, summary: dict, audit: dict, preflight: dict,
             summary_hash: str, preflight_hash: str, atomic_masses: np.ndarray,
             gates: dict) -> dict:
    q = np.asarray(plan_point["q_parallel_q_transverse_sqrt_amu_A"], dtype=float)
    if (summary.get("kind")
            != "bto_symmetry_restricted_qy_zero_variable_cell_local_candidate_not_PES_or_barrier"
            or summary.get("status") != "orthogonal_gradient_and_stress_converged_curvature_unchecked"
            or audit.get("status") != "verified_gradient_stationary_candidate_curvature_and_branches_unchecked"
            or preflight.get("no_dft_launched") is not True
            or summary.get("preflight_sha256") != preflight_hash
            or audit.get("source_sha256", {}).get("summary") != summary_hash
            or audit.get("source_sha256", {}).get("preflight") != preflight_hash
            or summary.get("phonon_supercell") != [1, 1, 1]
            or summary.get("electronic_kpoints") != [4, 4, 4]
            or summary.get("stress_target_kbar") != 2.0
            or summary.get("stress_target_passed") is not True
            or audit.get("third_soft_mode_restricted_at_zero") is not True
            or audit.get("n_individually_audited_DFT_points", 0) < 1
            or audit.get("n_individually_audited_DFT_points") != len(audit.get("evaluations", []))
            or abs(audit.get("final_third_soft_y_amplitude_sqrt_amu_A", np.inf)) > 1e-8
            or not np.allclose(summary.get("q_parallel_q_transverse_sqrt_amu_A"), q,
                               atol=1e-12, rtol=0)
            or not np.allclose(audit.get("q_parallel_q_transverse_sqrt_amu_A"), q,
                               atol=1e-12, rtol=0)
            or not np.allclose(preflight.get("q_parallel_q_transverse_sqrt_amu_A"), q,
                               atol=1e-12, rtol=0)):
        raise ValueError(f"raw-audited holdout contract does not match {q.tolist()}")
    energy = float(summary["energy_minus_c_eV_per_BTO"])
    gradient = float(summary["orthogonal_gradient_norm_eV_per_sqrt_amu_A"])
    stress = float(audit["final_maximum_absolute_stress_kbar"])
    if (abs(energy - audit["final_energy_minus_c_eV_per_BTO"]) > 1e-8
            or abs(gradient - audit["orthogonal_gradient_norm_eV_per_sqrt_amu_A"]) > 1e-8
            or abs(stress - summary["maximum_absolute_stress_kbar"]) > 1e-8):
        raise ValueError(f"summary/audit mismatch at {q.tolist()}")
    predicted = np.asarray(plan_point["predicted_coordinates_u_A_eta_voigt"], dtype=float)
    observed = np.asarray(summary["coordinates_u_A_eta_voigt"], dtype=float)
    strain_weights = np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float)
    weights = np.r_[np.repeat(atomic_masses, 3), strain_weights]
    if (predicted.shape != (21,) or observed.shape != (21,)
            or strain_weights.shape != (6,) or np.any(weights <= 0)
            or not np.all(np.isfinite(predicted)) or not np.all(np.isfinite(observed))):
        raise ValueError("invalid full atom-plus-strain coordinate metric")
    coordinate_error = float(np.sqrt(np.dot(weights, (observed - predicted) ** 2)))
    energy_error = 1000.0 * (energy - float(plan_point["predicted_energy_minus_c_eV_per_BTO"]))
    checks = {
        "energy": abs(energy_error) <= gates["energy_absolute_error_gate_meV_per_BTO"],
        "full_coordinate": coordinate_error <= gates[
            "full_atom_plus_strain_metric_coordinate_error_gate_sqrt_amu_A"],
        "gradient": gradient <= gates["gradient_gate_eV_per_sqrt_amu_A"],
        "stress": stress <= gates["maximum_absolute_stress_gate_kbar"],
    }
    return {
        "q_sqrt_amu_A": q.tolist(),
        "interpolation": plan_point["interpolation"],
        "predicted_energy_minus_C_eV_per_BTO": plan_point["predicted_energy_minus_c_eV_per_BTO"],
        "audited_energy_minus_C_eV_per_BTO": energy,
        "audited_minus_predicted_energy_meV_per_BTO": energy_error,
        "full_atom_plus_strain_metric_coordinate_error_sqrt_amu_A": coordinate_error,
        "audited_orthogonal_gradient_eV_per_sqrt_amu_A": gradient,
        "audited_maximum_absolute_stress_kbar": stress,
        "n_raw_DFT_evaluations_audited": audit["n_individually_audited_DFT_points"],
        "checks": checks,
        "all_gates_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--measured-patch", type=Path, required=True)
    parser.add_argument("--host", default="hf")
    parser.add_argument("--remote-root", required=True)
    args = parser.parse_args()
    plan_bytes = args.plan.read_bytes()
    plan = json.loads(plan_bytes)
    if (plan.get("kind") != "predeclared_three_BTO_Qy0_restricted_sheet_interpolation_holdouts_no_DFT"
            or plan.get("status") != "three_unmeasured_holdouts_predicted_before_DFT"
            or plan.get("n_holdouts") != 3 or len(plan.get("holdouts", [])) != 3
            or plan.get("source_measured_patch_sha256") != digest(args.measured_patch.read_bytes())):
        raise ValueError("the frozen holdout plan or measured patch changed")
    reference_bytes = fetch(args.host, args.remote_root, "inputs/cubic_CONTCAR")
    atoms = read(StringIO(reference_bytes.decode("utf-8")), format="vasp")
    if Counter(atoms.get_chemical_symbols()) != Counter({"Ba": 1, "Ti": 1, "O": 3}):
        raise ValueError("reference is not one BaTiO3 formula unit")
    masses = np.asarray(atoms.get_masses(), dtype=float)
    rows = []
    sources = []
    calculator_ids = set()
    evaluator_contracts = set()
    for plan_point, (label, summary_path, audit_path, preflight_path) in zip(plan["holdouts"], CASES):
        payloads = {key: fetch(args.host, args.remote_root, path) for key, path in (
            ("summary", summary_path), ("audit", audit_path), ("preflight", preflight_path))}
        data = {key: json.loads(value) for key, value in payloads.items()}
        rows.append(gate_one(plan_point, data["summary"], data["audit"],
                             data["preflight"], digest(payloads["summary"]),
                             digest(payloads["preflight"]), masses, plan))
        calculator_ids.add(data["preflight"]["calculator_id"])
        evaluator_contracts.add(data["summary"]["evaluator_contract_sha256"])
        sources.append({"label": label, "summary": summary_path, "audit": audit_path,
                        "preflight": preflight_path,
                        "sha256": {key: digest(value) for key, value in payloads.items()}})
    if len(calculator_ids) != 1 or len(evaluator_contracts) != 1:
        raise ValueError("holdouts were not run under one ABACUS calculator contract")
    report = {
        "kind": "BTO_Qy0_three_preregistered_holdouts_raw_audited_interpolation_gate",
        "status": "interpolation_gate_failed" if not all(row["all_gates_pass"] for row in rows)
                  else "local_interpolation_gates_passed_not_global_PES_certificate",
        "frozen_plan_sha256": digest(plan_bytes),
        "source_measured_patch_sha256": plan["source_measured_patch_sha256"],
        "reference_sha256": digest(reference_bytes),
        "gate_script_sha256": digest(Path(__file__).read_bytes()),
        "atomic_mass_order_amu": masses.tolist(),
        "calculator_id": next(iter(calculator_ids)),
        "evaluator_contract_sha256": next(iter(evaluator_contracts)),
        "gates": {key: plan[key] for key in (
            "energy_absolute_error_gate_meV_per_BTO",
            "full_atom_plus_strain_metric_coordinate_error_gate_sqrt_amu_A",
            "gradient_gate_eV_per_sqrt_amu_A",
            "maximum_absolute_stress_gate_kbar")},
        "points": rows,
        "source_files": sources,
        "limitations": "Passing local gates would not certify branch stability, global minima or a complete PES."
    }
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
