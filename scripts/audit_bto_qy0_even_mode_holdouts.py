"""Compare two new raw-audited BTO points to the frozen even-mode predictions.

This reader never submits DFT or edits a calculation. The full ABACUS
SCF/force/stress check is delegated to the independent per-point raw auditor.
"""

from __future__ import annotations

import argparse
from collections import Counter
from io import StringIO
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_bto_qy0_preregistered_holdouts import digest, fetch
from scripts.plan_bto_qy0_even_mode_validation import even_basis, plan as reconstruct_plan


CASES = (
    ("q045_q025", "run-model-holdout-q045_q025/q045_q025_model_holdout_result.json",
     "audit-model-holdout-q045_q025-27795230_0.json",
     "inputs/bto_qy0_q045_q025_model_preflight_2026-09-28.json"),
    ("q105_q025", "run-model-holdout-q105_q025/q105_q025_model_holdout_result.json",
     "audit-model-holdout-q105_q025-27795230_1.json",
     "inputs/bto_qy0_q105_q025_model_preflight_2026-09-28.json"),
)


def evaluate(frozen: dict, summary: dict, audit: dict, preflight: dict,
             summary_hash: str, preflight_hash: str, masses: np.ndarray,
             plan: dict) -> dict:
    q = np.asarray(frozen["q_sqrt_amu_A"], dtype=float)
    if (summary.get("status") != "orthogonal_gradient_and_stress_converged_curvature_unchecked"
            or summary.get("kind")
            != "bto_symmetry_restricted_qy_zero_variable_cell_local_candidate_not_PES_or_barrier"
            or audit.get("status") != "verified_gradient_stationary_candidate_curvature_and_branches_unchecked"
            or audit.get("third_soft_mode_restricted_at_zero") is not True
            or abs(audit.get("final_third_soft_y_amplitude_sqrt_amu_A", np.inf)) > 1e-8
            or audit.get("n_individually_audited_DFT_points", 0) < 1
            or audit.get("n_individually_audited_DFT_points") != len(audit.get("evaluations", []))
            or summary.get("preflight_sha256") != preflight_hash
            or audit.get("source_sha256", {}).get("summary") != summary_hash
            or audit.get("source_sha256", {}).get("preflight") != preflight_hash
            or summary.get("phonon_supercell") != [1, 1, 1]
            or summary.get("electronic_kpoints") != [4, 4, 4]
            or summary.get("stress_target_kbar") != 2.0
            or summary.get("stress_target_passed") is not True
            or preflight.get("no_dft_launched") is not True
            or not np.allclose(summary.get("q_parallel_q_transverse_sqrt_amu_A"), q,
                               atol=1e-12, rtol=0)
            or not np.allclose(audit.get("q_parallel_q_transverse_sqrt_amu_A"), q,
                               atol=1e-12, rtol=0)
            or not np.allclose(preflight.get("q_parallel_q_transverse_sqrt_amu_A"), q,
                               atol=1e-12, rtol=0)):
        raise ValueError(f"model holdout raw-audit contract mismatch at {q.tolist()}")
    energy = float(summary["energy_minus_c_eV_per_BTO"])
    gradient = float(summary["orthogonal_gradient_norm_eV_per_sqrt_amu_A"])
    stress = float(audit["final_maximum_absolute_stress_kbar"])
    if (abs(energy - audit["final_energy_minus_c_eV_per_BTO"]) > 1e-8
            or abs(gradient - audit["orthogonal_gradient_norm_eV_per_sqrt_amu_A"]) > 1e-8
            or abs(stress - summary["maximum_absolute_stress_kbar"]) > 1e-8
            or abs(float(even_basis(q)[0] @ np.asarray(plan["coefficients_eV_per_BTO"]))
                   - frozen["model_prediction_energy_minus_C_eV_per_BTO"]) > 1e-12):
        raise ValueError("model, ABACUS summary, and raw audit disagree")
    observed = np.asarray(summary["coordinates_u_A_eta_voigt"], dtype=float)
    predicted = np.asarray(frozen["predicted_structure_u_A_eta_voigt"], dtype=float)
    strain = np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float)
    weights = np.r_[np.repeat(masses, 3), strain]
    if (observed.shape != (21,) or predicted.shape != (21,) or strain.shape != (6,)
            or np.any(weights <= 0) or not np.all(np.isfinite(observed))):
        raise ValueError("invalid atom-plus-strain metric coordinates")
    metric_error = float(np.sqrt(np.dot(weights, (observed - predicted)**2)))
    energy_error = 1000 * (energy - frozen["model_prediction_energy_minus_C_eV_per_BTO"])
    checks = {
        "energy": abs(energy_error) <= plan["energy_absolute_error_gate_meV_per_BTO"],
        "full_structure": metric_error <= plan["full_structure_seed_error_gate_sqrt_amu_A"],
        "gradient": gradient <= plan["gradient_gate_eV_per_sqrt_amu_A"],
        "stress": stress <= plan["maximum_absolute_stress_gate_kbar"],
    }
    return {
        "q_sqrt_amu_A": q.tolist(),
        "frozen_model_prediction_energy_minus_C_eV_per_BTO": frozen[
            "model_prediction_energy_minus_C_eV_per_BTO"],
        "audited_energy_minus_C_eV_per_BTO": energy,
        "audited_minus_predicted_meV_per_BTO": energy_error,
        "full_structure_seed_error_sqrt_amu_A": metric_error,
        "audited_orthogonal_gradient_eV_per_sqrt_amu_A": gradient,
        "audited_maximum_absolute_stress_kbar": stress,
        "n_raw_DFT_evaluations_audited": audit["n_individually_audited_DFT_points"],
        "checks": checks,
        "all_gates_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen-plan", type=Path, required=True)
    parser.add_argument("--measured-patch", type=Path, required=True)
    parser.add_argument("--earlier-holdout-audit", type=Path, required=True)
    parser.add_argument("--force-constants", type=Path, required=True)
    parser.add_argument("--host", default="hf")
    parser.add_argument("--remote-root", required=True)
    args = parser.parse_args()
    plan_bytes = args.frozen_plan.read_bytes()
    frozen = json.loads(plan_bytes)
    reconstructed = reconstruct_plan(
        json.loads(args.measured_patch.read_text(encoding="utf-8")),
        json.loads(args.earlier_holdout_audit.read_text(encoding="utf-8")),
    )
    if (frozen.get("status") != "two_model_predictions_frozen_before_new_DFT"
            or len(frozen.get("new_holdouts", [])) != 2
            or any(frozen.get(key) != value for key, value in reconstructed.items())
            or frozen.get("source_sha256", {}).get("measured")
            != digest(args.measured_patch.read_bytes())
            or frozen.get("source_sha256", {}).get("linear_holdout_audit")
            != digest(args.earlier_holdout_audit.read_bytes())):
        raise ValueError("frozen model or its training sources changed")
    reference_bytes = fetch(args.host, args.remote_root, "inputs/cubic_CONTCAR")
    atoms = read(StringIO(reference_bytes.decode("utf-8")), format="vasp")
    if Counter(atoms.get_chemical_symbols()) != Counter({"Ba": 1, "Ti": 1, "O": 3}):
        raise ValueError("reference is not one BaTiO3 formula unit")
    force_constants_hash = digest(args.force_constants.read_bytes())
    with np.load(args.force_constants) as phonopy_source:
        masses = np.asarray(phonopy_source["masses_amu"], dtype=float)
    if (masses.shape != (5,) or np.any(masses <= 0)
            or not np.allclose(masses, atoms.get_masses(), rtol=0, atol=0.01)):
        raise ValueError("phonopy masses do not match the BaTiO3 reference")
    rows, sources = [], []
    calculator_ids, evaluator_contracts = set(), set()
    for point, (label, summary_path, audit_path, preflight_path) in zip(frozen["new_holdouts"], CASES):
        files = {key: fetch(args.host, args.remote_root, path) for key, path in (
            ("summary", summary_path), ("audit", audit_path), ("preflight", preflight_path))}
        data = {key: json.loads(value) for key, value in files.items()}
        if (data["preflight"].get("source_sha256", {}).get("force_constants")
                != force_constants_hash):
            raise ValueError("preflight used a different phonopy mass source")
        rows.append(evaluate(point, data["summary"], data["audit"], data["preflight"],
                             digest(files["summary"]), digest(files["preflight"]),
                             masses, frozen))
        calculator_ids.add(data["preflight"]["calculator_id"])
        evaluator_contracts.add(data["summary"]["evaluator_contract_sha256"])
        sources.append({"label": label, "summary": summary_path, "audit": audit_path,
                        "preflight": preflight_path,
                        "sha256": {key: digest(value) for key, value in files.items()}})
    if len(calculator_ids) != 1 or len(evaluator_contracts) != 1:
        raise ValueError("the two new points used different calculator contracts")
    report = {
        "kind": "BTO_Qy0_even_mode_two_prospective_holdout_raw_audit_gate",
        "status": ("two_prospective_local_gates_passed_not_global_PES_certificate"
                   if all(row["all_gates_pass"] for row in rows)
                   else "prospective_model_gate_failed"),
        "frozen_plan_sha256": digest(plan_bytes),
        "reference_sha256": digest(reference_bytes),
        "phonopy_force_constants_sha256": force_constants_hash,
        "metric_atomic_masses_amu": masses.tolist(),
        "gate_script_sha256": digest(Path(__file__).read_bytes()),
        "calculator_id": next(iter(calculator_ids)),
        "evaluator_contract_sha256": next(iter(evaluator_contracts)),
        "gates": {key: frozen[key] for key in (
            "energy_absolute_error_gate_meV_per_BTO",
            "full_structure_seed_error_gate_sqrt_amu_A",
            "gradient_gate_eV_per_sqrt_amu_A",
            "maximum_absolute_stress_gate_kbar")},
        "points": rows,
        "source_files": sources,
        "limitations": "Two pointwise holdouts do not certify positive restricted curvature, competing branches, global surface accuracy, or a physical barrier.",
    }
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
