"""Join two independently raw-audited BTO fixed-Q curvature screens.

This is a read-only evidence gate. Reproducible positive force-Hessian
eigenvalues do not by themselves establish a conditional minimum when the
energy/force integrability check is of comparable size.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def join_two_steps(folder: Path) -> dict:
    results = []
    reconstructions = []
    for label, step in (("0p05", 0.05), ("0p10", 0.10)):
        result_path = folder / f"curvature_{label}_result.json"
        audit_path = folder / f"curvature_{label}_raw_audit.json"
        reconstruction_path = folder / f"curvature_{label}_reconstruction.json"
        result = json.loads(result_path.read_text(encoding="utf-8"))
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        reconstruction = json.loads(reconstruction_path.read_text(encoding="utf-8"))
        sources = reconstruction.get("source_sha256", {})
        if not (
            result.get("status") == "orthogonal_curvature_screen_complete_not_branch_certification"
            and audit.get("status") == "verified_gradient_stationary_candidate_curvature_and_branches_unchecked"
            and reconstruction.get("status") == "no_negative_curvature_at_one_difference_step_only"
            and sources.get("curvature_result") == digest(result_path)
            and sources.get("all_points_audit") == digest(audit_path)
            and result.get("curvature_step_sqrt_amu_A") == step
            and reconstruction.get("step_sqrt_amu_A") == step
            and reconstruction.get("n_matched_signed_probes") == 32
            and len(reconstruction.get("eigenvalues_eV_per_amu_A2", [])) == 16
            and len(reconstruction.get("lowest_eigenvector_in_common_orthogonal_basis", [])) == 16
            and np.allclose(result.get("q1_q2_sqrt_amu_A"), (0.75, 0.15), atol=1e-12, rtol=0)
            and result.get("q1_q2_sqrt_amu_A") == audit.get("q1_q2_sqrt_amu_A")
            == reconstruction.get("q1_q2_sqrt_amu_A")
            and np.allclose(result["orthogonal_hessian_eigenvalues_eV_per_amu_A2"],
                            reconstruction["eigenvalues_eV_per_amu_A2"], atol=1e-12, rtol=0)
            and all(np.isfinite(reconstruction[key]) for key in (
                "hessian_antisymmetric_relative_defect",
                "energy_gradient_max_abs_difference_eV_per_sqrt_amu_A",
                "energy_hessian_diagonal_max_abs_difference_eV_per_amu_A2",
            ))
        ):
            raise ValueError(f"{label}: incomplete or inconsistent raw-audited curvature chain")
        results.append(result)
        reconstructions.append(reconstruction)
    if not (
        results[0]["start_from_result_sha256"] == results[1]["start_from_result_sha256"]
        and results[0]["evaluator_contract_sha256"] == results[1]["evaluator_contract_sha256"]
        and results[0]["runner_sha256"] == results[1]["runner_sha256"]
        and reconstructions[0]["source_sha256"]["refined_result"]
        == reconstructions[1]["source_sha256"]["refined_result"]
        and reconstructions[0]["source_sha256"]["branch_replay"]
        == reconstructions[1]["source_sha256"]["branch_replay"]
    ):
        raise ValueError("the two steps use different center, contract or branch")
    vectors = [np.asarray(r["lowest_eigenvector_in_common_orthogonal_basis"], dtype=float)
               for r in reconstructions]
    if any(not np.isclose(np.linalg.norm(vector), 1, atol=1e-8) for vector in vectors):
        raise ValueError("lowest eigenvectors are not metric-unit vectors")
    eigenvalues = [float(r["eigenvalues_eV_per_amu_A2"][0]) for r in reconstructions]
    overlap = float(abs(vectors[0] @ vectors[1]))
    if min(eigenvalues) <= 0 or overlap < 0.99:
        raise ValueError("the positive soft-direction screen does not reproduce")
    mismatches = [float(r["energy_hessian_diagonal_max_abs_difference_eV_per_amu_A2"])
                  for r in reconstructions]
    return {
        "kind": "bto_q075_q015_two_step_orthogonal_curvature_screen_not_PES",
        "status": "positive_force_hessian_sign_reproduced_energy_force_certification_open",
        "q1_q2_sqrt_amu_A": [0.75, 0.15],
        "steps_sqrt_amu_A": [0.05, 0.10],
        "minimum_eigenvalues_eV_per_amu_A2": eigenvalues,
        "minimum_eigenvalue_step_difference_eV_per_amu_A2": abs(eigenvalues[1] - eigenvalues[0]),
        "absolute_lowest_eigenvector_overlap": overlap,
        "energy_hessian_diagonal_max_abs_differences_eV_per_amu_A2": mismatches,
        "energy_gradient_max_abs_differences_eV_per_sqrt_amu_A": [
            float(r["energy_gradient_max_abs_difference_eV_per_sqrt_amu_A"])
            for r in reconstructions
        ],
        "hessian_antisymmetric_relative_defects": [
            float(r["hessian_antisymmetric_relative_defect"]) for r in reconstructions
        ],
        "n_raw_audited_DFT_points_by_step": [
            int(r["n_audited_DFT_points"]) for r in reconstructions
        ],
        "source_sha256": {
            path.name: digest(path)
            for label in ("0p05", "0p10")
            for path in (
                folder / f"curvature_{label}_result.json",
                folder / f"curvature_{label}_raw_audit.json",
                folder / f"curvature_{label}_reconstruction.json",
            )
        },
        "limitations": "The signed force Hessian is stable across two difference steps, but its weakest eigenvalue is smaller than the maximum 0.05-step energy/force diagonal-curvature disagreement. A direct energy-gradient check along the mixed soft eigenvector and branch/domain validation remain required before claiming a stable conditional surface.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folder", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    report = join_two_steps(args.folder)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "status", "minimum_eigenvalues_eV_per_amu_A2",
        "absolute_lowest_eigenvector_overlap",
        "energy_hessian_diagonal_max_abs_differences_eV_per_amu_A2",
    )}, indent=2))


if __name__ == "__main__":
    main()
