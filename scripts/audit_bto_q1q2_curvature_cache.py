"""Reconstruct a BTO fixed-Q orthogonal Hessian from audited cached DFT probes.

Read-only: no calculator is constructed. The 32 signed probes must match the
declared metric directions and a separately audited set of SCF/input logs.
Negative eigenvalues reject a conditional *minimum* at this Q point; they do
not by themselves identify a global lower branch or a T-to-C barrier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.bto_q1q2_reference import bto_transverse_soft_plane, load_bto_q1q2_reference, sha256
from vcneb.mode_surface import _orthogonal_directions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("workdir", "preflight", "report", "reference", "force_constants",
                 "phonopy_eigenpairs", "refined_result", "curvature_result", "all_points_audit", "output"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    parser.add_argument("--branch-replay", type=Path,
                        help="required audited-cache branch replay for the transverse-soft plane")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite a curvature audit: {args.output}")
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    refined = json.loads(args.refined_result.read_text(encoding="utf-8"))
    curvature = json.loads(args.curvature_result.read_text(encoding="utf-8"))
    all_points = json.loads(args.all_points_audit.read_text(encoding="utf-8"))
    physical_soft_plane = preflight.get("kind") == "bto_transverse_soft_conditional_preflight_no_dft"
    if (preflight.get("kind") not in {"bto_q1q2_conditional_preflight_no_dft",
                                       "bto_transverse_soft_conditional_preflight_no_dft"}
            or refined.get("status") != "orthogonal_gradient_and_stress_converged_curvature_unchecked"
            or curvature.get("status") != "orthogonal_curvature_screen_complete_not_branch_certification"
            or curvature.get("start_from_result_sha256") != sha256(args.refined_result)
            or curvature.get("evaluator_contract_sha256") != refined.get("evaluator_contract_sha256")
            or all_points.get("source_sha256", {}).get("summary") != sha256(args.refined_result)
            or all_points.get("status") not in {
                "verified_gradient_stationary_candidate_curvature_and_branches_unchecked",
                "verified_gradient_stationary_branch_candidate_local_curvature_unchecked",
            }):
        raise ValueError("BTO preflight/refined/curvature/SCF audit chain differs")
    if physical_soft_plane:
        if args.branch_replay is None:
            raise ValueError("transverse-soft curvature requires the audited three-branch replay")
        replay = json.loads(args.branch_replay.read_text(encoding="utf-8"))
        selected = refined["selected_start"]
        if (replay.get("status") != "all_three_branch_outcomes_reproduced_from_raw_audited_cache"
                or replay.get("source_sha256", {}).get("summary") != sha256(args.refined_result)
                or replay.get("selected_start") != selected
                or replay["branch_outcomes"][selected]["final_evaluation_directory"]
                != refined["final_evaluation_directory"]):
            raise ValueError("selected soft branch differs from audited replay")
    elif args.branch_replay is not None:
        raise ValueError("branch replay belongs only to the transverse-soft plane")
    for key, path in {
        "report": args.report, "reference": args.reference,
        "force_constants": args.force_constants, "phonopy_eigenpairs": args.phonopy_eigenpairs,
    }.items():
        if sha256(path) != preflight["source_sha256"][key]:
            raise ValueError(f"BTO mode source hash changed: {key}")
    loaded = load_bto_q1q2_reference(
        args.report, args.reference, args.force_constants, args.phonopy_eigenpairs,
        strain_metric_weights_amu_A2=np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float),
    )
    chart = loaded.chart
    plane = (bto_transverse_soft_plane(
        loaded, strain_metric_weights_amu_A2=np.asarray(
            preflight["strain_metric_weights_amu_A2"], dtype=float,
        ),
    ) if physical_soft_plane else loaded.plane)
    center = np.asarray(refined["coordinates_u_A_eta_voigt"], dtype=float)
    q = np.asarray(curvature["q1_q2_sqrt_amu_A"], dtype=float)
    if (center.shape != (chart.coordinate_count,)
            or not np.allclose(plane.project(center), q, atol=1e-8, rtol=0.0)
            or not np.allclose(refined["q1_q2_sqrt_amu_A"], q, atol=1e-12, rtol=0.0)):
        raise ValueError("curvature center or fixed Q differs from refined result")
    directions = _orthogonal_directions(plane, chart.rigid_translation_directions())
    count = directions.shape[1]
    step = float(curvature["curvature_step_sqrt_amu_A"])
    if count != 16 or not np.isfinite(step) or step <= 0 or curvature["n_gradient_evaluations"] != 2 * count:
        raise ValueError("unexpected BTO orthogonal dimension or finite-difference step")
    records = {}
    audited_names = {item["directory"] for item in all_points["evaluations"]}
    audited_count = all_points.get("n_individually_audited_DFT_points")
    if (type(audited_count) is not int or audited_count < 2 * count + 1
            or len(audited_names) != audited_count
            or len(all_points["evaluations"]) != audited_count):
        raise ValueError("SCF/input audit has missing or duplicate DFT point records")
    for directory in args.workdir.glob("eval-*-*"):
        result_path = directory / "result.json"
        if not result_path.is_file() or directory.name not in audited_names:
            raise ValueError(f"curvature cache contains an unaudited or incomplete point: {directory}")
        record = json.loads(result_path.read_text(encoding="utf-8"))
        digest = hashlib.sha256(np.asarray(record["coordinates"], dtype=float).tobytes()).hexdigest()
        if (record.get("contract_sha256") != refined["evaluator_contract_sha256"]
                or record.get("coordinate_sha256") != digest or digest in records):
            raise ValueError(f"corrupt or duplicate curvature cache point: {directory}")
        records[digest] = directory.name, record
    if len(records) != audited_count:
        raise ValueError("curvature cache differs from the independently audited DFT point set")
    center_digest = hashlib.sha256(center.tobytes()).hexdigest()
    if center_digest not in records or records[center_digest][0] != refined["final_evaluation_directory"]:
        raise ValueError("curvature center is not the audited selected DFT point")
    center_record = records[center_digest][1]
    center_energy = float(center_record["enthalpy_eV"])
    center_gradient = np.asarray(center_record["gradient"], dtype=float)
    if center_gradient.shape != center.shape or not np.all(np.isfinite(center_gradient)):
        raise ValueError("curvature center has no finite raw-audited gradient")
    hessian = np.empty((count, count), dtype=float)
    probe_names = []
    probe_energies = []
    for index in range(count):
        projected = []
        paired_energy = []
        for sign in (1.0, -1.0):
            expected = center + sign * step * directions[:, index]
            digest = hashlib.sha256(expected.tobytes()).hexdigest()
            if digest not in records:
                raise ValueError(f"missing exact signed DFT probe for orthogonal direction {index}, sign {sign}")
            name, record = records[digest]
            if (record["validation"]["slurm_job_id"] != str(curvature["slurm_job_id"])
                    or not np.array_equal(np.asarray(record["coordinates"]), expected)
                    or not np.allclose(plane.project(expected), q, atol=1e-8, rtol=0.0)):
                raise ValueError(f"curvature probe does not match its declared fixed-Q job: {name}")
            gradient = np.asarray(record["gradient"], dtype=float)
            if gradient.shape != center.shape or not np.all(np.isfinite(gradient)):
                raise ValueError(f"curvature gradient is invalid: {name}")
            projected.append(directions.T @ gradient)
            paired_energy.append(float(record["enthalpy_eV"]))
            probe_names.append(name)
        hessian[:, index] = (projected[0] - projected[1]) / (2.0 * step)
        probe_energies.append(paired_energy)
    if len(set(probe_names)) != 2 * count:
        raise ValueError("a DFT probe was reused for two curvature directions")
    antisymmetric_norm = float(np.linalg.norm(hessian - hessian.T))
    antisymmetric_relative = antisymmetric_norm / max(float(np.linalg.norm(hessian)), 1e-30)
    eigenvalues, eigenvectors = np.linalg.eigh(0.5 * (hessian + hessian.T))
    if not np.allclose(eigenvalues, curvature["orthogonal_hessian_eigenvalues_eV_per_amu_A2"],
                       rtol=0.0, atol=1e-9):
        raise ValueError("reconstructed eigenvalues disagree with the production curvature result")
    negative = []
    for index, eigenvalue in enumerate(eigenvalues):
        if eigenvalue >= 0:
            continue
        full_direction = directions @ eigenvectors[:, index]
        modal = plane.modal_amplitudes(full_direction)
        negative.append({
            "eigenvalue_eV_per_amu_A2": float(eigenvalue),
            "metric_normalized_chart_direction": full_direction.tolist(),
            "soft_gamma_triplet_fraction": float(np.sum(modal[:3] ** 2)),
            "second_gamma_triplet_fraction": float(np.sum(modal[6:9] ** 2)),
            "strain_fraction": float(np.sum(modal[-6:] ** 2)),
            "top_source_mode_indices": np.argsort(np.abs(modal))[::-1][:5].tolist(),
        })
    paired = np.asarray(probe_energies, dtype=float)
    energy_gradient = (paired[:, 0] - paired[:, 1]) / (2 * step)
    gradient_from_forces = directions.T @ center_gradient
    energy_diagonal_curvature = (paired[:, 0] + paired[:, 1] - 2 * center_energy) / (step * step)
    diagonal_from_gradient = np.diag(hessian)
    output = {
        "kind": ("bto_transverse_soft_curvature_reconstruction_from_audited_DFT_not_PES_or_barrier"
                 if physical_soft_plane else
                 "bto_fixed_q1q2_curvature_reconstruction_from_audited_DFT_not_PES_or_barrier"),
        "status": "negative_orthogonal_curvature_rejects_conditional_local_minimum"
                  if negative else "no_negative_curvature_at_one_difference_step_only",
        "q1_q2_sqrt_amu_A": q.tolist(),
        "axis_kind": preflight.get("axis_kind", "archived_soft_stable"),
        "step_sqrt_amu_A": step,
        "n_audited_DFT_points": len(records),
        "n_matched_signed_probes": len(probe_names),
        "eigenvalues_eV_per_amu_A2": eigenvalues.tolist(),
        "lowest_eigenvector_in_common_orthogonal_basis": eigenvectors[:, 0].tolist(),
        "negative_directions": negative,
        "hessian_antisymmetric_frobenius_norm_eV_per_amu_A2": antisymmetric_norm,
        "hessian_antisymmetric_relative_defect": antisymmetric_relative,
        "energy_gradient_max_abs_difference_eV_per_sqrt_amu_A": float(
            np.max(np.abs(energy_gradient - gradient_from_forces))
        ),
        "energy_hessian_diagonal_max_abs_difference_eV_per_amu_A2": float(
            np.max(np.abs(energy_diagonal_curvature - diagonal_from_gradient))
        ),
        "energy_gradient_per_open_direction_eV_per_sqrt_amu_A": energy_gradient.tolist(),
        "force_gradient_per_open_direction_eV_per_sqrt_amu_A": gradient_from_forces.tolist(),
        "energy_diagonal_curvature_eV_per_amu_A2": energy_diagonal_curvature.tolist(),
        "probe_evaluation_directories": probe_names,
        "probe_pair_energies_eV": probe_energies,
        "source_sha256": {
            "preflight": sha256(args.preflight), "refined_result": sha256(args.refined_result),
            "curvature_result": sha256(args.curvature_result),
            "all_points_audit": sha256(args.all_points_audit),
            "branch_replay": sha256(args.branch_replay) if physical_soft_plane else None,
        },
        "limitations": "one finite-difference step; a second step and competing-branch checks remain necessary",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": output["status"], "n_audited_DFT_points": len(records),
        "negative_eigenvalues_eV_per_amu_A2": [v["eigenvalue_eV_per_amu_A2"] for v in negative],
        "negative_soft_triplet_fractions": [v["soft_gamma_triplet_fraction"] for v in negative],
    }, indent=2))


if __name__ == "__main__":
    main()
