"""Reconstruct the BTO fixed-Q soft mixed direction from audited DFT probes.

Read-only until the final new JSON report. Neither this preflight nor positive
curvature along sampled basis axes certifies a conditional local minimum.
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

from examples.bto_q1q2_reference import (
    bto_transverse_soft_plane,
    load_bto_q1q2_reference,
    sha256,
)
from vcneb.mode_surface import _orthogonal_directions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("preflight", "report", "reference", "force_constants", "phonopy_eigenpairs",
                 "branch_result", "workdir", "curvature_small", "curvature_large",
                 "result_small", "result_large", "audit_small", "audit_large", "output"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    parser.add_argument("--branch-replay", type=Path,
                        help="required cache-only three-branch replay for a transverse-soft plane")
    return parser.parse_args()


def _reconstruct(
    *, reconstruction: dict, curvature_result: dict, audit: dict,
    reconstruction_path: Path, curvature_path: Path, audit_path: Path,
    workdir: Path, center: np.ndarray, q: np.ndarray, directions: np.ndarray,
    evaluator_contract: str, plane, expected_kind: str,
) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    step = float(reconstruction["step_sqrt_amu_A"])
    n = directions.shape[1]
    if (reconstruction.get("kind") != expected_kind
            or reconstruction.get("status") != "no_negative_curvature_at_one_difference_step_only"
            or reconstruction.get("source_sha256", {}).get("curvature_result") != sha256(curvature_path)
            or reconstruction.get("source_sha256", {}).get("all_points_audit") != sha256(audit_path)
            or reconstruction.get("n_matched_signed_probes") != 2 * n
            or curvature_result.get("curvature_step_sqrt_amu_A") != step
            or audit.get("n_individually_audited_DFT_points") != reconstruction.get("n_audited_DFT_points")
            or not np.allclose(reconstruction.get("q1_q2_sqrt_amu_A"), q, rtol=0.0, atol=1e-12)
            or not np.isfinite(step) or step <= 0):
        raise ValueError(f"curvature reconstruction/SCF audit chain differs: {reconstruction_path}")
    names = reconstruction["probe_evaluation_directories"]
    audited = {item["directory"] for item in audit["evaluations"]}
    if len(names) != 2 * n or len(set(names)) != len(names) or not set(names) <= audited:
        raise ValueError("curvature probes are missing or duplicate in the independent SCF audit")
    hessian = np.empty((n, n))
    records: list[dict] = []
    for axis in range(n):
        projected = []
        for offset, sign in enumerate((1.0, -1.0)):
            name = names[2 * axis + offset]
            path = workdir / name / "result.json"
            record = json.loads(path.read_text(encoding="utf-8"))
            expected = center + sign * step * directions[:, axis]
            coordinates = np.asarray(record["coordinates"], dtype=float)
            gradient = np.asarray(record["gradient"], dtype=float)
            digest = hashlib.sha256(expected.tobytes()).hexdigest()
            if (coordinates.shape != center.shape or gradient.shape != center.shape
                    or not np.array_equal(coordinates, expected)
                    or record.get("coordinate_sha256") != digest
                    or record.get("contract_sha256") != evaluator_contract
                    or record.get("validation", {}).get("slurm_job_id") != str(curvature_result["slurm_job_id"])
                    or not np.all(np.isfinite(gradient))
                    or not np.isfinite(float(record["enthalpy_eV"]))
                    or not np.allclose(plane.project(coordinates), q, rtol=0.0, atol=1e-8)):
                raise ValueError(f"signed DFT probe does not match fixed-Q contract: {path}")
            projected.append(directions.T @ gradient)
            records.append({"directory": name, "result_sha256": sha256(path)})
        hessian[:, axis] = (projected[0] - projected[1]) / (2.0 * step)
    eigenvalues, eigenvectors = np.linalg.eigh((hessian + hessian.T) * 0.5)
    if not np.allclose(eigenvalues, reconstruction["eigenvalues_eV_per_amu_A2"], rtol=0.0, atol=1e-9):
        raise ValueError("fresh cache reconstruction differs from the independently audited eigenvalues")
    full_direction = directions @ eigenvectors[:, 0]
    return eigenvalues, full_direction, records


def main() -> None:
    args = parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite soft-direction preflight: {args.output}")
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    branch = json.loads(args.branch_result.read_text(encoding="utf-8"))
    physical_soft_plane = preflight.get("kind") == "bto_transverse_soft_conditional_preflight_no_dft"
    if (preflight.get("kind") not in {
            "bto_q1q2_conditional_preflight_no_dft",
            "bto_transverse_soft_conditional_preflight_no_dft",
        }
            or branch.get("kind") != (
                "bto_transverse_soft_variable_cell_conditional_local_candidate_not_PES_or_barrier"
                if physical_soft_plane else
                "bto_fixed_q1q2_variable_cell_conditional_local_candidate_not_T_to_C_barrier")
            or branch.get("status") != "orthogonal_gradient_and_stress_converged_curvature_unchecked"
            or branch.get("preflight_sha256") != sha256(args.preflight)
            or not branch.get("stress_target_passed")):
        raise ValueError("branch center does not match the audited fixed-Q preflight")
    if physical_soft_plane:
        if args.branch_replay is None:
            raise ValueError("transverse-soft direct probes require the audited three-branch replay")
        replay = json.loads(args.branch_replay.read_text(encoding="utf-8"))
        selected = branch.get("selected_start")
        if (replay.get("status") != "all_three_branch_outcomes_reproduced_from_raw_audited_cache"
                or replay.get("source_sha256", {}).get("summary") != sha256(args.branch_result)
                or replay.get("selected_start") != selected
                or replay["branch_outcomes"][selected]["final_evaluation_directory"]
                != branch.get("final_evaluation_directory")):
            raise ValueError("selected soft branch differs from audited replay")
    elif args.branch_replay is not None:
        raise ValueError("branch replay belongs only to the transverse-soft plane")
    loaded = load_bto_q1q2_reference(
        args.report, args.reference, args.force_constants, args.phonopy_eigenpairs,
        strain_metric_weights_amu_A2=np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float),
    )
    if any(preflight.get("source_sha256", {}).get(key) != digest
           for key, digest in loaded.source_hashes.items()):
        raise ValueError("mode/reference sources differ from the audited preflight")
    chart = loaded.chart
    plane = (bto_transverse_soft_plane(
        loaded,
        strain_metric_weights_amu_A2=np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float),
    ) if physical_soft_plane else loaded.plane)
    center = np.asarray(branch["coordinates_u_A_eta_voigt"], dtype=float)
    q = np.asarray(branch["q1_q2_sqrt_amu_A"], dtype=float)
    if (center.shape != (chart.coordinate_count,) or not np.all(np.isfinite(center))
            or not np.allclose(plane.project(center), q, rtol=0.0, atol=1e-8)):
        raise ValueError("branch center is not a finite point at the declared fixed Q")
    directions = _orthogonal_directions(plane, chart.rigid_translation_directions())
    if directions.shape[1] != 16:
        raise ValueError("BTO fixed-Q open dimension changed from the audited 16")
    reconstructed = []
    for recon_path, result_path, audit_path in (
        (args.curvature_small, args.result_small, args.audit_small),
        (args.curvature_large, args.result_large, args.audit_large),
    ):
        recon = json.loads(recon_path.read_text(encoding="utf-8"))
        result = json.loads(result_path.read_text(encoding="utf-8"))
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if (recon.get("source_sha256", {}).get("preflight") != sha256(args.preflight)
                or recon.get("source_sha256", {}).get("refined_result") != sha256(args.branch_result)
                or result.get("evaluator_contract_sha256") != branch["evaluator_contract_sha256"]):
            raise ValueError("curvature source is not the same BTO branch/calculator contract")
        values, direction, records = _reconstruct(
            reconstruction=recon, curvature_result=result, audit=audit,
            reconstruction_path=recon_path, curvature_path=result_path, audit_path=audit_path,
            workdir=args.workdir, center=center, q=q, directions=directions,
            evaluator_contract=branch["evaluator_contract_sha256"], plane=plane,
            expected_kind=(
                "bto_transverse_soft_curvature_reconstruction_from_audited_DFT_not_PES_or_barrier"
                if physical_soft_plane else
                "bto_fixed_q1q2_curvature_reconstruction_from_audited_DFT_not_PES_or_barrier"
            ),
        )
        reconstructed.append((float(recon["step_sqrt_amu_A"]), values, direction, records))
    first, second = reconstructed
    if not first[0] < second[0]:
        raise ValueError("curvature inputs must be ordered by increasing step")
    overlap = float(abs(first[2] @ (plane.metric_weights * second[2])))
    if overlap < 0.95:
        raise ValueError("the soft eigenvector rotates substantially between finite-difference steps")
    chosen = first[2] * (1.0 if first[2][np.argmax(np.abs(first[2]))] > 0 else -1.0)
    norm = float(np.sqrt(chosen @ (plane.metric_weights * chosen)))
    if abs(norm - 1.0) > 1e-8 or not np.allclose(plane.project(center + chosen), q, rtol=0.0, atol=1e-8):
        raise ValueError("reconstructed soft direction does not preserve fixed Q or metric normalization")
    minimum_allowed = float(preflight["minimum_allowed_atomic_distance_A"])
    probes = []
    for step in (first[0], second[0]):
        for sign in (1, -1):
            coordinates = center + sign * step * chosen
            atoms = chart.to_atoms(coordinates)
            distances = atoms.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            minimum = float(np.min(distances))
            volume = float(atoms.get_volume())
            if minimum < minimum_allowed or volume <= 0 or not np.isfinite(volume):
                raise ValueError("proposed soft-direction probe fails the pre-DFT geometry guard")
            probes.append({"step_sqrt_amu_A": step, "sign": sign,
                           "coordinates_u_A_eta_voigt": coordinates.tolist(),
                           "minimum_atomic_distance_A": minimum, "volume_A3": volume})
    modal = plane.modal_amplitudes(chosen)
    output = {
        "kind": ("bto_transverse_soft_mixed_eigenvector_probe_preflight_no_dft"
                 if physical_soft_plane else
                 "bto_fixed_q1q2_soft_mixed_eigenvector_probe_preflight_no_dft"),
        "status": "four_soft_direction_geometries_passed_no_dft",
        "q1_q2_sqrt_amu_A": q.tolist(),
        "center_coordinates_u_A_eta_voigt": center.tolist(),
        "soft_direction_metric_unit_chart": chosen.tolist(),
        "soft_direction_mode_amplitudes": modal.tolist(),
        "two_step_soft_eigenvector_overlap_abs": overlap,
        "two_step_minimum_eigenvalues_eV_per_amu_A2": [float(first[1][0]), float(second[1][0])],
        "probes": probes,
        "minimum_allowed_atomic_distance_A": minimum_allowed,
        "cached_probe_result_sha256": {"small": first[3], "large": second[3]},
        "source_sha256": {key: sha256(path) for key, path in {
            "preflight": args.preflight, "report": args.report, "reference": args.reference,
            "force_constants": args.force_constants, "phonopy_eigenpairs": args.phonopy_eigenpairs,
            "branch_result": args.branch_result, "curvature_small": args.curvature_small,
            "curvature_large": args.curvature_large, "result_small": args.result_small,
            "result_large": args.result_large, "audit_small": args.audit_small,
            "audit_large": args.audit_large,
            **({"branch_replay": args.branch_replay} if physical_soft_plane else {}),
        }.items()},
        "limitations": "geometry and eigenvector consistency only; no new DFT energy curvature or local-minimum certification",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: output[key] for key in (
        "status", "two_step_soft_eigenvector_overlap_abs",
        "two_step_minimum_eigenvalues_eV_per_amu_A2",
    )}, indent=2))


if __name__ == "__main__":
    main()
