"""Project the final GaN VCNEB path tangent onto a local joint Hessian mode.

This is a local geometric diagnostic, not a transition-state or basin-link
certificate. The final chain, static center and Hessian are hash-checked.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.units import GPa

from vcneb.joint_curvature import JointCurvatureCoordinates, symmetric_strain_basis


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def joint_delta(reference: Atoms, target: Atoms, cell_scale_A: float) -> tuple[np.ndarray, float]:
    """Return local atomic/symmetric-strain coordinates and rotation residual.

    Fractional differences are placed in the nearest periodic image. Thus this
    chart is intended for neighboring images, not for a whole reconstructive
    path with potentially different atom mapping or lattice gauge.
    """

    if (reference.get_chemical_symbols() != target.get_chemical_symbols()
            or len(reference) != len(target)):
        raise ValueError("target atoms/order differ from the reference")
    cell0 = np.asarray(reference.cell.array, dtype=float)
    cell = np.asarray(target.cell.array, dtype=float)
    if (not np.isfinite(cell).all() or np.linalg.det(cell) <= 0
            or not np.isfinite(cell_scale_A) or cell_scale_A <= 0):
        raise ValueError("invalid cell or scale")
    dq = (target.get_scaled_positions(wrap=False)
          - reference.get_scaled_positions(wrap=False))
    dq -= np.rint(dq)
    atomic = dq @ cell0
    deformation = (np.linalg.solve(cell0, cell)).T
    symmetric = 0.5 * (deformation + deformation.T) - np.eye(3)
    rotation_residual = float(np.linalg.norm(
        0.5 * (deformation - deformation.T), ord="fro",
    ))
    strain = cell_scale_A * np.einsum("ij,aij->a", symmetric, symmetric_strain_basis())
    return np.concatenate((atomic.ravel(), strain)), rotation_residual


def analyze(
    trajectory: Path, source_summary: Path, center_poscar: Path,
    hessian_manifest: Path, hessian_audit: Path, hessian_npz: Path,
) -> dict:
    summary = json.loads(source_summary.read_text(encoding="utf-8"))
    manifest = json.loads(hessian_manifest.read_text(encoding="utf-8"))
    audit = json.loads(hessian_audit.read_text(encoding="utf-8"))
    if (summary.get("status") != "completed"
            or not summary.get("converged")
            or manifest.get("source_sha256", {}).get("center_POSCAR") != sha256(center_poscar)
            or audit.get("source_sha256", {}).get("manifest") != sha256(hessian_manifest)
            or audit.get("n_static_displacements") != 36
            or manifest.get("phonon_supercell_matrix") != np.eye(3, dtype=int).tolist()
            or audit.get("step_A") != manifest.get("step_A")):
        raise ValueError("source chain or audited 1x1x1 Hessian provenance changed")
    n_images = int(summary["n_images"])
    peak = int(summary["saddle_diagnostics"]["image_index"])
    if n_images < 3 or not 0 < peak < n_images - 1:
        raise ValueError("final chain has no interior peak for central tangent")
    with Trajectory(trajectory) as traj:
        if len(traj) < n_images or len(traj) % n_images:
            raise ValueError("trajectory does not contain whole VCNEB chains")
        frames = [traj[index] for index in range(len(traj) - n_images, len(traj))]
        n_frames = len(traj)
    pressure = float(audit["pressure_GPa"]) * GPa
    energies = np.asarray([
        frame.get_potential_energy() + pressure * frame.get_volume()
        for frame in frames
    ])
    if not np.allclose(energies, summary["image_enthalpies_eV"], atol=1e-7, rtol=0):
        raise ValueError("final trajectory does not reproduce source image enthalpies")
    if np.argmax(energies) != peak:
        raise ValueError("declared peak differs from final chain")
    center = read(center_poscar, format="vasp")
    scale = float(manifest["cell_scale_A"])
    coordinates = JointCurvatureCoordinates(center, scale)
    basis = coordinates.translation_free_basis()
    with np.load(hessian_npz) as npz:
        hessian = np.asarray(npz["hessian"], dtype=float)
        eigenvalues = np.asarray(npz["eigenvalues"], dtype=float)
        eigenvectors = np.asarray(npz["eigenvectors"], dtype=float)
    if (hessian.shape != (18, 18) or eigenvalues.shape != (15,)
            or eigenvectors.shape != (18, 15)
            or not np.allclose(eigenvalues, audit["translation_free_eigenvalues_eV_per_A2"], atol=1e-8, rtol=0)
            or not np.allclose(eigenvectors.T @ eigenvectors, np.eye(15), atol=1e-8, rtol=0)
            or not np.allclose(
                basis @ (basis.T @ hessian @ eigenvectors),
                eigenvectors * eigenvalues, atol=1e-7, rtol=0,
            )
            or not eigenvalues[0] < 0 < eigenvalues[1]):
        raise ValueError("audited Hessian eigensystem differs from archive")
    left, left_rotation = joint_delta(center, frames[peak - 1], scale)
    middle, middle_rotation = joint_delta(center, frames[peak], scale)
    right, right_rotation = joint_delta(center, frames[peak + 1], scale)
    chord = basis @ (basis.T @ (right - left))
    norm = float(np.linalg.norm(chord))
    if norm <= 1e-10:
        raise ValueError("neighboring images have no resolved joint-coordinate chord")
    tangent = chord / norm
    unstable = eigenvectors[:, 0]
    overlap = float(abs(np.dot(tangent, unstable)))
    return {
        "status": "GaN_local_path_tangent_projection_not_TS_or_basin_certificate",
        "n_trajectory_frames": n_frames,
        "n_images_final_chain": n_images,
        "peak_image_index": peak,
        "pressure_GPa": float(audit["pressure_GPa"]),
        "hessian_step_A": float(audit["step_A"]),
        "lowest_eigenvalue_eV_per_A2": float(eigenvalues[0]),
        "next_eigenvalue_eV_per_A2": float(eigenvalues[1]),
        "absolute_unstable_mode_path_tangent_overlap": overlap,
        "unstable_mode_atomic_squared_fraction": float(np.dot(unstable[:12], unstable[:12])),
        "unstable_mode_strain_squared_fraction": float(np.dot(unstable[12:], unstable[12:])),
        "path_tangent_atomic_squared_fraction": float(np.dot(tangent[:12], tangent[:12])),
        "path_tangent_strain_squared_fraction": float(np.dot(tangent[12:], tangent[12:])),
        "peak_to_refined_center_joint_distance_A": float(np.linalg.norm(basis.T @ middle)),
        "neighboring_chord_joint_length_A": norm,
        "near_peak_cell_rotation_residual_frobenius": {
            "left": left_rotation, "peak": middle_rotation, "right": right_rotation,
        },
        "source_sha256": {
            "trajectory": sha256(trajectory), "summary": sha256(source_summary),
            "center_POSCAR": sha256(center_poscar), "hessian_manifest": sha256(hessian_manifest),
            "hessian_audit": sha256(hessian_audit), "hessian_npz": sha256(hessian_npz),
        },
        "limitations": [
            "Only the final chain's peak-neighbor central chord is projected.",
            "This chord is a local tangent estimate, not a soft-mode eigenvector or basin link.",
            "The near-stationary candidate is not certified until step-size and both basin links pass.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trajectory", "source-summary", "center-poscar", "hessian-manifest",
                 "hessian-audit", "hessian-npz", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = analyze(args.trajectory, args.source_summary, args.center_poscar,
                     args.hessian_manifest, args.hessian_audit, args.hessian_npz)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "overlap": result["absolute_unstable_mode_path_tangent_overlap"],
        "rotation_residuals": result["near_peak_cell_rotation_residual_frobenius"],
    }))


if __name__ == "__main__":
    main()
