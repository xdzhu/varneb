"""Audit GaN image-15 VASP probes and a translation-free joint Hessian.

The center is a converged ordinary-NEB *maximum image*, not yet a stationary
point of H=E+PV. Therefore eigenvalues here diagnose local curvature only;
even one negative eigenvalue would not certify a variable-cell transition
state without stationarity, step-size stability and basin connections.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.prepare_gan_joint_curvature import same_geometry
from vcneb.core import deformation_from_cell
from vcneb.joint_curvature import (
    JointCurvatureCoordinates, central_difference_hessian, symmetric_strain_basis,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def completed_case(directory: Path, record: dict) -> tuple[object, dict]:
    hashes = json.loads((directory / "sha256.inputs.json").read_text(encoding="utf-8"))
    if hashes["POSCAR"] != record["POSCAR_sha256"]:
        raise ValueError(f"manifest/input hash mismatch: {directory}")
    for name, expected in hashes.items():
        if sha256(directory / name) != expected:
            raise ValueError(f"input changed after staging: {directory / name}")
    outcar = directory / "OUTCAR"
    if not outcar.is_file():
        raise FileNotFoundError(f"missing VASP result: {outcar}")
    output = outcar.read_text(encoding="utf-8", errors="replace")
    if ("General timing and accounting informations" not in output
            or "aborting loop because EDIFF is reached" not in output):
        raise ValueError(f"static VASP SCF/OUTCAR not complete: {outcar}")
    atoms = read(outcar)
    supplied = read(directory / "POSCAR", format="vasp")
    if not same_geometry(atoms, supplied, tolerance=2e-5):
        raise ValueError(f"VASP output geometry differs from input: {directory}")
    energy = float(atoms.get_potential_energy())
    forces = np.asarray(atoms.get_forces(), dtype=float)
    stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
    if (not np.isfinite(energy) or forces.shape != (4, 3) or stress.shape != (3, 3)
            or not np.isfinite(forces).all() or not np.isfinite(stress).all()):
        raise ValueError(f"invalid complete VASP outputs: {directory}")
    return atoms, {
        "case": record["name"], "input_sha256": hashes, "outcar_sha256": sha256(outcar),
        "energy_eV_per_cell": energy,
        "maximum_atomic_force_eV_per_A": float(np.max(np.linalg.norm(forces, axis=1))),
        "max_abs_stress_eV_per_A3": float(np.max(np.abs(stress))),
    }


def tangent_in_joint_coordinates(coordinates: JointCurvatureCoordinates, left, right) -> np.ndarray:
    center = coordinates.reference
    cell0 = center.cell.array
    q_left = left.get_scaled_positions(wrap=False)
    q_right = right.get_scaled_positions(wrap=False)
    delta_q = q_right - q_left
    # The archived chain is unwrapped. A large fractional jump means its
    # periodic gauge must be examined, not silently minimum-imaged here.
    if np.max(np.abs(delta_q)) > 0.5:
        raise ValueError("neighboring VCNEB images have an unexpected gauge jump")
    atomic = (delta_q @ cell0).ravel()
    deformation = (deformation_from_cell(right.cell.array, cell0)
                   - deformation_from_cell(left.cell.array, cell0))
    rotation = float(np.linalg.norm(0.5 * (deformation - deformation.T)))
    strain = coordinates.cell_scale_A * np.einsum(
        "ij,aij->a", deformation, symmetric_strain_basis(),
    )
    tangent = np.concatenate((atomic, strain))
    if not np.isfinite(tangent).all() or np.linalg.norm(tangent) <= 0:
        raise ValueError("invalid local VCNEB tangent")
    return tangent / np.linalg.norm(tangent), rotation


def audit(work_root: Path, trajectory: Path, summary: Path, output_dir: Path) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite previous Hessian audit: {output_dir}")
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("coordinate_count") != 18
            or manifest.get("physical_dimension_after_translation_removal") != 15
            or manifest.get("supercell_matrix") != np.eye(3, dtype=int).tolist()
            or len(manifest.get("cases", [])) != 36
            or manifest["source_sha256"]["trajectory"] != sha256(trajectory)
            or manifest["source_sha256"]["summary"] != sha256(summary)
            or manifest["source_sha256"]["preparer"]
            != sha256(Path(__file__).with_name("prepare_gan_joint_curvature.py"))
            or manifest["source_sha256"]["joint_curvature"]
            != sha256(Path(__file__).resolve().parents[1] / "vcneb/joint_curvature.py")):
        raise ValueError("GaN Hessian stage/source provenance gate failed")
    chain = read(trajectory, index="-29:")
    center = chain[15]
    if not same_geometry(center, read(work_root / "center_POSCAR", format="vasp")):
        raise ValueError("highest image differs from staged center")
    center_outcar = work_root / "center_OUTCAR"
    if sha256(center_outcar) != manifest["source_static_sha256"]["OUTCAR"]:
        raise ValueError("original image-15 center OUTCAR hash differs from manifest")
    original_text = center_outcar.read_text(encoding="utf-8", errors="replace")
    if ("General timing and accounting informations" not in original_text
            or "aborting loop because EDIFF is reached" not in original_text):
        raise ValueError("original image-15 static SCF was incomplete")
    baseline = read(center_outcar)
    if (not same_geometry(center, baseline, tolerance=2e-5)
            or not np.isclose(center.get_potential_energy(), baseline.get_potential_energy(), atol=1e-6)
            or not np.allclose(center.get_forces(), baseline.get_forces(), atol=1e-6)
            or not np.allclose(center.get_stress(voigt=False), baseline.get_stress(voigt=False), atol=1e-7)):
        raise ValueError("archived trajectory center disagrees with original VASP OUTCAR")
    metadata = json.loads(summary.read_text(encoding="utf-8"))
    if metadata.get("status") != "completed" or not metadata.get("converged"):
        raise ValueError("source VCNEB chain no longer passes completion gate")
    coordinates = JointCurvatureCoordinates(center, float(manifest["cell_scale_A"]))
    pressure = float(manifest["pressure_GPa"]) * GPa
    recorded = metadata["path_diagnostics"]["images"]
    if (len(recorded) != 29
            or not np.isclose(metadata["path_diagnostics"]["pressure_eV_per_A3"], pressure,
                              rtol=0, atol=1e-9)
            or any(item["image_index"] != index for index, item in enumerate(recorded))):
        raise ValueError("source VCNEB image/pressure summary disagrees with the staged chain")
    enthalpy_trace = np.array([image.get_potential_energy() + pressure * image.get_volume()
                               for image in chain])
    if (not np.allclose(enthalpy_trace, [item["enthalpy_eV"] for item in recorded],
                        rtol=0, atol=1e-6)
            or int(np.argmax(enthalpy_trace[1:-1])) + 1 != 15):
        raise ValueError("archived 29-image enthalpy trace disagrees with highest-image source")
    step = float(manifest["step_A"])
    expected_names = [f"axis{i:02d}_{side}" for i in range(18) for side in ("plus", "minus")]
    if [case["name"] for case in manifest["cases"]] != expected_names:
        raise ValueError("incomplete or reordered paired finite displacements")
    gradients = np.zeros((18, 2, 18))
    enthalpies = np.zeros((18, 2))
    details = []
    for record in manifest["cases"]:
        directory = work_root / "cases" / record["name"]
        atoms, detail = completed_case(directory, record)
        gradient = coordinates.enthalpy_gradient(
            atoms, atoms.get_forces(), atoms.get_stress(voigt=False), pressure,
        )
        side = 0 if record["sign"] == 1 else 1
        gradients[record["axis"], side] = gradient
        enthalpies[record["axis"], side] = atoms.get_potential_energy() + pressure * atoms.get_volume()
        details.append(detail)
    hessian, reciprocity_defect = central_difference_hessian(
        gradients[:, 0], gradients[:, 1], step,
    )
    translation_basis = coordinates.translation_free_basis()
    reduced = translation_basis.T @ hessian @ translation_basis
    eigenvalues, reduced_vectors = np.linalg.eigh(reduced)
    full_vectors = translation_basis @ reduced_vectors
    tangent, tangent_rotation = tangent_in_joint_coordinates(coordinates, chain[14], chain[16])
    tangent = translation_basis @ (translation_basis.T @ tangent)
    tangent /= np.linalg.norm(tangent)
    center_gradient = coordinates.enthalpy_gradient(
        baseline, baseline.get_forces(), baseline.get_stress(voigt=False), pressure,
    )
    center_enthalpy = baseline.get_potential_energy() + pressure * baseline.get_volume()
    energy_gradient = (enthalpies[:, 0] - enthalpies[:, 1]) / (2.0 * step)
    energy_diagonal_curvature = (
        enthalpies[:, 0] + enthalpies[:, 1] - 2.0 * center_enthalpy
    ) / (step * step)
    translation = np.eye(18) - translation_basis @ translation_basis.T
    acoustic_defect = float(np.linalg.norm(hessian @ translation) / max(np.linalg.norm(hessian), 1e-30))
    report = {
        "status": "GaN_image15_local_joint_curvature_at_nonstationary_NEB_candidate",
        "interpretation": "Local E+PV curvature only; not a TS Hessian-index certificate.",
        "highest_image_index": 15,
        "step_A": step,
        "pressure_GPa": manifest["pressure_GPa"],
        "n_static_displacements": len(details),
        "center_gradient_euclidean_eV_per_A": float(np.linalg.norm(center_gradient)),
        "center_gradient_translation_free_eV_per_A": float(np.linalg.norm(translation_basis.T @ center_gradient)),
        "center_max_abs_gradient_component_eV_per_A": float(np.max(np.abs(center_gradient))),
        "energy_gradient_max_abs_difference_eV_per_A": float(
            np.max(np.abs(energy_gradient - center_gradient))
        ),
        "energy_hessian_diagonal_max_abs_difference_eV_per_A2": float(
            np.max(np.abs(energy_diagonal_curvature - np.diag(hessian)))
        ),
        "raw_hessian_reciprocity_relative_defect": reciprocity_defect,
        "translation_null_relative_defect": acoustic_defect,
        "tangent_antisymmetric_deformation_norm": tangent_rotation,
        "translation_free_eigenvalues_eV_per_A2": eigenvalues.tolist(),
        "negative_counts_at_curvature_threshold_eV_per_A2": {
            str(cutoff): int(np.count_nonzero(eigenvalues < -cutoff))
            for cutoff in (0.01, 0.05, 0.1)
        },
        "tangent_squared_overlap_with_eigenvectors": np.square(full_vectors.T @ tangent).tolist(),
        "cases": details,
        "source_sha256": {
            "manifest": sha256(manifest_path), "trajectory": sha256(trajectory),
            "summary": sha256(summary), "auditor": sha256(Path(__file__)),
            "joint_curvature": sha256(Path(__file__).resolve().parents[1] / "vcneb/joint_curvature.py"),
            "center_OUTCAR": sha256(center_outcar),
        },
        "limitations": [
            "The center has a nonzero full enthalpy gradient.",
            "Only one displacement step is audited here; step-size stability needs an independent set.",
            "Rotations were excluded, translations projected; finite-q modes and NAC are absent.",
            "No local saddle refinement or two-sided basin connection has been performed.",
        ],
    }
    output_dir.mkdir(parents=True)
    np.savez_compressed(
        output_dir / "joint_hessian.npz", hessian=hessian,
        eigenvalues=eigenvalues, eigenvectors=full_vectors,
        center_gradient=center_gradient, tangent=tangent,
        paired_gradients=gradients, paired_enthalpies=enthalpies,
        energy_gradient=energy_gradient,
        energy_diagonal_curvature=energy_diagonal_curvature,
    )
    (output_dir / "audit.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.work_root, args.trajectory, args.summary, args.output_dir)
    print(json.dumps({
        "status": result["status"],
        "center_gradient": result["center_gradient_translation_free_eV_per_A"],
        "lowest_eigenvalues": result["translation_free_eigenvalues_eV_per_A2"][:3],
    }))


if __name__ == "__main__":
    main()
