"""Independently audit 1000 eV GaN joint curvature at the near-stationary candidate.

All 36 VASP outputs, signed geometries, forces, stresses, energies and input
hashes are checked. One step and one negative mode do not certify basin links.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256
from vcneb.joint_curvature import JointCurvatureCoordinates, central_difference_hessian


def audit(work_root: Path, center_dir: Path, center_manifest: Path,
          center_audit: Path, output_dir: Path) -> dict:
    if output_dir.exists():
        raise FileExistsError(output_dir)
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    source_manifest = json.loads(center_manifest.read_text(encoding="utf-8"))
    source_audit = json.loads(center_audit.read_text(encoding="utf-8"))
    hashes = manifest["source_sha256"]
    expected_names = [f"axis{axis:02d}_{label}"
                      for axis in range(18) for label in ("plus", "minus")]
    if (manifest.get("purpose") != "GaN_near_stationary_candidate_1000eV_joint_atomic_strain_Hessian_not_TS_certificate"
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("coordinate_count") != 18
            or manifest.get("physical_dimension_after_translations") != 15
            or manifest.get("phonon_supercell_matrix") != np.eye(3, dtype=int).tolist()
            or manifest.get("electronic_kpoints") != [8, 8, 6]
            or manifest.get("encut_eV") != 1000
            or not np.isclose(manifest.get("pressure_GPa"), 45.7)
            or [record["name"] for record in manifest["cases"]] != expected_names
            or hashes["center_manifest"] != sha256(center_manifest)
            or hashes["center_audit"] != sha256(center_audit)
            or hashes["center_POSCAR"] != sha256(center_dir / "POSCAR")
            or hashes["center_OUTCAR"] != sha256(center_dir / "OUTCAR")
            or hashes["preparer"] != sha256(Path(__file__).with_name("prepare_gan_ts_refined_hessian.py"))
            or hashes["joint_curvature"] != sha256(Path(__file__).resolve().parents[1] / "vcneb/joint_curvature.py")
            or source_audit["source_sha256"]["manifest"] != sha256(center_manifest)
            or source_manifest["cases"][0]["input_sha256"]["POSCAR"] != hashes["center_POSCAR"]):
        raise ValueError("near-TS Hessian stage or audited source provenance changed")
    for name in ("INCAR", "KPOINTS", "POTCAR"):
        if hashes[name] != sha256(center_dir / name):
            raise ValueError(f"source VASP static input changed: {name}")
    center = read(work_root / "center_POSCAR", format="vasp")
    reference = read(center_dir / "POSCAR", format="vasp")
    baseline = read(center_dir / "OUTCAR")
    if (sha256(work_root / "center_POSCAR") != hashes["center_POSCAR"]
            or not same_geometry(center, reference)
            or not same_geometry(center, baseline, tolerance=2e-5)):
        raise ValueError("center geometry differs from completed VASP source")
    coordinates = JointCurvatureCoordinates(center, float(manifest["cell_scale_A"]))
    pressure = float(manifest["pressure_GPa"]) * GPa
    center_gradient = coordinates.enthalpy_gradient(
        baseline, baseline.get_forces(), baseline.get_stress(voigt=False), pressure,
    )
    physical_basis = coordinates.translation_free_basis()
    center_norm = float(np.linalg.norm(physical_basis.T @ center_gradient))
    center_enthalpy = float(baseline.get_potential_energy() + pressure * baseline.get_volume())
    if (not np.allclose(center_gradient, manifest["center_gradient_components_eV_per_A"], atol=1e-8, rtol=0)
            or not np.isclose(center_norm, manifest["center_translation_free_gradient_eV_per_A"], atol=1e-8, rtol=0)
            or not np.isclose(center_enthalpy, manifest["center_enthalpy_eV_per_cell"], atol=1e-7, rtol=0)):
        raise ValueError("center raw energy, force or stress differs from staged manifest")
    step = float(manifest["step_A"])
    gradients = np.zeros((18, 2, 18))
    enthalpies = np.zeros((18, 2))
    details = []
    for record in manifest["cases"]:
        axis, sign = record["axis"], record["sign"]
        if record["name"] != f"axis{axis:02d}_{'plus' if sign == 1 else 'minus'}" or sign not in (1, -1):
            raise ValueError("signed curvature case index changed")
        directory = work_root / "cases" / record["name"]
        expected_shift = np.zeros(18)
        expected_shift[axis] = sign * step
        expected = coordinates.displaced(expected_shift)
        if not same_geometry(expected, read(directory / "POSCAR", format="vasp")):
            raise ValueError(f"signed probe geometry changed: {record['name']}")
        atoms, detail = completed_case(directory, record)
        index = 0 if sign == 1 else 1
        gradients[axis, index] = coordinates.enthalpy_gradient(
            atoms, atoms.get_forces(), atoms.get_stress(voigt=False), pressure,
        )
        enthalpies[axis, index] = atoms.get_potential_energy() + pressure * atoms.get_volume()
        details.append(detail)
    hessian, reciprocity = central_difference_hessian(
        gradients[:, 0], gradients[:, 1], step,
    )
    reduced = physical_basis.T @ hessian @ physical_basis
    eigenvalues, reduced_vectors = np.linalg.eigh(reduced)
    full_vectors = physical_basis @ reduced_vectors
    translation = np.eye(18) - physical_basis @ physical_basis.T
    acoustic_defect = float(np.linalg.norm(hessian @ translation)
                            / max(np.linalg.norm(hessian), 1e-30))
    energy_gradient = (enthalpies[:, 0] - enthalpies[:, 1]) / (2 * step)
    energy_diagonal_curvature = (
        enthalpies[:, 0] + enthalpies[:, 1] - 2 * center_enthalpy
    ) / (step * step)
    report = {
        "status": "GaN_near_stationary_1000eV_joint_curvature_one_step_audited_not_TS_certificate",
        "n_static_displacements": len(details), "step_A": step,
        "pressure_GPa": manifest["pressure_GPa"], "encut_eV": 1000,
        "center_gradient_translation_free_eV_per_A": center_norm,
        "center_enthalpy_eV_per_cell": center_enthalpy,
        "energy_gradient_max_abs_difference_eV_per_A": float(np.max(np.abs(energy_gradient - center_gradient))),
        "energy_hessian_diagonal_max_abs_difference_eV_per_A2": float(
            np.max(np.abs(energy_diagonal_curvature - np.diag(hessian)))
        ),
        "raw_hessian_reciprocity_relative_defect": reciprocity,
        "translation_null_relative_defect": acoustic_defect,
        "translation_free_eigenvalues_eV_per_A2": eigenvalues.tolist(),
        "negative_counts_at_curvature_threshold_eV_per_A2": {
            str(cutoff): int(np.count_nonzero(eigenvalues < -cutoff))
            for cutoff in (0.01, 0.05, 0.1)
        },
        "cases": details,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "center_audit": sha256(center_audit),
            "center_OUTCAR": sha256(center_dir / "OUTCAR"),
            "auditor": sha256(Path(__file__)),
            "joint_curvature": sha256(Path(__file__).resolve().parents[1] / "vcneb/joint_curvature.py"),
        },
        "limitations": [
            "Only one finite-difference step at the near-stationary 1000 eV candidate is audited.",
            "A second step-size test and both downhill basin connections remain required.",
            "Translations are projected and rotations excluded; finite-q phonons and NAC are not tested.",
        ],
    }
    output_dir.mkdir(parents=True)
    np.savez_compressed(
        output_dir / "joint_hessian.npz", hessian=hessian, eigenvalues=eigenvalues,
        eigenvectors=full_vectors, center_gradient=center_gradient,
        paired_gradients=gradients, paired_enthalpies=enthalpies,
        energy_gradient=energy_gradient,
        energy_diagonal_curvature=energy_diagonal_curvature,
    )
    (output_dir / "audit.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work-root", "center-dir", "center-manifest", "center-audit", "output-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.work_root, args.center_dir, args.center_manifest,
                   args.center_audit, args.output_dir)
    print(json.dumps({"status": result["status"],
                      "center_gradient": result["center_gradient_translation_free_eV_per_A"],
                      "lowest_eigenvalues": result["translation_free_eigenvalues_eV_per_A2"][:3]}))


if __name__ == "__main__":
    main()
