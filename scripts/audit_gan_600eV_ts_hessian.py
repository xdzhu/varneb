"""Audit the 600-eV GaN candidate's full atom--strain finite-difference Hessian."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case, tangent_in_joint_coordinates
from scripts.prepare_gan_600eV_ts_hessian import same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates, central_difference_hessian


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(work: Path, center_work: Path, center_audit: Path,
          trajectory: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    canary = json.loads(center_audit.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_600eV_candidate_joint_hessian_0p02_A"
            or manifest.get("step_A") != 0.02
            or manifest.get("coordinate_count") != 18
            or manifest.get("source_static_job") != 27792422
            or manifest.get("source_sha256", {}).get("center_OUTCAR") != sha256(center_work / "OUTCAR")
            or manifest.get("source_sha256", {}).get("center_audit") != sha256(center_audit)
            or manifest.get("source_sha256", {}).get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_ts_hessian.py"))
            or manifest.get("source_input_sha256") != PRODUCTION_INPUT_SHA256
            or canary.get("OUTCAR_sha256") != sha256(center_work / "OUTCAR")):
        raise ValueError("600-eV Hessian center/source provenance failed")
    names = [f"axis{i:02d}_{side}" for i in range(18) for side in ("plus", "minus")]
    if [item.get("name") for item in manifest.get("cases", [])] != names:
        raise ValueError("36 signed Hessian displacement cases are incomplete or reordered")
    baseline = read(center_work / "OUTCAR")
    center = read(work / "center_POSCAR", format="vasp")
    if (manifest["center_POSCAR_sha256"] != sha256(work / "center_POSCAR")
            or not same_geometry(center, baseline)):
        raise ValueError("Hessian center differs from audited static")
    chart = JointCurvatureCoordinates(center, float(manifest["cell_scale_A"]))
    basis = chart.translation_free_basis()
    pressure = float(manifest["pressure_GPa"]) * GPa
    center_gradient = chart.enthalpy_gradient(
        baseline, baseline.get_forces(), baseline.get_stress(voigt=False), pressure
    )
    center_norm = float(np.linalg.norm(basis.T @ center_gradient))
    if not np.isclose(center_norm, canary["candidate_gradient_translation_free_eV_per_A"],
                      rtol=0, atol=1e-6):
        # The canary reported this norm in image-15 reference coordinates;
        # the Hessian uses the refined center as its deformation reference.
        raise ValueError("center full gradient differs from the canary audit")
    gradients = np.zeros((18, 2, 18))
    enthalpies = np.zeros((18, 2))
    details = []
    for record in manifest["cases"]:
        directory = work / "cases" / record["name"]
        if (record.get("input_sha256", {}).get("POSCAR") != record["POSCAR_sha256"]
                or any(record.get("input_sha256", {}).get(name) != digest
                       for name, digest in PRODUCTION_INPUT_SHA256.items())):
            raise ValueError(f"case {record['name']} deviates from 600-eV contract")
        atoms, detail = completed_case(directory, record)
        gradient = chart.enthalpy_gradient(
            atoms, atoms.get_forces(), atoms.get_stress(voigt=False), pressure
        )
        side = 0 if record["sign"] == 1 else 1
        gradients[record["axis"], side] = gradient
        enthalpies[record["axis"], side] = (
            atoms.get_potential_energy() + pressure * atoms.get_volume()
        )
        details.append(detail)
    hessian, reciprocity_defect = central_difference_hessian(
        gradients[:, 0], gradients[:, 1], float(manifest["step_A"])
    )
    reduced = basis.T @ hessian @ basis
    eigenvalues, eigenvectors_reduced = np.linalg.eigh(reduced)
    eigenvectors = basis @ eigenvectors_reduced
    chain = read(trajectory, index=":")
    if len(chain) != 29 or chain[15].get_chemical_symbols() != center.get_chemical_symbols():
        raise ValueError("original 600-eV chain is not the expected 29 images")
    tangent, tangent_rotation = tangent_in_joint_coordinates(chart, chain[14], chain[16])
    tangent = basis @ (basis.T @ tangent)
    tangent /= np.linalg.norm(tangent)
    center_h = float(baseline.get_potential_energy() + pressure * baseline.get_volume())
    energy_gradient = (enthalpies[:, 0] - enthalpies[:, 1]) / (2 * manifest["step_A"])
    energy_diagonal = (enthalpies[:, 0] + enthalpies[:, 1] - 2 * center_h) / (manifest["step_A"] ** 2)
    translation = np.eye(18) - basis @ basis.T
    null_defect = float(np.linalg.norm(hessian @ translation)
                        / max(np.linalg.norm(hessian), 1e-30))
    result = {
        "status": "GaN_600eV_near_stationary_joint_hessian_one_step_not_TS_certificate",
        "n_static_displacements": 36,
        "pressure_GPa": manifest["pressure_GPa"],
        "step_A": manifest["step_A"],
        "center_gradient_translation_free_eV_per_A": center_norm,
        "energy_gradient_max_abs_difference_eV_per_A": float(np.max(np.abs(energy_gradient - center_gradient))),
        "energy_hessian_diagonal_max_abs_difference_eV_per_A2": float(np.max(np.abs(energy_diagonal - np.diag(hessian)))),
        "raw_reciprocity_relative_defect": reciprocity_defect,
        "translation_null_relative_defect": null_defect,
        "path_tangent_rotation_norm": tangent_rotation,
        "eigenvalues_eV_per_A2": eigenvalues.tolist(),
        "negative_count_below_minus_0p1_eV_per_A2": int(np.count_nonzero(eigenvalues < -0.1)),
        "path_tangent_squared_overlap_with_modes": np.square(eigenvectors.T @ tangent).tolist(),
        "cases": details,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "center_OUTCAR": sha256(center_work / "OUTCAR"),
            "center_audit": sha256(center_audit),
            "trajectory": sha256(trajectory),
            "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "One finite-difference step cannot establish step-size stability.",
            "The full energy/stress derivative consistency is a separate numerical gate.",
            "Basin links must be reverified at the same 600-eV electronic contract.",
        ],
    }
    output.mkdir(parents=True)
    np.savez_compressed(
        output / "joint_hessian.npz", hessian=hessian,
        eigenvalues=eigenvalues, eigenvectors=eigenvectors,
        center_gradient=center_gradient, tangent=tangent,
        paired_gradients=gradients, paired_enthalpies=enthalpies,
        energy_gradient=energy_gradient, energy_diagonal_curvature=energy_diagonal,
    )
    (output / "audit.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "center-work", "center-audit", "trajectory", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.work, args.center_work, args.center_audit,
                   args.trajectory, args.output)
    print(json.dumps({
        "status": result["status"],
        "lowest_two_eigenvalues": result["eigenvalues_eV_per_A2"][:2],
        "center_gradient": result["center_gradient_translation_free_eV_per_A"],
    }))


if __name__ == "__main__":
    main()
