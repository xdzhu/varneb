"""Audit the single guarded GaN Newton trial, including a possible VASP failure."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.audit_gan_ts_newton_probe import failed_bravais_case
from scripts.prepare_gan_ts_bracketed_trial import cross_axis_skew_angstrom
from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256
from vcneb.joint_curvature import JointCurvatureCoordinates


def audit(work_root: Path, trajectory: Path, previous_manifest: Path,
          previous_audit: Path, center_1000_outcar: Path) -> dict:
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    previous = json.loads(previous_manifest.read_text(encoding="utf-8"))
    prior_result = json.loads(previous_audit.read_text(encoding="utf-8"))
    if (manifest["purpose"] != "GaN_image15_bracketed_local_Newton_trial_not_TS_certificate"
            or manifest["status"] != "inputs_finalized_no_DFT"
            or manifest["evaluation_encut_eV"] != 1000
            or len(manifest["cases"]) != 1
            or manifest["source_sha256"]["trajectory"] != sha256(trajectory)
            or manifest["source_sha256"]["previous_manifest"] != sha256(previous_manifest)
            or manifest["source_sha256"]["previous_audit"] != sha256(previous_audit)
            or manifest["source_sha256"]["center_1000_OUTCAR"] != sha256(center_1000_outcar)
            or manifest["source_sha256"]["preparer"] != sha256(
                Path(__file__).with_name("prepare_gan_ts_bracketed_trial.py")
            )):
        raise ValueError("guarded-trial provenance or VASP contract failed")
    case = manifest["cases"][0]
    if (case["name"] != "newton_0p625"
            or case["fraction_from_original_center"] != 0.625
            or case["input_sha256"]["INCAR"] != previous["source_static_sha256"]["INCAR"]
            or case["input_sha256"]["KPOINTS"] != previous["source_static_sha256"]["KPOINTS"]
            or case["input_sha256"]["POTCAR"] != previous["source_static_sha256"]["POTCAR"]):
        raise ValueError("candidate name/fraction or static input contract changed")
    reference = read(trajectory, index=15)
    coordinates = JointCurvatureCoordinates(reference, float(manifest["cell_scale_A"]))
    correction = np.asarray(previous["correction_A"], dtype=float)
    candidate = coordinates.displaced(0.625 * correction)
    directory = work_root / "cases" / case["name"]
    supplied = read(directory / "POSCAR", format="vasp")
    if (not same_geometry(candidate, supplied)
            or not np.isclose(cross_axis_skew_angstrom(candidate.cell.array),
                              case["cross_axis_skew_A"], atol=1e-10, rtol=0)
            or not case["cross_axis_skew_A"]
            < manifest["empirical_bravais_risk"]["allowed_skew_A"]):
        raise ValueError("candidate geometry or empirical risk check changed")
    outcar = directory / "OUTCAR"
    output = outcar.read_text(encoding="utf-8", errors="replace")
    if "Inconsistent Bravais lattice types found for crystalline and" in output:
        result = failed_bravais_case(directory, case)
    else:
        atoms, result = completed_case(directory, case)
        if result["input_sha256"] != case["input_sha256"]:
            raise ValueError("completed-trial hashes differ from staged input")
        pressure = float(manifest["pressure_GPa"]) * GPa
        gradient = coordinates.enthalpy_gradient(
            atoms, atoms.get_forces(), atoms.get_stress(voigt=False), pressure,
        )
        norm = float(np.linalg.norm(coordinates.translation_free_basis().T @ gradient))
        previous_norm = float(prior_result["cases"][0]["translation_free_gradient_eV_per_A"])
        result.update({
            "status": "complete_converged_static_not_TS_certificate",
            "enthalpy_eV_per_cell": float(atoms.get_potential_energy() + pressure * atoms.get_volume()),
            "translation_free_gradient_eV_per_A": norm,
            "gradient_ratio_to_previous_half": norm / previous_norm,
            "gradient_ratio_to_original_center": norm / prior_result["center_gradient_eV_per_A"],
            "gradient_components_eV_per_A": gradient.tolist(),
        })
    return {
        "status": "GaN_bracketed_local_trial_audited_not_TS_certificate",
        "result": result,
        "empirical_risk": manifest["empirical_bravais_risk"],
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "trajectory": sha256(trajectory),
            "previous_manifest": sha256(previous_manifest),
            "previous_audit": sha256(previous_audit),
            "center_1000_OUTCAR": sha256(center_1000_outcar),
            "auditor": sha256(Path(__file__)),
        },
        "interpretation_boundary": (
            "The empirical lattice margin is case-local, not a guarantee against VASP "
            "initialization failure. A reduced candidate gradient is not a certified TS."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--previous-manifest", type=Path, required=True)
    parser.add_argument("--previous-audit", type=Path, required=True)
    parser.add_argument("--center-1000-outcar", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work_root, args.trajectory, args.previous_manifest,
                   args.previous_audit, args.center_1000_outcar)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "result": result["result"]["status"],
                      "gradient": result["result"]["translation_free_gradient_eV_per_A"]}))


if __name__ == "__main__":
    main()
