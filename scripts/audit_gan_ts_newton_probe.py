"""Audit half/full GaN local Newton probes without claiming a TS."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_ts_newton_probe import same_geometry
from vcneb.joint_curvature import JointCurvatureCoordinates


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def failed_bravais_case(directory: Path, record: dict) -> dict:
    """Accept only the observed pre-SCF VASP lattice-classification failure."""

    hashes = json.loads((directory / "sha256.inputs.json").read_text(encoding="utf-8"))
    if hashes != record["input_sha256"]:
        raise ValueError(f"failed-case input hash manifest mismatch: {directory}")
    for filename, expected in hashes.items():
        if sha256(directory / filename) != expected:
            raise ValueError(f"failed-case input changed: {directory / filename}")
    outcar = directory / "OUTCAR"
    stdout = directory / "vasp.stdout"
    oszicar = directory / "OSZICAR"
    if not (outcar.is_file() and stdout.is_file() and oszicar.is_file()):
        raise ValueError(f"missing failed-case evidence: {directory}")
    outcar_text = outcar.read_text(encoding="utf-8", errors="replace")
    stdout_text = stdout.read_text(encoding="utf-8", errors="replace")
    signature = "Inconsistent Bravais lattice types found for crystalline and"
    classes = ("Crystalline: simple monoclinic",
               "Reciprocal : base centered orthorhombic")
    if (signature not in outcar_text or signature not in stdout_text
            or any(item not in outcar_text or item not in stdout_text for item in classes)
            or "General timing and accounting informations" in outcar_text
            or "aborting loop because EDIFF is reached" in outcar_text
            or oszicar.stat().st_size != 0):
        raise ValueError(f"not the audited pre-SCF Bravais failure: {directory}")
    return {
        "case": record["name"],
        "status": "VASP_6p3p2_Bravais_classification_failure_before_SCF",
        "input_sha256": hashes,
        "outcar_sha256": sha256(outcar),
        "stdout_sha256": sha256(stdout),
        "oszicar_sha256": sha256(oszicar),
        "energy_eV_per_cell": None,
        "translation_free_gradient_eV_per_A": None,
        "gradient_ratio_to_center": None,
    }


def audit(work_root: Path, trajectory: Path, summary: Path, hessian_archive: Path,
          center_1000_outcar: Path) -> dict:
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_source = {
        "trajectory": sha256(trajectory), "summary": sha256(summary),
        "hessian_archive": sha256(hessian_archive),
        "center_1000_OUTCAR": sha256(center_1000_outcar),
        "preparer": sha256(Path(__file__).with_name("prepare_gan_ts_newton_probe.py")),
        "joint_curvature": sha256(Path(__file__).resolve().parents[1] / "vcneb/joint_curvature.py"),
    }
    if (manifest.get("purpose") != "GaN_image15_1000eV_two_step_local_Newton_stationarity_probe"
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("source_sha256") != expected_source
            or manifest.get("evaluation_encut_eV") != 1000
            or [record["name"] for record in manifest.get("cases", [])]
            != ["newton_half", "newton_full"]
            or manifest["source_static_sha256"]["OUTCAR"] != sha256(center_1000_outcar)):
        raise ValueError("Newton probe provenance or input contract failed")
    reference = read(trajectory, index=15)
    original = read(center_1000_outcar)
    if not same_geometry(reference, original, tolerance=2e-5):
        raise ValueError("1000 eV center geometry differs from source image 15")
    coordinates = JointCurvatureCoordinates(reference, float(manifest["cell_scale_A"]))
    pressure = float(manifest["pressure_GPa"]) * GPa
    center_gradient = coordinates.enthalpy_gradient(
        original, original.get_forces(), original.get_stress(voigt=False), pressure,
    )
    basis = coordinates.translation_free_basis()
    center_norm = float(np.linalg.norm(basis.T @ center_gradient))
    if not np.isclose(center_norm, manifest["center_gradient_translation_free_eV_per_A"],
                      rtol=0, atol=1e-7):
        raise ValueError("1000 eV center gradient changed from preflight")
    center_h = float(original.get_potential_energy() + pressure * original.get_volume())
    records = []
    for case in manifest["cases"]:
        directory = work_root / "cases" / case["name"]
        outcar_text = (directory / "OUTCAR").read_text(encoding="utf-8", errors="replace")
        if "Inconsistent Bravais lattice types found for crystalline and" in outcar_text:
            records.append({**failed_bravais_case(directory, case),
                            "step_fraction": case["step_fraction"]})
            continue
        atoms, details = completed_case(directory, case)
        gradient = coordinates.enthalpy_gradient(
            atoms, atoms.get_forces(), atoms.get_stress(voigt=False), pressure,
        )
        current_norm = float(np.linalg.norm(basis.T @ gradient))
        enthalpy = float(atoms.get_potential_energy() + pressure * atoms.get_volume())
        records.append({
            **details,
            "step_fraction": case["step_fraction"],
            "enthalpy_eV_per_cell": enthalpy,
            "enthalpy_change_eV_per_cell": enthalpy - center_h,
            "translation_free_gradient_eV_per_A": current_norm,
            "gradient_ratio_to_center": current_norm / center_norm,
            "gradient_components_eV_per_A": gradient.tolist(),
            "status": "complete_converged_static",
        })
    completed = [record for record in records if record["status"] == "complete_converged_static"]
    if not completed:
        raise ValueError("no completed Newton candidate")
    return {
        "status": "GaN_1000eV_local_Newton_probe_partially_completed_not_TS_certificate",
        "center_gradient_eV_per_A": center_norm,
        "center_gradient_components_eV_per_A": center_gradient.tolist(),
        "center_enthalpy_eV_per_cell": center_h,
        "best_trial": min(completed, key=lambda item: item["translation_free_gradient_eV_per_A"])["case"],
        "cases": records,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "trajectory": sha256(trajectory),
            "summary": sha256(summary),
            "hessian_archive": sha256(hessian_archive),
            "center_1000_OUTCAR": sha256(center_1000_outcar),
            "auditor": sha256(Path(__file__)),
        },
        "interpretation_boundary": (
            "The reduced half-step gradient is not TS convergence. The full step failed "
            "before SCF with a VASP direct/reciprocal Bravais-classification conflict; "
            "it has no energy or gradient. Never silently change symmetry settings or "
            "the crystal metric to make that candidate run. Rebuild a same-contract "
            "Hessian at a stationary point and verify both basin connections."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--hessian", type=Path, required=True)
    parser.add_argument("--center-1000-outcar", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.work_root, args.trajectory, args.summary,
                   args.hessian, args.center_1000_outcar)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"],
                      "best_trial": report["best_trial"],
                      "gradient_ratios": [case["gradient_ratio_to_center"]
                                          for case in report["cases"]]}))


if __name__ == "__main__":
    main()
