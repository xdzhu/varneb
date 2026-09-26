"""Audit isolated 1000 eV GaN strain probes against 600/800 eV evidence.

Only seven matched image-15 geometries are compared. This does not change
the published 600 eV VCNEB path or certify a stationary transition state.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from vcneb.joint_curvature import JointCurvatureCoordinates


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(work_root: Path, baseline_manifest_path: Path, baseline_audit_path: Path,
          baseline_archive_path: Path, report_800_path: Path) -> dict:
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    baseline = json.loads(baseline_manifest_path.read_text(encoding="utf-8"))
    prior = json.loads(baseline_audit_path.read_text(encoding="utf-8"))
    report_800 = json.loads(report_800_path.read_text(encoding="utf-8"))
    names = ["center"] + [f"axis{axis:02d}_{side}" for axis in (12, 13, 14)
                          for side in ("plus", "minus")]
    if (manifest.get("purpose") != "GaN_image15_VASP_strain_energy_stress_cutoff_sensitivity_only"
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("encut_eV") != 1000
            or manifest.get("baseline_encut_eV") != 600
            or manifest.get("step_A") != 0.02
            or manifest.get("source_600_manifest_sha256") != sha256(baseline_manifest_path)
            or manifest.get("preparer_sha256") != sha256(
                Path(__file__).with_name("prepare_gan_strain_cutoff_diagnostic.py")
            )
            or [record["name"] for record in manifest.get("cases", [])] != names
            or prior["source_sha256"]["manifest"] != sha256(baseline_manifest_path)
            or report_800["source_sha256"]["manifest_600"] != sha256(baseline_manifest_path)
            or report_800["source_sha256"]["audit_600"] != sha256(baseline_audit_path)
            or report_800["source_sha256"]["archive_600"] != sha256(baseline_archive_path)):
        raise ValueError("1000 eV stage or 600/800 eV evidence chain is inconsistent")
    baseline_cases = {case["name"]: case for case in baseline["cases"]}
    for record in manifest["cases"]:
        name = record["name"]
        expected_poscar = (baseline["center_POSCAR_sha256"] if name == "center"
                           else baseline_cases[name]["POSCAR_sha256"])
        if record["POSCAR_sha256"] != expected_poscar:
            raise ValueError(f"1000 eV geometry changed: {name}")
        for unchanged in ("KPOINTS", "POTCAR"):
            if record["input_sha256"][unchanged] != baseline["source_static_sha256"][unchanged]:
                raise ValueError(f"1000 eV input changed {unchanged}: {name}")
        incar = (work_root / "cases" / name / "INCAR").read_text(encoding="utf-8")
        if incar.count("ENCUT = 1000.000000") != 1:
            raise ValueError(f"1000 eV cutoff not isolated: {name}")
        restored = incar.replace("ENCUT = 1000.000000", "ENCUT = 600.000000")
        if hashlib.sha256(restored.encode("utf-8")).hexdigest() != baseline["source_static_sha256"]["INCAR"]:
            raise ValueError(f"1000 eV INCAR changed more than cutoff: {name}")
    evaluated = {}
    details = []
    for record in manifest["cases"]:
        atoms, info = completed_case(work_root / "cases" / record["name"], record)
        evaluated[record["name"]] = atoms
        details.append(info)
    reference = evaluated["center"]
    if len(reference) != 4 or not np.isfinite(reference.get_volume()):
        raise ValueError("invalid 1000 eV reference cell")
    coordinates = JointCurvatureCoordinates(reference, float(manifest["cell_scale_A"]))
    pressure = float(manifest["pressure_GPa"]) * GPa
    center_gradient = coordinates.enthalpy_gradient(
        reference, reference.get_forces(), reference.get_stress(voigt=False), pressure,
    )
    energy_gradient = []
    for axis in (12, 13, 14):
        plus, minus = (evaluated[f"axis{axis:02d}_{side}"] for side in ("plus", "minus"))
        h_plus = plus.get_potential_energy() + pressure * plus.get_volume()
        h_minus = minus.get_potential_energy() + pressure * minus.get_volume()
        energy_gradient.append((h_plus - h_minus) / (2 * float(manifest["step_A"])))
    mismatch_1000 = np.asarray(energy_gradient) - center_gradient[12:15]
    with np.load(baseline_archive_path, allow_pickle=False) as previous:
        mismatch_600 = previous["energy_gradient"][12:15] - previous["center_gradient"][12:15]
    mismatch_800 = np.asarray(
        report_800["energy_gradient_minus_stress_gradient_eV_per_A"]["800"], dtype=float,
    )
    if (mismatch_600.shape != (3,) or mismatch_800.shape != (3,)
            or not np.allclose(mismatch_600,
                               report_800["energy_gradient_minus_stress_gradient_eV_per_A"]["600"],
                               atol=1e-10, rtol=0)):
        raise ValueError("800 eV comparison does not reproduce the 600 eV baseline")
    factor = coordinates.cell_scale_A / reference.get_volume() * 1602.176634
    mismatches = {"600": mismatch_600, "800": mismatch_800, "1000": mismatch_1000}
    return {
        "status": "isolated_GaN_VASP_three_cutoff_stress_sensitivity_not_production_or_TS_certificate",
        "encut_eV": [600, 800, 1000],
        "axes": [12, 13, 14],
        "step_A": manifest["step_A"],
        "energy_gradient_minus_stress_gradient_eV_per_A": {
            key: value.tolist() for key, value in mismatches.items()
        },
        "equivalent_stress_difference_kbar": {
            key: (value * factor).tolist() for key, value in mismatches.items()
        },
        "mean_absolute_gradient_mismatch_eV_per_A": {
            key: float(np.mean(np.abs(value))) for key, value in mismatches.items()
        },
        "center_gradient_translation_free_eV_per_A": {
            "600": prior["center_gradient_translation_free_eV_per_A"],
            "1000": float(np.linalg.norm(
                coordinates.translation_free_basis().T @ center_gradient
            )),
        },
        "center_1000_stress_gradient_eV_per_A": center_gradient[12:15].tolist(),
        "center_1000_energy_eV_per_cell": float(reference.get_potential_energy()),
        "cases": details,
        "source_sha256": {
            "manifest_1000": sha256(manifest_path),
            "manifest_600": sha256(baseline_manifest_path),
            "audit_600": sha256(baseline_audit_path),
            "archive_600": sha256(baseline_archive_path),
            "report_800": sha256(report_800_path),
            "auditor": sha256(Path(__file__)),
        },
        "interpretation_boundary": (
            "Three local cutoff samples assess a finite-basis hypothesis, but do not "
            "replace full stress convergence of the production path or prove a TS."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--baseline-manifest", type=Path, required=True)
    parser.add_argument("--baseline-audit", type=Path, required=True)
    parser.add_argument("--baseline-archive", type=Path, required=True)
    parser.add_argument("--report-800", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work_root, args.baseline_manifest, args.baseline_audit,
                   args.baseline_archive, args.report_800)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "mean_absolute_gradient_mismatch_eV_per_A",
        "equivalent_stress_difference_kbar",
    )}))


if __name__ == "__main__":
    main()
