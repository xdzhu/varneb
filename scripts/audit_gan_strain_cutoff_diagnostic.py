"""Audit isolated 800 eV GaN strain statics against the 600 eV baseline."""

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


def audit(work_root: Path, baseline_manifest_path: Path,
          baseline_audit_path: Path, baseline_archive_path: Path) -> dict:
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    baseline = json.loads(baseline_manifest_path.read_text(encoding="utf-8"))
    prior = json.loads(baseline_audit_path.read_text(encoding="utf-8"))
    names = ["center"] + [f"axis{axis:02d}_{side}" for axis in (12, 13, 14)
                          for side in ("plus", "minus")]
    if (manifest.get("purpose") != "GaN_image15_VASP_strain_energy_stress_cutoff_sensitivity_only"
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("encut_eV") != 800
            or manifest.get("baseline_encut_eV") != 600
            or manifest.get("step_A") != 0.02
            or manifest.get("source_600_manifest_sha256") != sha256(baseline_manifest_path)
            or manifest.get("preparer_sha256") != sha256(
                Path(__file__).with_name("prepare_gan_strain_cutoff_diagnostic.py")
            )
            or [record["name"] for record in manifest.get("cases", [])] != names
            or baseline["source_static_sha256"]["OUTCAR"] != manifest["source_static_OUTCAR_sha256"]
            or prior["source_sha256"]["manifest"] != sha256(baseline_manifest_path)):
        raise ValueError("800 eV diagnostic or 600 eV source provenance gate failed")
    prior_cases = {record["name"]: record for record in baseline["cases"]}
    for record in manifest["cases"]:
        name = record["name"]
        expected_poscar = (baseline["center_POSCAR_sha256"] if name == "center"
                           else prior_cases[name]["POSCAR_sha256"])
        if record["POSCAR_sha256"] != expected_poscar:
            raise ValueError(f"800 eV probe geometry differs from 600 eV: {name}")
        for unchanged in ("KPOINTS", "POTCAR"):
            if record["input_sha256"][unchanged] != baseline["source_static_sha256"][unchanged]:
                raise ValueError(f"800 eV probe changed {unchanged}: {name}")
        incar = (work_root / "cases" / name / "INCAR").read_text(encoding="utf-8")
        if "ENCUT = 800.000000" not in incar or "ENCUT = 600.000000" in incar:
            raise ValueError(f"isolated VASP cutoff was not applied: {name}")
        restored = incar.replace("ENCUT = 800.000000", "ENCUT = 600.000000")
        if hashlib.sha256(restored.encode("utf-8")).hexdigest() != baseline["source_static_sha256"]["INCAR"]:
            raise ValueError(f"800 eV INCAR changed more than the cutoff: {name}")
    evaluated = {}
    details = []
    for record in manifest["cases"]:
        atoms, info = completed_case(work_root / "cases" / record["name"], record)
        evaluated[record["name"]] = atoms
        details.append(info)
    reference = evaluated["center"]
    coordinates = JointCurvatureCoordinates(reference, float(manifest["cell_scale_A"]))
    pressure = float(manifest["pressure_GPa"]) * GPa
    center_gradient = coordinates.enthalpy_gradient(
        reference, reference.get_forces(), reference.get_stress(voigt=False), pressure,
    )
    step = float(manifest["step_A"])
    energy_gradient_800 = []
    for axis in (12, 13, 14):
        plus, minus = (evaluated[f"axis{axis:02d}_{side}"] for side in ("plus", "minus"))
        h_plus = plus.get_potential_energy() + pressure * plus.get_volume()
        h_minus = minus.get_potential_energy() + pressure * minus.get_volume()
        energy_gradient_800.append((h_plus - h_minus) / (2 * step))
    energy_gradient_800 = np.array(energy_gradient_800)
    mismatch_800 = energy_gradient_800 - center_gradient[12:15]
    with np.load(baseline_archive_path, allow_pickle=False) as previous:
        mismatch_600 = previous["energy_gradient"][12:15] - previous["center_gradient"][12:15]
        center_gradient_600 = previous["center_gradient"][12:15].copy()
    volume = reference.get_volume()
    scale = coordinates.cell_scale_A
    kbar_per_eV_A3 = 1602.176634
    report = {
        "status": "isolated_GaN_VASP_800eV_strain_diagnostic_not_production_or_TS_certificate",
        "encut_eV": [600, 800],
        "step_A": step,
        "axes": [12, 13, 14],
        "energy_gradient_minus_stress_gradient_eV_per_A": {
            "600": mismatch_600.tolist(), "800": mismatch_800.tolist(),
        },
        "equivalent_stress_difference_kbar": {
            "600": (mismatch_600 * scale / volume * kbar_per_eV_A3).tolist(),
            "800": (mismatch_800 * scale / volume * kbar_per_eV_A3).tolist(),
        },
        "mean_absolute_gradient_mismatch_eV_per_A": {
            "600": float(np.mean(np.abs(mismatch_600))),
            "800": float(np.mean(np.abs(mismatch_800))),
        },
        "center_stress_gradient_eV_per_A": {
            "600": center_gradient_600.tolist(), "800": center_gradient[12:15].tolist(),
        },
        "center_800_energy_eV_per_cell": float(reference.get_potential_energy()),
        "center_800_volume_A3": volume,
        "cases": details,
        "source_sha256": {
            "manifest_800": sha256(manifest_path),
            "manifest_600": sha256(baseline_manifest_path),
            "audit_600": sha256(baseline_audit_path),
            "archive_600": sha256(baseline_archive_path),
            "auditor": sha256(Path(__file__)),
        },
        "interpretation_boundary": (
            "A single higher-cutoff comparison can support or weaken a finite-basis/Pulay "
            "hypothesis but cannot establish full stress convergence. The 600 eV VCNEB "
            "and BTO/ABACUS 100 Ry contracts remain untouched."
        ),
    }
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--baseline-manifest", type=Path, required=True)
    parser.add_argument("--baseline-audit", type=Path, required=True)
    parser.add_argument("--baseline-archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work_root, args.baseline_manifest, args.baseline_audit,
                   args.baseline_archive)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "mean_absolute_gradient_mismatch_eV_per_A",
        "equivalent_stress_difference_kbar",
    )}))


if __name__ == "__main__":
    main()
