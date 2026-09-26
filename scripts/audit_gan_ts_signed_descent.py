"""Audit both raw VASP statics along the GaN near-TS unstable mode.

Two lower-enthalpy signed probes establish only local downhill behavior.
They do not establish relaxation into distinct B4 and B1 basins.
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
from vcneb.joint_curvature import JointCurvatureCoordinates


def audit(work_root: Path, center_dir: Path, hessian_npz: Path,
          hessian_audit: Path, comparison: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    source = manifest.get("source_sha256", {})
    if (manifest.get("purpose") != "GaN_near_TS_1000eV_signed_unstable_mode_static_downhill_probe_not_basin_link"
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("phonon_supercell_matrix") != np.eye(3, dtype=int).tolist()
            or manifest.get("electronic_kpoints") != [8, 8, 6]
            or source.get("center_POSCAR") != sha256(center_dir / "POSCAR")
            or source.get("center_OUTCAR") != sha256(center_dir / "OUTCAR")
            or source.get("hessian_npz") != sha256(hessian_npz)
            or source.get("hessian_audit") != sha256(hessian_audit)
            or source.get("step_comparison") != sha256(comparison)
            or source.get("preparer") != sha256(Path(__file__).with_name("prepare_gan_ts_signed_descent.py"))
            or sha256(work_root / "center_POSCAR") != source["center_POSCAR"]
            or [r["name"] for r in manifest.get("cases", [])]
            != ["negative_mode_plus", "negative_mode_minus"]):
        raise ValueError("signed GaN source or 1x1x1 input provenance differs")
    for filename in ("INCAR", "KPOINTS", "POTCAR"):
        if source.get(filename) != sha256(center_dir / filename):
            raise ValueError(f"static VASP source changed: {filename}")
    center = read(center_dir / "POSCAR", format="vasp")
    center_output = read(center_dir / "OUTCAR")
    if not same_geometry(center, center_output, tolerance=2e-5):
        raise ValueError("center VASP output geometry differs")
    pressure = float(manifest["pressure_GPa"]) * GPa
    center_enthalpy = float(center_output.get_potential_energy() + pressure * center_output.get_volume())
    if not np.isclose(center_enthalpy, manifest["center_enthalpy_eV_per_cell"], atol=1e-8, rtol=0):
        raise ValueError("center enthalpy changed")
    with np.load(hessian_npz) as archive:
        vector = np.asarray(archive["eigenvectors"][:, 0], dtype=float)
    if vector.shape != (18,) or not np.isclose(np.linalg.norm(vector), 1.0, atol=1e-8):
        raise ValueError("source negative mode has invalid dimension or norm")
    chart = JointCurvatureCoordinates(center, float(manifest["cell_scale_A"]))
    details = []
    for record in manifest["cases"]:
        expected = chart.displaced(record["sign"] * float(manifest["amplitude_A"]) * vector)
        directory = work_root / "cases" / record["name"]
        if not same_geometry(expected, read(directory / "POSCAR", format="vasp")):
            raise ValueError(f"signed negative-mode geometry changed: {record['name']}")
        atoms, raw = completed_case(directory, record)
        enthalpy = float(atoms.get_potential_energy() + pressure * atoms.get_volume())
        gradient = chart.enthalpy_gradient(
            atoms, atoms.get_forces(), atoms.get_stress(voigt=False), pressure,
        )
        details.append({
            **raw, "sign": record["sign"], "enthalpy_eV_per_cell": enthalpy,
            "delta_enthalpy_from_center_eV_per_cell": enthalpy - center_enthalpy,
            "quadratic_prediction_eV_per_cell": record["quadratic_predicted_delta_enthalpy_eV"],
            "gradient_along_source_unstable_mode_eV_per_A": float(vector @ gradient),
        })
    result = {
        "status": "GaN_signed_unstable_mode_two_statics_audited_not_basin_links",
        "n_completed_raw_VASP_statics": 2,
        "amplitude_A": float(manifest["amplitude_A"]),
        "pressure_GPa": float(manifest["pressure_GPa"]),
        "center_enthalpy_eV_per_cell": center_enthalpy,
        "both_sides_downhill": all(row["delta_enthalpy_from_center_eV_per_cell"] < 0
                                   for row in details),
        "cases": details,
        "source_sha256": {
            "manifest": sha256(manifest_path), "center_OUTCAR": sha256(center_dir / "OUTCAR"),
            "hessian_audit": sha256(hessian_audit), "hessian_npz": sha256(hessian_npz),
            "step_comparison": sha256(comparison), "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "Signed static downhill checks do not prove B4/B1 basin connections.",
            "The original 600 eV VCNEB chain is unchanged; these are separate 1000 eV diagnostics.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"],
                      "both_sides_downhill": result["both_sides_downhill"],
                      "delta_enthalpies_eV": [row["delta_enthalpy_from_center_eV_per_cell"]
                                              for row in details]}))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work-root", "center-dir", "hessian-npz", "hessian-audit",
                 "comparison", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.work_root, args.center_dir, args.hessian_npz,
          args.hessian_audit, args.comparison, args.output)


if __name__ == "__main__":
    main()
