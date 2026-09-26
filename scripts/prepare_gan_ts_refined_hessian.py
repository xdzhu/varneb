"""Stage 1x1x1-cell GaN joint atomic/strain probes at an audited near-TS.

This is an isolated 1000 eV diagnostic, not a change to the original 600 eV
VCNEB chain. Every static VASP input is copied from the completed center.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256
from vcneb.joint_curvature import JointCurvatureCoordinates


def prepare(center_dir: Path, center_manifest: Path, center_audit: Path,
            original_center_poscar: Path, work_root: Path, step_A: float) -> dict:
    if work_root.exists():
        raise FileExistsError(work_root)
    if not np.isfinite(step_A) or not 0 < step_A <= 0.025:
        raise ValueError("finite-difference step must be in (0, 0.025] Å")
    manifest = json.loads(center_manifest.read_text(encoding="utf-8"))
    audit = json.loads(center_audit.read_text(encoding="utf-8"))
    if (manifest.get("purpose") != "GaN_image15_followup_1000eV_local_stationarity_probe_not_TS_certificate"
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("evaluation_encut_eV") != 1000
            or audit.get("status") != "GaN_followup_local_trial_audited_not_TS_certificate"
            or audit.get("source_sha256", {}).get("manifest") != sha256(center_manifest)
            or audit.get("result", {}).get("status") != "complete_converged_static_not_TS_certificate"
            or len(manifest.get("cases", [])) != 1):
        raise ValueError("the 1000 eV near-TS source has not passed its independent audit")
    record = manifest["cases"][0]
    source, parsed = completed_case(
        center_dir, {**record, "POSCAR_sha256": record["input_sha256"]["POSCAR"]},
    )
    if (sha256(center_dir / "OUTCAR") != audit["result"]["outcar_sha256"]
            or parsed["input_sha256"] != audit["result"]["input_sha256"]):
        raise ValueError("raw near-TS center differs from its audit")
    if sha256(original_center_poscar) != manifest["source_sha256"]["center_POSCAR"]:
        raise ValueError("original image-15 coordinate reference changed")
    center = read(center_dir / "POSCAR", format="vasp")
    if not same_geometry(source, center, tolerance=2e-5):
        raise ValueError("the VASP output geometry differs from its staged center")
    incar = (center_dir / "INCAR").read_text(encoding="utf-8")
    kpoints = (center_dir / "KPOINTS").read_text(encoding="utf-8")
    if ("ENCUT = 1000.000000" not in incar or "EDIFF = 1.00e-07" not in incar
            or "ISYM = -1" not in incar or "SYMPREC = 1.00e-04" not in incar
            or "8 8 6" not in kpoints):
        raise ValueError("1000 eV GaN static input contract changed")
    coordinates = JointCurvatureCoordinates(center, cell_scale_A=3.3982714330050063)
    original_coordinates = JointCurvatureCoordinates(
        read(original_center_poscar, format="vasp"), cell_scale_A=coordinates.cell_scale_A,
    )
    original_gradient = original_coordinates.enthalpy_gradient(
        source, source.get_forces(), source.get_stress(voigt=False), 45.7 * GPa,
    )
    original_norm = float(np.linalg.norm(
        original_coordinates.translation_free_basis().T @ original_gradient
    ))
    if (not np.allclose(original_gradient, audit["result"]["gradient_components_eV_per_A"], atol=1e-7, rtol=0)
            or not np.isclose(original_norm, audit["result"]["translation_free_gradient_eV_per_A"], atol=1e-7, rtol=0)):
        raise ValueError("raw center no longer reproduces its original-coordinate audit")
    gradient = coordinates.enthalpy_gradient(
        source, source.get_forces(), source.get_stress(voigt=False), 45.7 * GPa,
    )
    gradient_norm = float(np.linalg.norm(coordinates.translation_free_basis().T @ gradient))
    if gradient_norm > 0.005:
        raise ValueError("center is not a sufficiently stationary audited candidate")

    work_root.mkdir(parents=True)
    shutil.copy2(center_dir / "POSCAR", work_root / "center_POSCAR")
    cases = []
    for axis in range(coordinates.size):
        for sign, label in ((1, "plus"), (-1, "minus")):
            shift = np.zeros(coordinates.size)
            shift[axis] = sign * step_A
            displaced = coordinates.displaced(shift)
            distances = displaced.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            minimum = float(np.min(distances))
            if (len(displaced) != 4 or minimum < 1.4
                    or np.linalg.cond(displaced.cell.array) > 50):
                raise ValueError(f"unsafe joint probe axis {axis} {label}")
            name = f"axis{axis:02d}_{label}"
            directory = work_root / "cases" / name
            directory.mkdir(parents=True)
            write(directory / "POSCAR", displaced, format="vasp", direct=True, sort=False)
            if not same_geometry(displaced, read(directory / "POSCAR", format="vasp")):
                raise ValueError(f"POSCAR roundtrip altered axis {axis} {label}")
            for input_name in ("INCAR", "KPOINTS", "POTCAR"):
                shutil.copy2(center_dir / input_name, directory / input_name)
            hashes = {input_name: sha256(directory / input_name)
                      for input_name in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
            (directory / "sha256.inputs.json").write_text(
                json.dumps(hashes, indent=2) + "\n", encoding="utf-8",
            )
            cases.append({
                "name": name, "axis": axis, "sign": sign, "input_sha256": hashes,
                "POSCAR_sha256": hashes["POSCAR"], "volume_A3": float(displaced.get_volume()),
                "minimum_distance_A": minimum,
            })
    result = {
        "purpose": "GaN_near_stationary_candidate_1000eV_joint_atomic_strain_Hessian_not_TS_certificate",
        "status": "inputs_finalized_no_DFT", "n_atoms": 4,
        "coordinate_count": 18, "physical_dimension_after_translations": 15,
        "step_A": step_A, "pressure_GPa": 45.7, "encut_eV": 1000,
        "cell_scale_A": coordinates.cell_scale_A,
        "phonon_supercell_matrix": np.eye(3, dtype=int).tolist(),
        "electronic_kpoints": [8, 8, 6],
        "center_translation_free_gradient_eV_per_A": gradient_norm,
        "center_original_coordinate_gradient_eV_per_A": original_norm,
        "center_gradient_components_eV_per_A": gradient.tolist(),
        "center_energy_eV_per_cell": float(source.get_potential_energy()),
        "center_enthalpy_eV_per_cell": float(
            source.get_potential_energy() + 45.7 * GPa * source.get_volume()
        ),
        "cases": cases,
        "source_sha256": {
            "center_manifest": sha256(center_manifest),
            "center_audit": sha256(center_audit),
            "original_center_POSCAR": sha256(original_center_poscar),
            "center_POSCAR": sha256(center_dir / "POSCAR"),
            "center_OUTCAR": sha256(center_dir / "OUTCAR"),
            "INCAR": sha256(center_dir / "INCAR"),
            "KPOINTS": sha256(center_dir / "KPOINTS"),
            "POTCAR": sha256(center_dir / "POTCAR"),
            "preparer": sha256(Path(__file__)),
            "joint_curvature": sha256(Path(__file__).resolve().parents[1] / "vcneb/joint_curvature.py"),
        },
        "limitations": "First finite-difference step at a near-stationary candidate; double-step and basin connections still required.",
    }
    (work_root / "manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--center-dir", type=Path, required=True)
    parser.add_argument("--center-manifest", type=Path, required=True)
    parser.add_argument("--center-audit", type=Path, required=True)
    parser.add_argument("--original-center-poscar", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--step-A", type=float, default=0.02)
    args = parser.parse_args()
    result = prepare(args.center_dir, args.center_manifest, args.center_audit,
                     args.original_center_poscar, args.work_root, args.step_A)
    print(json.dumps({"status": result["status"], "n_cases": len(result["cases"]),
                      "center_gradient_eV_per_A": result["center_translation_free_gradient_eV_per_A"]}))


if __name__ == "__main__":
    main()
