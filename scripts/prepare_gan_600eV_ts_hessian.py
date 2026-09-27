"""Stage and finalize 600-eV GaN joint Hessian inputs at the static canary.

The 4-atom center is job 27792422's audited static result. This finite
difference experiment tests local curvature only; it does not certify a TS.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.prepare_gan_600eV_ts_newton_canary import (
    PRODUCTION_INPUT_SHA256,
    basal_projection_A,
)
from vcneb.joint_curvature import JointCurvatureCoordinates


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same_geometry(a, b, tol: float = 2e-5) -> bool:
    q = a.get_scaled_positions(wrap=False) - b.get_scaled_positions(wrap=False)
    q -= np.rint(q)
    return (a.get_chemical_symbols() == b.get_chemical_symbols()
            and np.allclose(a.cell.array, b.cell.array, atol=tol, rtol=0)
            and np.allclose(q, 0, atol=tol, rtol=0))


def geometry_preflight(atoms) -> dict:
    distances = atoms.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    volume = float(atoms.get_volume())
    minimum = float(np.min(distances))
    projection = basal_projection_A(atoms.cell.array)
    # This is a case-specific warning window from historical successful and
    # failed GaN inputs, not a prediction of VASP's internal classifier.
    ambiguous = bool(1e-4 <= projection <= 3.5e-4)
    if not np.isfinite([volume, minimum, projection]).all() or volume <= 0 or minimum < 1.4:
        raise ValueError("unsafe finite-displacement geometry")
    return {
        "volume_A3": volume, "minimum_distance_A": minimum,
        "basal_projection_A": projection,
        "empirical_near_symmetry_warning": ambiguous,
    }


def prepare(center_work: Path, center_audit: Path, output: Path,
            cell_scale_A: float, step_A: float) -> dict:
    if output.exists():
        raise FileExistsError(output)
    if not 0 < step_A <= 0.02:
        raise ValueError("step must be in (0, 0.02] A")
    center_report = json.loads(center_audit.read_text(encoding="utf-8"))
    if (center_report.get("status") != "GaN_600eV_one_step_static_audited_not_TS_certificate"
            or center_report.get("slurm_job") != 27792422
            or center_report.get("OUTCAR_sha256") != sha256(center_work / "OUTCAR")
            or center_report.get("input_contract", {}).get("ENCUT_eV") != 600
            or center_report.get("candidate_gradient_translation_free_eV_per_A", 1) >= 0.01
            or any(sha256(center_work / name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("600-eV center has not passed the required raw audit")
    center = read(center_work / "OUTCAR")
    chart = JointCurvatureCoordinates(center, cell_scale_A)
    if center.get_chemical_symbols() != ["Ga", "Ga", "N", "N"] or chart.size != 18:
        raise ValueError("unexpected joint-coordinate center")
    output.mkdir(parents=True)
    (output / "cases").mkdir()
    write(output / "center_POSCAR", center, format="vasp", direct=True, sort=False)
    if not same_geometry(center, read(output / "center_POSCAR", format="vasp")):
        raise ValueError("center POSCAR roundtrip failed")
    records = []
    for axis in range(18):
        for sign, side in ((1, "plus"), (-1, "minus")):
            delta = np.zeros(18)
            delta[axis] = sign * step_A
            displaced = chart.displaced(delta)
            measurements = geometry_preflight(displaced)
            if measurements["empirical_near_symmetry_warning"]:
                raise ValueError(f"near-symmetry risk for axis{axis:02d}_{side}")
            name = f"axis{axis:02d}_{side}"
            target = output / "cases" / name
            target.mkdir()
            write(target / "POSCAR", displaced, format="vasp", direct=True, sort=False)
            if not same_geometry(displaced, read(target / "POSCAR", format="vasp")):
                raise ValueError(f"POSCAR roundtrip failed for {name}")
            records.append({
                "name": name, "axis": axis, "sign": sign,
                **measurements, "POSCAR_sha256": sha256(target / "POSCAR"),
            })
    manifest = {
        "status": "geometry_preflight_passed_no_DFT",
        "purpose": "GaN_600eV_candidate_joint_hessian_0p02_A",
        "step_A": step_A, "cell_scale_A": cell_scale_A,
        "n_atoms": 4, "coordinate_count": 18, "pressure_GPa": 45.7,
        "source_static_job": 27792422,
        "source_input_sha256": PRODUCTION_INPUT_SHA256,
        "center_POSCAR_sha256": sha256(output / "center_POSCAR"),
        "source_sha256": {
            "center_OUTCAR": sha256(center_work / "OUTCAR"),
            "center_audit": sha256(center_audit),
            "preparer": sha256(Path(__file__)),
        },
        "cases": records,
        "limitations": [
            "Empirical geometry preflight does not guarantee VASP parser success.",
            "A Hessian at this candidate needs numerical/force consistency and basin checks before TS claims.",
        ],
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def finalize(work: Path, source_static: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "geometry_preflight_passed_no_DFT"
            or manifest.get("purpose") != "GaN_600eV_candidate_joint_hessian_0p02_A"
            or manifest.get("source_sha256", {}).get("center_OUTCAR") != sha256(source_static / "OUTCAR")
            or any(sha256(source_static / name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("center or 600-eV electronic contract changed before finalization")
    if not same_geometry(read(work / "center_POSCAR", format="vasp"),
                         read(source_static / "OUTCAR")):
        raise ValueError("Hessian center is not the audited 600-eV static")
    if any((work / "cases" / item["name"] / "INCAR").exists()
           or sha256(work / "cases" / item["name"] / "POSCAR") != item["POSCAR_sha256"]
           for item in manifest["cases"]):
        raise ValueError("staged case changed or already finalized")
    for item in manifest["cases"]:
        target = work / "cases" / item["name"]
        for name in PRODUCTION_INPUT_SHA256:
            shutil.copy2(source_static / name, target / name)
        hashes = {name: sha256(target / name) for name in ("POSCAR", *PRODUCTION_INPUT_SHA256)}
        (target / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
        )
        item["input_sha256"] = hashes
    manifest["status"] = "inputs_finalized_no_DFT"
    manifest["source_static_directory"] = str(source_static)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    stage = modes.add_parser("prepare")
    for name in ("center-work", "center-audit", "output"):
        stage.add_argument("--" + name, type=Path, required=True)
    stage.add_argument("--cell-scale-A", type=float, required=True)
    stage.add_argument("--step-A", type=float, default=0.02)
    complete = modes.add_parser("finalize")
    complete.add_argument("--work", type=Path, required=True)
    complete.add_argument("--source-static", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "prepare":
        result = prepare(args.center_work, args.center_audit, args.output,
                         args.cell_scale_A, args.step_A)
    else:
        result = finalize(args.work, args.source_static)
    print(json.dumps({"status": result["status"], "cases": len(result["cases"])}))


if __name__ == "__main__":
    main()
