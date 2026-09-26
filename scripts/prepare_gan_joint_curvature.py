"""Stage non-expanding GaN highest-image atomic/strain curvature probes.

The preparation mode writes only 4-atom POSCARs and a manifest. The separate
finalize mode runs on hf, validates the original VASP image-15 static source,
and copies its immutable electronic input into the staged displacement cases.
Neither mode runs VASP or changes the original NEB working directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

POTCAR_SHA256 = "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same_geometry(a, b, *, tolerance: float = 1e-7) -> bool:
    if a.get_chemical_symbols() != b.get_chemical_symbols():
        return False
    delta = (a.get_scaled_positions(wrap=False) - b.get_scaled_positions(wrap=False))
    delta -= np.rint(delta)
    return (
        np.allclose(a.cell.array, b.cell.array, rtol=0, atol=tolerance)
        and np.allclose(delta, 0, rtol=0, atol=tolerance)
    )


def prepare(trajectory: Path, summary: Path, work_root: Path, step_A: float) -> dict:
    from vcneb.joint_curvature import JointCurvatureCoordinates

    if work_root.exists():
        raise FileExistsError(f"refusing to overwrite existing work: {work_root}")
    if not np.isfinite(step_A) or not 0 < step_A <= 0.05:
        raise ValueError("joint displacement step must be in (0, 0.05] Å")
    metadata = json.loads(summary.read_text(encoding="utf-8"))
    if (metadata.get("status") != "completed" or not metadata.get("converged")
            or metadata.get("n_images") != 29
            or metadata.get("endpoint_evaluation_policy") != "fixed_cached_once"):
        raise ValueError("source VCNEB chain has not passed its existing completion gate")
    chain = read(trajectory, index="-29:")
    if len(chain) != 29:
        raise ValueError("the public trajectory does not contain a complete 29-image chain")
    candidate = chain[15]
    if candidate.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]:
        raise ValueError("unexpected highest-image atom identity")
    scale = float(metadata["path_diagnostics"]["cell_scale_A"])
    coordinates = JointCurvatureCoordinates(candidate, cell_scale_A=scale)
    if coordinates.size != 18:
        raise ValueError("unexpected GaN four-atom joint coordinate dimension")
    work_root.mkdir(parents=True)
    (work_root / "cases").mkdir()
    (work_root / "logs").mkdir()
    center = work_root / "center_POSCAR"
    write(center, candidate, format="vasp", direct=True, sort=False, vasp5=True)
    if not same_geometry(candidate, read(center, format="vasp")):
        raise ValueError("center POSCAR differs from the audited trajectory image")
    cases = []
    for axis in range(coordinates.size):
        for sign, label in ((1, "plus"), (-1, "minus")):
            shift = np.zeros(coordinates.size)
            shift[axis] = sign * step_A
            displaced = coordinates.displaced(shift)
            distances = displaced.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            minimum_distance = float(distances.min())
            if minimum_distance < 1.4:
                raise ValueError(f"too-close atoms in axis {axis} {label}: {minimum_distance}")
            name = f"axis{axis:02d}_{label}"
            target = work_root / "cases" / name
            target.mkdir()
            poscar = target / "POSCAR"
            write(poscar, displaced, format="vasp", direct=True, sort=False, vasp5=True)
            if not same_geometry(displaced, read(poscar, format="vasp")):
                raise ValueError(f"VASP POSCAR roundtrip changed {name}")
            cases.append({
                "name": name, "axis": axis, "sign": sign,
                "n_atoms": len(displaced), "volume_A3": displaced.get_volume(),
                "minimum_distance_A": minimum_distance, "POSCAR_sha256": sha256(poscar),
            })
    manifest = {
        "purpose": "GaN_B4_B1_highest_image_local_joint_atomic_symmetric_strain_curvature",
        "status": "geometry_staged_no_DFT",
        "highest_image_index": 15,
        "n_images_total": 29,
        "n_atoms": 4,
        "coordinate_count": 18,
        "physical_dimension_after_translation_removal": 15,
        "cell_scale_A": scale,
        "pressure_GPa": 45.7,
        "step_A": step_A,
        "supercell_matrix": np.eye(3, dtype=int).tolist(),
        "center_POSCAR_sha256": sha256(center),
        "source_sha256": {
            "trajectory": sha256(trajectory), "summary": sha256(summary),
            "preparer": sha256(Path(__file__)),
            "joint_curvature": sha256(Path(__file__).resolve().parents[1] / "vcneb/joint_curvature.py"),
        },
        "cases": cases,
        "limitations": (
            "Candidate is not stationary; this finite-difference Hessian is not a TS certificate. "
            "Rotation gauge is omitted, translations must be projected, and two-sided "
            "downhill connections remain untested."
        ),
    }
    (work_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8",
    )
    return manifest


def finalize_source(work_root: Path, source: Path) -> dict:
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "geometry_staged_no_DFT":
        raise ValueError("work root is not an unfinalized geometry stage")
    required = ("CONTCAR", "POSCAR.final", "INCAR", "KPOINTS", "POTCAR", "OUTCAR")
    if any(not (source / name).is_file() or (source / name).stat().st_size == 0 for name in required):
        raise ValueError("original image-15 static source is incomplete")
    if sha256(source / "POTCAR") != POTCAR_SHA256:
        raise ValueError("image-15 POTCAR differs from Ga_d+N audited production input")
    original = read(source / "POSCAR.final", format="vasp")
    center = read(work_root / "center_POSCAR", format="vasp")
    if not same_geometry(original, center) or not same_geometry(
        read(source / "CONTCAR", format="vasp"), center,
    ):
        raise ValueError("image-15 static geometry differs from the staged highest image")
    incar = (source / "INCAR").read_text(encoding="utf-8")
    parsed = {}
    for line in incar.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            parsed[key.strip().upper()] = value.strip().split("#", 1)[0].strip()
    for key, expected in {
        "ENCUT": 600, "EDIFF": 1e-7, "SYMPREC": 1e-4, "ISYM": -1,
        "IBRION": -1, "NSW": 0, "ISMEAR": 0, "SIGMA": 0.05,
    }.items():
        if key not in parsed or not np.isclose(float(parsed[key]), expected, atol=1e-12, rtol=1e-9):
            raise ValueError(f"image-15 VASP contract mismatch: {key}")
    if parsed.get("GGA") != "PE" or parsed.get("PREC", "").lower() != "accurate":
        raise ValueError("image-15 PBE/PREC contract mismatch")
    kpoints = (source / "KPOINTS").read_text(encoding="utf-8").splitlines()
    if (len(kpoints) < 5 or kpoints[2].strip().lower() != "gamma"
            or kpoints[3].split() != ["8", "8", "6"]
            or kpoints[4].split() != ["0", "0", "0"]):
        raise ValueError("image-15 electronic k-point contract mismatch")
    if "General timing and accounting informations" not in (source / "OUTCAR").read_text(
        encoding="utf-8", errors="replace",
    ):
        raise ValueError("image-15 source OUTCAR lacks completion marker")
    for item in manifest["cases"]:
        directory = work_root / "cases" / item["name"]
        if sha256(directory / "POSCAR") != item["POSCAR_sha256"]:
            raise ValueError(f"staged geometry changed: {item['name']}")
        for name in ("INCAR", "KPOINTS", "POTCAR"):
            if (directory / name).exists():
                raise FileExistsError(f"refusing to overwrite {directory / name}")
            shutil.copy2(source / name, directory / name)
        hashes = {name: sha256(directory / name) for name in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
        (directory / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8",
        )
    manifest["status"] = "inputs_finalized_no_DFT"
    manifest["source_static_directory"] = str(source)
    manifest["source_static_sha256"] = {name: sha256(source / name) for name in required}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_subparsers(dest="mode", required=True)
    stage_parser = mode.add_parser("prepare")
    stage_parser.add_argument("--trajectory", type=Path, required=True)
    stage_parser.add_argument("--summary", type=Path, required=True)
    stage_parser.add_argument("--work-root", type=Path, required=True)
    stage_parser.add_argument("--step-A", type=float, default=0.01)
    finalize_parser = mode.add_parser("finalize-source")
    finalize_parser.add_argument("--work-root", type=Path, required=True)
    finalize_parser.add_argument("--source", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "prepare":
        result = prepare(args.trajectory, args.summary, args.work_root, args.step_A)
    else:
        result = finalize_source(args.work_root, args.source)
    print(json.dumps({"status": result["status"], "cases": len(result["cases"])}))


if __name__ == "__main__":
    main()
