"""Stage two cautious 1000 eV Newton probes from GaN VCNEB image 15.

This is a local saddle-stationarity experiment, not a replacement path or a
transition-state certificate. The 600 eV joint Hessian supplies an initial
Jacobian; 1000 eV VASP statics decide whether its predicted correction works.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write
from ase.units import GPa


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same_geometry(a, b, tolerance: float = 1e-7) -> bool:
    delta = a.get_scaled_positions(wrap=False) - b.get_scaled_positions(wrap=False)
    delta -= np.rint(delta)
    return (a.get_chemical_symbols() == b.get_chemical_symbols()
            and np.allclose(a.cell.array, b.cell.array, rtol=0, atol=tolerance)
            and np.allclose(delta, 0, rtol=0, atol=tolerance))


def prepare(trajectory: Path, summary: Path, hessian_archive: Path,
            center_outcar: Path, work_root: Path) -> dict:
    from vcneb.joint_curvature import JointCurvatureCoordinates

    if work_root.exists():
        raise FileExistsError(work_root)
    metadata = json.loads(summary.read_text(encoding="utf-8"))
    if (metadata.get("status") != "completed" or not metadata.get("converged")
            or metadata.get("n_images") != 29):
        raise ValueError("source VCNEB chain is not the audited 29-image result")
    reference = read(trajectory, index=15)
    center = read(center_outcar)
    if not same_geometry(reference, center, tolerance=2e-5):
        raise ValueError("1000 eV center is not the same image-15 geometry")
    scale = float(metadata["path_diagnostics"]["cell_scale_A"])
    coordinates = JointCurvatureCoordinates(reference, scale)
    with np.load(hessian_archive, allow_pickle=False) as archived:
        eigenvalues = archived["eigenvalues"].copy()
        eigenvectors = archived["eigenvectors"].copy()
    if (eigenvalues.shape != (15,) or eigenvectors.shape != (18, 15)
            or int(np.count_nonzero(eigenvalues < -0.1)) != 1
            or eigenvalues[0] > -3.5 or eigenvalues[1] < 2.0
            or not np.allclose(eigenvectors.T @ eigenvectors, np.eye(15), atol=1e-8)):
        raise ValueError("source joint Hessian has no reliable isolated negative direction")
    pressure = 45.7 * GPa
    gradient = coordinates.enthalpy_gradient(
        center, center.get_forces(), center.get_stress(voigt=False), pressure,
    )
    correction = -eigenvectors @ ((eigenvectors.T @ gradient) / eigenvalues)
    atom_max = float(np.max(np.linalg.norm(correction[:12].reshape(4, 3), axis=1)))
    if (not np.isfinite(correction).all() or np.linalg.norm(correction) > 0.025
            or atom_max > 0.01 or np.linalg.norm(correction[12:]) > 0.02):
        raise ValueError("predicted Newton correction exceeds the prespecified trust radius")
    work_root.mkdir(parents=True)
    (work_root / "cases").mkdir()
    (work_root / "logs").mkdir()
    write(work_root / "center_POSCAR", reference, format="vasp", direct=True, sort=False)
    records = []
    for factor, name in ((0.5, "newton_half"), (1.0, "newton_full")):
        displaced = coordinates.displaced(factor * correction)
        distances = displaced.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        if len(displaced) != 4 or np.min(distances) < 1.4:
            raise ValueError(f"unsafe trial geometry: {name}")
        directory = work_root / "cases" / name
        directory.mkdir()
        write(directory / "POSCAR", displaced, format="vasp", direct=True, sort=False)
        if not same_geometry(displaced, read(directory / "POSCAR", format="vasp")):
            raise ValueError("Newton trial POSCAR roundtrip changed the geometry")
        records.append({
            "name": name, "step_fraction": factor,
            "POSCAR_sha256": sha256(directory / "POSCAR"),
            "volume_A3": displaced.get_volume(),
            "minimum_distance_A": float(np.min(distances)),
        })
    manifest = {
        "purpose": "GaN_image15_1000eV_two_step_local_Newton_stationarity_probe",
        "status": "geometry_staged_no_DFT",
        "n_atoms": 4,
        "pressure_GPa": 45.7,
        "cell_scale_A": scale,
        "source_hessian_encut_eV": 600,
        "evaluation_encut_eV": 1000,
        "center_gradient_translation_free_eV_per_A": float(np.linalg.norm(
            coordinates.translation_free_basis().T @ gradient
        )),
        "correction_joint_norm_A": float(np.linalg.norm(correction)),
        "correction_max_atom_A": atom_max,
        "correction_cell_scaled_norm_A": float(np.linalg.norm(correction[12:])),
        "predicted_full_step_gradient_eV_per_A": float(np.linalg.norm(
            coordinates.translation_free_basis().T @ (
                gradient + eigenvectors @ (eigenvalues * (eigenvectors.T @ correction))
            )
        )),
        "correction_A": correction.tolist(),
        "center_POSCAR_sha256": sha256(work_root / "center_POSCAR"),
        "source_sha256": {
            "trajectory": sha256(trajectory), "summary": sha256(summary),
            "hessian_archive": sha256(hessian_archive),
            "center_1000_OUTCAR": sha256(center_outcar),
            "preparer": sha256(Path(__file__)),
            "joint_curvature": sha256(Path(__file__).resolve().parents[1] / "vcneb/joint_curvature.py"),
        },
        "cases": records,
        "limitations": "A small Newton probe only; no 1000 eV Hessian or basin-connection certificate.",
    }
    (work_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8",
    )
    return manifest


def finalize_source(work_root: Path, source_center: Path) -> dict:
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "geometry_staged_no_DFT":
        raise ValueError("Newton probe is not an unfinalized geometry stage")
    if (sha256(source_center / "OUTCAR") != manifest["source_sha256"]["center_1000_OUTCAR"]
            or not same_geometry(read(source_center / "POSCAR", format="vasp"),
                                 read(work_root / "center_POSCAR", format="vasp"))):
        raise ValueError("1000 eV source differs from the audited center")
    incar = (source_center / "INCAR").read_text(encoding="utf-8")
    if "ENCUT = 1000.000000" not in incar or "ISYM = -1" not in incar:
        raise ValueError("source does not have the isolated 1000 eV VASP input")
    for case in manifest["cases"]:
        directory = work_root / "cases" / case["name"]
        if sha256(directory / "POSCAR") != case["POSCAR_sha256"]:
            raise ValueError("Newton trial geometry changed before VASP staging")
        for filename in ("INCAR", "KPOINTS", "POTCAR"):
            if (directory / filename).exists():
                raise FileExistsError(directory / filename)
            shutil.copy2(source_center / filename, directory / filename)
        hashes = {filename: sha256(directory / filename)
                  for filename in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
        (directory / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8",
        )
        case["input_sha256"] = hashes
    manifest["status"] = "inputs_finalized_no_DFT"
    manifest["source_static_directory"] = str(source_center)
    manifest["source_static_sha256"] = {
        filename: sha256(source_center / filename)
        for filename in ("POSCAR", "INCAR", "KPOINTS", "POTCAR", "OUTCAR")
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_subparsers(dest="mode", required=True)
    first = mode.add_parser("prepare")
    first.add_argument("--trajectory", type=Path, required=True)
    first.add_argument("--summary", type=Path, required=True)
    first.add_argument("--hessian", type=Path, required=True)
    first.add_argument("--center-outcar", type=Path, required=True)
    first.add_argument("--work-root", type=Path, required=True)
    final = mode.add_parser("finalize-source")
    final.add_argument("--work-root", type=Path, required=True)
    final.add_argument("--source-center", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "prepare":
        result = prepare(args.trajectory, args.summary, args.hessian,
                         args.center_outcar, args.work_root)
    else:
        result = finalize_source(args.work_root, args.source_center)
    print(json.dumps({"status": result["status"], "n_cases": len(result["cases"])}))


if __name__ == "__main__":
    main()
