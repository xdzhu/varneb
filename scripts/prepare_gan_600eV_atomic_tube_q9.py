"""Add four interleaved q values to each audited GaN central path station.

The original VASP/PBE 600-eV INCAR, KPOINTS and POTCAR are copied bytewise.
The existing 18x5 measured grid is reused.  All 72 new points are prospective
tests of a quartic q interpolant frozen before their DFT calculations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analyze_gan_600eV_atomic_tube_feasibility import (
    CELL_SCALE_A,
    atomic_normals,
)
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


Q_NEW = (-0.0375, -0.0125, 0.0125, 0.0375)
Q_FULL = (-0.05, -0.0375, -0.025, -0.0125, 0.0,
          0.0125, 0.025, 0.0375, 0.05)
IMAGE_INDICES = tuple(range(5, 23))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(args: argparse.Namespace) -> dict:
    if args.output.exists():
        raise FileExistsError(args.output)
    report = json.loads(args.dense_report.read_text(encoding="utf-8"))
    prior = json.loads(args.dense_manifest.read_text(encoding="utf-8"))
    feasibility = json.loads(args.feasibility.read_text(encoding="utf-8"))
    if (report.get("status") != "GaN_600eV_central_dense_atomic_tube_raw_audited"
            or report.get("central_image_indices") != list(IMAGE_INDICES)
            or report.get("q_atom_A") != [-0.05, -0.025, 0.0, 0.025, 0.05]
            or report.get("source_sha256", {}).get("manifest") != sha256(args.dense_manifest)
            or prior.get("status") != "inputs_finalized_no_DFT"
            or prior.get("n_new_statics") != 44
            or prior.get("source_sha256", {}).get("trajectory") != sha256(args.trajectory)
            or prior.get("source_sha256", {}).get("atomic_feasibility") != sha256(args.feasibility)
            or prior.get("source_sha256", {}).get("hessian_npz") != sha256(args.hessian_npz)
            or prior.get("source_sha256", {}).get("source_static_OUTCAR")
               != sha256(args.source_static / "OUTCAR")
            or feasibility.get("status") != "GaN_600eV_atomic_transverse_tube_geometry_only_no_DFT"
            or prior.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or any(sha256(args.source_static / name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("the 18x5 measured chart or 600-eV electronic contract changed")
    frames = read(args.trajectory, index=":")
    if len(frames) != 29:
        raise ValueError("expected the complete 29-image GaN path")
    with np.load(args.hessian_npz, allow_pickle=False) as archive:
        seed = np.asarray(archive["eigenvectors"][:, 1], dtype=float)
    normals, metrics = atomic_normals(frames, seed)
    if any(not np.isclose(value, feasibility["transport"][key], atol=1e-10, rtol=0)
           for key, value in metrics.items()):
        raise ValueError("path-adapted transverse direction changed")
    old_q = np.asarray(prior["q_atom_A"], dtype=float)
    observed = np.asarray(report["excess_enthalpy_meV_per_GaN"], dtype=float)
    if observed.shape != (18, 5) or not np.all(np.isfinite(observed)):
        raise ValueError("old measured 18x5 energy matrix is incomplete")
    if not np.array_equal(old_q, [-0.05, -0.025, 0.0, 0.025, 0.05]):
        raise ValueError("prior q coordinates changed")

    proposals = []
    for row, index in enumerate(IMAGE_INDICES):
        coefficient = np.polynomial.polynomial.polyfit(old_q, observed[row], 4)
        for q in Q_NEW:
            chart = JointCurvatureCoordinates(frames[index], CELL_SCALE_A)
            atoms = chart.displaced(np.r_[q * normals[index], np.zeros(6)])
            geometry = geometry_preflight(atoms)
            if not np.allclose(atoms.cell.array, frames[index].cell.array,
                               atol=1e-10, rtol=0):
                raise ValueError(f"off-path cell changed at image {index}")
            proposals.append((index, q, atoms, geometry, float(
                np.polynomial.polynomial.polyval(q, coefficient))))
    if len(proposals) != 72:
        raise RuntimeError("expected exactly 72 new interleaved statics")

    args.output.mkdir(parents=True)
    (args.output / "cases").mkdir()
    cases = []
    for index, q, atoms, geometry, prediction in proposals:
        side = "minus" if q < 0 else "plus"
        width = "0125" if abs(q) == 0.0125 else "0375"
        name = f"image_{index:02d}_q{width}_{side}"
        folder = args.output / "cases" / name
        folder.mkdir()
        write(folder / "POSCAR", atoms, format="vasp", direct=True, sort=False)
        if not same_geometry(atoms, read(folder / "POSCAR", format="vasp")):
            raise ValueError(f"POSCAR roundtrip failed for {name}")
        for filename in PRODUCTION_INPUT_SHA256:
            shutil.copy2(args.source_static / filename, folder / filename)
        hashes = {filename: sha256(folder / filename)
                  for filename in ("POSCAR", *PRODUCTION_INPUT_SHA256)}
        if any(hashes[filename] != digest
               for filename, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"electronic input changed for {name}")
        (folder / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
        )
        cases.append({
            "name": name, "image_index": index, "q_atom_A": q,
            "role": "prospective_interleaved_q_holdout",
            "prediction_excess_meV_per_GaN": prediction,
            "cell_unchanged_from_audited_path": True,
            "input_sha256": hashes, **geometry,
        })
    manifest = {
        "status": "72_interleaved_inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_central_18x9_frozen_atomic_transverse_cut",
        "claim_limit": "central images 5-22, fixed-cell off-path statics; not a relaxed or two-phonon PES",
        "pressure_GPa": 45.7, "formula_units_per_cell": 2,
        "central_image_indices": list(IMAGE_INDICES),
        "q_atom_A_full": list(Q_FULL), "q_atom_A_new": list(Q_NEW),
        "n_old_measured_points": 90, "n_new_statics": 72,
        "prospective_error_gate_meV_per_GaN": 1.0,
        "prospective_prediction_method": "per-image quartic interpolant of five prior measured q nodes",
        "input_contract_sha256": PRODUCTION_INPUT_SHA256,
        "source_sha256": {
            "trajectory": sha256(args.trajectory),
            "dense_report": sha256(args.dense_report),
            "dense_manifest": sha256(args.dense_manifest),
            "feasibility": sha256(args.feasibility),
            "hessian_npz": sha256(args.hessian_npz),
            "source_static_OUTCAR": sha256(args.source_static / "OUTCAR"),
            "preparer": sha256(Path(__file__)),
        },
        "cases": cases,
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trajectory", "dense-report", "dense-manifest", "feasibility",
                 "hessian-npz", "source-static", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args)
    print(json.dumps({
        "status": result["status"], "n_old_measured_points": 90,
        "n_new_statics": len(result["cases"]),
        "minimum_distance_A": min(row["minimum_distance_A"]
                                  for row in result["cases"]),
    }))


if __name__ == "__main__":
    main()
