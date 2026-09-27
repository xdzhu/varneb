"""Stage four fixed-path-cell GaN atomic-transverse canaries at 600 eV."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.analyze_gan_600eV_atomic_tube_feasibility import (
    AMPLITUDE_A,
    CELL_SCALE_A,
    atomic_normals,
)
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


ANCHORS = (5, 18)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(trajectory: Path, feasibility_path: Path, hessian_npz: Path,
            source_static: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    if (feasibility.get("status")
            != "GaN_600eV_atomic_transverse_tube_geometry_only_no_DFT"
            or feasibility["source_sha256"]["trajectory"] != sha256(trajectory)
            or feasibility["source_sha256"]["hessian_npz"] != sha256(hessian_npz)
            or feasibility["amplitude_A"] != AMPLITUDE_A
            or any(sha256(source_static / name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("atomic tube or 600-eV electronic source changed")
    frames = read(trajectory, index=":")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        seed = np.asarray(archive["eigenvectors"][:, 1], dtype=float)
    normals, transport = atomic_normals(frames, seed)
    if len(frames) != 29 or any(not np.isclose(value, feasibility["transport"][key], atol=1e-10, rtol=0)
                                for key, value in transport.items()):
        raise ValueError("atomic normal transport changed")
    proposed = []
    for index in ANCHORS:
        chart = JointCurvatureCoordinates(frames[index], CELL_SCALE_A)
        for sign, side in ((-1, "minus"), (1, "plus")):
            delta = np.r_[sign * AMPLITUDE_A * normals[index], np.zeros(6)]
            atoms = chart.displaced(delta)
            geometry = geometry_preflight(atoms)
            if (geometry["empirical_near_symmetry_warning"]
                    or not np.allclose(atoms.cell.array, frames[index].cell.array,
                                       atol=1e-10, rtol=0)):
                raise ValueError(f"atomic-only chart changed or risks cell at {index} {side}")
            proposed.append((index, sign, side, atoms, geometry))
    output.mkdir(parents=True)
    (output / "cases").mkdir()
    cases = []
    for index, sign, side, atoms, geometry in proposed:
        name = f"image_{index:02d}_q_{side}"
        target = output / "cases" / name
        target.mkdir()
        write(target / "POSCAR", atoms, format="vasp", direct=True, sort=False)
        if not same_geometry(atoms, read(target / "POSCAR", format="vasp")):
            raise ValueError(f"POSCAR roundtrip failed at {name}")
        for filename in PRODUCTION_INPUT_SHA256:
            shutil.copy2(source_static / filename, target / filename)
        hashes = {filename: sha256(target / filename)
                  for filename in ("POSCAR", *PRODUCTION_INPUT_SHA256)}
        if any(hashes[filename] != digest
               for filename, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"electronic contract changed at {name}")
        (target / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
        )
        cases.append({
            "name": name,
            "image_index": index,
            "q_atom_A": sign * AMPLITUDE_A,
            **geometry,
            "input_sha256": hashes,
        })
    manifest = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_fixed_cell_atomic_transverse_canary",
        "claim_limit": "Four atomic-only off-path statics; path reaction coordinate still variable-cell",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "anchors": list(ANCHORS),
        "q_atom_A": [-AMPLITUDE_A, AMPLITUDE_A],
        "input_contract_sha256": PRODUCTION_INPUT_SHA256,
        "source_sha256": {
            "trajectory": sha256(trajectory),
            "atomic_feasibility": sha256(feasibility_path),
            "hessian_npz": sha256(hessian_npz),
            "source_static_OUTCAR": sha256(source_static / "OUTCAR"),
            "preparer": sha256(Path(__file__)),
        },
        "cases": cases,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trajectory", "feasibility", "hessian-npz", "source-static", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.trajectory, args.feasibility, args.hessian_npz,
                     args.source_static, args.output)
    print(json.dumps({"status": result["status"], "cases": len(result["cases"])}))


if __name__ == "__main__":
    main()
