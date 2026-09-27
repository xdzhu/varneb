"""Stage a bounded atomic-only GaN central-path tube at q=+/-0.05 A.

The path coordinate retains the variable-cell chain. The transverse mode is
atomic-only, remains normal to the local atomic tangent, and keeps each
already-audited image cell exactly unchanged. This is not a global mode chart.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.analyze_gan_600eV_atomic_tube_feasibility import (
    CELL_SCALE_A,
    atomic_normals,
)
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


ANCHORS = (5, 8, 11, 14, 17, 18, 20, 22)
AMPLITUDE_A = 0.05


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(trajectory: Path, atomic_feasibility_path: Path,
            joint_feasibility_path: Path, transport_path: Path,
            hessian_npz: Path, source_static: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    atomic_feasibility = json.loads(atomic_feasibility_path.read_text(encoding="utf-8"))
    joint_feasibility = json.loads(joint_feasibility_path.read_text(encoding="utf-8"))
    transport = json.loads(transport_path.read_text(encoding="utf-8"))
    if (atomic_feasibility.get("status")
            != "GaN_600eV_atomic_transverse_tube_geometry_only_no_DFT"
            or joint_feasibility.get("status")
            != "GaN_600eV_path_adapted_joint_normal_tube_geometry_only"
            or transport.get("status")
            != "GaN_600eV_atomic_transverse_transport_offline_comparison"
            or transport["fixed_seed_projection"]["minimum_central_5_to_22_overlap"] <= 0.95
            or atomic_feasibility["source_sha256"]["trajectory"] != sha256(trajectory)
            or atomic_feasibility["source_sha256"]["hessian_npz"] != sha256(hessian_npz)
            or joint_feasibility["source_sha256"]["trajectory"] != sha256(trajectory)
            or joint_feasibility["source_sha256"]["hessian_npz"] != sha256(hessian_npz)
            or transport["source_sha256"]["trajectory"] != sha256(trajectory)
            or any(sha256(source_static / name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("central atomic chart or 600-eV electronic sources changed")
    frames = read(trajectory, index=":")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        seed = np.asarray(archive["eigenvectors"][:, 1], dtype=float)
    normals, metrics = atomic_normals(frames, seed)
    if len(frames) != 29 or any(
            not np.isclose(value, atomic_feasibility["transport"][key], atol=1e-10, rtol=0)
            for key, value in metrics.items()):
        raise ValueError("atomic normal definition changed")
    proposed = []
    for index in ANCHORS:
        chart = JointCurvatureCoordinates(frames[index], CELL_SCALE_A)
        for sign, side in ((-1, "minus"), (1, "plus")):
            delta = np.r_[sign * AMPLITUDE_A * normals[index], np.zeros(6)]
            atoms = chart.displaced(delta)
            geometry = geometry_preflight(atoms)
            if not np.allclose(atoms.cell.array, frames[index].cell.array,
                               atol=1e-10, rtol=0):
                raise ValueError(f"atomic-only chart changed cell at image {index} {side}")
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
            "arc_fraction_s": joint_feasibility["arc_fraction_s"][index],
            "path_center_enthalpy_eV_per_cell":
                joint_feasibility["path_enthalpy_eV_per_cell"][index],
            "cell_unchanged_from_audited_path": True,
            **geometry,
            "input_sha256": hashes,
        })
    manifest = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_atomic_only_central_path_tube_pilot",
        "claim_limit": "Central images 5-22 only; q is atomic-only and path s remains variable-cell",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "anchors": list(ANCHORS),
        "q_atom_A": [-AMPLITUDE_A, 0.0, AMPLITUDE_A],
        "predeclared_LOO_gate_meV_per_GaN": 1.0,
        "input_contract_sha256": PRODUCTION_INPUT_SHA256,
        "source_sha256": {
            "trajectory": sha256(trajectory),
            "atomic_feasibility": sha256(atomic_feasibility_path),
            "joint_feasibility": sha256(joint_feasibility_path),
            "transport_comparison": sha256(transport_path),
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
    for name in ("trajectory", "atomic-feasibility", "joint-feasibility",
                 "transport", "hessian-npz", "source-static", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.trajectory, args.atomic_feasibility,
                     args.joint_feasibility, args.transport, args.hessian_npz,
                     args.source_static, args.output)
    print(json.dumps({"status": result["status"], "cases": len(result["cases"])}))


if __name__ == "__main__":
    main()
