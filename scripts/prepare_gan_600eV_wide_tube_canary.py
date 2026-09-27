"""Stage a bounded wider GaN transverse-tube canary at unchanged 600 eV.

These eight frozen off-path points are exploratory geometry/DFT checks, not a
smooth 2D surface or a conditional relaxation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.analyze_gan_600eV_path_tube_feasibility import transported_normals
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


ANCHORS = (5, 15, 18, 20)
AMPLITUDE_A = 0.2
CELL_SCALE_A = 3.3982714330050063


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(trajectory: Path, feasibility_path: Path, hessian_npz: Path,
            source_static: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    if (feasibility.get("status")
            != "GaN_600eV_path_adapted_joint_normal_tube_geometry_only"
            or feasibility["source_sha256"]["trajectory"] != sha256(trajectory)
            or feasibility["source_sha256"]["hessian_npz"] != sha256(hessian_npz)
            or any(sha256(source_static / name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("600-eV wide-canary sources or electronic contract changed")
    frames = read(trajectory, index=":")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        seed = np.asarray(archive["eigenvectors"][:, 1], dtype=float)
    normals, transport = transported_normals(frames, seed, CELL_SCALE_A)
    if len(frames) != 29 or any(
            not np.isclose(value, feasibility["transport"][key], atol=1e-10, rtol=0)
            for key, value in transport.items()):
        raise ValueError("the original path-adapted coordinates changed")
    proposed = []
    for index in ANCHORS:
        chart = JointCurvatureCoordinates(frames[index], CELL_SCALE_A)
        for sign, side in ((-1, "minus"), (1, "plus")):
            atoms = chart.displaced(sign * AMPLITUDE_A * normals[index])
            geometry = geometry_preflight(atoms)
            if geometry["empirical_near_symmetry_warning"]:
                raise ValueError(f"wide geometry risk at image {index} {side}")
            proposed.append((index, sign, side, atoms, geometry))
    # The full batch passes geometry preflight before the output is created.
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
            raise ValueError(f"input contract changed at {name}")
        (target / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
        )
        cases.append({
            "name": name,
            "image_index": index,
            "q_perp_A": sign * AMPLITUDE_A,
            "arc_fraction_s": feasibility["arc_fraction_s"][index],
            "path_center_enthalpy_eV_per_cell": feasibility["path_enthalpy_eV_per_cell"][index],
            **geometry,
            "input_sha256": hashes,
        })
    manifest = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_path_tube_wide_amplitude_canary",
        "claim_limit": "Eight exploratory frozen statics only; no contour or MEP claim",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "anchors": list(ANCHORS),
        "q_perp_A": [-AMPLITUDE_A, AMPLITUDE_A],
        "input_contract_sha256": PRODUCTION_INPUT_SHA256,
        "source_sha256": {
            "trajectory": sha256(trajectory),
            "feasibility": sha256(feasibility_path),
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
