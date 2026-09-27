"""Screen same-contract GaN off-path anchor geometries without launching DFT.

The warning window is an empirical case-specific VASP Bravais-risk screen,
not a symmetry proof. Electronic inputs are deliberately absent here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.analyze_gan_600eV_path_tube_feasibility import transported_normals
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight
from vcneb.joint_curvature import JointCurvatureCoordinates


CELL_SCALE_A = 3.3982714330050063
AMPLITUDES_A = (0.015, 0.0125, 0.010)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def screen(trajectory: Path, hessian_npz: Path, feasibility_path: Path) -> dict:
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    if (feasibility.get("status")
            != "GaN_600eV_path_adapted_joint_normal_tube_geometry_only"
            or feasibility["source_sha256"]["trajectory"] != sha256(trajectory)
            or feasibility["source_sha256"]["hessian_npz"] != sha256(hessian_npz)):
        raise ValueError("600-eV path-tube sources changed")
    frames = read(trajectory, index=":")
    if len(frames) != 29:
        raise ValueError("expected the archived 29-image GaN chain")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        seed = np.asarray(archive["eigenvectors"][:, 1], dtype=float)
    normals, transport = transported_normals(frames, seed, CELL_SCALE_A)
    for key, value in transport.items():
        if not np.isclose(value, feasibility["transport"][key], atol=1e-10, rtol=0):
            raise ValueError("normal transport differs from the first 600-eV tube")
    rows = []
    for image_index in range(1, 28):
        chart = JointCurvatureCoordinates(frames[image_index], CELL_SCALE_A)
        for amplitude in AMPLITUDES_A:
            sides = {}
            for sign, side in ((-1, "minus"), (1, "plus")):
                atoms = chart.displaced(sign * amplitude * normals[image_index])
                sides[side] = geometry_preflight(atoms)
            rows.append({
                "image_index": image_index,
                "amplitude_A": amplitude,
                "both_sides_pass_empirical_screen": not any(
                    side["empirical_near_symmetry_warning"] for side in sides.values()
                ),
                "sides": sides,
            })
    return {
        "status": "GaN_600eV_offpath_refinement_geometry_only_no_DFT",
        "claim_limit": "Empirical geometry warning does not guarantee VASP parser success; no electronic parameters changed",
        "pressure_GPa": 45.7,
        "amplitudes_A": list(AMPLITUDES_A),
        "source_sha256": {
            "trajectory": sha256(trajectory),
            "hessian_npz": sha256(hessian_npz),
            "feasibility": sha256(feasibility_path),
            "screen_script": sha256(Path(__file__)),
        },
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--hessian-npz", type=Path, required=True)
    parser.add_argument("--feasibility", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = screen(args.trajectory, args.hessian_npz, args.feasibility)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "pass_by_amplitude": {
            str(amplitude): [row["image_index"] for row in report["rows"]
                             if row["amplitude_A"] == amplitude
                             and row["both_sides_pass_empirical_screen"]]
            for amplitude in AMPLITUDES_A
        },
    }))


if __name__ == "__main__":
    main()
