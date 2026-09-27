"""Define an atomic-only transverse chart along the variable-cell GaN path.

The reaction coordinate retains the full atomic+cell VCNEB arc. At each image,
the transverse displacement is atomic and orthogonal to the atomic part of
the local tangent, so the already-audited cell is held fixed for that point.
This is not the joint atom-strain mode used in the narrow-tube experiment.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.analyze_gan_ts_mode_path_overlap import joint_delta
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight
from vcneb.joint_curvature import JointCurvatureCoordinates


CELL_SCALE_A = 3.3982714330050063
AMPLITUDE_A = 0.2


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_normals(frames: list, joint_seed: np.ndarray) -> tuple[np.ndarray, dict]:
    if len(frames) != 29 or np.asarray(joint_seed).shape != (18,):
        raise ValueError("expected 29 GaN frames and one 18-D joint seed")
    seed = np.asarray(joint_seed[:12], dtype=float).reshape(4, 3)
    seed -= np.mean(seed, axis=0)
    if np.linalg.norm(seed) < 0.1:
        raise ValueError("joint stable mode has negligible atomic component")
    seed = seed.ravel() / np.linalg.norm(seed)
    result = []
    retained = []
    for index, frame in enumerate(frames):
        left = np.zeros(18) if index == 0 else joint_delta(frame, frames[index - 1], CELL_SCALE_A)[0]
        right = np.zeros(18) if index == 28 else joint_delta(frame, frames[index + 1], CELL_SCALE_A)[0]
        atomic_tangent = (right - left)[:12].reshape(4, 3)
        atomic_tangent -= np.mean(atomic_tangent, axis=0)
        tangent_norm = float(np.linalg.norm(atomic_tangent))
        if tangent_norm < 1e-8:
            raise ValueError(f"atomic path tangent vanishes at image {index}")
        tangent = atomic_tangent.ravel() / tangent_norm
        normal = seed - tangent * float(np.dot(seed, tangent))
        length = float(np.linalg.norm(normal))
        if length < 0.2:
            raise ValueError(f"atomic transverse seed collapses at image {index}")
        normal /= length
        if result and np.dot(normal, result[-1]) < 0:
            normal *= -1
        result.append(normal)
        retained.append(length)
    normals = np.asarray(result)
    if np.dot(normals[15], seed) < 0:
        normals *= -1
    return normals, {
        "minimum_seed_retention": float(min(retained)),
        "minimum_adjacent_normal_overlap": float(np.min(np.sum(normals[1:] * normals[:-1], axis=1))),
        "atomic_seed_norm_in_joint_eigenvector": float(np.linalg.norm(joint_seed[:12])),
    }


def analyze(trajectory: Path, feasibility_path: Path, hessian_npz: Path) -> dict:
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    if (feasibility.get("status")
            != "GaN_600eV_path_adapted_joint_normal_tube_geometry_only"
            or feasibility["source_sha256"]["trajectory"] != sha256(trajectory)
            or feasibility["source_sha256"]["hessian_npz"] != sha256(hessian_npz)):
        raise ValueError("same-600-eV path or mode source changed")
    frames = read(trajectory, index=":")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        seed = np.asarray(archive["eigenvectors"][:, 1], dtype=float)
    normals, transport = atomic_normals(frames, seed)
    cases = []
    for index, frame in enumerate(frames):
        chart = JointCurvatureCoordinates(frame, CELL_SCALE_A)
        for sign, side in ((-1, "minus"), (1, "plus")):
            delta = np.r_[sign * AMPLITUDE_A * normals[index], np.zeros(6)]
            atoms = chart.displaced(delta)
            if not np.allclose(atoms.cell.array, frame.cell.array, rtol=0, atol=1e-10):
                raise ValueError(f"atomic-only displacement changed cell at image {index}")
            geometry = geometry_preflight(atoms)
            cases.append({
                "image_index": index,
                "side": side,
                "q_atom_A": sign * AMPLITUDE_A,
                "cell_unchanged_from_audited_path": True,
                **geometry,
            })
    return {
        "status": "GaN_600eV_atomic_transverse_tube_geometry_only_no_DFT",
        "claim_limit": "Path-adapted atomic-only mode; s still includes VCNEB cell changes; not a joint transverse mode",
        "pressure_GPa": 45.7,
        "amplitude_A": AMPLITUDE_A,
        "n_images": len(frames),
        "transport": transport,
        "minimum_distance_A": min(row["minimum_distance_A"] for row in cases),
        "cases": cases,
        "source_sha256": {
            "trajectory": sha256(trajectory),
            "joint_tube_feasibility": sha256(feasibility_path),
            "hessian_npz": sha256(hessian_npz),
            "analyzer": sha256(Path(__file__)),
        },
        "limitation": (
            "Keeping the audited path cell avoids new cell-classification boundaries, "
            "but atomic geometry and SCF must still be tested with the original 600-eV settings."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trajectory", "feasibility", "hessian-npz", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = analyze(args.trajectory, args.feasibility, args.hessian_npz)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "transport": report["transport"],
        "minimum_distance_A": report["minimum_distance_A"],
    }))


if __name__ == "__main__":
    main()
