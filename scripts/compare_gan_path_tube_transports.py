"""Compare fixed-seed projection with recursive parallel transport for GaN.

Both are local normal-coordinate conventions on the same 600-eV path. This
read-only test uses archived path forces/stress to decide whether a costly
second off-path DFT grid is justified; no new electronic calculation occurs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.analyze_gan_600eV_path_tube_feasibility import transported_normals
from scripts.analyze_gan_ts_mode_path_overlap import joint_delta
from vcneb.joint_curvature import JointCurvatureCoordinates


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def local_tangents(frames: list, scale: float) -> np.ndarray:
    tangents = []
    for i, frame in enumerate(frames):
        chart = JointCurvatureCoordinates(frame, scale)
        left = np.zeros(18) if i == 0 else joint_delta(frame, frames[i - 1], scale)[0]
        right = np.zeros(18) if i == 28 else joint_delta(frame, frames[i + 1], scale)[0]
        free = chart.translation_free_basis()
        t = free @ (free.T @ (right - left))
        t /= np.linalg.norm(t)
        tangents.append(t)
    return np.asarray(tangents)


def parallel_transport(frames: list, seed: np.ndarray, scale: float) -> np.ndarray:
    tangent = local_tangents(frames, scale)
    normals = np.zeros((29, 18))
    center = seed - tangent[15] * float(np.dot(seed, tangent[15]))
    normals[15] = center / np.linalg.norm(center)
    for order in (range(16, 29), range(14, -1, -1)):
        previous = 15
        for i in order:
            value = normals[previous] - tangent[i] * float(np.dot(normals[previous], tangent[i]))
            norm = np.linalg.norm(value)
            if norm < 0.2:
                raise ValueError(f"recursive transverse mode collapses at image {i}")
            normals[i] = value / norm
            previous = i
    return normals


def summarize(frames: list, normals: np.ndarray, scale: float) -> dict:
    slopes = []
    for i, frame in enumerate(frames):
        chart = JointCurvatureCoordinates(frame, scale)
        gradient = chart.enthalpy_gradient(
            frame, frame.get_forces(), frame.get_stress(voigt=False), 45.7 * GPa
        )
        slopes.append(float(np.dot(gradient, normals[i])))
    overlaps = np.sum(normals[:-1] * normals[1:], axis=1)
    return {
        "transverse_slopes_eV_per_A": slopes,
        "max_adjacent_slope_jump_eV_per_A": float(np.max(np.abs(np.diff(slopes)))),
        "min_adjacent_normal_overlap": float(np.min(overlaps)),
        "overlap_03_04": float(overlaps[3]),
        "overlap_19_20": float(overlaps[19]),
        "slope_17_20": slopes[17:21],
    }


def compare(trajectory: Path, hessian_npz: Path) -> dict:
    frames = read(trajectory, index=":")
    if len(frames) != 29:
        raise ValueError("expected 29-image GaN path")
    with np.load(hessian_npz, allow_pickle=False) as data:
        seed = np.asarray(data["eigenvectors"][:, 1], dtype=float)
    scale = 3.3982714330050063
    projected, _ = transported_normals(frames, seed, scale)
    parallel = parallel_transport(frames, seed, scale)
    return {
        "status": "GaN_600eV_path_normal_transport_offline_comparison",
        "fixed_seed_projection": summarize(frames, projected, scale),
        "recursive_parallel_transport": summarize(frames, parallel, scale),
        "source_sha256": {
            "trajectory": sha256(trajectory),
            "hessian_npz": sha256(hessian_npz),
            "comparator": sha256(Path(__file__)),
        },
        "limitation": "Gradient smoothness alone does not determine the physical off-path enthalpy or guarantee a good tube.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trajectory", "hessian-npz", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = compare(args.trajectory, args.hessian_npz)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "old": result["fixed_seed_projection"]["max_adjacent_slope_jump_eV_per_A"],
        "new": result["recursive_parallel_transport"]["max_adjacent_slope_jump_eV_per_A"],
        "old_min_overlap": result["fixed_seed_projection"]["min_adjacent_normal_overlap"],
        "new_min_overlap": result["recursive_parallel_transport"]["min_adjacent_normal_overlap"],
    }))


if __name__ == "__main__":
    main()
