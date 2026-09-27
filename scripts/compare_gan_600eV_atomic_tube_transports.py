"""Compare two atomic-transverse transports on the archived GaN VCNEB chain.

This is an offline coordinate-gauge check. A smooth transported normal alone
does not establish a smooth potential surface or a unique phonon identity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.analyze_gan_600eV_atomic_tube_feasibility import (
    CELL_SCALE_A,
    atomic_normals,
)
from scripts.analyze_gan_ts_mode_path_overlap import joint_delta


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_tangents(frames: list) -> np.ndarray:
    tangents = []
    for index, frame in enumerate(frames):
        left = np.zeros(18) if index == 0 else joint_delta(frame, frames[index - 1], CELL_SCALE_A)[0]
        right = np.zeros(18) if index == 28 else joint_delta(frame, frames[index + 1], CELL_SCALE_A)[0]
        vector = (right - left)[:12].reshape(4, 3)
        vector -= np.mean(vector, axis=0)
        length = float(np.linalg.norm(vector))
        if length < 1e-8:
            raise ValueError(f"atomic tangent vanishes at image {index}")
        tangents.append(vector.ravel() / length)
    return np.asarray(tangents)


def parallel_transport_from_center(tangents: np.ndarray, center_normal: np.ndarray) -> np.ndarray:
    if tangents.shape != (29, 12) or center_normal.shape != (12,):
        raise ValueError("unexpected atomic transverse chart dimensions")
    result = np.zeros_like(tangents)
    result[15] = center_normal
    for indices in (range(16, 29), range(14, -1, -1)):
        for index in indices:
            previous = index - 1 if index > 15 else index + 1
            old_tangent = tangents[previous]
            new_tangent = tangents[index]
            tangent_overlap = float(np.dot(old_tangent, new_tangent))
            if tangent_overlap <= -0.95:
                raise ValueError(f"atomic tangent nearly reverses at image {index}")
            # Minimal rotation between neighboring tangent hyperplanes. This
            # preserves the normal norm even when simple projection collapses.
            normal = result[previous] - (
                float(np.dot(result[previous], new_tangent)) / (1 + tangent_overlap)
            ) * (old_tangent + new_tangent)
            length = float(np.linalg.norm(normal))
            if not np.isclose(length, 1, atol=1e-8):
                raise ValueError(f"minimal-rotation atomic normal changed norm at image {index}")
            result[index] = normal / length
    if not np.allclose(np.sum(result * tangents, axis=1), 0, atol=1e-10):
        raise ValueError("recursive normals lost tangent orthogonality")
    return result


def diagnostics(frames: list, normals: np.ndarray) -> dict:
    overlaps = np.sum(normals[:-1] * normals[1:], axis=1)
    slopes = [-float(np.dot(frame.get_forces().ravel(), normals[index]))
              for index, frame in enumerate(frames)]
    jumps = np.abs(np.diff(slopes))
    return {
        "minimum_neighbor_normal_overlap": float(np.min(overlaps)),
        "minimum_central_5_to_22_overlap": float(np.min(overlaps[5:22])),
        "maximum_adjacent_transverse_slope_jump_eV_per_A": float(np.max(jumps)),
        "maximum_central_5_to_22_slope_jump_eV_per_A": float(np.max(jumps[5:22])),
        "root_mean_square_adjacent_slope_jump_eV_per_A": float(np.sqrt(np.mean(jumps**2))),
        "transverse_slopes_eV_per_A": slopes,
        "neighbor_normal_overlaps": overlaps.tolist(),
    }


def analyze(trajectory: Path, feasibility_path: Path, hessian_npz: Path) -> dict:
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    if (feasibility.get("status")
            != "GaN_600eV_atomic_transverse_tube_geometry_only_no_DFT"
            or feasibility["source_sha256"].get("trajectory") != sha256(trajectory)
            or feasibility["source_sha256"].get("hessian_npz") != sha256(hessian_npz)):
        raise ValueError("atomic chart source changed")
    frames = read(trajectory, index=":")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        seed = np.asarray(archive["eigenvectors"][:, 1], dtype=float)
    fixed, pinned = atomic_normals(frames, seed)
    if any(not np.isclose(value, feasibility["transport"][key], atol=1e-10, rtol=0)
           for key, value in pinned.items()):
        raise ValueError("fixed-seed atomic normal differs from audited chart")
    recursive = parallel_transport_from_center(atomic_tangents(frames), fixed[15])
    agreement = np.sum(fixed * recursive, axis=1)
    return {
        "status": "GaN_600eV_atomic_transverse_transport_offline_comparison",
        "claim_limit": "Coordinate continuity diagnostic only; not a DFT PES or mode-identity proof",
        "fixed_seed_projection": diagnostics(frames, fixed),
        "recursive_parallel_projection": diagnostics(frames, recursive),
        "fixed_vs_recursive_normal_overlap_by_image": agreement.tolist(),
        "minimum_fixed_vs_recursive_normal_overlap": float(np.min(agreement)),
        "source_sha256": {
            "trajectory": sha256(trajectory),
            "atomic_feasibility": sha256(feasibility_path),
            "hessian_npz": sha256(hessian_npz),
            "analyzer": sha256(Path(__file__)),
        },
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
        "fixed_min_overlap": report["fixed_seed_projection"]["minimum_neighbor_normal_overlap"],
        "recursive_min_overlap": report["recursive_parallel_projection"]["minimum_neighbor_normal_overlap"],
        "fixed_max_slope_jump": report["fixed_seed_projection"]["maximum_adjacent_transverse_slope_jump_eV_per_A"],
        "recursive_max_slope_jump": report["recursive_parallel_projection"]["maximum_adjacent_transverse_slope_jump_eV_per_A"],
        "minimum_mode_agreement": report["minimum_fixed_vs_recursive_normal_overlap"],
    }))


if __name__ == "__main__":
    main()
