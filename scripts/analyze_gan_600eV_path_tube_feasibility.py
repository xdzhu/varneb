"""Preflight a GaN path-adapted 2D enthalpy tube at the 600-eV contract.

The reaction coordinate is the audited VCNEB joint arc length, not a phonon.
A same-protocol local stable joint Hessian direction is projected normal to
each image's path tangent. This script does geometry only; it launches no DFT.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.analyze_gan_ts_mode_path_overlap import joint_delta
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight
from vcneb.joint_curvature import JointCurvatureCoordinates


ANCHORS = (0, 5, 10, 13, 15, 17, 20, 25, 28)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def transported_normals(frames: list, seed: np.ndarray, scale: float) -> tuple[np.ndarray, dict]:
    """Project one joint direction normal to each local image tangent."""

    if len(frames) != 29 or np.asarray(seed).shape != (18,):
        raise ValueError("expected a 29-image, four-atom path and 18-D seed")
    normals = []
    overlaps = []
    rotations = []
    for i, frame in enumerate(frames):
        chart = JointCurvatureCoordinates(frame, scale)
        left = np.zeros(18) if i == 0 else joint_delta(frame, frames[i - 1], scale)[0]
        right = np.zeros(18) if i == 28 else joint_delta(frame, frames[i + 1], scale)[0]
        tangent = right - left
        free = chart.translation_free_basis()
        tangent = free @ (free.T @ tangent)
        tangent /= np.linalg.norm(tangent)
        normal = free @ (free.T @ seed)
        normal -= tangent * float(np.dot(normal, tangent))
        retained = float(np.linalg.norm(normal))
        if not np.isfinite(retained) or retained < 0.2:
            raise ValueError(f"transverse seed collapses into the tangent at image {i}")
        normal /= retained
        if normals and np.dot(normal, normals[-1]) < 0:
            normal *= -1
        normals.append(normal)
        overlaps.append(retained)
        if i < 28:
            rotations.append(joint_delta(frame, frames[i + 1], scale)[1])
    result = np.asarray(normals)
    if np.dot(result[15], seed) < 0:
        result *= -1
    return result, {
        "minimum_transverse_seed_retention": min(overlaps),
        "maximum_neighbor_rotation_frobenius": max(rotations),
        "minimum_neighbor_normal_overlap": float(np.min(np.sum(result[1:] * result[:-1], axis=1))),
    }


def analyze(trajectory: Path, chain_audit: Path, hessian_root: Path) -> dict:
    chain = json.loads(chain_audit.read_text(encoding="utf-8"))
    hessian_audit_path = hessian_root / "audit.json"
    hessian_npz = hessian_root / "joint_hessian.npz"
    curvature = json.loads(hessian_audit_path.read_text(encoding="utf-8"))
    if (chain.get("backends", {}).get("vasp", {}).get("compact_final_chain_sha256")
            != sha256(trajectory)
            or curvature.get("status")
            != "GaN_600eV_near_stationary_joint_hessian_one_step_not_TS_certificate"
            or curvature.get("negative_count_below_minus_0p1_eV_per_A2") != 1):
        raise ValueError("path or local mode provenance is not the audited 600-eV case")
    frames = read(trajectory, index=":")
    if (len(frames) != 29
            or any(frame.get_chemical_symbols() != ["Ga", "Ga", "N", "N"] for frame in frames)):
        raise ValueError("GaN path images or atom mapping changed")
    with np.load(hessian_npz, allow_pickle=False) as data:
        eigenvalues = np.asarray(data["eigenvalues"], dtype=float)
        modes = np.asarray(data["eigenvectors"], dtype=float)
    if (modes.shape != (18, 15)
            or not np.allclose(eigenvalues, curvature["eigenvalues_eV_per_A2"], atol=1e-8)
            or not eigenvalues[0] < 0 < eigenvalues[1]):
        raise ValueError("local 600-eV stable transverse mode changed")
    scale = 3.3982714330050063
    normals, transport = transported_normals(frames, modes[:, 1], scale)
    arc_steps = [float(np.linalg.norm(joint_delta(frames[i], frames[i + 1], scale)[0]))
                 for i in range(28)]
    arc = np.r_[0.0, np.cumsum(arc_steps)]
    arc /= arc[-1]
    pressure = 45.7 * GPa
    enthalpy = [float(frame.get_potential_energy() + pressure * frame.get_volume())
                for frame in frames]
    trials = []
    for amplitude in (0.02, 0.015, 0.01, 0.0075):
        cases = []
        for index in ANCHORS:
            chart = JointCurvatureCoordinates(frames[index], scale)
            for sign in (-1, 1):
                atoms = chart.displaced(sign * amplitude * normals[index])
                cases.append({"image_index": index, "sign": sign,
                              **geometry_preflight(atoms)})
        trials.append({
            "q_transverse_A": amplitude,
            "n_anchors": len(ANCHORS),
            "max_basal_projection_A": float(max(item["basal_projection_A"] for item in cases)),
            "min_distance_A": float(min(item["minimum_distance_A"] for item in cases)),
            "near_symmetry_warning_count": sum(item["empirical_near_symmetry_warning"] for item in cases),
            "cases": cases,
        })
    viable = next((trial for trial in trials if trial["near_symmetry_warning_count"] == 0), None)
    return {
        "status": "GaN_600eV_path_adapted_joint_normal_tube_geometry_only",
        "n_images": 29,
        "n_interior_images": 27,
        "n_offpath_static_cases_if_run": 2 * len(ANCHORS),
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "anchors": list(ANCHORS),
        "arc_fraction_s": arc.tolist(),
        "path_enthalpy_eV_per_cell": enthalpy,
        "transport": transport,
        "local_seed_eigenvalue_eV_per_A2": float(eigenvalues[1]),
        "amplitude_trials": trials,
        "largest_preflighted_amplitude_A": None if viable is None else viable["q_transverse_A"],
        "source_sha256": {
            "trajectory": sha256(trajectory), "chain_audit": sha256(chain_audit),
            "hessian_audit": sha256(hessian_audit_path), "hessian_npz": sha256(hessian_npz),
            "analyzer": sha256(Path(__file__)),
        },
        "limitations": [
            "A transported joint normal is not a fixed global phonon eigenmode.",
            "This report contains no off-path DFT energies or two-dimensional contour.",
            "Near-symmetry risk screen is empirical, not a guarantee of VASP parser success.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trajectory", "chain-audit", "hessian-root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = analyze(args.trajectory, args.chain_audit, args.hessian_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "largest_preflighted_amplitude_A": report["largest_preflighted_amplitude_A"],
        "transport": report["transport"],
    }))


if __name__ == "__main__":
    main()
