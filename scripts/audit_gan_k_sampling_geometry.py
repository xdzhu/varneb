"""Measure reciprocal-space sampling along the archived GaN VASP chain.

Integer k meshes need not be identical for different cells. This audit only
measures the geometric density of the *declared* mesh; it does not establish
k-point convergence or audit raw KPOINTS for every image.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sampling_metrics(cells: np.ndarray, meshes: np.ndarray) -> dict:
    """Return BZ volume per point and reciprocal-basis step lengths in Å units.

    Cell vectors are rows; reciprocal vectors include the conventional 2π.
    For a skewed cell, the three step lengths are not nearest-neighbor
    distances in the full grid and must not be interpreted as such.
    """
    cells = np.asarray(cells, dtype=float)
    meshes = np.asarray(meshes)
    if cells.ndim != 3 or cells.shape[1:] != (3, 3) or meshes.shape != (len(cells), 3):
        raise ValueError("expected N cells of shape 3x3 and N three-axis meshes")
    if (not np.isfinite(cells).all() or not np.isfinite(meshes).all()
            or not np.all(meshes == np.rint(meshes)) or np.any(meshes <= 0)):
        raise ValueError("cell and mesh entries must be finite, with positive integer meshes")
    volumes = np.linalg.det(cells)
    if np.any(volumes <= 0):
        raise ValueError("cell volumes must be positive")
    reciprocal = 2 * np.pi * np.linalg.inv(cells).transpose(0, 2, 1)
    steps = np.linalg.norm(reciprocal, axis=2) / meshes
    bz_per_point = (2 * np.pi) ** 3 / (volumes * np.prod(meshes, axis=1))
    return {
        "cell_volume_A3": volumes.tolist(),
        "bz_volume_per_kpoint_A_minus3": bz_per_point.tolist(),
        "reciprocal_basis_step_lengths_A_minus1": steps.tolist(),
    }


def audit(trajectory: Path, chain_audit: Path, mesh: tuple[int, int, int]) -> dict:
    if mesh != (8, 8, 6):
        raise ValueError("archived GaN VASP production contract declares Γ 8x8x6")
    archived = json.loads(chain_audit.read_text(encoding="utf-8"))
    expected = archived["backends"]["vasp"]["compact_final_chain_sha256"]
    if sha256(trajectory) != expected:
        raise ValueError("VASP trajectory differs from archived final-chain audit")
    atoms = read(trajectory, index=":")
    if len(atoms) != archived["n_images_total"] or len({len(a) for a in atoms}) != 1:
        raise ValueError("image count or atom count differs from archived chain")
    meshes = np.tile(np.asarray(mesh, dtype=int), (len(atoms), 1))
    data = sampling_metrics(np.array([a.cell.array for a in atoms]), meshes)
    density = np.array(data["bz_volume_per_kpoint_A_minus3"])
    steps = np.array(data["reciprocal_basis_step_lengths_A_minus1"])
    adjacent = np.abs(np.diff(density)) / density[:-1]
    worst = int(np.argmax(adjacent))
    return {
        "kind": "gan_45p7_vasp_chain_reciprocal_sampling_geometry_not_k_convergence",
        "n_images_total": len(atoms),
        "declared_gamma_centered_mesh": list(mesh),
        "mesh_provenance_limit": "Declared production mesh, not independently parsed from 29 raw KPOINTS files",
        "definition": "BZ volume per point = (2*pi)^3 / (cell volume * product(mesh)); reciprocal step lengths = |b_i|/mesh_i, including 2*pi",
        "bz_volume_per_kpoint_min_max_A_minus3": [float(density.min()), float(density.max())],
        "bz_volume_per_kpoint_endpoints_A_minus3": [float(density[0]), float(density[-1])],
        "max_adjacent_fractional_change": float(adjacent[worst]),
        "max_adjacent_image_pair_zero_based": [worst, worst + 1],
        "reciprocal_step_min_max_by_axis_A_minus1": [steps.min(axis=0).tolist(), steps.max(axis=0).tolist()],
        "per_image": data,
        "source_sha256": {"compact_final_trajectory": expected, "chain_audit": sha256(chain_audit)},
        "limitations": [
            "A smooth reciprocal-density sequence does not prove k-point convergence of energy, force, or stress.",
            "An integer mesh change may be acceptable when comparable density and property continuity are verified.",
            "FFT grid and plane-wave-basis changes are distinct from electronic k-point changes.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--chain-audit", type=Path, required=True)
    parser.add_argument("--mesh", type=int, nargs=3, default=(8, 8, 6))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.trajectory, args.chain_audit, tuple(args.mesh))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "bz_volume_per_kpoint_min_max_A_minus3", "max_adjacent_fractional_change",
        "max_adjacent_image_pair_zero_based")}))


if __name__ == "__main__":
    main()
