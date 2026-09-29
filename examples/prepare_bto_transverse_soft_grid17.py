"""Prepare only the 208 missing nodes of a nested 17x17 BTO frozen-mode cut.

The existing, independently audited 9x9 grid supplies the even/even nodes.
No calculator is constructed here.  The eight output shards are inert input
manifests for the already reviewed 100-Ry/10-au-DZP ABACUS static runner.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.bto_q1q2_reference import (
    bto_transverse_soft_plane,
    load_bto_q1q2_reference,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("report", "reference", "force-constants", "phonopy-eigenpairs",
                 "prior-points", "prior-audit", "output-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--shards", type=int, default=8)
    return parser.parse_args()


def prepare(args: argparse.Namespace) -> dict:
    if args.output_dir.exists():
        raise FileExistsError(args.output_dir)
    if args.shards < 1 or 208 % args.shards:
        raise ValueError("shards must be a positive divisor of 208")
    loaded = load_bto_q1q2_reference(
        args.report, args.reference, args.force_constants,
        args.phonopy_eigenpairs,
    )
    plane = bto_transverse_soft_plane(loaded)
    old_path = args.prior_points / "manifest.json"
    old = json.loads(old_path.read_text(encoding="utf-8"))
    audit = json.loads(args.prior_audit.read_text(encoding="utf-8"))
    if (old.get("kind") != "inert_frozen_cubic_cell_transverse_soft_mode_points_not_T_to_C_barrier"
            or old.get("n_points") != 81
            or old.get("input_sha256") != loaded.source_hashes
            or old.get("reference_id") != plane.reference_id
            or not np.allclose(old.get("axis_mode_weights"), plane.axis_weights,
                               rtol=0, atol=1e-12)
            or audit.get("n_total_DFT_points") != 81
            or audit.get("status") != "22_new_nodes_raw_energy_force_stress_audited"
            or audit.get("source_sha256", {}).get("dense_manifest") != sha256(old_path)
            or len(audit.get("samples", [])) != 81):
        raise ValueError("the audited 9x9 source or mode contract changed")
    old_by_index = {tuple(point["grid_index_q1_q2"]): point for point in old["points"]}
    old_energy = {(int(row["dense_i_q1"]), int(row["dense_j_q2"])): row
                  for row in audit["samples"]}
    if len(old_by_index) != 81 or len(old_energy) != 81:
        raise ValueError("9x9 source contains duplicate indices")

    values_q1 = [round(-1.2 + 0.15 * i, 12) for i in range(17)]
    values_q2 = [round(-0.6 + 0.075 * j, 12) for j in range(17)]
    new_points = []
    minimum_distance = float("inf")
    for j, q2 in enumerate(values_q2):
        for i, q1 in enumerate(values_q1):
            q = np.array([q1, q2], dtype=float)
            atomic = plane.frozen_coordinates(q)
            atoms = loaded.chart.to_atoms(np.r_[atomic, np.zeros(6)])
            if not np.allclose(plane.project(atomic), q, rtol=0, atol=1e-9):
                raise ValueError(f"mode-coordinate reconstruction failed at {(i, j)}")
            distances = atoms.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            separation = float(np.min(distances))
            minimum_distance = min(minimum_distance, separation)
            if separation < 1.6:
                raise ValueError(f"unsafe separation {separation:.6f} A at {(i, j)}")
            if i % 2 == 0 and j % 2 == 0:
                old_index = (i // 2, j // 2)
                prior = old_by_index[old_index]
                measured = old_energy[old_index]
                source = args.prior_points / prior["structure"]
                old_atoms = read(source)
                if (sha256(source) != prior["structure_sha256"]
                        or prior["structure_sha256"] != measured["source_structure_sha256"]
                        or abs(float(prior["q1"]) - q1) > 1e-12
                        or abs(float(prior["q2"]) - q2) > 1e-12
                        or not np.allclose(old_atoms.cell.array, atoms.cell.array,
                                           rtol=0, atol=1e-9)
                        or not np.allclose(old_atoms.get_positions(), atoms.get_positions(),
                                           rtol=0, atol=1e-9)
                        or not np.isfinite(float(measured["energy_eV_per_BTO"]))):
                    raise ValueError(f"old 9x9 node is not reusable at {(i, j)}")
                continue
            new_points.append((i, j, q1, q2, separation, atoms))
    if len(new_points) != 208:
        raise RuntimeError(f"expected 208 new points, found {len(new_points)}")

    args.output_dir.mkdir(parents=True)
    shard_rows: list[list[dict]] = [[] for _ in range(args.shards)]
    for index, (i, j, q1, q2, separation, atoms) in enumerate(new_points):
        shard = index % args.shards
        name = f"grid-q1-{i:02d}-q2-{j:02d}"
        folder = args.output_dir / f"shard-{shard:02d}" / name
        folder.mkdir(parents=True)
        structure = folder / "POSCAR"
        write(structure, atoms, format="vasp", direct=True, vasp5=True)
        reread = read(structure)
        if (not np.allclose(reread.cell.array, atoms.cell.array, atol=1e-9, rtol=0)
                or not np.allclose(reread.get_positions(), atoms.get_positions(),
                                   atol=1e-9, rtol=0)):
            raise ValueError(f"POSCAR roundtrip changed geometry at {(i, j)}")
        shard_rows[shard].append({
            "name": name, "path_image_index": None,
            "grid_index_q1_q2": [i, j], "q1": q1, "q2": q2,
            "minimum_distance_A": separation,
            "structure": f"{name}/POSCAR",
            "structure_sha256": sha256(structure),
        })
    hashes = {"old_manifest": sha256(old_path),
              "old_audit": sha256(args.prior_audit),
              **loaded.source_hashes}
    shard_hashes = {}
    for shard, rows in enumerate(shard_rows):
        folder = args.output_dir / f"shard-{shard:02d}"
        manifest = {
            "kind": "inert_frozen_cubic_cell_transverse_soft_mode_points_not_T_to_C_barrier",
            "status": "inputs_finalized_no_DFT",
            "n_points": len(rows), "n_atoms": loaded.chart.n_atoms,
            "amplitude_unit": plane.amplitude_unit,
            "axis_labels": list(plane.axis_labels),
            "axis_definition": "Ti_z-minus-Ba_z_and_Ti_x-minus-Ba_x_projected_to_Gamma_unstable_triplet",
            "axis_mode_weights": plane.axis_weights.tolist(),
            "reference_id": plane.reference_id,
            "boundary_condition": "all_unselected_atomic_modes_and_all_cell_strains_frozen_at_cubic_reference",
            "grid": None, "full_grid_shape": [17, 17],
            "new_only": True, "no_dft_launched": True,
            "source_sha256": hashes, "points": rows,
        }
        path = folder / "manifest.json"
        path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        shard_hashes[folder.name] = sha256(path)
    plan = {
        "status": "208_new_points_prepared_no_DFT",
        "purpose": "nested_17x17_BTO_frozen_soft_soft_surface",
        "claim_limit": "frozen cubic-cell cut; not a conditional PES or T-to-C barrier",
        "old_measured_points_reused": 81,
        "new_points": 208, "full_grid_shape": [17, 17],
        "q1": values_q1, "q2": values_q2,
        "minimum_distance_A": minimum_distance,
        "source_sha256": hashes,
        "shard_manifest_sha256": shard_hashes,
    }
    (args.output_dir / "plan.json").write_text(
        json.dumps(plan, indent=2) + "\n", encoding="utf-8"
    )
    return plan


def main() -> None:
    plan = prepare(parse_args())
    print(json.dumps({key: plan[key] for key in
                      ("status", "old_measured_points_reused", "new_points",
                       "minimum_distance_A", "shard_manifest_sha256")}, indent=2))


if __name__ == "__main__":
    main()
