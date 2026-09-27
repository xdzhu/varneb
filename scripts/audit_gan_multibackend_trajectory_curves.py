"""Rebuild five GaN 45.7-GPa figure curves from hashed full ASE trajectories.

This validates the saved *final evaluated chain* against plotted enthalpies.
It does not parse the underlying calculator logs or recompute VCNEB forces.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read, write
from ase.units import GPa


BACKENDS = ("abacus", "vasp", "qe", "abinit", "cp2k")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_pairs(values: list[str], *, label: str) -> dict[str, str]:
    parsed = {}
    for value in values:
        backend, delimiter, item = value.partition("=")
        if not delimiter or backend not in BACKENDS or not item or backend in parsed:
            raise ValueError(f"invalid {label}: {value}")
        parsed[backend] = item
    if set(parsed) != set(BACKENDS):
        raise ValueError(f"{label} must specify exactly {BACKENDS}")
    return parsed


def audit(evidence_path: Path, source_data_path: Path,
          chains: dict[str, Path], expected_hashes: dict[str, str],
          snapshot_dir: Path, output_path: Path) -> dict:
    if output_path.exists():
        raise FileExistsError(output_path)
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    with source_data_path.open(newline="", encoding="utf-8") as handle:
        source_rows = list(csv.DictReader(handle))
    n_images = int(evidence["contract"]["n_images_total"])
    pressure_gpa = float(evidence["contract"]["pressure_gpa"])
    if n_images != 29 or pressure_gpa != 45.7:
        raise ValueError("not the reviewed 29-image GaN 45.7-GPa contract")
    if set(chains) != set(BACKENDS) or set(expected_hashes) != set(BACKENDS):
        raise ValueError("all five reviewed backends are required")

    records = {}
    snapshots = {}
    for backend in BACKENDS:
        source = chains[backend]
        source_hash = sha256(source)
        if source_hash != expected_hashes[backend]:
            raise ValueError(f"{backend} full-trajectory SHA-256 changed")
        frames = read(source, index=":")
        if not frames or len(frames) % n_images:
            raise ValueError(f"{backend} has an incomplete trajectory block")
        images = frames[-n_images:]
        atom_order = images[0].get_chemical_symbols()
        if len(atom_order) != 4 or sorted(atom_order) != ["Ga", "Ga", "N", "N"]:
            raise ValueError(f"{backend} has the wrong cell composition")
        energies = []
        volumes = []
        for image in images:
            energy = float(image.get_potential_energy())
            volume = float(image.get_volume())
            forces = np.asarray(image.get_forces(), dtype=float)
            stress = np.asarray(image.get_stress(), dtype=float)
            if (image.get_chemical_symbols() != atom_order or forces.shape != (4, 3)
                    or stress.shape != (6,) or not np.isfinite(energy)
                    or not np.isfinite(volume) or volume <= 0
                    or not np.isfinite(forces).all() or not np.isfinite(stress).all()):
                raise ValueError(f"{backend} has an invalid final-chain image")
            energies.append(energy)
            volumes.append(volume)
        enthalpy = np.asarray(energies) + pressure_gpa * GPa * np.asarray(volumes)
        relative = enthalpy - enthalpy[0]
        record = evidence["backends"][backend]
        rows = [row for row in source_rows if row["backend"] == backend.upper()]
        if len(rows) != n_images or [int(row["image_index"]) for row in rows] != list(range(n_images)):
            raise ValueError(f"{backend} figure source does not contain one ordered final chain")
        plotted = np.asarray([float(row["relative_enthalpy_eV_per_GaN"]) * 2
                              for row in rows])
        max_curve_error = float(np.max(np.abs(relative - plotted)))
        if max_curve_error > 3e-9:
            raise ValueError(f"{backend} trajectory differs from the plotted chain")
        barrier = float(np.max(relative[1:-1]))
        reaction = float(relative[-1])
        if (record["status"] != "converged" or int(np.argmax(relative)) != record["highest_image_index"]
                or not np.isclose(barrier, record["barrier_eV_per_cell"], atol=1e-9, rtol=0)
                or not np.isclose(reaction, record["reaction_enthalpy_eV_per_cell"], atol=1e-9, rtol=0)):
            raise ValueError(f"{backend} does not match the compact evidence claim")
        if "relative_enthalpy_eV_per_cell" in record:
            np.testing.assert_allclose(relative, record["relative_enthalpy_eV_per_cell"], atol=1e-9, rtol=0)
        snapshot_path = snapshot_dir / f"gan_{backend}_45p7_final_chain.traj"
        if snapshot_path.exists():
            raise FileExistsError(snapshot_path)
        snapshots[backend] = (snapshot_path, images)
        records[backend] = {
            "job_id": record["job_id"],
            "full_trajectory_sha256": source_hash,
            "n_full_trajectory_frames": len(frames),
            "n_complete_chain_blocks": len(frames) // n_images,
            "n_final_images": n_images,
            "n_formula_units_per_image": 2,
            "peak_image_index": int(np.argmax(relative)),
            "barrier_eV_per_GaN": barrier / 2,
            "reaction_enthalpy_eV_per_GaN": reaction / 2,
            "maximum_figure_curve_difference_eV_per_cell": max_curve_error,
            "all_final_image_energies_forces_stresses_finite": True,
            "limitations": "ASE trajectory/figure audit, not raw calculator-log or NEB-force certification",
        }
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    for backend, (path, images) in snapshots.items():
        write(path, images)
        records[backend]["compact_final_chain_sha256"] = sha256(path)
        records[backend]["compact_final_chain_file"] = path.name
    output = {
        "kind": "gan_45p7_five_backend_final_trajectory_curve_audit_not_raw_calculator_audit",
        "status": "five_hashed_final_trajectories_match_plotted_enthalpy_curves",
        "pressure_gpa": pressure_gpa,
        "n_images_total": n_images,
        "backends": records,
        "source_sha256": {
            "compact_evidence": sha256(evidence_path),
            "figure_source_data": sha256(source_data_path),
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--source-data", type=Path, required=True)
    parser.add_argument("--chain", action="append", required=True,
                        help="backend=path to complete original ASE trajectory")
    parser.add_argument("--expected-source-sha256", action="append", required=True,
                        help="backend=SHA-256 pinned by the remote source index")
    parser.add_argument("--snapshot-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    chains = {name: Path(path) for name, path in parse_pairs(args.chain, label="chain").items()}
    expected = parse_pairs(args.expected_source_sha256, label="expected source hash")
    output = audit(args.evidence, args.source_data, chains, expected,
                   args.snapshot_dir, args.output)
    print(json.dumps({"status": output["status"],
                      "maximum_curve_error_eV_per_cell": max(
                          row["maximum_figure_curve_difference_eV_per_cell"]
                          for row in output["backends"].values())}, indent=2))


if __name__ == "__main__":
    main()
