"""Prepare eight fixed-contract local probes, not a TS or phonon calculation."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase import Atoms
from ase.io import read
from ase.stress import voigt_6_to_full_3x3_stress

from scripts.audit_hfo2_static_replica import INPUT_FILES, audited_results, sha256
from vcneb.joint_curvature import JointCurvatureCoordinates


def write_structure(path: Path, atoms: Atoms) -> None:
    """Write only geometry; retain the declared historical species/basis order."""
    if atoms.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8:
        raise ValueError("expected ordered Hf4O8")
    if not np.isfinite(atoms.positions).all() or atoms.get_volume() <= 0:
        raise ValueError("invalid probe geometry")
    lines = ["ATOMIC_SPECIES", "Hf 178.49 Hf.upf", "O 15.999 O.upf", "",
             "NUMERICAL_ORBITAL", "Hf_gga_10au_100Ry_4s2p2d1f.orb",
             "O_gga_10au_100Ry_2s2p1d.orb", "", "LATTICE_CONSTANT",
             "1.8897261258369282", "", "LATTICE_VECTORS"]
    lines.extend(" ".join(f"{value:.17g}" for value in row) for row in atoms.cell.array)
    lines.extend(["", "ATOMIC_POSITIONS", "Direct", ""])
    scaled = atoms.get_scaled_positions(wrap=True)
    for label, begin, end in (("Hf", 0, 4), ("O", 4, 12)):
        lines.extend([label, "0.0", str(end - begin)])
        lines.extend(" ".join(f"{value:.17g}" for value in row) + " 1 1 1" for row in scaled[begin:end])
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def probe_directions(center: Atoms, previous: Atoms, following: Atoms) -> np.ndarray:
    if any(image.get_chemical_symbols() != center.get_chemical_symbols()
           for image in (previous, following)):
        raise ValueError("neighbor atom identity/order differs")
    delta = following.get_scaled_positions() - previous.get_scaled_positions()
    delta -= np.rint(delta)
    atomic = delta @ center.cell.array
    atomic -= atomic.mean(axis=0)
    length = np.linalg.norm(atomic)
    if length < 1e-8:
        raise ValueError("atomic path secant is degenerate")
    directions = np.zeros((2, 3 * len(center) + 6))
    directions[0, :3 * len(center)] = atomic.ravel() / length
    directions[1, 3 * len(center):3 * len(center) + 2] = 1 / np.sqrt(2)
    return directions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--chain", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("refusing existing probe namespace")
    if sha256(args.source / "INPUT") != "dc6684ffa709bf3d953c647272589b05d844c1293477bd144226d495207cb9df":
        raise ValueError("historical 100-Ry INPUT contract changed")
    if sha256(args.source / "KPT") != "92b917107d9df11da28a465cc4900f3574956a057ddf9f7ad125e75c956ec508":
        raise ValueError("historical Gamma 2x2x2 KPT contract changed")
    center = read(args.source / "STRU", format="abacus")
    previous = read(args.chain / "POSCAR_01", format="vasp")
    following = read(args.chain / "POSCAR_03", format="vasp")
    chart = JointCurvatureCoordinates(center, center.get_volume() ** (1 / 3))
    directions = probe_directions(center, previous, following)
    raw = audited_results(args.source)
    center_gradient = chart.enthalpy_gradient(
        center, raw["forces"], voigt_6_to_full_3x3_stress(raw["stress"]), 0.0,
    )
    args.output.mkdir(parents=True)
    points = []
    for direction_id, direction in enumerate(directions):
        for step in (0.01, 0.02):
            for sign in (-1, 1):
                index = len(points)
                geometry = chart.displaced(sign * step * direction)
                distances = geometry.get_all_distances(mic=True)
                np.fill_diagonal(distances, np.inf)
                if distances.min() < 1.6:
                    raise ValueError("probe violates minimum-distance guard")
                directory = args.output / "points" / f"{index:02d}"
                directory.mkdir(parents=True)
                for name in INPUT_FILES:
                    if name != "STRU":
                        shutil.copyfile(args.source / name, directory / name)
                write_structure(directory / "STRU", geometry)
                # Roundtrip through the actual HF ASE-ABACUS reader before DFT.
                restored = read(directory / "STRU", format="abacus")
                if not np.allclose(restored.cell.array, geometry.cell.array, atol=1e-11, rtol=0):
                    raise ValueError("STRU lattice unit/precision roundtrip failed")
                delta = restored.get_scaled_positions() - geometry.get_scaled_positions()
                delta -= np.rint(delta)
                if not np.allclose(delta @ geometry.cell.array, 0, atol=1e-11, rtol=0):
                    raise ValueError("STRU position/unit roundtrip failed")
                points.append({"index": index, "direction": direction_id, "step_A": step,
                               "sign": sign, "minimum_distance_A": float(distances.min()),
                               "input_sha256": {name: sha256(directory / name) for name in INPUT_FILES}})
    manifest = {
        "schema_version": 1, "source_directory": str(args.source),
        "source_log_sha256": sha256(args.source / "OUT.ABACUS/running_scf.log"),
        "source_STRU_sha256": sha256(args.source / "STRU"),
        "pressure_eV_A3": 0, "cell_scale_A": chart.cell_scale_A,
        "center_energy_eV_cell": float(raw["energy"]), "center_gradient_eV_A": center_gradient.tolist(),
        "directions": directions.tolist(),
        "direction_semantics": ["atomic component of local chain secant, fixed cell; not a phonon",
                                "scaled symmetric xx+yy strain with fixed fractional atoms"],
        "predicted_center_directional_gradients_eV_A": (directions @ center_gradient).tolist(),
        "formula_units": 4, "n_new_static_points": len(points), "points": points,
        "interpretation": "nonstationary local work/curvature diagnostic; never a TS certificate",
    }
    (args.output / "probe_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"n_points": len(points), "directional_gradients_eV_A": (directions @ center_gradient).tolist()}))


if __name__ == "__main__":
    main()
