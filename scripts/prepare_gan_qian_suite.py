"""Build audited GaN initial paths for Qian et al. (2013), without DFT.

Hexagonal intermediates seed a branch, not a constraint on the converged band.
The B3 route deliberately retains the conventional-cell diagonal displacement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
from ase import Atoms
from ase.io import read, write
import spglib

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from vcneb import endpoint_structure_record, validate_path_geometry

FCC = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])


def from_record(record):
    return Atoms(record["species_order"], cell=record["cell_A"],
                 scaled_positions=record["fractional_positions_wrapped"], pbc=True)


def conventional_gan(a, offset):
    return Atoms("Ga4N4", cell=np.eye(3) * a,
                 scaled_positions=np.vstack((FCC, FCC + offset)), pbc=True)


def symmetry_audit(atoms):
    cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)
    records = []
    for tolerance in (1e-5, 1e-4, 1e-3):
        result = spglib.get_symmetry_dataset(cell, symprec=tolerance)
        records.append({"symprec_A": tolerance, "number": int(result.number),
                        "international": result.international})
    return records


def hexagonal_chain(initial, final, n_images=29, anchor_index=11):
    if len(initial) != 4 or initial.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]:
        raise ValueError("hexagonal route requires the audited ordered 2-formula-unit cell")
    q0, q1 = (a.get_scaled_positions(wrap=False) for a in (initial, final))
    q1 += np.rint(q0 - q1)
    # The signed N-Ga separation is a half period at the h-MgO anchor.
    anchor_q = q0.copy()
    for ga, nitrogen in ((0, 2), (1, 3)):
        offset = q0[nitrogen, 2] - q0[ga, 2]
        anchor_q[nitrogen, 2] = q0[ga, 2] + (.5 if offset > 0 else -.5)
    fraction = anchor_index / (n_images - 1)
    lengths = (1 - fraction) * initial.cell.lengths() + fraction * final.cell.lengths()
    anchor_cell = initial.cell.array * (lengths / initial.cell.lengths())[:, None]
    images = []
    for index in range(n_images):
        if index <= anchor_index:
            t = index / anchor_index
            cell = (1 - t) * initial.cell.array + t * anchor_cell
            q = (1 - t) * q0 + t * anchor_q
        else:
            t = (index - anchor_index) / (n_images - 1 - anchor_index)
            cell = (1 - t) * anchor_cell + t * final.cell.array
            q = (1 - t) * anchor_q + t * q1
        images.append(Atoms(initial.symbols, cell=cell, scaled_positions=q, pbc=True))
    return images


def body_diagonal_chain(initial, final, n_images=29):
    if initial.get_chemical_symbols() != ["Ga"] * 4 + ["N"] * 4:
        raise ValueError("B3/B1 route requires ordered conventional cells (4 formula units)")
    q0 = initial.get_scaled_positions(wrap=False)
    q1 = final.get_scaled_positions(wrap=False)
    q0 += np.rint(np.vstack((FCC, FCC + .25)) - q0)
    q1 += np.rint(np.vstack((FCC, FCC + .5)) - q1)
    if not np.allclose(q1 - q0, np.vstack((np.zeros((4, 3)), np.full((4, 3), .25))), atol=1e-6):
        raise ValueError("relaxed endpoints do not preserve the audited diagonal atom mapping")
    return [Atoms(initial.symbols, cell=(1 - t) * initial.cell.array + t * final.cell.array,
                  scaled_positions=(1 - t) * q0 + t * q1, pbc=True)
            for t in np.linspace(0, 1, n_images)]


def copy_contract(source, target):
    target.mkdir(parents=True, exist_ok=True)
    for name in ("INCAR", "KPOINTS", "POTCAR"):
        shutil.copy2(source / name, target / name)


def save_case(root, source, images, metadata):
    if (root / "initial.traj").exists():
        raise FileExistsError(f"refusing to replace a prepared path: {root}")
    geometry = validate_path_geometry(images, minimum_distance=1.4, maximum_deformation=.50)
    for label, atoms in (("initial", images[0]), ("final", images[-1])):
        directory = root / "input" / label
        copy_contract(source, directory)
        write(directory / "CONTCAR", atoms, format="vasp", direct=True, vasp5=True)
        reread = read(directory / "CONTCAR")
        if not np.allclose(reread.positions, atoms.positions, atol=1e-10):
            raise ValueError("POSCAR roundtrip changed atom mapping")
    write(root / "initial.traj", images)
    metadata.update({"n_images_total": len(images), "n_interiors": len(images) - 2,
                     "geometry": geometry, "coordinate_policy": "explicit_unwrapped_no_remapping",
                     "endpoint_structures": {label: endpoint_structure_record(atoms)
                         for label, atoms in (("initial", images[0]), ("final", images[-1]))},
                     "endpoint_symmetry": {label: symmetry_audit(atoms)
                         for label, atoms in (("initial", images[0]), ("final", images[-1]))},
                     "trajectory_sha256": hashlib.sha256((root / "initial.traj").read_bytes()).hexdigest()})
    (root / "preparation.json").write_text(json.dumps(metadata, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=("seeds", "hexagonal", "b3"))
    parser.add_argument("--parent-case", required=True, type=Path)
    parser.add_argument("--suite-root", required=True, type=Path)
    args = parser.parse_args()
    parent, suite = args.parent_case, args.suite_root
    baseline = json.loads((parent / "vcneb_hf_n29_w9_mpi32" / "vcneb_summary.json").read_text())
    if not baseline.get("converged") or baseline["calculator_parameters"]["symprec"] != 1e-4:
        raise ValueError("parent calculation must be converged with the agreed contract")
    initial, final = (from_record(baseline["endpoint_structures"][label]) for label in ("initial", "final"))
    source = parent / "vcneb_input" / "initial"
    if args.stage == "seeds":
        a_b3 = (2 * initial.get_volume()) ** (1 / 3)
        a_b1 = (2 * final.get_volume()) ** (1 / 3)
        for label, atoms, expected in (("b3", conventional_gan(a_b3, .25), 216),
                                       ("b1", conventional_gan(a_b1, .5), 225)):
            target = suite / "endpoint_seeds" / label
            if target.exists():
                raise FileExistsError(target)
            copy_contract(source, target)
            write(target / "POSCAR", atoms, format="vasp", direct=True, vasp5=True)
            symmetry = symmetry_audit(atoms)
            if any(record["number"] != expected for record in symmetry):
                raise ValueError(f"unexpected {label} seed symmetry")
            record = {"structure": endpoint_structure_record(atoms), "symmetry": symmetry,
                      "pressure_gpa_for_relaxation": 45.0, "seed_only_not_relaxed": True}
            (target / "seed.json").write_text(json.dumps(record, indent=2) + "\n")
    elif args.stage == "hexagonal":
        images = hexagonal_chain(initial, final)
        save_case(suite / "b4_b1_hexagonal", source, images,
                  {"route": "B4_to_B1_hexagonal", "pressure_gpa": 45.7,
                   "literature_barrier_eV_per_GaN": .39, "n_formula_units": 2,
                   "endpoints_reused_from_converged_tetragonal_case": str(parent),
                   "h_MgO_anchor_index_zero_based": 11,
                   "anchor_symmetry": symmetry_audit(images[11]),
                   "optimization_constraint": "none; branch identity must be checked after relaxation"})
    else:
        endpoints = []
        for label, expected in (("b3", 216), ("b1", 225)):
            directory = suite / "endpoint_relax" / label
            result = json.loads((directory / "endpoint_relax_summary.json").read_text())
            if not result["converged"] or result["external_pressure_gpa"] != 45.0:
                raise ValueError(f"{label} endpoint not converged at 45 GPa")
            atoms = read(directory / "CONTCAR")
            atoms.set_constraint()
            if symmetry_audit(atoms)[1]["number"] != expected:
                raise ValueError(f"{label} endpoint changed phase")
            endpoints.append(atoms)
        save_case(suite / "b3_b1_diagonal", source, body_diagonal_chain(*endpoints),
                  {"route": "B3_to_B1_conventional_body_diagonal", "pressure_gpa": 45.0,
                   "literature_barrier_eV_per_GaN": .57, "n_formula_units": 4,
                   "literature_expected_sequence": "B3-B1-B3-B1; three barriers",
                   "optimization_constraint": "none; intermediates and peak multiplicity must be audited"})
    print(f"[OK] prepared {args.stage} in {suite}")


if __name__ == "__main__":
    main()
