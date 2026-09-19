"""Prepare exact unimodular basal-basis diagnostic inputs; never symmetrize.

Only a complete, equal-sized basal Gamma mesh is supported here. This is a
diagnostic generator, not an automatic calculator recovery or production policy.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np
from ase.io import read, write


TRANSFORMS = {
    "identity": np.eye(3, dtype=int),
    "basal_sum_b": np.array([[1, 0, 0], [1, 1, 0], [0, 0, 1]]),
    "basal_sum_a": np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]]),
}


def mesh_equivalence(transform, mesh):
    mesh = np.asarray(mesh, dtype=int)
    if mesh.shape != (3,) or np.any(mesh <= 0):
        raise ValueError("positive three-dimensional mesh required")
    if mesh[0] != mesh[1]:
        raise ValueError("this basal probe requires equal basal mesh counts")
    grid = np.array(list(product(*(range(int(value)) for value in mesh))), dtype=float) / mesh
    # C'=T C, B'=T^-T B; express new k coordinates in the old reciprocal basis.
    original_coordinates = grid @ np.linalg.inv(transform).T
    integer_coordinates = original_coordinates * mesh
    if not np.allclose(integer_coordinates, np.rint(integer_coordinates), atol=1e-10, rtol=0):
        raise ValueError("basis transformation changes the physical Gamma mesh")
    mapped = np.rint(integer_coordinates).astype(int) % mesh
    if len({tuple(row) for row in mapped}) != len(grid):
        raise ValueError("reciprocal mesh mapping is not bijective")
    return len(grid)


def equivalent_atoms(atoms, transform, mesh):
    if not np.all(atoms.pbc) or atoms.get_volume() <= 0 or not np.all(np.isfinite(atoms.positions)):
        raise ValueError("finite, positive-volume, three-dimensional periodic input required")
    transform = np.asarray(transform)
    if transform.shape != (3, 3) or not np.array_equal(transform, np.rint(transform)) or not np.isclose(np.linalg.det(transform), 1):
        raise ValueError("orientation-preserving unimodular integer transformation required")
    mesh_equivalence(transform, mesh)
    result = atoms.copy()
    result.set_cell(transform @ atoms.cell.array, scale_atoms=False)
    if not np.array_equal(result.positions, atoms.positions):
        raise ValueError("Cartesian atom positions changed")
    if not np.isclose(result.get_volume(), atoms.get_volume(), rtol=0, atol=1e-10):
        raise ValueError("periodic cell volume changed")
    if not np.allclose(result.get_all_distances(mic=True), atoms.get_all_distances(mic=True), atol=1e-10, rtol=0):
        raise ValueError("periodic pair distances changed")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--mesh", nargs=3, type=int, default=[8, 8, 6])
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError("refusing to overwrite diagnostic inputs")
    atoms = read(args.structure)
    args.output_dir.mkdir(parents=True)
    payload = {"parent_sha256": sha256(args.structure.read_bytes()).hexdigest(),
               "input": str(args.structure.resolve()), "mesh": args.mesh,
               "units": {"lattice": "angstrom", "angles": "degree"},
               "coordinate_policy": "Cartesian unchanged; only lattice basis changes",
               "no_symmetry_projection_or_geometry_perturbation": True, "candidates": {}}
    for label, transform in TRANSFORMS.items():
        candidate = equivalent_atoms(atoms, transform, args.mesh)
        path = args.output_dir / f"{label}.vasp"
        write(path, candidate, format="vasp", direct=True, vasp5=True)
        reread = read(path)
        if not np.allclose(reread.positions, atoms.positions, atol=1e-10, rtol=0):
            raise ValueError("POSCAR roundtrip changed atom positions")
        payload["candidates"][label] = {"transform": transform.tolist(),
            "sha256": sha256(path.read_bytes()).hexdigest(), "lattice_A": candidate.cell.array.tolist(),
            "angles_deg": candidate.cell.angles().tolist(), "volume_A3": candidate.get_volume(),
            "physical_k_points_verified": mesh_equivalence(transform, args.mesh)}
    (args.output_dir / "basis_manifest.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
