"""Track a rotated T-distortion triplet without assuming phonon/irrep labels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ase.io import read
import numpy as np

from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_reference_variants import distortion, fluorite_scaffold
from vcneb.reference_variants import apply_parent_operation


def rotated_t_triplet(t):
    parent = fluorite_scaffold(t)
    # Unit cubic scaffold is used only to obtain exact cubic site permutations.
    # No structure, cell, or energy of a DFT calculation is changed by this step.
    cubic = parent.copy()
    cubic.set_cell(np.eye(3), scale_atoms=True)
    delta = t.get_scaled_positions() - parent.get_scaled_positions()
    delta -= np.rint(delta)
    rotations = (np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]),
                 np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]), np.eye(3))
    basis, mappings = [], []
    for rotation in rotations:
        permutation = apply_parent_operation(cubic, cubic, rotation, [0, 0, 0]).source_to_target
        rotated = np.empty_like(delta)
        rotated[np.array(permutation)] = delta @ rotation.T
        vector = rotated @ parent.cell.array
        vector -= vector.mean(axis=0)
        basis.append(vector / np.linalg.norm(vector))
        mappings.append(list(permutation))
    basis = np.asarray(basis)
    if not np.allclose(basis.reshape(3, 36) @ basis.reshape(3, 36).T, np.eye(3), atol=1e-10):
        raise ValueError("rotated T patterns fail orthonormality")
    return parent, basis, mappings


def analyze(variants, ordinary, guided, output):
    if output.exists():
        raise FileExistsError("refusing existing pattern analysis")
    t = read(variants / "T.vasp", format="vasp")
    parent, basis, mappings = rotated_t_triplet(t)
    sources = {"T": variants / "T.vasp", "PO": variants / "PO.vasp",
               "PO_minus_T_preserving_inversion": variants / "PO_minus_T_preserving_inversion.vasp",
               "PO_minus_T_reversing_inversion": variants / "PO_minus_T_reversing_inversion.vasp"}
    endpoints = {}
    for name, path in sources.items():
        vector = distortion(read(path, format="vasp"), parent)
        amplitudes = np.einsum("aij,ij->a", basis, vector)
        endpoints[name] = {"Q_A": amplitudes.tolist(), "source_sha256": sha256(path)}
    chains = {}
    for name, directory in (("ordinary", ordinary), ("guided", guided)):
        records = []
        for index in range(7):
            path = directory / f"POSCAR_{index:02d}"
            atoms = read(path, format="vasp")
            # This audit only accepts the local quarter-site chart, far from a
            # nearest-reference wrap discontinuity. A longer path needs an
            # explicitly continuous periodic gauge instead.
            fractional = atoms.get_scaled_positions() - parent.get_scaled_positions()
            fractional -= np.rint(fractional)
            if np.max(np.abs(fractional)) > .45:
                raise ValueError("path leaves the declared local parent-site chart")
            vector = distortion(atoms, parent)
            amplitudes = np.einsum("aij,ij->a", basis, vector)
            residual = vector - np.einsum("a,aij->ij", amplitudes, basis)
            records.append({"index": index, "Q_A": amplitudes.tolist(),
                            "orthogonal_residual_norm_A": float(np.linalg.norm(residual)),
                            "volume_A3": atoms.get_volume(), "source_sha256": sha256(path)})
        chains[name] = records
    report = {"schema_version": 1, "basis": basis.tolist(), "rotated_parent_permutations": mappings,
              "basis_order": ["rotated_T_pattern_x", "rotated_T_pattern_y", "T_pattern_z"],
              "reference_metric": "original T lattice, Cartesian Angstrom norm",
              "endpoints": endpoints, "chains": chains,
              "limitations": ["geometric pattern amplitudes; not energy contributions or phonon occupations",
                              "rotated T-pattern subspace; no unverified X2- irrep label",
                              "two historical paths share ordered endpoints; this is not a new DFT chain",
                              "orthogonal residual includes other atomic distortions; cell strain is recorded separately"]}
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("variants", "ordinary", "guided", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    report = analyze(args.variants, args.ordinary, args.guided, args.output)
    print(json.dumps({"endpoints": report["endpoints"], "chains_Q_A": {k: [x["Q_A"] for x in v] for k, v in report["chains"].items()}}))


if __name__ == "__main__":
    main()
