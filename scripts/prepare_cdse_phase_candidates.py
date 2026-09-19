"""Construct explicitly phase-constrained CdSe endpoint *candidates*.

The completed unconstrained-cell endpoints remain untouched.  These candidates
must pass a new full SCF and a raw force/virial gate before production.  This
is endpoint construction, never a transformation of an optimized NEB image.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
from ase import Atoms
from ase.io import read, write

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.prepare_cdse_sheppard_suite import symmetry


def phase_candidate(atoms, label):
    if atoms.get_chemical_symbols() != ["Cd"] * 4 + ["Se"] * 4:
        raise ValueError("ordered Cd4Se4 endpoints required")
    h = atoms.cell.array
    q = atoms.get_scaled_positions(wrap=False)
    if label == "rs":
        template = np.array([[2., 0, 0], [.5, .5, 0], [0, 0, 1.]])
        cell = template * np.sum(h * template) / np.sum(template**2)
        motif = np.array([[0,0,0],[.5,0,0],[.25,0,.5],[.75,0,.5],
                          [.25,0,0],[.75,0,0],[0,0,.5],[.5,0,.5]])
        expected = 225
    elif label == "wz":
        template = np.array([[2.,0,0],[.5,np.sqrt(3)/2,0],[0,0,0]])
        cell = template * np.sum(h * template) / np.sum(template**2)
        cell[2,2] = h[2,2]
        offsets = np.array([0,0,.5,.5])
        u = np.mean(q[4:,2] - offsets) - np.mean(q[:4,2] - offsets)
        if not .05 < u < .25:
            raise ValueError("WZ internal coordinate outside the declared motif")
        motif = np.array([[0,0,0],[.5,0,0],[1/6,1/3,.5],[2/3,1/3,.5],
                          [1/6,1/3,u],[2/3,1/3,u],[0,0,.5+u],[.5,0,.5+u]])
        expected = 186
    else:
        raise ValueError("unknown phase")
    # Retain the common translation, not a fitted species permutation.
    difference = q - motif
    difference -= np.rint(difference)
    motif += difference.mean(axis=0)
    candidate = Atoms(atoms.symbols, cell=cell, scaled_positions=motif, pbc=True)
    strain = float(np.max(np.abs(np.linalg.svd(np.linalg.solve(cell, h), compute_uv=False) - 1)))
    displacement = candidate.get_scaled_positions(wrap=False) - q
    displacement -= np.rint(displacement)
    displacement = float(np.linalg.norm(displacement @ cell, axis=1).max())
    # Fixed before use: this is a bounded reconstruction of numerical drift,
    # not permission to turn an arbitrary relaxed structure into the target.
    if strain > 1e-3 or displacement > 1e-3:
        raise ValueError("endpoint too far from the declared phase for bounded reconstruction")
    if any(record["number"] != expected for record in symmetry(candidate)):
        raise ValueError("candidate does not have the declared exact phase")
    return candidate, {"phase": label, "expected_spacegroup": expected,
                       "maximum_principal_stretch_change": strain,
                       "maximum_motif_displacement_A": displacement,
                       "original_symmetry": symmetry(atoms), "candidate_symmetry": symmetry(candidate),
                       "candidate_only_not_stationarity_verified": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-suite", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    records = {}
    candidates = {}
    for label in ("rs", "wz"):
        origin = args.original_suite / "endpoint_relax" / label
        summary = json.loads((origin / "endpoint_relax_summary.json").read_text())
        if not summary["converged"] or summary["external_pressure_gpa"] != 0:
            raise ValueError("completed zero-pressure relaxation required")
        candidates[label], records[label] = phase_candidate(read(origin / "CONTCAR"), label)
        records[label]["parent_sha256"] = {name: hashlib.sha256((origin/name).read_bytes()).hexdigest()
                                            for name in ("CONTCAR", "endpoint_relax_summary.json")}
    args.output.mkdir(parents=True, exist_ok=False)
    for label in ("rs", "wz"):
        directory = args.output / label
        directory.mkdir()
        write(directory / "CONTCAR", candidates[label], format="vasp", direct=True, vasp5=True)
        for name in ("INCAR", "KPOINTS", "POTCAR"):
            shutil.copy2(args.original_suite / "endpoint_seeds" / label / name, directory / name)
        (directory / "phase_candidate_audit.json").write_text(json.dumps(records[label], indent=2)+"\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
