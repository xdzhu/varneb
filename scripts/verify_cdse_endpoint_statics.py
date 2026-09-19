"""Gate constructed phase endpoints using full-SCF raw forces and stress."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.filters import FrechetCellFilter
from ase.io import read


def raw_stationarity(atoms, summary, threshold=.02):
    if summary["status"] != "completed":
        raise ValueError("completed static SCF required")
    forces = np.asarray(summary["forces_eV_per_A"], dtype=float)
    stress = np.asarray(summary["stress_eV_per_A3_voigt"], dtype=float)
    energy = float(summary["potential_energy_eV"])
    if forces.shape != (len(atoms), 3) or stress.shape != (6,) or not all(
            np.isfinite(value).all() for value in (forces, stress, energy)):
        raise ValueError("finite complete static results required")
    image = atoms.copy()
    image.set_constraint()
    image.calc = SinglePointCalculator(image, energy=energy, forces=forces, stress=stress)
    fmax = float(np.linalg.norm(FrechetCellFilter(image).get_forces(), axis=1).max())
    if fmax > threshold:
        raise ValueError(f"phase candidate is not stationary: raw generalized force {fmax} > {threshold}")
    return {"raw_maximum_generalized_force_eV_per_A": fmax, "threshold_eV_per_A": threshold,
            "symmetry_force_projection": False, "additional_DFT_evaluations": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-root", type=Path, required=True)
    args = parser.parse_args()
    output = args.case_root / "endpoint_static_gate.json"
    records = {}
    for label in ("initial", "final"):
        path = args.case_root / f"static_{label}" / "vasp_static_summary.json"
        atoms = read(args.case_root / "input" / label / "CONTCAR")
        records[label] = raw_stationarity(atoms, json.loads(path.read_text()))
        records[label]["static_summary_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    with output.open("x") as handle:
        json.dump({"all_passed": True, "records": records}, handle, indent=2)
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
