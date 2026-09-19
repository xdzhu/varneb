"""Create a separate diagnostic fractional optimizer step, never edit a run."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read, write
from ase.io.trajectory import Trajectory

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from vcneb.vasp_contract import validate_vasp_image_geometry
from vcneb.vasp_lattice import NativeVaspLatticeProbe


def fractional_step(previous, proposed, fraction):
    if not np.isfinite(fraction) or not 0 < fraction < 1:
        raise ValueError("fraction must be strictly between zero and one")
    if previous.get_chemical_symbols() != proposed.get_chemical_symbols():
        raise ValueError("atom order/species must be identical")
    for atoms in (previous, proposed):
        validate_vasp_image_geometry(atoms, minimum_distance=1.4)
    result = previous.copy()
    result.set_constraint()
    # VCNEB's deformation and unwrapped fractional coordinates are both linear
    # in its generalized coordinates for a fixed reference cell/cell_scale.
    result.set_cell((1-fraction)*previous.cell.array+fraction*proposed.cell.array)
    result.set_scaled_positions((1-fraction)*previous.get_scaled_positions(wrap=False)
                                +fraction*proposed.get_scaled_positions(wrap=False))
    validate_vasp_image_geometry(result, minimum_distance=1.4)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluated-trajectory", required=True, type=Path)
    parser.add_argument("--n-images", type=int, default=29)
    parser.add_argument("--proposed-structure", required=True, type=Path)
    parser.add_argument("--image-index", type=int, default=10)
    parser.add_argument("--fraction", type=float, default=0.5)
    parser.add_argument("--native-checker", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError("refusing to overwrite canary artifacts")
    if args.n_images < 3 or not 0 < args.image_index < args.n_images-1:
        raise ValueError("interior image index and valid total image count required")
    with Trajectory(args.evaluated_trajectory) as trajectory:
        if not len(trajectory) or len(trajectory)%args.n_images:
            raise ValueError("evaluated trajectory must contain only complete chains")
        previous=trajectory[len(trajectory)-args.n_images+args.image_index]
    if previous.calc is None or not {"energy","forces","stress"}.issubset(previous.calc.results):
        raise ValueError("reference must have a complete static result")
    proposed=read(args.proposed_structure)
    candidate=fractional_step(previous,proposed,args.fraction)
    probe=NativeVaspLatticeProbe(args.native_checker)
    classifications={name:probe.classify(atoms.cell.array).to_dict() for name,atoms in
                     [("previous",previous),("proposed",proposed),("fractional_step",candidate)]}
    if not classifications["previous"]["consistent"] or not classifications["fractional_step"]["consistent"]:
        raise ValueError("reference and canary must pass native classification")
    args.output_dir.mkdir(parents=True)
    artifacts={}
    for name,atoms in [("reference",previous),("fractional_step",candidate)]:
        path=args.output_dir/(name+".vasp")
        write(path,atoms,format="vasp",direct=True,vasp5=True)
        parsed=read(path)
        if (parsed.get_chemical_symbols()!=atoms.get_chemical_symbols()
            or not np.allclose(parsed.positions,atoms.positions,atol=1e-10,rtol=0)
            or not np.allclose(parsed.cell.array,atoms.cell.array,atol=1e-10,rtol=0)):
            raise ValueError("canary POSCAR roundtrip mismatch")
        artifacts[name]={"sha256":sha256(path.read_bytes()).hexdigest(),"classification":
                         classifications["previous" if name=="reference" else "fractional_step"]}
    manifest={"mode":"diagnostic_only_no_production_resume","fraction":args.fraction,
              "image_index_zero_based":args.image_index,"probe":probe.descriptor(),
              "parent_sha256":{"evaluated_trajectory":sha256(args.evaluated_trajectory.read_bytes()).hexdigest(),
                               "proposed_structure":sha256(args.proposed_structure.read_bytes()).hexdigest()},
              "reference_results":{"energy_eV":float(previous.calc.results["energy"]),
                                   "forces_eV_per_A":np.asarray(previous.calc.results["forces"]).tolist(),
                                   "stress_eV_per_A3":np.asarray(previous.calc.results["stress"]).tolist()},
              "classifications":classifications,"artifacts":artifacts,
              "coordinate_policy":"linear unwrapped fractional coordinates and cell; no MIC or projection",
              "units":{"cell":"angstrom","energy":"eV","forces":"eV/angstrom","stress":"eV/angstrom^3"}}
    (args.output_dir/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps(classifications,indent=2))


if __name__=="__main__":
    main()
