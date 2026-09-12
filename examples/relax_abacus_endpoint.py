"""Relax one ABACUS endpoint with ASE's cell filter before VC-NEB."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import BFGS, FIRE, LBFGS

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.abacus import attach_abacus_calculators, make_ase_abacus_factory
from vcneb.calculator import validate_image_calculators


DEFAULT_PSEUDO_DIR = "/home/zhuxd/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential"
DEFAULT_BASIS_DIR = "/home/zhuxd/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/selected_Orbs"
DEFAULT_COMMAND = os.environ.get("ABACUS_COMMAND", "abacus")


def parse_species_files(values: list[str], option: str) -> dict[str, str]:
    result = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"{option} expects SPECIES=FILE, got {value!r}")
        species, filename = value.split("=", 1)
        if not species or not filename:
            raise ValueError(f"{option} expects non-empty SPECIES and FILE, got {value!r}")
        result[species] = filename
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--command", default=DEFAULT_COMMAND)
    parser.add_argument("--pseudo-dir", default=os.environ.get("ABACUS_PP_PATH"))
    parser.add_argument("--basis-dir", default=os.environ.get("ABACUS_ORBITAL_PATH"))
    parser.add_argument("--pp", action="append", default=[], metavar="SPECIES=FILE")
    parser.add_argument("--basis", action="append", default=[], metavar="SPECIES=FILE")
    parser.add_argument("--ecutwfc", type=float, default=60.0)
    parser.add_argument("--scf-thr", type=float, default=1e-6)
    parser.add_argument("--scf-nmax", type=int, default=100)
    parser.add_argument("--mixing-beta", type=float, default=0.7)
    parser.add_argument("--kpts", type=int, nargs=3, default=[1, 1, 1], metavar=("NX", "NY", "NZ"))
    parser.add_argument("--fmax", type=float, default=0.05)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--optimizer", choices=["FIRE", "BFGS", "LBFGS"], default="FIRE")
    parser.add_argument(
        "--maxstep",
        type=float,
        default=None,
        help="Optional ASE optimizer maximum step in Angstrom-like coordinates",
    )
    parser.add_argument("--fixed-cell", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    structure = Path(args.structure).resolve()
    workdir = Path(args.workdir).resolve()
    atoms = read(structure)
    parameters = {
        "calculation": "scf",
        "basis_type": "lcao",
        "dft_functional": "pbe",
        "ecutwfc": args.ecutwfc,
        "scf_thr": args.scf_thr,
        "scf_nmax": args.scf_nmax,
        "mixing_type": "pulay",
        "mixing_beta": args.mixing_beta,
        "kpts": args.kpts,
        "cal_force": 1,
        "cal_stress": 1,
        "out_stru": 1,
    }
    if args.pseudo_dir:
        parameters["pseudo_dir"] = args.pseudo_dir
    if args.basis_dir:
        parameters["basis_dir"] = args.basis_dir
    pp = parse_species_files(args.pp, "--pp")
    basis = parse_species_files(args.basis, "--basis")
    if pp:
        parameters["pp"] = pp
    if basis:
        parameters["basis"] = basis

    factory = make_ase_abacus_factory(parameters=parameters, command=args.command)
    attach_abacus_calculators([atoms], workdir=workdir, factory=factory)
    validate_image_calculators(
        [atoms],
        require_stress=True,
        require_variable_cell=True,
        require_directory=True,
        require_unique_directories=True,
    )

    target = atoms if args.fixed_cell else FrechetCellFilter(atoms)
    optimizer_class = {"FIRE": FIRE, "BFGS": BFGS, "LBFGS": LBFGS}[args.optimizer]
    optimizer_kwargs = {
        "logfile": str(workdir / "relax.log"),
        "trajectory": str(workdir / "relax.traj"),
    }
    if args.maxstep is not None:
        optimizer_kwargs["maxstep"] = args.maxstep
    optimizer = optimizer_class(target, **optimizer_kwargs)
    optimizer.run(fmax=args.fmax, steps=args.steps)
    energy = float(atoms.get_potential_energy())
    forces = np.asarray(atoms.get_forces(), dtype=float)
    stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
    generalized_forces = np.asarray(target.get_forces(), dtype=float)
    write(workdir / "CONTCAR", atoms, format="vasp", direct=True, vasp5=True)
    write(workdir / "POSCAR.final", atoms, format="vasp", direct=True, vasp5=True)
    max_force = float(np.linalg.norm(forces, axis=1).max())
    max_generalized_force = float(np.linalg.norm(generalized_forces, axis=1).max())
    summary = {
        "structure": str(structure),
        "workdir": str(workdir),
        "natoms": len(atoms),
        "fixed_cell": args.fixed_cell,
        "optimizer": args.optimizer,
        "steps_requested": args.steps,
        "optimizer_steps": int(getattr(optimizer, "nsteps", -1)),
        "maxstep": args.maxstep,
        "fmax_target_eV_per_A": args.fmax,
        "energy_eV": energy,
        "max_force_eV_per_A": max_force,
        "max_generalized_force_eV_per_A": max_generalized_force,
        "max_abs_stress_eV_per_A3": float(np.abs(stress).max()),
        "stress_eV_per_A3": stress.tolist(),
        "volume_A3": float(atoms.get_volume()),
        "converged": bool(max_generalized_force < args.fmax),
    }
    (workdir / "relax_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
