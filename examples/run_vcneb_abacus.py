"""Run VC-NEB with ABACUS through a calculator factory.

The default path uses ``ase.calculators.abacus.Abacus`` when available.  If
your production environment uses a different ABACUS ASE adapter, replace
``make_ase_abacus_factory(...)`` with a custom factory that returns one
calculator per image directory.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import interpolate_vcneb, read_chain_trajectory, run_vcneb
from vcneb.abacus import attach_abacus_calculators, make_ase_abacus_factory


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--initial", required=True, help="ASE-readable initial endpoint")
    parser.add_argument("--final", required=True, help="ASE-readable final endpoint")
    parser.add_argument("--workdir", default="abacus_vcneb_run")
    parser.add_argument("--n-images", type=int, default=7)
    parser.add_argument("--fmax", type=float, default=0.05)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--k", type=float, default=0.10)
    parser.add_argument("--pressure-gpa", type=float, default=0.0)
    parser.add_argument("--command", default=None)
    parser.add_argument("--pseudo-dir", default=os.environ.get("ABACUS_PP_PATH"))
    parser.add_argument("--basis-dir", default=os.environ.get("ABACUS_ORBITAL_PATH"))
    parser.add_argument("--pp", action="append", default=[], metavar="SPECIES=FILE")
    parser.add_argument("--basis", action="append", default=[], metavar="SPECIES=FILE")
    parser.add_argument("--ecutwfc", type=float, default=None)
    parser.add_argument("--scf-thr", type=float, default=None)
    parser.add_argument("--scf-nmax", type=int, default=None)
    parser.add_argument("--mixing-beta", type=float, default=None)
    parser.add_argument("--kpts", type=int, nargs=3, default=None, metavar=("NX", "NY", "NZ"))
    parser.add_argument("--optimizer", choices=["FIRE", "BFGS", "LBFGS"], default="FIRE")
    parser.add_argument("--mic", action="store_true")
    parser.add_argument("--no-climb", action="store_true")
    parser.add_argument("--resume", action="store_true", help="Resume from the latest complete chain in vcneb.traj")
    parser.add_argument("--resume-trajectory", default=None, help="Trajectory to resume from; defaults to workdir/vcneb.traj")
    return parser.parse_args()


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


def main() -> None:
    args = parse_args()
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    initial = read(args.initial)
    final = read(args.final)
    traj_path = workdir / "vcneb.traj"
    resume_path = Path(args.resume_trajectory).resolve() if args.resume_trajectory else traj_path
    if args.resume:
        images = read_chain_trajectory(resume_path, n_images=args.n_images)
        print(f"[OK] resumed latest complete {args.n_images}-image chain from {resume_path}")
    else:
        images = interpolate_vcneb(initial, final, n_images=args.n_images, align_cells=True, mic=args.mic)
    write(workdir / "initial-vcneb.traj", images)

    parameters = {
        "calculation": "scf",
        "basis_type": "lcao",
        "dft_functional": "pbe",
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
    optional_parameters = {
        "ecutwfc": args.ecutwfc,
        "scf_thr": args.scf_thr,
        "scf_nmax": args.scf_nmax,
        "mixing_beta": args.mixing_beta,
        "kpts": args.kpts,
    }
    parameters.update({key: value for key, value in optional_parameters.items() if value is not None})
    factory = make_ase_abacus_factory(parameters=parameters, command=args.command)
    attach_abacus_calculators(images, workdir=workdir, factory=factory)

    chain, _ = run_vcneb(
        images,
        pressure_gpa=args.pressure_gpa,
        k=args.k,
        climb=not args.no_climb,
        optimizer=args.optimizer,
        fmax=args.fmax,
        steps=args.steps,
        logfile=workdir / "vcneb.opt.log",
        trajectory=traj_path,
        trajectory_mode="a" if args.resume and resume_path == traj_path and traj_path.exists() else "w",
        snapshot_dir=workdir / "snapshots",
    )
    chain.plot_band(workdir / "vcneb_barrier.png")
    barrier, delta = chain.barrier()
    print(f"[DONE] barrier={barrier:.6f} eV delta={delta:.6f} eV workdir={workdir}")


if __name__ == "__main__":
    main()
