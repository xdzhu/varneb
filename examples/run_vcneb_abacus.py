"""Run VC-NEB with ABACUS through a calculator factory.

The default path uses ``ase.calculators.abacus.Abacus`` when available.  If
your production environment uses a different ABACUS ASE adapter, replace
``make_ase_abacus_factory(...)`` with a custom factory that returns one
calculator per image directory.
"""

from __future__ import annotations

import argparse
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
    parser.add_argument("--resume", action="store_true", help="Resume from the latest complete chain in vcneb.traj")
    parser.add_argument("--resume-trajectory", default=None, help="Trajectory to resume from; defaults to workdir/vcneb.traj")
    return parser.parse_args()


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
        images = interpolate_vcneb(initial, final, n_images=args.n_images, align_cells=True)
    write(workdir / "initial-vcneb.traj", images)

    parameters = {
        "calculation": "scf",
        "basis_type": "lcao",
        "dft_functional": "pbe",
        "cal_force": 1,
        "cal_stress": 1,
        "out_stru": 1,
        # Fill in material-specific keys before a production run:
        # "ecutwfc": 100,
        # "kpts": [4, 4, 4],
        # "pp": {"P": "P.upf"},
        # "basis": {"P": "P_gga_10au_100Ry_2s2p1d.orb"},
        # "pseudo_dir": "...",
        # "basis_dir": "...",
    }
    factory = make_ase_abacus_factory(parameters=parameters, command=args.command)
    attach_abacus_calculators(images, workdir=workdir, factory=factory)

    chain, _ = run_vcneb(
        images,
        pressure_gpa=args.pressure_gpa,
        k=args.k,
        climb=True,
        optimizer="FIRE",
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
