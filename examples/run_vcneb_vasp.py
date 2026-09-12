"""Run VC-NEB with VASP through ASE.

This script is intentionally close to ``run_NEB/run_ase_cineb.sh`` but allows
the cell to vary along the band.  It reads relaxed endpoint structures and
copies INCAR/KPOINTS/POTCAR settings from the initial endpoint directory.
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
from vcneb.vasp import attach_vasp_calculators, default_vasp_command


def read_endpoint(directory: Path):
    for name in ["CONTCAR", "CONTCAR.vasp", "POSCAR"]:
        path = directory / name
        if path.exists():
            atoms = read(path)
            print(f"[OK] read {path}")
            return atoms
    raise FileNotFoundError(f"No endpoint structure found in {directory}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--initial", default=str(root / "initial_state" / "relax"))
    parser.add_argument("--final", default=str(root / "final_state" / "relax"))
    parser.add_argument("--workdir", default=str(root / "run_VCNEB" / "run"))
    parser.add_argument("--n-images", type=int, default=7, help="Total images including endpoints")
    parser.add_argument("--fmax", type=float, default=0.05)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--k", type=float, default=0.10)
    parser.add_argument("--pressure-gpa", type=float, default=0.0)
    parser.add_argument("--optimizer", choices=["FIRE", "BFGS", "LBFGS"], default="FIRE")
    parser.add_argument("--vasp-bin", default=os.environ.get("VASP_BIN", "vasp_std"))
    parser.add_argument("--ncores", type=int, default=int(os.environ.get("NP", "8")))
    parser.add_argument("--mic", action="store_true")
    parser.add_argument("--no-climb", action="store_true")
    parser.add_argument("--resume", action="store_true", help="Resume from the latest complete chain in vcneb.traj")
    parser.add_argument("--resume-trajectory", default=None, help="Trajectory to resume from; defaults to workdir/vcneb.traj")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    initial_dir = Path(args.initial).resolve()
    final_dir = Path(args.final).resolve()
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    initial = read_endpoint(initial_dir)
    final = read_endpoint(final_dir)

    traj_path = workdir / "vcneb.traj"
    resume_path = Path(args.resume_trajectory).resolve() if args.resume_trajectory else traj_path
    if args.resume:
        images = read_chain_trajectory(resume_path, n_images=args.n_images)
        print(f"[OK] resumed latest complete {args.n_images}-image chain from {resume_path}")
    else:
        images = interpolate_vcneb(initial, final, n_images=args.n_images, align_cells=True, mic=args.mic)
    write(workdir / "initial-vcneb.traj", images)

    command = default_vasp_command(args.ncores, args.vasp_bin)
    attach_vasp_calculators(
        images,
        source_dir=initial_dir,
        workdir=workdir,
        command=command,
        overrides={"xc": "PBE", "pp": "PBE"},
    )

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

    for image_index, image in enumerate(chain.images):
        write(workdir / f"{image_index:02d}" / "POSCAR.final", image, format="vasp", direct=True, vasp5=True)
    chain.plot_band(workdir / "vcneb_barrier.png")
    barrier, delta = chain.barrier()
    with open(workdir / "vcneb_summary.txt", "w", encoding="utf-8") as handle:
        handle.write(f"Forward barrier (enthalpy) = {barrier:.8f} eV\n")
        handle.write(f"Reaction enthalpy          = {delta:.8f} eV\n")
        handle.write("Image enthalpies (eV)      = " + " ".join(f"{e:.8f}" for e in chain.enthalpies) + "\n")
    print(f"[DONE] barrier={barrier:.6f} eV delta={delta:.6f} eV workdir={workdir}")


if __name__ == "__main__":
    main()
