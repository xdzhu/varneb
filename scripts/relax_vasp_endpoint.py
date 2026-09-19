"""Relax one VASP endpoint with BFGS in a specified external pressure field.

VASP remains a static energy/force/stress calculator (IBRION=-1, NSW=0,
ISIF=2).  ASE BFGS updates atoms and the cell through ``FrechetCellFilter``;
this keeps endpoint and VARNEB image energies on the same calculator surface.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import BFGS

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.vasp import attach_vasp_calculators, default_vasp_command, vasp_input_fingerprints

GPA_PER_EV_PER_A3 = 160.21766208


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", required=True)
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--pressure-gpa", type=float, required=True)
    parser.add_argument("--fmax", type=float, default=0.02)
    parser.add_argument("--steps", type=int, default=120)
    parser.add_argument("--maxstep", type=float, default=0.05)
    parser.add_argument("--symprec", type=float, default=1.0e-4)
    parser.add_argument("--ncores", type=int, default=int(os.environ.get("NP", "40")))
    parser.add_argument("--vasp-bin", default=os.environ.get("VASP_BIN", "vasp_std"))
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def git_revision() -> str | None:
    if declared := os.environ.get("VCNEB_GIT_REVISION"):
        return declared
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
        ).stdout.strip() or None
    except (OSError, subprocess.CalledProcessError):
        return None


def main() -> int:
    args = parse_args()
    source = Path(args.source_dir).resolve()
    workdir = Path(args.workdir).resolve()
    current = workdir / "CONTCAR.current"
    if workdir.exists() and not args.resume:
        raise FileExistsError(f"workdir already exists; use --resume: {workdir}")
    if args.resume and not current.is_file():
        raise FileNotFoundError(f"missing checkpoint: {current}")
    for name in ("INCAR", "KPOINTS", "POTCAR"):
        if not (source / name).is_file():
            raise FileNotFoundError(source / name)
    workdir.mkdir(parents=True, exist_ok=True)

    atoms = read(current if args.resume else Path(args.structure).resolve())
    atoms.set_constraint(FixSymmetry(atoms, symprec=args.symprec))
    attach_vasp_calculators(
        [atoms],
        source_dir=source,
        workdir=workdir / "vasp",
        command=default_vasp_command(args.ncores, args.vasp_bin),
        overrides={"xc": "PBE", "pp": "PBE"},
    )
    pressure_ev_a3 = args.pressure_gpa / GPA_PER_EV_PER_A3
    filtered = FrechetCellFilter(atoms, scalar_pressure=pressure_ev_a3)
    optimizer = BFGS(
        filtered,
        logfile=str(workdir / "relax.log"),
        trajectory=str(workdir / "relax.traj"),
        maxstep=args.maxstep,
    )

    def checkpoint() -> None:
        forces = np.asarray(atoms.get_forces(), dtype=float)
        energy = float(atoms.get_potential_energy())
        write(current, atoms, format="vasp", direct=True, vasp5=True)
        atomic_json(
            workdir / "endpoint_relax_state.json",
            {
                "status": "running",
                "optimizer": "BFGS(FrechetCellFilter)",
                "optimizer_steps": optimizer.get_number_of_steps(),
                "external_pressure_gpa": args.pressure_gpa,
                "potential_energy_eV": energy,
                "enthalpy_eV": energy + pressure_ev_a3 * float(atoms.get_volume()),
                "max_atomic_force_eV_per_A": float(np.linalg.norm(forces, axis=1).max()),
                "max_generalized_force_eV_per_A": float(np.linalg.norm(filtered.get_forces(), axis=1).max()),
                "stress_eV_per_A3_voigt": np.asarray(atoms.get_stress(), dtype=float).tolist(),
                "cell_A": np.asarray(atoms.cell, dtype=float).tolist(),
                "volume_A3": float(atoms.get_volume()),
            },
        )

    optimizer.attach(checkpoint, interval=1)
    converged = bool(optimizer.run(fmax=args.fmax, steps=args.steps))
    checkpoint()
    forces = np.asarray(atoms.get_forces(), dtype=float)
    energy = float(atoms.get_potential_energy())
    write(workdir / "CONTCAR", atoms, format="vasp", direct=True, vasp5=True)
    summary = {
        "status": "completed" if converged else "step_limit",
        "converged": converged,
        "git_revision": git_revision(),
        "optimizer": "BFGS(FrechetCellFilter)",
        "external_pressure_gpa": args.pressure_gpa,
        "fmax_target_eV_per_A": args.fmax,
        "optimizer_steps": optimizer.get_number_of_steps(),
        "potential_energy_eV": energy,
        "enthalpy_eV": energy + pressure_ev_a3 * float(atoms.get_volume()),
        "max_atomic_force_eV_per_A": float(np.linalg.norm(forces, axis=1).max()),
        "max_generalized_force_eV_per_A": float(np.linalg.norm(filtered.get_forces(), axis=1).max()),
        "forces_eV_per_A": forces.tolist(),
        "stress_eV_per_A3_voigt": np.asarray(atoms.get_stress(), dtype=float).tolist(),
        "cell_A": np.asarray(atoms.cell, dtype=float).tolist(),
        "volume_A3": float(atoms.get_volume()),
        "licensed_input_fingerprints": vasp_input_fingerprints(source),
    }
    atomic_json(workdir / "endpoint_relax_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if converged else 1


if __name__ == "__main__":
    raise SystemExit(main())
