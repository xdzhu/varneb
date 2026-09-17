"""Relax one physical VCA endpoint with externally coupled atomic/cell forces."""

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
from ase.optimize import FIRE

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.vasp import attach_vasp_calculators, default_vasp_command, vasp_input_fingerprints


def _write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def _git_revision() -> str | None:
    declared = os.environ.get("VCNEB_GIT_REVISION")
    if declared:
        return declared
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", required=True, help="physical five-site endpoint structure")
    parser.add_argument("--source-dir", required=True, help="directory containing INCAR/KPOINTS/POTCAR")
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--virtual-symbol", default=None, help="physical symbol representing one VCA site")
    parser.add_argument("--components", nargs="+", default=None, help="coincident VASP components for the VCA site")
    parser.add_argument("--symprec", type=float, default=1e-4)
    parser.add_argument("--fmax", type=float, default=0.03)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--maxstep", type=float, default=0.08)
    parser.add_argument("--ncores", type=int, default=int(os.environ.get("NP", "40")))
    parser.add_argument("--vasp-bin", default=os.environ.get("VASP_BIN", "vasp_std"))
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source_dir = Path(args.source_dir).resolve()
    workdir = Path(args.workdir).resolve()
    current = workdir / "CONTCAR.current"
    if workdir.exists() and not args.resume:
        raise FileExistsError(f"workdir already exists; pass --resume to continue: {workdir}")
    if args.resume and not current.is_file():
        raise FileNotFoundError(f"resume structure is missing: {current}")
    workdir.mkdir(parents=True, exist_ok=True)

    structure_path = current if args.resume else Path(args.structure).resolve()
    atoms = read(structure_path)
    atoms.set_constraint(FixSymmetry(atoms, symprec=args.symprec))
    command = default_vasp_command(args.ncores, args.vasp_bin)
    vca_options = {}
    if args.virtual_symbol is not None or args.components is not None:
        if args.virtual_symbol is None or not args.components:
            raise ValueError("pass both --virtual-symbol and --components for a VCA endpoint")
        vca_options = {"vca_virtual_symbol": args.virtual_symbol, "vca_components": args.components}
    attach_vasp_calculators(
        [atoms], source_dir=source_dir, workdir=workdir / "vasp", command=command,
        overrides={"xc": "PBE", "pp": "PBE"},
        **vca_options,
    )

    filtered = FrechetCellFilter(atoms, scalar_pressure=0.0)
    optimizer = FIRE(
        filtered,
        restart=str(workdir / "fire_restart.json"),
        logfile=str(workdir / "relax.log"),
        trajectory=str(workdir / "relax.traj"),
        maxstep=args.maxstep,
    )
    state = {"steps": 0}

    def checkpoint() -> None:
        state["steps"] = optimizer.get_number_of_steps()
        write(current, atoms, format="vasp", direct=True, vasp5=True)
        payload = {
            "status": "running",
            "optimizer_steps": state["steps"],
            "potential_energy_eV": float(atoms.get_potential_energy()),
            "max_atomic_force_eV_per_A": float(np.linalg.norm(atoms.get_forces(), axis=1).max()),
            "stress_eV_per_A3_voigt": np.asarray(atoms.get_stress(), dtype=float).tolist(),
            "cell_A": np.asarray(atoms.cell, dtype=float).tolist(),
            "volume_A3": float(atoms.get_volume()),
        }
        _write_json_atomic(workdir / "endpoint_relax_state.json", payload)

    optimizer.attach(checkpoint, interval=1)
    converged = bool(optimizer.run(fmax=args.fmax, steps=args.steps))
    checkpoint()
    forces = np.asarray(atoms.get_forces(), dtype=float)
    stress = np.asarray(atoms.get_stress(), dtype=float)
    write(workdir / "CONTCAR", atoms, format="vasp", direct=True, vasp5=True)
    summary = {
        "status": "completed" if converged else "step_limit",
        "converged": converged,
        "git_revision": _git_revision(),
        "optimizer": "FIRE(FrechetCellFilter)",
        "symmetry_constraint": {"name": "FixSymmetry", "symprec": args.symprec},
        "fmax_target_eV_per_A": args.fmax,
        "optimizer_steps": optimizer.get_number_of_steps(),
        "potential_energy_eV": float(atoms.get_potential_energy()),
        "max_atomic_force_eV_per_A": float(np.linalg.norm(forces, axis=1).max()),
        "forces_eV_per_A": forces.tolist(),
        "stress_eV_per_A3_voigt": stress.tolist(),
        "cell_A": np.asarray(atoms.cell, dtype=float).tolist(),
        "volume_A3": float(atoms.get_volume()),
        "vca": ({"virtual_symbol": args.virtual_symbol, "components": args.components} if args.virtual_symbol else None),
        "licensed_input_fingerprints": vasp_input_fingerprints(source_dir),
    }
    _write_json_atomic(workdir / "endpoint_relax_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if converged else 1


if __name__ == "__main__":
    raise SystemExit(main())
