"""Relax one variable-cell endpoint through any VARNEB ASE factory.

This is the calculator-independent preparation stage for the production
template.  It deliberately separates endpoint BFGS relaxation from VC-NEB:
the relaxed structure is written to ``CONTCAR`` and can then be passed through
``examples/run_vcneb_ase.py --static-only`` and the endpoint gate.
"""

from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
import sys

import numpy as np
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import BFGS

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import (  # noqa: E402
    attach_image_calculators,
    endpoint_structure_record,
    validate_image_calculators,
)

GPA_PER_EV_PER_A3 = 160.21766208
EV_A3_TO_KBAR = 1602.176634


def _stress_residual_kbar(atoms, pressure_gpa: float) -> float:
    """Maximum Cartesian stress residual relative to hydrostatic pressure."""

    stress = _finite_stress(atoms)
    target = -pressure_gpa / GPA_PER_EV_PER_A3
    residual = stress - np.diag([target, target, target])
    return float(np.max(np.abs(residual)) * EV_A3_TO_KBAR)


class EndpointBFGS(BFGS):
    """BFGS with separate physical force and pressure convergence tests.

    ASE cell filters express cell gradients in force-like units, so a single
    ``fmax`` can stop with a backend- and volume-dependent residual stress.
    The endpoint contract instead applies ``fmax`` to atomic forces and an
    explicit kbar tolerance to the stress tensor.  Cached calculator results
    are reused, so this check does not add electronic-structure evaluations.
    """

    def __init__(self, *args, endpoint_atoms, pressure_gpa, stress_kbar=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.endpoint_atoms = endpoint_atoms
        self.pressure_gpa = float(pressure_gpa)
        self.stress_kbar = None if stress_kbar is None else float(stress_kbar)

    def gradient_converged(self, gradient):
        if self.stress_kbar is None:
            return super().gradient_converged(gradient)
        forces = np.asarray(self.endpoint_atoms.get_forces(), dtype=float)
        force_max = float(np.linalg.norm(forces, axis=1).max())
        return force_max < self.fmax and _stress_residual_kbar(
            self.endpoint_atoms, self.pressure_gpa
        ) < self.stress_kbar


def _load_symbol(spec: str):
    if ":" not in spec:
        raise ValueError("symbol must use module:attribute form")
    module_name, attribute = spec.split(":", 1)
    if not module_name or not attribute:
        raise ValueError(f"invalid symbol specification {spec!r}")
    return getattr(importlib.import_module(module_name), attribute)


def _json_object(path: str | None) -> dict:
    if path is None:
        return {}
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"JSON object expected in {path}")
    return value


def _finite_stress(atoms) -> np.ndarray:
    stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
    if stress.shape != (3, 3) or not np.all(np.isfinite(stress)):
        raise RuntimeError("calculator did not return a finite 3x3 stress tensor")
    return stress


def _snapshot(
    atoms,
    filtered,
    optimizer,
    pressure_gpa: float,
    fmax: float,
    stress_kbar: float | None = None,
) -> dict:
    forces = np.asarray(atoms.get_forces(), dtype=float)
    if forces.shape != (len(atoms), 3) or not np.all(np.isfinite(forces)):
        raise RuntimeError("calculator did not return finite atomic forces")
    stress = _finite_stress(atoms)
    pressure_ev_a3 = pressure_gpa / GPA_PER_EV_PER_A3
    energy = float(atoms.get_potential_energy())
    generalized = np.asarray(filtered.get_forces(), dtype=float)
    return {
        "optimizer": "BFGS(FrechetCellFilter)",
        "optimizer_steps": optimizer.get_number_of_steps(),
        "external_pressure_gpa": pressure_gpa,
        "fmax_target_eV_per_A": fmax,
        "potential_energy_eV": energy,
        "enthalpy_eV": energy + pressure_ev_a3 * float(atoms.get_volume()),
        "max_atomic_force_eV_per_A": float(np.linalg.norm(forces, axis=1).max()),
        "max_generalized_force_eV_per_A": float(np.linalg.norm(generalized, axis=1).max()),
        "max_abs_stress_residual_kbar": _stress_residual_kbar(atoms, pressure_gpa),
        "stress_target_kbar": stress_kbar,
        "stress_eV_per_A3": stress.tolist(),
        "cell_A": np.asarray(atoms.cell, dtype=float).tolist(),
        "volume_A3": float(atoms.get_volume()),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--factory", required=True, help="module:attribute VARNEB ASE factory")
    parser.add_argument("--parameters", default=None, help="JSON calculator parameter object")
    parser.add_argument("--factory-kwargs", default=None, help="JSON factory-only keyword object")
    parser.add_argument("--command", default=None)
    parser.add_argument("--pressure-gpa", type=float, default=0.0)
    parser.add_argument("--fmax", type=float, default=0.10)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--maxstep", type=float, default=0.05)
    parser.add_argument(
        "--stress-kbar",
        type=float,
        default=None,
        help="variable-cell convergence threshold for residual stress",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--fixed-cell", action="store_true", help="relax atoms only, retaining the cell")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.fmax <= 0.0 or args.steps < 1 or args.maxstep <= 0.0:
        raise ValueError("fmax, steps and maxstep must be positive")
    if args.stress_kbar is not None and args.stress_kbar <= 0.0:
        raise ValueError("stress-kbar must be positive")
    if args.fixed_cell and args.stress_kbar is not None:
        raise ValueError("stress-kbar cannot be enforced during fixed-cell relaxation")
    workdir = Path(args.workdir).resolve()
    current = workdir / "CONTCAR.current"
    if workdir.exists() and not args.resume and any(workdir.iterdir()):
        raise FileExistsError(f"workdir is not empty; use --resume: {workdir}")
    if args.resume and not current.is_file():
        raise FileNotFoundError(f"missing endpoint checkpoint: {current}")
    workdir.mkdir(parents=True, exist_ok=True)

    atoms = read(current if args.resume else Path(args.structure).resolve())
    builder = _load_symbol(args.factory)
    factory_kwargs = _json_object(args.factory_kwargs)
    factory_kwargs["parameters"] = _json_object(args.parameters)
    if args.command is not None:
        factory_kwargs["command"] = args.command
    factory = builder(**factory_kwargs)
    attach_image_calculators([atoms], workdir=workdir / "calculator", factory=factory)
    validate_image_calculators(
        [atoms],
        require_stress=True,
        require_variable_cell=not args.fixed_cell,
        require_directory=True,
        require_unique_directories=True,
    )

    filtered = atoms if args.fixed_cell else FrechetCellFilter(
        atoms,
        scalar_pressure=args.pressure_gpa / GPA_PER_EV_PER_A3,
    )
    optimizer = EndpointBFGS(
        filtered,
        logfile=str(workdir / "relax.log"),
        trajectory=str(workdir / "relax.traj"),
        maxstep=args.maxstep,
        endpoint_atoms=atoms,
        pressure_gpa=args.pressure_gpa,
        stress_kbar=args.stress_kbar,
    )

    def checkpoint() -> None:
        payload = _snapshot(
            atoms,
            filtered,
            optimizer,
            args.pressure_gpa,
            args.fmax,
            args.stress_kbar,
        )
        payload.update({"status": "running", "endpoint": endpoint_structure_record(atoms)})
        write(workdir / "CONTCAR.current", atoms, format="vasp", direct=True, vasp5=True)
        (workdir / "endpoint_relax_state.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    optimizer.attach(checkpoint, interval=1)
    converged = bool(optimizer.run(fmax=args.fmax, steps=args.steps))
    checkpoint()
    payload = _snapshot(
        atoms,
        filtered,
        optimizer,
        args.pressure_gpa,
        args.fmax,
        args.stress_kbar,
    )
    payload.update(
        {
            "status": "completed" if converged else "step_limit",
            "converged": converged,
            "endpoint": endpoint_structure_record(atoms),
            "structure": str(workdir / "CONTCAR"),
        }
    )
    write(workdir / "CONTCAR", atoms, format="vasp", direct=True, vasp5=True)
    (workdir / "endpoint_relax_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if converged else 1


if __name__ == "__main__":
    raise SystemExit(main())
