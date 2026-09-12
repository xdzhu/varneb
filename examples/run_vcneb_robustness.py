"""Stress the analytic VCNEB path with deterministic random perturbations."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.run_toy_vcneb import ToyPhaseTransition
from vcneb import interpolate_vcneb, run_vcneb


def build_endpoints(reference_cell: np.ndarray) -> tuple[Atoms, Atoms]:
    initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final = Atoms("Ar", scaled_positions=[[0.75, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final_cell = reference_cell.copy()
    final_cell[0, 0] *= 1.25
    final.set_cell(final_cell, scale_atoms=True)
    return initial, final


def run_case(interpolation: str, seed: int, *, fmax: float, steps: int) -> dict:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial, final = build_endpoints(reference_cell)
    images = interpolate_vcneb(
        initial,
        final,
        n_images=7,
        align_cells=False,
        cell_interpolation=interpolation,
    )
    rng = np.random.default_rng(seed)
    for image in images[1:-1]:
        scaled = image.get_scaled_positions(wrap=False)
        scaled[0, 1:] += 0.08 * rng.normal(size=2)
        image.set_scaled_positions(scaled)
        cell = image.cell.array.copy()
        cell[1, 0] += 0.08 * rng.normal()
        cell[2, 0] += 0.08 * rng.normal()
        image.set_cell(cell, scale_atoms=False)
        image.calc = ToyPhaseTransition(reference_cell)
    initial.calc = ToyPhaseTransition(reference_cell)
    final.calc = ToyPhaseTransition(reference_cell)

    chain, optimizer = run_vcneb(
        images,
        k=0.15,
        climb=True,
        optimizer="FIRE",
        fmax=fmax,
        steps=steps,
        logfile=None,
        trajectory=None,
        snapshot_dir=None,
    )
    barrier, reaction = chain.barrier()
    final_force = chain.gradient_norm(-chain.get_forces())
    return {
        "cell_interpolation": interpolation,
        "seed": seed,
        "barrier_eV": barrier,
        "barrier_error_eV": abs(barrier - 0.25),
        "reaction_energy_eV": reaction,
        "final_max_generalized_force_eV_per_A": final_force,
        "optimizer_steps": int(getattr(optimizer, "nsteps", -1)),
        "converged": bool(final_force < fmax),
        "saddle": chain.saddle_diagnostics(),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(ROOT / "outputs" / "vcneb_robustness.json"))
    parser.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3])
    parser.add_argument("--fmax", type=float, default=0.002)
    parser.add_argument("--steps", type=int, default=800)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results = [
        run_case(interpolation, seed, fmax=args.fmax, steps=args.steps)
        for interpolation, seed in itertools.product(("linear", "log_strain"), args.seeds)
    ]
    report = {
        "calculator": "analytic ToyPhaseTransition",
        "expected_barrier_eV": 0.25,
        "settings": {
            "n_images": 7,
            "seeds": args.seeds,
            "fmax_eV_per_A": args.fmax,
            "steps": args.steps,
        },
        "case_count": len(results),
        "converged_cases": sum(result["converged"] for result in results),
        "max_barrier_error_eV": max(result["barrier_error_eV"] for result in results),
        "max_final_generalized_force_eV_per_A": max(
            result["final_max_generalized_force_eV_per_A"] for result in results
        ),
        "results": results,
    }
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
