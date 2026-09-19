"""Run a compact VC-NEB convergence matrix on the analytic toy calculator."""

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(ROOT / "outputs" / "vcneb_convergence_matrix.json"))
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--fmax", type=float, default=0.10)
    parser.add_argument("--images", type=int, nargs="+", default=[5, 7, 9])
    parser.add_argument("--springs", type=float, nargs="+", default=[0.05, 0.10, 0.20])
    parser.add_argument("--cell-scales", type=float, nargs="+", default=[4.0, 5.0, 6.0])
    parser.add_argument(
        "--optimizers",
        nargs="+",
        choices=["FIRE", "BlockFIRE", "SplitFIRE", "LBFGS"],
        default=["FIRE", "BlockFIRE", "SplitFIRE", "LBFGS"],
    )
    return parser.parse_args()


def build_endpoints(reference_cell: np.ndarray):
    initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final = Atoms("Ar", scaled_positions=[[0.75, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final_cell = reference_cell.copy()
    final_cell[0, 0] *= 1.25
    final.set_cell(final_cell, scale_atoms=True)
    return initial, final


def main() -> None:
    args = parse_args()
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    reference_cell = np.diag([5.0, 5.0, 5.0])
    results = []
    for n_images, spring, cell_scale, optimizer in itertools.product(
        args.images, args.springs, args.cell_scales, args.optimizers
    ):
        initial, final = build_endpoints(reference_cell)
        images = interpolate_vcneb(initial, final, n_images=n_images, align_cells=False)
        # A linear interpolation lies exactly on this toy MEP.  Add the same
        # smooth transverse perturbation at every resolution so the optimizer
        # scaling test measures relaxation work instead of a zero-step path.
        for index, image in enumerate(images[1:-1], start=1):
            scaled = image.get_scaled_positions(wrap=False)
            scaled[0, 1] += 0.04 * np.sin(np.pi * index / (n_images - 1))
            image.set_scaled_positions(scaled)
        for image in images:
            image.calc = ToyPhaseTransition(reference_cell)
        chain, opt = run_vcneb(
            images,
            k=spring,
            cell_scale=cell_scale,
            climb=False,
            optimizer=optimizer,
            fmax=args.fmax,
            steps=args.steps,
            logfile=None,
            trajectory=None,
            snapshot_dir=None,
        )
        final_force = chain.gradient_norm(-chain.get_forces())
        barrier, delta = chain.barrier()
        results.append(
            {
                "n_images": n_images,
                "spring_eV_per_A2": spring,
                "cell_scale_A": cell_scale,
                "optimizer": optimizer,
                "steps_requested": args.steps,
                "optimizer_steps": int(getattr(opt, "nsteps", -1)),
                "fmax_target_eV_per_A": args.fmax,
                "final_max_generalized_force_eV_per_A": final_force,
                "barrier_eV": barrier,
                "barrier_error_eV": abs(barrier - 0.25),
                "reaction_energy_eV": delta,
                "converged": bool(final_force < args.fmax),
            }
        )

    report = {
        "calculator": "analytic ToyPhaseTransition",
        "expected_barrier_eV": 0.25,
        "matrix_size": len(results),
        "settings": {
            "images": args.images,
            "springs_eV_per_A2": args.springs,
            "cell_scales_A": args.cell_scales,
            "optimizers": args.optimizers,
            "steps": args.steps,
            "fmax_eV_per_A": args.fmax,
        },
        "results": results,
    }
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
