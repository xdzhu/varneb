"""Compare unconstrained, mode-guided, and strict-mode VCNEB on the toy model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.run_toy_vcneb import ToyPhaseTransition
from vcneb import Mode, build_mode_basis, interpolate_vcneb, mode_guided_path, run_vcneb


def endpoints() -> tuple[Atoms, Atoms, np.ndarray]:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final = Atoms("Ar", scaled_positions=[[0.75, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final_cell = reference_cell.copy()
    final_cell[0, 0] *= 1.25
    final.set_cell(final_cell, scale_atoms=True)
    return initial, final, reference_cell


def attach(images: list[Atoms], reference_cell: np.ndarray) -> None:
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)


def optimize(
    images: list[Atoms],
    reference_cell: np.ndarray,
    *,
    mode_basis: np.ndarray | None = None,
    fmax: float = 0.002,
    steps: int = 600,
) -> dict:
    attach(images, reference_cell)
    chain, optimizer = run_vcneb(
        images,
        k=0.15,
        climb=True,
        mode_basis=mode_basis,
        optimizer="FIRE",
        fmax=fmax,
        steps=steps,
        logfile=None,
        trajectory=None,
        snapshot_dir=None,
    )
    barrier, reaction = chain.barrier()
    saddle = chain.saddle_diagnostics()
    return {
        "barrier_eV": barrier,
        "reaction_energy_eV": reaction,
        "final_max_generalized_force_eV_per_A": chain.gradient_norm(-chain.get_forces()),
        "optimizer_steps": int(getattr(optimizer, "nsteps", -1)),
        "saddle": saddle,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(ROOT / "outputs" / "mode_path_variants.json"))
    parser.add_argument("--fmax", type=float, default=0.002)
    parser.add_argument("--steps", type=int, default=600)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    initial, final, reference_cell = endpoints()

    direct = interpolate_vcneb(initial, final, n_images=7, align_cells=False)
    for index, image in enumerate(direct[1:-1], start=1):
        scaled = image.get_scaled_positions(wrap=False)
        scaled[0, 1] += 0.20 * np.sin(index)
        image.set_scaled_positions(scaled)

    guided = mode_guided_path(
        initial,
        final,
        n_images=7,
        mode=Mode([[0.0, 1.0, 0.0]]),
        amplitude=0.20,
        envelope="sin",
        align_cells=False,
    )

    strict = interpolate_vcneb(initial, final, n_images=7, align_cells=False)
    coupled_mode = Mode([[1.0, 0.0, 0.0]], cell=np.diag([0.1, 0.0, 0.0]))
    strict_basis = build_mode_basis(coupled_mode, initial, cell_scale=5.0)

    result = {
        "calculator": "analytic ToyPhaseTransition",
        "n_images": 7,
        "expected_barrier_eV": 0.25,
        "fmax_eV_per_A": args.fmax,
        "max_steps": args.steps,
        "variants": {
            "ordinary_unconstrained": optimize(direct, reference_cell, fmax=args.fmax, steps=args.steps),
            "mode_guided_then_released": optimize(guided, reference_cell, fmax=args.fmax, steps=args.steps),
            "strict_coupled_mode": optimize(
                strict,
                reference_cell,
                mode_basis=strict_basis,
                fmax=args.fmax,
                steps=args.steps,
            ),
        },
    }
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
