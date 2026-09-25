"""Compare fixed-cell VC-NEB with ASE's standard CINEB on one model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from ase.mep import NEB
from ase.optimize import FIRE

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.run_toy_vcneb import ToyPhaseTransition
from vcneb import interpolate_vcneb, run_vcneb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(ROOT / "outputs" / "fixed_cell_ase_comparison.json"))
    parser.add_argument("--images", type=int, default=7)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--fmax", type=float, default=0.10)
    return parser.parse_args()


def build_endpoints() -> tuple[Atoms, Atoms, np.ndarray]:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final = Atoms("Ar", scaled_positions=[[0.75, 0.5, 0.5]], cell=reference_cell, pbc=True)
    return initial, final, reference_cell


def make_images(initial: Atoms, final: Atoms, reference_cell: np.ndarray, n_images: int) -> list[Atoms]:
    images = interpolate_vcneb(initial, final, n_images=n_images, align_cells=False)
    for index, image in enumerate(images[1:-1], start=1):
        scaled = image.get_scaled_positions(wrap=False)
        scaled[0, 1] += 0.04 * np.sin(index)
        image.set_scaled_positions(scaled)
        image.calc = ToyPhaseTransition(reference_cell)
    images[0].calc = ToyPhaseTransition(reference_cell)
    images[-1].calc = ToyPhaseTransition(reference_cell)
    return images


def summarize(energies: list[float], steps: int) -> dict:
    return {
        "energies_eV": energies,
        "barrier_eV": float(max(energies) - energies[0]),
        "reaction_energy_eV": float(energies[-1] - energies[0]),
        "optimizer_steps": int(steps),
    }


def main() -> None:
    args = parse_args()
    if args.images < 3:
        raise ValueError("--images must be at least 3")
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    initial, final, reference_cell = build_endpoints()

    ase_images = make_images(initial, final, reference_cell, args.images)
    ase_neb = NEB(ase_images, k=0.15, climb=True, parallel=False, method="improvedtangent")
    ase_opt = FIRE(ase_neb, logfile=None)
    ase_opt.run(fmax=args.fmax, steps=args.steps)
    ase_energies = [float(image.get_potential_energy()) for image in ase_images]

    vc_images = make_images(initial, final, reference_cell, args.images)
    vc_chain, vc_opt = run_vcneb(
        vc_images,
        k=0.15,
        climb=True,
        cell_mask=np.zeros((3, 3)),
        optimizer="FIRE",
        fmax=args.fmax,
        steps=args.steps,
        logfile=None,
        trajectory=None,
        snapshot_dir=None,
    )
    vc_energies = [float(value) for value in vc_chain.enthalpies]

    result = {
        "calculator": "analytic ToyPhaseTransition",
        "n_images": args.images,
        "fmax_eV_per_A": args.fmax,
        "fixed_cell": True,
        "ase_cineb": summarize(ase_energies, getattr(ase_opt, "nsteps", -1)),
        "vcneb_cell_mask_zero": summarize(vc_energies, getattr(vc_opt, "nsteps", -1)),
        "barrier_absolute_difference_eV": abs(
            (max(ase_energies) - ase_energies[0]) - (max(vc_energies) - vc_energies[0])
        ),
        "max_cell_deviation_from_reference_A": float(max(
            np.max(np.abs(image.cell.array - reference_cell))
            for image in [*ase_images, *vc_images]
        )),
        "max_corresponding_image_position_difference_A": float(max(
            np.max(np.abs(ase_image.positions - vc_image.positions))
            for ase_image, vc_image in zip(ase_images, vc_images)
        )),
        "saddle_diagnostics": vc_chain.saddle_diagnostics(),
    }
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
