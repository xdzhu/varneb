"""Run a constrained VC-NEB path followed by an unconstrained refinement.

The example deliberately keeps a transverse displacement in the initial path.
The first stage projects updates onto the endpoint collective mode while the
transverse component is retained by ``constraint_mode='projected'``.  The
second stage releases the constraint and refines the same images in the full
extended space.  The calculator is analytic and is intended for workflow
validation, not for a material barrier.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.run_toy_vcneb import ToyPhaseTransition
from vcneb import Mode, build_mode_basis, mode_guided_path, run_vcneb


def attach_calculators(images: list[Atoms], reference_cell: np.ndarray) -> None:
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)


def main() -> None:
    output = ROOT / "outputs" / "release_and_refine_model"
    output.mkdir(parents=True, exist_ok=True)
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final = Atoms("Ar", scaled_positions=[[0.75, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final_cell = reference_cell.copy()
    final_cell[0, 0] *= 1.25
    final.set_cell(final_cell, scale_atoms=True)

    images = mode_guided_path(
        initial,
        final,
        n_images=7,
        mode=Mode([[0.0, 1.0, 0.0]]),
        amplitude=0.30,
        envelope="sin",
        align_cells=False,
    )
    attach_calculators(images, reference_cell)

    endpoint_mode = Mode([[2.5, 0.0, 0.0]], cell=np.diag([0.25, 0.0, 0.0]))
    basis = build_mode_basis(endpoint_mode, initial)
    constrained, constrained_opt = run_vcneb(
        images,
        k=0.15,
        climb=True,
        mode_basis=basis,
        constraint_mode="projected",
        optimizer="FIRE",
        fmax=0.01,
        steps=500,
        logfile=output / "constrained.log",
        trajectory=output / "constrained.traj",
        snapshot_dir=output / "constrained_snapshots",
    )

    refined_images = [image.copy() for image in constrained.images]
    attach_calculators(refined_images, reference_cell)
    refined, refined_opt = run_vcneb(
        refined_images,
        k=0.15,
        climb=True,
        optimizer="FIRE",
        fmax=0.01,
        steps=500,
        logfile=output / "refined.log",
        trajectory=output / "refined.traj",
        snapshot_dir=output / "refined_snapshots",
    )

    def result(chain, optimizer) -> dict:
        barrier, reaction = chain.barrier()
        return {
            "barrier_eV": barrier,
            "reaction_energy_eV": reaction,
            "max_generalized_force_eV_per_A": chain.gradient_norm(-chain.get_forces()),
            "optimizer_steps": int(getattr(optimizer, "nsteps", -1)),
            "saddle_diagnostics": chain.saddle_diagnostics(),
            "path_diagnostics": chain.path_diagnostics(),
        }

    report = {
        "calculator": "analytic ToyPhaseTransition",
        "workflow": "projected collective-mode VC-NEB -> unconstrained full-space refinement",
        "constrained": result(constrained, constrained_opt),
        "refined": result(refined, refined_opt),
    }
    (output / "release_and_refine.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
