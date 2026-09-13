"""Regression check for bounded explicit line-search recovery."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import vcneb.core as core
from vcneb.core import interpolate_vcneb, run_vcneb


class QuadraticCalculator(Calculator):
    implemented_properties = ["energy", "forces", "stress"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x = np.asarray(atoms.positions, dtype=float)
        delta = x[:, 0] - 1.0
        self.results = {
            "energy": float(np.dot(delta, delta)),
            "forces": np.column_stack((-2.0 * delta, np.zeros((len(x), 2)))),
            "stress": np.zeros((3, 3), dtype=float),
        }


def main() -> None:
    cell = np.diag([5.0, 5.0, 5.0])
    initial = Atoms("Ar", positions=[[0.0, 0.0, 0.0]], cell=cell, pbc=True)
    final = Atoms("Ar", positions=[[3.0, 0.0, 0.0]], cell=cell, pbc=True)
    images = interpolate_vcneb(initial, final, n_images=3, align_cells=False)
    images[1].positions[0, 1] = 0.5
    for image in images:
        image.calc = QuadraticCalculator()

    original_step = core.BFGSLineSearch.step
    calls = {"count": 0}

    def fail_once(self, *args, **kwargs):
        calls["count"] += 1
        if calls["count"] == 1:
            raise RuntimeError("LineSearch failed!")
        return None

    core.BFGSLineSearch.step = fail_once
    try:
        _, optimizer = run_vcneb(
            images,
            climb=False,
            optimizer="BFGSLineSearch",
            optimizer_kwargs={"maxstep": 0.2},
            line_search_retries=1,
            line_search_retry_factor=0.5,
            fmax=1e-12,
            steps=1,
            logfile=None,
            trajectory=None,
        )
    finally:
        core.BFGSLineSearch.step = original_step

    if calls["count"] != 2:
        raise SystemExit(f"line-search retry count mismatch: {calls}")
    if int(getattr(optimizer, "nsteps", 0)) != 1:
        raise SystemExit("retry did not preserve the global optimizer step count")
    if int(getattr(optimizer, "line_search_retries_used", 0)) != 1:
        raise SystemExit("retry provenance was not attached to the returned optimizer")

    # A calculator/optimizer error carrying the same words must not trigger a
    # retry when the selected optimizer has no explicit line search.
    guarded = interpolate_vcneb(initial, final, n_images=3, align_cells=False)
    guarded[1].positions[0, 1] = 0.5
    for image in guarded:
        image.calc = QuadraticCalculator()
    original_fire_step = core.FIRE.step
    fire_calls = {"count": 0}

    def fire_failure(self, *args, **kwargs):
        fire_calls["count"] += 1
        raise RuntimeError("LineSearch failed! (calculator wrapper)")

    core.FIRE.step = fire_failure
    try:
        try:
            run_vcneb(
                guarded,
                climb=False,
                optimizer="FIRE",
                line_search_retries=2,
                fmax=1e-12,
                steps=1,
                logfile=None,
                trajectory=None,
            )
        except RuntimeError:
            pass
        else:
            raise SystemExit("non-line-search optimizer failure unexpectedly succeeded")
    finally:
        core.FIRE.step = original_fire_step
    if fire_calls["count"] != 1:
        raise SystemExit("non-line-search optimizer failure was retried")
    print("line_search_retry_regression=ok")


if __name__ == "__main__":
    main()
