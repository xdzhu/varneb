"""Smoke test for the variable-cell NEB implementation.

The calculator below is analytic and cheap.  It has two minima in an extended
space made of one fractional coordinate and one cell strain component, with a
known saddle between them.  It is not a physical potential; it is a sign and
optimizer sanity check.
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import interpolate_vcneb, run_vcneb
from vcneb.core import deformation_from_cell


class ToyPhaseTransition(Calculator):
    implemented_properties = ["energy", "forces", "stress"]

    def __init__(self, reference_cell, *, barrier=0.25, k_perp=3.0, k_rest=1.0):
        super().__init__()
        self.reference_cell = np.asarray(reference_cell, dtype=float)
        self.barrier = float(barrier)
        self.k_perp = float(k_perp)
        self.k_rest = float(k_rest)
        self.a = np.array([0.25, 0.0])
        self.b = np.array([0.75, 0.25])
        self.mid = 0.5 * (self.a + self.b)
        self.direction = self.b - self.a
        self.length = float(np.linalg.norm(self.direction))
        self.unit = self.direction / self.length
        self.perp = np.array([-self.unit[1], self.unit[0]])

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        atoms = self.atoms
        cell = atoms.cell.array
        deform = deformation_from_cell(cell, self.reference_cell)
        q = atoms.get_scaled_positions(wrap=False)

        y = np.array([q[0, 0], deform[0, 0] - 1.0])
        rel = y - self.mid
        u = np.dot(rel, self.unit) / (0.5 * self.length)
        v = np.dot(rel, self.perp)

        energy = self.barrier * (u * u - 1.0) ** 2 + 0.5 * self.k_perp * v * v
        grad_y = (
            4.0 * self.barrier * u * (u * u - 1.0) * self.unit / (0.5 * self.length)
            + self.k_perp * v * self.perp
        )

        q_target = np.array([0.5, 0.5])
        rest = q[0, 1:] - q_target
        energy += 0.5 * self.k_rest * float(np.dot(rest, rest))

        grad_q = np.zeros_like(q)
        grad_q[0, 0] = grad_y[0]
        grad_q[0, 1:] = self.k_rest * rest

        grad_deform = np.zeros((3, 3))
        grad_deform[0, 0] = grad_y[1]
        off_target = deform - np.eye(3)
        off_target[0, 0] = 0.0
        energy += 0.5 * self.k_rest * float(np.sum(off_target * off_target))
        grad_deform += self.k_rest * off_target

        frac_force = -grad_q
        cart_forces = frac_force @ np.linalg.inv(cell.T)

        cell_force = -grad_deform
        volume = atoms.get_volume()
        virial = cell_force @ deform.T
        stress = -virial / volume

        self.results["energy"] = float(energy)
        self.results["forces"] = cart_forces
        self.results["stress"] = stress


def make_endpoint(qx: float, exx: float) -> Atoms:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    atoms = Atoms("Ar", scaled_positions=[[qx, 0.5, 0.5]], cell=reference_cell, pbc=True)
    deform = np.eye(3)
    deform[0, 0] = 1.0 + exx
    atoms.set_cell(reference_cell @ deform.T, scale_atoms=True)
    return atoms


def main() -> None:
    out = Path("toy_vcneb_run")
    out.mkdir(exist_ok=True)

    initial = make_endpoint(0.25, 0.0)
    final = make_endpoint(0.75, 0.25)
    images = interpolate_vcneb(initial, final, n_images=7, align_cells=False)

    # Start slightly off the exact MEP so the relaxation has work to do.
    for index, image in enumerate(images[1:-1], start=1):
        scaled = image.get_scaled_positions(wrap=False)
        scaled[0, 1] += 0.04 * np.sin(index)
        image.set_scaled_positions(scaled)

    reference_cell = initial.cell.array
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)

    chain, _ = run_vcneb(
        images,
        k=0.15,
        climb=True,
        optimizer="FIRE",
        fmax=0.01,
        steps=300,
        logfile=out / "toy-vcneb.log",
        trajectory=out / "toy-vcneb.traj",
        snapshot_dir=out / "snapshots",
    )
    chain.write_chain(out / "final-chain.traj")
    chain.plot_band(out / "toy-band.png")
    barrier, delta = chain.barrier()
    print(f"barrier_eV={barrier:.6f}")
    print(f"delta_eV={delta:.6f}")
    print("enthalpies_eV=" + " ".join(f"{e:.6f}" for e in chain.enthalpies))
    print(f"outputs={out.resolve()}")


if __name__ == "__main__":
    main()
