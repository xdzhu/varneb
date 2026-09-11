"""VC-NEB smoke test on mapped 12-atom HfO2 T -> PO endpoints.

This uses a synthetic endpoint double-well calculator.  It is not a physical
HfO2 potential; it validates the VC-NEB code on realistic atom counts, mapping,
and variable-cell geometry before spending DFT time.
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import interpolate_vcneb, run_vcneb
from vcneb.core import deformation_from_cell


FIXTURE = ROOT / "validation" / "hfo2_t_to_po"


def state_x(atoms, reference_cell, cell_scale):
    q = atoms.get_scaled_positions(wrap=False)
    deform = deformation_from_cell(atoms.cell.array, reference_cell)
    return np.concatenate([
        (q @ reference_cell).reshape(-1),
        (cell_scale * (deform - np.eye(3))).reshape(-1),
    ])


class EndpointDoubleWell(Calculator):
    implemented_properties = ["energy", "forces", "stress"]

    def __init__(self, initial, final, *, barrier=1.0, k_perp=1.0):
        super().__init__()
        self.reference_cell = initial.cell.array.copy()
        self.cell_scale = abs(np.linalg.det(self.reference_cell)) ** (1.0 / 3.0)
        self.x0 = state_x(initial, self.reference_cell, self.cell_scale)
        self.x1 = state_x(final, self.reference_cell, self.cell_scale)
        self.mid = 0.5 * (self.x0 + self.x1)
        self.direction = self.x1 - self.x0
        self.length = float(np.linalg.norm(self.direction))
        self.unit = self.direction / self.length
        self.barrier = float(barrier)
        self.k_perp = float(k_perp)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        atoms = self.atoms
        x = state_x(atoms, self.reference_cell, self.cell_scale)
        rel = x - self.mid
        u = float(np.dot(rel, self.unit) / (0.5 * self.length))
        parallel = np.dot(rel, self.unit) * self.unit
        perp = rel - parallel

        energy = self.barrier * (u * u - 1.0) ** 2 + 0.5 * self.k_perp * float(np.dot(perp, perp))
        grad_x = (
            4.0 * self.barrier * u * (u * u - 1.0) * self.unit / (0.5 * self.length)
            + self.k_perp * perp
        )
        force_x = -grad_x

        natoms = len(atoms)
        f_atoms_x = force_x[: 3 * natoms].reshape(natoms, 3)
        f_cell_x = force_x[3 * natoms :].reshape(3, 3)

        frac_force = f_atoms_x @ self.reference_cell.T
        cart_forces = frac_force @ np.linalg.inv(atoms.cell.array.T)

        deform = deformation_from_cell(atoms.cell.array, self.reference_cell)
        cell_force = f_cell_x * self.cell_scale
        virial = cell_force @ deform.T
        stress = -virial / atoms.get_volume()

        self.results["energy"] = float(energy)
        self.results["forces"] = cart_forces
        self.results["stress"] = stress


def main() -> None:
    out = FIXTURE / "model_vcneb_run"
    out.mkdir(parents=True, exist_ok=True)
    initial = read(FIXTURE / "T_HfO2_12.vasp")
    final = read(FIXTURE / "PO_HfO2_12_mapped.vasp")
    images = interpolate_vcneb(initial, final, n_images=7, align_cells=True, mic=True)
    final_unwrapped = images[-1].copy()

    # Push the initial guess slightly off the analytic MEP.
    for i, image in enumerate(images[1:-1], start=1):
        q = image.get_scaled_positions(wrap=False)
        q[4:, 2] += 0.01 * np.sin(i)
        image.set_scaled_positions(q)

    calc = EndpointDoubleWell(initial, final_unwrapped, barrier=0.8, k_perp=0.5)
    for image in images:
        image.calc = calc

    chain, _ = run_vcneb(
        images,
        k=0.2,
        climb=True,
        optimizer="FIRE",
        fmax=0.03,
        steps=500,
        logfile=out / "hfo2-model-vcneb.log",
        trajectory=out / "hfo2-model-vcneb.traj",
        snapshot_dir=out / "snapshots",
    )
    write(out / "final-chain.traj", chain.images)
    chain.plot_band(out / "hfo2-model-band.png")
    barrier, delta = chain.barrier()
    print(f"barrier_eV={barrier:.6f}")
    print(f"delta_eV={delta:.6f}")
    print("enthalpies_eV=" + " ".join(f"{e:.6f}" for e in chain.enthalpies))
    print(f"outputs={out.resolve()}")


if __name__ == "__main__":
    main()
