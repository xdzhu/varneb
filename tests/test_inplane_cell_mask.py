"""A same-cell endpoint pair need not have a fixed-cell minimum path."""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from vcneb import interpolate_vcneb, run_vcneb


class StrainCoupledSlidingToy(Calculator):
    """One sliding coordinate coupled to in-plane strain, with fixed vacuum."""

    implemented_properties = ["energy", "forces", "stress"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        a0 = 6.0
        a = float(atoms.cell[0, 0])
        q = float(atoms.positions[0, 0]) * a0 / a - 3.0
        eta = a / a0 - 1.0
        envelope = 1.0 - q * q
        offset = eta - 0.2 * envelope
        energy = envelope**2 + 10.0 * offset**2
        d_energy_d_q = -4.0 * q * envelope + 8.0 * q * offset
        d_energy_d_eta = 20.0 * offset
        forces = np.zeros((1, 3))
        forces[0, 0] = -d_energy_d_q * a0 / a
        stress = np.zeros((3, 3))
        stress[0, 0] = d_energy_d_eta * (a / a0) / atoms.get_volume()
        self.results = {"energy": float(energy), "forces": forces, "stress": stress}


def _endpoints():
    cell = np.diag([6.0, 6.0, 20.0])
    initial = Atoms("Ar", positions=[[2.0, 3.0, 10.0]], cell=cell, pbc=True)
    final = Atoms("Ar", positions=[[4.0, 3.0, 10.0]], cell=cell, pbc=True)
    return initial, final


def _run(cell_mask):
    initial, final = _endpoints()
    images = interpolate_vcneb(initial, final, 7, align_cells=False)
    for image in images:
        image.calc = StrainCoupledSlidingToy()
    chain, _ = run_vcneb(
        images, cell_mask=cell_mask, optimizer="FIRE", fmax=0.005,
        steps=300, logfile=None, trajectory=None, snapshot_dir=None,
    )
    return chain


def test_inplane_relaxation_changes_barrier_without_moving_vacuum() -> None:
    fixed = _run(np.zeros((3, 3)))
    inplane = _run(np.diag([1.0, 0.0, 0.0]))

    fixed_barrier = float(np.max(fixed.enthalpies) - fixed.enthalpies[0])
    inplane_barrier = float(np.max(inplane.enthalpies) - inplane.enthalpies[0])
    assert fixed_barrier == pytest.approx(1.4, abs=0.01)
    assert inplane_barrier == pytest.approx(1.0, abs=0.01)
    assert all(np.array_equal(image.cell.array, fixed.images[0].cell.array)
               for image in fixed.images)
    assert all(image.cell[2, 2] == 20.0 for image in inplane.images)
    assert inplane.images[3].cell[0, 0] == pytest.approx(7.2, abs=0.05)
    assert inplane_barrier < fixed_barrier - 0.3
