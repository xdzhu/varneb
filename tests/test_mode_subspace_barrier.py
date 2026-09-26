"""Strict two-mode NEB barrier versus fully released NEB on a coupled toy."""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from vcneb import (
    ModePlane, ReferenceCellCoordinates, interpolate_vcneb, run_vcneb,
    strict_mode_subspace_for_vcneb,
)


class CoupledModeDoubleWell(Calculator):
    implemented_properties = ["energy", "forces", "stress"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x, y, z = atoms.positions[0] - 3.0
        envelope = 1.0 - x**2
        offset = z - 0.3 * envelope
        energy = envelope**2 + 3.0 * offset**2 + 0.1 * y**2
        gradient = np.array([
            -4.0 * x * envelope + 3.6 * x * offset,
            0.2 * y,
            6.0 * offset,
        ])
        self.results = {
            "energy": float(energy),
            "forces": -gradient.reshape(1, 3),
            "stress": np.zeros((3, 3)),
        }


def _images(initial: Atoms, final: Atoms):
    images = interpolate_vcneb(initial, final, 7, align_cells=False)
    for image in images:
        x = image.positions[0, 0] - 3.0
        image.positions[0, 2] = 3.0 + 0.3 * (1.0 - x**2)
        image.calc = CoupledModeDoubleWell()
    return images


def test_strict_two_mode_barrier_is_higher_than_released_path() -> None:
    reference = Atoms("Ar", positions=[[3.0, 3.0, 3.0]], cell=np.diag([6.0] * 3), pbc=True)
    chart = ReferenceCellCoordinates(reference)
    weights = np.zeros((9, 2))
    weights[0, 0], weights[1, 1] = 1.0, 1.0
    plane = ModePlane(
        reference=np.zeros(9), basis=np.eye(9), metric_weights=np.ones(9),
        axis_weights=weights, axis_labels=("reaction", "transverse"),
        amplitude_unit="toy length", reference_id="toy/coupled-mode-neb",
    )
    initial = chart.to_atoms(plane.frozen_coordinates([-1.0, 0.0]))
    final = chart.to_atoms(plane.frozen_coordinates([1.0, 0.0]))
    subspace = strict_mode_subspace_for_vcneb(
        plane, chart, initial, final, cell_scale_A=6.0,
    )
    strict, _ = run_vcneb(
        _images(initial, final), mode_basis=subspace.mode_basis,
        constraint_mode="subspace", cell_mask=np.zeros((3, 3)),
        optimizer="FIRE", fmax=0.005, steps=200, logfile=None,
        trajectory=None, snapshot_dir=None,
    )
    released, _ = run_vcneb(
        _images(initial, final), cell_mask=np.zeros((3, 3)),
        optimizer="FIRE", fmax=0.005, steps=200, logfile=None,
        trajectory=None, snapshot_dir=None,
    )
    released_seed = [image.copy() for image in strict.images]
    for image in released_seed:
        image.calc = CoupledModeDoubleWell()
    released_from_strict, _ = run_vcneb(
        released_seed, cell_mask=np.zeros((3, 3)),
        optimizer="FIRE", fmax=0.005, steps=200, logfile=None,
        trajectory=None, snapshot_dir=None,
    )
    strict_barrier = float(np.max(strict.enthalpies) - strict.enthalpies[0])
    released_barrier = float(np.max(released.enthalpies) - released.enthalpies[0])
    released_from_strict_barrier = float(
        np.max(released_from_strict.enthalpies) - released_from_strict.enthalpies[0]
    )
    assert strict_barrier == pytest.approx(1.27, abs=0.01)
    assert released_barrier == pytest.approx(1.0, abs=0.03)
    assert released_from_strict_barrier == pytest.approx(released_barrier, abs=0.03)
    assert strict_barrier > released_barrier + 0.20
    assert all(abs(image.positions[0, 2] - 3.0) < 1e-8 for image in strict.images)
    diagnostics = strict.path_diagnostics()
    assert diagnostics["constraint_mode"] == "subspace"
    assert diagnostics["n_mode_directions"] == 2
    midpoint = diagnostics["images"][3]
    assert midpoint["neb_residual_generalized_force_eV_per_A"] < 0.01
    assert midpoint["raw_force_orthogonal_to_mode_subspace_max_vector_eV_per_A"] > 1.0
    assert "raw_force_orthogonal_to_mode_subspace_norm_eV_per_A" not in (
        released.path_diagnostics()["images"][3]
    )
