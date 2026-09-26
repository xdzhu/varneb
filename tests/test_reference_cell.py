"""Finite-difference checks for the backend-independent atomic/strain chart."""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.singlepoint import SinglePointCalculator

from vcneb.reference_cell import ReferenceCellCoordinates


class AffineQuadraticCalculator(Calculator):
    implemented_properties = ["energy", "forces", "stress"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        positions = atoms.positions
        cell = atoms.cell.array
        volume = atoms.get_volume()
        k, c = 0.7, 0.04
        energy = 0.5 * k * np.sum(positions**2) + 0.5 * c * np.sum(cell**2)
        forces = -k * positions
        # At fixed fractional coordinates, a spatial affine strain acts on
        # both atomic positions and lattice vectors from the right.
        stress = (k * positions.T @ positions + c * cell.T @ cell) / volume
        self.results = {"energy": energy, "forces": forces, "stress": stress}


def _chart() -> ReferenceCellCoordinates:
    reference = Atoms(
        "HeLi", scaled_positions=[[0.13, 0.21, 0.32], [0.61, 0.53, 0.44]],
        cell=[[4.0, 0.0, 0.0], [0.3, 5.0, 0.0], [0.1, 0.2, 6.0]], pbc=True,
    )
    return ReferenceCellCoordinates(reference)


def test_reference_cell_roundtrip_with_shear_and_periodic_wrap() -> None:
    chart = _chart()
    translations = chart.rigid_translation_directions()
    assert translations.shape == (chart.coordinate_count, 3)
    assert np.allclose(translations[:6], np.tile(np.eye(3), (2, 1)))
    assert np.allclose(translations[-6:], 0.0)
    coordinates = np.array([
        0.12, -0.08, 0.04, -0.11, 0.05, 0.07,
        0.02, -0.03, 0.06, 0.012, -0.015, 0.025,
    ])
    atoms = chart.to_atoms(coordinates)
    assert np.allclose(chart.from_atoms(atoms), coordinates, atol=1e-12)
    # Equivalent wrapped positions must represent the same reference-gauge u.
    scaled = atoms.get_scaled_positions(wrap=False)
    scaled[0, 0] += 1.0
    atoms.set_scaled_positions(scaled)
    assert np.allclose(chart.from_atoms(atoms), coordinates, atol=1e-12)
    centered, removed = chart.remove_mass_weighted_translation(coordinates, [4.0, 8.0])
    assert np.allclose(
        np.average(centered[:6].reshape(2, 3), axis=0, weights=[4.0, 8.0]),
        np.zeros(3), atol=1e-12,
    )
    assert np.allclose(centered[-6:], coordinates[-6:])
    assert np.allclose(centered[:6].reshape(2, 3) + removed, coordinates[:6].reshape(2, 3))


def test_energy_stress_and_force_gradient_matches_all_finite_differences() -> None:
    chart = _chart()
    coordinates = np.array([
        0.12, -0.08, 0.04, -0.11, 0.05, 0.07,
        0.02, -0.03, 0.06, 0.012, -0.015, 0.025,
    ])
    pressure = 0.003  # eV/Angstrom^3, positive compression

    def value(values: np.ndarray) -> float:
        atoms = chart.to_atoms(values)
        atoms.calc = AffineQuadraticCalculator()
        return atoms.get_potential_energy() + pressure * atoms.get_volume()

    atoms = chart.to_atoms(coordinates)
    atoms.calc = AffineQuadraticCalculator()
    energy, gradient = chart.energy_and_gradient(coordinates, atoms, pressure_eV_per_A3=pressure)
    assert energy == pytest.approx(value(coordinates))
    finite_difference = np.empty(chart.coordinate_count)
    step = 1e-6
    for index in range(chart.coordinate_count):
        plus, minus = coordinates.copy(), coordinates.copy()
        plus[index] += step
        minus[index] -= step
        finite_difference[index] = (value(plus) - value(minus)) / (2.0 * step)
    assert np.allclose(gradient, finite_difference, rtol=2e-7, atol=2e-7)


def test_invalid_cell_gauge_and_mismatched_calculator_structure_fail_early() -> None:
    chart = _chart()
    zero = np.zeros(chart.coordinate_count)
    with pytest.raises(ValueError, match="positive definite"):
        chart.to_atoms(np.r_[np.zeros(3 * chart.n_atoms), [-1.1, 0.0, 0.0, 0.0, 0.0, 0.0]])
    rotated = chart.to_atoms(zero)
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    rotated.set_cell(chart.reference_cell @ rotation, scale_atoms=True)
    with pytest.raises(ValueError, match="rotation or nonsymmetric"):
        chart.from_atoms(rotated)
    shifted = chart.to_atoms(zero)
    shifted.positions[0, 0] += 0.1
    shifted.calc = AffineQuadraticCalculator()
    with pytest.raises(ValueError, match="positions disagree"):
        chart.energy_and_gradient(zero, shifted)
    wrapped = chart.to_atoms(zero)
    wrapped.positions[0] += wrapped.cell.array[0]
    wrapped.calc = SinglePointCalculator(
        wrapped, energy=1.0, forces=np.zeros((2, 3)), stress=np.zeros(6),
    )
    energy, gradient = chart.energy_and_gradient(zero, wrapped)
    assert energy == pytest.approx(1.0)
    assert np.allclose(gradient, 0.0)
