"""Recoverable conditional evaluations at explicitly prescribed biaxial strain.

The physical map and stress pullback are the existing curvature chart's,
not a new deformation convention. The extra substrate coordinate is a
control, never an automatically relaxable internal degree of freedom.
"""
from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.stress import voigt_6_to_full_3x3_stress

from .biaxial_curvature import BiaxialClampedCurvatureCoordinates
from .clamped_mode_coordinates import ClampedModeCoordinates


class BiaxialClampedModeCoordinates(BiaxialClampedCurvatureCoordinates):
    """Atomic/open-cell Angstrom plus t=L_strain*(epsilon-epsilon_anchor).

    A fixed-Q solver must retain or explicitly freeze the controlled
    covector. CalculatorModeEvaluator forwards this requirement so the
    ordinary conditional optimizer rejects a releasing complement before
    any calculator call. No material mode, metric or strain range is chosen.
    """

    coordinate_convention = "reference_atomic_open_cell_A_plus_prescribed_biaxial_strain_A"

    def __init__(self, reference: Atoms, *, boundary, cell_scale_A: float,
                 anchor_strain: float, strain_scale_A: float):
        super().__init__(reference, boundary=boundary, cell_scale_A=cell_scale_A,
                         anchor_strain=anchor_strain, strain_scale_A=strain_scale_A)
        self._inverse = ClampedModeCoordinates(
            reference, boundary=boundary, cell_scale_A=cell_scale_A,
        )

    @property
    def coordinate_count(self) -> int:
        return self.size

    @property
    def n_atoms(self) -> int:
        return len(self.reference)

    @property
    def reference_atoms(self) -> Atoms:
        return self.reference

    @property
    def reference_cell(self) -> np.ndarray:
        return self.reference.cell.array.copy()

    @property
    def fixed_coordinate_covectors(self) -> np.ndarray:
        """Columns whose contractions with every released direction must vanish."""
        return self.controlled_direction()[:, None]

    @property
    def mechanical_contract(self) -> dict:
        contract = self._inverse.mechanical_contract
        contract.update(kind="controlled_biaxial_clamped_plane",
            anchor_strain=self.anchor_strain, strain_scale_A=self.strain_scale_A,
            controlled_coordinate_index=self.controlled_index,
            controlled_cell_rows=[0, 1],
            substrate_row_factor="(1+epsilon)/(1+anchor_strain)",
            fixed_coordinate_covectors=self.fixed_coordinate_covectors.tolist())
        return contract

    def rigid_translation_directions(self) -> np.ndarray:
        result = np.zeros((self.size, 3))
        result[:-1] = self._inverse.rigid_translation_directions()
        return result

    def to_atoms(self, coordinates) -> Atoms:
        return self.displaced(coordinates)

    def from_atoms(self, atoms: Atoms) -> np.ndarray:
        """Invert the declared ordered lift; do not wrap, rotate or rematch it."""
        if (not isinstance(atoms, Atoms) or atoms.constraints
                or atoms.get_chemical_symbols() != self.reference.get_chemical_symbols()
                or not np.array_equal(atoms.pbc, self.reference.pbc)):
            raise ValueError("ordered unconstrained atom/species/periodic identity required")
        rows = self.reference_cell[:2]
        factor = float(np.sum(atoms.cell.array[:2]*rows)/np.sum(rows*rows))
        if (not np.isfinite(factor) or factor <= 0.
                or not np.allclose(atoms.cell.array[:2], factor*rows, atol=1e-9, rtol=0)):
            raise ValueError("substrate rows must share the prescribed biaxial factor")
        restored = atoms.copy()
        cell = restored.cell.array.copy()
        cell[:2] = rows
        restored.set_cell(cell, scale_atoms=True)
        internal = self._inverse.from_atoms(restored)
        epsilon = (1+self.anchor_strain)*factor - 1
        return np.r_[internal, self.strain_scale_A*(epsilon-self.anchor_strain)]

    def energy_and_gradient(self, coordinates, atoms: Atoms, *, pressure_eV_per_A3: float = 0.):
        """Full d(E+PV)/dy, including substrate-reaction stress, in eV/Angstrom."""
        expected = self.to_atoms(coordinates)
        if (not isinstance(atoms, Atoms) or atoms.constraints
                or atoms.get_chemical_symbols() != expected.get_chemical_symbols()
                or not np.array_equal(atoms.pbc, expected.pbc)
                or not np.allclose(atoms.cell.array, expected.cell.array, atol=1e-9, rtol=0)
                or not np.allclose(atoms.positions, expected.positions, atol=1e-9, rtol=0)):
            raise ValueError("calculator geometry differs from the declared controlled point")
        pressure = float(pressure_eV_per_A3)
        if not np.isfinite(pressure):
            raise ValueError("pressure must be finite in eV/Angstrom^3")
        energy = float(atoms.get_potential_energy())
        forces = np.asarray(atoms.get_forces(apply_constraint=False), dtype=float)
        raw = np.asarray(atoms.calc.get_property("stress", atoms), dtype=float)
        stress = voigt_6_to_full_3x3_stress(raw) if raw.shape == (6,) else raw
        if (not np.isfinite(energy) or forces.shape != (self.n_atoms, 3)
                or not np.isfinite(forces).all() or stress.shape != (3, 3)
                or not np.isfinite(stress).all()
                or not np.allclose(stress, stress.T, atol=1e-8, rtol=0)):
            raise ValueError("finite energy/forces/symmetric raw ASE stress required")
        gradient = self.enthalpy_gradient(coordinates, atoms, forces, stress, pressure)
        return energy + pressure*atoms.get_volume(), gradient
