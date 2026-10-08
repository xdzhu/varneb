"""Controlled biaxial strain is an input parameter, not a relaxed cell DOF.

This sampling chart augments the existing clamped-plane internal coordinates
with one scaled substrate parameter. It is not an optimizer/filter and does
not certify a stationary branch or a barrier prediction. Every evaluated
geometry must match its declared coordinate point before gradients are used.
"""

from __future__ import annotations

from numbers import Real

import numpy as np
from ase import Atoms

from .active_curvature import ActiveJointCurvatureCoordinates
from .epitaxial_boundary import ClampedPlaneBoundary, cell_work_derivative


class BiaxialClampedCurvatureCoordinates:
    """Reference-atom Å, open-cell Å, and controlled t=L*(epsilon-epsilon0).

    The first two row-cell vectors are H0[:2]*(1+epsilon)/(1+epsilon0).
    The third vector uses the existing exact clamped-plane deformation basis.
    Fractional atom shifts are measured with H0 throughout. L is explicit so
    the final gradient is d(E+PV)/dt in eV/Å; d(E+PV)/dε is L times that entry.
    Only ``internal_relaxation_basis`` is eligible for orthogonal release.
    """

    def __init__(
        self, reference: Atoms, *, boundary: ClampedPlaneBoundary,
        cell_scale_A: float, anchor_strain: float, strain_scale_A: float,
    ) -> None:
        if not isinstance(boundary, ClampedPlaneBoundary):
            raise TypeError("boundary must be a ClampedPlaneBoundary")
        if (isinstance(anchor_strain, (bool, np.bool_))
                or not isinstance(anchor_strain, Real) or not np.isfinite(anchor_strain)
                or anchor_strain <= -1):
            raise ValueError("anchor strain must be finite and greater than -1")
        if (isinstance(strain_scale_A, (bool, np.bool_))
                or not isinstance(strain_scale_A, Real) or not np.isfinite(strain_scale_A)
                or strain_scale_A <= 0):
            raise ValueError("strain scale must be finite and positive in Angstrom")
        if (isinstance(cell_scale_A, (bool, np.bool_))
                or not isinstance(cell_scale_A, Real) or not np.isfinite(cell_scale_A)
                or cell_scale_A <= 0):
            raise ValueError("cell scale must be finite and positive in Angstrom")
        owned = reference.copy()  # no inherited calculator or energy cache
        boundary.validate_images([owned])
        self._internal = ActiveJointCurvatureCoordinates.for_clamped_plane(
            owned, cell_scale_A=cell_scale_A, boundary=boundary,
        )
        self._anchor_strain = float(anchor_strain)
        self._strain_scale_A = float(strain_scale_A)
        self._boundary = boundary
        self._cell0 = owned.cell.array.copy()
        self._cell0.flags.writeable = False
        derivatives = np.array([
            self._cell0 @ b.T / cell_scale_A
            for b in self._internal.deformation_basis
        ])
        # The geometry construction fixes these rows exactly, not approximately.
        derivatives[:, :2] = 0.
        prescribed = np.zeros((3, 3))
        prescribed[:2] = self._cell0[:2] / ((1 + self.anchor_strain) * self.strain_scale_A)
        self._cell_derivatives = np.concatenate((derivatives, prescribed[None]))
        self._cell_derivatives.flags.writeable = False

    @property
    def anchor_strain(self) -> float:
        return self._anchor_strain

    @property
    def reference(self) -> Atoms:
        return self._internal.reference.copy()

    @property
    def cell_scale_A(self) -> float:
        return float(self._internal.cell_scale_A)

    @property
    def strain_scale_A(self) -> float:
        return self._strain_scale_A

    @property
    def internal_size(self) -> int:
        return self._internal.size

    @property
    def size(self) -> int:
        return self.internal_size + 1

    @property
    def controlled_index(self) -> int:
        return self.size - 1

    def _coordinates(self, delta: np.ndarray) -> tuple[np.ndarray, float]:
        values = np.asarray(delta, dtype=float)
        if values.shape != (self.size,) or not np.isfinite(values).all():
            raise ValueError("finite coordinates of the declared augmented size required")
        strain = self.anchor_strain + values[-1] / self.strain_scale_A
        if not np.isfinite(strain) or strain <= -1:
            raise ValueError("controlled strain must be finite and greater than -1")
        return values, float(strain)

    def displaced(self, delta: np.ndarray) -> Atoms:
        """Generate a calculator-free point at its prescribed substrate strain."""
        values, strain = self._coordinates(delta)
        result = self._internal.displaced(values[:-1])
        cell = result.cell.array.copy()
        cell[:2] = self._cell0[:2] * ((1 + strain) / (1 + self.anchor_strain))
        result.set_cell(cell, scale_atoms=True)
        if (not np.isfinite(result.positions).all()
                or not np.isfinite(cell).all() or np.linalg.det(cell) <= 1e-12
                or np.linalg.cond(cell) > 1e12):
            raise ValueError("controlled probe has invalid geometry")
        result.calc = None
        return result

    def enthalpy_gradient(
        self, delta: np.ndarray, atoms: Atoms, forces_eV_per_A: np.ndarray,
        stress_eV_per_A3: np.ndarray, pressure_eV_per_A3: float,
    ) -> np.ndarray:
        """Physical full gradient at the declared point; no calculator calls.

        A strain component is a frozen partial derivative unless the internal
        branch is stationary. Stress must include substrate reactions: dropping
        them would remove the very response this controlled coordinate measures.
        """
        expected = self.displaced(delta)
        forces = np.asarray(forces_eV_per_A, dtype=float)
        if (atoms.get_chemical_symbols() != expected.get_chemical_symbols()
                or not np.array_equal(atoms.pbc, expected.pbc)
                or atoms.constraints or forces.shape != (len(expected), 3)
                or not np.isfinite(forces).all()):
            raise ValueError("ordered identity, periodicity and finite full atomic forces required")
        if (not np.allclose(atoms.cell.array, expected.cell.array, atol=1e-9, rtol=0)
                or not np.allclose(atoms.positions, expected.positions, atol=1e-9, rtol=0)):
            raise ValueError("evaluated geometry differs from the declared controlled point; no remapping")
        atom_gradient = -forces @ atoms.cell.array.T @ np.linalg.inv(self._cell0.T)
        cell_gradient = [
            cell_work_derivative(stress_eV_per_A3, atoms.cell.array, derivative,
                                 pressure=pressure_eV_per_A3)
            for derivative in self._cell_derivatives
        ]
        return np.concatenate((atom_gradient.ravel(), cell_gradient))

    def internal_relaxation_basis(self) -> np.ndarray:
        """Remove rigid translations, keeping the controlled coordinate fixed."""
        internal = self._internal.translation_free_basis()
        result = np.zeros((self.size, internal.shape[1]))
        result[:-1] = internal
        return result

    def controlled_direction(self) -> np.ndarray:
        """The prescribed input column; retain it, never eliminate/optimize it."""
        result = np.zeros(self.size)
        result[-1] = 1.
        return result
