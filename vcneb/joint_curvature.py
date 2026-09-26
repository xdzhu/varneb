"""Local atomic--symmetric-strain enthalpy coordinates for saddle diagnostics.

These coordinates deliberately exclude cell rotations and retain the same
fractional-atom/deformation convention as :mod:`vcneb.core`.  A Hessian at a
nonstationary image is only a *local curvature diagnostic*, never by itself a
transition-state certificate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from ase import Atoms

from .core import cell_force_from_arrays, fractional_force_from_arrays


def symmetric_strain_basis() -> np.ndarray:
    """Return six Frobenius-orthonormal 3x3 symmetric matrices."""

    basis = np.zeros((6, 3, 3), dtype=float)
    for axis in range(3):
        basis[axis, axis, axis] = 1.0
    for index, (i, j) in enumerate(((1, 2), (0, 2), (0, 1)), start=3):
        basis[index, i, j] = basis[index, j, i] = 1.0 / np.sqrt(2.0)
    return basis


@dataclass(frozen=True)
class JointCurvatureCoordinates:
    """Atomic reference positions (Å) plus six scaled symmetric strains (Å)."""

    reference: Atoms
    cell_scale_A: float

    def __post_init__(self) -> None:
        cell = np.asarray(self.reference.cell.array)
        if (not np.isfinite(cell).all() or np.linalg.det(cell) <= 0
                or not np.isfinite(self.cell_scale_A) or self.cell_scale_A <= 0):
            raise ValueError("reference cell and cell scale must be finite and positive")
        if len(self.reference) < 2:
            raise ValueError("at least two atoms are needed to remove translations")

    @property
    def size(self) -> int:
        return 3 * len(self.reference) + 6

    def displaced(self, delta: np.ndarray) -> Atoms:
        """Apply a joint perturbation without mutating the reference."""

        shift = np.asarray(delta, dtype=float)
        if shift.shape != (self.size,) or not np.isfinite(shift).all():
            raise ValueError("joint displacement has invalid shape or non-finite values")
        cell0 = self.reference.cell.array
        q = self.reference.get_scaled_positions(wrap=False).copy()
        q += shift[: 3 * len(self.reference)].reshape(-1, 3) @ np.linalg.inv(cell0)
        deform = np.eye(3) + np.einsum(
            "a,aij->ij", shift[3 * len(self.reference) :] / self.cell_scale_A,
            symmetric_strain_basis(),
        )
        result = self.reference.copy()
        result.set_cell(cell0 @ deform.T, scale_atoms=False)
        result.set_scaled_positions(q)
        if np.linalg.det(result.cell.array) <= 0:
            raise ValueError("joint displacement produced a non-positive cell volume")
        result.calc = None
        return result

    def enthalpy_gradient(
        self, atoms: Atoms, forces_eV_per_A: np.ndarray,
        stress_eV_per_A3: np.ndarray, pressure_eV_per_A3: float,
    ) -> np.ndarray:
        """Return d(E+PV)/dy, using ASE tensile-positive Cartesian stress."""

        n = len(self.reference)
        forces = np.asarray(forces_eV_per_A, dtype=float)
        stress = np.asarray(stress_eV_per_A3, dtype=float)
        if (len(atoms) != n or atoms.get_chemical_symbols() != self.reference.get_chemical_symbols()
                or forces.shape != (n, 3) or stress.shape != (3, 3)
                or not np.isfinite(forces).all() or not np.isfinite(stress).all()
                or not np.isfinite(pressure_eV_per_A3)):
            raise ValueError("invalid atomic force, stress, pressure or atom identity")
        cell0 = self.reference.cell.array
        fractional = fractional_force_from_arrays(forces, atoms.cell.array)
        atom_force = fractional @ np.linalg.inv(cell0.T)
        deformation_force = cell_force_from_arrays(
            stress, atoms, cell0, pressure=pressure_eV_per_A3,
        )
        strain_force = np.einsum(
            "ij,aij->a", deformation_force, symmetric_strain_basis(),
        ) / self.cell_scale_A
        return -np.concatenate((atom_force.ravel(), strain_force))

    def translation_free_basis(self) -> np.ndarray:
        """Orthonormal basis for the 3N+6 space with three translations removed."""

        n = len(self.reference)
        translation = np.zeros((self.size, 3))
        for atom in range(n):
            translation[3 * atom : 3 * atom + 3] = np.eye(3) / np.sqrt(n)
        _, _, vt = np.linalg.svd(translation.T, full_matrices=True)
        return vt[3:].T


def central_difference_hessian(
    plus_gradients: np.ndarray, minus_gradients: np.ndarray, step_A: float,
) -> tuple[np.ndarray, float]:
    """Return symmetric Hessian and pre-symmetrization reciprocity defect."""

    plus = np.asarray(plus_gradients, dtype=float)
    minus = np.asarray(minus_gradients, dtype=float)
    if (plus.ndim != 2 or plus.shape != minus.shape or plus.shape[0] != plus.shape[1]
            or not np.isfinite(plus).all() or not np.isfinite(minus).all()
            or not np.isfinite(step_A) or step_A <= 0):
        raise ValueError("paired gradients must be finite square arrays with positive step")
    raw = (plus - minus).T / (2.0 * step_A)
    defect = float(np.linalg.norm(raw - raw.T) / max(np.linalg.norm(raw), 1e-30))
    return 0.5 * (raw + raw.T), defect
