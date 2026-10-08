"""Joint curvature coordinates in the actual permitted deformation subspace.

The historical six-symmetric-strain implementation stays byte-identical for
its published evidence. This extension reuses its atomic gradient and gauge,
while constructing cell probes in a separately declared mechanical subspace.
Neither class certifies stationarity, identifies phonons or executes DFT.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np
from ase import Atoms

from .core import cell_force_from_arrays
from .epitaxial_boundary import ClampedPlaneBoundary
from .joint_curvature import JointCurvatureCoordinates, symmetric_strain_basis


@dataclass(frozen=True)
class ActiveJointCurvatureCoordinates(JointCurvatureCoordinates):
    """Reference atomic displacements and scaled deformation coefficients, in Å.

    An explicit Frobenius-orthonormal basis encodes the allowed cell variables.
    Zero directions mean fixed cell; the default matches the legacy six strains.
    ``for_clamped_plane`` shares the exact endpoint/NEB boundary, whose permitted
    tilts are generally nonsymmetric, not symmetric shear strains. All atoms
    remain movable: extra ASE atomic constraints are refused, not inherited.
    """

    deformation_basis: np.ndarray | None = field(default=None, kw_only=True)
    candidate_validator: Callable[[Sequence[Atoms]], None] | None = field(
        default=None, kw_only=True, repr=False,
    )

    def __post_init__(self) -> None:
        super().__post_init__()
        basis = np.array(
            symmetric_strain_basis() if self.deformation_basis is None else self.deformation_basis,
            dtype=float, copy=True,
        )
        if (basis.ndim != 3 or basis.shape[1:] != (3, 3)
                or not np.isfinite(basis).all()):
            raise ValueError("deformation basis must be finite matrices with shape (k, 3, 3)")
        flat = basis.reshape(len(basis), 9)
        if not np.allclose(flat @ flat.T, np.eye(len(basis)), rtol=0, atol=1e-10):
            raise ValueError("deformation basis must be Frobenius-orthonormal")
        if self.candidate_validator is not None and not callable(self.candidate_validator):
            raise TypeError("candidate_validator must be callable")
        basis.flags.writeable = False
        object.__setattr__(self, "deformation_basis", basis)
        self._validate_candidate(self.reference)

    @classmethod
    def for_clamped_plane(
        cls, reference: Atoms, *, cell_scale_A: float, boundary: ClampedPlaneBoundary,
    ) -> ActiveJointCurvatureCoordinates:
        """Use precisely the one/three cell directions of a validated substrate.

        Reference and later probes must already satisfy the boundary. They are
        not projected, remapped or recalculated. Changing the allowed subspace
        defines a different mechanical experiment.
        """
        if not isinstance(boundary, ClampedPlaneBoundary):
            raise TypeError("boundary must be a ClampedPlaneBoundary")
        boundary.validate_images([reference])
        count = 3 * len(reference)
        basis = boundary.mode_basis[count:, count:].T.reshape(-1, 3, 3)
        return cls(reference, cell_scale_A, deformation_basis=basis,
                   candidate_validator=boundary.validate_images)

    @property
    def cell_dofs(self) -> int:
        return len(self.deformation_basis)

    @property
    def size(self) -> int:
        return 3 * len(self.reference) + self.cell_dofs

    def _validate_candidate(self, atoms: Atoms) -> None:
        cell = np.asarray(atoms.cell.array)
        if (not np.isfinite(cell).all() or np.linalg.det(cell) <= 0
                or not np.isfinite(atoms.positions).all()):
            raise ValueError("joint displacement produced a non-positive or non-finite geometry")
        if atoms.constraints:
            raise ValueError("extra ASE constraints require an explicit atomic subspace")
        if self.candidate_validator is not None:
            self.candidate_validator([atoms])

    def displaced(self, delta: np.ndarray) -> Atoms:
        """Construct one calculator-free probe; do not mutate the reference."""
        shift = np.asarray(delta, dtype=float)
        if shift.shape != (self.size,) or not np.isfinite(shift).all():
            raise ValueError("joint displacement has invalid shape or non-finite values")
        cell0 = self.reference.cell.array
        q = self.reference.get_scaled_positions(wrap=False).copy()
        q += shift[:3 * len(self.reference)].reshape(-1, 3) @ np.linalg.inv(cell0)
        deform = np.eye(3) + np.einsum(
            "a,aij->ij", shift[3 * len(self.reference):] / self.cell_scale_A,
            self.deformation_basis,
        )
        result = self.reference.copy()
        result.set_cell(cell0 @ deform.T, scale_atoms=False)
        result.set_scaled_positions(q)
        result.calc = None
        self._validate_candidate(result)
        return result

    def enthalpy_gradient(
        self, atoms: Atoms, forces_eV_per_A: np.ndarray,
        stress_eV_per_A3: np.ndarray, pressure_eV_per_A3: float,
    ) -> np.ndarray:
        """Physical d(E+PV)/dy in the shared reference-atom/active-cell metric."""
        # Reuse legacy array/identity validation and the same atomic derivative.
        # Its symmetric-cell components are discarded, never sampled as probes.
        legacy = super().enthalpy_gradient(
            atoms, forces_eV_per_A, stress_eV_per_A3, pressure_eV_per_A3,
        )
        self._validate_candidate(atoms)
        force = cell_force_from_arrays(
            np.asarray(stress_eV_per_A3, dtype=float), atoms, self.reference.cell.array,
            pressure=pressure_eV_per_A3,
        )
        gradient = -np.einsum("ij,aij->a", force, self.deformation_basis) / self.cell_scale_A
        return np.concatenate((legacy[:3 * len(self.reference)], gradient))

    def translation_free_basis(self) -> np.ndarray:
        """Orthonormal atomic--active-cell basis with three translations removed."""
        return super().translation_free_basis()
