"""A fixed substrate plane, not a fixed vacuum or a stress-only mask.

ASE cell vectors are rows: ``H = H0 @ F.T``. Holding the first two
rows of H fixed permits exactly ``delta F = v outer n``, where n is
their unit normal. v is unrestricted when out-of-plane tilt is allowed;
otherwise it must be parallel to n. This affine subspace is exact, even
for an oblique or rotated reference cell. No global cell rotation can
hold two non-collinear substrate vectors fixed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from ase import Atoms
from ase.filters import UnitCellFilter


def _cell(value: np.ndarray, name: str) -> np.ndarray:
    result = np.array(value, dtype=float, copy=True)
    if result.shape != (3, 3) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite 3x3 matrix")
    if np.linalg.det(result) <= 1e-12 or np.linalg.cond(result) > 1e12:
        raise ValueError(f"{name} must have positive volume and be nonsingular")
    return result


@dataclass(frozen=True)
class ClampedPlaneBoundary:
    """Validated optimizer subspace for one prescribed substrate.

    Obtain this object with :func:`clamped_plane_vcneb_boundary`, then use
    ``VCNEB(images, **boundary.vcneb_kwargs(images))``. The pre-construction
    check is important: VCNEB's general subspace interface can intentionally
    project interior seeds. This boundary rejects incompatible seeds instead.
    Arrays are read-only, atomic freedoms are unrestricted, and the all-one
    deformation mask must be paired with the supplied subspace.
    """

    n_atoms: int
    reference_cell: np.ndarray
    normal: np.ndarray
    allow_tilt: bool
    cell_mask: np.ndarray
    mode_basis: np.ndarray

    @property
    def cell_dofs(self) -> int:
        return 3 if self.allow_tilt else 1

    def validate_images(self, images: Sequence[Atoms]) -> None:
        """Check substrate and tilt constraints without modifying any image."""
        if not images:
            raise ValueError("clamped-plane boundary needs at least one image")
        symbols = images[0].get_chemical_symbols()
        for index, image in enumerate(images):
            if len(image) != self.n_atoms or image.get_chemical_symbols() != symbols:
                raise ValueError(f"image {index}: inconsistent atom count or order")
            cell = _cell(image.cell.array, f"image {index} cell")
            if not np.allclose(cell[:2], self.reference_cell[:2], rtol=0.0, atol=1e-9):
                raise ValueError(f"image {index}: substrate vectors differ from the prescribed plane")
            if not self.allow_tilt:
                delta = cell[2] - self.reference_cell[2]
                in_plane = delta - np.dot(delta, self.normal) * self.normal
                if np.linalg.norm(in_plane) > 1e-9:
                    raise ValueError(f"image {index}: out-of-plane tilt is not allowed")

    def vcneb_kwargs(self, images: Sequence[Atoms]) -> dict:
        """Validate the *unmodified* chain before returning VCNEB keywords.

        Additional geometry validators can be composed explicitly with
        ``validate_images``; do not silently replace this candidate check.
        The same fixed plane must be used when relaxing the endpoints.
        """
        self.validate_images(images)
        return {
            "cell_mask": self.cell_mask,
            "mode_basis": self.mode_basis,
            "constraint_mode": "subspace",
            "candidate_validator": self.validate_images,
        }


def clamped_plane_vcneb_boundary(
    n_atoms: int, reference_cell: np.ndarray, *, allow_tilt: bool,
) -> ClampedPlaneBoundary:
    """Fix ASE cell rows 0/1; release row 2's length and optional two tilts.

    ``allow_tilt`` is required because normal-only and shear-released
    epitaxy are different mechanical ensembles. The mask is all ones:
    elementwise masking cannot describe a general rotated substrate. The
    orthonormal deformation basis enforces the actual constraint.
    """
    if isinstance(n_atoms, bool) or not isinstance(n_atoms, (int, np.integer)) or n_atoms < 1:
        raise ValueError("n_atoms must be a positive integer")
    if not isinstance(allow_tilt, bool):
        raise ValueError("allow_tilt must be an explicit boolean")
    cell = _cell(reference_cell, "reference_cell")
    normal = np.cross(cell[0], cell[1])
    normal /= np.linalg.norm(normal)
    directions = np.eye(3) if allow_tilt else normal.reshape(1, 3)
    atomic_count = 3 * int(n_atoms)
    basis = np.zeros((atomic_count + 9, atomic_count + len(directions)))
    basis[:atomic_count, :atomic_count] = np.eye(atomic_count)
    for index, direction in enumerate(directions):
        basis[atomic_count:, atomic_count + index] = np.outer(direction, normal).reshape(-1)
    mask = np.ones((3, 3))
    for array in (cell, normal, basis, mask):
        array.flags.writeable = False
    return ClampedPlaneBoundary(int(n_atoms), cell, normal, allow_tilt, mask, basis)


def cell_work_derivative(
    stress: np.ndarray, current_cell: np.ndarray, cell_derivative: np.ndarray,
    *, pressure: float = 0.0,
) -> float:
    """Return d(E+PV)/d(parameter) for a prescribed row-cell derivative.

    Stress is ASE tensile-positive in eV/A^3, cell and cell_derivative in
    A and A per parameter respectively, and pressure is positive for
    compression in eV/A^3. Fractional positions are held fixed in this
    partial derivative. It is an envelope derivative only on a stationary
    branch with consistent open-coordinate constraints; a NEB peak image
    or a frozen slice does not establish that condition.
    """
    cell = _cell(current_cell, "current_cell")
    sigma = np.asarray(stress, dtype=float)
    derivative = np.asarray(cell_derivative, dtype=float)
    if sigma.shape != (3, 3) or not np.all(np.isfinite(sigma)):
        raise ValueError("stress must be a finite 3x3 matrix")
    if not np.allclose(sigma, sigma.T, atol=1e-12, rtol=1e-10):
        raise ValueError("stress must be symmetric in the ASE convention")
    if derivative.shape != (3, 3) or not np.all(np.isfinite(derivative)):
        raise ValueError("cell_derivative must be a finite 3x3 matrix")
    if not np.isfinite(pressure):
        raise ValueError("pressure must be finite")
    velocity_gradient = np.linalg.solve(cell, derivative).T
    return float(np.linalg.det(cell) * np.sum((sigma + pressure * np.eye(3)) * velocity_gradient))


class ClampedPlaneFilter(UnitCellFilter):
    """ASE BFGS-compatible endpoint relaxation in the same NEB ensemble.

    All atoms and only the boundary's open cell directions are optimized.
    Raw reaction stresses on the substrate are not a convergence failure.
    Read ``atoms.get_stress()`` for the physical stress, rather than treating
    the filter's generalized virial as an independently relaxed stress tensor.
    Extra ASE atom/cell constraints are refused: this interface promises all
    atomic freedoms and cannot silently let another constraint alter the plane.
    """

    def __init__(
        self, atoms: Atoms, boundary: ClampedPlaneBoundary, *,
        cell_scale_A: float, pressure_eV_per_A3: float = 0.0,
    ) -> None:
        if not isinstance(boundary, ClampedPlaneBoundary):
            raise TypeError("boundary must be a ClampedPlaneBoundary")
        if atoms.constraints:
            raise ValueError("extra ASE constraints must be removed or composed explicitly")
        if not np.isfinite(cell_scale_A) or cell_scale_A <= 0:
            raise ValueError("cell_scale_A must be finite and positive")
        if not np.isfinite(pressure_eV_per_A3):
            raise ValueError("pressure_eV_per_A3 must be finite")
        boundary.validate_images([atoms])
        self.boundary = boundary
        count = 3 * len(atoms)
        self._cell_basis = boundary.mode_basis[count:, count:]
        super().__init__(atoms, mask=np.ones((3, 3)), cell_factor=cell_scale_A,
                         scalar_pressure=pressure_eV_per_A3)

    def set_positions(self, new: np.ndarray, **kwargs) -> None:
        candidate = np.array(new, dtype=float, copy=True)
        if candidate.shape != (len(self.atoms) + 3, 3) or not np.all(np.isfinite(candidate)):
            raise ValueError("filter positions must have the expected finite shape")
        delta = (candidate[-3:] / self.cell_factor - np.eye(3)).reshape(-1)
        projected = self._cell_basis @ (self._cell_basis.T @ delta)
        deformation = np.eye(3) + projected.reshape(3, 3)
        candidate[-3:] = self.cell_factor * deformation
        preview = self.atoms.copy()
        preview.set_cell(self.orig_cell @ deformation.T, scale_atoms=False)
        self.boundary.validate_images([preview])  # reject before mutating atoms
        super().set_positions(candidate, **kwargs)

    def get_forces(self, **kwargs) -> np.ndarray:
        forces = super().get_forces(**kwargs)
        cell_force = forces[-3:].reshape(-1)
        forces[-3:] = (self._cell_basis @ (self._cell_basis.T @ cell_force)).reshape(3, 3)
        return forces


__all__ = ["ClampedPlaneBoundary", "ClampedPlaneFilter", "clamped_plane_vcneb_boundary",
           "cell_work_derivative"]
