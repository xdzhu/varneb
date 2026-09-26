"""Explicit optimizer subspace for a flat two-dimensional slab.

The VCNEB cell mask alone leaves ``F_xy`` and ``F_yx`` independent. For an
in-plane symmetric-strain comparison, that admits an antisymmetric rotation.
This helper pairs the mask with an orthonormal optimizer basis that includes
all atomic freedoms but only ``F_xx``, ``F_yy`` and symmetric ``F_xy=F_yx``.
It assumes a Cartesian z-directed vacuum vector and does not define a 2D
stress or an effective slab thickness.
"""

from __future__ import annotations

import numpy as np


def symmetric_inplane_vcneb_boundary(
    n_atoms: int, reference_cell: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(cell_mask, mode_basis)`` for ``run_vcneb``/``VCNEB``.

    Pass both arrays with ``constraint_mode='subspace'``. All atoms remain
    active. Only the x/y symmetric deformation components are active; the
    vacuum vector and xz/yz shear remain fixed. The reference cell may be
    oblique in the xy plane, but its vacuum vector must be parallel to z.
    """

    if isinstance(n_atoms, bool) or not isinstance(n_atoms, (int, np.integer)) or n_atoms < 1:
        raise ValueError("n_atoms must be a positive integer")
    cell = np.asarray(reference_cell, dtype=float)
    if cell.shape != (3, 3) or not np.all(np.isfinite(cell)):
        raise ValueError("reference_cell must be a finite 3x3 matrix")
    if np.linalg.det(cell) <= 0:
        raise ValueError("reference_cell must have positive volume")
    scale = float(np.max(np.abs(cell)))
    if not np.allclose(cell[:2, 2], 0.0, rtol=0.0, atol=1e-10 * scale) or not np.allclose(
        cell[2, :2], 0.0, rtol=0.0, atol=1e-10 * scale,
    ):
        raise ValueError("reference_cell must have an xy plane and a z-directed vacuum vector")

    mask = np.zeros((3, 3), dtype=float)
    mask[:2, :2] = 1.0
    atomic_count = 3 * int(n_atoms)
    basis = np.zeros((atomic_count + 9, atomic_count + 3), dtype=float)
    basis[:atomic_count, :atomic_count] = np.eye(atomic_count)
    basis[atomic_count + 0, atomic_count + 0] = 1.0  # F_xx
    basis[atomic_count + 4, atomic_count + 1] = 1.0  # F_yy
    basis[atomic_count + 1, atomic_count + 2] = 1.0 / np.sqrt(2.0)  # F_xy
    basis[atomic_count + 3, atomic_count + 2] = 1.0 / np.sqrt(2.0)  # F_yx
    mask.flags.writeable = False
    basis.flags.writeable = False
    return mask, basis


__all__ = ["symmetric_inplane_vcneb_boundary"]
