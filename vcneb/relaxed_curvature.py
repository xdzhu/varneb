"""Local curvature after *stable* orthogonal degrees of freedom relax.

The Schur complement is a standard harmonic result, not a new physical model.
All variables must use one declared coordinate/unit convention. This module
does not certify stationarity, identify phonons or call a calculator. Unstable
or unresolved eliminated directions are rejected instead of pseudo-inverted.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class RelaxedCurvature:
    frozen: np.ndarray
    relaxed: np.ndarray
    softening: np.ndarray
    orthogonal_response: np.ndarray
    minimum_eliminated_curvature: float | None
    eliminated_condition_number: float | None


def relax_orthogonal_curvature(
    hessian: np.ndarray,
    retained_basis: np.ndarray,
    eliminated_basis: np.ndarray,
    *,
    stability_floor: float,
    symmetry_tolerance: float = 1e-8,
) -> RelaxedCurvature:
    """Condense selected stable directions of a symmetric local Hessian.

    Bases are orthonormal columns in the *same* full coordinate space. They
    need not span the whole space: omitted directions are clamped. At fixed
    retained displacement q, dr/dq = -H_rr^-1 H_rq. The relaxed curvature is
    H_qq - H_qr H_rr^-1 H_rq. Retained directions may be unstable (e.g. a
    saddle), but every eliminated eigenvalue must exceed stability_floor,
    which the caller must choose from measured curvature resolution.
    """

    matrix = np.asarray(hessian, dtype=float)
    q = np.asarray(retained_basis, dtype=float)
    r = np.asarray(eliminated_basis, dtype=float)
    if (matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]
            or not matrix.size or not np.isfinite(matrix).all()):
        raise ValueError("hessian must be a nonempty finite square matrix")
    if (not np.isfinite(stability_floor) or stability_floor < 0
            or not np.isfinite(symmetry_tolerance) or symmetry_tolerance <= 0):
        raise ValueError("stability floor must be nonnegative and symmetry tolerance positive")
    n = matrix.shape[0]
    if (q.ndim != 2 or r.ndim != 2 or q.shape[0] != n or r.shape[0] != n
            or q.shape[1] == 0 or not np.isfinite(q).all() or not np.isfinite(r).all()):
        raise ValueError("bases must be finite columns of the Hessian coordinate space")
    combined = np.concatenate((q, r), axis=1)
    if not np.allclose(combined.T @ combined, np.eye(combined.shape[1]), rtol=0, atol=1e-10):
        raise ValueError("retained and eliminated bases must be mutually orthonormal")
    defect = np.linalg.norm(matrix - matrix.T)
    if defect > symmetry_tolerance * max(np.linalg.norm(matrix), 1e-30):
        raise ValueError("hessian reciprocity defect exceeds symmetry tolerance")
    matrix = (matrix + matrix.T) * 0.5
    frozen = q.T @ matrix @ q
    if r.shape[1] == 0:
        return RelaxedCurvature(
            frozen.copy(), frozen.copy(), np.zeros_like(frozen),
            np.empty((0, q.shape[1])), None, None,
        )
    hrr = r.T @ matrix @ r
    hrq = r.T @ matrix @ q
    eigenvalues = np.linalg.eigvalsh(hrr)
    minimum = float(eigenvalues[0])
    if minimum <= stability_floor:
        raise ValueError(
            "eliminated curvature is unstable or unresolved: "
            f"minimum={minimum:.8g}, floor={stability_floor:.8g}; "
            "retain the soft direction or resolve a separate branch"
        )
    response = -np.linalg.solve(hrr, hrq)
    softening = -(hrq.T @ response)
    relaxed = frozen - softening
    return RelaxedCurvature(
        frozen, (relaxed + relaxed.T) * 0.5, (softening + softening.T) * 0.5,
        response, minimum, float(eigenvalues[-1] / minimum),
    )
