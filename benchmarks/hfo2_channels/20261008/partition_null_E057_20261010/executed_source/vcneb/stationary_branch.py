"""Fixed-control stationary branches of an explicitly restricted quadratic.

This is standard harmonic continuation, not a material saddle certificate,
branch search or barrier predictor. Stable omitted directions may be released;
the external control is never relaxed or counted in the internal Morse index.
Unsampled/frozen directions, trust-region coverage and anharmonicity remain
explicit caller responsibilities. No calculator is invoked.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real

import numpy as np

from .quadratic_reduction import ConditionalQuadraticModel, condition_quadratic_energy


def _owned(values):
    result = np.array(values, dtype=float, copy=True)
    result.setflags(write=False)
    return result


def _scalar(value, name):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
            or not np.isfinite(value)):
        raise ValueError(f"{name} must be a finite real scalar")
    return float(value)


@dataclass(frozen=True)
class StationaryQuadraticPoint:
    """Model stationary in its declared internal space, not the full crystal."""

    control_value: float
    retained_internal_coordinates: np.ndarray
    eliminated_coordinates: np.ndarray
    full_displacement: np.ndarray
    energy_change: float
    control_gradient: float
    declared_internal_gradient_norm: float
    unrepresented_internal_gradient_norm: float
    clamped_gradient_norm: float


@dataclass(frozen=True)
class StationaryQuadraticBranch:
    """Continuation with fixed index in Q plus the stable released R space.

    All coordinates share a declared metric and units. The reported eigenvalues
    are those of the condensed Q block, not phonon frequencies or eigenvalues
    of a full physical Hessian. No sampled trust region is inferred.
    """

    conditional: ConditionalQuadraticModel
    hessian: np.ndarray
    gradient: np.ndarray
    admissible_internal_basis: np.ndarray
    retained_internal_eigenvalues: np.ndarray
    internal_offset: np.ndarray
    internal_response: np.ndarray
    control_curvature: float
    expected_index: int

    def evaluate(self, control_value: float) -> StationaryQuadraticPoint:
        t = _scalar(control_value, "control value")
        q = self.internal_offset + self.internal_response * t
        retained = np.concatenate((q, [t]))
        r = self.conditional.eliminated_coordinates(retained)
        displacement = self.conditional.full_displacement(retained)
        physical = self.gradient + self.hessian @ displacement
        basis = self.conditional.retained_basis
        released = self.conditional.eliminated_basis
        internal_gradient = np.concatenate((basis[:, :-1].T @ physical, released.T @ physical))
        represented = basis[:, :-1] @ (basis[:, :-1].T @ physical) + released @ (released.T @ physical)
        admissible = self.admissible_internal_basis @ (self.admissible_internal_basis.T @ physical)
        unrepresented = admissible - represented
        clamped = physical - admissible - basis[:, -1] * (basis[:, -1] @ physical)
        energy = self.conditional.energy_change(retained)
        control_gradient = float(basis[:, -1] @ physical)
        norms = np.array([np.linalg.norm(internal_gradient), np.linalg.norm(unrepresented), np.linalg.norm(clamped)])
        if (not all(np.isfinite(v).all() for v in (q, r, displacement, physical, norms))
                or not np.isfinite(energy) or not np.isfinite(control_gradient)):
            raise ValueError("stationary quadratic evaluation is nonfinite; no extrapolation repair")
        return StationaryQuadraticPoint(
            t, _owned(q), _owned(r), _owned(displacement), float(energy),
            control_gradient, float(norms[0]), float(norms[1]), float(norms[2]),
        )


def stationary_quadratic_branch(
    hessian: np.ndarray, gradient: np.ndarray, reference_energy: float,
    retained_internal_basis: np.ndarray, eliminated_basis: np.ndarray,
    control_direction: np.ndarray, *, admissible_internal_basis: np.ndarray,
    expected_index: int,
    stability_floor: float, internal_curvature_floor: float,
    symmetry_tolerance: float = 1e-8,
) -> StationaryQuadraticBranch:
    """Continue a minimum (index0) or first-order saddle (index1) at fixed t.

    Q, R and control c are mutually orthonormal columns in one chart. The
    controlled coordinate must be supplied separately; it cannot be eliminated.
    Q and R must lie in the supplied actual mechanical internal space, which
    must be orthogonal to the control (e.g. the biaxial chart's translation-free
    internal basis). Mechanically forbidden reactions are distinct from
    unrepresented, physically movable internal gradients.
    The R block must exceed the measured stability floor. Every eigenvalue of
    the condensed Q block must be resolved above the internal floor in absolute
    value, with exactly the requested negative count. No pseudoinverse is used.

    Nonzero reference gradients are retained: q*(t)=-Kqq^-1(gq+Kqt*t).
    Eliminated offsets, reference relaxation energy and control work all follow
    the unchanged same-centre conditional model. Curvature along this moving
    stationary branch need not soften: Ktt-Ktq*Kqq^-1*Kqt includes negative Q.
    """
    floor = _scalar(internal_curvature_floor, "internal curvature floor")
    release_floor = _scalar(stability_floor, "stability floor")
    tolerance = _scalar(symmetry_tolerance, "symmetry tolerance")
    if floor < 0 or release_floor < 0 or tolerance <= 0:
        raise ValueError("curvature floors must be nonnegative and symmetry tolerance positive")
    if (isinstance(expected_index, (bool, np.bool_)) or not isinstance(expected_index, Integral)
            or expected_index not in (0, 1)):
        raise ValueError("expected index must be integer0 (minimum) or1 (first-order saddle)")
    if any(np.iscomplexobj(v) for v in (hessian, gradient, retained_internal_basis,
                                      eliminated_basis, control_direction, admissible_internal_basis)):
        raise ValueError("a declared real coordinate chart is required; do not discard complex phases")
    h = np.asarray(hessian, dtype=float)
    q = np.asarray(retained_internal_basis, dtype=float)
    c = np.asarray(control_direction, dtype=float)
    allowed = np.asarray(admissible_internal_basis, dtype=float)
    if (h.ndim != 2 or h.shape[0] != h.shape[1] or not h.size
            or q.ndim != 2 or q.shape[0] != h.shape[0] or q.shape[1] == 0
            or c.shape != (h.shape[0],)):
        raise ValueError("nonempty retained internal columns and one matching control vector required")
    r = np.asarray(eliminated_basis, dtype=float)
    if (allowed.ndim != 2 or allowed.shape[0] != h.shape[0] or allowed.shape[1] == 0
            or not np.isfinite(allowed).all()
            or not np.allclose(allowed.T @ allowed, np.eye(allowed.shape[1]), atol=1e-10, rtol=0)
            or not np.allclose(allowed.T @ c, 0., atol=1e-10, rtol=0)):
        raise ValueError("admissible internal basis must be finite, orthonormal and exclude the control")
    if r.ndim != 2 or r.shape[0] != h.shape[0]:
        raise ValueError("matching eliminated columns required")
    represented = np.column_stack((q, r))
    if not np.allclose(allowed @ (allowed.T @ represented), represented, atol=1e-10, rtol=0):
        raise ValueError("retained/released directions leave the admissible mechanical internal space")
    model = condition_quadratic_energy(
        h, gradient, _scalar(reference_energy, "reference energy"),
        np.column_stack((q, c)), eliminated_basis,
        stability_floor=release_floor, symmetry_tolerance=tolerance,
    )
    k = model.curvature.relaxed
    internal = k[:-1, :-1]
    eigenvalues = np.linalg.eigvalsh(internal)
    if not np.isfinite(eigenvalues).all() or np.any(np.abs(eigenvalues) <= floor):
        raise ValueError("retained internal curvature is singular or unresolved; resolve or retain another branch")
    if np.count_nonzero(eigenvalues < 0) != expected_index:
        raise ValueError("retained internal curvature does not have the expected stationary index")
    offset = -np.linalg.solve(internal, model.retained_gradient_at_zero[:-1])
    response = -np.linalg.solve(internal, k[:-1, -1])
    control_curvature = float(k[-1, -1] + k[-1, :-1] @ response)
    if (not np.isfinite(offset).all() or not np.isfinite(response).all()
            or not np.isfinite(control_curvature)):
        raise ValueError("stationary response is nonfinite; no pseudoinverse repair")
    return StationaryQuadraticBranch(
        model, _owned(0.5 * (h + h.T)), _owned(gradient), _owned(allowed), _owned(eigenvalues),
        _owned(offset), _owned(response), control_curvature, int(expected_index),
    )
