"""Stable quadratic minimization with shared linear coordinate constraints.

Standard equality-constrained quadratic algebra, not a new physical model.
Omitted directions remain clamped. Positive curvature is checked only in
the allowed constraint-null space; the reference need not be stationary.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def _real(values):
    if np.iscomplexobj(values):
        raise ValueError('real coordinates required; no imaginary-part removal')
    return np.array(values, dtype=float, copy=True)


def _readonly(values):
    result = np.array(values, dtype=float, copy=True)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class AffineQuadraticModel:
    """Minimize E0+g.x+x.H.x/2 over x in span(U), with A.T@x=q.

    Only allowed-space stationarity is checked. No DFT stability, trust
    domain, minimum-distance or whole-path prediction is certified here.
    """
    hessian: np.ndarray
    gradient: np.ndarray
    allowed_basis: np.ndarray
    observables: np.ndarray
    offset: np.ndarray
    response: np.ndarray
    multiplier_offset: np.ndarray
    multiplier_response: np.ndarray
    minimum_free_curvature: float | None
    constraint_singular_values: np.ndarray

    def _q(self, q):
        values = _real(q)
        if values.shape != (self.observables.shape[1],) or not np.isfinite(values).all():
            raise ValueError('finite coordinates of the declared constraint dimension required')
        return values

    def coordinates(self, q):
        values = self._q(q)
        return self.offset + self.response @ values

    def energy_change(self, q):
        x = self.coordinates(q)
        return float(self.gradient @ x + .5 * x @ self.hessian @ x)

    def retained_gradient(self, q):
        values = self._q(q)
        return -(self.multiplier_offset + self.multiplier_response @ values)

    def residuals(self, q):
        values = self._q(q)
        x = self.coordinates(values)
        multipliers = -self.retained_gradient(values)
        c = self.observables.T @ self.allowed_basis
        return dict(constraint_norm=float(np.linalg.norm(self.observables.T @ x - values)),
            allowed_stationarity_norm=float(np.linalg.norm(
                self.allowed_basis.T @ (self.gradient + self.hessian @ x) + c.T @ multipliers)))


def constrain_quadratic(hessian, gradient, allowed_basis, observables, *, stability_floor,
                        rank_tolerance=1e-10, symmetry_tolerance=1e-8):
    """Solve the KKT system without clamping/clipping unstable free directions.

    All arrays share one caller-declared physical coordinate metric and
    energy unit. U is orthonormal; A need not be, and is not renormalized.
    A.T@U must have full row rank. Its null-space Hessian must lie strictly
    above the explicit floor. A zero floor is only a formal raw-sign test,
    not measured physical resolution. Empty free space is not a stability
    measurement. No pseudoinverse, fallback branch or calculator is used.
    """
    h, g, u, a = map(_real, (hessian, gradient, allowed_basis, observables))
    if (h.ndim != 2 or h.shape[0] != h.shape[1] or not h.size
            or g.shape != (h.shape[0],) or u.ndim != 2 or a.ndim != 2
            or u.shape[0] != h.shape[0] or a.shape[0] != h.shape[0]
            or not 1 <= a.shape[1] <= u.shape[1] <= h.shape[0]
            or not all(np.isfinite(v).all() for v in (h, g, u, a))):
        raise ValueError('finite compatible Hessian/gradient/allowed/constraint arrays required')
    if (not np.isfinite(stability_floor) or stability_floor < 0
            or not np.isfinite(rank_tolerance) or rank_tolerance <= 0
            or not np.isfinite(symmetry_tolerance) or symmetry_tolerance <= 0):
        raise ValueError('nonnegative stability floor and positive numerical tolerances required')
    if not np.allclose(u.T @ u, np.eye(u.shape[1]), rtol=0, atol=1e-10):
        raise ValueError('orthonormal allowed basis required; do not silently change its metric')
    if np.linalg.norm(h - h.T) > symmetry_tolerance * max(np.linalg.norm(h), 1e-30):
        raise ValueError('Hessian reciprocity defect exceeds declared tolerance')
    h = (h + h.T) * .5
    c = a.T @ u
    _, singular, vh = np.linalg.svd(c, full_matrices=True)
    if singular[-1] <= rank_tolerance * max(float(singular[0]), 1.):
        raise ValueError('allowed reference cannot resolve all declared shared coordinates')
    free = vh[a.shape[1]:].T
    hu = u.T @ h @ u
    minimum = None
    if free.shape[1]:
        minimum = float(np.linalg.eigvalsh(free.T @ hu @ free)[0])
        if minimum <= stability_floor:
            raise ValueError('constraint-null curvature is unstable or unresolved; do not invert it')
    k, m = u.shape[1], a.shape[1]
    kkt = np.block([[hu, c.T], [c, np.zeros((m, m))]])
    rhs = np.zeros((k + m, 1 + m))
    rhs[:k, 0] = -(u.T @ g)
    rhs[k:, 1:] = np.eye(m)
    solution = np.linalg.solve(kkt, rhs)
    if not np.isfinite(solution).all():
        raise ValueError('nonfinite constrained quadratic solution')
    return AffineQuadraticModel(hessian=_readonly(h), gradient=_readonly(g),
        allowed_basis=_readonly(u), observables=_readonly(a),
        offset=_readonly(u @ solution[:k, 0]), response=_readonly(u @ solution[:k, 1:]),
        multiplier_offset=_readonly(solution[k:, 0]), multiplier_response=_readonly(solution[k:, 1:]),
        minimum_free_curvature=minimum, constraint_singular_values=_readonly(singular))
