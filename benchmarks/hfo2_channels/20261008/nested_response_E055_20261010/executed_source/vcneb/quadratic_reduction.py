"""Restricted harmonic release about a possibly nonstationary reference.

This is a local quadratic model, not a conditional DFT surface or a new
Schur-complement theory. Coordinates, gradients and Hessian must share one
metric, reference, Hamiltonian and mechanical boundary. Directions absent
from both bases stay clamped. The caller must check sampling coverage and
anharmonicity before interpreting any model minimum as physical evidence.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .relaxed_curvature import RelaxedCurvature, relax_orthogonal_curvature


def _readonly(values: np.ndarray) -> np.ndarray:
    result = np.array(values, dtype=float, copy=True)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class ConditionalQuadraticModel:
    """Analytic minimum in the selected stable eliminated coordinates only.

    For E=E0+g.y+y.H.y/2 and y=Qq+Rr, r*(q)=offset+response.q.
    Energy changes are relative to E0, not to the released reference energy.
    Evaluation deliberately does not certify a trust region or full stability.
    """

    reference_energy: float
    reference_energy_change: float
    retained_gradient_at_zero: np.ndarray
    eliminated_offset: np.ndarray
    retained_basis: np.ndarray
    eliminated_basis: np.ndarray
    curvature: RelaxedCurvature

    def _coordinates(self, q: np.ndarray) -> np.ndarray:
        values = np.asarray(q, dtype=float)
        if (values.shape != (self.retained_basis.shape[1],)
                or not np.isfinite(values).all()):
            raise ValueError("retained coordinates must be a finite vector of the declared size")
        return values

    def eliminated_coordinates(self, q: np.ndarray) -> np.ndarray:
        values = self._coordinates(q)
        return self.eliminated_offset + self.curvature.orthogonal_response @ values

    def full_displacement(self, q: np.ndarray) -> np.ndarray:
        values = self._coordinates(q)
        return self.retained_basis @ values + self.eliminated_basis @ self.eliminated_coordinates(values)

    def energy_change(self, q: np.ndarray) -> float:
        """Return E_model(q,r*(q))-E0 without subtracting large total energies."""
        values = self._coordinates(q)
        return float(self.reference_energy_change + self.retained_gradient_at_zero @ values
                     + 0.5 * values @ self.curvature.relaxed @ values)

    def retained_gradient(self, q: np.ndarray) -> np.ndarray:
        values = self._coordinates(q)
        return self.retained_gradient_at_zero + self.curvature.relaxed @ values


def condition_quadratic_energy(
    hessian: np.ndarray,
    gradient: np.ndarray,
    reference_energy: float,
    retained_basis: np.ndarray,
    eliminated_basis: np.ndarray,
    *,
    stability_floor: float,
    symmetry_tolerance: float = 1e-8,
) -> ConditionalQuadraticModel:
    """Include g_r-induced offset and energy lowering in stable condensation.

    Input units are consistent caller-declared energy and coordinate units,
    e.g. eV, Angstrom, eV/Angstrom and eV/Angstrom**2. No mass conversion or
    silent coordinate rescaling is performed. Retained curvature need not
    be positive. The eliminated block must pass the existing strict stable
    curvature gate; an unstable/unresolved direction is never pseudoinverted.
    """
    curvature = relax_orthogonal_curvature(
        hessian, retained_basis, eliminated_basis, stability_floor=stability_floor,
        symmetry_tolerance=symmetry_tolerance,
    )
    matrix = np.asarray(hessian, dtype=float)
    matrix = 0.5 * (matrix + matrix.T)
    g = np.asarray(gradient, dtype=float)
    q = np.asarray(retained_basis, dtype=float)
    r = np.asarray(eliminated_basis, dtype=float)
    if (g.shape != (matrix.shape[0],) or not np.isfinite(g).all()
            or not np.isfinite(reference_energy)):
        raise ValueError("reference energy and matching full gradient must be finite")
    gr = r.T @ g
    offset = -np.linalg.solve(r.T @ matrix @ r, gr) if r.shape[1] else np.empty(0)
    effective_gradient = q.T @ g + q.T @ matrix @ r @ offset
    immutable_curvature = RelaxedCurvature(
        _readonly(curvature.frozen), _readonly(curvature.relaxed),
        _readonly(curvature.softening), _readonly(curvature.orthogonal_response),
        curvature.minimum_eliminated_curvature, curvature.eliminated_condition_number,
    )
    return ConditionalQuadraticModel(
        float(reference_energy), float(0.5 * gr @ offset), _readonly(effective_gradient),
        _readonly(offset), _readonly(q), _readonly(r), immutable_curvature,
    )
