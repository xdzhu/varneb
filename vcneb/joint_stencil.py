"""Paired joint-coordinate probes and measured Hessian actions.

All directions use the chart's declared Angstrom metric. Calculators and
physical input parameters stay outside this module. A projected Hessian is
not the missing full-space Hessian, and a local Hessian is not a TS certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np
from ase import Atoms

from .biaxial_curvature import BiaxialClampedCurvatureCoordinates
from .joint_curvature import JointCurvatureCoordinates, central_difference_hessian


def _basis(values: np.ndarray, size: int | None = None) -> np.ndarray:
    basis = np.array(values, dtype=float, copy=True)
    if (basis.ndim != 2 or min(basis.shape) == 0 or not np.isfinite(basis).all()
            or (size is not None and basis.shape[0] != size)
            or not np.allclose(basis.T @ basis, np.eye(basis.shape[1]), rtol=0, atol=1e-10)):
        raise ValueError("directions must be finite orthonormal columns of the chart space")
    return basis


def _step(value: float) -> float:
    if isinstance(value, (bool, np.bool_)) or not np.isfinite(value) or value <= 0:
        raise ValueError("step_A must be finite and positive")
    return float(value)


def _readonly(values: np.ndarray) -> np.ndarray:
    result = np.array(values, copy=True)
    result.flags.writeable = False
    return result


@dataclass(frozen=True)
class JointCurvatureProbe:
    """One inert point tagged by direction and sign; atoms remain ASE-editable."""

    direction_index: int
    sign: int
    step_A: float
    delta_A: np.ndarray
    atoms: Atoms


def joint_curvature_probes(
    coordinates: JointCurvatureCoordinates | BiaxialClampedCurvatureCoordinates,
    directions: np.ndarray, *, step_A: float,
    center_delta_A: np.ndarray | None = None,
    candidate_validator: Callable[[Sequence[Atoms]], None] | None = None,
) -> tuple[JointCurvatureProbe, ...]:
    """Generate exactly 2*k probes, ordered direction then (+,-), without SCF.

    Directions must already be orthonormal in the chart metric. No modes,
    translations, integer lifts or mechanical freedoms are silently adjusted.
    Optional geometric guards reject candidates; they must not modify them.
    The full stencil is built and validated before returning any point.
    """
    if not isinstance(coordinates, (JointCurvatureCoordinates, BiaxialClampedCurvatureCoordinates)):
        raise TypeError("coordinates must be a joint curvature chart")
    if coordinates.reference.constraints:
        raise ValueError("extra ASE constraints require an explicit atomic subspace")
    basis = _basis(directions, coordinates.size)
    h = _step(step_A)
    center = np.zeros(coordinates.size) if center_delta_A is None else np.array(center_delta_A, dtype=float)
    if center.shape != (coordinates.size,) or not np.isfinite(center).all():
        raise ValueError("center_delta_A must be a finite chart coordinate")
    if candidate_validator is not None and not callable(candidate_validator):
        raise TypeError("candidate_validator must be callable")
    # Check the declared center as well as its paired probes. No calculator is
    # inherited even when the caller's reference already contains a result.
    origin = coordinates.displaced(center)
    probes = []
    for index in range(basis.shape[1]):
        for sign in (1, -1):
            delta = center + sign * h * basis[:, index]
            probes.append(JointCurvatureProbe(index, sign, h, _readonly(delta), coordinates.displaced(delta)))
    points = (origin, *(probe.atoms for probe in probes))
    if any(not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell.array).all()
           or np.linalg.det(atoms.cell.array) <= 0 for atoms in points):
        raise ValueError("probe geometry must be finite with positive volume")
    if candidate_validator is not None:
        before = [(atoms.numbers.copy(), atoms.pbc.copy(), atoms.positions.copy(), atoms.cell.array.copy())
                  for atoms in points]
        candidate_validator(points)
        if any(atoms.calc is not None or atoms.constraints
               or any(not np.array_equal(actual, expected) for actual, expected in
                      zip((atoms.numbers, atoms.pbc, atoms.positions, atoms.cell.array), original))
               for atoms, original in zip(points, before)):
            raise ValueError("candidate_validator must reject only, not modify or attach calculators")
    return tuple(probes)


@dataclass(frozen=True)
class DirectionalJointCurvature:
    """Measured H*B, B.T*H*B and response outside B, all in eV/Angstrom².

    Transverse response is coupling to unsampled directions, not necessarily
    a numerical error. It does not determine their self-curvature or justify
    eliminating them. The raw projected matrix is retained before symmetry.
    """

    directions: np.ndarray
    step_A: float
    full_hessian_action: np.ndarray
    raw_projected: np.ndarray
    symmetric_projected: np.ndarray
    transverse_action: np.ndarray
    reciprocity_relative_defect: float
    symmetrization_operator_change_eV_A2: float


def assemble_joint_directional_curvature(
    directions: np.ndarray, plus_gradients_eV_A: np.ndarray,
    minus_gradients_eV_A: np.ndarray, *, step_A: float,
) -> DirectionalJointCurvature:
    """Assemble from complete full-chart gradients at center +/- h*B_j.

    Each gradient array has shape (k directions, d chart coordinates), not
    (k,k). Use the chart's physical enthalpy gradient, never projected NEB or
    spring forces. The caller must verify paired geometries, pressure, common
    center/basis and the unchanged calculator contract from raw provenance.
    No acceptance threshold, error bar or stability certificate is invented.
    """
    basis = _basis(directions)
    h = _step(step_A)
    plus = np.asarray(plus_gradients_eV_A, dtype=float)
    minus = np.asarray(minus_gradients_eV_A, dtype=float)
    expected = (basis.shape[1], basis.shape[0])
    if (plus.shape != expected or minus.shape != expected
            or not np.isfinite(plus).all() or not np.isfinite(minus).all()):
        raise ValueError("paired gradients must have finite shape (directions, full chart coordinates)")
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            action = (plus - minus).T / (2 * h)
            raw = basis.T @ action
            symmetric, defect = central_difference_hessian(plus @ basis, minus @ basis, h)
            transverse = action - basis @ raw
            change = float(np.linalg.norm(raw - symmetric, ord=2))
    except FloatingPointError as exc:
        raise ValueError("curvature difference is not numerically representable at step_A") from exc
    if (not all(np.isfinite(array).all() for array in (action, raw, symmetric, transverse))
            or not np.isfinite(defect) or not np.isfinite(change)):
        raise ValueError("curvature difference is not numerically representable at step_A")
    return DirectionalJointCurvature(
        _readonly(basis), h, _readonly(action), _readonly(raw), _readonly(symmetric),
        _readonly(transverse), defect, change,
    )
