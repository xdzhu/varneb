"""Evidence gate for interpolating a *measured* conditional mode surface.

Passing this screen permits an explicitly labelled exploratory interpolation.
It never certifies the global lower envelope, a transition state, or a barrier.
No calculator is invoked here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .mode_surface import (
    ConditionalModeSurface,
    ModePlane,
    audit_conditional_surface_continuity,
)


@dataclass(frozen=True)
class ConditionalPointEvidence:
    """Independent evidence for a sampled, fixed-Q conditional point."""

    calculator_id: str
    energy_reference_id: str
    raw_audit_sha256: str
    branch_audit_sha256: str
    curvature_audit_sha256: str
    curvature_resolution_sha256: str
    minimum_orthogonal_curvature: float | None
    curvature_uncertainty: float | None


@dataclass(frozen=True)
class ConditionalHoldout:
    """An independently calculated interior point of one grid cell."""

    cell_i: int
    cell_j: int
    q: tuple[float, float]
    energy: float
    coordinates: np.ndarray
    evidence: ConditionalPointEvidence
    orthogonal_gradient_norm: float


@dataclass(frozen=True)
class ConditionalInterpolationScreen:
    """A conservative screen, not a conditional-PES certificate."""

    ready_for_exploratory_interpolation: bool
    blockers: tuple[str, ...]
    maximum_holdout_energy_error: float | None
    maximum_holdout_coordinate_error: float | None
    maximum_neighbor_orthogonal_jump: float
    n_cells: int
    n_independent_holdouts: int


def screen_conditional_interpolation(
    plane: ModePlane,
    surface: ConditionalModeSurface,
    point_evidence: Sequence[Sequence[ConditionalPointEvidence]],
    holdouts: Sequence[ConditionalHoldout],
    *,
    energy_zero_id: str,
    gradient_tolerance: float,
    maximum_neighbor_orthogonal_jump: float,
    maximum_holdout_energy_error: float,
    maximum_holdout_coordinate_error: float,
) -> ConditionalInterpolationScreen:
    """Require raw audits, resolved local curvature, continuity and cell holdouts.

    Every cell needs an independent interior holdout. The cell prediction is
    bilinear in the *measured* four corners; the screen does not fit a PES.
    Curvature must exceed its independently estimated uncertainty. A pass does
    not rule out unsampled lower branches or guarantee contour fidelity away
    from the holdout points.
    """

    if not energy_zero_id.strip():
        raise ValueError("an explicit common energy-zero identifier is required")
    tolerances = (
        gradient_tolerance, maximum_neighbor_orthogonal_jump,
        maximum_holdout_energy_error, maximum_holdout_coordinate_error,
    )
    if any(not np.isfinite(value) or value <= 0 for value in tolerances):
        raise ValueError("all screening tolerances must be finite and positive")
    shape = surface.energies.shape
    if len(point_evidence) != shape[0] or any(len(row) != shape[1] for row in point_evidence):
        raise ValueError("point evidence must match the measured energy grid")
    continuity = audit_conditional_surface_continuity(
        plane, surface, maximum_allowed_orthogonal_jump=maximum_neighbor_orthogonal_jump,
    )
    blockers: list[str] = []
    calculator_ids = set()

    def check_evidence(evidence: ConditionalPointEvidence, label: str) -> None:
        if not isinstance(evidence, ConditionalPointEvidence):
            raise TypeError(f"{label}: expected ConditionalPointEvidence")
        if not evidence.calculator_id or evidence.energy_reference_id != energy_zero_id:
            blockers.append(f"{label}: incompatible calculator or energy zero")
        calculator_ids.add(evidence.calculator_id)
        digests = (
            evidence.raw_audit_sha256, evidence.branch_audit_sha256,
            evidence.curvature_audit_sha256, evidence.curvature_resolution_sha256,
        )
        if any(len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest)
               for digest in digests):
            blockers.append(f"{label}: missing raw/branch/curvature/resolution audit hashes")
        curvature = evidence.minimum_orthogonal_curvature
        uncertainty = evidence.curvature_uncertainty
        if (curvature is None or uncertainty is None
                or not np.isfinite(curvature) or not np.isfinite(uncertainty)
                or uncertainty < 0 or curvature <= uncertainty):
            blockers.append(f"{label}: positive orthogonal curvature is unresolved")

    for j, row in enumerate(point_evidence):
        for i, evidence in enumerate(row):
            check_evidence(evidence, f"grid({i},{j})")
    if np.any(surface.orthogonal_gradient_norms > gradient_tolerance):
        blockers.append("one or more grid points fail the orthogonal-gradient threshold")
    if not continuity.within_bound:
        blockers.append("an off-plane branch jump exceeds the neighbor bound")

    covered_cells: set[tuple[int, int]] = set()
    energy_errors: list[float] = []
    coordinate_errors: list[float] = []
    for index, holdout in enumerate(holdouts):
        label = f"holdout({index})"
        i, j = holdout.cell_i, holdout.cell_j
        if not 0 <= i < shape[1] - 1 or not 0 <= j < shape[0] - 1:
            raise ValueError(f"{label}: cell index is outside the measured grid")
        q = np.asarray(holdout.q, dtype=float)
        coordinates = np.asarray(holdout.coordinates, dtype=float)
        if (q.shape != (2,) or not np.all(np.isfinite(q))
                or coordinates.shape != (plane.reference.size,)
                or not np.all(np.isfinite(coordinates))
                or not np.isfinite(holdout.energy)
                or not np.isfinite(holdout.orthogonal_gradient_norm)):
            raise ValueError(f"{label}: nonfinite or malformed measured holdout")
        x0, x1 = surface.q1[i:i + 2]
        y0, y1 = surface.q2[j:j + 2]
        if not x0 < q[0] < x1 or not y0 < q[1] < y1:
            raise ValueError(f"{label}: holdout must lie strictly inside its cell")
        if not np.allclose(plane.project(coordinates), q, rtol=0, atol=1e-8):
            blockers.append(f"{label}: relaxed coordinates drifted from fixed Q")
        check_evidence(holdout.evidence, label)
        if holdout.orthogonal_gradient_norm < 0 or holdout.orthogonal_gradient_norm > gradient_tolerance:
            blockers.append(f"{label}: orthogonal gradient has not converged")
        tx, ty = (q[0] - x0) / (x1 - x0), (q[1] - y0) / (y1 - y0)
        weights = np.array([[(1 - tx) * (1 - ty), tx * (1 - ty)],
                            [(1 - tx) * ty, tx * ty]])
        prediction = float(np.sum(weights * surface.energies[j:j + 2, i:i + 2]))
        predicted_coordinates = np.einsum(
            "ji,jik->k", weights, surface.coordinates[j:j + 2, i:i + 2],
        )
        energy_errors.append(abs(holdout.energy - prediction))
        coordinate_errors.append(float(np.linalg.norm(
            np.sqrt(plane.metric_weights) * (coordinates - predicted_coordinates),
        )))
        covered_cells.add((i, j))
    n_cells = (shape[0] - 1) * (shape[1] - 1)
    if len(covered_cells) != n_cells:
        blockers.append("not every grid cell has an independent interior holdout")
    if len(calculator_ids) != 1:
        blockers.append("grid and holdout points use inconsistent calculator contracts")
    if energy_errors and max(energy_errors) > maximum_holdout_energy_error:
        blockers.append("holdout energy interpolation error exceeds its threshold")
    if coordinate_errors and max(coordinate_errors) > maximum_holdout_coordinate_error:
        blockers.append("holdout branch-coordinate error exceeds its threshold")
    return ConditionalInterpolationScreen(
        ready_for_exploratory_interpolation=not blockers,
        blockers=tuple(blockers),
        maximum_holdout_energy_error=max(energy_errors) if energy_errors else None,
        maximum_holdout_coordinate_error=max(coordinate_errors) if coordinate_errors else None,
        maximum_neighbor_orthogonal_jump=continuity.maximum_orthogonal_jump,
        n_cells=n_cells, n_independent_holdouts=len(holdouts),
    )


__all__ = [
    "ConditionalPointEvidence", "ConditionalHoldout", "ConditionalInterpolationScreen",
    "screen_conditional_interpolation",
]
