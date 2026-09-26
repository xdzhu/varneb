"""Auditable two-order-parameter coordinates and frozen energy surfaces.

This module does not label arbitrary path images as phonons and never launches
an electronic-structure job by itself.  The caller supplies a fixed reference,
an explicitly normalized mode basis, and (for a surface) an energy evaluator.
Cell degrees of freedom can be included only when the caller supplies a
documented metric and coordinate convention for them.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

from .phonons import GammaModes


Array = np.ndarray


@dataclass(frozen=True)
class ModePlane:
    """Two signed, orthonormal order parameters made from one or more modes.

    ``reference`` and each column of ``basis`` use the same physical
    coordinate convention. ``metric_weights`` are the positive diagonal of
    the mode-projection metric: ``basis.T @ diag(metric_weights) @ basis = I``.
    ``axis_weights`` has two unit, mutually orthogonal columns in mode space.
    Thus either axis may represent one mode or a fixed signed combination.

    For Phonopy Gamma eigenvectors, use :meth:`from_gamma_modes`, which takes
    Cartesian atomic displacements in Angstrom and a mass metric.  A norm of a
    degenerate group is *not* a signed axis and is deliberately not accepted as
    an ``axis_weights`` column.
    """

    reference: Array
    basis: Array
    metric_weights: Array
    axis_weights: Array
    axis_labels: tuple[str, str]
    amplitude_unit: str
    reference_id: str

    def __post_init__(self) -> None:
        reference = np.asarray(self.reference, dtype=float).reshape(-1)
        basis = np.asarray(self.basis, dtype=float)
        metric = np.asarray(self.metric_weights, dtype=float).reshape(-1)
        weights = np.asarray(self.axis_weights, dtype=float)
        if basis.ndim != 2 or basis.shape[0] != reference.size or basis.shape[1] < 2:
            raise ValueError("basis must have shape (coordinate_count, mode_count >= 2)")
        if metric.shape != reference.shape or np.any(metric <= 0.0):
            raise ValueError("metric_weights must contain one positive value per coordinate")
        if weights.shape != (basis.shape[1], 2):
            raise ValueError("axis_weights must have shape (mode_count, 2)")
        if not all(np.all(np.isfinite(item)) for item in (reference, basis, metric, weights)):
            raise ValueError("mode-plane arrays must be finite")
        if len(self.axis_labels) != 2 or any(not label.strip() for label in self.axis_labels):
            raise ValueError("two nonempty axis labels are required")
        if not self.amplitude_unit.strip() or not self.reference_id.strip():
            raise ValueError("amplitude_unit and reference_id must be nonempty")
        gram = basis.T @ (metric[:, None] * basis)
        if not np.allclose(gram, np.eye(basis.shape[1]), rtol=1e-7, atol=1e-7):
            raise ValueError("basis is not orthonormal in the declared metric")
        if not np.allclose(weights.T @ weights, np.eye(2), rtol=1e-10, atol=1e-10):
            raise ValueError("axis combinations must be normalized and mutually orthogonal")
        for name, value in (("reference", reference), ("basis", basis),
                            ("metric_weights", metric), ("axis_weights", weights)):
            frozen = value.copy()
            frozen.flags.writeable = False
            object.__setattr__(self, name, frozen)

    @classmethod
    def from_gamma_modes(
        cls,
        modes: GammaModes,
        axis_weights: Array,
        *,
        axis_labels: tuple[str, str],
        reference_id: str,
    ) -> "ModePlane":
        """Construct atomic-only coordinates from endpoint Gamma eigenvectors.

        Input structures must first be aligned and represented as Cartesian
        displacements in the reference cell.  Cell strain is intentionally not
        projected onto these atomic-only axes.
        """

        return cls(
            reference=np.zeros(3 * modes.n_atoms),
            basis=modes.cartesian_eigenvectors(),
            metric_weights=np.repeat(modes.masses_amu, 3),
            axis_weights=axis_weights,
            axis_labels=axis_labels,
            amplitude_unit="sqrt(amu)*angstrom",
            reference_id=reference_id,
        )

    @classmethod
    def from_gamma_modes_with_strain(
        cls,
        modes: GammaModes,
        axis_weights: Array,
        *,
        strain_metric_weights_amu_A2: Array,
        axis_labels: tuple[str, str],
        reference_id: str,
    ) -> "ModePlane":
        """Append six symmetric-strain freedoms to atomic Gamma modes.

        Coordinate order is atomic reference-cell displacements followed by
        ASE Voigt strain ``(xx,yy,zz,yz,xz,xy)``. The six positive metric
        weights must be supplied in ``amu*angstrom^2``; there is no universal
        physical choice of atom/strain relative metric. The two selected axes
        remain purely atomic, and the added strain directions belong to the
        orthogonal complement available to conditional relaxation.
        """

        strain_weights = np.asarray(strain_metric_weights_amu_A2, dtype=float).reshape(-1)
        if strain_weights.shape != (6,) or not np.all(np.isfinite(strain_weights)) or np.any(strain_weights <= 0.0):
            raise ValueError("six positive finite strain metric weights in amu*angstrom^2 are required")
        atomic_basis = modes.cartesian_eigenvectors()
        n_atomic, n_modes = atomic_basis.shape
        weights = np.asarray(axis_weights, dtype=float)
        if weights.shape != (n_modes, 2):
            raise ValueError("atomic axis_weights must have shape (Gamma mode_count, 2)")
        basis = np.zeros((n_atomic + 6, n_modes + 6))
        basis[:n_atomic, :n_modes] = atomic_basis
        basis[n_atomic:, n_modes:] = np.diag(1.0 / np.sqrt(strain_weights))
        combined_weights = np.zeros((n_modes + 6, 2))
        combined_weights[:n_modes, :] = weights
        return cls(
            reference=np.zeros(n_atomic + 6), basis=basis,
            metric_weights=np.r_[np.repeat(modes.masses_amu, 3), strain_weights],
            axis_weights=combined_weights, axis_labels=axis_labels,
            amplitude_unit="sqrt(amu)*angstrom", reference_id=reference_id,
        )

    @property
    def axis_vectors(self) -> Array:
        """Return the two displacement directions in physical coordinates."""

        return self.basis @ self.axis_weights

    def _coordinates(self, coordinates: Array) -> tuple[Array, bool]:
        values = np.asarray(coordinates, dtype=float)
        single = values.ndim == 1
        if single:
            values = values[None, :]
        if values.ndim != 2 or values.shape[1] != self.reference.size:
            raise ValueError(f"coordinates must have shape ({self.reference.size},) or (n, {self.reference.size})")
        if not np.all(np.isfinite(values)):
            raise ValueError("coordinates contain non-finite values")
        return values, single

    def modal_amplitudes(self, coordinates: Array) -> Array:
        """Project coordinates onto every source mode; units are declared."""

        values, single = self._coordinates(coordinates)
        amplitudes = ((values - self.reference) * self.metric_weights) @ self.basis
        return amplitudes[0] if single else amplitudes

    def project(self, coordinates: Array) -> Array:
        """Return signed ``(Q1, Q2)`` for one or many configurations."""

        amplitudes = self.modal_amplitudes(coordinates)
        return amplitudes @ self.axis_weights

    def frozen_coordinates(self, q: Array) -> Array:
        """Reconstruct the frozen two-axis slice with all other modes zero."""

        values = np.asarray(q, dtype=float)
        single = values.ndim == 1
        if single:
            values = values[None, :]
        if values.ndim != 2 or values.shape[1] != 2 or not np.all(np.isfinite(values)):
            raise ValueError("q must have shape (2,) or (n, 2), with finite values")
        result = self.reference + values @ self.axis_vectors.T
        return result[0] if single else result

    def residual_norm(self, coordinates: Array) -> Array | float:
        """Metric norm of displacement outside the two selected axes."""

        values, single = self._coordinates(coordinates)
        residual = values - self.frozen_coordinates(self.project(values))
        norms = np.sqrt(np.sum(self.metric_weights * residual**2, axis=1))
        return float(norms[0]) if single else norms

    def require_representable_endpoints(
        self, initial: Array, final: Array, *, tolerance: float = 1e-8
    ) -> tuple[Array, Array]:
        """Reject a strict two-axis pathway that cannot connect its endpoints.

        Conditional relaxation is a different construction: its off-plane
        coordinates may be nonzero at either endpoint. This gate applies to
        frozen/strict-two-axis paths only.
        """

        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError("endpoint representability tolerance must be finite and nonnegative")
        if np.asarray(initial).ndim != 1 or np.asarray(final).ndim != 1:
            raise ValueError("initial and final endpoints must each be one coordinate vector")
        endpoints = np.vstack([self._coordinates(initial)[0][0], self._coordinates(final)[0][0]])
        residuals = self.residual_norm(endpoints)
        if np.any(residuals > tolerance):
            raise ValueError(
                "strict Q1/Q2 plane does not contain both endpoints: "
                f"residuals {residuals.tolist()} {self.amplitude_unit}; "
                "add the missing mode/strain directions or use conditional relaxation"
            )
        return self.project(endpoints[0]), self.project(endpoints[1])


@dataclass(frozen=True)
class ModeEnergyGrid:
    """Raw, unfitted ``Q1``–``Q2`` energies with explicit energy semantics."""

    q1: Array
    q2: Array
    energies: Array
    axis_labels: tuple[str, str]
    amplitude_unit: str
    reference_id: str
    energy_unit: str
    potential_kind: str
    boundary_condition: str

    def __post_init__(self) -> None:
        q1 = np.asarray(self.q1, dtype=float).reshape(-1)
        q2 = np.asarray(self.q2, dtype=float).reshape(-1)
        energies = np.asarray(self.energies, dtype=float)
        if q1.size < 2 or q2.size < 2 or not np.all(np.diff(q1) > 0.0) or not np.all(np.diff(q2) > 0.0):
            raise ValueError("surface axes must be strictly increasing")
        if energies.shape != (len(q2), len(q1)):
            raise ValueError("surface energies must have shape (len(q2), len(q1))")
        if not all(np.all(np.isfinite(item)) for item in (q1, q2, energies)):
            raise ValueError("surface arrays must be finite")
        if len(self.axis_labels) != 2 or any(not label.strip() for label in self.axis_labels):
            raise ValueError("surface requires two nonempty axis labels")
        if not self.amplitude_unit.strip() or not self.reference_id.strip():
            raise ValueError("surface requires amplitude_unit and reference_id")
        if not self.energy_unit.strip() or self.potential_kind not in ("E", "H"):
            raise ValueError("surface requires an energy_unit and potential_kind E or H")
        if not self.boundary_condition.strip():
            raise ValueError("surface requires an explicit boundary_condition")
        for name, value in (("q1", q1), ("q2", q2), ("energies", energies)):
            frozen = value.copy()
            frozen.flags.writeable = False
            object.__setattr__(self, name, frozen)


@dataclass(frozen=True)
class FrozenModeSurface(ModeEnergyGrid):
    """Raw energies with every coordinate outside the two axes frozen."""


@dataclass(frozen=True)
class ConditionalModeSurface(ModeEnergyGrid):
    """Raw locally stationary candidates at fixed ``Q1`` and ``Q2``.

    Each grid point is solved independently from the same declared starts.
    This is not a certified global lower envelope or a continuous branch.
    """

    orthogonal_gradient_norms: Array
    n_evaluations: Array
    selected_starts: Array
    coordinates: Array

    def __post_init__(self) -> None:
        super().__post_init__()
        grid_shape = self.energies.shape
        arrays = {
            "orthogonal_gradient_norms": np.asarray(self.orthogonal_gradient_norms, dtype=float),
            "n_evaluations": np.asarray(self.n_evaluations, dtype=int),
            "selected_starts": np.asarray(self.selected_starts, dtype=int),
            "coordinates": np.asarray(self.coordinates, dtype=float),
        }
        if any(arrays[name].shape != grid_shape for name in
               ("orthogonal_gradient_norms", "n_evaluations", "selected_starts")):
            raise ValueError("conditional diagnostics must match the energy grid shape")
        if arrays["coordinates"].shape[:2] != grid_shape or arrays["coordinates"].ndim != 3:
            raise ValueError("conditional coordinates must have shape (len(q2), len(q1), coordinate_count)")
        if not all(np.all(np.isfinite(value)) for value in arrays.values()):
            raise ValueError("conditional diagnostics must be finite")
        if np.any(arrays["orthogonal_gradient_norms"] < 0.0) or np.any(arrays["n_evaluations"] < 1):
            raise ValueError("conditional gradient norms and evaluation counts are invalid")
        if np.any(arrays["selected_starts"] < 0):
            raise ValueError("conditional selected-start indices are invalid")
        for name, value in arrays.items():
            frozen = value.copy()
            frozen.flags.writeable = False
            object.__setattr__(self, name, frozen)


@dataclass(frozen=True)
class ConditionalSurfaceContinuityAudit:
    """Neighboring off-plane jumps; a screen, not proof of one smooth branch."""

    q1_neighbor_jumps: Array
    q2_neighbor_jumps: Array
    maximum_orthogonal_jump: float
    allowed_jump: float
    within_bound: bool
    selected_start_switches: int


def audit_conditional_surface_continuity(
    plane: ModePlane,
    surface: ConditionalModeSurface,
    *,
    maximum_allowed_orthogonal_jump: float,
) -> ConditionalSurfaceContinuityAudit:
    """Screen all grid-neighbor changes outside the declared Q1/Q2 plane.

    The metric norm removes the prescribed Q step before measuring the
    remaining structural change. A small jump cannot exclude two nearly
    identical branches; a large jump requires additional branch analysis.
    """

    limit = float(maximum_allowed_orthogonal_jump)
    if not np.isfinite(limit) or limit <= 0.0:
        raise ValueError("maximum_allowed_orthogonal_jump must be finite and positive")
    if (surface.reference_id != plane.reference_id
            or surface.axis_labels != plane.axis_labels
            or surface.amplitude_unit != plane.amplitude_unit
            or surface.coordinates.shape[-1] != plane.reference.size):
        raise ValueError("conditional surface and mode plane have different coordinate contracts")
    projected = plane.project(surface.coordinates.reshape(-1, plane.reference.size))
    projected = projected.reshape((*surface.energies.shape, 2))
    q1_expected, q2_expected = np.meshgrid(surface.q1, surface.q2)
    if not np.allclose(projected, np.stack([q1_expected, q2_expected], axis=-1), rtol=0.0, atol=1e-8):
        raise ValueError("conditional surface coordinates drifted from their declared Q grid")

    def jumps(axis: int, q_step: Array, direction: Array) -> Array:
        delta = np.diff(surface.coordinates, axis=axis)
        if axis == 1:
            prescribed = q_step[None, :, None] * direction[None, None, :]
        else:
            prescribed = q_step[:, None, None] * direction[None, None, :]
        residual = delta - prescribed
        return np.sqrt(np.sum(plane.metric_weights * residual**2, axis=-1))

    q1_jumps = jumps(1, np.diff(surface.q1), plane.axis_vectors[:, 0])
    q2_jumps = jumps(0, np.diff(surface.q2), plane.axis_vectors[:, 1])
    maximum = float(max(np.max(q1_jumps), np.max(q2_jumps)))
    switches = int(np.count_nonzero(np.diff(surface.selected_starts, axis=1))
                   + np.count_nonzero(np.diff(surface.selected_starts, axis=0)))
    q1_jumps.flags.writeable = False
    q2_jumps.flags.writeable = False
    return ConditionalSurfaceContinuityAudit(
        q1_neighbor_jumps=q1_jumps, q2_neighbor_jumps=q2_jumps,
        maximum_orthogonal_jump=maximum, allowed_jump=limit,
        within_bound=maximum <= limit, selected_start_switches=switches,
    )


@dataclass(frozen=True)
class ConditionalBranchOutcome:
    """Final state of one fixed-Q orthogonal-relaxation start.

    An unconverged branch is retained as diagnostic evidence, never as a
    conditional-surface point or a competing local minimum.
    """

    start_index: int
    coordinates: Array
    energy: float
    orthogonal_gradient_norm: float
    n_evaluations: int
    converged: bool


@dataclass(frozen=True)
class ConditionalModePoint:
    """One locally stationary candidate at fixed signed ``(Q1, Q2)``.

    Gradient convergence alone does not certify positive orthogonal curvature
    or a global conditional minimum; these are separate scientific checks.
    """

    q: Array
    coordinates: Array
    energy: float
    orthogonal_gradient_norm: float
    n_evaluations: int
    n_starts: int
    selected_start: int
    branch_outcomes: tuple[ConditionalBranchOutcome, ...] = ()


@dataclass(frozen=True)
class OrthogonalCurvatureAudit:
    """Finite-difference curvature in the declared metric-orthogonal space.

    A positive minimum eigenvalue establishes only *local* stability for the
    sampled coordinates and finite-difference step. It does not prove a
    global conditional minimum or continuity with adjacent grid points.
    """

    eigenvalues: Array
    step_amplitude: float
    n_gradient_evaluations: int

    @property
    def minimum_eigenvalue(self) -> float:
        return float(np.min(self.eigenvalues)) if self.eigenvalues.size else float("inf")


@dataclass(frozen=True)
class PairedEnergyCurvatureAudit:
    """Central energy curvatures along the *sampled* orthogonal basis vectors.

    This is an independent energy-side check of gradient-derived Hessians.
    Positive basis diagonals do not establish a positive-definite Hessian:
    mixed directions can still have negative curvature.
    """

    steps: tuple[float, float]
    curvature_by_step: Array
    absolute_step_disagreement: Array

    @property
    def minimum_sampled_curvature(self) -> float:
        return float(np.min(self.curvature_by_step))

    @property
    def maximum_step_disagreement(self) -> float:
        return float(np.max(self.absolute_step_disagreement))


@dataclass(frozen=True)
class DirectionalCurvatureConsistencyAudit:
    """Energy-versus-gradient finite-difference diagnostics, not certification.

    The signed energies and directional gradients must come from the same
    physical coordinate line and calculator contract. Nonzero disagreement
    can reflect finite-step anharmonicity, electronic noise, basis/Pulay
    effects, or a wrong energy-gradient adapter; this object cannot infer
    which cause applies.
    """

    steps: tuple[float, float]
    energy_curvatures: Array
    gradient_curvatures: Array
    absolute_energy_gradient_disagreement: Array


@dataclass(frozen=True)
class DirectionalWorkConsistencyAudit:
    """Energy-minus-integrated-gradient residuals on a five-point line.

    The two three-point Simpson integrals use the center and both signed
    probe radii. A cubic directional gradient is integrated exactly, so a
    persistent residual cannot be dismissed as ordinary harmonic curvature.
    It may still reflect higher-order anharmonicity or calculator noise; this
    diagnostic does not assign a cause or certify a local minimum.
    """

    abscissas: Array
    energy_increments: Array
    simpson_gradient_work: Array
    energy_minus_gradient_work: Array

    @property
    def maximum_absolute_work_residual(self) -> float:
        return float(np.max(np.abs(self.energy_minus_gradient_work)))


@dataclass(frozen=True)
class CurvatureEnergyResolutionScreen:
    """Compare sampled curvature signals with observed origin sensitivity.

    This is deliberately a *screen*, not a numerical error bound. A rigid
    translation samples a different direction from a curvature probe; its
    energy span cannot be subtracted from or assigned as the uncertainty of
    the curvature. Ratios above one only flag that the sampled origin effect
    is at least as large as the energy scale being interpreted.
    """

    steps: tuple[float, float]
    energy_double_differences_eV: Array
    energy_gradient_disagreement_scales_eV: Array
    observed_origin_energy_span_eV: float
    origin_span_over_energy_signal: Array
    origin_span_over_energy_gradient_disagreement: Array


def screen_curvature_energy_resolution(
    curvature: DirectionalCurvatureConsistencyAudit,
    observed_origin_energy_span_eV: float,
) -> CurvatureEnergyResolutionScreen:
    """Expose whether a measured translation span rivals two curvature signals.

    The energy curvature signal is ``|E(+h)+E(-h)-2E(0)|=|K_E|h²``.
    The disagreement scale is ``|K_E-K_gradient|h²``. The caller must verify
    common calculator settings and physical provenance of both datasets.
    A small ratio does **not** certify curvature or PES accuracy.
    """

    if not isinstance(curvature, DirectionalCurvatureConsistencyAudit):
        raise TypeError("curvature must be a directional energy-gradient audit")
    steps = np.asarray(curvature.steps, dtype=float)
    energy_curvatures = np.asarray(curvature.energy_curvatures, dtype=float)
    disagreement = np.asarray(curvature.absolute_energy_gradient_disagreement, dtype=float)
    span = float(observed_origin_energy_span_eV)
    if (steps.shape != (2,) or np.any(~np.isfinite(steps)) or np.any(steps <= 0.0)
            or energy_curvatures.shape != (2,) or disagreement.shape != (2,)
            or np.any(~np.isfinite(energy_curvatures)) or np.any(~np.isfinite(disagreement))
            or np.any(disagreement < 0.0) or not np.isfinite(span) or span < 0.0):
        raise ValueError("two finite curvature steps and a nonnegative origin span are required")
    signal = np.abs(energy_curvatures) * steps**2
    disagreement_scale = disagreement * steps**2

    def ratio(values: Array) -> Array:
        result = np.zeros_like(values)
        np.divide(span, values, out=result, where=values > 0.0)
        result[(values == 0.0) & (span > 0.0)] = np.inf
        result.flags.writeable = False
        return result

    for array in (signal, disagreement_scale):
        array.flags.writeable = False
    return CurvatureEnergyResolutionScreen(
        tuple(float(value) for value in steps), signal, disagreement_scale,
        span, ratio(signal), ratio(disagreement_scale),
    )


def audit_directional_work_consistency(
    abscissas: Array, energies: Array, directional_gradients: Array,
) -> DirectionalWorkConsistencyAudit:
    """Compare E differences with dE/dq over [-2h,0] and [0,2h].

    Inputs must be ordered at ``(-2h,-h,0,+h,+2h)`` with the same geometry
    line and calculator settings. Provenance and numerical convergence are
    caller responsibilities; this function introduces no acceptance threshold.
    """

    x = np.asarray(abscissas, dtype=float).reshape(-1)
    values = np.asarray(energies, dtype=float).reshape(-1)
    gradients = np.asarray(directional_gradients, dtype=float).reshape(-1)
    if (x.shape != (5,) or values.shape != (5,) or gradients.shape != (5,)
            or not all(np.all(np.isfinite(array)) for array in (x, values, gradients))
            or not np.all(np.diff(x) > 0.0)):
        raise ValueError("five finite, strictly ordered line points are required")
    h = float(x[3])
    if h <= 0.0 or not np.allclose(x, [-2*h, -h, 0.0, h, 2*h], rtol=0.0, atol=1e-10):
        raise ValueError("line points must be equally spaced at (-2h,-h,0,+h,+2h)")
    increments = np.array([values[2] - values[0], values[4] - values[2]])
    work = h / 3.0 * np.array([
        gradients[0] + 4.0 * gradients[1] + gradients[2],
        gradients[2] + 4.0 * gradients[3] + gradients[4],
    ])
    residuals = increments - work
    for array in (x, increments, work, residuals):
        array.flags.writeable = False
    return DirectionalWorkConsistencyAudit(x, increments, work, residuals)


def audit_directional_curvature_consistency(
    center_energy: float,
    pair_energies: Array,
    pair_directional_gradients: Array,
    steps: tuple[float, float],
) -> DirectionalCurvatureConsistencyAudit:
    """Cross-check one metric-unit direction at two steps using ``(+,-)`` data.

    Arrays have shape ``(2 steps, 2 signs)``. The direction and its sign must
    be identical across steps. Physical provenance and SCF convergence are
    caller responsibilities; no tolerance is silently imposed.
    """

    energies = np.asarray(pair_energies, dtype=float)
    gradients = np.asarray(pair_directional_gradients, dtype=float)
    if (energies.shape != (2, 2) or gradients.shape != (2, 2)
            or not np.all(np.isfinite(gradients))):
        raise ValueError("two finite (+,-) energy and directional-gradient pairs are required")
    energy_audit = audit_paired_energy_curvature(
        center_energy, energies[:1], steps[0], energies[1:], steps[1],
    )
    derivative = (gradients[:, 0] - gradients[:, 1]) / (2.0 * np.asarray(steps, dtype=float))
    energy_values = energy_audit.curvature_by_step[:, 0].copy()
    disagreement = np.abs(energy_values - derivative)
    for values in (energy_values, derivative, disagreement):
        values.flags.writeable = False
    return DirectionalCurvatureConsistencyAudit(
        energy_audit.steps, energy_values, derivative, disagreement,
    )


def audit_paired_energy_curvature(
    center_energy: float,
    first_pair_energies: Array,
    first_step: float,
    second_pair_energies: Array,
    second_step: float,
) -> PairedEnergyCurvatureAudit:
    """Compare two central energy differences using the same center and axes.

    Each pair array has shape ``(n_directions, 2)`` with ``(+,-)`` energies.
    The caller must independently verify that the probes really are those
    fixed-Q, metric-unit directions and share one calculator contract.
    """

    center = float(center_energy)
    steps = (float(first_step), float(second_step))
    pairs = tuple(np.asarray(values, dtype=float) for values in
                  (first_pair_energies, second_pair_energies))
    if not np.isfinite(center) or any(not np.isfinite(step) or step <= 0 for step in steps):
        raise ValueError("center energy and positive finite steps are required")
    if steps[0] == steps[1]:
        raise ValueError("the two finite-difference steps must differ")
    if (pairs[0].ndim != 2 or pairs[0].shape[1] != 2 or pairs[0].shape[0] < 1
            or pairs[1].shape != pairs[0].shape
            or any(not np.all(np.isfinite(values)) for values in pairs)):
        raise ValueError("both energy-pair arrays must be finite and have matching (n, 2) shapes")
    curvatures = np.stack([
        (values[:, 0] + values[:, 1] - 2.0 * center) / step**2
        for values, step in zip(pairs, steps)
    ])
    disagreement = np.abs(curvatures[0] - curvatures[1])
    curvatures.flags.writeable = False
    disagreement.flags.writeable = False
    return PairedEnergyCurvatureAudit(steps, curvatures, disagreement)


def _orthogonal_directions(
    plane: ModePlane, frozen_directions: Array | None = None,
) -> Array:
    """Metric-orthonormal directions outside the axes and fixed gauge space."""

    weighted_axes = np.sqrt(plane.metric_weights)[:, None] * plane.axis_vectors
    if frozen_directions is None:
        weighted_constraints = weighted_axes
    else:
        directions = np.asarray(frozen_directions, dtype=float)
        if (directions.ndim != 2 or directions.shape[0] != plane.reference.size
                or directions.shape[1] < 1 or not np.all(np.isfinite(directions))):
            raise ValueError("frozen_directions must be finite columns in physical coordinate space")
        weighted_frozen = np.sqrt(plane.metric_weights)[:, None] * directions
        lengths = np.linalg.norm(weighted_frozen, axis=0)
        if np.any(lengths <= 0.0):
            raise ValueError("frozen_directions cannot contain zero vectors")
        weighted_frozen /= lengths
        if np.max(np.abs(weighted_axes.T @ weighted_frozen)) > 1e-8:
            raise ValueError("frozen_directions must be metric-orthogonal to both Q axes")
        weighted_constraints = np.column_stack([weighted_axes, weighted_frozen])
        if weighted_constraints.shape[1] > plane.reference.size:
            raise ValueError("more constrained directions than physical coordinates")
    left, singular, _ = np.linalg.svd(weighted_constraints, full_matrices=True)
    if singular.size and np.min(singular) < 1e-8:
        raise ValueError("frozen_directions must be linearly independent")
    return left[:, weighted_constraints.shape[1]:] / np.sqrt(plane.metric_weights)[:, None]


def audit_orthogonal_curvature(
    plane: ModePlane,
    coordinates: Array,
    energy_and_gradient: Callable[[Array], tuple[float, Array]],
    *,
    step_amplitude: float,
    frozen_directions: Array | None = None,
) -> OrthogonalCurvatureAudit:
    """Probe the complete local orthogonal Hessian with raw physical gradients.

    Costs two evaluator calls per orthogonal coordinate; therefore run this
    deliberately on selected material points, not silently on every DFT grid
    point. Directions are unit vectors in the declared coordinate metric.
    Each displaced probe preserves Q1 and Q2. A negative eigenvalue rejects a
    stationary candidate as a conditional local minimum.
    """

    values, single = plane._coordinates(coordinates)
    if not single or not callable(energy_and_gradient):
        raise ValueError("one coordinate vector and an energy/gradient evaluator are required")
    step = float(step_amplitude)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("step_amplitude must be finite and positive")
    directions = _orthogonal_directions(plane, frozen_directions)
    n_orthogonal = directions.shape[1]
    hessian = np.empty((n_orthogonal, n_orthogonal), dtype=float)
    for index in range(n_orthogonal):
        projected = []
        for sign in (1.0, -1.0):
            trial = values[0] + sign * step * directions[:, index]
            energy, gradient = energy_and_gradient(trial)
            gradient = np.asarray(gradient, dtype=float).reshape(-1)
            if (not np.isfinite(float(energy)) or gradient.shape != values[0].shape
                    or not np.all(np.isfinite(gradient))):
                raise ValueError("curvature evaluator must return finite energy and a full gradient")
            projected.append(directions.T @ gradient)
        hessian[:, index] = (projected[0] - projected[1]) / (2.0 * step)
    eigenvalues = np.linalg.eigvalsh(0.5 * (hessian + hessian.T))
    eigenvalues.flags.writeable = False
    return OrthogonalCurvatureAudit(
        eigenvalues=eigenvalues, step_amplitude=step,
        n_gradient_evaluations=2 * n_orthogonal,
    )


def orthogonal_branch_seeds(
    plane: ModePlane,
    coordinates: Array,
    direction: Array,
    *,
    amplitude: float,
    frozen_directions: Array | None = None,
    orthogonal_amplitude_bound: float | None = None,
) -> tuple[Array, Array]:
    """Create signed fixed-Q starts along one audited unstable direction.

    This performs geometry/gauge checks only. The caller must separately
    establish that ``direction`` comes from a trustworthy negative-curvature
    calculation; the seeds are not energies or certified minima.
    """

    values, single = plane._coordinates(coordinates)
    vector = np.asarray(direction, dtype=float).reshape(-1)
    step = float(amplitude)
    if not single or vector.shape != values[0].shape or not np.all(np.isfinite(vector)):
        raise ValueError("one finite center and a matching finite branch direction are required")
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("branch amplitude must be finite and positive")
    if orthogonal_amplitude_bound is not None and (
        not np.isfinite(orthogonal_amplitude_bound) or orthogonal_amplitude_bound <= 0.0
    ):
        raise ValueError("orthogonal amplitude bound must be finite and positive")
    orthogonal = _orthogonal_directions(plane, frozen_directions)
    coefficients = orthogonal.T @ (plane.metric_weights * vector)
    represented = orthogonal @ coefficients
    metric_norm = float(np.sqrt(np.dot(plane.metric_weights * vector, vector)))
    residual_norm = float(np.sqrt(np.dot(plane.metric_weights * (vector - represented), vector - represented)))
    if abs(metric_norm - 1.0) > 1e-6 or residual_norm > 1e-6:
        raise ValueError("branch direction must be unit-normalized inside the open fixed-Q subspace")
    target_q = plane.project(values[0])
    seeds = (values[0] + step * vector, values[0] - step * vector)
    if orthogonal_amplitude_bound is not None:
        base = plane.frozen_coordinates(target_q)
        if any(np.any(np.abs(orthogonal.T @ (plane.metric_weights * (seed - base)))
                      > orthogonal_amplitude_bound + 1e-10) for seed in seeds):
            raise ValueError("branch seed exceeds the orthogonal amplitude bound")
    if any(not np.allclose(plane.project(seed), target_q, atol=1e-8, rtol=0.0) for seed in seeds):
        raise ValueError("branch seed changes its fixed order parameters")
    return seeds


def _safeguarded_bfgs(
    evaluate: Callable[[Array], tuple[float, Array]],
    valid: Callable[[Array], bool],
    start: Array,
    *,
    gradient_tolerance: float,
    max_iterations: int,
    initial_trust_radius: float,
) -> tuple[Array, int]:
    """BFGS with geometry-screened Armijo backtracking in open coordinates."""

    z = np.asarray(start, dtype=float).copy()
    if not valid(z):
        raise ValueError("supplied safeguarded-BFGS start violates the geometry/bound guard")
    energy, gradient = evaluate(z)
    calls = 1
    inverse_hessian = np.eye(z.size)
    radius = initial_trust_radius
    for _ in range(max_iterations):
        if np.linalg.norm(gradient) <= gradient_tolerance:
            break
        direction = -inverse_hessian @ gradient
        if np.dot(direction, gradient) >= -1e-12 * np.linalg.norm(direction) * np.linalg.norm(gradient):
            direction = -gradient
            inverse_hessian = np.eye(z.size)
        norm = np.linalg.norm(direction)
        if norm <= 0.0:
            break
        direction *= min(1.0, radius / norm)
        trial = None
        step = 1.0
        for _ in range(20):
            proposed = z + step * direction
            if valid(proposed):
                proposed_energy, proposed_gradient = evaluate(proposed)
                calls += 1
                if proposed_energy <= energy + 1e-4 * step * np.dot(gradient, direction):
                    trial = (proposed, proposed_energy, proposed_gradient)
                    break
            step *= 0.5
        if trial is None:
            break
        next_z, next_energy, next_gradient = trial
        displacement = next_z - z
        gradient_change = next_gradient - gradient
        curvature = float(np.dot(displacement, gradient_change))
        if curvature > 1e-10 * np.linalg.norm(displacement) * np.linalg.norm(gradient_change):
            rho = 1.0 / curvature
            transform = np.eye(z.size) - rho * np.outer(displacement, gradient_change)
            inverse_hessian = (transform @ inverse_hessian @ transform.T
                               + rho * np.outer(displacement, displacement))
        z, energy, gradient = next_z, next_energy, next_gradient
        if step >= 0.75:
            radius = min(2.0 * initial_trust_radius, 1.25 * radius)
    return z, calls


def relax_orthogonal_at_q(
    plane: ModePlane,
    q: Array,
    energy_and_gradient: Callable[[Array], tuple[float, Array]],
    *,
    starts: Sequence[Array] = (),
    include_frozen_start: bool = True,
    frozen_directions: Array | None = None,
    orthogonal_amplitude_bound: float | None = None,
    gradient_tolerance: float = 1e-6,
    max_iterations: int = 200,
    optimizer: str = "L-BFGS-B",
    trial_validator: Callable[[Array], bool] | None = None,
    initial_trust_radius: float = 0.25,
) -> ConditionalModePoint:
    """Minimize all coordinates orthogonal to two fixed order parameters.

    ``energy_and_gradient`` returns physical energy and its gradient with
    respect to the *physical coordinates* used by ``plane``. A backend adapter
    must consistently transform raw atomic forces and cell stress into this
    dual gradient; NEB spring/climbing forces are not appropriate. The result
    is the lowest-energy gradient-converged candidate among supplied starting
    points. ``frozen_directions`` additionally fixes explicit gauge or other
    directions at their reference values. It is not certified as a minimum
    rather than a saddle, nor as the global conditional minimum. Set
    ``include_frozen_start=False`` only when deliberately continuing from one
    or more audited prior candidates; the caller then assumes responsibility
    for checking competing branches. SciPy is an optional extra.

    ``safeguarded_bfgs`` is an opt-in backend-independent optimizer for
    geometry-limited DFT workflows: it screens each proposed step with
    ``trial_validator`` *before* calling the energy evaluator, then Armijo
    backtracks within an initial metric trust radius. A rejected trial is
    neither a DFT point nor a changed physical Hamiltonian. The default
    remains L-BFGS-B for compatibility; neither optimizer certifies a local
    minimum without an independent curvature check.
    """

    try:
        from scipy.optimize import minimize
    except ImportError as exc:
        raise ImportError("conditional mode relaxation requires the 'mode' extra: pip install varneb[mode]") from exc
    point = np.asarray(q, dtype=float).reshape(-1)
    if point.shape != (2,) or not np.all(np.isfinite(point)):
        raise ValueError("q must contain two finite order-parameter values")
    if not callable(energy_and_gradient):
        raise TypeError("energy_and_gradient must be callable")
    if not np.isfinite(gradient_tolerance) or gradient_tolerance <= 0.0:
        raise ValueError("gradient_tolerance must be finite and positive")
    if isinstance(max_iterations, bool) or int(max_iterations) != max_iterations or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")
    if orthogonal_amplitude_bound is not None and (
        not np.isfinite(orthogonal_amplitude_bound) or orthogonal_amplitude_bound <= 0.0
    ):
        raise ValueError("orthogonal_amplitude_bound must be finite and positive")
    if optimizer not in ("L-BFGS-B", "safeguarded_bfgs"):
        raise ValueError("optimizer must be L-BFGS-B or safeguarded_bfgs")
    if optimizer == "safeguarded_bfgs" and not callable(trial_validator):
        raise ValueError("safeguarded_bfgs requires a callable trial_validator")
    if not np.isfinite(initial_trust_radius) or initial_trust_radius <= 0.0:
        raise ValueError("initial_trust_radius must be finite and positive")
    if not isinstance(include_frozen_start, bool) or (not include_frozen_start and len(starts) == 0):
        raise ValueError("include_frozen_start=False requires at least one supplied start")
    base = plane.frozen_coordinates(point)
    orthogonal = _orthogonal_directions(plane, frozen_directions)
    supplied = ([base] if include_frozen_start else []) + list(starts)
    candidates: list[ConditionalModePoint] = []
    branch_outcomes: list[ConditionalBranchOutcome] = []

    def evaluate(z: Array) -> tuple[float, Array]:
        coordinates = base + orthogonal @ z
        value, gradient = energy_and_gradient(coordinates)
        energy_value = float(value)
        gradient = np.asarray(gradient, dtype=float).reshape(-1)
        if gradient.shape != base.shape or not np.isfinite(energy_value) or not np.all(np.isfinite(gradient)):
            raise ValueError("energy_and_gradient must return finite energy and a full physical gradient")
        return energy_value, orthogonal.T @ gradient

    def valid(z: Array) -> bool:
        if orthogonal_amplitude_bound is not None and np.any(
            np.abs(z) > orthogonal_amplitude_bound + 1e-10
        ):
            return False
        return trial_validator is None or bool(trial_validator(base + orthogonal @ z))

    for index, start in enumerate(supplied):
        coordinates = np.asarray(start, dtype=float).reshape(-1)
        if coordinates.shape != base.shape or not np.all(np.isfinite(coordinates)):
            raise ValueError("every start must be a finite full coordinate vector")
        z0 = orthogonal.T @ (plane.metric_weights * (coordinates - base))
        if (orthogonal_amplitude_bound is not None
                and np.any(np.abs(z0) > orthogonal_amplitude_bound + 1e-10)):
            raise ValueError("a supplied start exceeds the orthogonal amplitude bound")
        if orthogonal.shape[1] and optimizer == "safeguarded_bfgs":
            z, n_evaluations = _safeguarded_bfgs(
                evaluate, valid, z0, gradient_tolerance=gradient_tolerance,
                max_iterations=max_iterations, initial_trust_radius=initial_trust_radius,
            )
            n_evaluations += 1  # final independent verification
        elif orthogonal.shape[1]:
            result = minimize(
                evaluate, z0, method="L-BFGS-B", jac=True,
                bounds=(None if orthogonal_amplitude_bound is None else
                        [(-orthogonal_amplitude_bound, orthogonal_amplitude_bound)] * orthogonal.shape[1]),
                # SciPy's gtol is an infinity-norm criterion, whereas the
                # public contract and final gate use the Euclidean norm.
                # Tighten the component criterion so an early SciPy success
                # cannot still violate our requested full-vector tolerance.
                options={"gtol": gradient_tolerance / np.sqrt(orthogonal.shape[1]),
                         "maxiter": int(max_iterations)},
            )
            z = result.x
            n_evaluations = int(result.nfev) + 1  # final independent verification
        else:
            z = z0
            n_evaluations = 1
        value, projected_gradient = evaluate(z)
        residual = float(np.linalg.norm(projected_gradient))
        branch_outcomes.append(ConditionalBranchOutcome(
            start_index=index, coordinates=base + orthogonal @ z,
            energy=value, orthogonal_gradient_norm=residual,
            n_evaluations=n_evaluations, converged=residual <= gradient_tolerance,
        ))
        if residual <= gradient_tolerance:
            candidate = ConditionalModePoint(
                q=point.copy(), coordinates=base + orthogonal @ z,
                energy=value, orthogonal_gradient_norm=residual,
                n_evaluations=n_evaluations, n_starts=len(supplied), selected_start=index,
            )
            candidates.append(candidate)
    if not candidates:
        raise RuntimeError("no starting point reached the orthogonal-gradient tolerance")
    best = min(candidates, key=lambda candidate: candidate.energy)
    if not np.allclose(plane.project(best.coordinates), point, rtol=0.0, atol=1e-8):
        raise RuntimeError("conditional optimizer drifted away from fixed Q1/Q2")
    return ConditionalModePoint(
        q=best.q, coordinates=best.coordinates, energy=best.energy,
        orthogonal_gradient_norm=best.orthogonal_gradient_norm,
        n_evaluations=best.n_evaluations, n_starts=best.n_starts,
        selected_start=best.selected_start,
        branch_outcomes=tuple(branch_outcomes),
    )


def sample_frozen_mode_surface(
    plane: ModePlane,
    q1: Array,
    q2: Array,
    energy: Callable[[Array], float],
    *,
    energy_unit: str,
    potential_kind: str,
    boundary_condition: str,
) -> FrozenModeSurface:
    """Evaluate a user-supplied energy function on a frozen two-mode plane.

    The returned array is indexed ``energies[j_q2, i_q1]``. The evaluator gets
    one *physical coordinate vector* per point and may wrap any ASE calculator.
    It must use the same Hamiltonian and boundary condition as the compared
    VCNEB path. The required metadata labels the returned grid but cannot
    independently certify that the evaluator obeyed the declared contract.
    """

    axes = []
    for label, values in (("q1", q1), ("q2", q2)):
        grid = np.asarray(values, dtype=float).reshape(-1)
        if grid.size < 2 or not np.all(np.isfinite(grid)) or not np.all(np.diff(grid) > 0.0):
            raise ValueError(f"{label} must contain at least two strictly increasing finite values")
        axes.append(grid)
    if not callable(energy):
        raise TypeError("energy must be callable")
    if not energy_unit.strip() or potential_kind not in ("E", "H") or not boundary_condition.strip():
        raise ValueError("energy_unit, potential_kind E/H, and boundary_condition are required")
    values = np.empty((len(axes[1]), len(axes[0])), dtype=float)
    for j, second in enumerate(axes[1]):
        for i, first in enumerate(axes[0]):
            coordinate = plane.frozen_coordinates(np.array([first, second]))
            result = float(energy(coordinate))
            if not np.isfinite(result):
                raise ValueError(f"energy evaluator returned a non-finite value at grid ({i}, {j})")
            values[j, i] = result
    return FrozenModeSurface(
        q1=axes[0].copy(), q2=axes[1].copy(), energies=values,
        axis_labels=plane.axis_labels, amplitude_unit=plane.amplitude_unit,
        reference_id=plane.reference_id, energy_unit=energy_unit,
        potential_kind=potential_kind, boundary_condition=boundary_condition,
    )


def sample_conditional_mode_surface(
    plane: ModePlane,
    q1: Array,
    q2: Array,
    energy_and_gradient: Callable[[Array], tuple[float, Array]],
    *,
    energy_unit: str,
    potential_kind: str,
    boundary_condition: str,
    starts: Sequence[Array] = (),
    frozen_directions: Array | None = None,
    orthogonal_amplitude_bound: float | None = None,
    gradient_tolerance: float = 1e-6,
    max_iterations: int = 200,
) -> ConditionalModeSurface:
    """Evaluate independent fixed-``Q1,Q2`` orthogonal relaxations on a grid.

    This low-level serial function does not launch or schedule DFT jobs. It
    fails rather than silently filling a nonconverged point. The caller must
    verify positive orthogonal curvature, competing branches and the physical
    identity of every material structure before calling the result a PES.
    """

    # Validate both axes and physical labels before the first evaluation.
    axes = []
    for label, values in (("q1", q1), ("q2", q2)):
        grid = np.asarray(values, dtype=float).reshape(-1)
        if grid.size < 2 or not np.all(np.isfinite(grid)) or not np.all(np.diff(grid) > 0.0):
            raise ValueError(f"{label} must contain at least two strictly increasing finite values")
        axes.append(grid)
    if not energy_unit.strip() or potential_kind not in ("E", "H") or not boundary_condition.strip():
        raise ValueError("energy_unit, potential_kind E/H, and boundary_condition are required")
    shape = (len(axes[1]), len(axes[0]))
    energies = np.empty(shape)
    norms = np.empty(shape)
    evaluations = np.empty(shape, dtype=int)
    selected = np.empty(shape, dtype=int)
    coordinates = np.empty((*shape, plane.reference.size))
    for j, second in enumerate(axes[1]):
        for i, first in enumerate(axes[0]):
            result = relax_orthogonal_at_q(
                plane, np.array([first, second]), energy_and_gradient,
                starts=starts, frozen_directions=frozen_directions,
                orthogonal_amplitude_bound=orthogonal_amplitude_bound,
                gradient_tolerance=gradient_tolerance,
                max_iterations=max_iterations,
            )
            energies[j, i] = result.energy
            norms[j, i] = result.orthogonal_gradient_norm
            evaluations[j, i] = result.n_evaluations
            selected[j, i] = result.selected_start
            coordinates[j, i] = result.coordinates
    return ConditionalModeSurface(
        q1=axes[0], q2=axes[1], energies=energies,
        axis_labels=plane.axis_labels, amplitude_unit=plane.amplitude_unit,
        reference_id=plane.reference_id, energy_unit=energy_unit,
        potential_kind=potential_kind, boundary_condition=boundary_condition,
        orthogonal_gradient_norms=norms, n_evaluations=evaluations,
        selected_starts=selected, coordinates=coordinates,
    )


__all__ = [
    "ModePlane", "ModeEnergyGrid", "FrozenModeSurface", "ConditionalBranchOutcome", "ConditionalModePoint",
    "ConditionalModeSurface", "ConditionalSurfaceContinuityAudit",
    "audit_conditional_surface_continuity", "OrthogonalCurvatureAudit", "audit_orthogonal_curvature",
    "PairedEnergyCurvatureAudit", "audit_paired_energy_curvature",
    "DirectionalCurvatureConsistencyAudit", "audit_directional_curvature_consistency",
    "DirectionalWorkConsistencyAudit", "audit_directional_work_consistency",
    "CurvatureEnergyResolutionScreen", "screen_curvature_energy_resolution",
    "sample_frozen_mode_surface",
    "sample_conditional_mode_surface", "relax_orthogonal_at_q",
]
