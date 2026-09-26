"""Mode-plane projection and frozen-surface numerical contracts."""

from __future__ import annotations

import numpy as np
import pytest

from vcneb import (
    ConditionalModeSurface, GammaModes, ModePlane,
    anchored_gamma_axis_weights,
    audit_conditional_surface_continuity, audit_orthogonal_curvature, relax_orthogonal_at_q,
    orthogonal_branch_seeds, sample_conditional_mode_surface, sample_frozen_mode_surface,
)
from vcneb.mode_surface import (
    audit_directional_curvature_consistency, audit_directional_work_consistency,
    audit_paired_energy_curvature, screen_curvature_energy_resolution,
)


def test_directional_work_integrates_quartic_energy_and_detects_wrong_gradients() -> None:
    x = np.array([-0.2, -0.1, 0.0, 0.1, 0.2])
    energy = 0.8 + 0.3*x + 0.7*x*x + 0.4*x**3 + 0.2*x**4
    gradient = 0.3 + 1.4*x + 1.2*x*x + 0.8*x**3
    audit = audit_directional_work_consistency(x, energy, gradient)
    assert audit.maximum_absolute_work_residual < 1e-15
    wrong = audit_directional_work_consistency(x, energy, gradient + 0.05*x)
    assert wrong.maximum_absolute_work_residual > 1e-4
    with pytest.raises(ValueError, match="equally spaced"):
        audit_directional_work_consistency(x + [0, 0, 0, 0.01, 0], energy, gradient)


def _identity_plane() -> ModePlane:
    weights = np.array([
        [1.0 / np.sqrt(2.0), 0.0],
        [1.0 / np.sqrt(2.0), 0.0],
        [0.0, 1.0],
    ])
    return ModePlane(
        reference=np.array([0.0, 0.0, 0.0]),
        basis=np.eye(3),
        metric_weights=np.ones(3),
        axis_weights=weights,
        axis_labels=("Q1: combined x+y", "Q2: z"),
        amplitude_unit="angstrom",
        reference_id="toy/reference/zero",
    )


def test_single_and_composite_axes_project_and_reconstruct() -> None:
    plane = _identity_plane()
    coordinates = np.array([2.0, 2.0, 3.0])
    assert np.allclose(plane.modal_amplitudes(coordinates), coordinates)
    assert np.allclose(plane.project(coordinates), [2.0 * np.sqrt(2.0), 3.0])
    assert np.allclose(plane.frozen_coordinates(plane.project(coordinates)), coordinates)
    assert plane.residual_norm(coordinates) < 1e-12
    with_residual = np.array([3.0, 1.0, 3.0])
    assert np.allclose(plane.project(with_residual), plane.project(coordinates))
    assert plane.residual_norm(with_residual) == pytest.approx(np.sqrt(2.0))
    assert np.allclose(plane.project(np.vstack([coordinates, with_residual]))[0],
                       plane.project(coordinates))
    initial_q, final_q = plane.require_representable_endpoints(np.zeros(3), coordinates)
    assert np.allclose(initial_q, [0.0, 0.0])
    assert np.allclose(final_q, plane.project(coordinates))
    with pytest.raises(ValueError, match="does not contain both endpoints"):
        plane.require_representable_endpoints(np.zeros(3), with_residual)


def test_invalid_axis_or_mode_metric_fails_before_energy_evaluation() -> None:
    plane = _identity_plane()
    arguments = {
        "reference": plane.reference,
        "basis": plane.basis,
        "metric_weights": plane.metric_weights,
        "axis_weights": plane.axis_weights,
        "axis_labels": plane.axis_labels,
        "amplitude_unit": plane.amplitude_unit,
        "reference_id": plane.reference_id,
    }
    with pytest.raises(ValueError, match="orthogonal"):
        ModePlane(**{**arguments, "axis_weights": np.array([[1.0, 1.0], [0.0, 0.0], [0.0, 0.0]])})
    with pytest.raises(ValueError, match="orthonormal"):
        ModePlane(**{**arguments, "basis": 2.0 * np.eye(3)})
    with pytest.raises(ValueError, match="positive"):
        ModePlane(**{**arguments, "metric_weights": np.array([1.0, 0.0, 1.0])})
    with pytest.raises(ValueError, match="amplitude_unit"):
        ModePlane(**{**arguments, "amplitude_unit": ""})


def test_phonopy_gamma_metric_and_mass_weighted_units() -> None:
    modes = GammaModes(
        masses_amu=np.array([4.0]),
        eigenvalues_eV_per_A2_amu=np.ones(3),
        frequencies_cm1=np.ones(3),
        eigenvectors=np.eye(3),
        translations_projected=False,
    )
    plane = ModePlane.from_gamma_modes(
        modes,
        np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]),
        axis_labels=("Qx", "Qy"),
        reference_id="toy/gamma",
    )
    assert plane.amplitude_unit == "sqrt(amu)*angstrom"
    assert np.allclose(plane.project(np.array([0.5, -0.25, 0.0])), [1.0, -0.5])
    extended = ModePlane.from_gamma_modes_with_strain(
        modes, np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]),
        strain_metric_weights_amu_A2=np.full(6, 25.0),
        axis_labels=("Qx", "Qy"), reference_id="toy/gamma+strain",
    )
    values = np.r_[[0.5, -0.25, 0.0], [0.02, 0.0, 0.0, 0.0, 0.0, 0.0]]
    assert np.allclose(extended.project(values), [1.0, -0.5])
    assert extended.residual_norm(values) == pytest.approx(0.1)
    with pytest.raises(ValueError, match="strain metric"):
        ModePlane.from_gamma_modes_with_strain(
            modes, np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]),
            strain_metric_weights_amu_A2=np.zeros(6),
            axis_labels=("Qx", "Qy"), reference_id="toy/invalid",
        )


def test_endpoint_defined_axis_is_invariant_to_degenerate_basis_rotation() -> None:
    endpoint_delta = np.array([2.0, 1.0, 0.0])
    probe = np.array([0.7, -0.3, 0.5])
    angle = 0.63
    rotation = np.array([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])
    projections = []
    for basis in (np.eye(3), rotation):
        endpoint_amplitudes = basis.T @ endpoint_delta
        weights = np.column_stack([
            np.array([*endpoint_amplitudes[:2] / np.linalg.norm(endpoint_amplitudes[:2]), 0.0]),
            np.array([0.0, 0.0, 1.0]),
        ])
        plane = ModePlane(
            reference=np.zeros(3), basis=basis, metric_weights=np.ones(3),
            axis_weights=weights, axis_labels=("endpoint-aligned doublet", "other"),
            amplitude_unit="toy length", reference_id="toy/degenerate-doublet",
        )
        projections.append(plane.project(probe))
    assert np.allclose(projections[0], projections[1], atol=1e-12)


def test_anchored_soft_subspace_axes_survive_eigenvector_rotation() -> None:
    angle = 0.47
    rotation = np.array([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])
    anchors = np.array([[0.0, 1.0], [0.0, 0.0], [1.0, 0.0]])
    axes = []
    for eigenvectors in (np.eye(3), rotation):
        modes = GammaModes(
            masses_amu=np.array([4.0]),
            eigenvalues_eV_per_A2_amu=np.full(3, -1.0),
            frequencies_cm1=np.full(3, -1.0),
            eigenvectors=eigenvectors,
            translations_projected=False,
        )
        weights = anchored_gamma_axis_weights(modes, [0, 1, 2], anchors)
        axes.append(modes.cartesian_eigenvectors() @ weights)
    assert np.allclose(axes[0], axes[1], atol=1e-12)
    assert np.allclose(axes[0], anchors / 2.0)


def test_anchored_soft_subspace_rejects_dependent_or_missing_anchors() -> None:
    modes = GammaModes(
        masses_amu=np.ones(1),
        eigenvalues_eV_per_A2_amu=np.ones(3),
        frequencies_cm1=np.ones(3),
        eigenvectors=np.eye(3),
        translations_projected=False,
    )
    with pytest.raises(ValueError, match="linearly dependent"):
        anchored_gamma_axis_weights(modes, [0, 1], np.array([[1.0, 2.0], [0.0, 0.0], [0.0, 0.0]]))
    with pytest.raises(ValueError, match="negligible projection"):
        anchored_gamma_axis_weights(modes, [0, 1], np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 1.0]]))


def test_frozen_surface_is_raw_grid_and_rejects_bad_evaluations() -> None:
    plane = _identity_plane()
    q1 = np.array([-1.0, 0.0, 1.0])
    q2 = np.array([-0.5, 0.0, 0.5])
    conditions = {"energy_unit": "eV", "potential_kind": "E", "boundary_condition": "fixed_cell"}
    surface = sample_frozen_mode_surface(
        plane, q1, q2,
        lambda x: (plane.project(x)[0] ** 2 - 1.0) ** 2
        + 2.0 * (plane.project(x)[1] - 0.5 * (1.0 - plane.project(x)[0] ** 2)) ** 2,
        **conditions,
    )
    assert surface.energies.shape == (3, 3)
    assert (surface.energy_unit, surface.potential_kind, surface.boundary_condition) == ("eV", "E", "fixed_cell")
    assert surface.energies[1, 1] == pytest.approx(1.5)
    assert surface.energies[2, 1] == pytest.approx(1.0)
    assert surface.energies[1, 0] == pytest.approx(0.0)
    with pytest.raises(ValueError, match="strictly increasing"):
        sample_frozen_mode_surface(plane, [0.0, 0.0], q2, lambda _: 0.0, **conditions)
    with pytest.raises(ValueError, match="non-finite"):
        sample_frozen_mode_surface(plane, q1, q2, lambda _: np.nan, **conditions)
    with pytest.raises(ValueError, match="boundary_condition"):
        sample_frozen_mode_surface(plane, q1, q2, lambda _: 0.0,
                                   energy_unit="eV", potential_kind="E", boundary_condition="")


def test_conditional_relaxation_fixes_composite_axes_and_lowers_frozen_energy() -> None:
    plane = _identity_plane()
    q = np.array([0.5, 0.2])

    def energy_gradient(x: np.ndarray) -> tuple[float, np.ndarray]:
        q1, q2 = plane.project(x)
        orthogonal = (x[0] - x[1]) / np.sqrt(2.0)
        preferred = 0.3 * q1 * q2
        energy = (q1**2 - 1.0)**2 + 2.0 * (q2 - 0.5 * (1.0 - q1**2))**2
        energy += 3.0 * (orthogonal - preferred)**2
        dq1 = 4.0 * q1 * (q1**2 - 1.0) + 4.0 * (q2 - 0.5 * (1.0 - q1**2)) * q1
        dq1 -= 6.0 * (orthogonal - preferred) * 0.3 * q2
        dq2 = 4.0 * (q2 - 0.5 * (1.0 - q1**2))
        dq2 -= 6.0 * (orthogonal - preferred) * 0.3 * q1
        dz = 6.0 * (orthogonal - preferred)
        gradient = dq1 * plane.axis_vectors[:, 0] + dq2 * plane.axis_vectors[:, 1]
        gradient += dz * np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
        return float(energy), gradient

    frozen = plane.frozen_coordinates(q)
    frozen_energy, _ = energy_gradient(frozen)
    minimum = relax_orthogonal_at_q(plane, q, energy_gradient, gradient_tolerance=1e-9)
    assert np.allclose(plane.project(minimum.coordinates), q, atol=1e-12)
    assert minimum.energy < frozen_energy
    assert minimum.energy == pytest.approx(frozen_energy - 3.0 * (0.3 * q[0] * q[1]) ** 2, abs=1e-12)
    assert minimum.orthogonal_gradient_norm < 1e-9
    assert minimum.n_starts == 1


def test_conditional_relaxation_respects_nonidentity_metric_and_physical_gradient() -> None:
    plane = ModePlane(
        reference=np.zeros(3), basis=np.diag([0.5, 1.0, 1.0]),
        metric_weights=np.array([4.0, 1.0, 1.0]),
        axis_weights=np.array([
            [1.0 / np.sqrt(2.0), 0.0],
            [1.0 / np.sqrt(2.0), 0.0],
            [0.0, 1.0],
        ]),
        axis_labels=("mass-weighted combined", "third mode"),
        amplitude_unit="toy sqrt(mass)*length", reference_id="toy/mass-metric",
    )

    def value_and_gradient(x: np.ndarray) -> tuple[float, np.ndarray]:
        a1, a2, q2 = np.array([2.0 * x[0], x[1], x[2]])
        q1 = (a1 + a2) / np.sqrt(2.0)
        orthogonal = (a1 - a2) / np.sqrt(2.0)
        error = orthogonal - 0.3 * q1
        value = q1**2 + q2**2 + 2.0 * error**2
        d_q1 = 2.0 * q1 - 1.2 * error
        d_orthogonal = 4.0 * error
        d_a1 = (d_q1 + d_orthogonal) / np.sqrt(2.0)
        d_a2 = (d_q1 - d_orthogonal) / np.sqrt(2.0)
        return float(value), np.array([2.0 * d_a1, d_a2, 2.0 * q2])

    result = relax_orthogonal_at_q(
        plane, [0.8, -0.4], value_and_gradient, gradient_tolerance=1e-10,
    )
    assert np.allclose(plane.project(result.coordinates), [0.8, -0.4], atol=1e-12)
    assert result.energy == pytest.approx(0.8**2 + 0.4**2, abs=1e-12)
    assert result.orthogonal_gradient_norm < 1e-10
    guarded = relax_orthogonal_at_q(
        plane, [0.8, -0.4], value_and_gradient,
        optimizer="safeguarded_bfgs", trial_validator=lambda _: True,
        gradient_tolerance=1e-9, initial_trust_radius=0.2,
    )
    assert guarded.energy == pytest.approx(result.energy, abs=1e-10)
    assert guarded.orthogonal_gradient_norm < 1e-9


def test_conditional_warm_continuation_uses_only_audited_start() -> None:
    plane = _identity_plane()
    q = np.array([0.5, 0.2])
    direction = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    start = plane.frozen_coordinates(q) + 0.4 * direction
    calls: list[np.ndarray] = []

    def value_and_gradient(x: np.ndarray) -> tuple[float, np.ndarray]:
        calls.append(x.copy())
        displacement = x - (plane.frozen_coordinates(q) + 0.2 * direction)
        return float(displacement @ displacement), 2.0 * displacement

    result = relax_orthogonal_at_q(
        plane, q, value_and_gradient, starts=[start],
        include_frozen_start=False, gradient_tolerance=1e-10,
    )
    assert result.n_starts == 1
    assert result.selected_start == 0
    assert np.allclose(calls[0], start, atol=1e-12)
    assert np.allclose(result.coordinates, plane.frozen_coordinates(q) + 0.2 * direction, atol=1e-10)
    with pytest.raises(ValueError, match="requires at least one supplied start"):
        relax_orthogonal_at_q(plane, q, value_and_gradient, include_frozen_start=False)


def test_conditional_point_preserves_all_branch_outcomes() -> None:
    plane = _identity_plane()
    q = np.array([0.6, 0.0])
    base = plane.frozen_coordinates(q)
    direction = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)

    def double_well(x: np.ndarray) -> tuple[float, np.ndarray]:
        z = float(direction @ (x - base))
        return (z * z - 1.0)**2, 4.0 * z * (z * z - 1.0) * direction

    point = relax_orthogonal_at_q(
        plane, q, double_well, starts=[base + 0.8 * direction, base - 0.8 * direction],
        gradient_tolerance=1e-8,
    )
    outcomes = point.branch_outcomes
    assert [row.start_index for row in outcomes] == [0, 1, 2]
    assert all(row.converged for row in outcomes)
    assert outcomes[0].energy == pytest.approx(1.0)
    assert outcomes[1].energy == pytest.approx(0.0, abs=1e-12)
    assert outcomes[2].energy == pytest.approx(0.0, abs=1e-12)
    assert direction @ (outcomes[1].coordinates - base) > 0.0
    assert direction @ (outcomes[2].coordinates - base) < 0.0
    assert point.selected_start in (1, 2)
    assert point.energy == pytest.approx(0.0, abs=1e-12)


def test_orthogonal_branch_seeds_preserve_fixed_axes_and_reject_wrong_gauge() -> None:
    plane = _identity_plane()
    center = np.array([0.3, 0.1, -0.4])
    direction = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    plus, minus = orthogonal_branch_seeds(plane, center, direction, amplitude=0.25)
    assert np.allclose(plus - center, 0.25 * direction)
    assert np.allclose(minus - center, -0.25 * direction)
    assert np.allclose(plane.project([plus, minus]), plane.project(center))
    with pytest.raises(ValueError, match="unit-normalized inside"):
        orthogonal_branch_seeds(plane, center, plane.axis_vectors[:, 0], amplitude=0.25)
    with pytest.raises(ValueError, match="finite and positive"):
        orthogonal_branch_seeds(plane, center, direction, amplitude=0.0)
    with pytest.raises(ValueError, match="exceeds the orthogonal amplitude bound"):
        orthogonal_branch_seeds(
            plane, center, direction, amplitude=0.25,
            orthogonal_amplitude_bound=0.1,
        )


def test_safeguarded_bfgs_backtracks_before_invalid_geometry_evaluation() -> None:
    plane = _identity_plane()
    q = np.array([0.3, -0.1])
    base = plane.frozen_coordinates(q)
    direction = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    evaluated: list[float] = []
    rejected: list[float] = []

    def coordinate(x: np.ndarray) -> float:
        return float(np.dot(direction, x - base))

    def guard(x: np.ndarray) -> bool:
        value = coordinate(x)
        if abs(value) > 1.4:
            rejected.append(value)
            return False
        return True

    def energy_gradient(x: np.ndarray) -> tuple[float, np.ndarray]:
        value = coordinate(x)
        evaluated.append(value)
        assert abs(value) <= 1.4
        return 10.0 * (value**2 - 1.0)**2, 40.0 * value * (value**2 - 1.0) * direction

    result = relax_orthogonal_at_q(
        plane, q, energy_gradient, starts=[base + 0.2 * direction],
        include_frozen_start=False, optimizer="safeguarded_bfgs",
        trial_validator=guard, initial_trust_radius=2.0,
        gradient_tolerance=1e-8, max_iterations=80,
    )
    assert rejected
    assert all(abs(value) <= 1.4 for value in evaluated)
    assert coordinate(result.coordinates) == pytest.approx(1.0, abs=1e-8)
    assert result.energy == pytest.approx(0.0, abs=1e-10)


def test_conditional_surface_continuity_screens_off_plane_branch_jump() -> None:
    plane = ModePlane(
        reference=np.zeros(3), basis=np.eye(3), metric_weights=np.ones(3),
        axis_weights=np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]),
        axis_labels=("Q1", "Q2"), amplitude_unit="toy length",
        reference_id="toy/continuity",
    )

    def energy_gradient(x: np.ndarray) -> tuple[float, np.ndarray]:
        target = 0.2 * x[0] + 0.1 * x[1]
        error = x[2] - target
        return error**2, np.array([-0.4 * error, -0.2 * error, 2.0 * error])

    smooth = sample_conditional_mode_surface(
        plane, [-1.0, 0.0, 1.0], [-0.5, 0.5], energy_gradient,
        energy_unit="toy", potential_kind="E", boundary_condition="fixed Q1,Q2",
    )
    audit = audit_conditional_surface_continuity(
        plane, smooth, maximum_allowed_orthogonal_jump=0.25,
    )
    assert audit.within_bound
    assert np.allclose(audit.q1_neighbor_jumps, 0.2, atol=1e-8)
    assert np.allclose(audit.q2_neighbor_jumps, 0.1, atol=1e-8)
    coordinates = smooth.coordinates.copy()
    coordinates[1, 1, 2] += 1.0
    selected = smooth.selected_starts.copy()
    selected[1, 1] = 1
    jump = ConditionalModeSurface(
        q1=smooth.q1, q2=smooth.q2, energies=smooth.energies,
        axis_labels=smooth.axis_labels, amplitude_unit=smooth.amplitude_unit,
        reference_id=smooth.reference_id, energy_unit=smooth.energy_unit,
        potential_kind=smooth.potential_kind, boundary_condition=smooth.boundary_condition,
        orthogonal_gradient_norms=smooth.orthogonal_gradient_norms,
        n_evaluations=smooth.n_evaluations, selected_starts=selected,
        coordinates=coordinates,
    )
    flagged = audit_conditional_surface_continuity(
        plane, jump, maximum_allowed_orthogonal_jump=0.25,
    )
    assert not flagged.within_bound
    assert flagged.maximum_orthogonal_jump > 0.7
    assert flagged.selected_start_switches > 0


def test_conditional_relaxation_reports_nonconvergence() -> None:
    plane = _identity_plane()
    with pytest.raises(RuntimeError, match="orthogonal-gradient tolerance"):
        relax_orthogonal_at_q(
            plane, [0.0, 0.0],
            lambda x: (float(x[0]), np.array([1.0, -1.0, 0.0])),
            gradient_tolerance=1e-10, max_iterations=2,
        )


def test_optimizer_component_tolerance_cannot_fake_vector_convergence() -> None:
    plane = ModePlane(
        reference=np.zeros(4), basis=np.eye(4), metric_weights=np.ones(4),
        axis_weights=np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0], [0.0, 0.0]]),
        axis_labels=("Q1", "Q2"), amplitude_unit="toy length",
        reference_id="toy/l2-versus-linf",
    )

    def energy_gradient(x: np.ndarray) -> tuple[float, np.ndarray]:
        delta = x[2:] - 0.01
        return float(np.dot(delta, delta)), np.r_[0.0, 0.0, 2.0 * delta]

    # At the frozen start, each gradient component is 0.02 (<0.025) but
    # the full vector norm is sqrt(2)*0.02 (>0.025). L-BFGS-B must not stop.
    result = relax_orthogonal_at_q(
        plane, [0.0, 0.0], energy_gradient, gradient_tolerance=0.025,
    )
    assert result.n_evaluations > 2
    assert result.orthogonal_gradient_norm < 0.025
    assert np.allclose(result.coordinates[2:], [0.01, 0.01], atol=1e-10)


def test_orthogonal_curvature_rejects_stationary_saddle_and_accepts_local_wells() -> None:
    plane = _identity_plane()
    orthogonal = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    calls: list[np.ndarray] = []

    def double_well(x: np.ndarray) -> tuple[float, np.ndarray]:
        calls.append(x.copy())
        r = float(orthogonal @ x)
        return r**4 - r**2, (4.0 * r**3 - 2.0 * r) * orthogonal

    saddle = audit_orthogonal_curvature(
        plane, np.zeros(3), double_well, step_amplitude=1e-3,
    )
    assert saddle.minimum_eigenvalue == pytest.approx(-2.0, abs=1e-5)
    assert saddle.n_gradient_evaluations == 2
    assert all(np.allclose(plane.project(x), [0.0, 0.0]) for x in calls)
    well = audit_orthogonal_curvature(
        plane, orthogonal / np.sqrt(2.0), double_well,
        step_amplitude=1e-3,
    )
    assert well.minimum_eigenvalue == pytest.approx(4.0, abs=1e-5)
    with pytest.raises(ValueError, match="step_amplitude"):
        audit_orthogonal_curvature(plane, np.zeros(3), double_well, step_amplitude=0.0)


def test_two_step_energy_curvature_checks_only_sampled_axes() -> None:
    # The first axis is quartic; the second is harmonic. The mixed term makes
    # the Hessian indefinite despite positive diagonal energy curvatures.
    def energy(x: float, y: float) -> float:
        return x*x + y*y - 3.0*x*y + 0.5*x**4

    def pairs(step: float) -> np.ndarray:
        return np.array([[energy(step, 0), energy(-step, 0)],
                         [energy(0, step), energy(0, -step)]])

    audit = audit_paired_energy_curvature(0.0, pairs(0.1), 0.1, pairs(0.2), 0.2)
    assert audit.curvature_by_step.shape == (2, 2)
    assert np.allclose(audit.curvature_by_step[:, 1], 2.0)
    assert audit.curvature_by_step[0, 0] == pytest.approx(2.01)
    assert audit.curvature_by_step[1, 0] == pytest.approx(2.04)
    assert audit.maximum_step_disagreement == pytest.approx(0.03)
    assert audit.minimum_sampled_curvature > 0
    assert np.linalg.eigvalsh(np.array([[2.0, -3.0], [-3.0, 2.0]])).min() < 0
    assert not audit.curvature_by_step.flags.writeable
    with pytest.raises(ValueError, match="must differ"):
        audit_paired_energy_curvature(0.0, pairs(0.1), 0.1, pairs(0.1), 0.1)
    with pytest.raises(ValueError, match="matching"):
        audit_paired_energy_curvature(0.0, pairs(0.1), 0.1, pairs(0.2)[:1], 0.2)


def test_two_step_energy_gradient_consistency_reports_anharmonic_difference() -> None:
    stiffness, quartic = 2.0, 3.0
    steps = (0.1, 0.2)
    energies = np.array([[0.5 * stiffness * h**2 + 0.25 * quartic * h**4] * 2
                         for h in steps])
    gradients = np.array([[stiffness * h + quartic * h**3,
                           -stiffness * h - quartic * h**3] for h in steps])
    audit = audit_directional_curvature_consistency(0.0, energies, gradients, steps)
    assert np.allclose(audit.energy_curvatures,
                       [stiffness + 0.5 * quartic * h**2 for h in steps])
    assert np.allclose(audit.gradient_curvatures,
                       [stiffness + quartic * h**2 for h in steps])
    assert np.allclose(audit.absolute_energy_gradient_disagreement,
                       [0.5 * quartic * h**2 for h in steps])
    assert not audit.energy_curvatures.flags.writeable
    with pytest.raises(ValueError, match="directional-gradient pairs"):
        audit_directional_curvature_consistency(0.0, energies, gradients[:1], steps)


def test_curvature_energy_resolution_screen_is_explicit_and_not_a_certificate() -> None:
    steps = (0.05, 0.1)
    energy_curvature, gradient_curvature = 0.0034, 0.0022
    energies = np.array([[0.5 * energy_curvature * h*h] * 2 for h in steps])
    gradients = np.array([[gradient_curvature * h, -gradient_curvature * h]
                          for h in steps])
    curvature = audit_directional_curvature_consistency(0.0, energies, gradients, steps)
    observed_origin_span = 2.1e-5
    screen = screen_curvature_energy_resolution(curvature, observed_origin_span)
    assert np.allclose(screen.energy_double_differences_eV,
                       energy_curvature * np.square(steps))
    assert np.allclose(screen.energy_gradient_disagreement_scales_eV,
                       abs(energy_curvature - gradient_curvature) * np.square(steps))
    assert screen.origin_span_over_energy_signal[0] > 1.0
    assert screen.origin_span_over_energy_signal[1] < 1.0
    assert np.all(screen.origin_span_over_energy_gradient_disagreement > 1.0)
    assert not screen.energy_double_differences_eV.flags.writeable
    with pytest.raises(ValueError, match="nonnegative origin span"):
        screen_curvature_energy_resolution(curvature, -1.0)
    with pytest.raises(TypeError, match="directional energy-gradient audit"):
        screen_curvature_energy_resolution(None, observed_origin_span)


def test_conditional_relaxation_can_freeze_an_explicit_zero_mode_gauge() -> None:
    plane = ModePlane(
        reference=np.zeros(4), basis=np.eye(4), metric_weights=np.ones(4),
        axis_weights=np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0], [0.0, 0.0]]),
        axis_labels=("Q1", "Q2"), amplitude_unit="toy length",
        reference_id="toy/frozen-gauge",
    )

    def energy_gradient(x: np.ndarray) -> tuple[float, np.ndarray]:
        target = np.array([0.2, -0.1, 0.3, 2.0])
        difference = x - target
        return float(np.dot(difference, difference)), 2.0 * difference

    gauge = np.array([[0.0], [0.0], [0.0], [1.0]])
    free = relax_orthogonal_at_q(plane, [0.2, -0.1], energy_gradient)
    fixed = relax_orthogonal_at_q(
        plane, [0.2, -0.1], energy_gradient,
        starts=[np.array([0.2, -0.1, 0.9, 2.0])],
        frozen_directions=gauge, gradient_tolerance=1e-10,
        orthogonal_amplitude_bound=1.0,
    )
    assert free.energy == pytest.approx(0.0, abs=1e-12)
    assert fixed.energy == pytest.approx(4.0, abs=1e-12)
    assert np.allclose(fixed.coordinates, [0.2, -0.1, 0.3, 0.0], atol=1e-10)
    audit = audit_orthogonal_curvature(
        plane, fixed.coordinates, energy_gradient,
        step_amplitude=1e-3, frozen_directions=gauge,
    )
    assert audit.n_gradient_evaluations == 2
    assert audit.minimum_eigenvalue == pytest.approx(2.0, abs=1e-10)
    with pytest.raises(ValueError, match="metric-orthogonal"):
        relax_orthogonal_at_q(
            plane, [0.0, 0.0], energy_gradient,
            frozen_directions=np.array([[1.0], [0.0], [0.0], [0.0]]),
        )
    with pytest.raises(ValueError, match="linearly independent"):
        relax_orthogonal_at_q(
            plane, [0.0, 0.0], energy_gradient,
            frozen_directions=np.column_stack([gauge, gauge]),
        )
    with pytest.raises(ValueError, match="orthogonal_amplitude_bound"):
        relax_orthogonal_at_q(
            plane, [0.0, 0.0], energy_gradient,
            orthogonal_amplitude_bound=0.0,
        )


def test_analytic_barrier_hierarchy_is_not_a_single_mode_contribution() -> None:
    """A coupled double well separates frozen, strict-plane and relaxed barriers."""
    plane = ModePlane(
        reference=np.zeros(3), basis=np.eye(3), metric_weights=np.ones(3),
        axis_weights=np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]),
        axis_labels=("primary", "secondary"), amplitude_unit="toy length",
        reference_id="analytic/coupled-double-well",
    )
    k_secondary, k_orthogonal = 2.0, 3.0
    alpha, beta = 0.4, 0.3

    def value_and_gradient(x: np.ndarray) -> tuple[float, np.ndarray]:
        primary, secondary, orthogonal = x
        envelope = 1.0 - primary**2
        s = secondary - alpha * envelope
        o = orthogonal - beta * envelope
        value = envelope**2 + k_secondary * s**2 + k_orthogonal * o**2
        gradient = np.array([
            -4.0 * primary * envelope + 4.0 * primary * (k_secondary * alpha * s + k_orthogonal * beta * o),
            2.0 * k_secondary * s,
            2.0 * k_orthogonal * o,
        ])
        return float(value), gradient

    assert value_and_gradient(np.array([-1.0, 0.0, 0.0]))[0] == pytest.approx(0.0)
    assert value_and_gradient(np.array([1.0, 0.0, 0.0]))[0] == pytest.approx(0.0)
    # Any connecting path crosses primary=0. These three constrained minima
    # at that section are exact barriers for their respective path spaces.
    frozen_line = value_and_gradient(np.array([0.0, 0.0, 0.0]))[0]
    strict_plane = value_and_gradient(np.array([0.0, alpha, 0.0]))[0]
    conditional = relax_orthogonal_at_q(
        plane, [0.0, alpha], value_and_gradient, gradient_tolerance=1e-10,
    )
    assert frozen_line == pytest.approx(1.0 + k_secondary * alpha**2 + k_orthogonal * beta**2)
    assert strict_plane == pytest.approx(1.0 + k_orthogonal * beta**2)
    assert conditional.energy == pytest.approx(1.0)
    assert np.allclose(conditional.coordinates, [0.0, alpha, beta], atol=1e-10)
    assert frozen_line > strict_plane > conditional.energy
    grid = sample_conditional_mode_surface(
        plane, [-1.0, 0.0, 1.0], [0.0, alpha], value_and_gradient,
        energy_unit="toy", potential_kind="E", boundary_condition="all_axes_free_except_Q1_Q2",
        gradient_tolerance=1e-10,
    )
    assert grid.energies.shape == (2, 3)
    assert grid.energies[1, 1] == pytest.approx(1.0)
    assert grid.orthogonal_gradient_norms.max() < 1e-10
    assert np.allclose(plane.project(grid.coordinates[1, 1]), [0.0, alpha])
    assert np.all(grid.n_evaluations >= 1)
    frozen_grid = sample_frozen_mode_surface(
        plane, grid.q1, grid.q2, lambda x: value_and_gradient(x)[0],
        energy_unit="toy", potential_kind="E", boundary_condition="all_axes_frozen_except_Q1_Q2",
    )
    assert np.all(grid.energies <= frozen_grid.energies + 1e-12)
    assert frozen_grid.energies[1, 1] == pytest.approx(strict_plane)
