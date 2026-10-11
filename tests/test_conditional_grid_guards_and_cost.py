"""Guarded grid dispatch and actual callback costs, not material predictions."""
import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT

from vcneb import (CalculatorModeEvaluator, ClampedModeCoordinates, ConditionalModePoint,
                   ModePlane, clamped_plane_vcneb_boundary, relax_orthogonal_at_q,
                   sample_conditional_mode_surface)

LABELS = dict(energy_unit="eV", potential_kind="E", boundary_condition="declared-fixed-chart")


def plane(size=3):
    return ModePlane(reference=np.zeros(size), basis=np.eye(size), metric_weights=np.ones(size),
        axis_weights=np.eye(size)[:, :2], axis_labels=("Q1", "Q2"),
        amplitude_unit="Angstrom", reference_id="analytic-not-HfO2")


def test_grid_counts_every_supplied_start_not_only_the_selected_one():
    chart = plane(2)
    calls = []
    def evaluate(x):
        calls.append(x.copy())
        return float(x @ x), 2*x
    result = sample_conditional_mode_surface(chart, [-.1, .1], [-.1, .1], evaluate,
        starts=[np.zeros(2), np.ones(2)], **LABELS)
    assert len(calls) == 12
    assert int(result.n_evaluations.sum()) == len(calls)
    np.testing.assert_array_equal(result.n_evaluations, np.full((2, 2), 3))


def test_grid_can_forward_guard_before_any_prescribed_static_call():
    called = []
    with pytest.raises(ValueError, match="geometry/bound guard"):
        sample_conditional_mode_surface(plane(2), [-.1, .1], [-.1, .1],
            lambda x: called.append(x), optimizer="safeguarded_bfgs",
            trial_validator=lambda x: False, initial_trust_radius=.03, **LABELS)
    assert called == []


def test_nonlinear_multi_branch_grid_preserves_signed_starts_and_total_calls():
    calls = []
    def evaluate(x):
        calls.append(x.copy())
        energy = x[0]**2 + x[1]**2 + (x[2]**2-1)**2
        return energy, np.array([2*x[0], 2*x[1], 4*x[2]*(x[2]**2-1)])
    result = sample_conditional_mode_surface(plane(), [-.1, .1], [-.2, .2], evaluate,
        starts=[np.array([0., 0., .7]), np.array([0., 0., -.7])],
        include_frozen_start=False, optimizer="safeguarded_bfgs", initial_trust_radius=.1,
        trial_validator=lambda x: np.max(np.abs(x)) < 1.5,
        gradient_tolerance=1e-7, max_iterations=80, **LABELS)
    np.testing.assert_allclose(result.energies, .05, atol=1e-12, rtol=0)
    np.testing.assert_allclose(np.abs(result.coordinates[..., 2]), 1., atol=1e-7, rtol=0)
    assert int(result.n_evaluations.sum()) == len(calls)
    assert np.any(np.array(calls)[:, 2] < 0) and np.any(np.array(calls)[:, 2] > 0)


@pytest.mark.parametrize("options", [
    dict(optimizer="unknown"),
    dict(optimizer="safeguarded_bfgs"),
    dict(include_frozen_start=False),
    dict(initial_trust_radius=0.),
])
def test_invalid_grid_options_fail_before_energy_callback(options):
    calls = []
    with pytest.raises(ValueError):
        sample_conditional_mode_surface(plane(), [-.1, .1], [-.1, .1],
            lambda x: calls.append(x), **options, **LABELS)
    assert calls == []


def test_point_total_keeps_unconverged_start_cost_and_legacy_selected_count():
    calls = []
    def evaluate(x):
        calls.append(x.copy())
        return float(x @ x), 2*x
    result = relax_orthogonal_at_q(plane(), [0., 0.], evaluate,
        starts=[np.array([0., 0., 5.])], max_iterations=1,
        gradient_tolerance=1e-10, optimizer="safeguarded_bfgs", initial_trust_radius=.01,
        trial_validator=lambda x: True)
    assert [row.converged for row in result.branch_outcomes] == [True, False]
    assert result.n_evaluations == result.branch_outcomes[result.selected_start].n_evaluations
    assert result.total_n_evaluations == len(calls)
    assert result.total_n_evaluations > result.n_evaluations


def test_legacy_point_with_missing_multi_start_evidence_has_unknown_total():
    point = ConditionalModePoint(q=np.zeros(2), coordinates=np.zeros(3), energy=0.,
        orthogonal_gradient_norm=0., n_evaluations=3, n_starts=2, selected_start=0)
    assert point.total_n_evaluations is None
    single = ConditionalModePoint(q=np.zeros(2), coordinates=np.zeros(3), energy=0.,
        orthogonal_gradient_norm=0., n_evaluations=3, n_starts=1, selected_start=0)
    assert single.total_n_evaluations == 3


def test_guarded_grid_reuses_clamped_joint_calculator_flow_and_exact_cache(tmp_path):
    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    boundary = clamped_plane_vcneb_boundary(len(atoms), atoms.cell.array, allow_tilt=False)
    coordinates = ClampedModeCoordinates(atoms, boundary=boundary, cell_scale_A=3.6)
    calls = []
    def factory(index, candidate, directory):
        boundary.validate_images([candidate])
        calls.append(index)
        return EMT()
    evaluate = CalculatorModeEvaluator(coordinates, factory, workdir=tmp_path/"grid",
        calculator_id="Cu-EMT-grid-not-HfO2-G3")
    axes = np.zeros((coordinates.coordinate_count, 2))
    axes[[0, 3], 0] = np.array([1., -1.])/np.sqrt(2)
    axes[[1, 4], 1] = np.array([1., -1.])/np.sqrt(2)
    modes = ModePlane(reference=np.zeros(coordinates.coordinate_count),
        basis=np.eye(coordinates.coordinate_count), metric_weights=np.ones(coordinates.coordinate_count),
        axis_weights=axes, axis_labels=("atomic Q1", "atomic Q2"), amplitude_unit="Angstrom",
        reference_id="Cu-clamped-grid-not-HfO2")
    accepted = []
    def guard(y):
        try:
            candidate = coordinates.to_atoms(y)
            boundary.validate_images([candidate])
            accepted.append(y.copy())
            return True
        except (ValueError, np.linalg.LinAlgError):
            return False
    callback_calls = []
    def counted(y):
        callback_calls.append(y.copy())
        return evaluate(y)
    result = sample_conditional_mode_surface(modes, [-.002, .002], [-.003, .003], counted,
        frozen_directions=coordinates.rigid_translation_directions(), optimizer="safeguarded_bfgs",
        trial_validator=guard, initial_trust_radius=.03, gradient_tolerance=2e-4,
        max_iterations=60, orthogonal_amplitude_bound=.3, **LABELS)
    assert int(result.n_evaluations.sum()) == len(callback_calls)
    assert 0 < len(calls) == evaluate.n_new_evaluations < len(callback_calls)
    for j, q2 in enumerate(result.q2):
        for i, q1 in enumerate(result.q1):
            y = result.coordinates[j, i]
            np.testing.assert_allclose(modes.project(y), [q1, q2], atol=1e-12, rtol=0)
            np.testing.assert_allclose(coordinates.to_atoms(y).cell[:2], atoms.cell[:2], atol=1e-13, rtol=0)
    assert len(accepted) >= 4
