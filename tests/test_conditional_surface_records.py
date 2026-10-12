"""Grid branch evidence survives selection and prospective registration."""
from dataclasses import replace
import json
from pathlib import Path
import runpy

import numpy as np
import pytest

from vcneb import (
    ConditionalHoldout, ConditionalPointEvidence, ConditionalProbeRecord,
    ConditionalProbeSpec, ModePlane, audit_conditional_probe_coverage,
    register_conditional_probe_plan, sample_conditional_mode_surface,
)
from vcneb.conditional_coverage import _snapshot


def fixture(n=2, *, starts=()):
    plane = ModePlane(
        reference=np.zeros(3), basis=np.eye(3), metric_weights=np.ones(3),
        axis_weights=np.eye(3)[:, :2], axis_labels=("Q1", "Q2"),
        amplitude_unit="toy", reference_id="synthetic-records-not-material",
    )
    calls = []

    def energy_gradient(x):
        calls.append(x.copy())
        return float(x[0] + 2*x[1] + (x[2]**2-1)**2), np.array([
            1., 2., 4*x[2]*(x[2]**2-1),
        ])

    surface = sample_conditional_mode_surface(
        plane, np.linspace(0., 1., n), np.linspace(0., 1., n), energy_gradient,
        starts=starts, gradient_tolerance=1e-8,
        energy_unit="toy", potential_kind="E", boundary_condition="fixed Q; toy only",
    )
    return plane, surface, calls


def test_sampler_keeps_all_point_records_and_signed_competing_starts():
    plane, surface, calls = fixture(starts=(np.array([0., 0., 1.]), np.array([0., 0., -1.])))
    assert len(getattr(surface, "point_records", ())) == 2
    assert sum(surface.n_evaluations.flat) == len(calls)
    for j, row in enumerate(surface.point_records):
        assert len(row) == 2
        for i, point in enumerate(row):
            np.testing.assert_array_equal(point.q, [surface.q1[i], surface.q2[j]])
            np.testing.assert_array_equal(point.coordinates, surface.coordinates[j, i])
            assert point.energy == surface.energies[j, i]
            assert point.total_n_evaluations == surface.n_evaluations[j, i]
            assert point.n_starts == 3 and len(point.branch_outcomes) == 3
            assert point.selected_start in (1, 2)
            assert point.branch_outcomes[1].coordinates[2] > 0
            assert point.branch_outcomes[2].coordinates[2] < 0
            # Stationary frozen start is retained, not certified as a minimum.
            assert point.branch_outcomes[0].converged
            assert point.branch_outcomes[0].energy > point.energy
            np.testing.assert_allclose(plane.project(point.coordinates), point.q)


def test_grid_keeps_unconverged_start_and_its_callbacks():
    plane, _, _ = fixture()
    calls = []

    def evaluate(x):
        calls.append(x.copy())
        return float(x @ x), 2*x

    grid = sample_conditional_mode_surface(
        plane, [0., 1.], [0., 1.], evaluate, starts=[np.array([0., 0., 5.])],
        optimizer="safeguarded_bfgs", trial_validator=lambda x: True,
        initial_trust_radius=.01, gradient_tolerance=1e-10, max_iterations=1,
        energy_unit="toy", potential_kind="E", boundary_condition="toy",
    )
    assert int(grid.n_evaluations.sum()) == len(calls)
    for row in grid.point_records:
        for point in row:
            assert [item.converged for item in point.branch_outcomes] == [True, False]
            assert point.selected_start == 0
            assert point.total_n_evaluations > point.n_evaluations
            assert point.branch_outcomes[1].orthogonal_gradient_norm > 1e-10


def with_first(surface, point):
    rows = [list(row) for row in surface.point_records]
    rows[0][0] = point
    return replace(surface, point_records=rows)


def test_records_are_defensively_copied_and_nested_arrays_read_only():
    _, surface, _ = fixture()
    point = surface.point_records[0][0]
    supplied = replace(point, q=point.q.copy(), coordinates=point.coordinates.copy(),
                       branch_outcomes=tuple(replace(row, coordinates=row.coordinates.copy())
                                             for row in point.branch_outcomes))
    owned = with_first(surface, supplied).point_records[0][0]
    supplied.q[:] = 99
    supplied.coordinates[:] = 99
    supplied.branch_outcomes[0].coordinates[:] = 99
    np.testing.assert_array_equal(owned.q, [0., 0.])
    np.testing.assert_array_equal(owned.coordinates, [0., 0., 0.])
    for array in (owned.q, owned.coordinates, owned.branch_outcomes[0].coordinates):
        with pytest.raises(ValueError):
            array[0] = 99


@pytest.mark.parametrize("change", [
    {"q": np.array([.01, 0.])}, {"coordinates": np.array([.1, 0., 0.])},
    {"energy": .1}, {"orthogonal_gradient_norm": .1},
    {"n_evaluations": 999}, {"n_starts": 2}, {"selected_start": 9},
    {"n_evaluations": True}, {"energy": np.nan},
])
def test_inconsistent_point_record_is_rejected(change):
    _, surface, _ = fixture()
    with pytest.raises(ValueError, match="record"):
        with_first(surface, replace(surface.point_records[0][0], **change))


@pytest.mark.parametrize("change", [
    {"n_evaluations": 0}, {"start_index": True}, {"start_index": 2},
    {"converged": "yes"}, {"converged": False}, {"energy": np.inf},
    {"coordinates": np.array([0., 0., np.nan])}, {"orthogonal_gradient_norm": -1.},
])
def test_invalid_or_incomplete_branch_record_is_rejected(change):
    _, surface, _ = fixture()
    point = surface.point_records[0][0]
    row = replace(point.branch_outcomes[0], **change)
    with pytest.raises(ValueError, match="record"):
        with_first(surface, replace(point, branch_outcomes=(row,)))


@pytest.mark.parametrize("records", [((None, None),), ((None,), (None,)), ((None, None), (None, None))])
def test_wrong_shape_or_type_is_rejected(records):
    _, surface, _ = fixture()
    with pytest.raises(ValueError, match="record"):
        replace(surface, point_records=records)


def test_existing_legacy_plan_signatures_remain_byte_semantically_identical():
    example = runpy.run_path(str(Path(__file__).with_name("test_conditional_coverage.py")))["example"]
    _, surface, _, _, _, plan, _ = example(2)
    assert surface.point_records == ()
    assert plan.training_snapshot_sha256 == "0c1c8d8ab56a623e424357202cf6833a38ace08c1c5254bc8285f046332f1c0c"
    assert plan.to_dict()["registration_sha256"] == "20fa5ea0beb468b7a76ce5a8ae47a9b55aca952a1e91def5416c119d7f34a41f"


def registered_fixture():
    plane, surface, _ = fixture(3, starts=[np.array([0., 0., 1.]), np.array([0., 0., -1.])])
    evidence = ConditionalPointEvidence(
        calculator_id="toy", energy_reference_id="toy-zero", raw_audit_sha256="a"*64,
        branch_audit_sha256="b"*64, curvature_audit_sha256="c"*64,
        curvature_resolution_sha256="d"*64, minimum_orthogonal_curvature=8., curvature_uncertainty=.1,
    )
    rows = [[evidence]*3 for _ in range(3)]
    plan = register_conditional_probe_plan(
        plane, surface, rows, [ConditionalProbeSpec("probe", 1, 1, (.75, .75))],
        training_evaluation_ids=[f"toy-{i}" for i in range(9)], protocol_id="toy",
        start_policy_id="toy-signed", selection_rule_id="toy-last-cell", probe_labels_unread=True,
        energy_zero_id="toy-zero", gradient_tolerance=.001, maximum_neighbor_orthogonal_jump=.1,
        maximum_holdout_energy_error=.001, maximum_holdout_coordinate_error=.001,
    )
    observations = {"probe": ConditionalProbeRecord("toy-new", point=ConditionalHoldout(
        cell_i=1, cell_j=1, q=(.75, .75), energy=2.25,
        coordinates=np.array([.75, .75, 1.]), evidence=evidence, orthogonal_gradient_norm=0.,
    ))}
    return plane, surface, rows, plan, observations


def test_complete_records_serialize_and_survive_internal_probe_patch_slicing():
    plane, surface, rows, plan, observations = registered_fixture()
    json.dumps(plan.to_dict(), allow_nan=False)
    report = audit_conditional_probe_coverage(plane, surface, rows, plan, observations)
    assert report.n_passed == 1 and report.n_cells == 4
    assert not report.whole_grid_screen.ready_for_exploratory_interpolation
    assert not report.to_dict()["continuous_region_error_bound_established"]


@pytest.mark.parametrize("remove", [False, True])
def test_removing_or_changing_unselected_record_invalidates_registration(remove):
    plane, surface, rows, plan, observations = registered_fixture()
    if remove:
        altered = replace(surface, point_records=())
    else:
        point = surface.point_records[0][0]
        outcomes = list(point.branch_outcomes)
        outcomes[0] = replace(outcomes[0], energy=outcomes[0].energy+.25)
        altered = with_first(surface, replace(point, branch_outcomes=tuple(outcomes)))
    assert _snapshot(plane, altered, rows) != plan.training_snapshot_sha256
    with pytest.raises(ValueError, match="changed after"):
        audit_conditional_probe_coverage(plane, altered, rows, plan, observations)


def test_fully_prescribed_static_points_retain_all_starts_without_extra_calls():
    plane = ModePlane(reference=np.zeros(2), basis=np.eye(2), metric_weights=np.ones(2),
                      axis_weights=np.eye(2), axis_labels=("Q1", "Q2"),
                      amplitude_unit="toy", reference_id="toy-static")
    calls = []

    def evaluate(x):
        calls.append(x.copy())
        return float(x @ x), 2*x

    surface = sample_conditional_mode_surface(
        plane, [0., 1.], [0., 1.], evaluate, starts=[np.ones(2), np.zeros(2)],
        energy_unit="toy", potential_kind="E", boundary_condition="toy-static",
    )
    assert len(calls) == 12 == int(surface.n_evaluations.sum())
    assert all(len(point.branch_outcomes) == 3 for row in surface.point_records for point in row)
