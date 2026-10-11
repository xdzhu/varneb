"""Finite independent probes must not acquire whole-grid/area guarantees."""
from dataclasses import replace
import json

import numpy as np
import pytest

from vcneb import (
    ConditionalHoldout, ConditionalModeSurface, ConditionalPointEvidence,
    ConditionalProbePlan, ConditionalProbeRecord, ConditionalProbeSpec, ModePlane,
    audit_conditional_probe_coverage, register_conditional_probe_plan,
)


def example(n=9, positions=None):
    # Non-unit mass metric makes the coordinate error's units explicit.
    plane = ModePlane(
        reference=np.zeros(3), basis=np.diag([0.5, 1., 1./3.]),
        metric_weights=np.array([4., 1., 9.]),
        axis_weights=np.eye(3)[:, :2], axis_labels=("Q1", "Q2"),
        amplitude_unit="sqrt(amu)*angstrom", reference_id="synthetic-common-reference",
    )
    axes = np.linspace(0., 1., n)
    x, y = np.meshgrid(axes, axes)
    coordinates = np.stack([x/2., y, 0.03*x + 0.02*y], axis=-1)
    surface = ConditionalModeSurface(
        q1=axes, q2=axes, energies=2*x + 3*y + 0.5*x*y,
        axis_labels=plane.axis_labels, amplitude_unit=plane.amplitude_unit,
        reference_id=plane.reference_id, energy_unit="eV/cell", potential_kind="E",
        boundary_condition="synthetic fixed-Q local branch; not material data",
        orthogonal_gradient_norms=np.zeros((n, n)), n_evaluations=np.ones((n, n), dtype=int),
        selected_starts=np.zeros((n, n), dtype=int), coordinates=coordinates,
    )
    evidence = ConditionalPointEvidence(
        calculator_id="synthetic-only-not-HfO2", energy_reference_id="common-zero",
        raw_audit_sha256="a"*64, branch_audit_sha256="b"*64,
        curvature_audit_sha256="c"*64, curvature_resolution_sha256="d"*64,
        minimum_orthogonal_curvature=0.3, curvature_uncertainty=0.05,
    )
    grid_evidence = [[evidence for _ in axes] for _ in axes]
    if positions is None:
        positions = [(i, i, 0.5, 0.5) for i in range(n - 1)]
    probes = tuple(ConditionalProbeSpec(
        f"p{k}", i, j, (axes[i] + tx*(axes[i + 1]-axes[i]), axes[j] + ty*(axes[j + 1]-axes[j]))
    ) for k, (i, j, tx, ty) in enumerate(positions))
    settings = dict(energy_zero_id="common-zero", gradient_tolerance=0.01,
                    maximum_neighbor_orthogonal_jump=0.1,
                    maximum_holdout_energy_error=0.02, maximum_holdout_coordinate_error=0.05)
    registration = dict(training_evaluation_ids=[f"training-{k}" for k in range(n*n)],
                        protocol_id="synthetic-protocol", start_policy_id="same-starts-v1",
                        selection_rule_id="declared-locations-not-labels", probe_labels_unread=True,
                        **settings)
    plan = register_conditional_probe_plan(plane, surface, grid_evidence, probes, **registration)
    observations = {}
    for probe in probes:
        a, b = probe.q
        observations[probe.probe_id] = ConditionalProbeRecord(
            evaluation_id=f"test-{probe.probe_id}",
            point=ConditionalHoldout(
                cell_i=probe.cell_i, cell_j=probe.cell_j, q=probe.q,
                energy=2*a + 3*b + 0.5*a*b,
                coordinates=np.array([a/2., b, 0.03*a + 0.02*b]),
                evidence=evidence, orthogonal_gradient_norm=0.,
            ),
        )
    return plane, surface, grid_evidence, probes, registration, plan, observations


def test_eight_checks_of_64_cells_remain_finite_pointwise_evidence():
    plane, surface, evidence, _, _, plan, observations = example()
    report = audit_conditional_probe_coverage(plane, surface, evidence, plan, observations)
    assert (report.n_cells, report.n_registered, report.n_measured, report.n_passed) == (64, 8, 8, 8)
    assert len(report.cells_with_passed_point_checks) == 8
    assert report.maximum_observed_absolute_energy_error < 1e-14
    assert report.maximum_observed_metric_coordinate_error < 1e-14
    assert not report.whole_grid_screen.ready_for_exploratory_interpolation
    assert any("not every grid cell" in item for item in report.whole_grid_screen.blockers)
    data = json.loads(json.dumps(report.to_dict(), allow_nan=False))
    assert not data["continuous_region_error_bound_established"]
    assert not data["prospective_blinding_independently_proven"]
    assert data["new_calculator_calls"] == 0
    assert data["coordinate_error_unit"] == "sqrt(amu)*angstrom"


def test_one_cell_screen_can_pass_but_still_has_no_area_error_bound():
    plane, surface, evidence, _, _, plan, observations = example(2)
    report = audit_conditional_probe_coverage(plane, surface, evidence, plan, observations)
    assert report.n_passed == report.n_cells == 1
    assert report.whole_grid_screen.ready_for_exploratory_interpolation
    assert not report.to_dict()["continuous_region_error_bound_established"]


def test_multiple_probes_in_same_cell_do_not_inflate_cell_coverage():
    plane, surface, evidence, _, _, plan, observations = example(
        positions=[(0, 0, .2, .3), (0, 0, .7, .8)],
    )
    report = audit_conditional_probe_coverage(plane, surface, evidence, plan, observations)
    assert report.n_passed == 2
    assert report.cells_with_measurements == report.cells_with_passed_point_checks == ((0, 0),)


def test_signed_failure_metric_error_and_pending_points_are_not_dropped():
    plane, surface, evidence, _, _, plan, observations = example()
    point = observations["p0"].point
    coordinates = point.coordinates.copy()
    coordinates[2] += .04
    observations["p0"] = replace(observations["p0"], point=replace(
        point, energy=point.energy-.07, coordinates=coordinates,
    ))
    observations["p1"] = ConditionalProbeRecord("failed-DFT-id", failure="SCF did not converge")
    del observations["p2"]
    report = audit_conditional_probe_coverage(plane, surface, evidence, plan, observations)
    assert (report.n_passed, report.n_failed_checks, report.n_failed_calculations, report.n_pending) == (5, 1, 1, 1)
    assert report.n_measured == 6
    assert report.maximum_observed_absolute_energy_error == pytest.approx(.07)
    assert report.maximum_observed_metric_coordinate_error == pytest.approx(.12)
    assert report.point_results[0].measured_minus_predicted_energy == pytest.approx(-.07)
    assert report.point_results[1].status == "calculation_failed"
    assert report.point_results[1].measured_minus_predicted_energy is None
    assert report.point_results[2].status == "pending"
    assert (0, 0) not in report.cells_with_passed_point_checks


@pytest.mark.parametrize("change", ["energy", "coordinate", "reference", "metric", "evidence", "pressure_label"])
def test_post_registration_refit_or_contract_changes_are_rejected(change):
    plane, surface, evidence, _, _, plan, observations = example()
    if change == "energy":
        surface = replace(surface, energies=surface.energies + .001)
    elif change == "coordinate":
        c = surface.coordinates.copy(); c[:, :, 2] += .001
        surface = replace(surface, coordinates=c)
    elif change == "reference":
        plane = replace(plane, reference=np.array([0., 0., .01]))
    elif change == "metric":
        plane = replace(plane, metric_weights=plane.metric_weights*4, basis=plane.basis/2)
    elif change == "evidence":
        evidence[0][0] = replace(evidence[0][0], raw_audit_sha256="e"*64)
    else:
        surface = replace(surface, potential_kind="H", boundary_condition="a different declared pressure")
    with pytest.raises(ValueError, match="prediction changed"):
        audit_conditional_probe_coverage(plane, surface, evidence, plan, observations)


def test_manifest_roundtrip_is_exact_and_tampering_detected():
    *_, plan, _ = example()
    data = json.loads(json.dumps(plan.to_dict(), allow_nan=False))
    assert ConditionalProbePlan.from_dict(data).to_dict() == data
    data["maximum_holdout_energy_error"] *= 10
    with pytest.raises(ValueError, match="checksum mismatch"):
        ConditionalProbePlan.from_dict(data)


@pytest.mark.parametrize("field,value", [("probe_labels_unread", False), ("probe_labels_unread", 1),
                                        ("training_evaluation_ids", ["repeated"]*81),
                                        ("training_evaluation_ids", ["too-few"]),
                                        ("start_policy_id", ""), ("selection_rule_id", "")])
def test_registration_requires_blinding_attestation_ids_and_policies(field, value):
    plane, surface, evidence, probes, settings, _, _ = example()
    settings[field] = value
    with pytest.raises(ValueError):
        register_conditional_probe_plan(plane, surface, evidence, probes, **settings)


@pytest.mark.parametrize("probe", [ConditionalProbeSpec("bad", 8, 0, (1.1, .05)),
                                  ConditionalProbeSpec("bad", 0, 0, (0., .05)),
                                  ConditionalProbeSpec("bad", 0, 0, (.5, .05))])
def test_outside_or_boundary_locations_do_not_register(probe):
    plane, surface, evidence, _, settings, _, _ = example()
    with pytest.raises(ValueError, match="outside|strictly inside"):
        register_conditional_probe_plan(plane, surface, evidence, [probe], **settings)


@pytest.mark.parametrize("q", [(np.nan, .5), (1.,), (1., 2., 3.)])
def test_malformed_probe_coordinates_rejected(q):
    with pytest.raises(ValueError):
        ConditionalProbeSpec("bad", 0, 0, q)


@pytest.mark.parametrize("index", [True, 0.0, -1])
def test_probe_cell_indices_are_strict_integers(index):
    with pytest.raises(ValueError):
        ConditionalProbeSpec("bad", index, 0, (.5, .5))


def test_duplicate_probe_positions_and_identifiers_rejected():
    plane, surface, evidence, probes, settings, _, _ = example()
    for extra in (probes[0], replace(probes[0], probe_id="another")):
        with pytest.raises(ValueError, match="must be unique"):
            register_conditional_probe_plan(plane, surface, evidence, [probes[0], extra], **settings)


@pytest.mark.parametrize("reuse", ["training-0", "test-p1"])
def test_declared_training_and_test_calculation_reuse_rejected(reuse):
    plane, surface, evidence, _, _, plan, observations = example()
    observations["p0"] = replace(observations["p0"], evaluation_id=reuse)
    with pytest.raises(ValueError, match="repeats a training or test"):
        audit_conditional_probe_coverage(plane, surface, evidence, plan, observations)


def test_unregistered_or_moved_observation_rejected():
    plane, surface, evidence, _, _, plan, observations = example()
    observations["extra"] = observations["p0"]
    with pytest.raises(ValueError, match="unregistered"):
        audit_conditional_probe_coverage(plane, surface, evidence, plan, observations)
    del observations["extra"]
    observations["p0"] = replace(observations["p0"], point=replace(observations["p0"].point, q=(.07, .07)))
    with pytest.raises(ValueError, match="registered cell/Q"):
        audit_conditional_probe_coverage(plane, surface, evidence, plan, observations)


@pytest.mark.parametrize("bad_evidence", [dict(curvature_uncertainty=.4), dict(branch_audit_sha256=""),
                                        dict(energy_reference_id="other-zero"), dict(calculator_id="other-inputs")])
def test_existing_raw_curvature_branch_and_contract_gates_remain(bad_evidence):
    plane, surface, evidence, _, _, plan, observations = example()
    original = observations["p0"].point
    observations["p0"] = replace(observations["p0"], point=replace(
        original, evidence=replace(original.evidence, **bad_evidence),
    ))
    report = audit_conditional_probe_coverage(plane, surface, evidence, plan, observations)
    assert report.point_results[0].status == "point_checks_failed"
    assert report.n_failed_checks == 1
    assert not report.whole_grid_screen.ready_for_exploratory_interpolation


def test_no_observations_returns_pending_not_zero_error_or_all_passed():
    plane, surface, evidence, _, _, plan, _ = example()
    report = audit_conditional_probe_coverage(plane, surface, evidence, plan, {})
    assert report.n_pending == 8 and report.n_measured == report.n_passed == 0
    assert report.maximum_observed_absolute_energy_error is None
    assert report.maximum_observed_metric_coordinate_error is None


def test_record_without_measurement_or_failure_cannot_fake_a_result():
    with pytest.raises(ValueError, match="explicit calculation failure"):
        ConditionalProbeRecord("empty")
    *_, observations = example()
    with pytest.raises(ValueError, match="also be a failed calculation"):
        replace(observations["p0"], failure="failed")


def test_loaded_or_directly_constructed_plan_still_checks_locations_and_ids():
    plane, surface, evidence, _, _, plan, observations = example()
    for altered in (replace(plan, training_evaluation_ids=("one",)),
                    replace(plan, probes=(ConditionalProbeSpec("p0", 8, 0, (1.1, .05)),))):
        with pytest.raises(ValueError, match="per measured grid|outside"):
            audit_conditional_probe_coverage(plane, surface, evidence, altered, observations)
