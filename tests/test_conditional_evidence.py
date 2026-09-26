"""An audited conditional grid must not be conflated with a certified PES."""

from dataclasses import replace

import numpy as np
import pytest

from vcneb import (
    ConditionalHoldout,
    ConditionalModeSurface,
    ConditionalPointEvidence,
    ModePlane,
    screen_conditional_interpolation,
)


def _example():
    plane = ModePlane(
        reference=np.zeros(3), basis=np.eye(3), metric_weights=np.ones(3),
        axis_weights=np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]),
        axis_labels=("Q1", "Q2"), amplitude_unit="toy",
        reference_id="same-cubic-reference",
    )
    coordinates = np.array([
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.2]],
        [[0.0, 1.0, 0.1], [1.0, 1.0, 0.3]],
    ])
    surface = ConditionalModeSurface(
        q1=np.array([0.0, 1.0]), q2=np.array([0.0, 1.0]),
        energies=np.array([[0.0, 1.0], [2.0, 3.0]]),
        axis_labels=plane.axis_labels, amplitude_unit=plane.amplitude_unit,
        reference_id=plane.reference_id, energy_unit="eV", potential_kind="E",
        boundary_condition="Q1,Q2 fixed; other directions relaxed",
        orthogonal_gradient_norms=np.zeros((2, 2)),
        n_evaluations=np.ones((2, 2), dtype=int),
        selected_starts=np.zeros((2, 2), dtype=int), coordinates=coordinates,
    )
    evidence = ConditionalPointEvidence(
        calculator_id="same-input-and-pseudopotentials",
        energy_reference_id="same-cubic-energy-zero",
        raw_audit_sha256="a" * 64, branch_audit_sha256="b" * 64,
        curvature_audit_sha256="c" * 64,
        curvature_resolution_sha256="d" * 64,
        minimum_orthogonal_curvature=0.3, curvature_uncertainty=0.05,
    )
    point_evidence = [[evidence, evidence], [evidence, evidence]]
    holdout = ConditionalHoldout(
        cell_i=0, cell_j=0, q=(0.5, 0.5), energy=1.5,
        coordinates=np.array([0.5, 0.5, 0.15]), evidence=evidence,
        orthogonal_gradient_norm=0.0,
    )
    settings = dict(
        energy_zero_id="same-cubic-energy-zero",
        gradient_tolerance=0.01, maximum_neighbor_orthogonal_jump=0.25,
        maximum_holdout_energy_error=0.02,
        maximum_holdout_coordinate_error=0.05,
    )
    return plane, surface, point_evidence, holdout, settings


def test_measured_patch_passes_only_exploratory_interpolation_screen():
    plane, surface, evidence, holdout, settings = _example()
    screen = screen_conditional_interpolation(
        plane, surface, evidence, [holdout], **settings,
    )
    assert screen.ready_for_exploratory_interpolation
    assert screen.blockers == ()
    assert screen.n_cells == screen.n_independent_holdouts == 1
    assert screen.maximum_holdout_energy_error == pytest.approx(0.0)
    assert screen.maximum_holdout_coordinate_error == pytest.approx(0.0)


def test_missing_curvature_branch_audit_or_holdout_blocks_interpolation():
    plane, surface, evidence, holdout, settings = _example()
    uncertain = replace(
        evidence[1][0], minimum_orthogonal_curvature=0.03,
        curvature_uncertainty=0.05, branch_audit_sha256="",
    )
    evidence[1][0] = uncertain
    screen = screen_conditional_interpolation(
        plane, surface, evidence, [], **settings,
    )
    assert not screen.ready_for_exploratory_interpolation
    assert any("curvature" in reason for reason in screen.blockers)
    assert any("raw/branch/curvature/resolution audit" in reason for reason in screen.blockers)
    assert any("holdout" in reason for reason in screen.blockers)


def test_holdout_energy_or_branch_mismatch_blocks_interpolation():
    plane, surface, evidence, holdout, settings = _example()
    bad = replace(
        holdout, energy=1.6, coordinates=np.array([0.5, 0.5, 0.45]),
        evidence=replace(holdout.evidence, calculator_id="different-input"),
    )
    screen = screen_conditional_interpolation(
        plane, surface, evidence, [bad], **settings,
    )
    assert not screen.ready_for_exploratory_interpolation
    assert any("energy interpolation" in reason for reason in screen.blockers)
    assert any("branch-coordinate" in reason for reason in screen.blockers)
    assert any("calculator contracts" in reason for reason in screen.blockers)


def test_holdout_must_be_independent_interior_and_fixed_q():
    plane, surface, evidence, holdout, settings = _example()
    with pytest.raises(ValueError, match="strictly inside"):
        screen_conditional_interpolation(
            plane, surface, evidence, [replace(holdout, q=(0.0, 0.5))], **settings,
        )
    drifted = replace(holdout, coordinates=np.array([0.6, 0.5, 0.15]))
    screen = screen_conditional_interpolation(
        plane, surface, evidence, [drifted], **settings,
    )
    assert not screen.ready_for_exploratory_interpolation
    assert any("drifted from fixed Q" in reason for reason in screen.blockers)


def test_energy_zero_gradient_and_neighbor_jump_each_block_a_smooth_claim():
    plane, surface, evidence, holdout, settings = _example()
    evidence[0][0] = replace(evidence[0][0], energy_reference_id="other-zero")
    norms = surface.orthogonal_gradient_norms.copy()
    norms[0, 1] = 0.02
    coordinates = surface.coordinates.copy()
    coordinates[1, 1, 2] += 0.8
    altered = replace(
        surface, orthogonal_gradient_norms=norms, coordinates=coordinates,
    )
    screen = screen_conditional_interpolation(
        plane, altered, evidence, [holdout], **settings,
    )
    assert not screen.ready_for_exploratory_interpolation
    assert any("energy zero" in reason for reason in screen.blockers)
    assert any("orthogonal-gradient" in reason for reason in screen.blockers)
    assert any("branch jump" in reason for reason in screen.blockers)
