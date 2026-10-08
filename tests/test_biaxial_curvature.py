"""Controlled strain derivatives are not confused with open-cell relaxation."""

from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator

from scripts.relax_clamped_ase_endpoint import load_seed
from vcneb import (ActiveJointCurvatureCoordinates, BiaxialClampedCurvatureCoordinates,
                   assemble_joint_directional_curvature, clamped_plane_vcneb_boundary,
                   joint_curvature_probes)
from vcneb.quadratic_reduction import condition_quadratic_energy


def chart_fixture(*, allow_tilt=True, rotated=False, anchor=.01, scale=4.2):
    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = atoms.cell.array.copy()
    cell[0] += [0., .10, .15]
    cell[1] += [.05, 0., .08]
    cell[2] += [.11, -.07, 0.]
    if rotated:
        rotation, _ = np.linalg.qr(np.random.default_rng(81).normal(size=(3, 3)))
        cell = cell @ rotation.T
    atoms.set_cell(cell, scale_atoms=True)
    boundary = clamped_plane_vcneb_boundary(len(atoms), cell, allow_tilt=allow_tilt)
    chart = BiaxialClampedCurvatureCoordinates(
        atoms, boundary=boundary, cell_scale_A=3.6,
        anchor_strain=anchor, strain_scale_A=scale,
    )
    return atoms, boundary, chart


@pytest.mark.parametrize("allow_tilt", [False, True])
@pytest.mark.parametrize("rotated", [False, True])
@pytest.mark.parametrize("deformed", [False, True])
def test_all_gradients_match_physical_energy_at_general_point(allow_tilt, rotated, deformed):
    reference, _, chart = chart_fixture(allow_tilt=allow_tilt, rotated=rotated)
    delta = np.zeros(chart.size)
    if deformed:
        delta = np.random.default_rng(17).normal(size=chart.size) * .007
    atoms = chart.displaced(delta)
    atoms.calc = EMT()
    pressure = .03  # only an implementation check, not HfO2 production pressure
    analytic = chart.enthalpy_gradient(delta, atoms, atoms.get_forces(),
                                      atoms.get_stress(voigt=False), pressure)
    numerical = []
    h = 1e-5
    for i in range(chart.size):
        pair = []
        for sign in (1, -1):
            offset = delta.copy()
            offset[i] += sign * h
            probe = chart.displaced(offset)
            eps = chart.anchor_strain + offset[-1] / chart.strain_scale_A
            np.testing.assert_allclose(probe.cell[:2], reference.cell[:2]
                                       * (1+eps)/(1+chart.anchor_strain), atol=1e-13)
            probe.calc = EMT()
            pair.append(probe.get_potential_energy() + pressure * probe.get_volume())
        numerical.append((pair[0]-pair[1])/(2*h))
    np.testing.assert_allclose(analytic, numerical, atol=3e-5, rtol=0)


def test_zero_control_coordinate_reproduces_existing_internal_chart():
    reference, boundary, chart = chart_fixture(rotated=True)
    original = ActiveJointCurvatureCoordinates.for_clamped_plane(
        reference, cell_scale_A=3.6, boundary=boundary,
    )
    delta = np.random.default_rng(21).normal(size=chart.size) * .01
    delta[-1] = 0.
    expected, current = original.displaced(delta[:-1]), chart.displaced(delta)
    np.testing.assert_allclose(current.cell.array, expected.cell.array, atol=1e-14)
    np.testing.assert_allclose(current.positions, expected.positions, atol=1e-14)
    current.calc = EMT()
    forces, stress = current.get_forces(), current.get_stress(voigt=False)
    np.testing.assert_allclose(chart.enthalpy_gradient(delta, current, forces, stress, 0.)[:-1],
                               original.enthalpy_gradient(expected, forces, stress, 0.), atol=1e-12)


def test_substrate_reaction_is_not_an_open_residual_but_drives_external_work():
    reference, boundary, chart = chart_fixture()
    delta = np.zeros(chart.size)
    sigma = (np.eye(3) - np.outer(boundary.normal, boundary.normal)) * .04
    np.testing.assert_allclose(boundary.open_traction(sigma), 0, atol=1e-15)
    gradient = chart.enthalpy_gradient(delta, reference, np.zeros((len(reference), 3)), sigma, 0.)
    np.testing.assert_allclose(gradient[:-1], 0, atol=1e-15)
    assert gradient[-1] * chart.strain_scale_A == pytest.approx(
        2 * reference.get_volume() * .04 / (1+chart.anchor_strain), abs=1e-12,
    )


def test_scale_changes_coordinate_not_physical_strain_derivative():
    ref, _, a = chart_fixture(scale=3.)
    _, _, b = chart_fixture(scale=6.)
    qa, qb = np.zeros(a.size), np.zeros(b.size)
    qa[-1], qb[-1] = 3.*.002, 6.*.002
    aa, ab = a.displaced(qa), b.displaced(qb)
    np.testing.assert_allclose(aa.cell.array, ab.cell.array, atol=1e-14)
    aa.calc = EMT()
    force, stress = aa.get_forces(), aa.get_stress(voigt=False)
    ga, gb = a.enthalpy_gradient(qa, aa, force, stress, 0.), b.enthalpy_gradient(qb, ab, force, stress, 0.)
    np.testing.assert_allclose(ga[:-1], gb[:-1], atol=1e-13)
    assert ga[-1]*a.strain_scale_A == pytest.approx(gb[-1]*b.strain_scale_A, abs=1e-12)


def test_release_basis_excludes_control_parameter_and_removes_translations():
    reference, _, chart = chart_fixture()
    released = chart.internal_relaxation_basis()
    retained = chart.controlled_direction()[:, None]
    assert released.shape == (chart.size, chart.internal_size-3)
    combined = np.column_stack((retained, released))
    np.testing.assert_allclose(combined.T @ combined, np.eye(chart.size-3), atol=1e-12)
    np.testing.assert_array_equal(released[-1], 0.)
    translations = np.zeros((chart.size, 3))
    translations[:3*len(reference)] = np.tile(np.eye(3), (len(reference), 1))
    np.testing.assert_allclose(translations.T @ released, 0, atol=1e-12)
    # A resolved synthetic quadratic tests the existing release machinery.
    hessian = 2.*np.eye(chart.size) + .1*np.ones((chart.size, chart.size))
    model = condition_quadratic_energy(hessian, np.ones(chart.size)*.03, 0.,
                                      retained, released, stability_floor=.01)
    displacement = model.full_displacement(np.array([.02]))
    assert displacement[-1] == pytest.approx(.02)
    np.testing.assert_allclose(released.T @ (np.ones(chart.size)*.03 + hessian@displacement), 0, atol=1e-12)


def test_reference_copy_drops_cache_and_cannot_be_changed_by_caller():
    reference, _, chart = chart_fixture()
    baseline = chart.displaced(np.zeros(chart.size))
    reference.calc = SinglePointCalculator(reference, energy=100.)
    reference.positions += 1.
    reference.cell[0] *= 1.01
    repeated = chart.displaced(np.zeros(chart.size))
    assert repeated.calc is None
    np.testing.assert_array_equal(repeated.positions, baseline.positions)
    np.testing.assert_array_equal(repeated.cell.array, baseline.cell.array)


@pytest.mark.parametrize("fault", ["cell", "lift", "species", "periodicity"])
def test_wrong_declared_geometry_is_rejected_without_repair(fault):
    reference, _, chart = chart_fixture()
    delta = np.zeros(chart.size)
    wrong = reference.copy()
    if fault == "cell":
        wrong.cell[0, 0] += .001
    elif fault == "lift":
        wrong.positions[0] += wrong.cell[0]
    elif fault == "species":
        wrong.symbols[0] = "Au"
    else:
        wrong.pbc[0] = False
    before = wrong.positions.copy()
    with pytest.raises(ValueError):
        chart.enthalpy_gradient(delta, wrong, np.zeros((len(wrong), 3)), np.zeros((3, 3)), 0.)
    np.testing.assert_array_equal(wrong.positions, before)


@pytest.mark.parametrize("anchor,scale", [(True, 4.), (-1., 4.), (np.nan, 4.),
                                        (0., 0.), (0., np.inf), (0., True),
                                        ("0.0", 4.), (0., "4.0")])
def test_invalid_control_definitions_are_rejected(anchor, scale):
    with pytest.raises(ValueError):
        chart_fixture(anchor=anchor, scale=scale)


def test_prescribed_internal_mixed_curvature_retains_full_gradient_and_reciprocity():
    _, _, chart = chart_fixture(rotated=True)
    directions = np.zeros((chart.size, 2))
    directions[:, 0] = chart.controlled_direction()
    directions[0, 1], directions[3, 1] = 1/np.sqrt(2), -1/np.sqrt(2)
    center = np.random.default_rng(9).normal(size=chart.size)*.005
    projections = []
    for step in (1e-4, 2e-4):
        probes = joint_curvature_probes(chart, directions, step_A=step, center_delta_A=center)
        assert len(probes) == 4 and all(p.atoms.calc is None for p in probes)
        gradients = []
        for probe in probes:
            probe.atoms.calc = EMT()
            gradients.append(chart.enthalpy_gradient(
                probe.delta_A, probe.atoms, probe.atoms.get_forces(),
                probe.atoms.get_stress(voigt=False), .03,
            ))
        block = assemble_joint_directional_curvature(
            directions, np.array(gradients)[::2], np.array(gradients)[1::2], step_A=step,
        )
        assert block.full_hessian_action.shape == (chart.size, 2)
        assert block.reciprocity_relative_defect < 1e-6
        assert np.linalg.norm(block.transverse_action) > 1e-3
        projections.append(block.raw_projected)
    np.testing.assert_allclose(projections[0], projections[1], atol=2e-5, rtol=0)


def test_collapse_bad_arrays_and_nonsymmetric_stress_are_rejected():
    reference, _, chart = chart_fixture()
    for delta in (np.zeros(chart.size-1), np.full(chart.size, np.nan),
                  np.r_[np.zeros(chart.internal_size), -5.*chart.strain_scale_A]):
        with pytest.raises(ValueError):
            chart.displaced(delta)
    stress = np.eye(3)
    stress[0, 1] = .1
    with pytest.raises(ValueError, match="symmetric"):
        chart.enthalpy_gradient(np.zeros(chart.size), reference,
                                np.zeros((len(reference), 3)), stress, 0.)


def test_actual_hafnia_training_starters_share_the_registered_controlled_plane():
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds"
    for condition, epsilon in (("strain_0000", 0.), ("strain_p0100", .01)):
        for phase in ("T", "PO_plus", "PO_minus_T_preserving", "PO_minus_T_reversing", "M"):
            atoms, boundary, seed = load_seed(root/condition/phase/"endpoint_seed.json")
            chart = BiaxialClampedCurvatureCoordinates(
                atoms, boundary=boundary, cell_scale_A=seed["cell_scale_A"],
                anchor_strain=epsilon, strain_scale_A=seed["cell_scale_A"],
            )
            assert chart.size == 40 and chart.internal_size == 39
            point = np.zeros(chart.size)
            point[-1] = .0001*chart.strain_scale_A
            changed = chart.displaced(point)
            np.testing.assert_allclose(changed.cell[:2], atoms.cell[:2]
                                       * (1+epsilon+.0001)/(1+epsilon), atol=1e-13)
            assert changed.calc is None
            assert changed.get_chemical_symbols() == atoms.get_chemical_symbols()
