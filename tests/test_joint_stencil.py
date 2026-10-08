"""Measured actions retain out-of-plane coupling and the real mechanical chart."""

import json
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms

from scripts.relax_clamped_ase_endpoint import load_seed
from vcneb import (
    ActiveJointCurvatureCoordinates, assemble_joint_directional_curvature,
    clamped_plane_vcneb_boundary, joint_curvature_probes,
)
from vcneb.joint_curvature import central_difference_hessian
from vcneb.joint_curvature import JointCurvatureCoordinates as LegacyCoordinates
from vcneb.mode_surface import audit_directional_curvature_consistency


def _harmonic(basis, hessian, h=.01, offset=None):
    offset = np.zeros(len(hessian)) if offset is None else offset
    plus = np.array([offset + hessian @ (h * column) for column in basis.T])
    minus = np.array([offset - hessian @ (h * column) for column in basis.T])
    return assemble_joint_directional_curvature(basis, plus, minus, step_A=h)


def test_selected_directions_retain_coupling_instead_of_claiming_full_hessian():
    matrix = np.array([[2., .3, .7], [.3, -1., .5], [.7, .5, 4.]])
    basis = np.eye(3)[:, :2]
    report = _harmonic(basis, matrix, offset=np.array([3., -1., 2.]))
    np.testing.assert_allclose(report.full_hessian_action, matrix @ basis, atol=1e-13)
    np.testing.assert_allclose(report.raw_projected, matrix[:2, :2], atol=1e-13)
    np.testing.assert_allclose(report.transverse_action[2], [.7, .5], atol=1e-13)
    assert report.full_hessian_action.shape == (3, 2)
    assert report.reciprocity_relative_defect < 1e-13
    # H[2,2] was not sampled: any value gives the same two-column evidence.
    changed = matrix.copy()
    changed[2, 2] = -10.
    np.testing.assert_array_equal(_harmonic(basis, changed).full_hessian_action,
                                  _harmonic(basis, matrix).full_hessian_action)


def test_full_direction_set_matches_legacy_assembly_and_basis_rotation():
    rng = np.random.default_rng(82)
    factor = rng.normal(size=(7, 7))
    matrix = (factor + factor.T) / 2
    basis, _ = np.linalg.qr(rng.normal(size=(7, 7)))
    report = _harmonic(basis, matrix)
    np.testing.assert_allclose(report.symmetric_projected, basis.T @ matrix @ basis, atol=1e-12)
    np.testing.assert_allclose(report.transverse_action, 0, atol=1e-12)
    plus = .01 * matrix.T
    legacy, defect = central_difference_hessian(plus, -plus, .01)
    identity = _harmonic(np.eye(7), matrix)
    np.testing.assert_array_equal(identity.symmetric_projected, legacy)
    assert identity.reciprocity_relative_defect == defect
    subspace = basis[:, :3]
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    first, second = _harmonic(subspace, matrix), _harmonic(subspace @ rotation, matrix)
    np.testing.assert_allclose(second.symmetric_projected,
                               rotation.T @ first.symmetric_projected @ rotation, atol=1e-12)
    np.testing.assert_allclose(second.transverse_action, first.transverse_action @ rotation, atol=1e-12)


def test_reciprocity_failure_is_exposed_and_arrays_do_not_alias_inputs():
    basis = np.eye(3)[:, :2]
    plus = np.array([[.02, .001, .004], [0., .01, 0.]])
    minus = -plus
    report = assemble_joint_directional_curvature(basis, plus, minus, step_A=.01)
    assert report.reciprocity_relative_defect > 0
    assert report.symmetrization_operator_change_eV_A2 > 0
    assert report.raw_projected[1, 0] == pytest.approx(.1)
    assert report.symmetric_projected[1, 0] == pytest.approx(.05)
    basis[:] = 0
    plus[:] = 0
    np.testing.assert_array_equal(report.directions, np.eye(3)[:, :2])
    for value in (report.directions, report.full_hessian_action, report.raw_projected,
                  report.symmetric_projected, report.transverse_action):
        assert not value.flags.writeable


@pytest.mark.parametrize("basis", [np.empty((3, 0)), np.ones((3, 2)),
                                  np.array([1., 0., 0.]), np.full((3, 2), np.nan)])
def test_invalid_direction_basis_is_rejected_without_normalizing(basis):
    with pytest.raises(ValueError, match="orthonormal"):
        assemble_joint_directional_curvature(basis, np.zeros((2, 3)), np.zeros((2, 3)), step_A=.01)


@pytest.mark.parametrize("step", [0., -1., np.nan, np.inf, True])
def test_invalid_probe_step_is_rejected(step):
    with pytest.raises(ValueError, match="step_A"):
        assemble_joint_directional_curvature(np.eye(3), np.eye(3), -np.eye(3), step_A=step)


@pytest.mark.parametrize("fault", ["square_subspace", "missing_pair", "nonfinite"])
def test_incomplete_or_wrong_metric_gradients_do_not_become_curvature(fault):
    plus, minus = np.zeros((2, 3)), np.zeros((2, 3))
    if fault == "square_subspace":
        plus = minus = np.zeros((2, 2))
    elif fault == "missing_pair":
        minus = minus[:1]
    else:
        plus[0, 0] = np.nan
    with pytest.raises(ValueError, match="full chart"):
        assemble_joint_directional_curvature(np.eye(3)[:, :2], plus, minus, step_A=.01)


def test_finite_inputs_that_overflow_differences_are_not_reported_as_curvature():
    with pytest.raises(ValueError, match="numerically representable"):
        assemble_joint_directional_curvature(np.eye(2), np.eye(2) * 1e308,
                                              -np.eye(2) * 1e308, step_A=.01)


def test_probe_order_center_and_guard_preserve_reference_and_call_no_calculator():
    reference = bulk("Cu", "fcc", a=3.6, cubic=True)
    reference.calc = EMT()
    original = reference.positions.copy()
    boundary = clamped_plane_vcneb_boundary(len(reference), reference.cell.array, allow_tilt=True)
    chart = ActiveJointCurvatureCoordinates.for_clamped_plane(reference, cell_scale_A=3.6, boundary=boundary)
    basis = chart.translation_free_basis()[:, :3]
    center = np.ones(chart.size) * .001
    checked = []
    probes = joint_curvature_probes(chart, basis, step_A=.005, center_delta_A=center,
                                   candidate_validator=lambda atoms: checked.extend(atoms))
    assert len(probes) == 6 and len(checked) == 7
    assert [(point.direction_index, point.sign) for point in probes] == [(i, s) for i in range(3) for s in (1, -1)]
    for point in probes:
        np.testing.assert_allclose(point.delta_A, center + point.sign * .005 * basis[:, point.direction_index])
        np.testing.assert_array_equal(point.atoms.positions, chart.displaced(point.delta_A).positions)
        assert point.atoms.calc is None and not point.delta_A.flags.writeable
        boundary.validate_images([point.atoms])
    np.testing.assert_array_equal(reference.positions, original)
    assert reference.calc.results == {}
    with pytest.raises(ValueError, match="center_delta"):
        joint_curvature_probes(chart, basis, step_A=.005, center_delta_A=np.zeros(2))
    with pytest.raises(ValueError, match="guard refused"):
        joint_curvature_probes(chart, basis, step_A=.005,
                               candidate_validator=lambda _: (_ for _ in ()).throw(ValueError("guard refused")))
    assert reference.calc.results == {}


@pytest.mark.parametrize("mutation", ["geometry", "calculator", "constraints"])
def test_guard_cannot_repair_a_probe_or_make_it_non_inert(mutation):
    reference = bulk("Cu", "fcc", a=3.6, cubic=True)
    chart = ActiveJointCurvatureCoordinates(reference, 3.6)

    def bad_guard(points):
        if mutation == "geometry":
            points[1].positions[0, 0] += .001
        elif mutation == "calculator":
            points[1].calc = EMT()
        else:
            points[1].set_constraint(FixAtoms(indices=[0]))

    before = reference.positions.copy()
    with pytest.raises(ValueError, match="reject only"):
        joint_curvature_probes(chart, chart.translation_free_basis()[:, :2], step_A=.005,
                               candidate_validator=bad_guard)
    np.testing.assert_array_equal(reference.positions, before)


def test_legacy_chart_constraints_and_nonfinite_positions_are_refused():
    reference = bulk("Cu", "fcc", a=3.6, cubic=True)
    chart = LegacyCoordinates(reference, 3.6)
    reference.set_constraint(FixAtoms(indices=[0]))
    with pytest.raises(ValueError, match="explicit atomic"):
        joint_curvature_probes(chart, chart.translation_free_basis()[:, :1], step_A=.005)
    reference.set_constraint()
    reference.positions[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        joint_curvature_probes(chart, chart.translation_free_basis()[:, :1], step_A=.005)


@pytest.mark.parametrize("tilt", [False, True])
@pytest.mark.parametrize("rotated", [False, True])
def test_EMT_joint_hessian_energy_gradient_check_at_two_steps(tilt, rotated):
    reference = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = reference.cell.array.copy()
    cell[0] += [0., .10, .15]
    cell[1] += [.05, 0., .08]
    cell[2] += [.11, -.07, 0.]
    if rotated:
        rotation, _ = np.linalg.qr(np.random.default_rng(81).normal(size=(3, 3)))
        cell = cell @ rotation.T
    reference.set_cell(cell, scale_atoms=True)
    boundary = clamped_plane_vcneb_boundary(len(reference), cell, allow_tilt=tilt)
    chart = ActiveJointCurvatureCoordinates.for_clamped_plane(reference, cell_scale_A=3.6, boundary=boundary)
    basis = chart.translation_free_basis()
    center = np.zeros(chart.size)
    center[0] = .003
    center[-1] = -.004
    origin = chart.displaced(center)
    origin.calc = EMT()
    pressure = .03  # Cu analytic implementation check, not an HfO2 setting
    energy0 = origin.get_potential_energy() + pressure * origin.get_volume()
    pairs, matrices = [], []
    steps = (.001, .002)
    for h in steps:
        gradients, energies = [], []
        for probe in joint_curvature_probes(chart, basis, step_A=h, center_delta_A=center):
            atoms = probe.atoms
            atoms.calc = EMT()
            gradients.append(chart.enthalpy_gradient(atoms, atoms.get_forces(), atoms.get_stress(voigt=False), pressure))
            energies.append(atoms.get_potential_energy() + pressure * atoms.get_volume())
        gradients = np.array(gradients).reshape(basis.shape[1], 2, chart.size)
        pairs.append((np.array(energies).reshape(-1, 2), gradients))
        matrices.append(assemble_joint_directional_curvature(basis, gradients[:, 0], gradients[:, 1], step_A=h))
    for index in range(basis.shape[1]):
        audit = audit_directional_curvature_consistency(
            energy0, np.array([pair[0][index] for pair in pairs]),
            np.array([pair[1][index] @ basis[:, index] for pair in pairs]), steps,
        )
        np.testing.assert_allclose(audit.gradient_curvatures,
                                   [matrix.raw_projected[index, index] for matrix in matrices], atol=1e-12)
        assert audit.absolute_energy_gradient_disagreement.max() < .002
    assert max(matrix.reciprocity_relative_defect for matrix in matrices) < 1e-4
    assert np.linalg.norm(matrices[0].symmetric_projected - matrices[1].symmetric_projected, ord=2) < .005


def test_all_ten_hafnia_starters_stage_72_translation_free_joint_probes_without_DFT():
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds"
    records = json.loads((root / "clamped_seed_manifest.json").read_text())["seeds"]
    for item in records:
        atoms, boundary, record = load_seed(root / item["manifest_file"])
        chart = ActiveJointCurvatureCoordinates.for_clamped_plane(atoms, cell_scale_A=record["cell_scale_A"], boundary=boundary)
        probes = joint_curvature_probes(chart, chart.translation_free_basis(), step_A=.005)
        assert len(probes) == 72
        for probe in probes:
            boundary.validate_images([probe.atoms])
            assert probe.atoms.calc is None and probe.atoms.get_chemical_symbols() == atoms.get_chemical_symbols()
            assert np.linalg.norm(probe.delta_A) == pytest.approx(.005)


def test_standalone_receipt_separates_implementation_checks_from_material_claims():
    from scripts.verify_joint_curvature_stencil import verify

    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds"
    result = verify(root)
    assert result["total_inert_HfO2_geometry_probes"] == 720
    assert result["EMT_evaluations"] == 180 and len(result["EMT_two_step_checks"]) == 4
    assert result["new_DFT_calls"] == 0
    assert not result["saddle_certified"] and not result["material_joint_Hessian_measured"]
    assert not result["material_barrier_prediction_validated"]
    assert result["synthetic_unsampled_direction_self_curvature_not_measured"]
    assert all(item["status"] == "geometry_seed_not_relaxed" for item in result["geometry_starters"])
    assert all(len(digest) == 64 for digest in result["module_sha256"].values())
    json.dumps(result)


def test_standalone_verifier_refuses_existing_receipt_before_work(tmp_path, monkeypatch):
    import scripts.verify_joint_curvature_stencil as checker

    target = tmp_path / "existing.json"
    target.write_text("original")
    monkeypatch.setattr("sys.argv", ["verify", "--report", str(target)])
    monkeypatch.setattr(checker, "verify", lambda _: (_ for _ in ()).throw(RuntimeError("must not run")))
    with pytest.raises(FileExistsError, match="existing"):
        checker.main()
    assert target.read_text() == "original"
