"""A recoverable calculator mode flow must preserve the actual substrate."""
import json
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms

from vcneb import (
    CalculatorModeEvaluator, ClampedModeCoordinates, ModePlane,
    clamped_plane_vcneb_boundary, relax_orthogonal_at_q,
)
from vcneb.active_curvature import ActiveJointCurvatureCoordinates
from scripts.analyze_hfo2_clamped_reference import terminal_images


def reference(rotated=False):
    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = atoms.cell.array.copy()
    cell += [[0., .12, .09], [.04, 0., .07], [.19, -.13, .06]]
    if rotated:
        rotation, _ = np.linalg.qr(np.random.default_rng(82).normal(size=(3, 3)))
        cell = cell @ rotation.T
    atoms.set_cell(cell, scale_atoms=True)
    return atoms


def make_chart(atoms, allow_tilt=True, scale=3.6):
    boundary = clamped_plane_vcneb_boundary(len(atoms), atoms.cell.array, allow_tilt=allow_tilt)
    return ClampedModeCoordinates(atoms, boundary=boundary, cell_scale_A=scale), boundary


@pytest.mark.parametrize("rotated", [False, True])
@pytest.mark.parametrize("allow_tilt", [False, True])
@pytest.mark.parametrize("pressure", [0., .03])
def test_recoverable_mode_gradient_matches_independent_energy_differences(tmp_path, rotated, allow_tilt, pressure):
    atoms = reference(rotated)
    chart, boundary = make_chart(atoms, allow_tilt)
    calls = []

    def factory(index, candidate, directory):
        boundary.validate_images([candidate])
        calls.append(index)
        return EMT()

    evaluator = CalculatorModeEvaluator(chart, factory, workdir=tmp_path / "mode", calculator_id="Cu-EMT-check",
        pressure_eV_per_A3=pressure)
    y = np.random.default_rng(83).normal(size=chart.coordinate_count)*.01
    energy, gradient = evaluator(y)
    h = 1e-5
    numerical = []
    # Independent direct ASE evaluations, not the adapter's gradient/energies.
    for i in range(chart.coordinate_count):
        energies = []
        for sign in (1., -1.):
            delta = y.copy(); delta[i] += sign*h
            candidate = chart.to_atoms(delta)
            candidate.calc = EMT()
            energies.append(candidate.get_potential_energy()+pressure*candidate.get_volume())
        numerical.append((energies[0]-energies[1])/(2*h))
    np.testing.assert_allclose(gradient, numerical, atol=3e-5, rtol=0)
    assert len(calls) == evaluator.n_new_evaluations == 1
    cached_energy, cached_gradient = evaluator(y)
    assert cached_energy == energy
    np.testing.assert_array_equal(cached_gradient, gradient)
    assert len(calls) == 1
    np.testing.assert_allclose(chart.to_atoms(y).cell[:2], atoms.cell[:2], atol=1e-13, rtol=0)
    assert evaluator.contract["mechanical_boundary"]["allow_tilt"] is allow_tilt
    assert evaluator.contract["mechanical_boundary"]["prescribed_substrate_relaxable"] is False


def test_contract_resume_preserves_cell_basis_scale_and_reference(tmp_path):
    atoms = reference()
    chart, boundary = make_chart(atoms)
    calls = []

    def factory(index, image, directory):
        calls.append(index); return EMT()

    work = tmp_path / "mode"
    evaluator = CalculatorModeEvaluator(chart, factory, workdir=work, calculator_id="Cu-EMT-check")
    y = np.zeros(chart.coordinate_count)
    original = evaluator(y)
    resumed = CalculatorModeEvaluator(chart, factory, workdir=work, calculator_id="Cu-EMT-check", resume=True)
    assert resumed(y)[0] == original[0] and calls == [0]
    np.testing.assert_array_equal(resumed(y)[1], original[1])
    for altered, _ in (make_chart(atoms, False), make_chart(atoms, True, scale=4.1)):
        with pytest.raises(ValueError, match="resume contract changed"):
            CalculatorModeEvaluator(altered, factory, workdir=work, calculator_id="Cu-EMT-check", resume=True)
    altered = reference(); altered.positions[0, 0] += .001
    with pytest.raises(ValueError, match="resume contract changed"):
        CalculatorModeEvaluator(make_chart(altered)[0], factory, workdir=work,
                                calculator_id="Cu-EMT-check", resume=True)
    assert calls == [0]
    # Callers cannot mutate chart state through original/reference/contract copies.
    cell = chart.reference_cell.copy()
    atoms.cell[0, 0] += .2
    returned = chart.reference_atoms; returned.cell[0, 0] += .1
    contract = chart.mechanical_contract; contract["deformation_basis"][0][0][0] += 1.
    np.testing.assert_array_equal(chart.reference_cell, cell)
    np.testing.assert_array_equal(chart.to_atoms(y).cell.array, cell)
    assert chart.reference_atoms.calc is None


@pytest.mark.parametrize("failure", ["cell", "positions", "order", "constraints", "pbc", "stress"])
def test_calculator_mutation_or_bad_stress_not_saved_as_valid_result(tmp_path, failure):
    atoms = reference()
    chart, _ = make_chart(atoms)

    def factory(index, candidate, directory):
        if failure == "cell": candidate.cell[0, 0] += .01
        if failure == "positions": candidate.positions[0, 0] += .01
        if failure == "order": candidate.numbers[0] = 8
        if failure == "constraints": candidate.set_constraint(FixAtoms(indices=[0]))
        if failure == "pbc": candidate.pbc[0] = False
        stress = np.eye(3)
        if failure == "stress": stress[0, 1] = .1
        return SinglePointCalculator(candidate, energy=1., forces=np.zeros((len(candidate), 3)), stress=stress)

    evaluator = CalculatorModeEvaluator(chart, factory, workdir=tmp_path / "mode", calculator_id="injected-invalid")
    with pytest.raises(ValueError):
        evaluator(np.zeros(chart.coordinate_count))
    assert not list((tmp_path / "mode").glob("eval-*/result.json"))
    assert len(list((tmp_path / "mode").glob("eval-*/failure.json"))) == 1
    assert evaluator.n_new_evaluations == 0


def test_invalid_trial_rejects_before_calculator_and_keeps_legacy_chart_guard(tmp_path):
    atoms = reference()
    chart, _ = make_chart(atoms)
    called = []
    evaluator = CalculatorModeEvaluator(chart, lambda *args: called.append(True),
        workdir=tmp_path / "invalid", calculator_id="never-called")
    with pytest.raises(ValueError): evaluator(np.full(chart.coordinate_count, np.nan))
    # Reverse the open-row height, leaving prescribed rows unchanged.
    y = np.zeros(chart.coordinate_count)
    y[-1] = -100.
    with pytest.raises(ValueError): evaluator(y)
    assert called == [] and not list((tmp_path / "invalid").glob("eval-*"))
    with pytest.raises(TypeError):
        CalculatorModeEvaluator(object(), lambda *args: EMT(), workdir=tmp_path / "unknown", calculator_id="unknown")


@pytest.mark.parametrize("scale", [True, 0., -1., np.nan])
def test_bad_chart_scale_is_not_an_implicit_metric(scale):
    with pytest.raises(ValueError, match="cell_scale_A"):
        make_chart(reference(), scale=scale)


def test_conditional_optimizer_runs_through_clamped_calculator_flow(tmp_path):
    atoms = reference()
    chart, boundary = make_chart(atoms)
    evaluator = CalculatorModeEvaluator(chart, lambda *args: EMT(), workdir=tmp_path / "conditional",
        calculator_id="Cu-EMT-conditional-check")
    axes = np.zeros((chart.coordinate_count, 2))
    axes[[0, 3], 0] = [1., -1.]/np.sqrt(2)
    axes[[1, 4], 1] = [1., -1.]/np.sqrt(2)
    plane = ModePlane(reference=np.zeros(chart.coordinate_count), basis=np.eye(chart.coordinate_count),
        metric_weights=np.ones(chart.coordinate_count), axis_weights=axes,
        axis_labels=("atomic Q1", "atomic Q2"), amplitude_unit="Angstrom",
        reference_id="Cu-check-not-HfO2-G3")
    q = np.array([.004, -.003])
    def trial_valid(values):
        try:
            chart.to_atoms(values)
            return True
        except (ValueError, np.linalg.LinAlgError):
            return False
    result = relax_orthogonal_at_q(plane, q, evaluator, frozen_directions=chart.rigid_translation_directions(),
        gradient_tolerance=1e-4, max_iterations=60, optimizer="safeguarded_bfgs",
        initial_trust_radius=.03, orthogonal_amplitude_bound=.3, trial_validator=trial_valid)
    np.testing.assert_allclose(plane.project(result.coordinates), q, atol=1e-12, rtol=0)
    boundary.validate_images([chart.to_atoms(result.coordinates)])
    assert result.orthogonal_gradient_norm <= 1e-4
    assert result.n_starts == len(result.branch_outcomes) == 1
    assert evaluator.n_new_evaluations > 1
    for path in (tmp_path / "conditional").glob("eval-*/result.json"):
        values = np.array(json.loads(path.read_text())["coordinates"])
        boundary.validate_images([chart.to_atoms(values)])
        np.testing.assert_allclose(plane.project(values), q, atol=1e-12, rtol=0)


@pytest.mark.parametrize("allow_tilt", [False, True])
def test_inverse_keeps_registered_long_periodic_lift_not_nearest_images(allow_tilt):
    chart, _ = make_chart(reference(True), allow_tilt)
    y = np.random.default_rng(89).normal(size=chart.coordinate_count)*.012
    # More than a full cell is a valid already registered displacement.
    y[:3] += 1.3*chart.reference_cell[0]
    image = chart.to_atoms(y)
    np.testing.assert_allclose(chart.from_atoms(image), y, atol=1e-12, rtol=0)
    wrapped = image.copy(); wrapped.wrap()
    assert not np.allclose(chart.from_atoms(wrapped), y)
    np.testing.assert_allclose(chart.to_atoms(chart.from_atoms(wrapped)).positions,
                               wrapped.positions, atol=1e-12, rtol=0)
    changed = image.copy(); changed.cell[0, 0] += .001
    with pytest.raises(ValueError, match="substrate"):
        chart.from_atoms(changed)


@pytest.mark.parametrize("case,digest", [
    ("clamped_T_terminal_E079_20261011", "c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e"),
    ("clamped_M_terminal_E068_20261011", "56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c"),
])
def test_actual_native_frames_use_same_active_chart_and_raw_energy(case, digest):
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"
    images, observation, _, _ = terminal_images(root / case, digest)
    boundary = clamped_plane_vcneb_boundary(len(images[0]),
        np.array(observation["mechanical_boundary"]["reference_cell_A"]), allow_tilt=True)
    common_reference = terminal_images(root / "clamped_T_terminal_E079_20261011",
        "c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e")[0][-1]
    common = ClampedModeCoordinates(common_reference, boundary=boundary, cell_scale_A=observation["cell_scale_A"])
    for native in images:
        chart = ClampedModeCoordinates(native, boundary=boundary, cell_scale_A=observation["cell_scale_A"])
        energy, gradient = chart.energy_and_gradient(np.zeros(chart.coordinate_count), native)
        legacy = ActiveJointCurvatureCoordinates.for_clamped_plane(native,
            boundary=boundary, cell_scale_A=observation["cell_scale_A"])
        expected = legacy.enthalpy_gradient(native, native.get_forces(), native.get_stress(voigt=False), 0.)
        np.testing.assert_array_equal(gradient, expected)
        assert energy == native.get_potential_energy()
        assert chart.coordinate_count == 39 and chart.cell_dofs == 3
        q = common.from_atoms(native)
        reconstructed = common.to_atoms(q)
        np.testing.assert_allclose(reconstructed.positions, native.positions, atol=1e-12, rtol=0)
        np.testing.assert_allclose(reconstructed.cell.array, native.cell.array, atol=1e-12, rtol=0)
        energy, gradient = common.energy_and_gradient(q, native)
        common_primitive = ActiveJointCurvatureCoordinates.for_clamped_plane(common_reference,
            boundary=boundary, cell_scale_A=observation["cell_scale_A"])
        np.testing.assert_array_equal(gradient, common_primitive.enthalpy_gradient(native,
            native.get_forces(), native.get_stress(voigt=False), 0.))
        assert energy == native.get_potential_energy()
