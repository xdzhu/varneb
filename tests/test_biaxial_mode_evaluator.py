"""Controlled strain must survive conditional relaxation, caching and raw EFS."""
import json
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms

from vcneb import (
    BiaxialClampedModeCoordinates, CalculatorModeEvaluator, ModePlane,
    audit_orthogonal_curvature, clamped_plane_vcneb_boundary,
    relax_orthogonal_at_q,
)
from scripts.analyze_hfo2_clamped_reference import terminal_images


def chart_fixture(*, tilted=True, rotated=False, anchor=.01, scale=4.2, cell_scale=3.6):
    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = atoms.cell.array + [[0., .10, .15], [.05, 0., .08], [.11, -.07, 0.]]
    if rotated:
        rotation, _ = np.linalg.qr(np.random.default_rng(181).normal(size=(3, 3)))
        cell = cell @ rotation.T
    atoms.set_cell(cell, scale_atoms=True)
    boundary = clamped_plane_vcneb_boundary(len(atoms), cell, allow_tilt=tilted)
    return BiaxialClampedModeCoordinates(atoms, boundary=boundary, cell_scale_A=cell_scale,
        anchor_strain=anchor, strain_scale_A=scale)


def make_plane(chart, axes):
    return ModePlane(reference=np.zeros(chart.size), basis=np.eye(chart.size),
        metric_weights=np.ones(chart.size), axis_weights=axes,
        axis_labels=("Cu implementation check Q", "scaled prescribed strain"),
        amplitude_unit="Angstrom", reference_id="Cu-test-not-HfO2-G3")


@pytest.mark.parametrize("rotated", [False, True])
@pytest.mark.parametrize("tilted", [False, True])
@pytest.mark.parametrize("pressure", [0., .03])
def test_evaluator_full_gradient_matches_direct_energy_differences(tmp_path, rotated, tilted, pressure):
    chart = chart_fixture(rotated=rotated, tilted=tilted)
    calls = []
    def factory(index, atoms, directory):
        calls.append(index)
        return EMT()
    evaluate = CalculatorModeEvaluator(chart, factory, calculator_id="Cu-EMT",
        workdir=tmp_path/"gradient", pressure_eV_per_A3=pressure)
    y = np.random.default_rng(183).normal(size=chart.size)*.007
    energy, gradient = evaluate(y)
    numerical = []
    for i in range(chart.size):
        pair = []
        for sign in (1., -1.):
            delta = y.copy(); delta[i] += sign*1e-5
            atoms = chart.to_atoms(delta); atoms.calc = EMT()
            pair.append(atoms.get_potential_energy()+pressure*atoms.get_volume())
        numerical.append((pair[0]-pair[1])/2e-5)
    np.testing.assert_allclose(gradient, numerical, atol=3e-5, rtol=0)
    assert evaluate(y)[0] == energy and calls == [0]
    np.testing.assert_array_equal(evaluate(y)[1], gradient)
    assert evaluate.contract["mechanical_boundary"]["anchor_strain"] == .01
    assert evaluate.contract["mechanical_boundary"]["prescribed_substrate_relaxable"] is False


@pytest.mark.parametrize("rotated", [False, True])
@pytest.mark.parametrize("tilted", [False, True])
def test_inverse_preserves_open_cell_and_long_periodic_lifts(rotated, tilted):
    chart = chart_fixture(rotated=rotated, tilted=tilted)
    y = np.random.default_rng(184).normal(size=chart.size)*.01
    y[:3] += 1.3*chart.reference_cell[0]
    atoms = chart.to_atoms(y)
    np.testing.assert_allclose(chart.from_atoms(atoms), y, atol=1e-12, rtol=0)
    wrapped = atoms.copy(); wrapped.wrap()
    assert not np.allclose(chart.from_atoms(wrapped), y)
    reconstructed = chart.to_atoms(chart.from_atoms(wrapped))
    np.testing.assert_allclose(reconstructed.positions, wrapped.positions, atol=1e-12)
    assert reconstructed.calc is None


@pytest.mark.parametrize("changed", [dict(anchor=.02), dict(scale=5.), dict(cell_scale=4.), dict(tilted=False)])
def test_resume_keeps_physical_chart_and_control_scales(tmp_path, changed):
    chart = chart_fixture()
    calls = []
    def factory(index, atoms, directory):
        calls.append(index); return EMT()
    work = tmp_path/"resume"
    evaluate = CalculatorModeEvaluator(chart, factory, calculator_id="Cu-EMT", workdir=work,
        max_evaluation_attempts=1)
    y = np.zeros(chart.size); y[-1] = .001*chart.strain_scale_A
    expected = evaluate(y)
    resumed = CalculatorModeEvaluator(chart, factory, calculator_id="Cu-EMT", workdir=work,
        max_evaluation_attempts=1, resume=True)
    assert resumed(y)[0] == expected[0] and calls == [0]
    with pytest.raises(ValueError, match="resume contract changed"):
        CalculatorModeEvaluator(chart_fixture(**changed), factory, calculator_id="Cu-EMT",
            workdir=work, max_evaluation_attempts=1, resume=True)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        resumed(y + .00001)
    assert resumed.n_evaluation_attempts == 1 and calls == [0]


@pytest.mark.parametrize("fault", ["cell", "positions", "wrap", "species", "pbc", "constraint", "stress"])
def test_bad_backend_data_never_enter_cache(tmp_path, fault):
    chart = chart_fixture()
    def factory(index, atoms, directory):
        if fault == "cell": atoms.cell[0, 0] += .01
        if fault == "positions": atoms.positions[0, 0] += .01
        if fault == "wrap": atoms.positions[0] += atoms.cell[0]
        if fault == "species": atoms.numbers[0] = 8
        if fault == "pbc": atoms.pbc[0] = False
        if fault == "constraint": atoms.set_constraint(FixAtoms(indices=[0]))
        stress = np.eye(3)
        if fault == "stress": stress[0, 1] = .1
        return SinglePointCalculator(atoms, energy=1., forces=np.zeros((len(atoms), 3)), stress=stress)
    evaluate = CalculatorModeEvaluator(chart, factory, calculator_id="invalid", workdir=tmp_path/"bad")
    with pytest.raises(ValueError): evaluate(np.zeros(chart.size))
    assert evaluate.n_new_evaluations == 0
    assert len(list((tmp_path/"bad").glob("eval-*/failure.json"))) == 1
    assert not list((tmp_path/"bad").glob("eval-*/result.json"))


@pytest.mark.parametrize("fault", ["unequal_scale", "rotation", "constraint", "collapse", "nan"])
def test_inverse_rejects_incompatible_geometry_without_repair(fault):
    chart = chart_fixture()
    atoms = chart.reference_atoms
    if fault == "unequal_scale": atoms.cell[0] *= 1.001
    if fault == "rotation": atoms.cell[0, 1] += .001
    if fault == "constraint": atoms.set_constraint(FixAtoms(indices=[0]))
    if fault == "collapse": atoms.cell[2] *= -1.
    if fault == "nan": atoms.positions[0, 0] = np.nan
    before = atoms.cell.array.copy()
    with pytest.raises(ValueError): chart.from_atoms(atoms)
    np.testing.assert_array_equal(atoms.cell.array, before)


@pytest.mark.parametrize("operation", ["relax", "curvature"])
@pytest.mark.parametrize("mixed_control_axis", [False, True])
def test_releasing_prescribed_strain_is_rejected_before_factory(tmp_path, operation, mixed_control_axis):
    chart = chart_fixture()
    calls = []
    evaluate = CalculatorModeEvaluator(chart, lambda *args: calls.append(1),
        calculator_id="must-not-run", workdir=tmp_path/"guard")
    axes = np.zeros((chart.size, 2)); axes[0, 0] = 1.; axes[1, 1] = 1.
    if mixed_control_axis:
        axes[0, 0] = 1/np.sqrt(2); axes[-1, 0] = 1/np.sqrt(2)
    plane = make_plane(chart, axes)
    with pytest.raises(ValueError, match="release a prescribed control"):
        if operation == "relax": relax_orthogonal_at_q(plane, [0., 0.], evaluate)
        else: audit_orthogonal_curvature(plane, np.zeros(chart.size), evaluate, step_amplitude=.01)
    assert calls == [] and evaluate.n_evaluation_attempts == 0


@pytest.mark.parametrize("optimizer", ["L-BFGS-B", "safeguarded_bfgs"])
def test_fixed_control_is_preserved_in_every_actual_conditional_evaluation(tmp_path, optimizer):
    chart = chart_fixture(rotated=True)
    axes = np.zeros((chart.size, 2))
    axes[[0, 3], 0] = np.array([1., -1.])/np.sqrt(2)
    axes[-1, 1] = 1.
    plane = make_plane(chart, axes)
    evaluate = CalculatorModeEvaluator(chart, lambda *args: EMT(),
        calculator_id="Cu-EMT", workdir=tmp_path/"conditional", max_evaluation_attempts=120)
    q = np.array([.004, .001*chart.strain_scale_A])
    result = relax_orthogonal_at_q(plane, q, evaluate,
        frozen_directions=chart.rigid_translation_directions(), optimizer=optimizer,
        trial_validator=lambda y: np.linalg.det(chart.to_atoms(y).cell.array) > 0,
        gradient_tolerance=1e-4, max_iterations=60, initial_trust_radius=.03,
        orthogonal_amplitude_bound=.3)
    assert result.orthogonal_gradient_norm <= 1e-4
    np.testing.assert_allclose(plane.project(result.coordinates), q, atol=1e-12)
    assert evaluate.n_new_evaluations > 1
    for path in (tmp_path/"conditional").glob("eval-*/result.json"):
        y = np.array(json.loads(path.read_text())["coordinates"])
        assert y[-1] == pytest.approx(q[1], abs=1e-12)
        np.testing.assert_allclose(chart.to_atoms(y).cell.array[:2],
            chart.reference_cell[:2]*(1+chart.anchor_strain+.001)/(1+chart.anchor_strain), atol=1e-12)


def test_explicitly_frozen_control_is_allowed_and_metadata_are_defensive(tmp_path):
    chart = chart_fixture()
    axes = np.zeros((chart.size, 2)); axes[0, 0] = 1.; axes[1, 1] = 1.
    plane = make_plane(chart, axes)
    evaluate = CalculatorModeEvaluator(chart, lambda *args: EMT(),
        calculator_id="Cu-EMT", workdir=tmp_path/"frozen")
    # Freeze the entire complement: one raw point, not a stability certificate.
    frozen = np.eye(chart.size)[:, 2:]
    result = relax_orthogonal_at_q(plane, [.001, -.001], evaluate, frozen_directions=frozen)
    assert result.coordinates[-1] == 0 and evaluate.n_new_evaluations == 1
    controls = evaluate.fixed_coordinate_covectors; controls[-1] = 0.
    assert evaluate.fixed_coordinate_covectors[-1, 0] == 1.
    contract = chart.mechanical_contract; contract["anchor_strain"] = 90.
    assert chart.mechanical_contract["anchor_strain"] == .01


@pytest.mark.parametrize("case,digest,anchor", [
    ("clamped_T_terminal_E079_20261011", "c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e", 0.),
    ("clamped_plus1_reversing_terminal_E148_20261014/completed_HF", "450dca48d9c8fec24dd649df868c47b122b340a88eaf89df941811855ef590e7", .01),
])
def test_native_hafnia_frames_replay_at_their_actual_strain_without_new_DFT(case, digest, anchor):
    root = Path(__file__).resolve().parents[1]/"benchmarks/hfo2_channels/20261008"
    images, observation, _, _ = terminal_images(root/case, digest)
    boundary = clamped_plane_vcneb_boundary(len(images[0]),
        np.array(observation["mechanical_boundary"]["reference_cell_A"]), allow_tilt=True)
    chart = BiaxialClampedModeCoordinates(images[0], boundary=boundary,
        cell_scale_A=observation["cell_scale_A"], anchor_strain=anchor,
        strain_scale_A=observation["cell_scale_A"])
    assert chart.size == 40
    for native in images:
        y = chart.from_atoms(native)
        assert y[-1] == pytest.approx(0., abs=1e-12)
        restored = chart.to_atoms(y)
        np.testing.assert_allclose(restored.positions, native.positions, atol=1e-12)
        np.testing.assert_allclose(restored.cell.array, native.cell.array, atol=1e-12)
        energy, gradient = chart.energy_and_gradient(y, native)
        assert energy == native.get_potential_energy()
        assert gradient.shape == (40,) and np.isfinite(gradient).all()
        np.testing.assert_allclose(gradient, chart.enthalpy_gradient(y, native,
            native.get_forces(), native.get_stress(voigt=False), 0.), atol=1e-13)


def test_one_native_common_chart_represents_both_observed_training_conditions():
    root = Path(__file__).resolve().parents[1]/"benchmarks/hfo2_channels/20261008"
    zero, observation, _, _ = terminal_images(root/"clamped_T_terminal_E079_20261011",
        "c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e")
    plus, _, _, _ = terminal_images(root/"clamped_plus1_reversing_terminal_E148_20261014/completed_HF",
        "450dca48d9c8fec24dd649df868c47b122b340a88eaf89df941811855ef590e7")
    boundary = clamped_plane_vcneb_boundary(len(zero[0]),
        np.array(observation["mechanical_boundary"]["reference_cell_A"]), allow_tilt=True)
    chart = BiaxialClampedModeCoordinates(zero[0], boundary=boundary,
        cell_scale_A=observation["cell_scale_A"], anchor_strain=0.,
        strain_scale_A=observation["cell_scale_A"])
    for epsilon, images in ((0., zero), (.01, plus)):
        for native in images:
            y = chart.from_atoms(native)
            assert y[-1]/chart.strain_scale_A == pytest.approx(epsilon, abs=1e-12)
            reconstructed = chart.to_atoms(y)
            np.testing.assert_allclose(reconstructed.cell.array, native.cell.array, atol=1e-12)
            np.testing.assert_allclose(reconstructed.positions, native.positions, atol=1e-12)
            energy, gradient = chart.energy_and_gradient(y, native)
            assert energy == native.get_potential_energy() and np.isfinite(gradient).all()


@pytest.mark.parametrize("controls", [np.array([1., 0., 0.]), np.zeros((3, 1)),
                                     np.full((3, 1), np.nan), np.ones((4, 1))])
def test_invalid_fixed_control_metadata_is_rejected_before_any_callback(controls):
    plane = ModePlane(reference=np.zeros(3), basis=np.eye(3), metric_weights=np.ones(3),
        axis_weights=np.eye(3)[:, :2], axis_labels=("Q1", "Q2"),
        amplitude_unit="test", reference_id="not-material")
    class Evaluator:
        fixed_coordinate_covectors = controls
        calls = 0
        def __call__(self, y):
            self.calls += 1
            return 0., np.zeros(3)
    evaluator = Evaluator()
    with pytest.raises(ValueError, match="fixed_coordinate_covectors"):
        relax_orthogonal_at_q(plane, [0., 0.], evaluator)
    assert evaluator.calls == 0
