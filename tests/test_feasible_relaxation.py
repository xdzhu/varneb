"""Analytic feasibility tests are not material acceleration evidence."""
import numpy as np
import pytest
from vcneb import ModePlane
from vcneb.feasible_relaxation import relax_orthogonal_in_ball


def plane():
    return ModePlane(
        np.zeros(4),
        np.eye(4),
        np.ones(4),
        np.eye(4)[:, :2],
        ("Q1", "Q2"),
        "test amplitude",
        "analytic/ball",
    )


def test_interior_seeded_quadratic_keeps_Q_radius_guard_and_all_callback_costs():
    p = plane()
    calls = []
    target = np.array([0.1, -0.2, 0.2, -0.1])

    def evaluate(x):
        assert np.linalg.norm(x[2:]) <= 0.5
        np.testing.assert_allclose(x[:2], target[:2], atol=1e-14)
        calls.append(x.copy())
        return float(0.5 * np.sum((x - target) ** 2)), x - target

    result = relax_orthogonal_in_ball(
        p,
        target[:2],
        evaluate,
        start=np.r_[target[:2], 0.0, 0.0],
        orthogonal_radius_bound=0.5,
        trial_validator=lambda x: True,
        initial_hessian=np.eye(4),
        initial_trust_radius=1.0,
        gradient_tolerance=1e-12,
    )
    np.testing.assert_allclose(result.coordinates, target, atol=1e-14)
    assert result.n_evaluations == result.total_n_evaluations == len(calls) == 3


def test_boundary_KKT_is_not_full_gradient_convergence():
    p = plane()
    seen = []

    def evaluate(x):
        seen.append(x.copy())
        g = x - np.array([0.0, 0.0, 1.0, 0.0])
        return float(0.5 * g @ g), g

    with pytest.raises(RuntimeError, match="full released-gradient"):
        relax_orthogonal_in_ball(
            p,
            [0, 0],
            evaluate,
            start=np.zeros(4),
            orthogonal_radius_bound=0.5,
            trial_validator=lambda x: True,
            initial_trust_radius=1.0,
            gradient_tolerance=0.03,
            max_iterations=30,
        )
    assert seen and all(np.linalg.norm(x[2:]) <= 0.5 for x in seen)


@pytest.mark.parametrize(
    "fault",
    [
        "wrong_Q",
        "outside_radius",
        "guard",
        "negative_seed",
        "complex_seed",
        "nonfinite_seed",
        "bad_radius",
        "bad_target",
        "bad_iterations",
    ],
)
def test_invalid_contract_refuses_before_any_evaluation(fault):
    p = plane()
    x = np.zeros(4)
    calls = []
    kw = dict(start=x, orthogonal_radius_bound=0.5, trial_validator=lambda x: True)
    if fault == "wrong_Q":
        x[0] = 0.1
    elif fault == "outside_radius":
        x[2] = 0.6
    elif fault == "guard":
        kw["trial_validator"] = lambda x: False
    elif fault == "negative_seed":
        kw["initial_hessian"] = -np.eye(4)
    elif fault == "complex_seed":
        kw["initial_hessian"] = np.eye(4).astype(complex)
    elif fault == "nonfinite_seed":
        kw["initial_hessian"] = np.eye(4) * np.nan
    elif fault == "bad_radius":
        kw["orthogonal_radius_bound"] = np.nan
    elif fault == "bad_target":
        kw["gradient_tolerance"] = -1.0
    else:
        kw["max_iterations"] = True
    with pytest.raises(ValueError):
        relax_orthogonal_in_ball(
            p, [0, 0], lambda x: (calls.append(x.copy()) or (0.0, x)), **kw
        )
    assert not calls


def test_nonunit_metric_and_frozen_coordinate_are_not_projected_away():
    w = np.array([2.0, 3.0, 4.0, 5.0, 6.0])
    basis = np.eye(5) / np.sqrt(w)[:, None]
    p = ModePlane(
        np.zeros(5),
        basis,
        w,
        np.eye(5)[:, :2],
        ("Q1", "Q2"),
        "test",
        "analytic/metric-ball",
    )
    target = basis @ np.array([0.1, -0.2, 0.0, 0.15, -0.2])
    h = np.diag(w)

    def evaluate(x):
        return float(0.5 * (x - target) @ h @ (x - target)), h @ (x - target)

    result = relax_orthogonal_in_ball(
        p,
        [0.1, -0.2],
        evaluate,
        start=p.frozen_coordinates([0.1, -0.2]),
        frozen_directions=basis[:, 2:3],
        orthogonal_radius_bound=0.5,
        initial_hessian=h,
        trial_validator=lambda x: True,
        initial_trust_radius=1.0,
        gradient_tolerance=1e-11,
    )
    np.testing.assert_allclose(result.coordinates, target, atol=1e-13)


def test_fixed_control_cannot_be_released_by_the_new_solver():
    class Evaluator:
        fixed_coordinate_covectors = np.eye(4)[:, 3:]

        def __call__(self, x):
            pytest.fail("must reject release of fixed control first")

    with pytest.raises(ValueError, match="prescribed control"):
        relax_orthogonal_in_ball(
            plane(),
            [0, 0],
            Evaluator(),
            start=np.zeros(4),
            orthogonal_radius_bound=0.5,
            trial_validator=lambda x: True,
        )


def test_callback_budget_exception_propagates_and_is_never_retried():
    calls = []

    def evaluate(x):
        calls.append(x.copy())
        if len(calls) == 2:
            raise RuntimeError("explicit analytic budget exhausted")
        g = x - np.array([0.0, 0.0, 0.2, 0.1])
        return float(0.5 * g @ g), g

    with pytest.raises(RuntimeError, match="budget exhausted"):
        relax_orthogonal_in_ball(
            plane(),
            [0, 0],
            evaluate,
            start=np.zeros(4),
            orthogonal_radius_bound=0.5,
            trial_validator=lambda x: True,
            initial_trust_radius=0.01,
        )
    assert len(calls) == 2


def test_ASE_EFS_adapter_retains_its_conservative_factory_budget(tmp_path):
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointCalculator
    from vcneb.mode_evaluator import CalculatorModeEvaluator
    from vcneb.reference_cell import ReferenceCellCoordinates

    atoms = Atoms(
        "HeLi",
        scaled_positions=[[0.1, 0.2, 0.3], [0.6, 0.5, 0.4]],
        cell=np.diag([4.0, 5.0, 6.0]),
        pbc=True,
    )
    chart = ReferenceCellCoordinates(atoms)
    n = chart.coordinate_count
    factories = []

    def factory(index, geometry, directory):
        factories.append(index)
        p = geometry.positions
        h = geometry.cell.array
        return SinglePointCalculator(
            geometry,
            energy=float(0.5 * (np.sum(p * p) + 0.04 * np.sum(h * h))),
            forces=-p,
            stress=(p.T @ p + 0.04 * h.T @ h) / geometry.get_volume(),
        )

    evaluator = CalculatorModeEvaluator(
        chart,
        factory,
        workdir=tmp_path / "capped",
        calculator_id="analytic:feasible-budget:v1",
        max_evaluation_attempts=1,
    )
    p = ModePlane(
        np.zeros(n),
        np.eye(n),
        np.ones(n),
        np.eye(n)[:, :2],
        ("Q1", "Q2"),
        "test amplitude",
        "analytic/feasible-ASE-budget",
    )
    with pytest.raises(RuntimeError, match="budget exhausted"):
        relax_orthogonal_in_ball(
            p,
            [0, 0],
            evaluator,
            start=np.zeros(n),
            orthogonal_radius_bound=0.1,
            trial_validator=lambda x: True,
            initial_trust_radius=0.01,
            initial_hessian=np.eye(n),
        )
    assert (
        factories == [0]
        and evaluator.n_evaluation_attempts == evaluator.n_new_evaluations == 1
    )
