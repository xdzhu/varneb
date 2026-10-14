"""Optimizer seeding contracts, not physical stability or DFT acceleration."""
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from vcneb import ModePlane, relax_orthogonal_at_q, sample_conditional_mode_surface
from vcneb.mode_evaluator import CalculatorModeEvaluator
from vcneb.reference_cell import ReferenceCellCoordinates


def system():
    rng = np.random.default_rng(157)
    metric = np.arange(1., 9.)
    rotation = np.linalg.qr(rng.normal(size=(8, 8)))[0]
    basis = rotation/np.sqrt(metric)[:, None]
    plane = ModePlane(np.zeros(8), basis, metric, np.eye(8)[:, :2],
                      ('Q1', 'Q2'), 'test amplitude', 'analytic/seed-contract')
    reduced = np.diag([-2., -1., -3., 1., 2., 5., 10., 30.])
    reduced[0, 5] = reduced[5, 0] = .2
    sqrt_metric = np.sqrt(metric)
    hessian = sqrt_metric[:, None]*(rotation@reduced@rotation.T)*sqrt_metric[None, :]
    linear = rng.normal(size=8)*.03
    frozen = basis[:, 2:3]
    return plane, hessian, linear, frozen


def oracle(plane, hessian, linear, frozen, q):
    constraints = np.column_stack((plane.axis_vectors, frozen)).T*plane.metric_weights[None, :]
    n, k = hessian.shape[0], constraints.shape[0]
    matrix = np.block([[hessian, constraints.T], [constraints, np.zeros((k, k))]])
    return np.linalg.solve(matrix, np.r_[-linear, q, np.zeros(k-2)])[:n]


def test_full_physical_hessian_projects_correctly_in_nonunit_metric_and_keeps_input():
    plane, h, linear, frozen = system(); original = h.copy(); calls = []
    def evaluate(x):
        calls.append(x.copy())
        return float(.5*x@h@x+linear@x), h@x+linear
    q = np.array([.08, -.03])
    result = relax_orthogonal_at_q(plane, q, evaluate, frozen_directions=frozen,
        optimizer='safeguarded_bfgs', trial_validator=lambda x: True,
        initial_hessian=h, initial_trust_radius=5., gradient_tolerance=1e-10)
    np.testing.assert_allclose(result.coordinates, oracle(plane,h,linear,frozen,q), atol=2e-14)
    np.testing.assert_array_equal(h, original)
    assert len(calls) == result.total_n_evaluations == 3
    np.testing.assert_allclose(plane.project(result.coordinates), q, atol=1e-14)


def test_grid_forwards_seed_without_changing_fixed_q_or_all_start_costs():
    plane, h, linear, frozen = system()
    def evaluate(x): return float(.5*x@h@x+linear@x), h@x+linear
    result = sample_conditional_mode_surface(plane, [-.1,.1], [-.1,.1], evaluate,
        frozen_directions=frozen, optimizer='safeguarded_bfgs', trial_validator=lambda x: True,
        initial_hessian=h, initial_trust_radius=5., gradient_tolerance=1e-10,
        energy_unit='test', potential_kind='E', boundary_condition='analytic constraints')
    assert np.all(result.n_evaluations == 3)
    for j, second in enumerate(result.q2):
        for i, first in enumerate(result.q1):
            np.testing.assert_allclose(result.coordinates[j,i],
                oracle(plane,h,linear,frozen,[first,second]), atol=2e-14)


@pytest.mark.parametrize('fault', ['shape','complex','nonfinite','asymmetric','negative_R','zero_R',
                                 'wrong_optimizer','empty_R','unrepresentable_inverse'])
def test_invalid_seed_refuses_before_physical_calls(fault):
    plane,h,linear,frozen = system(); calls = []
    kwargs = dict(optimizer='safeguarded_bfgs', trial_validator=lambda x: True, frozen_directions=frozen)
    if fault == 'shape': h = h[:-1,:-1]
    elif fault == 'complex': h = h.astype(complex)
    elif fault == 'nonfinite': h[0,0] = np.nan
    elif fault == 'asymmetric': h[0,1] += .1
    elif fault == 'negative_R': h = -np.eye(8)
    elif fault == 'zero_R': h = np.zeros((8,8))
    elif fault == 'unrepresentable_inverse': h = np.eye(8)*1e-320
    elif fault == 'wrong_optimizer': kwargs['optimizer'] = 'L-BFGS-B'
    elif fault == 'empty_R': kwargs['frozen_directions'] = plane.basis[:,2:]
    with pytest.raises(ValueError):
        relax_orthogonal_at_q(plane,[0.,0.],lambda x: (calls.append(x.copy()) or (0.,np.zeros(8))),
                              initial_hessian=h,**kwargs)
    assert not calls


def test_none_is_exact_legacy_trace_not_a_new_default():
    plane,h,linear,frozen = system()
    logs = []
    for kwargs in ({}, {'initial_hessian':None}):
        trace = []
        def evaluate(x):
            trace.append(x.copy()); return float(.5*x@h@x+linear@x), h@x+linear
        result = relax_orthogonal_at_q(plane,[.01,.02],evaluate,frozen_directions=frozen,
            optimizer='safeguarded_bfgs',trial_validator=lambda x: True,gradient_tolerance=1e-7,**kwargs)
        logs.append(trace)
    np.testing.assert_array_equal(logs[0],logs[1])


def test_optimistic_seed_still_checks_geometry_before_energy():
    plane = ModePlane(np.zeros(4),np.eye(4),np.ones(4),np.eye(4)[:,:2],
                      ('Q1','Q2'),'test','analytic/guard')
    evaluated, rejected = [], []
    def evaluate(x):
        assert abs(x[2]) <= .45
        evaluated.append(x.copy()); g=np.zeros(4);g[2]=.03*(x[2]-.4)
        return float(.015*(x[2]-.4)**2),g
    def guard(x):
        if abs(x[2])>.45: rejected.append(x.copy());return False
        return True
    result = relax_orthogonal_at_q(plane,[.1,-.2],evaluate,frozen_directions=np.eye(4)[:,3:],
        optimizer='safeguarded_bfgs',trial_validator=guard,initial_trust_radius=1.,
        initial_hessian=np.diag([-1.,-2.,1e-6,-3.]),gradient_tolerance=1e-10)
    assert rejected and result.coordinates[2] == pytest.approx(.4,abs=1e-9)
    assert result.total_n_evaluations == len(evaluated)


def test_seed_does_not_bypass_conservative_calculator_attempt_budget(tmp_path):
    atoms = Atoms('HeLi',scaled_positions=[[.1,.2,.3],[.6,.5,.4]],cell=np.diag([4.,5.,6.]),pbc=True)
    chart = ReferenceCellCoordinates(atoms); n=chart.coordinate_count; factories=[]
    def factory(index,geometry,directory):
        factories.append(index);p=geometry.positions;h=geometry.cell.array
        return SinglePointCalculator(geometry, energy=float(.5*(np.sum(p*p)+.04*np.sum(h*h))),
            forces=-p,stress=(p.T@p+.04*h.T@h)/geometry.get_volume())
    evaluator = CalculatorModeEvaluator(chart,factory,workdir=tmp_path/'capped',
        calculator_id='analytic:seed-budget:v1',max_evaluation_attempts=1)
    plane = ModePlane(np.zeros(n),np.eye(n),np.ones(n),np.eye(n)[:,:2],
                      ('Q1','Q2'),'test','analytic/budget')
    with pytest.raises(RuntimeError,match='budget exhausted'):
        relax_orthogonal_at_q(plane,[0.,0.],evaluator,frozen_directions=np.eye(n)[:,3:],
            optimizer='safeguarded_bfgs',trial_validator=lambda x: True,
            initial_trust_radius=.01,initial_hessian=np.eye(n))
    assert factories == [0] and evaluator.n_evaluation_attempts == evaluator.n_new_evaluations == 1


def test_measured_H_replay_is_explicit_quadratic_validation_not_material_speed_claim():
    from scripts.check_hfo2_conditional_seed import check
    result = check()
    assert result['new_DFT_or_scheduler_operations'] == 0
    assert result['DFT_speedup_or_novel_algorithm_claimed'] is False
    assert result['physical_conditional_axes_domains_or_forecasts_selected'] is False
    assert len(result['rows']) == 4
    for channel in ('T','preserving'):
        baseline,seeded = [r for r in result['rows'] if r['channel']==channel]
        assert baseline['status'] == seeded['status'] == 'converged'
        assert seeded['callback_evaluations'] < baseline['callback_evaluations']
        for row in (baseline,seeded):
            assert row['coordinate_oracle_error_A'] < 1e-7
            assert row['energy_oracle_error_eV'] < 1e-12
