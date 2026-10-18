"""Bounded signed-branch execution plumbing, not a DFT speed benchmark."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

BASE = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008'


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


worker = load(BASE/'preserving_pattern_branches_E272_20261018/worker.py', 'P_worker_test')
study = load(BASE/'matched_pattern_release_E270_20261018/study.py', 'P_worker_native')
helper = load(BASE/'matched_pattern_seeded_E271_20261018/repair.py', 'P_worker_space')


@pytest.fixture(scope='module')
def prepared():
    seed = BASE/'preserving_pattern_branches_E272_20261018/prepared/seeds/0'
    data = json.loads((seed/'seed.json').read_bytes())
    p = json.loads((seed/'branch_preparation.json').read_bytes())
    chart = study.chart_for(seed); e, g = study.cached_frozen(seed, chart)
    return chart, data, p, e, g


@pytest.mark.parametrize('sign', [-1., 1.])
def test_signed_A_known_direction_and_cached_J_keep_physical_controls(prepared, sign):
    chart, data, p, e, g = prepared
    h = np.asarray(p['numerical_hessian_39']); v = np.asarray(p['signed_starts'][1])/.1
    hard = h-np.outer(v,v); g = g-v*(v@g); g[36:] = 0
    k = .4474490104
    def model(x):
        a = v@x
        return float(e+g@x+.5*x@hard@x-.5*k*a*a+.25*k*a**4/.3**2), g+hard@x+(-k*a+k*a**3/.3**2)*v
    cache = study.ExactCache(model, e, g)
    a = worker.arm(study, helper, chart, data, cache, h, sign*.1*v, 'A')
    assert a['candidate'] is not None and a['new_attempts'] <= 12
    x = np.asarray(a['candidate']['coordinates'])
    assert .26 < sign*(v@x) < .34 and not np.any(x[36:])
    np.testing.assert_allclose(np.asarray(data['covectors_39_by_4']).T@x, 0, atol=1e-10, rtol=0)
    check = worker.directional_check(study, chart, cache, a['candidate'], v)
    assert check['new_attempts'] == 2 and check['status'] == 'known_direction_positive_at_one_step'
    assert not check['whole_stability_certified']
    j = worker.arm(study, helper, chart, data, cache, h, x, 'J')
    assert j['candidate'] is not None and j['new_attempts'] == 0
    assert cache.fresh == a['new_attempts']+2
    assert not j['candidate']['stationary_TS_or_stability_certified']


def test_missing_A_blocks_checks_and_J_without_evaluation(prepared):
    chart, data, p, e, g = prepared
    cache = study.ExactCache(lambda x: pytest.fail('no calls without A'), e, g)
    row = worker.arm(study, helper, chart, data, cache, np.asarray(p['numerical_hessian_39']), None, 'J')
    check = worker.directional_check(study, chart, cache, None, np.zeros(39))
    assert row['candidate'] is None and row['new_attempts'] == check['new_attempts'] == cache.fresh == 0


def test_direction_check_respects_existing_geometric_bound(prepared):
    chart, _, p, e, g = prepared
    v = np.asarray(p['signed_starts'][1])/.1
    cache = study.ExactCache(lambda x: pytest.fail('guarded check cannot call DFT'), e, g)
    row = worker.directional_check(study, chart, cache, dict(coordinates=(.499*v).tolist()), v)
    assert row['status'] == 'original_geometry_guard_prevented_checks'
    assert row['new_attempts'] == cache.fresh == 0


def test_nonstationary_linear_model_remains_a_bounded_failure(prepared):
    chart, data, p, e, _ = prepared
    v = np.asarray(p['signed_starts'][1])/.1; g = -.2*v
    cache = study.ExactCache(lambda x: (float(e+g@x), g), e, g)
    row = worker.arm(study, helper, chart, data, cache, np.asarray(p['numerical_hessian_39']), .1*v, 'A')
    assert row['candidate'] is None and row['new_attempts'] <= 12
    assert row['new_attempts'] == cache.fresh
