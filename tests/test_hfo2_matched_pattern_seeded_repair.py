"""A finite numerical repair, not a new physical-input or prediction search."""
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


repair = load(BASE/'matched_pattern_seeded_E271_20261018/repair.py', 'E271_test')
study = load(BASE/'matched_pattern_release_E270_20261018/study.py', 'E270_repair_test')
numeric = load(BASE/'matched_pattern_release_E270_20261018/numerical_seed.py', 'E270_numeric_test')


@pytest.fixture(scope='module')
def source():
    case = BASE/'matched_pattern_release_E270_20261018'
    seed = case/'prepared/seeds/0'
    chart = study.chart_for(seed)
    return chart, json.loads((seed/'seed.json').read_bytes()), numeric.seed_for(case, 0)[0]


def test_preparation_freezes_first_probe_and_does_not_run_DFT(tmp_path, monkeypatch):
    from examples.hfo2_fixed_input_factory import FixedHfo2Calculator
    monkeypatch.setattr(FixedHfo2Calculator, 'calculate', lambda *a, **k: pytest.fail('no DFT in preparation'))
    output = tmp_path/'prepared.json'
    receipt = repair.prepare(BASE/'matched_pattern_seeded_E271_20261018', output)
    data = json.loads(output.read_bytes())
    assert receipt['new_DFT_calls'] == data['new_DFT_calls'] == 0
    assert data['base_bundle_pin'] == repair.BASE_PIN and data['base_job'] == '28947338'
    assert len(data['matrices']) == 3
    assert all(r['measured_probe_step_A'] == .01 for r in data['metadata'])
    assert not data['policy']['eigenvalue_clipping']
    with pytest.raises(FileExistsError):
        repair.prepare(BASE/'matched_pattern_seeded_E271_20261018', output)


def originals():
    unavailable = dict(treatment='A', status='attempt_budget_exhausted', candidate=None, new_attempts=12)
    skip = dict(treatment='J', status='not_attempted_atomic_candidate_unavailable', candidate=None, new_attempts=0)
    available = dict(treatment='A', status='gradient_candidate_stability_unchecked', candidate={}, new_attempts=5)
    failed_j = dict(treatment='J', status='attempt_budget_exhausted', candidate=None, new_attempts=12)
    return [[dict(unavailable), dict(skip)], [dict(available), dict(failed_j)],
            [dict(available), dict(available, treatment='J')]]


def test_plan_repairs_only_unavailable_arms_and_never_erases_baseline():
    rows = originals(); before = json.dumps(rows)
    assert repair.baseline_plan(rows) == [dict(index=0, repair_A=True, repair_J=True),
        dict(index=1, repair_A=False, repair_J=True), dict(index=2, repair_A=False, repair_J=False)]
    assert json.dumps(rows) == before


@pytest.mark.parametrize('fault', ['count', 'J_before_A', 'attempts', 'status', 'false_skip'])
def test_invalid_or_incomplete_baseline_refused(fault):
    rows = originals()
    if fault == 'count': rows.pop()
    elif fault == 'J_before_A': rows[0][1] = dict(rows[2][1])
    elif fault == 'attempts': rows[0][0]['new_attempts'] = 13
    elif fault == 'status': rows[1][1]['status'] = 'calculator_error'
    else: rows[1][1] = dict(rows[0][1])
    with pytest.raises(ValueError): repair.baseline_plan(rows)


def test_seeded_quadratic_keeps_physical_patterns_and_exact_A_cell(source):
    chart, data, h = source
    _, g = study.cached_frozen(BASE/'matched_pattern_release_E270_20261018/prepared/seeds/0', chart)
    def quadratic(x):
        return float(g@x+.5*x@h@x), g+h@x
    cache = study.ExactCache(quadratic, 0., g)
    row = repair.seeded_arm(study, chart, data, cache, h, 'A')
    assert row['candidate'] is not None and 0 < row['new_attempts'] <= 12
    x = np.asarray(row['candidate']['coordinates'])
    assert not np.any(x[36:])
    np.testing.assert_allclose(np.asarray(data['covectors_39_by_4']).T@x, 0., atol=1e-10, rtol=0)
    assert np.array_equal(chart.to_atoms(x).cell.array, chart.reference_cell)
    assert row['candidate']['released_gradient_norm_eV_A'] <= .03
    assert not row['candidate']['stationary_TS_or_stability_certified']
    assert row['origin'] == 'E271_seeded_repair'


def test_no_A_means_no_J_calculator_attempt(source):
    chart, data, h = source
    cache = study.ExactCache(lambda x: pytest.fail('unavailable A cannot start J'), 0., np.zeros(39))
    row = repair.seeded_arm(study, chart, data, cache, h, 'J')
    assert row['candidate'] is None and row['new_attempts'] == cache.fresh == 0


def test_seed_never_bypasses_cap_or_changes_failure_to_candidate(source):
    chart, data, h = source
    plane, frozen, _, _ = repair.release_space(data, True, np.zeros(39))
    direction = np.r_[repair._orthogonal_directions(plane, frozen)[:, 0], np.zeros(3)]
    g = .2*direction
    cache = study.ExactCache(lambda x: (float(g@x), g), 0., g)
    row = repair.seeded_arm(study, chart, data, cache, h, 'A')
    assert row['candidate'] is None and row['new_attempts'] <= 12
    assert cache.fresh == row['new_attempts']


def test_actual_seeded_zero_A_repairs_one_budget_failure_not_the_original_label():
    case = BASE/'matched_pattern_seeded_E271_20261018'
    reader = load(BASE/'matched_pattern_release_E270_20261018/review_arm.py', 'E271_native_reader_test')
    r = reader.review(BASE/'matched_pattern_release_E270_20261018', case/'zero_A_prefix', 0, 'A')
    assert r['total_arm_attempts_reported'] == 2
    assert r['released_gradient_norm_eV_A'] == pytest.approx(.02101509741313)
    assert r['energy_change_meV_fu'] == pytest.approx(-.057641912008)
    assert r['physical_retained_defect_A'] < 1.5e-18 and r['cell_coordinate_norm_A'] == 0
    assert r['complete_allowed_gradient_norm_eV_A'] == pytest.approx(.10574104875965)
    assert not r['energy_partition_or_stationary_TS_or_minimum_certified']
