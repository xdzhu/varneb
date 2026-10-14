"""Signed-treatment arithmetic; injected atomic candidates are synthetic, not DFT."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'benchmarks/hfo2_channels/20261008'
CASE = BASE / 'G3_treatment_comparison_E178_20261015'
spec = importlib.util.spec_from_file_location('treatment_comparison_tests', CASE / 'analyze.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
EMPTY = BASE / 'G3_atomic_audit_E176_20261015/native_prefix.json'
EMPTY_SHA = 'b7027bd8db935fbc9e87f8b56a92999dbedeab50454e75538c2e5442c3eccd45'
CONFIG = json.loads((BASE / 'G3_atomic_release_E175_20261015/config.json').read_bytes())


def values(f, x, j):
    return dict(frozen=f, atomic=x, joint=j)


def test_missing_atomic_keeps_known_joint_frozen_absolute_and_shape():
    out = a.changes(values(.2, None, .1), values(0., None, -.05))
    assert out['joint_minus_frozen']['absolute_eV_cell'] == pytest.approx(-.1)
    assert out['joint_minus_frozen']['shape_eV_cell'] == pytest.approx(-.05)
    assert out['atomic_minus_frozen'] == dict(absolute_eV_cell=None, shape_eV_cell=None)
    assert out['joint_minus_atomic'] == dict(absolute_eV_cell=None, shape_eV_cell=None)


def test_missing_atomic_origin_retains_absolute_and_does_not_hide_joint_shape():
    out = a.changes(values(0., .2, .3), values(0., None, .1))
    assert out['atomic_minus_frozen']['absolute_eV_cell'] == .2
    assert out['atomic_minus_frozen']['shape_eV_cell'] is None
    assert out['joint_minus_atomic']['absolute_eV_cell'] == pytest.approx(.1)
    assert out['joint_minus_frozen']['shape_eV_cell'] == pytest.approx(.2)


@pytest.mark.parametrize('seed', range(6))
def test_all_three_signed_error_identities_and_common_zero_invariance(seed):
    v = np.random.default_rng(seed).normal(size=(4, 3))
    inputs = [values(*row) for row in v]
    result = a.compare(*inputs)
    shifted = a.compare(*[{k: x + 20 for k, x in row.items()} for row in inputs])
    for group in result:
        for name in result[group]:
            for field, x in result[group][name].items():
                assert shifted[group][name][field] == pytest.approx(x, abs=2e-12)
    for field in ('absolute_eV_cell', 'shape_eV_cell'):
        miss = result['missed_treatment_change']
        assert miss['atomic_minus_frozen'][field] + miss['joint_minus_atomic'][field] == pytest.approx(
            miss['joint_minus_frozen'][field], abs=2e-12)


@pytest.mark.parametrize('bad', [float('nan'), float('inf'), True, '0.'])
def test_unknown_is_none_not_nan_bool_or_string(bad):
    with pytest.raises(ValueError, match='finite numeric'):
        a.changes(values(0., bad, 1.), values(0., 0., 0.))


@pytest.fixture(scope='module')
def actual():
    return a.analyze(EMPTY, EMPTY_SHA)


def test_actual_twenty_rows_preserve_known_observations_and_abstain_atomic(actual):
    assert len(actual['rows']) == 4
    assert [(r['channel'], r['probe_step_A']) for r in actual['rows']] == [
        ('T', .01), ('T', .02), ('preserving', .01), ('preserving', .02)]
    for row in actual['rows']:
        assert [p['id'] for p in row['points']] == list(range(5))
        for p in row['points']:
            assert p['observed']['atomic'] is None and p['atomic_status'] == 'not_processed'
            assert p['observed_changes']['joint_minus_frozen']['absolute_eV_cell'] is not None
            assert p['model_minus_observed']['atomic']['absolute_eV_cell'] is None
    assert actual['costs']['atomic_completed_nodes'] == actual['costs']['atomic_new_SCFs'] == 0
    assert actual['atomic_audit_mode'] == 'prefix' and '|PENDING|' in actual['atomic_source_sacct']
    assert 'excludes in-flight/unaudited' in actual['costs']['atomic_cost_scope']
    assert actual['costs']['frozen_completed_call_proxy_core_hours'] == pytest.approx(8.13681406982243, abs=2e-12, rel=0)
    assert actual['costs']['original_curvature_SCFs_additional'] == 296
    assert actual['no_step_or_point_selection'] and actual['formula_units'] == 4
    assert not actual['new_raw_SCF_reparse']
    assert actual['independent_slots_used'] == actual['new_DFT_or_scheduler_writes'] == 0
    assert not actual['unique_mode_energy_TS_conditional_stability_or_predictive_advantage_certified']


def test_statistics_never_call_missing_atomic_zero_error_and_exclude_origin(actual):
    for row in actual['rows']:
        stats = row['non_origin_coverage_statistics']
        for field in ('absolute_eV_cell', 'shape_eV_cell'):
            joint = stats['model_minus_observed']['joint'][field]
            assert joint['ids'] == [1, 2, 3, 4] and joint['complete_four_axial_points']
            atom = stats['model_minus_observed']['atomic'][field]
            assert atom['ids'] == [] and atom['max_abs_meV_cell'] is None
            assert not atom['complete_four_axial_points'] and atom['denominator'] == 4


def test_partial_statistics_keep_coverage_without_ranking_the_subset(actual):
    rows = deepcopy(actual['rows'][0]['points'])
    rows[2]['model_minus_observed']['atomic']['absolute_eV_cell'] = -.001
    stat = a.summarize(rows)['model_minus_observed']['atomic']['absolute_eV_cell']
    assert stat == dict(ids=[2], denominator=4, complete_four_axial_points=False, max_abs_meV_cell=1.)


def synthetic_prefix():
    report = json.loads(EMPTY.read_bytes())
    report['sacct'] = '28867177|RUNNING|0:0|10|32|node999|hfacnormal01'
    report['completed_nodes_audited'] = report['gradient_candidates'] = 1
    report['native_new_SCFs_reparsed'] = 1
    report['native_sources'] = [dict(explicitly_synthetic=True)]
    report['old_native_cache_preflight'] = dict(explicitly_synthetic=True)
    report['pending'] = dict(channel='T', id=1, reason='point_not_present')
    report['candidate_rows'] = [dict(channel='T', id=0, q_A=[0., 0.], reason=None,
        new_completed_SCFs=1, sum_completed_call_seconds=1.,
        replay=dict(exact_request_order_optimizer_replay=True),
        candidate=dict(energy_eV_cell=-1., atomic_released_gradient_norm_eV_A=.02,
            joint_released_gradient_norm_eV_A=1., released_stability_unchecked=True))]
    return report


def test_synthetic_atomic_candidate_is_not_rejected_for_large_joint_gradient():
    report = synthetic_prefix()
    out = a.validate_atomic(report, CONFIG)
    assert out[('T', 0)]['candidate']['joint_released_gradient_norm_eV_A'] == 1.


@pytest.mark.parametrize('fault', ['job', 'Q', 'bool_Q', 'duplicate', 'budget', 'gradient',
    'stability', 'replay', 'native_count', 'pending', 'complete'])
def test_mismatched_audit_or_false_candidate_rejected(fault):
    report = synthetic_prefix()
    row = report['candidate_rows'][0]
    if fault == 'job': report['job_id'] = 'other'
    elif fault == 'Q': row['q_A'] = [.1, 0.]
    elif fault == 'bool_Q': row['q_A'] = [False, 0.]
    elif fault == 'duplicate':
        report['candidate_rows'].append(deepcopy(row))
        report['completed_nodes_audited'] = 2
    elif fault == 'budget': row['new_completed_SCFs'] = 13
    elif fault == 'gradient': row['candidate']['atomic_released_gradient_norm_eV_A'] = .031
    elif fault == 'stability': row['candidate']['released_stability_unchecked'] = False
    elif fault == 'replay': row['replay']['exact_request_order_optimizer_replay'] = False
    elif fault == 'native_count': report['native_new_SCFs_reparsed'] = 2
    elif fault == 'pending': report['pending']['id'] = 2
    else: report['audit_mode'] = 'complete'
    with pytest.raises(ValueError): a.validate_atomic(report, CONFIG)


def test_changed_atomic_hash_or_bound_original_source_refused(monkeypatch):
    with pytest.raises(ValueError, match='atomic-audit byte hash'):
        a.analyze(EMPTY, '0' * 64)
    monkeypatch.setattr(a, 'PINS', {**a.PINS, 'G3_atomic_audit_E176_20261015/forecast.json': '0' * 64})
    with pytest.raises(ValueError, match='original source bytes'):
        a.analyze(EMPTY, EMPTY_SHA)


def test_existing_output_refused_before_any_analysis(tmp_path, monkeypatch):
    out = tmp_path / 'old.json'
    out.write_text('preserve')
    monkeypatch.setattr(a.sys, 'argv', ['analyze', '--atomic-audit', str(EMPTY),
        '--atomic-audit-sha256', EMPTY_SHA, '--output', str(out)])
    monkeypatch.setattr(a, 'analyze', lambda *args: pytest.fail('must refuse before reading sources'))
    with pytest.raises(FileExistsError): a.main()
    assert out.read_text() == 'preserve'
