"""Actual same-Q treatment results; the supplied prefix is not a full batch."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'benchmarks/hfo2_channels/20261008'
CASE = BASE / 'G3_atomic_treatment_E192_20261015'
PREFIX_SHA = '673714cb35951ee60b856a1bfb7775aba91ab6ee8a3c2cd60bbe7d8c50dd2e96'
REPORT_SHA = '16b9dde68203a9888a642959c2a6ce979cd448f9ceb926a0cec911d06d61d5b5'
spec = importlib.util.spec_from_file_location(
    'actual_atomic_comparison', BASE / 'G3_treatment_comparison_E178_20261015/analyze.py')
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


@pytest.fixture(scope='module')
def actual():
    return json.loads((CASE / 'analysis.json').read_bytes())


@pytest.fixture(scope='module')
def native():
    return json.loads((CASE / 'atomic_native_prefix.json').read_bytes())


def test_exact_native_prefix_and_reproduced_report(actual):
    assert comparison.sha(CASE / 'atomic_native_prefix.json') == PREFIX_SHA
    assert comparison.sha(CASE / 'analysis.json') == REPORT_SHA
    assert comparison.analyze(CASE / 'atomic_native_prefix.json', PREFIX_SHA) == actual
    forecast = BASE / 'G3_atomic_audit_E176_20261015/forecast.json'
    assert hashlib.sha256(forecast.read_bytes()).hexdigest() == (
        '321ba7007c7f950166aab4eaa2dbe579773164968f04a94034a8a942bee1318d')


def test_actual_audit_retains_order_replay_and_incomplete_job(native):
    assert native['audit_mode'] == 'prefix' and '|RUNNING|' in native['sacct']
    assert not native['all10points_covered']
    assert native['native_new_SCFs_reparsed'] == len(native['native_sources']) == 20
    assert native['gradient_candidates'] == native['completed_nodes_audited'] == 7
    assert [(r['channel'], r['id']) for r in native['candidate_rows']] == (
        [('T', i) for i in range(5)] + [('preserving', 0), ('preserving', 1)])
    assert native['pending'] == dict(channel='preserving', id=2, reason='point_not_present')
    assert all(r['replay']['exact_request_order_optimizer_replay']
               for r in native['candidate_rows'])


def test_atomic_candidates_do_not_pass_joint_stationarity(native):
    for row in native['candidate_rows']:
        p = row['candidate']
        assert 0 <= p['atomic_released_gradient_norm_eV_A'] <= .03
        assert p['joint_released_gradient_norm_eV_A'] > .03
        assert p['released_stability_unchecked']


@pytest.mark.parametrize('channel,atomic,joint', [
    ('T', -.27139817029819824, -.6521419600176159),
    ('preserving', -48.348690623242874, -49.4849300084752)])
def test_origin_treatment_differences_are_signed_not_unique_partitions(actual, channel, atomic, joint):
    row = next(r for r in actual['rows'] if r['channel'] == channel)
    origins = row['observed_origins']
    assert origins['frozen'] == 0
    assert origins['atomic'] * 1000 == pytest.approx(atomic, rel=0, abs=2e-9)
    assert origins['joint'] * 1000 == pytest.approx(joint, rel=0, abs=2e-9)
    origin = row['points'][0]['observed_changes']
    a = origin['atomic_minus_frozen']['absolute_eV_cell']
    b = origin['joint_minus_atomic']['absolute_eV_cell']
    c = origin['joint_minus_frozen']['absolute_eV_cell']
    assert a + b == pytest.approx(c, rel=0, abs=2e-12)
    assert a < 0 and b < 0


@pytest.mark.parametrize('step,expected', [
    (.01, (.4162635209878788, .5018789359126777, 1.3271772093247154)),
    (.02, (.4031595155926308, .5157039519863094, 1.3460706421835107))])
def test_same_four_T_axial_points_have_measured_three_treatment_errors(actual, step, expected):
    row = next(r for r in actual['rows'] if r['channel'] == 'T' and r['probe_step_A'] == step)
    stats = row['non_origin_coverage_statistics']['model_minus_observed']
    for treatment, value in zip(('frozen', 'atomic', 'joint'), expected):
        s = stats[treatment]['shape_eV_cell']
        assert s['ids'] == [1, 2, 3, 4] and s['denominator'] == 4
        assert s['complete_four_axial_points']
        assert s['max_abs_meV_cell'] == pytest.approx(value, rel=0, abs=2e-9)


def test_preserving_pending_labels_abstain_without_hiding_joint_data(actual):
    for row in actual['rows']:
        if row['channel'] != 'preserving':
            continue
        stats = row['non_origin_coverage_statistics']['model_minus_observed']
        assert stats['atomic']['shape_eV_cell']['ids'] == [1]
        assert not stats['atomic']['shape_eV_cell']['complete_four_axial_points']
        assert stats['joint']['shape_eV_cell']['ids'] == [1, 2, 3, 4]
        for p in row['points'][2:]:
            assert p['atomic_status'] == 'not_processed'
            assert p['observed']['atomic'] is None
            assert p['model_minus_observed']['atomic']['shape_eV_cell'] is None
            assert p['observed_changes']['joint_minus_frozen']['shape_eV_cell'] is not None


def test_all_frozen_model_predictions_and_pointwise_error_identities_remain(actual):
    previous = json.loads((BASE / 'G3_treatment_comparison_E178_20261015/analysis.json').read_bytes())
    assert len(actual['rows']) == 4
    for row, old in zip(actual['rows'], previous['rows']):
        assert len(row['points']) == 5
        for p, before in zip(row['points'], old['points']):
            assert (p['id'], p['q_A'], p['model']) == (before['id'], before['q_A'], before['model'])
            for field in ('absolute_eV_cell', 'shape_eV_cell'):
                missed = p['missed_treatment_change']
                a, b, c = [missed[k][field] for k in (
                    'atomic_minus_frozen', 'joint_minus_atomic', 'joint_minus_frozen')]
                if a is not None and b is not None:
                    assert a + b == pytest.approx(c, rel=0, abs=2e-12)


def test_audited_cost_scope_is_not_later_inflight_work_or_speedup(actual, native):
    assert sum(r['new_completed_SCFs'] for r in native['candidate_rows'] if r['channel'] == 'T') == 9
    assert sum(r['new_completed_SCFs'] for r in native['candidate_rows'] if r['channel'] == 'preserving') == 11
    costs = actual['costs']
    assert costs['atomic_new_SCFs'] == 20 and costs['atomic_completed_nodes'] == 7
    assert costs['atomic_completed_call_proxy_core_hours'] == pytest.approx(22.331631473484965)
    assert costs['original_curvature_SCFs_additional'] == 296
    assert costs['not_total_campaign_cost_or_speedup']
    assert actual['new_DFT_or_scheduler_writes'] == actual['independent_slots_used'] == 0
    assert actual['fitted_coefficients'] == 0 and actual['no_step_or_point_selection']
    assert not actual['physical_reserved_geometry_generated_or_read']
    assert not actual['unique_mode_energy_TS_conditional_stability_or_predictive_advantage_certified']
