"""Complete actual atomic treatment cohort; no stable-branch or TS promotion."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'benchmarks/hfo2_channels/20261008'
CASE = BASE / 'G3_atomic_complete_E194_20261015'
PREFIX = BASE / 'G3_atomic_treatment_E192_20261015'
AUDIT_SHA = '8a2fdf2a2beec306bb2f39aed8b347019a4f049b2d8cbc42cd79de1fce922b92'
REPORT_SHA = '3cc150ec3a6d0dbfa5f5879f7215a3692524309add3e5bb7f5abed384a197444'
spec = importlib.util.spec_from_file_location(
    'complete_atomic_comparison', BASE / 'G3_treatment_comparison_E178_20261015/analyze.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


@pytest.fixture(scope='module')
def actual():
    return json.loads((CASE / 'analysis.json').read_bytes())


@pytest.fixture(scope='module')
def native():
    return json.loads((CASE / 'atomic_native_complete.json').read_bytes())


def test_complete_audit_hash_and_existing_postprocessor_reproduce(actual, native):
    assert a.sha(CASE / 'atomic_native_complete.json') == AUDIT_SHA
    assert a.sha(CASE / 'analysis.json') == REPORT_SHA
    assert a.analyze(CASE / 'atomic_native_complete.json', AUDIT_SHA) == actual
    assert native['audit_mode'] == 'complete'
    assert native['sacct'] == '28867177|COMPLETED|0:0|4577|32|node6|hfacnormal01'
    assert native['completed_nodes_audited'] == native['gradient_candidates'] == 10
    assert native['native_new_SCFs_reparsed'] == len(native['native_sources']) == 36
    assert native['all10points_covered'] and native['pending'] is None


def test_previous_complete_nodes_and_native_results_are_preserved(native):
    old = json.loads((PREFIX / 'atomic_native_prefix.json').read_bytes())
    assert native['candidate_rows'][:7] == old['candidate_rows']
    assert native['native_sources'][:20] == old['native_sources']
    assert [(r['channel'], r['id']) for r in native['candidate_rows']] == (
        [(c, i) for c in ('T', 'preserving') for i in range(5)])
    assert all(r['replay']['exact_request_order_optimizer_replay']
               for r in native['candidate_rows'])


def test_all_atomic_gradient_candidates_are_not_joint_stationary(native):
    for row in native['candidate_rows']:
        p = row['candidate']
        assert 0 <= p['atomic_released_gradient_norm_eV_A'] <= .03
        assert p['joint_released_gradient_norm_eV_A'] > .03
        assert p['released_stability_unchecked']
    assert not native['conditional_stability_TS_acceleration_or_independent_prediction_certified']


@pytest.mark.parametrize('step,expected', [
    (.01, (.075880169, 4.349823830, 4.671403795)),
    (.02, (.116535481, 4.399770177, 4.720758460))])
def test_full_preserving_axial_cohort_discrepancies(actual, step, expected):
    row = next(r for r in actual['rows'] if r['channel'] == 'preserving' and r['probe_step_A'] == step)
    stats = row['non_origin_coverage_statistics']['model_minus_observed']
    for treatment, value in zip(('frozen', 'atomic', 'joint'), expected):
        s = stats[treatment]['shape_eV_cell']
        assert s['ids'] == [1, 2, 3, 4] and s['denominator'] == 4
        assert s['complete_four_axial_points']
        assert s['max_abs_meV_cell'] == pytest.approx(value, rel=0, abs=2e-9)


def test_all_original_forecasts_and_known_F_J_values_unchanged(actual):
    old = json.loads((PREFIX / 'analysis.json').read_bytes())
    assert len(actual['rows']) == 4
    for row, before in zip(actual['rows'], old['rows']):
        assert (row['channel'], row['probe_step_A']) == (before['channel'], before['probe_step_A'])
        for p, previous in zip(row['points'], before['points']):
            assert (p['id'], p['q_A'], p['model']) == (previous['id'], previous['q_A'], previous['model'])
            assert p['observed']['frozen'] == previous['observed']['frozen']
            assert p['observed']['joint'] == previous['observed']['joint']
            assert p['observed']['atomic'] is not None
            assert p['atomic_status'] == 'gradient_candidate_stability_unchecked'


def test_five_preserving_release_values_and_pointwise_decompositions(actual):
    row = next(r for r in actual['rows'] if r['channel'] == 'preserving')
    expected = [(-48.348691, -1.136239), (-44.073849, -.804122),
                (-47.241146, -1.006188), (-47.238928, -1.006405), (-44.063824, -.816911)]
    for p, (atomic, joint) in zip(row['points'], expected):
        d = p['observed_changes']
        assert 1000*d['atomic_minus_frozen']['absolute_eV_cell'] == pytest.approx(atomic, rel=0, abs=5.1e-7)
        assert 1000*d['joint_minus_atomic']['absolute_eV_cell'] == pytest.approx(joint, rel=0, abs=5.1e-7)
    for row in actual['rows']:
        for p in row['points']:
            for group in ('observed_changes', 'model_changes', 'missed_treatment_change'):
                for field in ('absolute_eV_cell', 'shape_eV_cell'):
                    d = p[group]
                    assert d['atomic_minus_frozen'][field] + d['joint_minus_atomic'][field] == (
                        pytest.approx(d['joint_minus_frozen'][field], rel=0, abs=2e-12))


def test_complete_cost_and_prediction_status_are_not_acceleration_or_stability(actual, native):
    assert [sum(r['new_completed_SCFs'] for r in native['candidate_rows'] if r['channel'] == c)
            for c in ('T', 'preserving')] == [9, 27]
    costs = actual['costs']
    assert costs['atomic_new_SCFs'] == 36 and costs['atomic_completed_nodes'] == 10
    assert costs['atomic_completed_call_proxy_core_hours'] == pytest.approx(39.74955780385269)
    assert costs['original_curvature_SCFs_additional'] == 296
    assert costs['not_total_campaign_cost_or_speedup']
    assert actual['formula_units'] == 4 and actual['atomic_audit_mode'] == 'complete'
    assert actual['fitted_coefficients'] == actual['independent_slots_used'] == 0
    assert actual['new_DFT_or_scheduler_writes'] == 0
    assert not actual['physical_reserved_geometry_generated_or_read']
    assert not actual['unique_mode_energy_TS_conditional_stability_or_predictive_advantage_certified']
