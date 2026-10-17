"""Complete archived training replay; not independent labels or a new DFT run."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
TRAINING = ROOT/'benchmarks/hfo2_channels/20261008/G3_preserving_branches_E260_20261017'
TESTS = ROOT/'benchmarks/hfo2_channels/20261008/G3_P_branch_independent_E263_20261018'
AUDIT_SHA = '30971f61ab7d12d8e9bfa729a517603566c20cefa6bf442b06a4100264bd8811'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


@pytest.fixture(scope='module')
def native():
    raw = (TRAINING/'native_complete_E264.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest() == AUDIT_SHA
    return json.loads(raw)


@pytest.fixture(scope='module')
def analyzed():
    return load('E264_actual_training_replay_test', TRAINING/'analyze_training_E264.py').analyze()


def test_all_original_ten_calls_caps_and_real_allocation_cost(native):
    assert [r['id'] for r in native['rows']] == list(range(10))
    assert [r['attempted_SCFs'] for r in native['rows']] == [12, 12, 6, 6, 7, 7, 12, 7, 6, 6]
    assert native['new_native_SCFs_reparsed'] == len(native['native_records']) == 81
    assert native['gradient_candidates'] == 7
    assert native['allocation_core_hours'] == 85.6
    assert all(r['validation']['genuine_MPI_ranks'] == 32 for r in native['native_records'])
    assert all(r['minimum_distance_A'] >= 1.6 for r in native['native_records'])


def test_all_native_records_consumed_by_unchanged_optimizer_not_a_nearest_cache(analyzed):
    assert analyzed['native_SCFs_consumed'] == 81
    assert all(r['maximum_replay_coordinate_arithmetic_difference_A'] <= 2e-12 for r in analyzed['rows'])
    assert analyzed['new_DFT_or_scheduler_operations'] == 0
    assert not analyzed['stable_branches_or_TS_or_predictive_gain_certified']


def test_different_capped_failure_mechanisms_are_not_conflated(analyzed):
    first, second, other = [analyzed['rows'][i] for i in (0, 1, 6)]
    assert first['radius_margin_A'] < 1e-5 and second['radius_margin_A'] < 1e-6
    assert first['tangent_gradient_norm_eV_A'] > .29 and second['tangent_gradient_norm_eV_A'] > .29
    assert other['radius_margin_A'] > .1 and other['released_gradient_norm_eV_A'] < .03
    assert abs(other['zeta_gradient_eV_A']) > 1e-4
    assert other['energy_rejected_native_calls'] == 5
    assert [r['id'] for r in analyzed['rows'] if r['replay_status'] == 'budget_exception_replayed'] == [0, 1, 6]


def test_three_two_sided_pairs_do_not_establish_distinct_stable_branches(native):
    paired = [r for r in native['pairs'] if r['both_gradient_candidates']]
    assert [r['Q_id'] for r in paired] == [1, 2, 4]
    assert all(abs(r['plus_minus_energy_meV_per_fu']) < .001 for r in paired)
    assert all(not r['stable_distinct_branches_certified'] for r in native['pairs'])


def test_eight_real_training_only_starts_respect_fixed_Q_without_test_labels():
    module = load('E264_actual_start_preparation_test', TESTS/'prepare_starts.py')
    result = module.prepare(TRAINING/'native_complete_E264.json', AUDIT_SHA)
    saved = json.loads((TESTS/'starts_E264.json').read_bytes())
    assert result == saved
    assert [r['training_parent_id'] for r in result['points']] == [2, 4, 7, 8, 2, 2, 4, 7]
    assert all(r['independent_DFT_label'] is None for r in result['points'])
    assert all(np.linalg.norm(r['starting_internal_coordinates_A']) <= .5 for r in result['points'])
    assert result['maximum_new_P_attempts'] == 96 and result['remaining_P_tests'] == 8
    assert not result['final_models_domains_and_predictions_frozen']
    assert not result['authorizes_independent_DFT']


def test_wrong_native_audit_pin_refuses_before_materialization():
    module = load('E264_wrong_start_audit_pin', TESTS/'prepare_starts.py')
    with pytest.raises(ValueError, match='externally pinned'):
        module.prepare(TRAINING/'native_complete_E264.json', '0'*64)
