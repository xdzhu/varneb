"""Numerical and information-boundary checks for the separate sampled target."""
import copy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT/'benchmarks/hfo2_channels/20261008/physical_holdout_E268_20261018/study.py'
spec = importlib.util.spec_from_file_location('sampled_holdout', PATH)
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)


def network():
    return json.loads((ROOT/study.NETWORK).read_bytes())


def test_actual_training_only_pair_and_halfway_values():
    n = network()
    assert study.sha(ROOT/study.NETWORK) == study.NETWORK_SHA
    assert study.sampled_selection(n) == {'switching': 'preserving', 'decay': 'T'}
    rows = study.training_rows(n)
    for name, expected in [('preserving', 111.1263070804398), ('T', 90.59013713795138)]:
        values = study.simple_values(rows[0., name]['forward_meV_fu']/1000,
                                     rows[.01, name]['forward_meV_fu']/1000)
        assert values['B0']['prediction_eV_fu']*1000 == pytest.approx(expected)
        assert values['B1_fixed']['prediction_eV_fu'] is None
        assert values['B0']['physical_error_bound_eV_fu'] is None


def test_units_signed_abstention_and_common_initial_cancellation():
    a = study.simple_values(.10, .14, .02, .03, reaction=.01)
    b = study.simple_values(.08, .10, .02, -.01, reaction=-.02)
    assert a['B1_fixed']['prediction_eV_fu'] == pytest.approx(.08)
    assert a['B1_follow']['prediction_eV_fu'] == pytest.approx(.11)
    assert (a['B1_fixed']['prediction_eV_fu']-b['B1_fixed']['prediction_eV_fu']) == pytest.approx(.02)
    invalid = study.simple_values(.01, .02, .03, -.02, reaction=.005)
    assert invalid['B1_fixed']['raw_signed_prediction_eV_fu'] == pytest.approx(-.02)
    assert invalid['B1_fixed']['prediction_eV_fu'] is None
    assert invalid['B1_fixed']['abstention_reason'] == 'below_endpoint_inclusive_peak_limit'


@pytest.mark.parametrize('change', ['duplicate', 'force', 'profile', 'bound'])
def test_incomplete_or_fabricated_training_rejected(change):
    n = copy.deepcopy(network())
    if change == 'duplicate': n['rows'][-1] = n['rows'][0]
    if change == 'force': n['rows'][0]['ordinary_fmax_eV_A'] = .2
    if change == 'profile': n['rows'][0]['forward_meV_fu'] += 1
    if change == 'bound': n['rows'][0]['barrier_error_bound_meV_fu'] = .1
    with pytest.raises(ValueError): study.training_rows(n)


def test_registration_no_holdout_geometry_and_no_overwrite(tmp_path):
    dest = tmp_path/'registration.json'
    r = study.register(ROOT, dest)
    assert not r['target_structures_generated'] and r['target_DFT_calls'] == 0
    assert not r['physical_ordering_certified']
    assert len(r['unprovided_physical_transfers']) == 4
    assert all(len(x['representations']) == 8 for x in r['unprovided_physical_transfers'].values())
    assert r['original_2D_preserving_refusal_retained']
    original = dest.read_bytes()
    with pytest.raises(FileExistsError): study.register(ROOT, dest)
    assert dest.read_bytes() == original
