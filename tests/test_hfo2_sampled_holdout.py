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


def test_prepared_endpoints_keep_one_substrate_and_no_cached_results():
    from scripts.relax_clamped_ase_endpoint import load_seed
    import numpy as np
    case = ROOT/study.CASE
    m = json.loads((case/'endpoint_seeds/seeds.json').read_bytes())
    assert len(m['entries']) == 3 and not m['target_interiors_generated']
    planes = []
    for row in m['entries']:
        path = case/'endpoint_seeds'/row['phase']/'endpoint_seed.json'
        assert study.sha(path) == row['manifest_sha256']
        atoms, boundary, seed = load_seed(path)
        assert seed['strain'] == .005 and not seed['old_endpoint_results_reused']
        assert atoms.calc is None and atoms.get_chemical_symbols() == ['Hf']*4+['O']*8
        boundary.validate_images([atoms])
        planes.append(atoms.cell.array[:2])
    assert all(np.array_equal(planes[0], p) for p in planes[1:])
    original = study.training_rows(network())[0., 'T']['fixed_substrate_rows_A']
    np.testing.assert_allclose(planes[0], np.asarray(original)*1.005, atol=1e-12, rtol=0)


def test_submission_guard_checks_real_allocation_and_exact_handle():
    path = PATH.parent/'dispatch_endpoints.py'
    s = importlib.util.spec_from_file_location('endpoint_dispatch', path)
    module = importlib.util.module_from_spec(s); s.loader.exec_module(module)
    valid = ('JobId=123 JobState=PENDING Reason=JobHeldUser Partition=hfacnormal01 '
             'NumCPUs=32 NumTasks=1 CPUs/Task=32 NumNodes=1-1 TimeLimit=08:00:00 '
             'JobName=hfo2-E268-endpoints-0')
    module.validate_held(valid, 0, '123')
    for bad in (valid.replace('NumCPUs=32', 'NumCPUs=1'),
                valid.replace('hfacnormal01', 'other'), valid.replace('JobId=123', 'JobId=124'),
                valid.replace('NumNodes=1-1', 'NumNodes=2')):
        with pytest.raises(ValueError): module.validate_held(bad, 0, '123')
