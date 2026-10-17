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


def continuation_module():
    path = PATH.parent/'forecast_and_chain.py'
    s = importlib.util.spec_from_file_location('sampled_forecast', path)
    module = importlib.util.module_from_spec(s); s.loader.exec_module(module)
    return module


def test_all_registered_references_retained_when_physical_transfer_unavailable(tmp_path):
    module = continuation_module()
    definitions = module.model_definitions()
    assert len(definitions) == 69
    assert all(any(f'Cmma_{i}' in name for name in definitions) for i in (103, 106, 114, 127))
    assert len([name for name in definitions if name.startswith('B5_retained_anharmonic')]) == 2
    # Synthetic endpoint features test arithmetic/storage only; not native DFT evidence.
    f = dict(status='three_native_visible_endpoint_features_audited',
        registration_sha256=module.REGISTRATION_PIN, target_interiors_read=False,
        physical_input_sha256=study.CONTRACT, new_SCF_calls=3, SCF_core_hours=1.,
        allocation_core_hours=2., endpoints={p: dict(anchor_energy_eV_cell=-100.,
        target_energy_eV_cell=-99.96 if p != 'T' else -99.94) for p in study.PHASES})
    source = tmp_path/'synthetic_features.json'; study.save(source, f)
    result = module.freeze(ROOT, source, study.sha(source), tmp_path/'forecast')
    assert result['models'] == 69 and result['forecast_records'] == 138
    assert result['numeric_records'] == 6
    assert not result['target_interiors_read']
    assert result['baselines']['PO_to_T']['B1_fixed']['prediction_eV_fu'] == pytest.approx(.08441662388122495-.01)
    from vcneb.forecast_freeze import load_frozen_forecast_panel
    rows = load_frozen_forecast_panel(tmp_path/'forecast/frozen', expected_manifest_sha256=result['manifest_sha256'])
    assert len([r for r in rows if r.prediction_eV_fu is None and r.abstention_reason]) == 132
    assert all(r.abstention_reason is None for r in rows if r.model_id == 'B0')


def test_geometry_only_polygon_preserves_order_lift_endpoints_and_plane():
    from ase import Atoms
    from vcneb import clamped_plane_vcneb_boundary
    import numpy as np
    module = continuation_module()
    positions = [[x, y, z] for x in (1., 3., 5.) for y in (1., 3.) for z in (1., 3.)]
    start = Atoms(['Hf']*4+['O']*8, cell=np.eye(3)*8, positions=positions, pbc=True)
    images = [start.copy() for _ in range(9)]
    for i, atoms in enumerate(images): atoms.positions[4:, 0] += .01*i
    reference = start.cell.array.copy(); reference[:2] *= 1.02
    boundary = clamped_plane_vcneb_boundary(12, reference, allow_tilt=True)
    ends = [images[0].copy(), images[-1].copy()]
    for a in ends: a.set_cell(reference, scale_atoms=True)
    ends[1].positions[4:, 2] += .02
    result = module.deform_zero_band(images, ends, boundary)
    assert all(a.calc is None and a.get_chemical_symbols() == start.get_chemical_symbols() for a in result)
    assert all(np.array_equal(a.cell.array[:2], reference[:2]) for a in result)
    np.testing.assert_allclose(result[0].positions, ends[0].positions, atol=1e-12)
    np.testing.assert_allclose(result[-1].positions, ends[-1].positions, atol=1e-12)
    np.testing.assert_allclose(images[0].positions, start.positions, atol=1e-12)


def test_sampled_scores_keep_caps_and_model_abstentions_out_of_accuracy():
    path = PATH.parent/'audit_paths.py'
    s = importlib.util.spec_from_file_location('holdout_score', path)
    module = importlib.util.module_from_spec(s); s.loader.exec_module(module)
    from vcneb.prediction_evaluation import FrozenChannelForecast
    channels = module.CHANNELS
    forecasts = [FrozenChannelForecast('B0', c, p, None, 'a'*64, 'b'*64)
                 for c, p in zip(channels, (.11, .09))]
    forecasts += [FrozenChannelForecast('B5', c, None, 'physical transfer unavailable', 'a'*64, 'c'*64)
                  for c in channels]
    refs = {c: dict(ordinary_converged=True, sampled_peak_eV_fu=e)
            for c, e in zip(channels, (.111, .088))}
    rows = {r['model']: r for r in module.point_scores(forecasts, refs)}
    assert rows['B0']['maximum_absolute_error_meV_fu'] == pytest.approx(2.)
    assert rows['B0']['T_minus_preserving_error_meV_fu'] == pytest.approx(3.)
    assert rows['B5']['scored'] == 0 and rows['B5']['maximum_absolute_error_meV_fu'] is None
    refs[channels[1]]['ordinary_converged'] = False
    rows = {r['model']: r for r in module.point_scores(forecasts, refs)}
    assert rows['B0']['scored'] == 1 and rows['B0']['T_minus_preserving_error_meV_fu'] is None
    assert all(not r['physically_resolved_advantage_certified'] for r in rows.values())
