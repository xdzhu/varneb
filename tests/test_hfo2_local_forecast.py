"""Check prospective elementary forecasts without any independent DFT label."""
from pathlib import Path
import importlib.util
import itertools
import json

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'benchmarks/hfo2_channels/20261008/G3_local_forecast_E185_20261015'
spec = importlib.util.spec_from_file_location('E185_freeze', CASE / 'freeze.py')
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)


@pytest.fixture(scope='module')
def forecast():
    return f.freeze()


def test_geometric_targets_are_eight_new_points_not_training_or_origin():
    targets = f.targets()
    assert len(targets) == 8 and len({tuple(p['q_A']) for p in targets}) == 8
    old = {tuple(p['q_A']) for p in f.faces.ordered_points()}
    assert not old.intersection(tuple(p['q_A']) for p in targets)
    assert sum(0 not in p['q_A'] for p in targets) == 4


def test_central_nine_have_identified_quadratic_and_only_audited_labels(forecast):
    t, p = forecast['channels']
    for row in (t, p):
        assert {tuple(x['original_index']) for x in row['training']} == set(itertools.product((-1,0,1), repeat=2))
        assert len(row['refitted_total_quadratic_coefficients']) == 6
    assert [x['historical_face_SCFs_at_nine_training_targets'] for x in (t,p)] == [28,36]
    assert [x['original_whole_training_face_SCFs'] for x in (t,p)] == [118,172]


def test_all_five_forecasts_both_steps_and_calibration_are_retained(forecast):
    assert forecast['scalar_predictions'] == 80
    for channel in forecast['channels']:
        for point in channel['targets']:
            assert point['independent_DFT_label'] is None
            assert set(point['models']) == set(f.MODELS)
            for model in point['models'].values():
                assert np.isfinite(model['value'])
                assert model['seen_origin_calibrated_energy_change_eV_cell'] == pytest.approx(
                    model['value'] - model['origin_bias_eV_cell'], abs=1e-14)
                assert model['relative_to_seen_origin_eV_cell'] == pytest.approx(
                    model['value'] - model['origin_model_energy_change_eV_cell'], abs=1e-14)


def test_start_uses_training_only_exact_ties_origin_and_original_guard(forecast):
    for row in forecast['channels']:
        data = f.faces.load_inputs(row['channel'])
        for point in row['targets']:
            assert point['training_parent_id'] == 0
            start = np.array(point['starting_internal_coordinates_A'])
            np.testing.assert_allclose(data['axes'].T @ start, point['q_A'], atol=1e-13, rtol=0)
            np.testing.assert_allclose(data['basis'] @ start, point['starting_full_coordinates_A'], atol=1e-14, rtol=0)
            assert f.faces.canary.geometry_guard(data['chart'], data['basis'], start)


def test_hessian_predictions_match_preexisting_complete_reference_polynomial(forecast):
    # Independent scalar oracle: E180 full-reference energies on the original
    # central nine nodes determine this quadratic exactly. No new labels.
    source = f.elementary.pinned(f.elementary.SOURCE, f.elementary.PIN)
    for row in forecast['channels']:
        for step in (.01, .02):
            original = next(r for r in source['rows'] if r['channel'] == row['channel']
                            and r['reference'] == 'T_all_optical' and r['step_A'] == step)
            labels = {tuple(2 * i for i in t['original_index']):
                original['points'][t['id']]['model_energy_change_eV_cell'] for t in row['training']}
            for point in row['targets']:
                oracle = f.elementary.interpolate(point['half_grid_index'], labels, 'tensor_quadratic')['value']
                assert point['models'][f'original_H_{step:.2f}']['value'] == pytest.approx(oracle, abs=1e-12)


def test_interpolant_prediction_uses_central_training_scale_not_old_outer_grid(forecast):
    for row in forecast['channels']:
        labels = {tuple(t['original_index']):t['energy_change_eV_cell'] for t in row['training']}
        negative_x = next(p for p in row['targets'] if p['q_A'] == [-.05,0.])
        linear = .5 * (labels[(-1,0)] + labels[(0,0)])
        tensor = .375 * labels[(-1,0)] + .75 * labels[(0,0)] - .125 * labels[(1,0)]
        assert negative_x['models']['bilinear']['value'] == pytest.approx(linear, abs=1e-14)
        assert negative_x['models']['tensor_quadratic']['value'] == pytest.approx(tensor, abs=1e-14)


def test_unsafe_starts_do_not_repair_or_use_test_labels(monkeypatch):
    monkeypatch.setattr(f.faces.canary, 'geometry_guard', lambda *args: False)
    data = dict(axes=np.eye(3)[:,:2], chart=None, basis=None)
    training = [dict(id=0, original_index=[0,0], coordinates_internal_A=[0.,0.,0.])]
    with pytest.raises(ValueError, match='no registered safe'):
        f.fixed_start(data, training, f.targets()[0])


def test_model_policies_do_not_claim_physics_or_full_validation(forecast):
    assert forecast['local_test_labels_seen'] is False
    assert forecast['physical_reserved_plus005_geometry_generated_or_read'] is False
    assert forecast['strong_T_Cmma_full_B0_B5_physical_strain_predictions_and_domain_freeze_complete'] is False
    assert forecast['whole_path_error_gain_stability_TS_or_JCTC_validation_claimed'] is False
    policy = forecast['test_execution_policy']
    assert policy['maximum_new_attempts_two_faces'] == 192
    assert policy['independent_points_computed_per_face'] == 0
    assert policy['maximum_internal_displacement_A'] == .5
    assert policy['gradient_tolerance_eV_A'] == .03
    assert forecast['scoring']['reporting_target_meV_12atom_cell'] / 4 == forecast['scoring']['reporting_target_meV_formula_unit']
    assert forecast['scoring']['coverage_denominator_per_face'] == 8


def test_frozen_report_matches_computation(forecast):
    frozen = json.loads((CASE / 'forecast_v2.json').read_bytes())
    assert forecast == frozen


def test_original_numbers_and_policies_remain_unchanged(forecast):
    first = f.elementary.pinned(CASE / 'forecast.json', f.FIRST_SHA)
    assert forecast['channels'] == first['channels']
    for key in ('scoring', 'test_execution_policy', 'models', 'scalar_predictions',
                'modeling_training_commitments_unchanged', 'branch_start_policy', 'selection'):
        assert forecast[key] == first[key]


def test_only_text_code_hash_normalizes_windows_line_endings(tmp_path):
    lf, crlf = tmp_path/'lf.py', tmp_path/'crlf.py'
    lf.write_bytes(b'print(1)\n')
    crlf.write_bytes(b'print(1)\r\n')
    assert f.lf_text_sha(lf) == f.lf_text_sha(crlf)
    # Native JSON/DFT result and binary pins continue to be exact byte hashes.
    with pytest.raises(ValueError, match='original source bytes differ'):
        f.elementary.pinned(crlf, f.lf_text_sha(crlf))
