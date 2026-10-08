"""Synthetic mathematical controls, not strained HfO2 material results."""

import copy
import hashlib
import json
from pathlib import Path

from ase import Atoms
import numpy as np
import pytest

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.freeze_hfo2_prediction_controls import REQUIRED, freeze
from vcneb.channel_competition import ChannelPath, summarize_competing_paths
from vcneb.prediction_controls import simple_prediction_controls


def network(strain, *, swap=False, overlapping=False, pressure=0.):
    names = ["PO_flip_T_pattern_preserving", "PO_flip_T_pattern_reversing", "PO_to_M", "PO_to_T"]
    barriers = [.20, .30, .40, .50] if not swap else [.35, .15, .55, .25]
    initial = -100. + 10*strain
    paths = []
    for i, (name, barrier) in enumerate(zip(names, barriers)):
        images = [Atoms('Ar', scaled_positions=[[.2+.05*j, .3, .4]], cell=np.eye(3)*5, pbc=True) for j in range(3)]
        paths.append(ChannelPath(name=name, role=REQUIRED[name], images=images,
                                 energies_eV_cell=[initial, initial+4*barrier, initial-.4],
                                 physical_contract=CONTRACT, mechanical_family='same_substrate_tilt_open',
                                 mechanical_parameters={'biaxial_strain': strain, 'external_field_V_A': 0., 'allow_tilt': True},
                                 pressure_eV_A3=pressure, neb_fmax_eV_A=.08, source_id=f'synthetic_training_{strain}_{i}',
                                 barrier_error_eV_cell=.4 if overlapping else .004,
                                 sampling_audit_sha256='a'*64, barrier_error_audit_sha256='b'*64))
    return summarize_competing_paths(paths, formula_units=4, required_channels=REQUIRED)


def features(lower):
    def record(anchor, delta):
        return {'anchor_energy_eV_cell': anchor, 'target_energy_eV_cell': anchor+delta,
                'anchor_audit_sha256': 'c'*64, 'target_audit_sha256': 'd'*64,
                'target_endpoint_DFT_calls': 3}
    energy = lower['common_initial_raw_energy_eV_cell']
    return {'physical_contract': CONTRACT.copy(), 'mechanical_family': lower['mechanical_family'],
            'formula_units': 4, 'pressure_eV_A3': 0.,
            'anchor_parameters': lower['mechanical_parameters'].copy(),
            'target_parameters': dict(lower['mechanical_parameters'], biaxial_strain=.005),
            'initial': record(energy, .08),
            'final_by_channel': {r['name']: record(energy+4*r['reaction_energy_from_common_initial_eV_fu'], .04)
                                 for r in lower['channels']}}


def predict(a, b, f=None):
    return simple_prediction_controls(a, b, parameter='biaxial_strain', target=.005, endpoint_features=f)


def test_training_only_selection_keeps_lower_representative_if_upper_swaps():
    result = predict(network(0.), network(.01, swap=True))
    switching, decay = result['predictions']
    assert switching['selected_channel'] == 'PO_flip_T_pattern_preserving'
    assert switching['upper_training_selection']['selected'] == 'PO_flip_T_pattern_reversing'
    assert not switching['same_representative_at_both_conditions']
    assert decay['selected_channel'] == 'PO_to_M'
    assert decay['upper_training_selection']['selected'] == 'PO_to_T'
    assert switching['B0_direct_barrier_interpolation']['barrier_prediction_eV_fu'] == pytest.approx(.275)
    assert decay['B0_direct_barrier_interpolation']['barrier_prediction_eV_fu'] == pytest.approx(.475)
    assert switching['B1_fixed_absolute_bottleneck']['status'] == 'unavailable_endpoint_features_missing'
    assert result['visible_target_endpoint_DFT_calls'] == 0


def test_reference_interval_is_not_a_model_error_bound_and_no_input_mutates():
    a, b = network(0.), network(.01, swap=True)
    original = copy.deepcopy([a, b])
    result = predict(a, b)
    assert [a, b] == original
    b0 = result['predictions'][0]['B0_direct_barrier_interpolation']
    assert b0['propagated_training_reference_interval_eV_fu'] == pytest.approx([.274, .276])
    assert b0['prediction_error_bound_eV_fu'] is None
    assert not result['target_path_labels_read_by_this_API']
    assert result['new_DFT_calls_by_this_analysis'] == 0


def test_endpoint_response_sign_formula_normalization_and_visible_cost():
    a, b = network(0.), network(.01)
    f = features(a)
    originals = copy.deepcopy([a, b, f])
    result = predict(a, b, f)
    row = result['predictions'][0]
    assert row['B1_fixed_absolute_bottleneck']['barrier_prediction_eV_fu'] == pytest.approx(.18)
    assert row['B1_bottleneck_follows_final']['barrier_prediction_eV_fu'] == pytest.approx(.19)
    assert result['visible_target_endpoint_DFT_calls'] == 15  # not only selected features' cost
    assert [a, b, f] == originals
    result['endpoint_features']['initial']['target_energy_eV_cell'] = 0
    assert [a, b, f] == originals


def test_partial_final_features_do_not_manufacture_a_missing_baseline():
    a, b = network(0.), network(.01)
    f = features(a); f['final_by_channel'] = {}
    result = predict(a, b, f)
    assert result['visible_target_endpoint_DFT_calls'] == 3
    row = result['predictions'][0]
    assert row['B1_fixed_absolute_bottleneck']['barrier_prediction_eV_fu'] == pytest.approx(.18)
    assert row['B1_bottleneck_follows_final']['status'] == 'unavailable_endpoint_features_missing'


def test_finite_pressure_B0_uses_H_barriers_B1_cannot_silently_omit_P_deltaV():
    # Synthetic images have constant volume, so their H barriers match E.
    # The unseen endpoint volume is not supplied: E shifts cannot stand for H.
    a,b=network(0.,pressure=.001),network(.01,pressure=.001)
    result=predict(a,b)
    assert result['pressure_eV_A3']==.001
    assert result['predictions'][0]['B0_direct_barrier_interpolation']['barrier_prediction_eV_fu']==pytest.approx(.20)
    f=features(a); f['pressure_eV_A3']=.001
    with pytest.raises(ValueError,match='finite pressure needs audited volumes/enthalpy'):
        predict(a,b,f)


def test_interval_overlap_uses_declared_lexical_representative_not_best_posthoc_channel():
    result = predict(network(0., overlapping=True), network(.01))
    s = result['predictions'][0]['lower_training_selection']
    assert len(s['possible_lowest_within_training_bounds']) == 2
    assert not s['unique_lowest_resolved']
    assert s['selected'] == sorted(s['possible_lowest_within_training_bounds'])[0]


@pytest.mark.parametrize('fault', ['unconverged','missing_sampling','missing_error','bad_profile','bad_interval','duplicate','mixed_contract','mixed_family','different_field','empty_family','missing_residual'])
def test_bad_training_cannot_be_used_as_a_prediction(fault):
    a, b = network(0.), network(.01)
    if fault == 'unconverged': b['channels'][0]['source_NEB_fmax_eV_A'] = .11
    elif fault == 'missing_sampling': b['channels'][0]['sampling_audit_sha256'] = None
    elif fault == 'missing_error': b['channels'][0]['barrier_interval_eV_fu'] = None
    elif fault == 'bad_profile': b['channels'][0]['relative_profile_from_common_initial_eV_fu'][1] += .01
    elif fault == 'bad_interval': b['channels'][0]['barrier_interval_eV_fu'] = [.30, .40]
    elif fault == 'duplicate': b['channels'][1] = b['channels'][0]
    elif fault == 'mixed_contract': b['physical_contract']['INPUT'] = 'e'*64
    elif fault == 'mixed_family': b['mechanical_family'] = 'free_cell'
    elif fault == 'empty_family': b['mechanical_family'] = ''
    elif fault == 'missing_residual': del b['channels'][0]['source_NEB_fmax_eV_A']
    else: b['mechanical_parameters']['external_field_V_A'] = .1
    with pytest.raises(ValueError): predict(a, b)


@pytest.mark.parametrize('fault', ['path_label','wrong_condition','wrong_anchor','wrong_final','missing_audit','zero_cost','nan','wrong_pressure'])
def test_endpoint_features_are_not_target_path_labels_or_free_information(fault):
    a, b = network(0.), network(.01)
    f = features(a)
    if fault == 'path_label': f['holdout_path_barrier_eV_fu'] = .19
    elif fault == 'wrong_condition': f['target_parameters']['biaxial_strain'] = .006
    elif fault == 'wrong_anchor': f['initial']['anchor_energy_eV_cell'] += 1.
    elif fault == 'wrong_final': f['final_by_channel']['PO_to_M']['anchor_energy_eV_cell'] += 1.
    elif fault == 'missing_audit': f['initial']['target_audit_sha256'] = None
    elif fault == 'zero_cost': f['initial']['target_endpoint_DFT_calls'] = 0
    elif fault == 'nan': f['initial']['target_energy_eV_cell'] = float('nan')
    else: f['pressure_eV_A3'] = False
    with pytest.raises(ValueError): predict(a, b, f)


def test_unphysical_endpoint_ordering_abstains_without_clipping_or_erasing_prediction():
    a, b = network(0.), network(.01)
    f = features(a)
    f['initial']['target_energy_eV_cell'] += 2.
    row = predict(a, b, f)['predictions'][0]['B1_fixed_absolute_bottleneck']
    assert row['raw_signed_prediction_eV_fu'] == pytest.approx(-.32)
    assert row['barrier_prediction_eV_fu'] is None
    assert row['status'] == 'abstain_unphysical_endpoint_ordering'
    f = features(a)
    f['final_by_channel']['PO_flip_T_pattern_preserving']['target_energy_eV_cell'] += 2.
    row = predict(a, b, f)['predictions'][0]['B0_direct_barrier_interpolation']
    assert row['raw_signed_prediction_eV_fu'] == pytest.approx(.20)
    assert row['barrier_prediction_eV_fu'] is None


@pytest.mark.parametrize('target', [0.,.01,.02,float('nan'),True])
def test_no_extrapolation_or_seen_condition_claimed_as_unseen(target):
    with pytest.raises(ValueError):
        simple_prediction_controls(network(0.),network(.01),parameter='biaxial_strain',target=target)


def test_freeze_is_fresh_hashed_not_a_claim_full_prediction_batch_is_ready(tmp_path):
    files = [tmp_path/'zero.json',tmp_path/'one.json']
    for p,x in zip(files,(0.,.01)): p.write_text(json.dumps(network(x)))
    with pytest.raises(ValueError): freeze(*files,tmp_path/'denied')
    assert not (tmp_path/'denied').exists()
    result = freeze(*files,tmp_path/'prediction',holdout_path_labels_unread=True)
    assert not result['complete_forecast_batch_ready_for_holdout_paths']
    assert not result['independent_label_blinding_proven']
    receipt=json.loads((tmp_path/'prediction/freeze_receipt.json').read_text())
    assert receipt['prediction_sha256']==hashlib.sha256((tmp_path/'prediction/prediction_controls.json').read_bytes()).hexdigest()
    assert result['training_source_sha256']==[hashlib.sha256(p.read_bytes()).hexdigest() for p in files]
    with pytest.raises(FileExistsError): freeze(*files,tmp_path/'prediction',holdout_path_labels_unread=True)


def test_actual_G1_observations_are_rejected_before_output_and_not_synthetic_training(tmp_path):
    case=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/network_update_20261009/analysis.json'
    report=json.loads(case.read_text())
    with pytest.raises(ValueError,match='incomplete network'): predict(report,report)
    with pytest.raises(ValueError,match='G1 free cell is not training'):
        freeze(case,case,tmp_path/'not_a_forecast',holdout_path_labels_unread=True)
    assert not (tmp_path/'not_a_forecast').exists()
