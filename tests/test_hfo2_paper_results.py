"""Actual archived table replay and safe local export, no calculators."""
import csv
from copy import deepcopy
from pathlib import Path
import pytest
from examples.hfo2_paper_results import export, replay, load_sources, preserving_tables, SOURCES, ROOT


def test_all_paths_common_units_and_negative_material_signature_retained():
    result = replay()
    assert len(result['paths']) == 8
    values = {(r['channel'], r['substrate_strain']): r['sampled_forward_meV_fu'] for r in result['paths']}
    assert values['T', 0.] == pytest.approx(84.41662388122495)
    assert values['preserving', .01]-values['preserving', 0.] == pytest.approx(28.5962481552815)
    assert all(r['barrier_error_bound_meV_fu'] is None for r in result['paths'])
    assert not result['stable_TS_or_error_bound_or_method_advantage_certified']


def test_all_nineteen_models_raw_and_aligned_errors_and_abstentions_recomputed():
    result = replay(); models = {r['model']: r for r in result['local_models']}
    assert len(models) == 19 and len(result['local_points']) == 19*8
    assert models['original_H_0.01']['uncalibrated_max_error_meV_fu'] == pytest.approx(.07606453515501332)
    assert models['original_H_0.01']['calibrated_max_error_meV_fu'] == pytest.approx(.12976860140672872)
    assert any(r['supported_predictions'] < r['denominator'] for r in models.values())
    assert result['local_test_native_SCFs'] == 18
    assert not result['physical_plus005_results_included']


def test_fresh_export_is_portable_and_contains_all_rows(tmp_path):
    output = tmp_path/'result'; export(output)
    assert {p.name for p in output.iterdir()} == {'channel_comparison.csv', 'T_local_model_errors.csv', 'T_local_all_predictions.csv',
        'P_local_candidates.csv', 'P_local_model_errors.csv', 'P_local_all_predictions.csv', 'replay.json'}
    with (output/'T_local_all_predictions.csv').open(encoding='utf-8', newline='') as handle:
        assert len(list(csv.DictReader(handle))) == 152
    with (output/'P_local_all_predictions.csv').open(encoding='utf-8', newline='') as handle:
        assert len(list(csv.DictReader(handle))) == 184
    before = (output/'replay.json').read_bytes()
    with pytest.raises(FileExistsError): export(output)
    assert (output/'replay.json').read_bytes() == before


def test_changed_report_bytes_refused_before_output_creation(tmp_path):
    relative = SOURCES['network'][0]; path = tmp_path/relative; path.parent.mkdir(parents=True)
    path.write_bytes(b'{}\n')
    with pytest.raises(ValueError, match='source bytes changed'): load_sources(tmp_path)
    with pytest.raises(ValueError): export(tmp_path/'must-not-exist', root=tmp_path)
    assert not (tmp_path/'must-not-exist').exists()


def test_all_revised_P8_candidates_models_refusals_and_costs_replayed():
    result=replay();models={r['model']:r for r in result['preserving_models']}
    assert len(result['preserving_candidates'])==8 and len(models)==23
    assert len(result['preserving_points'])==184
    assert result['preserving_native_SCFs']==24
    assert result['preserving_allocation_core_hours']==pytest.approx(27.946666666666665)
    assert models['full_space_harmonic_reference@0.02']['supported_predictions']==7
    assert models['anharmonic_retained_residual@0.01']['raw_max_error_meV_fu']==pytest.approx(.03562514530131011)
    assert result['preserving_all23_common_ids']==[] and result['preserving_original_2D_refusal_preserved']


def test_changed_frozen_P8_prediction_is_refused_not_silently_rescored():
    sources=deepcopy(load_sources(ROOT))
    sources['P8_errors']['arms'][0]['points'][0]['prediction_raw_eV_cell']+=.01
    with pytest.raises(ValueError,match='pre-label panel'):preserving_tables(sources)


def test_changed_P8_native_label_cannot_be_shared_as_a_new_calibration():
    sources=deepcopy(load_sources(ROOT))
    sources['P8_native']['rows'][0]['label']['energy_change_eV_cell']+=.01
    with pytest.raises(ValueError,match='same native label'):preserving_tables(sources)
