"""Verify frozen scalar controls without private mode arrays or test labels."""
from pathlib import Path
import hashlib
import importlib.util
import json

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'benchmarks/hfo2_channels/20261008/G3_reference_forecast_E188_20261015'
spec = importlib.util.spec_from_file_location('E188_reference_forecast', CASE/'freeze.py')
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)


@pytest.fixture(scope='module')
def report():
    return json.loads((CASE/'forecast.json').read_bytes())


def test_all_seven_references_frames_steps_and_locked_targets_retained(report):
    assert len(report['rows']) == 28 and report['scalar_model_points'] == 224
    locked = f.locked_forecast()
    for channel in ('T','preserving'):
        targets = next(r['targets'] for r in locked['channels'] if r['channel'] == channel)
        for step in (.01,.02):
            rows = [r for r in report['rows'] if r['channel'] == channel and r['step_A'] == step]
            assert [(r['reference'],r['nominal_atomic_rank']) for r in rows] == f.control.NAMES
            for row in rows:
                assert row['allowed_dimension'] == row['nominal_atomic_rank'] + 3
                assert row['free_dimension'] == row['allowed_dimension'] - 2
                assert row['omitted_directions_clamped'] == 36 - row['allowed_dimension']
                assert [p['q_A'] for p in row['points']] == [p['q_A'] for p in targets]


def test_coverage_and_empty_common_set_do_not_imply_bad_accuracy(report):
    assert f.summarize_coverage(report['rows']) == report['coverage']
    for row in report['coverage']:
        counts = row['supported_per_reference']
        assert counts['T_all_optical'] == counts['T_two_complete_optical_doublets'] == 8
        if row['channel'] == 'T':
            assert row['common_all_reference_ids'] == [1,2,3,4]
            assert counts['T_geometric_triplet'] == 6
            assert all(counts[f'Cmma_frame_{i}'] == 6 for i in (103,106,114,127))
        else:
            assert row['common_all_reference_ids'] == [] and counts['T_geometric_triplet'] == 0
            assert all(counts[f'Cmma_frame_{i}'] == 2 for i in (103,106,114,127))
        assert row['denominator'] == 8 and row['independent_labels'] == 0
    assert sum(p['prediction_calibrated_eV_cell'] is not None for r in report['rows'] for p in r['points']) == 140


def test_outside_guard_values_are_diagnostics_and_never_filled(report):
    for row in report['rows']:
        for point in row['points']:
            assert point['independent_DFT_label'] is None
            assert np.isfinite(point['formal_energy_change_eV_cell'])
            valid = point['passes_original_geometry_guard'] and row['origin_geometry_eligible']
            assert (point['prediction_calibrated_eV_cell'] is not None) is valid
            if not valid:
                assert point['reason'] is not None
            else:
                assert point['prediction_calibrated_eV_cell'] == pytest.approx(
                    point['formal_energy_change_eV_cell'] - row['origin_bias_eV_cell'], abs=1e-13)


def test_full_reference_reproduces_already_locked_complete_space_predictions(report):
    locked = f.locked_forecast()
    for row in [r for r in report['rows'] if r['reference'] == 'T_all_optical']:
        old = next(r for r in locked['channels'] if r['channel'] == row['channel'])
        key = f"original_H_{row['step_A']:.2f}"
        for point, target in zip(row['points'], old['targets']):
            assert point['prediction_raw_eV_cell'] == pytest.approx(target['models'][key]['value'], abs=1e-13)
            assert point['prediction_calibrated_eV_cell'] == pytest.approx(
                target['models'][key]['seen_origin_calibrated_energy_change_eV_cell'], abs=1e-13)
    assert max(r['max_energy_difference_from_E185_complete_H_eV_cell'] for r in report['complete_space_equalities']) < 1e-14


def test_all224_formal_energies_and_norm_squared_match_independent_old_scalar_polynomials(report):
    # x(Q) is affine, so E(Q) and ||x(Q)||^2 are quadratic. Reconstruct
    # independently from nine OLD E180 nodes, without private basis access.
    path = f.BASE/'G3_shared_Q_controls_E180_20261015/analysis.json'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == '39e31d790592e7a2f640a8428f00efecc28fd14775b7fb4b2b699b7da54c6402'
    source = json.loads(path.read_bytes())
    def features(q):
        x,y = np.asarray(q)/.1
        return [1.,x,y,x*x,x*y,y*y]
    for row in report['rows']:
        old = next(r for r in source['rows'] if all(r[k] == row[k] for k in ('channel','reference','step_A')))
        central = [p for p in old['points'] if max(map(abs,p['q_A'])) <= .1]
        design = np.array([features(p['q_A']) for p in central])
        beta,_,rank,_ = np.linalg.lstsq(design,np.array([
            [p['model_energy_change_eV_cell'],p['model_total_internal_norm_A']**2] for p in central]),rcond=None)
        assert rank == 6 and len(central) == 9
        for point in row['points']:
            energy, norm2 = np.asarray(features(point['q_A'])) @ beta
            assert point['formal_energy_change_eV_cell'] == pytest.approx(energy, abs=1e-9)
            assert point['model_internal_norm_A']**2 == pytest.approx(norm2, abs=1e-8)


def test_scoring_budget_information_and_novelty_boundaries_are_preserved(report):
    assert report['inherited_scoring_policy'] == f.locked_forecast()['scoring']
    assert report['independent_slots_reserved_computed_each'] == [8,0]
    assert report['training_commitments_unchanged'] == {'T':34,'preserving':35}
    assert report['restricted_shared_Q_controls_not_complete_native_Cmma_chart_models']
    assert report['same_H_g_Q_and_mechanical_boundary']
    assert not report['full_B0_B5_final_physical_prediction_and_domain_freeze_complete']
    assert not report['physical_plus005_geometry_or_label_access']
    assert report['new_DFT_or_scheduler_writes'] == 0
    assert not report['independent_gain_stability_TS_or_JCTC_readiness_claimed']


def test_private_raw_and_transformed_vectors_are_not_in_report(report):
    forbidden = {'allowed_basis','source_basis','local_basis','coordinates','coordinates_full_A',
                 'response','hessian','raw_mapped_columns','source_to_local'}
    def visit(value):
        if isinstance(value, dict):
            assert not forbidden.intersection(value)
            for child in value.values():visit(child)
        elif isinstance(value,list):
            for child in value:visit(child)
    visit(report)
    assert not report['private_raw_or_transformed_arrays_exported']
    assert not report['private_sources_written']
    assert not report['LDA_source_frequencies_used_as_PBE_curvature']


def test_existing_output_rejected_before_private_input_access(tmp_path,monkeypatch):
    output = tmp_path/'old.json'
    output.write_text('keep')
    monkeypatch.setattr(f.sys,'argv',['freeze','--source-root','unused','--output',str(output)])
    monkeypatch.setattr(f,'freeze',lambda *args:pytest.fail('must refuse first'))
    with pytest.raises(FileExistsError):f.main()
    assert output.read_text() == 'keep'


def test_private_root_cannot_be_output_destination(tmp_path,monkeypatch):
    monkeypatch.setattr(f.sys,'argv',['freeze','--source-root',str(tmp_path),'--output',str(tmp_path/'new.json')])
    monkeypatch.setattr(f,'freeze',lambda *args:pytest.fail('private root must remain read-only'))
    with pytest.raises(ValueError,match='private sources'):f.main()
