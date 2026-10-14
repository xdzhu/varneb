"""Saved scalar material evidence; no private reference arrays or new DFT in tests."""
from pathlib import Path
import importlib.util
import json
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'benchmarks/hfo2_channels/20261008/G3_shared_Q_controls_E180_20261015'
s = importlib.util.spec_from_file_location('shared_Q_controls_tests', CASE / 'analyze.py')
a = importlib.util.module_from_spec(s)
s.loader.exec_module(a)


@pytest.fixture(scope='module')
def report():
    return json.loads((CASE / 'analysis.json').read_bytes())


def test_all_references_frames_steps_points_and_holes_are_retained(report):
    assert len(report['rows']) == 32 and report['model_point_rows'] == 800
    for channel in ('T', 'preserving'):
        for step in (.01, .02):
            rows = [r for r in report['rows'] if r['channel'] == channel and r['step_A'] == step]
            assert [r['reference'] for r in rows] == [n for n, _ in a.NAMES] + ['adaptive_two_directions_plus_complete_complement']
            for r in rows:
                assert [p['id'] for p in r['points']] == list(range(25))
                assert r['shared_observable_count'] == 2
                assert r['unmodeled_dimensions_clamped'] == 36 - r['allowed_dimension']
                assert r['free_dimension'] == r['allowed_dimension'] - 2


def test_complete_fixed_reference_equals_schur_not_new_physical_accuracy(report):
    for r in report['complete_space_equivalence']:
        assert r['max_complete_T_vs_adaptive_Schur_energy_difference_eV_cell'] < 1e-14
        assert r['max_complete_T_vs_adaptive_Schur_coordinate_difference_A'] < 1e-12
        assert r['equality_is_standard_same_problem_algebra_not_an_accuracy_improvement']
    for channel in ('T', 'preserving'):
        for step in (.01, .02):
            selected = [r for r in report['rows'] if r['channel'] == channel and r['step_A'] == step]
            full = next(r for r in selected if r['reference'] == 'T_all_optical')
            adaptive = selected[-1]
            assert full['nominal_atomic_reference_rank'] == 33
            assert full['allowed_dimension'] == adaptive['allowed_dimension'] == 36
            np.testing.assert_allclose([p['model_energy_change_eV_cell'] for p in full['points']],
                [p['model_energy_change_eV_cell'] for p in adaptive['points']], atol=1e-13, rtol=0)


def test_no_common_nonorigin_cohort_means_no_all_reference_accuracy_ranking(report):
    for r in report['complete_space_equivalence']:
        assert r['common_all_reference_non_origin_ids'] == [] and r['no_non_origin_common_comparison']
    for r in report['rows']:
        stat = r['common_all_reference_geometry_eligible_statistics']
        assert stat['count'] == 0 and stat['ids'] == []
        assert stat['max_absolute_energy_discrepancy_meV_cell'] is None
        assert stat['max_absolute_seen_origin_aligned_shape_discrepancy_meV_cell'] is None


def test_small_span_outside_domain_values_are_diagnostics_not_bad_material_predictions(report):
    for r in report['rows']:
        for p in r['points']:
            assert p['outside_guard_values_diagnostic_only_not_material_predictions'] is not p['passes_original_geometry_guard']
        eligible = r['own_geometry_eligible_statistics']['ids']
        assert all(r['points'][i]['passes_original_geometry_guard'] for i in eligible)
    triplet = next(r for r in report['rows'] if r['channel'] == 'preserving'
                   and r['step_A'] == .01 and r['reference'] == 'T_geometric_triplet')
    assert triplet['own_geometry_eligible_statistics']['count'] == 0
    assert triplet['own_geometry_eligible_statistics']['max_absolute_energy_discrepancy_meV_cell'] is None


@pytest.mark.parametrize('channel,expected', [('T', 25), ('preserving', 17)])
def test_missing_DFT_candidates_not_replaced_by_last_trials(report, channel, expected):
    for r in [r for r in report['rows'] if r['channel'] == channel]:
        assert len(r['observed_DFT_candidate_ids']) == expected
        assert sum(p['DFT_energy_change_eV_cell'] is not None for p in r['points']) == expected
        for p in r['points']:
            if p['DFT_energy_change_eV_cell'] is None:
                assert p['model_minus_DFT_energy_eV_cell'] is None
                assert p['seen_origin_aligned_shape_discrepancy_eV_cell'] is None


def test_no_new_data_or_private_vectors_or_independent_success_exported(report):
    assert report['fitted_coefficients'] == report['new_DFT_or_scheduler_writes'] == report['independent_material_slots_used'] == 0
    assert report['restricted_controls_are_not_full_outside_reference_release_B4']
    assert not report['raw_transformed_author_basis_or_response_arrays_exported']
    assert not report['physical_reserved_geometry_generated_or_read']
    assert not report['stable_conditional_branch_TS_physical_error_or_independent_advantage_certified']
    forbidden = {'allowed_basis', 'response', 'source_basis', 'local_basis', 'full_displacement', 'coordinates_full_A'}
    def visit(value):
        if isinstance(value, dict):
            assert not forbidden.intersection(value)
            for child in value.values(): visit(child)
        elif isinstance(value, list):
            for child in value: visit(child)
    visit(report)
    assert report['source_registration_sha256'] == a.REGISTRATION_SHA256


def test_changed_face_pin_refuses_before_any_observation_arithmetic(monkeypatch):
    monkeypatch.setattr(a, 'FACE_PINS', {**a.FACE_PINS, 'T': (a.FACE_PINS['T'][0], '0'*64)})
    with pytest.raises(ValueError, match='native face bytes'): a.observed_points('T', {})


def test_statistics_refuse_missing_declared_DFT_and_do_not_insert_zero():
    with pytest.raises(ValueError, match='all declared DFT'):
        a.statistics([dict(id=1, DFT_energy_change_eV_cell=None)], {1})
    empty = a.statistics([], set())
    assert empty['count'] == 0 and empty['max_absolute_energy_discrepancy_meV_cell'] is None


def test_existing_output_refused_before_private_source_access(tmp_path, monkeypatch):
    out = tmp_path / 'old.json'
    out.write_text('preserve')
    monkeypatch.setattr(a.sys, 'argv', ['compare', '--source-root', 'unused', '--output', str(out)])
    monkeypatch.setattr(a, 'analyze', lambda *args: pytest.fail('must refuse before reading private inputs'))
    with pytest.raises(FileExistsError): a.main()
    assert out.read_text() == 'preserve'
