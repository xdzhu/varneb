"""Actual scalar diagnostics and basis-invariant algebra; not forecasts."""
import json
from pathlib import Path

import numpy as np
import pytest

from scripts import analyze_hfo2_G3_reference_complements as audit
from scripts import verify_hfo2_reference_transport as reference

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/G3_reference_complements_E156_20261014'


def fixture_arrays():
    """Explicit arithmetic fixture, not a material measurement."""
    rng = np.random.default_rng(156)
    b = np.zeros((40, 36)); b[:33, :33] = np.eye(33); b[36:39, 33:] = np.eye(3)
    h = np.diag(np.r_[-2., -.4, np.linspace(.5, 7., 34)])
    h[0, -1] = h[-1, 0] = .3
    h[2, -2] = h[-2, 2] = .2
    q = b[:, :4].copy()
    return h, b, q, rng.normal(size=40)


@pytest.mark.parametrize('channel', ['T', 'preserving'])
def test_actual_all_seven_references_keep_raw_signs_rank_and_strong_optical_control(channel):
    result = json.loads((CASE/'scalar_report.json').read_text())
    rows = [r for r in result['rows'] if r['channel'] == channel]
    expected = [('T_geometric_triplet', 3), ('T_two_complete_optical_doublets', 4), ('T_all_optical', 33),
                *((f'Cmma_frame_{i}', 4) for i in (103, 106, 114, 127))]
    assert [(r['reference'], r['retained_rank']) for r in rows] == expected
    assert all(r['step_A'] == .01 and r['full_internal_raw_negative_count'] == 2 for r in rows)
    for row in rows:
        assert row['signs_are_raw_not_resolution_certified']
        assert not row['released_branch_or_stationary_index_or_TS_certified']
        assert row['atomic_release']['released_dimension'] == 33-row['retained_rank']
        assert row['joint_release']['released_dimension'] == 36-row['retained_rank']
        if row['reference'] == 'T_all_optical':
            atomic, joint = row['atomic_release'], row['joint_release']
            assert atomic['released_dimension'] == 0 and atomic['raw_minimum_eigenvalue_eV_A2'] is None
            assert atomic['empty_release_is_not_a_physical_stability_test']
            assert joint['raw_minimum_eigenvalue_eV_A2'] > 2.9
            assert joint['raw_negative_count'] == 0
            assert joint['raw_positive_block_Schur_retained_negative_count'] == 2
        else:
            for name in ('atomic_release', 'joint_release'):
                assert row[name]['raw_minimum_eigenvalue_eV_A2'] < 0
                assert row[name]['raw_negative_count'] >= 1
                assert row[name]['raw_positive_block_Schur_retained_negative_count'] is None
    assert result['raw_or_transformed_author_arrays_exported'] is False
    assert result['conditional_axes_or_forecasts_selected'] is False
    assert result['source_frequencies_used_as_PBE_curvature'] is False
    assert result['new_DFT_or_scheduler_operations'] == 0
    assert result['new_raw_SCF_reaudit'] is False


def test_orthogonal_internal_chart_and_retained_basis_changes_preserve_restrictions():
    h, b, q, g = fixture_arrays()
    original = audit.release_diagnostics(h, b, q, g)
    rng = np.random.default_rng(157)
    rotation = np.linalg.qr(rng.normal(size=(36, 36)))[0]
    qr = np.linalg.qr(rng.normal(size=(4, 4)))[0]
    changed = audit.release_diagnostics(rotation.T @ h @ rotation, b @ rotation, q @ qr, g)
    for key in ('full_internal_raw_negative_count', 'retained_projected_raw_negative_count'):
        assert changed[key] == original[key]
    assert changed['mean_raw_negative_subspace_squared_coverage'] == pytest.approx(
        original['mean_raw_negative_subspace_squared_coverage'], abs=1e-14)
    for name in ('atomic_release', 'joint_release'):
        for key in ('raw_minimum_eigenvalue_eV_A2', 'released_gradient_norm_eV_A'):
            assert changed[name][key] == pytest.approx(original[name][key], abs=2e-14)
        assert changed[name]['raw_positive_block_Schur_retained_negative_count'] == 2


def test_negative_release_is_exposed_not_pseudoinverted_or_confused_with_external_control():
    h, b, _, g = fixture_arrays()
    q = b[:, 1:5]  # omit the first true negative direction
    baseline = audit.release_diagnostics(h, b, q, g)
    g[-1] = -1e7  # prescribed reaction gradient cannot change internal restrictions
    changed = audit.release_diagnostics(h, b, q, g)
    assert baseline == changed
    for name in ('atomic_release', 'joint_release'):
        assert baseline[name]['raw_negative_count'] == 1
        assert baseline[name]['raw_positive_block_Schur_retained_negative_count'] is None


@pytest.mark.parametrize('fault', ['complex', 'missing_measured', 'mixed_cell_reference', 'mixed_control_basis',
                                 'rank_loss', 'wrong_atomic_space', 'asymmetric', 'nonfinite'])
def test_bad_or_incomplete_spaces_are_refused(fault):
    h, b, q, g = fixture_arrays()
    if fault == 'complex': q = q.astype(complex)
    elif fault == 'missing_measured': q[35, 0] = 1.; q[0, 0] = 0.
    elif fault == 'mixed_cell_reference': q[:, 0] = b[:, -1]
    elif fault == 'mixed_control_basis': b[-1, 0] = .1
    elif fault == 'rank_loss': q[:, 0] = q[:, 1]
    elif fault == 'wrong_atomic_space': b[33, -1] = 1.; b[38, -1] = 0.
    elif fault == 'asymmetric': h[0, 1] = .1
    elif fault == 'nonfinite': h[0, 0] = np.nan
    with pytest.raises(ValueError): audit.release_diagnostics(h, b, q, g)


@pytest.mark.parametrize('existing', [True, False])
def test_author_source_and_existing_outputs_are_protected(tmp_path, monkeypatch, existing):
    source = tmp_path/'private'; source.mkdir(); output = source/'report.json'
    if existing: output.write_text('untouched')
    monkeypatch.setattr('sys.argv', ['analyze', '--source-root', str(source), '--T-report', 'unused',
        '--T-sha256', 'unused', '--preserving-report', 'unused', '--preserving-sha256', 'unused',
        '--output', str(output)])
    with pytest.raises((FileExistsError, ValueError)): audit.main()
    if existing: assert output.read_text() == 'untouched'
    else: assert not output.exists()


def test_original_reference_loader_keeps_same_source_gate_before_private_arrays(tmp_path):
    with pytest.raises(FileNotFoundError): reference.load_registered_reference_bank(tmp_path)
