"""A mixed quadratic-block observation is not a new physical stability result."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

CASE = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/matched_pattern_release_E270_20261018'
spec = importlib.util.spec_from_file_location('E270_diagnostic', CASE/'diagnose.py')
diagnostic = importlib.util.module_from_spec(spec); spec.loader.exec_module(diagnostic)


def test_synthetic_positive_blocks_can_have_negative_mixed_direction():
    matrix = np.diag([2., 3., 1., 1., 1.]); matrix[0, 2] = matrix[2, 0] = 2.
    result = diagnostic.block_screen(matrix, 2)
    assert result['separate_blocks_positive'] and result['removing_mixed_block_removes_raw_negative']
    assert result['coupled_negative_count'] == 1
    assert result['cell_Schur_eigenvalues_eV_A2'] == pytest.approx([-1., 1., 1.])
    assert not result['stationary_physical_stability_or_error_bound_certified']
    matrix[0, 2] = matrix[2, 0] = 0.
    assert diagnostic.block_screen(matrix, 2)['coupled_negative_count'] == 0


def test_nonpositive_atomic_block_refuses_conditional_elimination():
    matrix = np.diag([-1., 3., 1., 1., 1.])
    result = diagnostic.block_screen(matrix, 2)
    assert result['cell_Schur_eigenvalues_eV_A2'] is None
    assert result['atomic_condition_number'] is None and not result['separate_blocks_positive']


@pytest.mark.parametrize('fault', ['shape', 'nan', 'asymmetric'])
def test_invalid_raw_matrix_refused(fault):
    matrix = np.eye(5)
    if fault == 'shape': matrix = matrix[:4]
    elif fault == 'nan': matrix[0, 0] = np.nan
    else: matrix[0, 1] = 1.
    with pytest.raises(ValueError): diagnostic.block_screen(matrix, 2)


def test_actual_two_steps_and_first_native_trial_keep_their_limits():
    result = diagnostic.diagnose(CASE)
    assert [r['probe_step_A'] for r in result['rows']] == [.01, .02]
    assert all(r['separate_blocks_positive'] and r['coupled_negative_count'] == 1 for r in result['rows'])
    assert [r['coupled_eigenvalue_range_eV_A2'][0] for r in result['rows']] == pytest.approx(
        [-.00367328239635, -.00621625978468])
    assert all(r['atomic_condition_number'] > 6000 for r in result['rows'])
    for row in result['rows']:
        assert sum(row[f'weak_direction_{k}_squared_norm'] for k in ('Hf','oxygen','scaled_cell')) == pytest.approx(1.)
    trial = result['first_trial']
    assert trial['step_norm_A'] == pytest.approx(.05)
    assert trial['kappa_energy_eV_A2'] == pytest.approx(29.511371404586573)
    assert trial['kappa_force_eV_A2'] == pytest.approx(29.417598178620924)
    assert all(28 < v < 31 for v in trial['original_measured_Rayleigh_eV_A2'])
    assert result['new_DFT_calls'] == result['new_curvature_probes'] == 0
    assert not result['conditional_stability_or_TS_or_prediction_advantage_certified']
    assert not result['live_source_or_parameters_changed']
    assert not result['numerical_fallback_execution_admitted']
    assert result['native_STRU_identity_checked']
    assert np.any(result['constant_reference_integer_offset'])


def test_periodic_reference_offset_is_constant_for_atomic_and_cell_increments():
    from scripts.hfo2_G3_conditional_faces import load_inputs
    spec = importlib.util.spec_from_file_location('E270_reference_gauge', CASE/'study.py')
    study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)
    original = load_inputs('T')['chart']
    native = study.chart_for(CASE/'prepared/seeds/0')
    offset = native.reference_atoms.get_scaled_positions(wrap=False)-original.reference_atoms.get_scaled_positions(wrap=False)
    rng = np.random.default_rng(927)
    for x in rng.normal(scale=.002, size=(5, 39)):
        first, second = native.to_atoms(x), original.to_atoms(np.r_[x, 0.])
        np.testing.assert_allclose(first.cell.array, second.cell.array, atol=1e-12, rtol=0)
        np.testing.assert_allclose(first.get_scaled_positions(wrap=False)-second.get_scaled_positions(wrap=False),
                                  offset, atol=1e-12, rtol=0)
