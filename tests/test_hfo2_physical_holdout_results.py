"""Actual independent physical-strain results, units, coverage and costs."""
from pathlib import Path

import numpy as np
import pytest

from scripts.analyze_hfo2_physical_holdout import analyze, independent_figure, FORECAST_PIN

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def data():
    return analyze(ROOT)


def test_two_registered_outcomes_and_common_units(data):
    rows = data['holdout_outcomes']
    assert [r['channel'] for r in rows] == ['preserving', 'T']
    assert [r['sampled_peak_meV_fu'] for r in rows] == pytest.approx([98.53087691635665, 90.41817833576715])
    assert [r['response_from_zero_meV_fu'] for r in rows] == pytest.approx([1.7026939135576, 6.0015544545422])
    assert [r['updates'] for r in rows] == [0, 2]
    assert [r['native_SCFs'] for r in rows] == [7, 21]
    assert all(r['force_eV_A'] < .10 for r in rows)
    assert not data['full_four_channel_network_at_half_strain']


def test_frozen_panel_keeps_every_abstention_and_correct_errors(data):
    assert data['forecast_manifest_sha256'] == FORECAST_PIN
    assert (data['numeric_forecasts'], data['unavailable_forecasts']) == (6, 132)
    assert len(data['holdout_scores']) == 69
    models = {r['model']: r for r in data['holdout_scores']}
    assert models['B0']['maximum_absolute_error_meV_fu'] == pytest.approx(12.5954301640832)
    assert models['B1_follow']['maximum_absolute_error_meV_fu'] == pytest.approx(5.0046391475007)
    for row in models.values():
        if row['model'] not in ('B0', 'B1_fixed', 'B1_follow'):
            assert row['scored'] == 0 and row['maximum_absolute_error_meV_fu'] is None
            assert all(p['abstention_reason'] and p['signed_error_meV_fu'] is None for p in row['points'].values())
        else:
            errors = [p['signed_error_meV_fu'] for p in row['points'].values()]
            assert row['T_minus_preserving_error_meV_fu'] == pytest.approx(errors[1]-errors[0])
        assert not row['physically_resolved_advantage_certified']


def test_pair_improvement_is_still_negative_and_not_linear(data):
    gaps = [r['T_minus_preserving_meV_fu'] for r in data['selected_pair_gaps']]
    assert gaps == pytest.approx([-12.4115591215741, -8.1126985805895, -28.6607807634027])
    assert all(g < 0 for g in gaps)
    assert data['holdout_gap_change_from_zero_meV_fu'] == pytest.approx(4.2988605409846)
    assert abs(gaps[1]-(gaps[0]+gaps[2])/2) > 12


def test_original_ordered_mode_frame_and_strain_retained(data):
    assert len(data['bands']) == 10 and len(data['geometry']) == 90
    for band in [b for b in data['bands'] if b['strain'] == .005]:
        rows = band['rows']
        np.testing.assert_allclose([r['F_xx'] for r in rows], 1.005, atol=1e-10)
        np.testing.assert_allclose([r['F_zz'] for r in rows], 1.005, atol=1e-10)
        assert np.all(np.diff([r['s'] for r in rows]) > 0)
        if band['channel'] == 'preserving':
            assert all(r['Q_x_A'] < 0 for r in rows)
            assert rows[0]['polar_dz_A'] > 0 > rows[-1]['polar_dz_A']
        else:
            assert abs(rows[-1]['Q_x_A']) < 1e-6 and rows[-1]['Q_z_A'] > .8


def test_actual_cost_does_not_erase_warm_start_or_preparation(data):
    cost = data['cost']
    assert cost['endpoint_native_SCFs']+cost['path_native_SCFs'] == 39
    assert cost['endpoint_and_path_SCF_core_hours'] == pytest.approx(47.28988831847906)
    assert cost['endpoint_and_path_allocation_core_hours'] == pytest.approx(48.80888888888889)
    assert cost['all_earlier_training_and_model_preparation_additional']
    assert data['no_new_model_advantage']


def test_new_main_figure_retains_all_samples_and_errors(data):
    import matplotlib.pyplot as plt
    fig = independent_figure(data)
    try:
        assert len(fig.axes) == 4
        for ax in fig.axes[:2]:
            assert [len(line.get_xdata()) for line in ax.lines] == [9, 9, 9]
        bars = fig.axes[3].patches
        assert len(bars) == 6
        assert [p.get_height() for p in bars] == pytest.approx([
            12.595430164, -2.933378260, -1.667842943, .171958802, -7.232238800, -5.004639148])
    finally:
        plt.close(fig)
