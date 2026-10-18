"""Actual seven-call +0.5% terminal, not a synthetic target or forecast label."""
import importlib.util
from pathlib import Path
import json
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019'


def load(name):
    spec = importlib.util.spec_from_file_location('E275_completed_'+name, CASE/(name+'.py'))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def test_real_midpoint_terminal_all_native_calls_and_zero_updates(tmp_path):
    b, cost = load('analyze_terminal').read_terminal('strain_0050', tmp_path)
    assert cost['job_id'] == '28951226' and cost['calls'] == 7 and cost['step'] == 0
    assert cost['ordinary_converged'] and cost['all_new_native_EFS_reparsed_locally']
    assert not cost['full_orbital_bytes_rechecked_locally']
    assert cost['allocated_core_hours'] == pytest.approx(6.1777777778)
    assert b['sampled_polygon_peak_meV_fu'] == pytest.approx(87.71854473298)
    assert b['fmax_eV_A'] == pytest.approx(.09231530152813)
    assert b['rows'][4]['polar_proxy_y_A'] == pytest.approx(.12108229235712)


def test_saved_comparison_keeps_all_four_midpoint_initializations_and_pending_high_strain():
    data = json.loads((CASE/'completed_analysis_0050/analysis.json').read_bytes())
    assert data['status'] == 'partial_terminal_connection_comparison'
    assert not data['complete_selected_target_pair_available']
    assert data['pending_terminal_archives'] == ['strain_0100']
    assert data['old_labels_and_scores_unchanged'] and not data['independent_prediction_test']
    assert len([b for b in data['bands'] if b['strain'] == .005]) == 4
    row = next(r for r in data['outcomes'] if r['strain'] == .005)
    assert row['old_T_minus_preserving_meV_fu'] == pytest.approx(-8.1126985806)
    assert row['new_T_minus_preserving_meV_fu'] == pytest.approx(2.6996336028)
    assert data['physical_error_bound'] is None


def test_real_figure_rejects_missing_target_and_keeps_endpoint_forces_undefined(tmp_path):
    p = load('plot_comparison'); source = CASE/'completed_analysis_0050/analysis.json'
    with pytest.raises(ValueError, match='actual completed'): p.draw(source, .01, tmp_path/'missing')
    assert not (tmp_path/'missing').exists()
    qa = json.loads((CASE/'figures/strain_0050/qa.json').read_bytes())
    assert qa['source_rows'] == 36 and qa['all_axes_and_legend_bounds_inside']
    assert not qa['smoothing'] and not qa['inferred_energy_error_bars']
    data = json.loads(source.read_bytes())
    for b in data['bands']:
        assert b['rows'][0]['NEB_residual_eV_A'] is None
        assert b['rows'][-1]['NEB_residual_eV_A'] is None
