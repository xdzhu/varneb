"""Completed original +1% handle and complete three-strain native mechanism data."""
import hashlib
import importlib.util
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019'


def load(name):
    spec = importlib.util.spec_from_file_location('E275_full_'+name, CASE/(name+'.py'))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def test_actual_high_strain_terminal_replays_every_new_EFS_call(tmp_path):
    band, cost = load('analyze_terminal').read_terminal('strain_0100', tmp_path)
    assert cost['job_id'] == '28951227' and cost['calls'] == 35 and cost['step'] == 4
    assert cost['ordinary_converged'] and cost['all_new_native_EFS_reparsed_locally']
    assert cost['allocated_core_hours'] == pytest.approx(29.8311111111)
    assert band['sampled_polygon_peak_meV_fu'] == pytest.approx(91.6989365251)
    assert band['fmax_eV_A'] == pytest.approx(.09093136383856)
    assert band['rows'][4]['polar_proxy_y_A'] == pytest.approx(.12267478821244)


def test_all_actual_strains_preserve_original_labels_and_report_bounded_candidate_response():
    data = json.loads((CASE/'completed_analysis_full/analysis.json').read_bytes())
    assert data['complete_selected_target_pair_available'] and data['pending_terminal_archives'] == []
    assert len(data['bands']) == 10 and sum(c['calls'] for c in data['terminal_cohorts']) == 42
    assert data['old_labels_and_scores_unchanged'] and not data['independent_prediction_test']
    assert data['physical_error_bound'] is None and not data['global_MEP_or_TS_certified']
    assert [r['released_seed_meV_fu'] for r in data['outcomes']] == pytest.approx([
        83.1662809401, 87.7185447330, 91.6989365251])
    assert [r['new_T_minus_preserving_meV_fu'] for r in data['outcomes']] == pytest.approx([
        1.2503429411, 2.6996336028, 5.0647138696])
    old = data['outcomes'][-1]['original_preserving_meV_fu']-data['outcomes'][0]['original_preserving_meV_fu']
    new = data['outcomes'][-1]['released_seed_meV_fu']-data['outcomes'][0]['released_seed_meV_fu']
    assert old == pytest.approx(28.5962481553) and new == pytest.approx(8.5326555850)


def test_main_summary_contains_nine_real_rows_and_retains_all_symmetry_tolerances():
    source = CASE/'completed_analysis_full/analysis.json';data = json.loads(source.read_bytes())
    rows = load('plot_summary').summary_rows(data)
    assert len(rows) == 9 and all(r['ordinary_residual_pass'] for r in rows)
    qa = json.loads((CASE/'figures/summary_three_strains/qa.json').read_bytes())
    assert qa['source_analysis_sha256'] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert qa['all_axes_legend_bounds_inside'] and not qa['interpolation_or_fit']
    preserving = [b for b in data['bands'] if b['kind'] in ('registered_preserving', 'released_seed')]
    assert len(preserving) == 6
    for b in preserving:
        audit = b['peak_structure_audit']
        assert [r['symprec_A'] for r in audit['symmetry_sweep']] == [.001, .01, .05]
        assert not audit['standardized_wrapped_reordered_or_relaxed']
        assert b['peak_structure_not_a_stable_phase_or_electronic_polarization']


def test_manuscript_uses_complete_native_outcomes_without_replacing_the_failed_local_control():
    text = (ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    main = text.split('## 3.', 1)[1].split('## 4.', 1)[0]
    data = json.loads((CASE/'completed_analysis_full/analysis.json').read_bytes())
    for row in data['outcomes']:
        assert f"{row['released_seed_meV_fu']:.6f}" in main
        assert f"{row['new_T_minus_preserving_meV_fu']:.6f}" in main
    assert 'Four updates and 35 fresh interior SCFs' in main
    assert 'so no optimizer\nupdate was performed' in main
    assert 'Figure 9. Completed same-strain comparison' in main
    assert 'summary_three_strains/strain_mechanism_summary.pdf' in main
    assert 'their J energies remain unavailable' in main
    assert 'It is not a substitute for the missing local A/J' in text
    assert '8.532656 rather than 28.596248' in text
