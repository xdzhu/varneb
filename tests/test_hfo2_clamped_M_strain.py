"""Actual matched M strains; no external author arrays or live DFT in tests."""
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from scripts import plot_hfo2_clamped_M_strain as plot
from scripts import plot_hfo2_clamped_comparison as paired
from scripts.plot_hfo2_clamped_T_strain import build_strain_data
from scripts.analyze_hfo2_clamped_external_work import ADMITTED_TRAINING_AUDITS, analyze

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/plot.REPORT.parent


def test_actual_eight_completed_training_receipts_only_and_M_raw_stress_replay():
    assert len(ADMITTED_TRAINING_AUDITS) == 8 and plot.AUDIT_SHA in ADMITTED_TRAINING_AUDITS
    r = analyze(ROOT/plot.TERMINAL, plot.AUDIT_SHA)
    saved = json.loads((CASE/'external_work.json').read_text())
    # Preserve the old production-independent LF/CRLF checkout evidence.
    # Only this exact source-metadata variant differs in a clean archive;
    # native inputs, output bytes and every scientific field stay strict.
    module = 'vcneb/biaxial_curvature.py'
    source = (ROOT/module).read_bytes()
    lf = 'd63eb858a07e58b176f21c5f8203b6a116b0a9c7134e204818e0441eb92bc0c7'
    crlf = 'd3d2d32dd2aa4736ae9ec6b0cf52d56c2d00f56c4af7e2b296fc06e58ee672ce'
    assert hashlib.sha256(source.replace(b'\r\n', b'\n')).hexdigest() == lf
    assert saved['core_source_sha256'][module] == lf
    assert r['core_source_sha256'][module] == hashlib.sha256(source).hexdigest()
    assert r['core_source_sha256'][module] in {lf, crlf}
    r['core_source_sha256'][module] = lf
    # E148 appends one native-audited training receipt, without changing any
    # computation. Keep the old report's actual executed-source provenance;
    # compare every scientific field strictly, not a rewritten old report.
    script = ROOT/'scripts/analyze_hfo2_clamped_external_work.py'
    assert hashlib.sha256(script.read_bytes().replace(b'\r\n', b'\n')).hexdigest() == '6b5f75cca07ff5baa18a03074d8a5ce0c72d42100998aa0ea3d728d566a6d0e3'
    assert r['script_sha256'] == hashlib.sha256(script.read_bytes()).hexdigest()
    assert saved['script_sha256'] == 'ad954a7710fa6b169cff7d4027f972e50fbb85daf6620c4a8791b972fee0a508'
    r['script_sha256'] = saved['script_sha256']
    assert r == saved
    assert r['source_job_id'] == '28727368' and r['strain_training_condition'] == .01
    assert r['native_portable_frames_reparsed'] == 9 and not r['complete_six_physical_bytes_rechecked_locally']
    assert r['sampled_peak_image'] == 3
    assert r['sampled_peak']['fixed_open_row_relative_PO_partial_meV_fu_per_unit_strain'] == pytest.approx(-396.1197570751782)
    assert r['new_DFT_calls'] == r['scheduler_mutations'] == 0
    assert not r['G3_selection_or_holdout_access']


def test_actual_paired_M_training_not_same_initial_cache_or_stationary_response():
    d = plot.build_data(ROOT)
    assert d['training_strains'] == [0., .01] and d['energies_relative_to_each_matched_PO']
    assert not d['different_strains_reuse_same_initial_SCF']
    assert d['observed_sampled_forward_peak_increment_meV_fu'] == pytest.approx(-1.2304034162298194)
    assert d['observed_PO_energy_increment_eV_cell'] == pytest.approx(.027226009508012794)
    assert not d['G3_selection_or_holdout_access'] and not d['independent_prediction']
    assert not d['TS_certified'] and not d['model_accuracy_certified']
    assert [b['channel'] for b in d['bands']] == ['M', 'M']
    assert [b['highest_image_index'] for b in d['bands']] == [3, 3]
    z, s = [b['rows'] for b in d['bands']]
    peak = s[3]
    assert peak['T_geometric_triplet_rank3_squared_fraction'] == pytest.approx(.5492512791711277)
    assert peak['T_geometric_triplet_retained_squared_gradient_fraction'] == pytest.approx(.039008425653634)
    assert peak['T_two_lowest_complete_optical_doublets_retained_squared_gradient_fraction'] == pytest.approx(.3177927318953525)
    assert peak['Cmma_frame_103_retained_squared_gradient_fraction'] == pytest.approx(.9158959031480227)
    assert peak['T_geometric_triplet_internal_gradient_norm'] == pytest.approx(.035561530520180615)
    assert s[-1]['energy_relative_PO_meV_fu'] - z[-1]['energy_relative_PO_meV_fu'] == pytest.approx(-22.206478717287314)
    for band in d['bands']:
        for row in band['rows']:
            assert -135. < row['energy_relative_PO_meV_fu'] < 125.
            for span in paired.SPANS:
                assert row[f'{span}_internal_gradient_norm']**2 == pytest.approx(
                    row[f'{span}_retained_gradient_norm']**2 + row[f'{span}_omitted_gradient_norm']**2, abs=1e-12)


@pytest.mark.parametrize('strains,channel', [([0., .005], 'M'), ([0., .01], 'other'), ([.01, 0.], 'M')])
def test_generic_pair_builder_rejects_holdout_or_unregistered_pair_before_labels(strains, channel):
    with pytest.raises(ValueError, match='before label access'):
        build_strain_data(None, [(s,) for s in strains], channel=channel)


@pytest.mark.parametrize('relative', [plot.REPORT, plot.GRADIENT, plot.OBSERVATION, plot.AUDIT])
def test_changed_plus1_M_material_rejects_before_export(relative, tmp_path):
    from scripts.plot_hfo2_clamped_M import MODE_REPORT, OBSERVATION, AUDIT
    for path in (MODE_REPORT, paired.M_GRADIENT, OBSERVATION, AUDIT,
                 plot.REPORT, plot.GRADIENT, plot.OBSERVATION, plot.AUDIT):
        dest = tmp_path/path
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/path, dest)
    value = json.loads((tmp_path/relative).read_text())
    value['tampered_training_input'] = True
    (tmp_path/relative).write_text(json.dumps(value))
    with pytest.raises(ValueError, match='frozen material|pinned completed'):
        plot.build_data(tmp_path)


def test_actual_M_signed_work_keeps_omitted_atoms_cell_and_all_strong_frames():
    r = json.loads((CASE/'path_work.json').read_text())
    assert r['source_job_id'] == '28727368' and r['sampled_peak_image'] == 3
    assert r['terminal_audit_sha256'] == plot.AUDIT_SHA and r['formula_units'] == 4
    assert list(r['representations']) == paired.SPANS
    for name, s in r['representations'].items():
        assert s['retained_rank'] == (3 if name == paired.SPANS[0] else 4)
        assert not s['physical_error_bound_measured']
        assert not s['source_modes_or_forces_ASR_corrected']
        assert s['independent_native_work_max_defect_eV_cell'] < 1e-12
        np.testing.assert_allclose(np.sum(list(s['segment_partition_work'].values()), axis=0),
                                   s['segment_full_work'], atol=1e-12, rtol=0)
        assert s['cumulative_full_work'][3]*250 == pytest.approx(75.9225790735665)
        assert s['cumulative_energy_minus_work'][3]*250 == pytest.approx(21.30513946383038)
        assert s['cumulative_energy_minus_work'][-1]*250 == pytest.approx(2.658882212207794)


def test_manuscript_M_rows_replay_actual_data_without_future_labels():
    text = (ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    d = plot.build_data(ROOT)
    band = d['bands'][1]
    peak = band['rows'][band['highest_image_index']]
    forward = peak['energy_relative_PO_meV_fu']
    well = band['rows'][-1]['energy_relative_PO_meV_fu']
    assert f'| PO-to-M, +1% substrate strain | {forward:.6f} | {forward-well:.6f} | {well:.6f} |' in text
    assert '| PO-to-M, +1% | T geometric triplet | 3 | 54.9251 | 3.9008 |' in text
    work = json.loads((CASE/'path_work.json').read_text())['representations'][paired.SPANS[0]]
    values = [forward, work['cumulative_full_work'][3]*250,
              work['cumulative_energy_minus_work'][3]*250, work['cumulative_energy_minus_work'][-1]*250]
    assert '| PO-to-M, +1% strain | ' + ' | '.join(f'{v:.6f}' for v in values) + ' |' in text
    section = text.split('### 3.1.', 1)[1].split('### 3.2.', 1)[0]
    assert 'not independent predictive wins' in section
    assert 'ordinary residuals do not bound sampled-barrier error' in section
    assert 'same four registered correspondences' in section
    work_section = text.split('### 3.4.', 1)[1].split('### 3.5.', 1)[0]
    assert 'not a measured barrier uncertainty' in work_section
    assert 'Eight of eight registered G2 training bands' in text


def test_same_grid_style_no_endpoint_clipping_and_legends_do_not_overlap():
    import matplotlib.pyplot as plt
    fig = paired.make_figure(plot.build_data(ROOT))
    try:
        positions = np.array([ax.get_position().bounds for ax in fig.axes]).reshape(3, 2, 4)
        np.testing.assert_allclose(positions[1:, :, [0, 2]], np.repeat(positions[0:1, :, [0, 2]], 2, axis=0))
        np.testing.assert_allclose(positions[:, 1, [1, 3]], positions[:, 0, [1, 3]])
        renderer = fig.canvas.get_renderer()
        legends = fig.legends + [a.get_legend() for a in fig.axes if a.get_legend() is not None]
        for i, ax in enumerate(fig.axes):
            assert not ax.get_title() and all(s.get_visible() for s in ax.spines.values())
            assert ax.texts[0].get_text() == f'({chr(97+i)})' and ax.texts[0].get_weight() == 'normal'
            for legend in legends:
                assert not legend.get_window_extent(renderer).overlaps(ax.get_window_extent(renderer))
                assert not legend.get_window_extent(renderer).overlaps(ax.texts[0].get_window_extent(renderer))
                assert legend.get_frame().get_alpha() == .85
        assert fig.axes[0].get_ylim() == fig.axes[1].get_ylim() == (-135., 125.)
        assert len(fig.axes[3].lines) == 7  # every tied Cmma frame is retained
    finally:
        plt.close(fig)
