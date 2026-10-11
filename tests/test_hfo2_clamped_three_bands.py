"""Pinned real training records and honest layout; no DFT or model labels."""
import csv
import json
from pathlib import Path

import numpy as np
import pytest

from scripts import plot_hfo2_clamped_three_bands as plot

ROOT = Path(__file__).resolve().parents[1]


def test_actual_three_bands_share_physical_zero_without_new_audit_or_forecast():
    data = plot.build_data(ROOT)
    assert [b["channel"] for b in data["bands"]] == ["T", "M", "reversing"]
    assert [b["source_job_id"] for b in data["bands"]] == ["28692775", "28722320", "28692776"]
    assert data["image_records"] == 27
    assert data["distinct_manifest_native_log_STRU_identities"] == 25
    assert data["new_DFT_calls"] == data["scheduler_mutations"] == data["independently_audited_new_SCFs"] == 0
    assert not data["G3_selection_or_holdout_access"] and not data["TS_certified"]
    assert not data["native_data_reaudit_claimed"] and not data["smooth_MEP_interpolation_used"]
    assert data["common_initial_energy_eV_cell"] == pytest.approx(-9782.973830941437)
    expected = [(4, 84.41662388122495), (3, 98.45812195362669), (4, 338.9689653767456)]
    for band, (peak, energy) in zip(data["bands"], expected):
        assert band["highest_image_index"] == peak and band["ordinary_fmax_eV_A"] <= .10
        assert band["rows"][peak]["energy_relative_PO_meV_fu"] == pytest.approx(energy, abs=1e-8)
        assert len(band["rows"]) == 9
        for row in band["rows"]:
            for span in plot.SPANS:
                norm = row[f"{span}_internal_gradient_norm"]
                kept = row[f"{span}_retained_gradient_norm"]
                omitted = row[f"{span}_omitted_gradient_norm"]
                assert norm*norm == pytest.approx(kept*kept+omitted*omitted, abs=1e-12)
    peak = data["bands"][2]["rows"][4]
    assert peak["T_geometric_triplet_rank3_squared_fraction"] == pytest.approx(.262915, abs=1e-6)
    assert peak["T_geometric_triplet_retained_squared_gradient_fraction"] == pytest.approx(.000320, abs=1e-6)
    assert peak["Cmma_frame_103_retained_squared_gradient_fraction"] == pytest.approx(.651046, abs=1e-6)


@pytest.mark.parametrize('pin', ['R_REPORT_SHA', 'R_GRADIENT_SHA', 'R_OBSERVATION_SHA', 'R_AUDIT_SHA'])
def test_changed_frozen_source_rejects_before_drawing(monkeypatch, pin):
    monkeypatch.setattr(plot, pin, '0'*64)
    with pytest.raises(ValueError, match='frozen'):
        plot.build_data(ROOT)


@pytest.mark.parametrize('dependency', list(plot.DEPENDENCY_SHA256))
def test_historical_adapter_code_must_remain_byte_pinned(monkeypatch, dependency):
    changed = plot.DEPENDENCY_SHA256.copy(); changed[dependency] = '0'*64
    monkeypatch.setattr(plot, 'DEPENDENCY_SHA256', changed)
    with pytest.raises(ValueError, match='historical plotting dependency'):
        plot.build_data(ROOT)


def test_shared_initial_identity_is_checked_not_assumed_from_normalized_arcs(monkeypatch):
    original = plot.build_TM
    def wrong_identity(repository):
        data = original(repository)
        data['common_initial_record_canonical_sha256'] = '0'*64
        return data
    monkeypatch.setattr(plot, 'build_TM', wrong_identity)
    with pytest.raises(ValueError, match='common native record'):
        plot.build_data(ROOT)


def test_nine_aligned_panels_legends_labels_spines_and_common_energy_scale():
    import matplotlib.pyplot as plt
    data = plot.build_data(ROOT)
    fig = plot.make_figure(data)
    try:
        assert len(fig.axes) == 9
        np.testing.assert_allclose(fig.get_size_inches()*25.4, [183, 207])
        positions = np.array([a.get_position().bounds for a in fig.axes]).reshape(3, 3, 4)
        np.testing.assert_allclose(positions[1:, :, [0, 2]], np.repeat(positions[0:1, :, [0, 2]], 2, axis=0))
        np.testing.assert_allclose(positions[:, 1:, [1, 3]], np.repeat(positions[:, 0:1, [1, 3]], 2, axis=1))
        renderer = fig.canvas.get_renderer()
        legends = fig.legends+[a.get_legend() for a in fig.axes if a.get_legend() is not None]
        assert len(legends) == 4
        for i, ax in enumerate(fig.axes):
            assert not ax.get_title() and all(s.get_visible() for s in ax.spines.values())
            assert ax.texts[0].get_text() == f'({chr(97+i)})'
            assert ax.texts[0].get_weight() == 'normal' and ax.texts[0].get_fontsize() == 13
            assert ax.texts[0].get_position() == (-.22, 1.025)
            assert ax.xaxis.get_major_ticks()[0].tick2line.get_visible()
            assert ax.yaxis.get_major_ticks()[0].tick2line.get_visible()
            for legend in legends:
                box = legend.get_window_extent(renderer)
                assert not box.overlaps(ax.get_window_extent(renderer))
                assert not box.overlaps(ax.texts[0].get_window_extent(renderer))
                assert legend.get_frame().get_alpha() == .85
        for i, band in enumerate(data['bands']):
            energy_ax = fig.axes[i]
            assert energy_ax.get_ylim() == (-115., 370.)
            np.testing.assert_array_equal(energy_ax.lines[0].get_ydata(),
                [r['energy_relative_PO_meV_fu'] for r in band['rows']])
            np.testing.assert_array_equal(energy_ax.lines[0].get_xdata(),
                [r['normalized_generalized_arc'] for r in band['rows']])
            for offset in (3, 6):
                assert len(fig.axes[offset+i].lines) == 7  # six curves/all four frames + peak guide
        assert np.isnan(fig.axes[3].lines[0].get_ydata()[-1])
        assert np.isnan(fig.axes[3].lines[2].get_ydata()[-1])
        assert np.isfinite(fig.axes[5].lines[0].get_ydata()[-1])
    finally:
        plt.close(fig)


def test_export_preserves_all_scalar_rows_and_editable_text_and_refuses_overwrite(tmp_path):
    from PIL import Image
    output = tmp_path/'new-figure'
    data = plot.build_data(ROOT)
    qa = plot.export(data, output)
    rows = list(csv.DictReader((output/'source_data.csv').open(encoding='utf-8', newline='')))
    assert len(rows) == 27 and [r['channel'] for r in rows[18:]] == ['reversing']*9
    assert qa['image_records'] == 27 and qa['grid_aligned']
    assert qa['all_four_Cmma_frames_retained'] and qa['all_spines_visible'] and qa['panel_titles_empty']
    assert not qa['Cmma_spread_is_statistical_uncertainty']
    assert qa['energy_ylim_meV_fu'] == [-115, 370] and qa['visual_review'] == 'pending'
    svg = (output/(plot.STEM+'.svg')).read_text(encoding='utf-8')
    assert '<text' in svg and '(a)' in svg and '(i)' in svg
    assert {p.suffix for p in output.glob(plot.STEM+'.*')} == {'.pdf', '.svg', '.png', '.tiff'}
    with Image.open(output/(plot.STEM+'.tiff')) as im:
        assert im.info['dpi'] == pytest.approx((600., 600.), abs=.1)
    for name, digest in qa['file_sha256'].items():
        assert plot.sha256(output/name) == digest
    with pytest.raises(FileExistsError, match='fresh dated'):
        plot.export(data, output)


def test_main_text_describes_training_records_and_all_nine_panels_without_TS_claim():
    text = (ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    assert 'hfo2_clamped_TMR_E099_20261011/hfo2_clamped_TMR.pdf' in text
    assert 'The twenty-seven records include cached endpoints' in text
    assert 'Figure 5. Three ordinary-converged zero-strain training candidates' in text
    assert '(a--c)' in text and '(d--f)' in text and '(g--i)' in text
    assert 'sampled peaks, not certified TSs' in text
