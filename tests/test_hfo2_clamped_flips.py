"""Real admitted flip observations and layout, without DFT or unseen labels."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from scripts import plot_hfo2_clamped_flips as plot

ROOT = Path(__file__).resolve().parents[1]


def test_actual_flips_and_absolute_gradient_context():
    data = plot.build_data(ROOT)
    assert [b["channel"] for b in data["bands"]] == ["preserving", "reversing"]
    assert [b["source_job_id"] for b in data["bands"]] == ["28727382", "28692776"]
    assert data["image_records"] == 18
    assert data["new_DFT_calls"] == data["scheduler_mutations"] == 0
    assert not data["TS_certified"] and not data["G3_selection_or_holdout_access"]
    assert data["ordered_minus_endpoints_are_not_merged"]
    assert not data["smooth_MEP_interpolation_used"]
    assert data["distinct_manifest_native_log_STRU_identities"] == 17
    for band, energy in zip(data["bands"], (96.8281830028, 338.9689653767)):
        assert band["highest_image_index"] == 4
        assert band["rows"][4]["energy_relative_PO_meV_fu"] == pytest.approx(energy, abs=1e-8)
        assert band["ordinary_fmax_eV_A"] <= .10
        for row in band["rows"]:
            for span in plot.comparison.SPANS:
                assert row[f"{span}_internal_gradient_norm"]**2 == pytest.approx(
                    row[f"{span}_retained_gradient_norm"]**2 +
                    row[f"{span}_omitted_gradient_norm"]**2, abs=1e-12)
    peak = data["bands"][0]["rows"][4]
    assert peak["T_geometric_triplet_rank3_squared_fraction"] == pytest.approx(.826300716489)
    assert peak["T_geometric_triplet_retained_squared_gradient_fraction"] == pytest.approx(.129582585287)
    assert peak["T_geometric_triplet_internal_gradient_norm"] == pytest.approx(.0701964521014)


def test_preserving_external_work_replays_actual_nine_native_frames():
    from scripts.analyze_hfo2_clamped_external_work import analyze
    actual = analyze(ROOT/plot.TERMINAL, plot.AUDIT_SHA)
    stored = json.loads((ROOT/plot.CASE/"external_work.json").read_text())
    actual = json.loads(json.dumps(actual, allow_nan=False))
    # The existing checkout/archive retain LF/CRLF serializations of this
    # historical curvature module. Reports keep the bytes actually used;
    # fresh replay identifies its own bytes. Require exact pinned LF content
    # and known raw-byte variants, then compare every scientific field.
    curvature = "vcneb/biaxial_curvature.py"
    source = (ROOT/curvature).read_bytes()
    lf_pin = "d63eb858a07e58b176f21c5f8203b6a116b0a9c7134e204818e0441eb92bc0c7"
    crlf_pin = "d3d2d32dd2aa4736ae9ec6b0cf52d56c2d00f56c4af7e2b296fc06e58ee672ce"
    old_raw = "d63eb858a07e58b176f21c5f8203b6a116b0a9c7134e204818e0441eb92bc0c7"
    assert hashlib.sha256(source.replace(b"\r\n", b"\n")).hexdigest() == lf_pin
    assert stored["core_source_sha256"][curvature] == old_raw
    assert actual["core_source_sha256"][curvature] == hashlib.sha256(source).hexdigest()
    assert actual["core_source_sha256"][curvature] in {lf_pin, crlf_pin}
    actual["core_source_sha256"][curvature] = old_raw
    assert actual == stored
    assert actual["new_DFT_calls"] == actual["scheduler_mutations"] == 0
    assert actual["G3_selection_or_holdout_access"] is False


@pytest.mark.parametrize("pin", ["MODE_SHA", "GRADIENT_SHA", "OBSERVATION_SHA", "AUDIT_SHA"])
def test_changed_material_pin_rejected(monkeypatch, pin):
    monkeypatch.setattr(plot, pin, "0"*64)
    with pytest.raises(ValueError, match="frozen"):
        plot.build_data(ROOT)


@pytest.mark.parametrize("dependency", list(plot.DEPENDENCIES))
def test_old_adapters_remain_byte_pinned(monkeypatch, dependency):
    pins = plot.DEPENDENCIES.copy(); pins[dependency] = "0"*64
    monkeypatch.setattr(plot, "DEPENDENCIES", pins)
    with pytest.raises(ValueError, match="historical plotting dependency"):
        plot.build_data(ROOT)


def test_common_native_initial_record_not_inferred_from_arc(monkeypatch):
    original = plot.historical.build_data
    def changed(repository):
        data = original(repository)
        data["common_initial_record_canonical_sha256"] = "0"*64
        return data
    monkeypatch.setattr(plot.historical, "build_data", changed)
    with pytest.raises(ValueError, match="common initial native record"):
        plot.build_data(ROOT)


def test_six_panel_layout_all_frames_and_no_legend_overlap():
    import matplotlib.pyplot as plt
    data = plot.build_data(ROOT)
    fig = plot.make_figure(data)
    try:
        np.testing.assert_allclose(fig.get_size_inches()*25.4, [183, 195])
        positions = np.array([a.get_position().bounds for a in fig.axes]).reshape(3, 2, 4)
        np.testing.assert_allclose(positions[1:, :, [0, 2]], np.repeat(positions[:1, :, [0, 2]], 2, axis=0))
        np.testing.assert_allclose(positions[:, 1:, [1, 3]], positions[:, :1, [1, 3]])
        renderer = fig.canvas.get_renderer()
        legends = fig.legends + [a.get_legend() for a in fig.axes if a.get_legend() is not None]
        for i, ax in enumerate(fig.axes):
            assert all(s.get_visible() for s in ax.spines.values()) and not ax.get_title()
            label = ax.texts[0]
            assert label.get_text() == f"({chr(97+i)})" and label.get_weight() == "normal"
            assert label.get_fontsize() == 13
            assert ax.xaxis.get_major_ticks()[0].tick2line.get_visible()
            assert ax.yaxis.get_major_ticks()[0].tick2line.get_visible()
            for legend in legends:
                box = legend.get_window_extent(renderer)
                assert not box.overlaps(ax.get_window_extent(renderer))
                assert not box.overlaps(label.get_window_extent(renderer))
                assert legend.get_frame().get_alpha() == .85
        for i, band in enumerate(data["bands"]):
            assert fig.axes[i].get_ylim() == (-10., 370.)
            np.testing.assert_array_equal(fig.axes[i].lines[0].get_ydata(),
                [r["energy_relative_PO_meV_fu"] for r in band["rows"]])
            for offset in (2, 4):
                assert len(fig.axes[offset+i].lines) == 7
    finally:
        plt.close(fig)


def test_export_actual_rows_formats_and_no_overwrite(tmp_path):
    output = tmp_path / "fresh"
    data = plot.build_data(ROOT)
    qa = plot.export(data, output)
    with (output/"source_data.csv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 18 and [r["channel"] for r in rows[:9]] == ["preserving"]*9
    assert qa["grid_aligned"] and qa["all_four_Cmma_frames_retained"]
    assert qa["visual_review"] == "pending" and qa["independent_replicates"] == 0
    assert not qa["Cmma_spread_is_statistical_uncertainty"]
    for name, digest in qa["file_sha256"].items():
        assert plot.sha256(output/name) == digest
    assert {p.suffix for p in output.glob(plot.STEM+".*")} == {".pdf", ".svg", ".png", ".tiff"}
    assert "<text" in (output/(plot.STEM+".svg")).read_text(encoding="utf-8")
    with pytest.raises(FileExistsError, match="fresh dated"):
        plot.export(data, output)


def test_new_table5_and6_values_are_real_samples_not_TS_or_energy_accuracy():
    text = (ROOT/"paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md").read_text(encoding="utf-8")
    assert "| Preserving-pattern switching, same clamping | 96.828183 | 96.863442 | -0.035259 |" in text
    assert "| Preserving-pattern switching | T geometric triplet | 3 | 82.6301 | 12.9583 |" in text
    assert "| Preserving-pattern switching | Registered Cmma span | 4 | 14.8336–14.8337 | 52.7203 |" in text
    data = plot.build_data(ROOT)
    peak = data["bands"][0]["rows"][4]
    assert 100*peak["T_optical_doublets_rank4_squared_fraction"] < .0001
    assert 100*peak["T_two_lowest_complete_optical_doublets_retained_squared_gradient_fraction"] < .0001
    section = text.split("### 3.9. ", 1)[1].split("## 4.", 1)[0]
    assert "not a certified index-one saddle" in section
    assert "the globally easiest homogeneous switching mechanism" in section
    assert "Five of eight registered G2 training bands" in section
    assert "Eighteen frame records share the initial" in section
    for suffix in ("hfo2_clamped_flips.pdf", "hfo2_clamped_flips.png"):
        assert (ROOT/"paper/VARNEB_JCTC/figures/hfo2_clamped_flips_E101_20261011"/suffix).is_file()
