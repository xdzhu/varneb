from pathlib import Path
import json
import shutil

import numpy as np
import pytest

from scripts import plot_hfo2_sampling_bridge as plot

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "benchmarks/hfo2_channels/20261008/preserving_sampling_bridge_20261009/completed_HF"
OBSERVATION = ROOT / "benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255/preserving_step69"


def test_actual_all_six_statics_and_original_samples_replay_without_DFT():
    data = plot.build_data(OBSERVATION, CASE)
    assert len(data["rows"]) == 15
    assert data["new_static_SCF_count"] == 6 and data["new_DFT_calls_for_plot"] == 0
    assert data["old_max_meV_fu"] == pytest.approx(32.806023216835456)
    assert data["new_max_meV_fu"] > data["old_max_meV_fu"]
    assert not data["stationary_TS_or_relaxed_MEP_certified"]
    assert data["reconstructed_band"]["total_images"] == 15
    assert data["reconstructed_band"]["moving_images"] == 13
    assert data["reconstructed_band"]["new_SCF_calls"] == 0
    assert data["reconstructed_band"]["optimizer_steps"] == 0
    assert np.isfinite(data["reconstructed_band"]["replayed_ordinary_fmax_eV_A"])
    new = data["rows"][9:]
    assert [r["segment"] for r in new] == ["2->3"]*3+["5->6"]*3
    assert [r["fraction"] for r in new] == [.25, .50, .75]*2
    np.testing.assert_allclose(data["rows"][0]["s_normalized"], 0.)
    np.testing.assert_allclose(data["rows"][8]["s_normalized"], 1.)


@pytest.mark.parametrize("key,value", [("climb", True), ("pressure_GPa", 1.), ("formula_units", 12)])
def test_physical_or_TS_interpretation_changed_is_rejected(tmp_path, key, value):
    copied = tmp_path / "copied"
    shutil.copytree(CASE, copied)
    m = json.loads((copied / "manifest.json").read_text())
    m[key] = value
    (copied / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(ValueError, match="contract changed"):
        plot.build_data(OBSERVATION, copied)


def test_changed_raw_log_rejected_even_if_parsed_results_unchanged(tmp_path):
    copied = tmp_path / "copied"
    shutil.copytree(CASE, copied)
    log = copied / "calculations/00/scf_000000/OUT.ABACUS/running_scf.log"
    log.write_bytes(log.read_bytes()+b"\nchanged provenance\n")
    with pytest.raises(ValueError, match="provenance changed"):
        plot.build_data(OBSERVATION, copied)


def test_actual_figure_has_aligned_axes_plain_labels_and_framed_legends():
    import matplotlib.pyplot as plt
    data = plot.build_data(OBSERVATION, CASE)
    fig = plot.make_figure(data)
    axes = fig.axes
    assert [a.texts[0].get_text() for a in axes] == ["(a)", "(b)"]
    assert all(a.texts[0].get_fontweight() == "normal" for a in axes)
    assert all(not a.get_title() for a in axes)
    assert all(s.get_visible() for a in axes for s in a.spines.values())
    assert axes[0].get_position().y0 == axes[1].get_position().y0
    assert axes[0].get_position().y1 == axes[1].get_position().y1
    assert all(a.get_legend().get_frame().get_alpha() == .85 for a in axes)
    assert len(axes[0].collections[0].get_offsets()) == 6
    assert all(len(line.get_xdata()) == 5 for line in axes[1].lines[:2])
    plt.close(fig)


def test_export_retains_all_actual_source_rows_and_no_smooth_MEP_claim(tmp_path):
    output = tmp_path / "figure"
    qa = plot.export(plot.build_data(OBSERVATION, CASE), output)
    assert qa["new_DFT_calls_for_plot"] == 0
    assert "not smooth MEP" in qa["line_interpretation"]
    assert qa["statistics"].startswith("9 reused band samples and 6 new SCFs")
    assert len((output / "source_data.csv").read_text().splitlines()) == 16
    assert all((output / ("hfo2_G1_sampling_bridge."+ext)).exists() for ext in ("svg", "pdf", "png", "tiff"))
    with pytest.raises(FileExistsError):
        plot.export({}, output)
