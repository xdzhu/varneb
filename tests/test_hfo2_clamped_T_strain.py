"""Actual completed +0.01 training, unchanged reference and no holdout draw."""
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from scripts import plot_hfo2_clamped_T_strain as plot
from scripts import plot_hfo2_clamped_comparison as paired
from scripts.analyze_hfo2_clamped_external_work import analyze
from scripts.plot_hfo2_clamped_M import build_band_data

ROOT = Path(__file__).resolve().parents[1]


def test_actual_strained_receipt_and_raw_stress_replay():
    assert hashlib.sha256((ROOT/plot.AUDIT).read_bytes()).hexdigest() == plot.AUDIT_SHA
    result = analyze(ROOT/plot.AUDIT.parent, plot.AUDIT_SHA)
    assert result["source_job_id"] == "28709788"
    assert result["strain_training_condition"] == .01
    assert result["ordinary_fmax_eV_A"] == pytest.approx(.0949479822061884)
    assert result["sampled_peak_image"] == 4
    assert result["sampled_peak"]["fixed_open_row_relative_PO_partial_meV_fu_per_unit_strain"] == pytest.approx(1366.1298301198235)
    assert result["new_DFT_calls"] == result["scheduler_mutations"] == 0
    assert not result["G3_selection_or_holdout_access"]


def test_same_fixed_reference_not_same_physical_PO_cache():
    data = plot.build_data(ROOT)
    assert data["training_strains"] == [0., .01]
    assert not data["different_strains_reuse_same_initial_SCF"]
    assert data["energies_relative_to_each_matched_PO"]
    assert data["observed_sampled_forward_peak_increment_meV_fu"] == pytest.approx(12.347026513452874)
    assert data["observed_PO_energy_increment_eV_cell"] == pytest.approx(.027226009508)
    assert data["unit_normal_serialization_max_difference"] < 1e-16
    assert not data["G3_selection_or_holdout_access"]
    assert not data["TS_certified"] and not data["model_accuracy_certified"]
    assert not data["independent_prediction"] and not data["native_data_reaudit_claimed"]
    assert data["new_DFT_calls"] == data["scheduler_mutations"] == 0
    zero, one = data["bands"]
    z, s = zero["rows"][4], one["rows"][4]
    assert z["T_geometric_triplet_rank3_squared_fraction"] == pytest.approx(.9909572569571496)
    assert s["T_geometric_triplet_rank3_squared_fraction"] == pytest.approx(.990044681252212)
    assert s["T_geometric_triplet_retained_squared_gradient_fraction"] == pytest.approx(.6972449231451816)
    assert s["Cmma_frame_103_retained_squared_gradient_fraction"] == pytest.approx(.1057660783721135)
    assert s["T_geometric_triplet_internal_gradient_norm"] == pytest.approx(.038221368665941294)
    assert s["Green_E_00_original_T"] == pytest.approx(.01005)
    assert s["Green_E_22_original_T"] == pytest.approx(.01005)
    assert zero["zero_displacement"][-1] is True
    assert not any(one["zero_displacement"])
    for band in data["bands"]:
        assert len(band["rows"]) == 9
        for row in band["rows"]:
            for span in paired.SPANS:
                assert row[f"{span}_internal_gradient_norm"]**2 == pytest.approx(
                    row[f"{span}_retained_gradient_norm"]**2 + row[f"{span}_omitted_gradient_norm"]**2, abs=1e-12)


@pytest.mark.parametrize("strain", [.005, -.01, .02])
def test_reserved_or_unregistered_strain_rejected_before_reading(strain, tmp_path):
    with pytest.raises(ValueError, match="registered zero"):
        build_band_data(tmp_path, Path("missing"), Path("missing"), Path("missing"),
            mode_digest="x", observation_digest="x", audit_digest="x", job="x", step=0,
            expected_strain=strain)


@pytest.mark.parametrize("relative", [plot.REPORT, plot.GRADIENT, plot.OBSERVATION, plot.AUDIT])
def test_changed_material_inputs_reject_before_export(relative, tmp_path):
    # Enough frozen inputs to load the first known band and the next source.
    paths = (paired.T_REPORT, paired.T_GRADIENT, paired.T_OBSERVATION, paired.T_AUDIT,
             plot.REPORT, plot.GRADIENT, plot.OBSERVATION, plot.AUDIT)
    for path in paths:
        dest = tmp_path/path
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/path, dest)
    value = json.loads((tmp_path/relative).read_text())
    value["tampered_training_input"] = True
    (tmp_path/relative).write_text(json.dumps(value))
    with pytest.raises(ValueError, match="frozen material|pinned completed"):
        plot.build_data(tmp_path)


def test_same_grid_style_and_distinct_defined_endpoint():
    import matplotlib.pyplot as plt
    fig = paired.make_figure(plot.build_data(ROOT))
    try:
        positions = np.array([ax.get_position().bounds for ax in fig.axes]).reshape(3, 2, 4)
        np.testing.assert_allclose(positions[1:, :, [0, 2]], np.repeat(positions[0:1, :, [0, 2]], 2, axis=0))
        np.testing.assert_allclose(positions[:, 1, [1, 3]], positions[:, 0, [1, 3]])
        renderer = fig.canvas.get_renderer()
        legends = fig.legends + [ax.get_legend() for ax in fig.axes if ax.get_legend() is not None]
        for i, ax in enumerate(fig.axes):
            assert not ax.get_title() and all(s.get_visible() for s in ax.spines.values())
            assert ax.texts[0].get_text() == f"({chr(97+i)})" and ax.texts[0].get_weight() == "normal"
            for legend in legends:
                box = legend.get_window_extent(renderer)
                assert not box.overlaps(ax.get_window_extent(renderer))
                assert not box.overlaps(ax.texts[0].get_window_extent(renderer))
                assert legend.get_frame().get_alpha() == .85
        assert "0\\%" in fig.axes[0].get_legend().texts[0].get_text()
        assert "1\\%" in fig.axes[1].get_legend().texts[0].get_text()
        assert np.isnan(fig.axes[2].lines[0].get_ydata()[-1])
        assert np.isfinite(fig.axes[3].lines[0].get_ydata()[-1])
        assert len(fig.axes[3].lines) == 7  # all four registered Cmma frames
    finally:
        plt.close(fig)


def test_export_refuses_existing_directory_and_path_stem(tmp_path):
    with pytest.raises(FileExistsError):
        paired.export(plot.build_data(ROOT), tmp_path, stem_name="hfo2_clamped_T_strain")
    with pytest.raises(ValueError, match="plain extension-free"):
        paired.export({}, tmp_path/"new", stem_name="../elsewhere")
