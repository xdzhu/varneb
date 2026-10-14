"""Real pinned paired observations, no author-mode redistribution or new DFT."""
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from scripts import plot_hfo2_clamped_comparison as plot
from scripts.analyze_hfo2_clamped_external_work import ADMITTED_TRAINING_AUDITS, analyze

ROOT = Path(__file__).resolve().parents[1]


def test_actual_T_is_admitted_and_its_native_work_replays():
    reversing_sha = "632e0e9b2b821b3bfc7d66b361c44c7db734988fb5d018ea2f31441fab14ed41"
    reversing_receipt = ROOT / "benchmarks/hfo2_channels/20261008/clamped_reversing_terminal_E083_20261011/audit_receipt.json"
    assert hashlib.sha256(reversing_receipt.read_bytes()).hexdigest() == reversing_sha
    strained_sha = "32a6e49a618e6e62632ae71a5e0ed093095254f6858221be2ade3d83cc194a28"
    strained_receipt = ROOT / "benchmarks/hfo2_channels/20261008/clamped_plus1_T_terminal_E087_20261011/audit_receipt.json"
    assert hashlib.sha256(strained_receipt.read_bytes()).hexdigest() == strained_sha
    preserving_sha = "46145da14f24b2b8bd46aa4656f0982650d6b2dc661282f58d809aab8a6d3202"
    preserving_receipt = ROOT / "benchmarks/hfo2_channels/20261008/clamped_held_flips_E100_20261011/zero_strain_preserving/audit_receipt.json"
    assert hashlib.sha256(preserving_receipt.read_bytes()).hexdigest() == preserving_sha
    plus1_M_sha = "b86cdb72eab9e3d8ff2125e16e6fbb2d0954b34d54c0f05099a581c066108389"
    plus1_M_receipt = ROOT / "benchmarks/hfo2_channels/20261008/clamped_reversing_release_E104_20261012/plus1_M_terminal/audit_receipt.json"
    assert hashlib.sha256(plus1_M_receipt.read_bytes()).hexdigest() == plus1_M_sha
    plus1_P_sha = '0c8d91532bc9051decdfd600eae12ce409fd84a86c677158077948fa352c9d8a'
    plus1_P_receipt = ROOT/'benchmarks/hfo2_channels/20261008/clamped_plus1_preserving_terminal_E145_20261012/completed_HF/audit_receipt.json'
    assert hashlib.sha256(plus1_P_receipt.read_bytes()).hexdigest() == plus1_P_sha
    # E145 admits only the seventh actual completed training band.
    # Keep an exact byte-pinned set, not arbitrary or holdout permission.
    assert ADMITTED_TRAINING_AUDITS == {plot.AUDIT_SHA, plot.T_AUDIT_SHA, reversing_sha, strained_sha, preserving_sha, plus1_M_sha, plus1_P_sha}
    result = analyze(ROOT / plot.T_AUDIT.parent, plot.T_AUDIT_SHA)
    assert result["source_job_id"] == "28692775"
    assert result["sampled_peak_image"] == 4
    assert result["ordinary_fmax_eV_A"] == pytest.approx(.09358976572425506)
    assert result["native_portable_frames_reparsed"] == 9
    assert result["new_DFT_calls"] == result["scheduler_mutations"] == 0
    assert result["G3_selection_or_holdout_access"] is False


def test_actual_same_reference_bands_have_distinct_displacement_and_gradient_coverage():
    data = plot.build_data(ROOT)
    assert data["common_initial_energy_eV_cell"] == pytest.approx(-9782.973830941437)
    assert data["new_DFT_calls"] == data["scheduler_mutations"] == 0
    assert not data["TS_certified"] and not data["G3_selection_or_holdout_access"]
    assert not data["native_data_reaudit_claimed"]
    t, m = data["bands"]
    assert t["channel"] == "T" and m["channel"] == "M"
    assert t["highest_image_index"] == 4 and m["highest_image_index"] == 3
    pt, pm = t["rows"][4], m["rows"][3]
    assert pt["energy_relative_PO_meV_fu"] == pytest.approx(84.41662388122495)
    assert pm["energy_relative_PO_meV_fu"] == pytest.approx(98.45812195362669)
    assert pt["T_geometric_triplet_rank3_squared_fraction"] == pytest.approx(.9909572569571496)
    assert pt["T_geometric_triplet_retained_squared_gradient_fraction"] == pytest.approx(.14301312220712437)
    assert pt["T_two_lowest_complete_optical_doublets_retained_squared_gradient_fraction"] == pytest.approx(.2039633700619102)
    assert pt["Cmma_frame_103_retained_squared_gradient_fraction"] == pytest.approx(.5407648515179279)
    assert pm["Cmma_frame_103_retained_squared_gradient_fraction"] == pytest.approx(.9350846598522286)
    assert pt["T_geometric_triplet_internal_gradient_norm"] == pytest.approx(.010196606055786928)
    for band in (t, m):
        assert len(band["rows"]) == 9
        for row in band["rows"]:
            for span in plot.SPANS:
                internal = row[f"{span}_internal_gradient_norm"]
                retained = row[f"{span}_retained_gradient_norm"]
                omitted = row[f"{span}_omitted_gradient_norm"]
                assert internal**2 == pytest.approx(retained**2 + omitted**2, abs=1e-12)
    assert t["rows"][-1]["zero_T_displacement_denominator"] is True
    assert t["rows"][4]["Green_E_11_original_T"] == pytest.approx(-.00391575926852783)


def test_changed_T_report_rejects_before_plotting(tmp_path):
    relatives = (plot.T_REPORT, plot.T_GRADIENT, plot.T_OBSERVATION, plot.T_AUDIT)
    for relative in relatives:
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, target)
    changed = json.loads((tmp_path / plot.T_REPORT).read_text())
    changed["T_triplet_rank3_fraction"][4] = .9
    (tmp_path / plot.T_REPORT).write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="frozen material sources"):
        plot.build_data(tmp_path)


def test_grid_labels_legends_and_undefined_endpoint_are_honest():
    import matplotlib.pyplot as plt
    fig = plot.make_figure(plot.build_data(ROOT))
    try:
        assert len(fig.axes) == 6
        np.testing.assert_allclose(fig.get_size_inches()*25.4, [183, 202])
        positions = np.array([ax.get_position().bounds for ax in fig.axes]).reshape(3, 2, 4)
        np.testing.assert_allclose(positions[1:, :, [0, 2]], np.repeat(positions[0:1, :, [0, 2]], 2, axis=0))
        np.testing.assert_allclose(positions[:, 1, [1, 3]], positions[:, 0, [1, 3]])
        renderer = fig.canvas.get_renderer()
        legends = fig.legends + [ax.get_legend() for ax in fig.axes if ax.get_legend() is not None]
        for i, ax in enumerate(fig.axes):
            assert all(s.get_visible() for s in ax.spines.values()) and not ax.get_title()
            assert ax.texts[0].get_text() == f"({chr(97+i)})"
            assert ax.texts[0].get_weight() == "normal"
            assert ax.get_ylim() == fig.axes[i//2*2].get_ylim()
            label_box = ax.texts[0].get_window_extent(renderer)
            for legend in legends:
                box = legend.get_window_extent(renderer)
                assert not box.overlaps(ax.get_window_extent(renderer))
                assert not box.overlaps(label_box)
                assert legend.get_frame().get_alpha() == .85
        # Three representative curves plus all four Cmma frames and peak guide.
        assert len(fig.axes[2].lines) == len(fig.axes[4].lines) == 7
        assert np.isnan(fig.axes[2].lines[0].get_ydata()[-1])
        assert np.isnan(fig.axes[2].lines[2].get_ydata()[-1])
        assert np.isfinite(fig.axes[4].lines[0].get_ydata()[-1])
    finally:
        plt.close(fig)


def test_no_overwrite_export(tmp_path):
    with pytest.raises(FileExistsError, match="fresh dated"):
        plot.export(plot.build_data(ROOT), tmp_path)
