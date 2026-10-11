"""Native terminal samples, original frames and honest cross-boundary semantics."""
from pathlib import Path

import numpy as np
import pytest

from scripts import plot_hfo2_preserving_boundary_role as plot

ROOT = Path(__file__).resolve().parents[1]


def test_actual_two_centres_and_not_a_strain_derivative():
    data = plot.build_data(ROOT)
    assert data["image_records"] == 18
    assert data["new_DFT_calls"] == data["scheduler_mutations"] == 0
    assert data["same_electronic_contract"] and not data["same_mechanical_boundary"]
    assert not data["common_initial_well_used"]
    assert not data["G3_selection_or_holdout_access"]
    assert not data["TS_or_stable_intermediate_certified"]
    free, clamped = data["bands"]
    assert free["central_sample"]["energy_relative_own_PO_meV_fu"] == pytest.approx(-14.8382007842)
    assert clamped["central_sample"]["energy_relative_own_PO_meV_fu"] == pytest.approx(96.8281830028)
    assert free["centre_minus_mean_adjacent_meV_fu"] < 0
    assert clamped["centre_minus_mean_adjacent_meV_fu"] > 0
    assert free["central_sample"]["T_F_xx"] == pytest.approx(1.12202198270)
    assert clamped["central_sample"]["T_F_xx"] == pytest.approx(1)
    assert data["free_union_sampling_peak_meV_fu"] == pytest.approx(38.5084480845)
    assert free["sampled_peak_meV_fu"] == pytest.approx(32.8060232168)
    assert free["original_cell_scale_A"] == pytest.approx(5.131328250206249)
    assert clamped["original_cell_scale_A"] == pytest.approx(5.12968067458423)
    for band in data["bands"]:
        assert band["ordinary_fmax_eV_A"] <= .10
        assert len(band["rows"]) == 9
        assert band["rows"][0]["energy_relative_own_PO_meV_fu"] == 0
        assert band["rows"][0]["s_normalized"] == 0
        assert band["rows"][-1]["s_normalized"] == 1
        centre = band["central_sample"]
        assert [centre[k] for k in ("symprec_0p001", "symprec_0p01", "symprec_0p05")] == ["Pbcn"]*3
        for row in band["rows"]:
            audit = row["structure_audit"]
            assert not audit["standardized_wrapped_reordered_or_relaxed"]
            assert not audit["oxidation_states_guessed"]
            assert audit["minimum_periodic_distance_A"] >= 1.6
            assert [s["angle_tolerance_deg"] for s in audit["symmetry_sweep"]] == [1.]*3
    np.testing.assert_allclose([r["T_F_xx"] for r in clamped["rows"]], 1, atol=1e-12, rtol=0)
    np.testing.assert_allclose([r["T_F_zz"] for r in clamped["rows"]], 1, atol=1e-12, rtol=0)


@pytest.mark.parametrize("name", list(plot.PINS))
def test_wrong_source_receipt_fails_closed(monkeypatch, name):
    monkeypatch.setitem(plot.PINS, name, "0"*64)
    with pytest.raises(ValueError, match="digest"):
        plot.build_data(ROOT)


def test_figure_alignment_and_outside_legend():
    import matplotlib.pyplot as plt
    fig = plot.make_figure(plot.build_data(ROOT))
    try:
        np.testing.assert_allclose(fig.get_size_inches()*25.4, [183, 88])
        a, b = fig.axes
        np.testing.assert_allclose(a.get_position().bounds[1::2], b.get_position().bounds[1::2])
        assert [ax.texts[0].get_text() for ax in fig.axes] == ["(a)", "(b)"]
        for ax in fig.axes:
            assert all(s.get_visible() for s in ax.spines.values())
            assert ax.texts[0].get_fontweight() == "normal"
            assert ax.get_title() == ""
            assert len(ax.lines) >= 4
        renderer = fig.canvas.get_renderer()
        legend = fig.legends[0].get_window_extent(renderer)
        assert all(not legend.overlaps(ax.get_window_extent(renderer)) for ax in fig.axes)
        assert fig.legends[0].get_frame().get_alpha() == .85
    finally:
        plt.close(fig)


def test_fresh_output_only(tmp_path):
    with pytest.raises(FileExistsError, match="fresh"):
        plot.export({}, tmp_path)


def test_manuscript_does_not_certify_the_two_centres():
    text = (ROOT/"paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md").read_text(encoding="utf-8")
    section = text.split("### 3.10.")[1].split("## 4.")[0]
    for value in ("-14.838201", "+96.828183", "1.122022", "38.508448", "Figure 8."):
        assert value in section
    for limit in ("different, boundary-relaxed PO+ well", "not present on this clamped",
                  "stationary minimum or index-one saddle", "not a structural matching",
                  "independent error-bounded forecasts remain unmeasured"):
        assert limit in section
