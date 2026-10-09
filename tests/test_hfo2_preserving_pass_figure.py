"""Actual frozen material and fail-closed captions for the dated G1 figure."""

import copy
import json
from pathlib import Path

import pytest

from scripts import plot_hfo2_G1_preserving_pass as plot

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255"


@pytest.fixture(scope="module")
def data():
    return plot.build_data(CASE / "network_specification.json",
                           CASE / "network_progress_local.json", CASE / "analysis_local.json")


def test_actual_stage_has_three_passes_not_a_TS_certificate(data):
    assert [c["snapshot_step"] for c in data["channels"]] == [6, 39, 69, 32]
    assert [c["ordinary_residual_passed"] for c in data["channels"]] == [True, True, True, False]
    assert data["existing_image_records_including_cached_duplicates"] == 37
    assert data["new_DFT_calls"] == 0
    assert not data["full_G1_passed"] and not data["TS_certified"]
    assert not data["smooth_MEP_interpolation_used"]
    preserving = data["preserving"]
    assert preserving["replayed_fmax_eV_A"] == pytest.approx(.09937245841139401)
    assert preserving["images"][4]["relative_energy_meV_fu"] == pytest.approx(-14.838200784197397)
    assert preserving["images"][3]["true_tangential_euclidean_eV_A"] == pytest.approx(.24248563953626004)


def test_csv_keeps_both_peaks_center_and_no_endpoint_residual(data):
    rows = [r for r in data["rows"] if r["channel"] == "PO_flip_T_pattern_preserving"]
    assert rows[3]["energy_relative_PO_meV_fu"] == pytest.approx(32.806023216835456)
    assert rows[3]["space_group_symprec_0.001_A"] == "Pca2_1"
    assert rows[4]["space_group_symprec_0.001_A"] == "Pbcn"
    assert rows[5]["space_group_symprec_0.001_A"] == "Pca2_1"
    assert rows[5]["preserving_true_tangent_eV_A"] < -.24
    for r in (rows[0], rows[-1]):
        assert r["preserving_true_tangent_eV_A"] is None
        assert r["NEB_max_vector_eV_A"] is None
    assert all(r["preserving_true_tangent_eV_A"] is None
               for r in data["rows"] if r["channel"] != "PO_flip_T_pattern_preserving")


@pytest.mark.parametrize("bad", ["step", "pass", "gate"])
def test_new_stage_or_unmeasured_gate_cannot_keep_this_caption(monkeypatch, data, bad):
    report = json.loads((CASE / "network_progress_local.json").read_text())
    if bad == "step":
        report["channels"][2]["snapshot_step"] = 70
    elif bad == "pass":
        report["channels"][3]["ordinary_residual_passed"] = True
    else:
        report["TS_certified"] = True
    monkeypatch.setattr(plot, "analyze_progress", lambda _: copy.deepcopy(report))
    monkeypatch.setattr(plot, "verify_material_replay", lambda *_: [])
    with pytest.raises(ValueError, match="stage|unmeasured"):
        plot.build_data(CASE / "network_specification.json",
                        CASE / "network_progress_local.json", CASE / "analysis_local.json")


def test_changed_chain_report_rejected(tmp_path, monkeypatch, data):
    report = json.loads((CASE / "analysis_local.json").read_text())
    report["analysis_script_sha256"] = "wrong"
    file = tmp_path / "stale.json"
    file.write_text(json.dumps(report))
    monkeypatch.setattr(plot, "analyze_progress", lambda _: json.loads((CASE / "network_progress_local.json").read_text()))
    monkeypatch.setattr(plot, "verify_material_replay", lambda *_: [])
    with pytest.raises(ValueError, match="hash"):
        plot.build_data(CASE / "network_specification.json", CASE / "network_progress_local.json", file)


def test_layout_is_aligned_and_labels_are_normal_not_titles(data):
    import matplotlib.pyplot as plt
    fig = plot.make_figure(data)
    try:
        positions = [a.get_position().bounds for a in fig.axes]
        for i in (0, 2):
            assert positions[i][1] == positions[i+1][1]
        for i in (0, 1):
            assert positions[i][0] == positions[i+2][0]
        for i, axis in enumerate(fig.axes):
            assert not axis.get_title()
            assert all(s.get_visible() for s in axis.spines.values())
            label = next(t for t in axis.texts if t.get_text() == f"({chr(97+i)})")
            assert label.get_fontweight() == "normal"
        assert fig.legends[0].get_frame().get_alpha() == .85
        assert len(fig.axes[0].lines) == 4  # three profiles and zero baseline
        assert len(fig.axes[1].lines) == 2  # reversing profile and zero baseline
    finally:
        plt.close(fig)


def test_existing_export_refused_before_input_read(tmp_path, monkeypatch):
    output = tmp_path / "keep"
    output.mkdir()
    marker = output / "original.txt"
    marker.write_text("retain")
    monkeypatch.setattr("sys.argv", ["plot", "--specification", "missing", "--material-report", "missing",
                                    "--chain-report", "missing", "--output", str(output)])
    with pytest.raises(FileExistsError):
        plot.main()
    assert marker.read_text() == "retain"
