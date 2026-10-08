"""Publication plotting must not manufacture a saddle or a complete mode plane."""

import copy
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from scripts.plot_hfo2_converged_band import (
    assert_replay_equal, build_data, export, make_figure, optical_group_weights,
)
from scripts.audit_hfo2_static_replica import sha256

ROOT = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"


@pytest.fixture(scope="module")
def real_data():
    return build_data(ROOT / "converged_gap", ROOT / "reference_variants",
                      ROOT / "gamma_analysis/T_d0.01.npz")


def test_real_figure_data_replay_and_descriptive_limits(real_data):
    rows = real_data["rows"]
    assert len(rows) == 10 and real_data["peak_image_index"] == 3
    assert real_data["largest_residual_image_index"] == 2
    assert rows[3]["energy_relative_T_meV_fu"] == pytest.approx(33.88897141439884)
    assert real_data["replayed_fmax_eV_A"] == pytest.approx(.05988161569025288)
    assert real_data["Gamma_selected_groups"][0]["mode_indices_zero_based"] == [13, 14]
    assert real_data["Gamma_selected_groups"][1]["mode_indices_zero_based"] == [19]
    assert real_data["peak_triplet_captured_fraction"] == pytest.approx(.9553920684674048)
    assert real_data["PO_triplet_captured_fraction"] == pytest.approx(.31433542834208184)
    assert not real_data["TS_certified"] and not real_data["energy_partition_or_prediction_claimed"]
    assert real_data["new_DFT_calls"] == 0 and not real_data["parameters_changed"]


def test_reference_zero_and_fixed_endpoint_residuals_are_absent_not_zero(real_data):
    rows = real_data["rows"]
    assert rows[0]["Gamma_optical_squared_norm_amu_A2"] == 0
    assert all(rows[0][key] is None for key in ("Gamma_group1_fraction", "Gamma_group2_fraction", "Gamma_other_fraction"))
    for i in (0, -1):
        assert rows[i]["NEB_atomic_max_vector_eV_A"] is None
        assert rows[i]["NEB_cell_max_vector_eV_A"] is None
    for row in rows[1:]:
        assert sum(row[key] for key in ("Gamma_group1_fraction", "Gamma_group2_fraction", "Gamma_other_fraction")) == pytest.approx(1)


def test_complete_doublet_weights_are_rotation_invariant():
    analysis = json.loads((ROOT / "converged_gap/analysis.json").read_text())
    obs = analysis["observations"][0]
    q = np.array([r["T_Gamma_Q_sqrt_amu_A"] for r in obs["images"]])
    args = (analysis["T_Gamma_frequencies_THz"], q, obs["Gamma_acoustic_mode_indices"], 3)
    groups, weights, total = optical_group_weights(*args)
    rotated = q.copy()
    angle = .371
    rotated[:, [13, 14]] = q[:, [13, 14]] @ np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    other_groups, other_weights, other_total = optical_group_weights(args[0], rotated, args[2], 3)
    assert groups == other_groups
    assert np.allclose(weights, other_weights, equal_nan=True, atol=1e-14)
    assert np.allclose(total, other_total, atol=1e-13)


@pytest.mark.parametrize("change", ("number", "nan", "hash", "missing"))
def test_stale_or_nonfinite_descriptors_are_rejected(change):
    original = {"curve": [0., 1.], "source": "original"}
    changed = copy.deepcopy(original)
    if change == "number":
        changed["curve"][1] += .01
    elif change == "nan":
        changed["curve"][1] = float("nan")
    elif change == "hash":
        changed["source"] = "other"
    else:
        changed.pop("source")
    with pytest.raises(ValueError, match="stale"):
        assert_replay_equal(original, changed)


def test_user_panel_style_and_grid_alignment(real_data):
    import matplotlib.pyplot as plt
    fig = make_figure(real_data)
    try:
        assert len(fig.axes) == 6
        assert all(not axis.get_title() for axis in fig.axes)
        assert all(spine.get_visible() for axis in fig.axes for spine in axis.spines.values())
        positions = [axis.get_position().bounds for axis in fig.axes]
        assert all(positions[i][1] == pytest.approx(positions[i+1][1]) for i in (0, 2, 4))
        assert all(positions[i][0] == pytest.approx(positions[i+2][0]) for i in (0, 1, 2, 3))
        for i, axis in enumerate(fig.axes):
            label = next(text for text in axis.texts if text.get_text() == f"({chr(97+i)})")
            assert label.get_fontweight() == "normal" and label.get_fontsize() == 13
            legend = axis.get_legend()
            if legend:
                assert legend.get_frame().get_alpha() == .85
        assert all(axis.xaxis.get_label().get_fontsize() == 11 for axis in fig.axes[-2:])
    finally:
        plt.close(fig)


def test_export_refuses_overwriting_existing_directory(tmp_path):
    with pytest.raises(FileExistsError, match="existing"):
        export({}, tmp_path)


def test_checked_export_bundle_keeps_hashes_editable_text_and_limits():
    folder = Path(__file__).resolve().parents[1] / "paper/VARNEB_JCTC/figures/hfo2_T_PO_ordinary_20261008"
    qa = json.loads((folder / "qa.json").read_text())
    assert sha256(folder / "source_data.csv") == qa["source_data_sha256"]
    assert sha256(Path(__file__).resolve().parents[1] / "scripts/plot_hfo2_converged_band.py") == qa["plot_script_sha256"]
    for name, digest in qa["file_sha256"].items():
        assert sha256(folder / name) == digest
    assert qa["row_alignment_verified"] and qa["column_alignment_verified"]
    assert "inspected" in qa["visual_review"]
    assert not qa["TS_certified"] and not qa["energy_partition_or_prediction_claimed"]
    root = ET.parse(folder / "hfo2_T_PO_ordinary_modes.svg").getroot()
    assert len(root.findall(".//{http://www.w3.org/2000/svg}text")) == 61
