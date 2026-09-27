"""The alternate BTO frozen plane must stay distinct from a conditional PES."""

import hashlib
import json

import pytest

from scripts.plot_bto_paper_path_adapted_79 import (
    FIGURES, GRID, PATH, SOURCE_QA, interpolation_diagnostics, read_data,
)


def test_snapshot_has_79_real_dft_nodes_and_seven_path_images():
    nodes, images, source_qa = read_data()
    assert len(nodes) == 79
    assert len(images) == 7
    assert sum(row["sample_group"] == "nonuniform_refinement_16_DFT" for row in nodes) == 16
    assert all(len(row["log_sha256"]) == 64 for row in nodes)
    assert source_qa["contour_claim_status"] == "exploratory_not_globally_validated"
    assert float(images[0]["Q1_sqrt_amu_A"]) - max(
        float(row["Q1_sqrt_amu_A"]) for row in nodes
    ) == pytest.approx(0.0043298776588884)
    assert source_qa["maximum_atomic_projection_residual_sqrt_amu_A"] == pytest.approx(
        0.060278288479626345
    )


def test_interpolated_pixels_are_not_counted_as_dft_and_holdout_is_recomputed():
    nodes, _, _ = read_data()
    _, _, _, diagnostics = interpolation_diagnostics(nodes)
    assert diagnostics["n_real_DFT_points"] == 79
    assert diagnostics["n_display_pixels_not_DFT"] == 401 * 201
    assert diagnostics["preregistered_16_point_cubic_max_abs_error_meV_per_BTO"] == pytest.approx(
        7.191791931629894, abs=1e-8
    )
    assert diagnostics["new_16_point_leave_one_out_cubic_max_abs_error_meV_per_BTO"] == pytest.approx(
        6.952572958074326, abs=1e-8
    )


def test_figure_qa_binds_sources_and_warns_against_mixed_plane_claim():
    report = json.loads((FIGURES / "bto_frozen_path_adapted_79_qa.json").read_text(encoding="utf-8"))
    assert report["axis_mode_groups_zero_based"] == [[0, 1, 2], [6, 7, 8]]
    assert report["status"] == "exploratory_frozen_soft_stable_plane_not_conditional_PES_or_barrier"
    assert "100 Ry" in report["calculator"] and "10 au DZP" in report["calculator"]
    assert report["T_endpoint_Q1_outside_sampled_range_sqrt_amu_A"] == pytest.approx(0.0043298776588884)
    for path in (GRID, PATH, SOURCE_QA):
        assert report["source_sha256"][path.name] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert any("different from the two-soft-mode" in note for note in report["limitations"])
