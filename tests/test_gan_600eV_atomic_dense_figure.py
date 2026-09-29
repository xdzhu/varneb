"""The GaN contour uses only a validated bounded interpolant."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from matplotlib.collections import PathCollection
import matplotlib.pyplot as plt

from scripts.plot_gan_600eV_atomic_dense_surface import load, model_surface
from scripts.plot_gan_600eV_atomic_transverse_landscape import transverse_surface
from scripts.plot_gan_600eV_atomic_transverse_q9 import draw as draw_q9, load as load_q9


ROOT = Path(__file__).resolve().parents[1]


def test_dense_figure_refuses_unaudited_report(tmp_path) -> None:
    report = tmp_path / "report.json"
    trajectory = tmp_path / "path.traj"
    report.write_text(json.dumps({"status": "inputs_finalized_no_DFT"}), encoding="utf-8")
    trajectory.write_bytes(b"not a trajectory")
    with pytest.raises(ValueError, match="raw audit or prospective contour gate"):
        load(report, trajectory)


def test_dense_model_preserves_path_and_outer_transverse_points() -> None:
    arc = np.linspace(0.15, 0.80, 18)
    h0 = -100 + 0.3 * np.sin(np.linspace(0, np.pi, 18))
    width = 0.05
    slope = np.linspace(-5, 5, 18)
    halfcurv = np.linspace(50, 100, 18)
    matrix = np.array([
        [slope[index] * q + halfcurv[index] * q**2
         for q in (-width, -0.025, 0.0, 0.025, width)]
        for index in range(18)
    ])
    report = {
        "arc_fraction_s": arc.tolist(),
        "path_enthalpy_eV_per_cell": h0.tolist(),
        "excess_enthalpy_meV_per_GaN": matrix.tolist(),
    }
    full_h = np.full(29, -100.0)
    s_mesh, q_mesh, surface = model_surface(report, full_h)
    _, _, transverse = transverse_surface(report, full_h)
    assert s_mesh.shape == q_mesh.shape == surface.shape == (201, 401)
    assert transverse.shape == surface.shape
    assert np.allclose(transverse[100], 0.0, atol=1e-9)
    for position, index in ((0, 0), (-1, 17)):
        assert surface[100, position] == pytest.approx((h0[index] - full_h[0]) * 500)
        assert surface[0, position] == pytest.approx(
            (h0[index] - full_h[0]) * 500 + matrix[index, 0]
        )
        assert surface[-1, position] == pytest.approx(
            (h0[index] - full_h[0]) * 500 + matrix[index, -1]
        )
        assert transverse[0, position] == pytest.approx(matrix[index, 0])
        assert transverse[-1, position] == pytest.approx(matrix[index, -1])


def test_frozen_dense_figure_has_ninety_traceable_coordinates() -> None:
    report_path = ROOT / "benchmarks/numerical_integrity/gan_600eV_atomic_tube_dense_20260928.json"
    trajectory = (
        ROOT / "paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927"
        / "gan_vasp_45p7_final_chain.traj"
    )
    report, frames = load(report_path, trajectory)
    assert len(frames) == 29
    assert report["maximum_prospective_s_holdout_error_meV_per_GaN"] < 0.17
    assert report["maximum_q_halfstep_error_meV_per_GaN"] < 0.16
    assert len({case["outcar_sha256"] for case in report["new_cases"]}) == 44
    csv_path = ROOT / "paper/VARNEB_CPC/figures/gan_600eV_atomic_dense_surface_source_data.csv"
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 90
    assert sum(row["source_kind"] == "path_trajectory" for row in rows) == 18
    assert sum(row["source_kind"] == "prior_static" for row in rows) == 28
    assert sum(row["source_kind"] == "new_static" for row in rows) == 44


def test_main_transverse_figure_tracks_162_dft_points_without_claiming_lower_mep() -> None:
    figures = ROOT / "paper/VARNEB_CPC/figures"
    stem = "gan_600eV_atomic_transverse_162_20260929_v2"
    qa = json.loads(
        (figures / f"{stem}_qa.json").read_text(encoding="utf-8")
    )
    manuscript = (ROOT / "paper/VARNEB_CPC/varneb_CPC.tex").read_text(encoding="utf-8")
    assert f"{{figures/{stem}.pdf}}" in manuscript
    assert r"{figures/gan_600eV_atomic_dense_surface.pdf}" not in manuscript
    assert qa["status"] == "GaN_600eV_central_18x9_frozen_transverse_figure"
    assert qa["n_audited_DFT_points"] == 162
    assert qa["n_measured_lower_than_centerline_by_0p05_meV"] == 2
    assert "not a whole-path" in qa["claim_limit"]
    assert qa["source_sha256"]["plotter"] == hashlib.sha256(
        (ROOT / "scripts/plot_gan_600eV_atomic_transverse_q9.py").read_bytes()
    ).hexdigest()
    assert qa["source_sha256"]["audit"] == hashlib.sha256(
        (ROOT / "benchmarks/numerical_integrity/gan_600eV_atomic_tube_q9_20260929.json").read_bytes()
    ).hexdigest()
    csv_path = figures / f"{stem}_source_data.csv"
    assert qa["source_sha256"]["source_data"] == hashlib.sha256(
        csv_path.read_bytes()
    ).hexdigest()

    report, full_arc, full_h = load_q9(
        ROOT / "benchmarks/numerical_integrity/gan_600eV_atomic_tube_q9_20260929.json",
        ROOT / "benchmarks/numerical_integrity/gan_600eV_atomic_tube_refinement_20260928.json",
        ROOT / "paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927/gan_vasp_45p7_final_chain.traj",
    )
    measured = np.asarray(report["excess_enthalpy_meV_per_GaN"], dtype=float)
    assert report["pressure_GPa"] == 45.7
    assert report["n_reused_measured_points"] == 90
    assert report["n_new_raw_audited_statics"] == 72
    assert measured.shape == (18, 9)
    assert np.allclose(measured[:, 4], 0.0, atol=1e-8)
    assert int(np.count_nonzero(measured < -0.05)) == 2
    assert measured.min() == pytest.approx(
        qa["minimum_measured_excess_meV_per_GaN"], abs=1e-9
    )
    assert report["prospective_gate_pass"] is True
    assert report["prospective_max_abs_error_meV_per_GaN"] < 0.01
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 162
    assert {row["kind"] for row in rows} == {"audited_DFT"}
    assert np.allclose(
        [float(row["excess_enthalpy_meV_per_GaN"]) for row in rows],
        measured.ravel(), atol=1e-9, rtol=0,
    )

    fig = draw_q9(report, full_arc, full_h)
    try:
        panel = fig.axes[0]
        sampled_grid = next(
            item for item in panel.collections
            if isinstance(item, PathCollection) and len(item.get_offsets()) == 162
        )
        assert sampled_grid.get_zorder() > 2  # above the filled contour
        lower_markers = next(
            item for item in panel.collections
            if isinstance(item, PathCollection) and len(item.get_offsets()) == 2
        )
        assert lower_markers.get_zorder() > max(
            line.get_zorder() for line in panel.lines
        )
    finally:
        plt.close(fig)
