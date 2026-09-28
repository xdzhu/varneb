"""The GaN contour uses only a validated bounded interpolant."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from ase.units import GPa
from matplotlib.collections import PathCollection
import matplotlib.pyplot as plt

from scripts.plot_gan_600eV_atomic_dense_surface import load, model_surface
from scripts.plot_gan_600eV_atomic_transverse_landscape import draw, transverse_surface


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


def test_main_figure_is_transverse_and_does_not_imply_lower_dft_route() -> None:
    figures = ROOT / "paper/VARNEB_CPC/figures"
    qa = json.loads(
        (figures / "gan_600eV_atomic_transverse_landscape_qa.json").read_text(encoding="utf-8")
    )
    manuscript = (ROOT / "paper/VARNEB_CPC/varneb_CPC.tex").read_text(encoding="utf-8")
    assert r"{figures/gan_600eV_atomic_transverse_landscape.pdf}" in manuscript
    assert r"{figures/gan_600eV_atomic_dense_surface.pdf}" not in manuscript
    assert qa["status"] == "GaN_600eV_central_atomic_transverse_landscape"
    assert qa["source_sha256"]["plotter"] == hashlib.sha256(
        (ROOT / "scripts/plot_gan_600eV_atomic_transverse_landscape.py").read_bytes()
    ).hexdigest()

    report, frames = load(
        ROOT / "benchmarks/numerical_integrity/gan_600eV_atomic_tube_dense_20260928.json",
        ROOT / "paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927/gan_vasp_45p7_final_chain.traj",
    )
    measured = np.asarray(report["excess_enthalpy_meV_per_GaN"], dtype=float)
    assert np.allclose(measured[:, 2], 0.0)
    assert measured[:, [0, 1, 3, 4]].min() > 0.0
    full_h = np.asarray([
        frame.get_potential_energy() + 45.7 * GPa * frame.get_volume()
        for frame in frames
    ])
    _, _, interpolated = transverse_surface(report, full_h)
    assert interpolated.min() == pytest.approx(
        qa["transverse_energy_min_max_meV_per_GaN"][0], abs=1e-9
    )
    assert -0.3 < interpolated.min() < -0.2

    fig = draw(report, frames)
    try:
        panel = fig.axes[0]
        sampled_grid = next(
            item for item in panel.collections
            if isinstance(item, PathCollection) and len(item.get_offsets()) == 90
        )
        assert sampled_grid.get_zorder() > max(line.get_zorder() for line in panel.lines)
    finally:
        plt.close(fig)
