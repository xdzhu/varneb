"""The GaN contour uses only a validated bounded interpolant."""

from __future__ import annotations

import json

import numpy as np
import pytest

from scripts.plot_gan_600eV_atomic_dense_surface import load, model_surface


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
    assert s_mesh.shape == q_mesh.shape == surface.shape == (201, 401)
    for position, index in ((0, 0), (-1, 17)):
        assert surface[100, position] == pytest.approx((h0[index] - full_h[0]) * 500)
        assert surface[0, position] == pytest.approx(
            (h0[index] - full_h[0]) * 500 + matrix[index, 0]
        )
        assert surface[-1, position] == pytest.approx(
            (h0[index] - full_h[0]) * 500 + matrix[index, -1]
        )
