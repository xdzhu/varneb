"""Keep the GaN figure stars tied to their audited, non-endpoint references."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from scripts.plot_gan_600eV_atomic_transverse_q9 import draw as draw_transverse
from scripts.plot_gan_600eV_atomic_transverse_q9 import load as load_transverse
from scripts.plot_gan_600eV_local_joint_cut17 import draw as draw_joint
from scripts.plot_gan_600eV_local_joint_cut17 import load_rows


ROOT = Path(__file__).resolve().parents[1]
NUMERIC = ROOT / "benchmarks" / "numerical_integrity"


def test_joint_cut_star_is_refined_image15_reference():
    pilot = json.loads((NUMERIC / "gan_600eV_ts_2d_pilot_20260928.json")
                       .read_text(encoding="utf-8"))
    newton = json.loads((NUMERIC / "gan_600eV_ts_newton_canary_20260928.json")
                        .read_text(encoding="utf-8"))
    assert pilot["source_sha256"]["center_OUTCAR"] == newton["OUTCAR_sha256"]
    rows, report = load_rows(
        NUMERIC / "gan_600eV_ts_2d_pilot_20260928.json",
        NUMERIC / "gan_600eV_ts_2d_5x5_refinement_20260928.json",
        NUMERIC / "gan_600eV_ts_2d_9x9_refinement_20260928.json",
        NUMERIC / "gan_600eV_ts_2d_17x17_refinement_20260930.json",
    )
    fig = draw_joint(rows, report)
    try:
        assert [label.get_text() for label in fig.legends[0].get_texts()] == [
            "One-step-refined image 15"
        ]
        np.testing.assert_allclose(fig.axes[0].collections[-1].get_offsets(), [[0, 0]])
    finally:
        plt.close(fig)


def test_transverse_stars_follow_audited_peak_image():
    report, arc, enthalpy = load_transverse(
        NUMERIC / "gan_600eV_atomic_tube_q9_20260929.json",
        NUMERIC / "gan_600eV_atomic_tube_refinement_20260928.json",
        ROOT / "paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927/gan_vasp_45p7_final_chain.traj",
    )
    assert int(np.argmax(enthalpy)) == 15
    fig = draw_transverse(report, arc, enthalpy)
    try:
        labels = [label.get_text() for legend in fig.legends
                  for label in legend.get_texts()]
        assert "Peak image 15" in labels
        peak_row = report["central_image_indices"].index(15)
        np.testing.assert_allclose(fig.axes[0].collections[-1].get_offsets(),
                                   [[report["arc_fraction_s"][peak_row], 0]])
        np.testing.assert_allclose(fig.axes[1].collections[-1].get_offsets(),
                                   [[arc[15], (enthalpy[15] - enthalpy[0]) * 500]])
    finally:
        plt.close(fig)

    wrong_peak = enthalpy.copy()
    wrong_peak[14] = wrong_peak[15] + 1.0
    with pytest.raises(ValueError, match="peaks at central image 15"):
        draw_transverse(report, arc, wrong_peak)
