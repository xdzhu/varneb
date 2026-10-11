"""Real frozen scalar-data/figure checks, not another native SCF audit."""
import copy
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from scripts.plot_hfo2_clamped_M import (
    AUDIT, MODE_REPORT, OBSERVATION, build_data, canonical_sha, export, make_figure,
)

ROOT = Path(__file__).resolve().parents[1]


def test_real_profile_projections_and_green_strains():
    data = build_data(ROOT)
    assert len(data["rows"]) == 9 and data["highest_image_index"] == 3
    assert data["new_DFT_calls"] == 0 and not data["TS_certified"]
    peak, final = data["rows"][3], data["rows"][-1]
    assert peak["energy_relative_PO_meV_fu"] == pytest.approx(98.45812195362669)
    assert final["energy_relative_PO_meV_fu"] == pytest.approx(-92.54049439232404)
    assert peak["T_geometric_triplet_rank3_squared_fraction"] == pytest.approx(.5792373469688469)
    assert final["Green_E_01_original_T"] == pytest.approx(-.0752091109578331)
    assert final["Green_E_11_original_T"] == pytest.approx(.015509423995985316)
    assert all(r["Green_E_00_original_T"] == 0 for r in data["rows"])
    assert not data["native_data_reaudit_claimed"]


def test_checkout_newlines_do_not_change_mode_content_pin(tmp_path):
    for relative in (AUDIT, MODE_REPORT, OBSERVATION):
        dest = tmp_path / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, dest)
    path = tmp_path / MODE_REPORT
    text = path.read_text()
    path.write_bytes(text.replace("\n", "\r\n").encode())
    assert build_data(tmp_path)["rows"] == build_data(ROOT)["rows"]
    value = json.loads(text)
    other = copy.deepcopy(value)
    other["T_triplet_rank3_fraction"][0] += .01
    assert canonical_sha(value) != canonical_sha(other)


@pytest.mark.parametrize("relative", [MODE_REPORT, OBSERVATION, AUDIT])
def test_changed_frozen_source_is_rejected(tmp_path, relative):
    for src in (AUDIT, MODE_REPORT, OBSERVATION):
        dest = tmp_path / src
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / src, dest)
    value = json.loads((tmp_path / relative).read_text())
    value["source_job_id"] = "changed"
    (tmp_path / relative).write_text(json.dumps(value))
    with pytest.raises(ValueError, match="frozen material sources"):
        build_data(tmp_path)


def test_aligned_axes_and_requested_plot_contract():
    import matplotlib.pyplot as plt
    fig = make_figure(build_data(ROOT))
    try:
        assert len(fig.axes) == 3
        np.testing.assert_allclose(fig.get_size_inches()*25.4, [183, 207])
        positions = [ax.get_position().bounds for ax in fig.axes]
        for ax, pos, letter in zip(fig.axes, positions, ("(a)", "(b)", "(c)")):
            np.testing.assert_allclose([pos[0], pos[2]], [positions[0][0], positions[0][2]])
            assert not ax.get_title() and all(s.get_visible() for s in ax.spines.values())
            assert ax.texts[0].get_text() == letter and ax.texts[0].get_weight() == "normal"
            assert ax.get_legend().get_frame().get_alpha() == .85
        assert len(fig.axes[1].lines) == 7  # 2 T + all 4 Cmma + sampled-peak guide
    finally:
        plt.close(fig)


def test_export_is_fresh_only(tmp_path):
    with pytest.raises(FileExistsError, match="fresh dated"):
        export(build_data(ROOT), tmp_path)
