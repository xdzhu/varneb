import csv
from pathlib import Path

import pytest

from scripts import plot_hfo2_G1_terminal_network as figure


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def data():
    return figure.build_data(ROOT)


def test_forty_actual_terminal_records_and_refinement_identity(data):
    assert len(data["rows"]) == 40 and data["new_DFT_calls"] == 0
    assert [c["snapshot_step"] for c in data["channels"]] == [6, 1, 69, 45]
    assert [len(c["rows"]) for c in data["channels"]] == [10, 12, 9, 9]
    assert [c["peak_energy_relative_PO_meV_fu"] for c in data["channels"]] == pytest.approx(
        [115.21016416145358, 82.67995730602706, 32.806023216835456, 392.8229051971357], abs=1e-9)
    assert all(c["ordinary_residual_passed"] for c in data["channels"])
    assert not data["full_G1_passed"] and not data["TS_certified"]
    assert abs(data["rows"][22+3]["preserving_true_tangent_eV_A"]) > .24


def test_terminal_export_source_alignment_and_no_clipped_energy(data, tmp_path):
    qa = figure.export(data, tmp_path/"figure")
    assert qa["row_alignment_verified"] and qa["column_alignment_verified"]
    assert qa["panel_label_weight"] == "normal" and qa["all_spines_visible"]
    assert qa["panel_titles_empty"] and qa["shared_legend_frame_alpha"] == .85
    with (tmp_path/"figure/source_data.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 40 and max(float(r["energy_relative_PO_meV_fu"]) for r in rows) > 392
    assert not qa["smooth_MEP_interpolation_used"]
    with pytest.raises(FileExistsError):
        figure.export(data, tmp_path/"figure")
