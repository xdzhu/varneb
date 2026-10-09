import importlib
import json
from pathlib import Path
import shutil

import pytest

from scripts import plot_hfo2_path_Berry as figure


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/figure.CASE
checker = importlib.import_module("benchmarks.hfo2_channels.20261008.switching_path_polarization_20261009.check_export")


def test_actual_56_raw_logs_band_tables_replay_identically_to_HF():
    result = checker.audit(ROOT, CASE/"completed_HF")
    assert result == json.loads((CASE/"offline_replay_HF.json").read_text())
    assert result["raw_SCF_logs"] == 14 and result["raw_NSCF_logs"] == 42
    assert result["native_integer_reported_DFT_seconds"] == 3155
    for c in result["channels"]:
        assert c["maximum_224_to_228_change_C_m2"] < .00084
        assert c["sampled_branch"]["delta_reduced"] == pytest.approx(1.6236913046100134)
        assert not c["sampled_branch"]["continuous_path_certified"]


@pytest.mark.parametrize("path", ["scf/OUT.ABACUS/running_scf.log", "nscf_222/OUT.ABACUS/BANDS_1.dat", "scf/INPUT"])
def test_any_raw_or_nongap_band_byte_mutation_is_refused(tmp_path, path):
    target = tmp_path/"export"
    shutil.copytree(CASE/"completed_HF", target)
    file = target/"calculations/preserving_01"/path
    file.write_bytes(file.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="inventory"):
        checker.audit(ROOT, target)


def test_figure_eighteen_rows_keep_gauge_quantum_and_measured_limit(tmp_path):
    data = figure.build_data(ROOT)
    assert len(data["rows"]) == 18 and not data["spontaneous_P_selected"]
    assert data["native_reduced_period"] == 2 and data["initial_gauge_integer"] == 0
    assert max(r["224_to_228_change_mC_m2"] for r in data["rows"]) < .84
    qa = figure.export(data, tmp_path/"figure")
    assert qa["row_alignment_verified"] and qa["all_spines_visible"]
    assert qa["panel_titles_empty"] and qa["panel_label_weight"] == "normal"
    assert not qa["continuous_path_full_G1_or_TS_certified"]
    with pytest.raises(FileExistsError):
        figure.export(data, tmp_path/"figure")
