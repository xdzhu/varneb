import importlib
import json
from pathlib import Path

import pytest

from scripts.audit_hfo2_static_replica import sha256
from vcneb.polarization import sampled_band_gap

CASE = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/reversing_peak_sampling_20261009"
DRIVER = importlib.import_module("benchmarks.hfo2_channels.20261008.reversing_peak_sampling_20261009.check_terminal_network")


@pytest.fixture(scope="module")
def terminal():
    return DRIVER.audit(CASE / "network_specification.json")


def test_all_four_actual_terminals_common_well_not_G1_certificate(terminal):
    assert len(terminal["channels"]) == 4 and terminal["existing_image_records"] == 40
    assert all(c["ordinary_residual_passed"] for c in terminal["channels"])
    assert terminal["common_PO_energy_eV_cell"] == pytest.approx(-9783.249675811956, rel=0, abs=1e-10)
    assert all(c["rows"][0]["energy_relative_PO_meV_fu"] == pytest.approx(0., rel=0, abs=1e-10)
               for c in terminal["channels"])
    assert terminal["new_DFT_calls"] == 0 and not terminal["full_G1_passed"] and not terminal["TS_certified"]


def test_preserving_nonpolar_basin_is_not_hidden(terminal):
    channel = next(c for c in terminal["channels"] if "preserving" in c["name"])
    center = channel["rows"][4]
    assert center["energy_relative_PO_meV_fu"] == pytest.approx(-14.838200784197397, rel=0, abs=1e-8)
    assert [s["symbol"] for s in center["structure_audit"]["symmetry_sweep"]] == ["Pbcn"]*3
    assert all(not r["structure_audit"]["standardized_wrapped_reordered_or_relaxed"] for r in channel["rows"])


def test_reversing_center_is_not_misidentified_from_zero_T_projection(terminal):
    channel = next(c for c in terminal["channels"] if "reversing" in c["name"])
    peak = channel["rows"][4]
    assert peak["energy_relative_PO_meV_fu"] == pytest.approx(392.8229051971357, rel=0, abs=1e-8)
    assert [s["symbol"] for s in peak["structure_audit"]["symmetry_sweep"]] == ["Pa-3"]*3
    assert peak["triplet_captured_squared_norm_fraction"] < 1e-20


def test_all40_original_band_files_independently_replay():
    root = CASE / "raw_gap_audit_HF"
    report = json.loads((root / "report.json").read_text())
    assert report["records"] == 40 and report["new_DFT_calls"] == 0 and not report["full_G1_passed"]
    assert len(report["channels"]) == 4
    for channel in report["channels"]:
        gaps = []
        for row in channel["rows"]:
            p = root / channel["name"] / f"istate_{row['source_image_index']:02d}.info"
            assert sha256(p) == row["istate_sha256"]
            gap = sampled_band_gap(p.read_text(), 48)
            assert gap == row["gap"]
            gaps.append(gap["sampled_indirect_gap_eV"])
        assert min(gaps) == channel["minimum_sampled_gap_eV"] and min(gaps) > 4.


def test_raw_gap_records_bound_to_current_observations(terminal):
    report = json.loads((CASE / "raw_gap_audit_HF/report.json").read_text())
    for channel, frozen in zip(report["channels"], terminal["channels"]):
        assert channel["name"] == frozen["name"]
        assert channel["observation_sha256"] == frozen["source_observation_sha256"]
        assert len(channel["rows"]) == len(frozen["rows"])
        by_index = {row["source_image_index"]: row for row in channel["rows"]}
        for row in frozen["rows"]:
            assert by_index[row["source_image_index"]]["raw_log_sha256"] == row["source_log_sha256"]
