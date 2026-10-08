import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from scripts.analyze_hfo2_channel_network import analyze, read_evaluated_observation


ROOT = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/channel_network"


def test_actual_four_channels_share_PO_initial_and_replay_without_DFT(tmp_path):
    expected = json.loads((ROOT / "analysis.json").read_text())
    actual = analyze(ROOT / "network_specification.json", tmp_path / "analysis.json")
    assert actual["existing_image_evaluations"] == 37 and actual["new_DFT_calls"] == 0
    assert not actual["physical_parameters_changed"] and not actual["missing_required_channels"]
    assert not actual["ready_for_bounded_discrete_comparison"]
    assert actual["H1_selectivity_conclusion"].startswith("not_evaluated")
    assert not actual["distinct_switching_MEPs_certified"]
    assert all(v is None for v in actual["minimum_barrier_intervals_eV_fu"].values())
    assert actual["selectivity_interval_eV_fu"] is None
    for candidate, reference in zip(actual["channels"], expected["channels"]):
        assert candidate["name"] == reference["name"]
        assert candidate["source_NEB_fmax_eV_A"] == pytest.approx(reference["source_NEB_fmax_eV_A"], abs=1e-10)
        assert candidate["barrier_from_common_initial_eV_fu"] == pytest.approx(reference["barrier_from_common_initial_eV_fu"], abs=1e-10)
        assert not candidate["NEB_residual_passed"]
    t = actual["channels"][0]
    assert t["source_direction_used"] == "reverse"
    assert t["barrier_from_common_initial_eV_fu"] == pytest.approx(.11627697978838114, abs=1e-10)
    assert t["source_forward_barrier_eV_fu"] == pytest.approx(.0349557870413264, abs=1e-10)
    modes = actual["reference_mode_observations_in_source_direction"]
    preserving, reversing = modes[-2:]
    assert preserving["residual_dominant_block"] == reversing["residual_dominant_block"] == "cell"
    assert .54 < preserving["images"][4]["parent_pattern_captured_squared_norm_fraction"] < .55
    assert reversing["images"][4]["parent_pattern_captured_squared_norm_fraction"] < 1e-20
    # Full-basis reconstruction is an identity, not a harmonic/TS claim.
    assert all(np.isfinite(im["T_Gamma_full_basis_reconstruction_residual_sqrt_amu_A"])
               for mode in modes for im in mode["images"])
    with pytest.raises(FileExistsError):
        analyze(ROOT / "network_specification.json", tmp_path / "analysis.json")


@pytest.mark.parametrize("fault", ["raw_energy", "forces", "stress", "input", "POSCAR", "trajectory", "optimizer_row"])
def test_changed_frozen_evidence_is_not_a_network_result(tmp_path, fault):
    folder = tmp_path / "observation"
    shutil.copytree(ROOT / "PO_flip_preserving_step1", folder)
    report_path = folder / "observation.json"
    report = json.loads(report_path.read_text())
    if fault in ("POSCAR", "trajectory"):
        with (folder / ("POSCAR_03" if fault == "POSCAR" else "evaluated_chain.traj")).open("ab") as f:
            f.write(b"tampered")
    else:
        point = report["raw_image_evaluations"][3]
        if fault == "raw_energy": point["energy_eV_cell"] += .001
        if fault == "forces": point["forces_eV_A"][0][0] += .001
        if fault == "stress": point["stress_ASE_voigt_eV_A3"][0] += .001
        if fault == "input": point["input_sha256"]["INPUT"] = "a" * 64
        if fault == "optimizer_row": report["optimizer_log_row"]["fmax_eV_A"] += .001
        report_path.write_text(json.dumps(report))
    with pytest.raises(ValueError): read_evaluated_observation(folder)


def test_terminal_switch_health_observation_does_not_promote_a_barrier_ranking(tmp_path):
    actual = analyze(ROOT / "terminal_switch_step10_specification.json", tmp_path / "latest.json")
    assert actual["existing_image_evaluations"] == 37 and actual["new_DFT_calls"] == 0
    assert not actual["ready_for_bounded_discrete_comparison"]
    assert actual["selectivity_interval_eV_fu"] is None
    assert all(interval is None for interval in actual["minimum_barrier_intervals_eV_fu"].values())
    assert [source["snapshot_step"] for source in actual["sources"]] == [10, 10, 10, 10]
    assert [source["source_job_id"] for source in actual["sources"][-2:]] == ["28288045", "28288063"]
    preserving, reversing = actual["channels"][-2:]
    assert preserving["barrier_from_common_initial_eV_fu"] == pytest.approx(.09720146879044478, abs=1e-10)
    assert reversing["barrier_from_common_initial_eV_fu"] == pytest.approx(.4107436409649381, abs=1e-10)
    modes = actual["reference_mode_observations_in_source_direction"][-2:]
    # The largest residual can move away from the energy maximum. Both energy
    # maxima are cell-dominated; the reversing chain's largest residual is now
    # atomic. Do not force one mechanism label onto every image/iteration.
    assert modes[0]["residual_dominant_block"] == "cell"
    assert modes[1]["residual_dominant_block"] == "atomic"
    assert all(mode["images"][4]["NEB_dominant_block"] == "cell" for mode in modes)
    assert .57 < modes[0]["images"][4]["parent_pattern_captured_squared_norm_fraction"] < .58
    assert modes[1]["images"][4]["parent_pattern_captured_squared_norm_fraction"] < 1e-20
