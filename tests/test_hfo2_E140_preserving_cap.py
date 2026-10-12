"""Native successful child40 cap: valid seed, not ordinary convergence."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
CASE = REPO / "benchmarks/hfo2_channels/20261008/clamped_preserving_cap_E140_20261012"
ROOT = CASE / "completed_HF"
AUDIT_SHA = "df34573384f76650e123691dbb508a52d7d8ee529e6e50dca0fce5e4e6835c02"
SEED_SHA = "1a222cc63b89344816f0829bbf079793fa06306d3d4d13f5cc3c2b401f5fbf3f"
FMAX = 0.11414663193850914


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def receipt():
    assert sha(ROOT / "audit_receipt.json") == AUDIT_SHA
    return json.loads((ROOT / "audit_receipt.json").read_text())


def test_real_cap_has_all_call_receipt_but_is_not_converged():
    from examples.hfo2_fixed_input_factory import CONTRACT
    from scripts.audit_hfo2_clamped_terminal import validate_terminal

    record = receipt()
    assert record["job_id"] == "28737465" and record["terminal_step"] == 40
    assert not record["ordinary_converged"]
    summary = json.loads((ROOT / "vcneb_summary.json").read_text())
    assert sha(ROOT / "vcneb_summary.json") == record["terminal_summary_sha256"]
    validate_terminal("28737465", record["scheduler_row"], summary, False, outcome="step-cap")
    with pytest.raises(ValueError):
        validate_terminal("28737465", record["scheduler_row"], summary, False)
    assert record["fresh_interior_SCFs"] == record["expected_fresh_interior_SCFs"] == 280
    assert record["initial_fresh_interior_SCFs"] == 0
    assert len(record["calls"]) == 280 and len(record["runtime_code_sha256"]) == 338
    assert {(Path(row["path"]).parts[-3], Path(row["path"]).parts[-2])
            for row in record["calls"]} == {
                (f"image_{i:04d}", f"scf_{j:06d}") for i in range(1, 8) for j in range(40)}
    assert all({name: row["input_sha256"][name] for name in CONTRACT} == CONTRACT
               for row in record["calls"])
    assert record["all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked"]
    assert record["current_frame_exact_caches_checked"] == 9
    assert record["new_DFT_calls"] == record["scheduler_mutations"] == record["new_independent_chains"] == 0
    assert not record["physical_inputs_or_running_source_changed"]
    assert not record["holdout_generated_or_read"] and not record["FIRE_state_restored"]


def test_nine_retained_raw_records_and_residual_are_replayed():
    from ase.io import read
    from examples.hfo2_fixed_input_factory import read_fixed_hfo2_stru, same_ordered_geometry
    from scripts.analyze_hfo2_clamped_residual import analyze
    from scripts.audit_hfo2_static_replica import audited_results

    record = receipt()
    observation = ROOT / "observations/step_0040"
    report = json.loads((observation / "observation.json").read_text())
    assert report == record["terminal_observation"]
    result = analyze(observation)
    assert result["fmax_eV_A"] == pytest.approx(FMAX, abs=1e-12)
    assert not result["ordinary_residual_pass"]
    images = read(observation / "evaluated_chain.traj", index=":")
    assert len(images) == 9
    for image, row in zip(images, report["raw_image_evaluations"]):
        raw = observation / "raw" / f"image_{row['image_index']:04d}"
        assert sha(raw / "OUT.ABACUS/running_scf.log") == row["raw_log_sha256"]
        assert sha(raw / "source_audit.json") == row["audit_sha256"]
        assert all(sha(raw / name) == row["input_sha256"][name] for name in ("INPUT", "KPT", "STRU"))
        assert same_ordered_geometry(image, read_fixed_hfo2_stru(raw / "STRU"))
        parsed = audited_results(raw)
        for name, expected in (("energy", row["energy_eV_cell"]),
                               ("forces", row["forces_eV_A"]),
                               ("stress", row["stress_ASE_voigt_eV_A3"])):
            np.testing.assert_allclose(parsed[name], expected, atol=1e-12, rtol=0)


def test_latest_seed_is_exact_geometry_and_cache_not_optimizer_history():
    from ase.io import read
    from examples.hfo2_fixed_input_factory import same_ordered_geometry

    record = receipt()
    seed = ROOT / "seed"
    assert sha(seed / "manifest.json") == record["seed_manifest_sha256"] == SEED_SHA
    manifest = json.loads((seed / "manifest.json").read_text())
    assert manifest["source_job_id"] == "28737465" and manifest["source_step"] == 40
    assert manifest["current_frame_exact_caches"] == 9 and not manifest["optimizer_state_restored"]
    assert manifest["new_DFT_calls"] == manifest["new_independent_chains"] == 0
    assert not manifest["physical_inputs_changed"] and not manifest["holdout_generated_or_read"]
    assert all(sha(seed / name) == pin for name, pin in manifest["files_sha256"].items())
    observation = ROOT / "observations/step_0040"
    assert sha(observation / "observation.json") == manifest["source_observation_sha256"]
    old = read(observation / "evaluated_chain.traj", index=":")
    new = read(seed / "seed.traj", index=":")
    assert len(old) == len(new) == 9
    for a, b in zip(old, new):
        assert same_ordered_geometry(a, b)
        for key in ("energy", "forces", "stress"):
            np.testing.assert_allclose(a.calc.results[key], b.calc.results[key], atol=1e-12, rtol=0)


def test_cost_and_unchanged_executed_driver_provenance():
    record = receipt()
    assert record["fresh_SCF_wall_seconds_sum"] == pytest.approx(
        sum(row["elapsed_seconds"] for row in record["calls"]), abs=1e-8)
    assert record["fresh_SCF_core_hour_proxy_32"] == pytest.approx(257.81047677434566, abs=1e-9)
    driver = json.loads((ROOT / "driver_receipt.json").read_text())
    assert driver["job_id"] == "28737465" and driver["original_parent"] == "28709790"
    assert driver["terminal_audit_receipt_sha256"] == AUDIT_SHA
    helper = REPO / "benchmarks/hfo2_channels/20261008/clamped_flip_continue_E103_20261012/continue_flips_HF.py"
    assert driver["executed_driver_sha256"] == sha(helper)
    assert driver["new_DFT_calls"] == driver["scheduler_writes"] == driver["new_independent_chains"] == 0
