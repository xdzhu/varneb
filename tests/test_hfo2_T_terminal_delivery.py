"""Actual E079 data and native admission, not synthetic scheduler evidence."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts.analyze_hfo2_clamped_reference import terminal_images
from scripts.audit_hfo2_static_replica import sha256
from scripts.audit_hfo2_clamped_terminal import completed_step, fresh_interior_call_count
from scripts.export_hfo2_clamped_observation import PILOT_SCRIPT_SHA256

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "benchmarks/hfo2_channels/20261008"
CASE = BASE / "clamped_T_terminal_E079_20261011"
AUDIT = "c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e"


@pytest.fixture(scope="module")
def actual_terminal():
    return terminal_images(CASE, AUDIT)


def test_actual_full_native_audit_and_fresh_starter_cost(actual_terminal):
    _, report, _, audit = actual_terminal
    assert audit["job_id"] == "28692775"
    assert audit["scheduler_row"] == "28692775|COMPLETED|0:0|32|hfacnormal01|iai806"
    assert completed_step(CASE / "vcneb.opt.log") == audit["terminal_step"] == 10
    assert audit["initial_fresh_interior_SCFs"] == 7
    assert audit["fresh_interior_SCFs"] == audit["expected_fresh_interior_SCFs"] == 77
    assert fresh_interior_call_count(10, json.loads(
        (CASE / "observations/step_0010/runtime_preflight.json").read_text()),
        PILOT_SCRIPT_SHA256) == (7, 77)
    assert len(audit["runtime_code_sha256"]) == 338
    expected_paths = {(f"image_{i:04d}", f"scf_{j:06d}")
                      for i in range(1, 8) for j in range(11)}
    calls = audit["calls"]
    assert len(calls) == 77
    assert {(Path(r["path"]).parts[-3], Path(r["path"]).parts[-2])
            for r in calls} == expected_paths
    assert all({k: r["input_sha256"][k] for k in CONTRACT} == CONTRACT for r in calls)
    assert all(r["elapsed_seconds"] > 0 for r in calls)
    wall_seconds = sum(r["elapsed_seconds"] for r in calls)
    assert wall_seconds == pytest.approx(audit["fresh_SCF_wall_seconds_sum"], abs=1e-10)
    # Seconds * 32 allocated ranks / 3600 seconds per hour is a transport-cost
    # proxy, not total campaign core-hours or an optimizer speedup benchmark.
    assert wall_seconds * 32 / 3600 == pytest.approx(92.21180905897585, abs=1e-10)
    assert audit["new_DFT_calls"] == audit["scheduler_mutations"] == 0
    assert audit["continuous_TS_or_barrier_error_bound_certified"] is False
    assert report["n_total_images"] == 9 and report["n_active_images"] == 7


def test_actual_nine_frame_energy_units_boundary_and_residual(actual_terminal):
    images, report, residual, _ = actual_terminal
    assert residual["ordinary_residual_pass"] is True
    assert residual["fmax_eV_A"] == pytest.approx(.09358976572425506, abs=1e-12)
    assert report["formula_units"] == 4 and report["pressure_GPa"] == 0
    assert report["climb"] is False and report["k_eV_A2"] == .2
    assert report["mechanical_boundary"]["fixed_cell_rows"] == [0, 1]
    assert report["mechanical_boundary"]["allow_tilt"] is True
    substrate_rows = np.asarray(report["mechanical_boundary"]["reference_cell_A"])[:2]
    assert all(np.allclose(image.cell.array[:2], substrate_rows, atol=1e-12, rtol=0)
               for image in images)
    energy_eV_cell = np.array([image.get_potential_energy() for image in images])
    # Explicit cell-to-formula-unit normalization, followed by eV -> meV.
    energy_meV_fu = (energy_eV_cell - energy_eV_cell[0]) / report["formula_units"] * 1000
    assert np.argmax(energy_meV_fu) == 4
    assert max(energy_meV_fu) == pytest.approx(84.41662388122495, abs=1e-9)
    assert max(energy_meV_fu) - energy_meV_fu[-1] == pytest.approx(72.05664876710216, abs=1e-9)
    assert energy_meV_fu[-1] == pytest.approx(12.35997511412279, abs=1e-9)


def test_comparison_reuses_same_initial_well_and_mechanical_ensemble(actual_terminal):
    images, report, _, _ = actual_terminal
    m_images, m_report, _, _ = terminal_images(BASE / "clamped_M_terminal_E068_20261011",
        "56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c")
    assert same_ordered_geometry(images[0], m_images[0])
    assert images[0].get_potential_energy() == m_images[0].get_potential_energy()
    for key in ("cell_scale_A", "pressure_GPa", "climb", "k_eV_A2", "formula_units"):
        assert report[key] == m_report[key]
    for key in ("reference_cell_A", "fixed_cell_rows", "allow_tilt", "raw_chain_policy"):
        assert report["mechanical_boundary"][key] == m_report["mechanical_boundary"][key]
    assert report["raw_image_evaluations"][0]["raw_log_sha256"] == m_report["raw_image_evaluations"][0]["raw_log_sha256"]


def test_actual_release_journal_and_separate_native_startup(actual_terminal):
    audit = actual_terminal[-1]
    helper_path = BASE / "clamped_T_terminal_E075_20261011/release_existing_reversing_once.py"
    spec = importlib.util.spec_from_file_location("actual_E079_gate", helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    helper.verified_audit(audit)  # Only a pure check, never call release_once.
    admission = json.loads((CASE / "admission_receipt.json").read_text())
    rows = [json.loads(line) for line in (CASE / "admission_journal.jsonl").read_text().splitlines()]
    assert [r["event"] for r in rows] == ["admission_verified", "release_intent", "release_accepted"]
    assert rows[0]["parent_audit_sha256"] == AUDIT
    assert rows[0]["active_jobs"] == ["28723655"]
    assert rows[0]["held_row"]["Dependency"] == "(null)"
    assert rows[0]["held_row"]["JobId"] == "28692776"
    pending = rows[0]["pending_seed_preflight"]
    assert pending["geometry_images_checked"] == 9 and pending["fresh_interiors"] == 7
    assert len(pending["cached_endpoints"]) == 2 and pending["new_DFT_calls"] == 0
    assert all(r["original_six_physical_bytes_native_DSIZE32_full_EFS_checked"] for r in pending["cached_endpoints"])
    assert admission["executed_helper_sha256"] == sha256(helper_path)
    assert admission["journal_sha256"] == sha256(CASE / "admission_journal.jsonl")
    assert admission["release_writes"] == 1 and admission["max_study_simultaneous_allocations"] == 2
    assert admission["new_submissions"] == admission["dependency_writes"] == admission["new_independent_chains"] == 0
    assert admission["startup_claimed"] is False
    startup = json.loads((CASE / "startup_observation.json").read_text())
    assert startup["target_job"] == "28692776" and startup["state"] == "RUNNING"
    assert startup["CPUs"] == 32 and startup["partition"] == "hfacnormal01"
    assert "DSIZE = 32" in startup["native_first_fresh_SCF_DSIZE_line"]
    assert startup["current_running_study_allocations"] == 2
    assert startup["four_plus1percent_states"] == "PENDING/JobHeldUser"
    assert startup["stderr_bytes_at_capture"] == 0
    assert startup["first_complete_step_zero_or_fresh_SCF_completion_claimed"] is False


def test_actual_rebound_was_followed_by_convergence_not_stop(actual_terminal):
    rows = [line.split() for line in (CASE / "vcneb.opt.log").read_text().splitlines()
            if line.startswith("CheckedFIRE:")]
    forces_eV_A = [float(row[-1]) for row in rows]
    assert len(rows) == 11 and forces_eV_A[5] == .103118
    assert forces_eV_A[6] > forces_eV_A[5] and forces_eV_A[7] > forces_eV_A[6]
    assert forces_eV_A[-1] < .10
