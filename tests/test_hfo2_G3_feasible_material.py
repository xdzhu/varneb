"""Matched prospective source/cache contract; injected delegates are not DFT."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "benchmarks/hfo2_channels/20261008/G3_feasible_material_E170_20261015"
NATIVE = (
    ROOT
    / "benchmarks/hfo2_channels/20261008/G3_bound_diagnostic_E168_20261015/preserving_native_prefix_and_holes.json"
)
CACHE = CASE / "original_point13_cache.zip"
SOLVER = ROOT / "vcneb/feasible_relaxation.py"
spec = importlib.util.spec_from_file_location(
    "matched_feasible_case", CASE / "matched.py"
)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


@pytest.fixture(scope="module")
def sources():
    return m.load_sources(NATIVE, CACHE)


def test_actual_configuration_freezes_one_matched_training_arm_not_an_independent_point(
    sources,
):
    config = json.loads((CASE / "config.json").read_bytes())
    assert config == m.prepare(NATIVE, CACHE, SOLVER)
    m.validate(config, NATIVE, CACHE, SOLVER)
    point, records = sources
    assert point["id"] == 13 and point["point"] is None and len(records) == 12
    assert (
        config["original_requested_start_internal_A"]
        == point["starting_internal_coordinates_A"]
    )
    assert config["original_first_evaluated_full_A"] == records[0]["coordinates"]
    assert config["policy"]["preserving_training_remaining"] == 51
    assert config["policy"]["T_training_remaining"] == 52
    assert config["policy"]["independent_slots_untouched_per_face"] == 8
    assert config["policy"]["maximum_new_attempts"] == 12
    assert config["policy"]["original_E166_controls_have_priority"]
    assert (
        config["old_baseline_existing_SCFs"] == 12
        and config["old_baseline_new_SCFs"] == 0
    )
    assert config["released_radius_A"] ** 2 + sum(
        q * q for q in config["policy"]["q_A"]
    ) == pytest.approx(0.25, abs=2e-16)
    assert not config["physical_reserved_geometry_generated_or_read"]


@pytest.mark.parametrize(
    "fault",
    ["target", "radius", "trust", "budget", "holdout", "source", "start", "solver"],
)
def test_changed_physical_start_domain_target_budget_or_source_is_rejected(fault):
    c = deepcopy(m.prepare(NATIVE, CACHE, SOLVER))
    if fault == "target":
        c["policy"]["gradient_target_eV_A"] = 0.1
    elif fault == "radius":
        c["policy"]["original_total_internal_radius_A"] = 0.6
    elif fault == "trust":
        c["policy"]["initial_trust_radius_A"] = 0.1
    elif fault == "budget":
        c["policy"]["maximum_new_attempts"] = 24
    elif fault == "holdout":
        c["policy"]["independent_slots_untouched_per_face"] = 7
    elif fault == "source":
        c["original_six_input_sha256"]["INPUT"] = "0" * 64
    elif fault == "start":
        c["original_requested_start_internal_A"][2] += 0.01
    else:
        c["new_solver_text_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="contract differs"):
        m.validate(c, NATIVE, CACHE, SOLVER)


def test_all_original_records_are_exact_cache_hits_with_defensive_gradient_copies(
    sources,
):
    _, records = sources
    oracle = m.SharedExactCache(records)
    for r in records:
        e, g = oracle(np.array(r["coordinates"]))
        assert e == r["enthalpy_eV"]
        np.testing.assert_array_equal(g, r["gradient"])
        g[0] += 1
        np.testing.assert_array_equal(
            oracle(np.array(r["coordinates"]))[1], r["gradient"]
        )
    assert oracle.old_cache_hits == 24 and len(oracle.requests) == 24


def test_zero_DFT_probe_never_satisfies_a_missing_coordinate_from_a_fake_result(
    sources,
):
    _, records = sources
    oracle = m.SharedExactCache(records)
    x = np.array(records[0]["coordinates"])
    x[0] += 1e-9
    with pytest.raises(m.NewCoordinateRequired) as caught:
        oracle(x)
    np.testing.assert_array_equal(caught.value.coordinates, x)
    assert not oracle.requests[-1]["original_cache_hit"]


def test_new_delegate_must_share_original_contract_and_propagates_its_budget_error(
    sources,
):
    _, records = sources
    called = []

    def capped(x):
        called.append(x.copy())
        raise RuntimeError("explicit toy delegate cap")

    class Delegate:
        contract = {"contract_sha256": records[0]["contract_sha256"]}
        fixed_coordinate_covectors = np.eye(40)[:, -1:]

        def __call__(self, x):
            return capped(x)

    delegate = Delegate()
    oracle = m.SharedExactCache(records, delegate)
    x = np.array(records[0]["coordinates"])
    x[0] += 1e-9
    with pytest.raises(RuntimeError, match="toy delegate cap"):
        oracle(x)
    assert len(called) == 1
    bad = SimpleNamespace(
        contract={"contract_sha256": "different"},
        fixed_coordinate_covectors=np.eye(40)[:, -1:],
    )
    with pytest.raises(ValueError, match="physical contract"):
        m.SharedExactCache(records, bad)


def test_exact_control_and_duplicate_cache_rows_are_refused(sources):
    _, records = sources
    oracle = m.SharedExactCache(records)
    x = np.array(records[0]["coordinates"])
    x[-1] = 1e-10
    with pytest.raises(ValueError, match="prescribed control"):
        oracle(x)
    with pytest.raises(ValueError, match="unique original"):
        m.SharedExactCache(records + [records[0]])


def test_existing_config_is_refused_before_preparation(tmp_path, monkeypatch):
    p = tmp_path / "old.json"
    p.write_text("preserve")
    monkeypatch.setattr(
        "sys.argv",
        [
            "matched",
            "prepare",
            "--native",
            "unused",
            "--cache",
            "unused",
            "--solver",
            "unused",
            "--config",
            str(p),
        ],
    )
    monkeypatch.setattr(m, "prepare", lambda *a: pytest.fail("must refuse first"))
    with pytest.raises(FileExistsError):
        m.main()
    assert p.read_text() == "preserve"


def test_preparation_exposes_no_run_or_submission_mode(monkeypatch):
    monkeypatch.setattr("sys.argv", ["matched", "run"])
    with pytest.raises(SystemExit):
        m.main()


def test_frozen_configuration_bytes_and_actual_hf_probe_preserve_start_and_old_contract():
    assert (
        hashlib.sha256((CASE / "config.json").read_bytes()).hexdigest()
        == "f937a2e31ce30ec9b15852fa6a8e06cb73f4abbd4aa9b7f1f26bd29d561b4fd3"
    )
    raw = (CASE / "native_matched_probe.json").read_bytes()
    assert (
        hashlib.sha256(raw).hexdigest()
        == "2c22c66c58511c51d768271c17f0b6cadec7dfb6781af5106cfbd0822d5ce7c7"
    )
    r = json.loads(raw)
    assert r["status"] == "matched_old_start_reused_first_novel_proposal_no_DFT"
    assert (
        r["first_original_start_static_exact"]
        and r["old_cache_hits_before_first_novel"] == 7
    )
    assert (
        r["original_runtime_Python_files_rechecked"] == 380
        and r["original_native_SCFs_reparsed"] == 12
    )
    assert r["original_baseline_exact_cache_replay"] == dict(
        original_order_cache_only_replay_passed=True, new_DFT=0
    )
    assert r["first_novel_internal_norm_A"] <= 0.5 and r["first_novel_full_A"][-1] == 0
    assert r["inert_first_novel_minimum_distance_A"] >= 1.6
    assert all(
        r["inert_first_novel_input_sha256"][n] == h
        for n, h in m.faces.study.CONTRACT.items()
    )
    assert (
        r["new_DFT_or_scheduler_writes"] == 0
        and not r["original_worker_source_or_parameters_changed"]
    )
    assert not r["physical_reserved_geometry_generated_or_read"]
