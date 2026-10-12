"""Finite evaluator attempts, never a count of raw backend SCFs."""
import json
from pathlib import Path
import runpy

import numpy as np
import pytest
from ase.calculators.singlepoint import SinglePointCalculator

from vcneb.mode_evaluator import CalculatorModeEvaluator
from vcneb import ModePlane, sample_conditional_mode_surface


HELPERS = runpy.run_path(str(Path(__file__).with_name("test_mode_evaluator.py")))


def evaluator(tmp_path, limit, calls, *, resume=False, factory=None, minimum=0.):
    chart = HELPERS["_chart"]()
    return CalculatorModeEvaluator(
        chart, HELPERS["_factory"](calls) if factory is None else factory,
        workdir=tmp_path/"run", calculator_id="analytic:test:v1",
        max_evaluation_attempts=limit, resume=resume, minimum_distance_A=minimum,
    )


def point(e, displacement=0.):
    x = np.zeros(e.chart.coordinate_count)
    x[0] = displacement
    return x


def test_completed_and_cached_trials_obey_cumulative_attempt_budget(tmp_path):
    calls = []
    e = evaluator(tmp_path, 2, calls)
    first = e(point(e))
    e(point(e, .01))
    for _ in range(3):
        cached = e(point(e))
        assert cached[0] == first[0]
        np.testing.assert_array_equal(cached[1], first[1])
    assert calls == [0, 1]
    assert e.n_new_evaluations == e.n_evaluation_attempts == 2
    with pytest.raises(RuntimeError, match="attempt budget exhausted"):
        e(point(e, .02))
    assert calls == [0, 1]
    assert len(list(e.workdir.glob("attempt-*.json"))) == 2
    assert len(list(e.workdir.glob("eval-*-*"))) == 2


def test_resume_keeps_budget_consumption_and_cached_results(tmp_path):
    calls = []
    first = evaluator(tmp_path, 1, calls)
    expected = first(point(first))
    resumed = evaluator(tmp_path, 1, calls, resume=True)
    assert resumed.n_new_evaluations == 0 and resumed.n_evaluation_attempts == 1
    actual = resumed(point(resumed))
    assert actual[0] == expected[0]
    np.testing.assert_array_equal(actual[1], expected[1])
    with pytest.raises(RuntimeError, match="budget exhausted"):
        resumed(point(resumed, .01))
    assert calls == [0]


def test_factory_failure_consumes_attempt_after_resume(tmp_path):
    calls = []

    def fail(index, atoms, directory):
        calls.append(index)
        assert (directory.parent/f"attempt-{index:06d}.json").is_file()
        raise RuntimeError("injected calculator failure")

    e = evaluator(tmp_path, 1, calls, factory=fail)
    with pytest.raises(RuntimeError, match="injected"):
        e(point(e))
    assert e.n_new_evaluations == 0 and e.n_evaluation_attempts == 1
    resumed = evaluator(tmp_path, 1, calls, resume=True)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        resumed(point(resumed))
    assert calls == [0]
    assert len(list(e.workdir.glob("eval-*/failure.json"))) == 1


def test_invalid_geometry_consumes_no_attempt(tmp_path):
    calls = []
    e = evaluator(tmp_path, 1, calls, minimum=10.)
    with pytest.raises(ValueError, match="calculator was not called"):
        e(point(e))
    assert e.n_evaluation_attempts == 0 and calls == []
    assert not list(e.workdir.glob("attempt-*.json"))


def test_zero_budget_rejects_first_fresh_attempt(tmp_path):
    calls = []
    e = evaluator(tmp_path, 0, calls)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        e(point(e))
    assert e.n_evaluation_attempts == 0 and calls == []
    assert not list(e.workdir.glob("eval-*-*"))


@pytest.mark.parametrize("bad", [True, np.bool_(False), -1, 1.5, 1., "2", np.nan])
def test_invalid_budget_rejected_before_namespace_creation(tmp_path, bad):
    with pytest.raises(ValueError, match="nonnegative integer"):
        evaluator(tmp_path, bad, [])
    assert not (tmp_path/"run").exists()


@pytest.mark.parametrize("changed", [None, 0, 3])
def test_resume_cannot_remove_or_change_budget(tmp_path, changed):
    e = evaluator(tmp_path, 2, [])
    e(point(e))
    with pytest.raises(ValueError, match="budget changed"):
        evaluator(tmp_path, changed, [], resume=True)


def test_budget_is_separate_from_unchanged_legacy_physical_contract(tmp_path):
    plain = evaluator(tmp_path/"plain", None, [])
    capped = evaluator(tmp_path/"capped", 2, [])
    assert plain.contract == capped.contract
    assert plain.contract["contract_sha256"] == "1308f3467e845dcfdb138871212baca9e0a6a05fcc78c281abf9b81e4212226f"
    plain(point(plain))
    assert plain.n_evaluation_attempts == 1
    assert not list(plain.workdir.glob("attempt-*.json"))
    assert not (plain.workdir/"evaluation-budget.json").exists()
    with pytest.raises(ValueError, match="budget changed or missing"):
        evaluator(tmp_path/"plain", 2, [], resume=True)


def test_reservation_without_directory_is_conservatively_counted(tmp_path, monkeypatch):
    calls = []
    e = evaluator(tmp_path, 1, calls)
    original_mkdir = Path.mkdir

    def fail(path, *args, **kwargs):
        if path.name.startswith("eval-"):
            raise OSError("injected directory creation interruption")
        return original_mkdir(path, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(Path, "mkdir", fail)
        with pytest.raises(OSError, match="interruption"):
            e(point(e))
    assert calls == [] and e.n_evaluation_attempts == 1
    assert not list(e.workdir.glob("eval-*-*"))
    resumed = evaluator(tmp_path, 1, calls, resume=True)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        resumed(point(resumed))
    assert calls == []


def test_pending_reservation_slot_is_not_overwritten_when_budget_remains(tmp_path, monkeypatch):
    calls = []
    e = evaluator(tmp_path, 2, calls)
    original_mkdir = Path.mkdir

    def fail(path, *args, **kwargs):
        if path.name.startswith("eval-"):
            raise OSError("interruption")
        return original_mkdir(path, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(Path, "mkdir", fail)
        with pytest.raises(OSError):
            e(point(e))
    marker = (e.workdir/"attempt-000000.json").read_bytes()
    resumed = evaluator(tmp_path, 2, calls, resume=True)
    resumed(point(resumed, .01))
    assert calls == [1] and resumed.n_evaluation_attempts == 2
    assert (e.workdir/"attempt-000000.json").read_bytes() == marker


@pytest.mark.parametrize("field,value", [
    ("schema_version", 2), ("evaluation_index", True), ("evaluation_index", 99),
    ("budget_sha256", "a"*64), ("coordinate_sha256", "a"*64),
    ("coordinates", [0.]),
])
def test_corrupt_marker_refused_before_factory(tmp_path, field, value):
    e = evaluator(tmp_path, 2, [])
    e(point(e))
    marker = e.workdir/"attempt-000000.json"
    data = json.loads(marker.read_text())
    data[field] = value
    marker.write_text(json.dumps(data))
    calls = []
    with pytest.raises(ValueError, match="reservation"):
        evaluator(tmp_path, 2, calls, resume=True)
    assert calls == []


def test_missing_marker_is_not_invented_for_completed_data(tmp_path):
    e = evaluator(tmp_path, 1, [])
    e(point(e))
    # A known temporary fixture file only, never material output.
    (e.workdir/"attempt-000000.json").unlink()
    with pytest.raises(ValueError, match="lacks its budget reservation"):
        evaluator(tmp_path, 1, [], resume=True)


def test_completed_coordinate_must_match_its_reserved_slot(tmp_path):
    import hashlib
    e = evaluator(tmp_path, 1, [])
    e(point(e))
    result = next(e.workdir.glob("eval-*/result.json"))
    data = json.loads(result.read_text())
    changed = point(e, .01)
    data["coordinates"] = changed.tolist()
    data["coordinate_sha256"] = hashlib.sha256(changed.tobytes()).hexdigest()
    result.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="corrupt cached"):
        evaluator(tmp_path, 1, [], resume=True)


def test_stale_writer_cannot_reserve_same_index_for_different_geometry(tmp_path):
    calls = []
    first = evaluator(tmp_path, 2, calls)
    stale = evaluator(tmp_path, 2, calls, resume=True)
    first(point(first))
    marker = (first.workdir/"attempt-000000.json").read_bytes()
    with pytest.raises(FileExistsError):
        stale(point(stale, .01))
    assert calls == [0] and len(list(first.workdir.glob("eval-*-*"))) == 1
    assert (first.workdir/"attempt-000000.json").read_bytes() == marker


def test_conditional_grid_budget_is_checked_before_next_new_calculator(tmp_path):
    calls = []

    def flat(index, atoms, directory):
        calls.append(index)
        return SinglePointCalculator(atoms, energy=0., forces=np.zeros((len(atoms), 3)), stress=np.zeros(6))

    e = evaluator(tmp_path, 2, calls, factory=flat)
    size = e.chart.coordinate_count
    plane = ModePlane(reference=np.zeros(size), basis=np.eye(size), metric_weights=np.ones(size),
                      axis_weights=np.eye(size)[:, :2], axis_labels=("toy Q1", "toy Q2"),
                      amplitude_unit="toy", reference_id="toy-budget-not-HfO2")
    with pytest.raises(RuntimeError, match="budget exhausted"):
        sample_conditional_mode_surface(plane, [-.001, .001], [-.002, .002], e,
                                        energy_unit="toy", potential_kind="E", boundary_condition="toy")
    resumed = evaluator(tmp_path, 2, calls, resume=True, factory=flat)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        sample_conditional_mode_surface(plane, [-.001, .001], [-.002, .002], resumed,
                                        energy_unit="toy", potential_kind="E", boundary_condition="toy")
    assert calls == [0, 1] and resumed.n_evaluation_attempts == 2


@pytest.mark.parametrize("value", [True, 1., 2, "1"])
def test_budget_policy_type_or_value_drift_is_refused(tmp_path, value):
    e = evaluator(tmp_path, 1, [])
    policy = e.workdir/"evaluation-budget.json"
    data = json.loads(policy.read_text())
    data["max_evaluation_attempts"] = value
    policy.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="budget changed"):
        evaluator(tmp_path, 1, [], resume=True)


def test_truncated_reservation_cannot_trigger_retry_or_refund(tmp_path):
    e = evaluator(tmp_path, 1, [])
    e(point(e))
    (e.workdir/"attempt-000000.json").write_text("{")
    calls = []
    with pytest.raises(ValueError):
        evaluator(tmp_path, 1, calls, resume=True)
    assert calls == []


def test_zero_DFT_verifier_executes_real_budget_and_recovery_paths(tmp_path):
    from scripts.verify_mode_evaluation_budget import verify
    result = verify(tmp_path/"verify")
    assert result["factory_invocations"] == [0, 1]
    assert result["cumulative_reserved_attempts"] == 2
    assert result["completed_result_files"] == result["failed_result_files"] == 1
    assert result["physical_contract_sha256"] == "1308f3467e845dcfdb138871212baca9e0a6a05fcc78c281abf9b81e4212226f"
    assert result["new_material_DFT_calls"] == result["scheduler_mutations"] == 0
    with pytest.raises(FileExistsError, match="fresh workdir"):
        verify(tmp_path/"verify")


@pytest.mark.parametrize("invalid", [[], None, 1])
def test_non_object_marker_is_rejected(tmp_path, invalid):
    e = evaluator(tmp_path, 1, [])
    e(point(e))
    (e.workdir/"attempt-000000.json").write_text(json.dumps(invalid))
    with pytest.raises(ValueError, match="corrupt.*reservation"):
        evaluator(tmp_path, 1, [], resume=True)
