"""Offline transaction tests: no native scheduler calls in pytest."""
import importlib.util
import hashlib
import json
from pathlib import Path
import pytest

SOURCE = Path(__file__).resolve().parents[1]/"benchmarks/hfo2_channels/20261008/clamped_capacity_E074_20261011/hold_pending_once.py"
spec = importlib.util.spec_from_file_location("capacity_E074", SOURCE)
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)


def body(job, held=False):
    channel = guard.CHANNELS[job]
    first = guard.DEPS[job]
    second = "28692776" if first == "28692775" else "28709789"
    return (f"JobId={job} JobName=hfo2-G2-E060-{channel} UserId=iai806(16284) "
            f"JobState=PENDING Partition=hfacnormal01 NumCPUs=32 NumTasks=1 CPUs/Task=32 "
            f"Command={guard.PILOT} WorkDir=/public/home/iai806 "
            f"StdOut={guard.ROOT}/{channel}.slurm.out StdErr={guard.ROOT}/{channel}.slurm.err "
            f"Dependency=afterok:{first}(unfulfilled),afterok:{second}(unfulfilled) "
            f"Reason={'JobHeldUser' if held else 'Dependency'} Priority={0 if held else 1163}")


def test_exact_transaction_and_exclusive_replay(tmp_path):
    held = set()
    commands = []
    def run(command):
        commands.append(command)
        job = command.split()[-1]
        if command.startswith("scontrol hold "):
            held.add(job)
            return ""
        return body(job, held=job in held)
    out = tmp_path/"new"
    receipt = guard.hold_once(out, run=run)
    assert held == set(guard.CHANNELS)
    assert [c for c in commands if c.startswith("scontrol hold ")] == [
        "scontrol hold "+j for j in guard.CHANNELS]
    assert receipt["dependency_writes"] == receipt["new_submissions"] == 0
    with pytest.raises(FileExistsError):
        guard.hold_once(out, run=run)
    assert len(commands) == 16


@pytest.mark.parametrize("old,new", [("PENDING", "RUNNING"), ("32", "64"),
    ("iai806(", "another("), ("hfacnormal01", "other"),
    ("(unfulfilled)", "(fulfilled)"), ("Priority=1163", "Priority=0")])
def test_wrong_native_identity_rejected_before_any_write(tmp_path, old, new):
    commands = []
    def run(command):
        commands.append(command)
        return body(command.split()[-1]).replace(old, new)
    with pytest.raises(ValueError):
        guard.hold_once(tmp_path/"new", run=run)
    assert not any(c.startswith("scontrol hold ") for c in commands)


def test_ambiguous_hold_is_not_repeated(tmp_path):
    commands = []
    def run(command):
        commands.append(command)
        if command.startswith("scontrol hold "):
            raise TimeoutError("transport lost after possible write")
        return body(command.split()[-1])
    with pytest.raises(TimeoutError):
        guard.hold_once(tmp_path/"new", run=run)
    assert sum(c.startswith("scontrol hold ") for c in commands) == 1
    assert "reconcile_do_not_rerun" in (tmp_path/"new/journal.jsonl").read_text()


def test_actual_native_receipt_retains_dependencies_and_exact_four_targets():
    case = SOURCE.parent/"native_transaction"
    events = [json.loads(line) for line in (case/"journal.jsonl").read_text().splitlines()]
    assert events[0]["helper_sha256"] == hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    before = {row["job_id"]: guard.verify(row["raw"], row["job_id"])
              for row in events if row["event"] == "preflight"}
    after = {row["job_id"]: guard.verify(row["raw"], row["job_id"], held=True)
             for row in events if row["event"] == "hold_verified"}
    assert set(before) == set(after) == set(guard.CHANNELS)
    assert sum(row["event"] == "hold_intent" for row in events) == 4
    assert all(before[j]["Dependency"] == after[j]["Dependency"] for j in before)
    assert not any(row["event"] == "reconcile_do_not_rerun" for row in events)
    receipt = json.loads((case/"receipt.json").read_text())
    assert receipt["new_submissions"] == receipt["dependency_writes"] == 0
    assert receipt["running_jobs_changed"] is False
