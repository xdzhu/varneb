import importlib.util
import json
from pathlib import Path

import pytest


def helper():
    path = Path(__file__).resolve().parents[1]/"benchmarks/hfo2_channels/20261008/clamped_M_terminal_E068_20261011/release_existing_PO_to_T.py"
    spec = importlib.util.spec_from_file_location("rolling_slot_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def receipt():
    return dict(job_id="28722320", ordinary_converged=True,
        status="audited_ordinary_converged_not_TS_or_sampling_certificate",
        new_DFT_calls=0, scheduler_mutations=0, physical_inputs_or_running_source_changed=False,
        holdout_generated_or_read=False,
        all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked=True,
        terminal_observation=dict(ordinary_residual_pass=True, replayed_fmax_eV_A=.09))


def test_terminal_gate():
    helper().verify_terminal_receipt(receipt())


@pytest.mark.parametrize("field,value", [("ordinary_converged",False), ("job_id","28661020"),
    ("new_DFT_calls",1), ("scheduler_mutations",1), ("physical_inputs_or_running_source_changed",True),
    ("holdout_generated_or_read",True), ("status","max_steps_reached"),
    ("all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked",False)])
def test_no_unconverged_or_changed_evidence(field, value):
    report = receipt()
    report[field] = value
    with pytest.raises(ValueError, match="converged terminal"):
        helper().verify_terminal_receipt(report)


@pytest.mark.parametrize("force", [.11, -.1, True, None, float("nan"), float("inf")])
def test_no_force_flag_only_release(force):
    report = receipt()
    report["terminal_observation"]["replayed_fmax_eV_A"] = force
    with pytest.raises(ValueError):
        helper().verify_terminal_receipt(report)


def test_single_spare_slot_ignores_unrelated_projects():
    body = "28722810|hfo2-G2-flip-E067|RUNNING|hfacnormal01|32\n28692775|hfo2-G2-E058-PO_to_T|PENDING|hfacnormal01|32\n999|pc-user-job|RUNNING|other|128"
    assert helper().single_spare_slot(body) == ["28722810"]


def test_empty_queue_only_with_positive_terminal_flip_accounting():
    assert helper().single_spare_slot("", terminal_flip=True) == []
    with pytest.raises(ValueError):
        helper().single_spare_slot("", terminal_flip=False)


@pytest.mark.parametrize("body", ["", "28722320|hfo2-G2-M-E062|COMPLETING|hfacnormal01|32\n28722810|hfo2-G2-flip-E067|RUNNING|hfacnormal01|32",
    "28722810|hfo2-G2-flip-E067|RUNNING|other|32", "28722810|hfo2-G2-flip-E067|RUNNING|hfacnormal01|1",
    "999|hfo2-unregistered|RUNNING|hfacnormal01|32", "broken"])
def test_unknown_busy_or_absent_slot_refuses(body):
    with pytest.raises(ValueError):
        helper().single_spare_slot(body)


def test_previous_intent_never_repeated(tmp_path, monkeypatch):
    module = helper()
    monkeypatch.setattr(module, "ROOT", tmp_path)
    (tmp_path/"rolling_handoff_journal.jsonl").write_text('{}\n')
    monkeypatch.setattr(module, "load_rules", lambda: pytest.fail("no scheduler interaction after an existing intent"))
    with pytest.raises(FileExistsError, match="never repeat"):
        module.release_once("a"*64)


def test_does_not_submit_or_release_other_waiter():
    source = Path(helper().__file__).read_text()
    assert "sbatch" not in source
    assert "['scontrol','release',TARGET]" in source
    assert "['scontrol','release',OTHER_HELD]" not in source
    assert "'Dependency='+dependency" in source
    assert "dependency = 'afterok:'+DONE" in source
    assert "rules.run(['scontrol','update'" in source


@pytest.fixture
def reviewed_workflow(tmp_path, monkeypatch):
    module = helper()
    root = tmp_path/"audit"
    root.mkdir()
    wait = tmp_path/"wait"
    wait.mkdir()
    monkeypatch.setattr(module, "ROOT", root)
    monkeypatch.setattr(module, "WAIT_ROOT", wait)
    work = tmp_path/"old/PO_to_M/band"
    work.mkdir(parents=True)
    (work/"vcneb_summary.json").write_text('{"synthetic":true}\n')
    report = {**receipt(), "workdir": str(work), "terminal_summary_sha256": module.sha(work/"vcneb_summary.json")}
    (root/"audit_receipt.json").write_text(json.dumps(report))
    source = tmp_path/"source"
    pilot = source/"cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm"
    pilot.parent.mkdir(parents=True)
    pilot.write_text("# synthetic tested recipe\n")
    binary = tmp_path/"binary"
    binary.write_text("synthetic nonexecutable\n")
    seed = tmp_path/"seed"
    seed.mkdir()
    (seed/"manifest.json").write_text('{}\n')
    (seed/"seed.traj").write_bytes(b"synthetic guard fixture only")
    preflight = dict(source_root=str(source), pilot_sha256=module.sha(pilot),
        ABACUS_binary_sha256=module.sha(binary), source_code_sha256={}, reports=[dict(
        channel="PO_to_T", seed_root=str(seed), manifest_sha256=module.sha(seed/"manifest.json"),
        seed_files_sha256={"seed.traj": module.sha(seed/"seed.traj")})])
    (wait/"preflight.json").write_text(json.dumps(preflight))
    monkeypatch.setattr(module, "PREFLIGHT_SHA", module.sha(wait/"preflight.json"))
    real_rule_path = Path(__file__).resolve().parents[1]/"benchmarks/hfo2_channels/20261008/clamped_flip_continue_E067_20261011/submit_held_once.py"
    spec = importlib.util.spec_from_file_location("E067_synthetic_workflow_rules", real_rule_path)
    rules = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rules)
    monkeypatch.setattr(rules, "BINARY", binary)
    monkeypatch.setattr(module, "load_rules", lambda: rules)
    monkeypatch.setattr(module, "pending_seed_preflight", lambda path: dict(synthetic=True, new_DFT_calls=0))
    state = dict(updated=False, writes=[], all_calls=[], stale=False, ambiguous_release=False, terminal_flip=False)
    def body(job):
        channel = rules.WAITERS[job]
        dependency = "(null)" if state["terminal_flip"] else ("afterok:28722320" if job == module.TARGET and state["updated"] and not state["stale"] else "afterok:28722320,afterok:28722810")
        values = dict(JobId=job, JobName="hfo2-G2-E058-"+channel, UserId="iai806(16284)",
            JobState="PENDING", Partition="hfacnormal01", NumCPUs="32", Priority="0", Reason="JobHeldUser",
            Command=rules.PILOT, StdOut=rules.WAIT_ROOT+"/"+channel+".slurm.out",
            StdErr=rules.WAIT_ROOT+"/"+channel+".slurm.err", Dependency=dependency)
        return " ".join(f"{k}={v}" for k,v in values.items())
    def run(argv):
        state["all_calls"].append(argv)
        if argv[0] == "sacct":
            if argv[argv.index("-j")+1] == module.LIVE:
                status = "COMPLETED" if state["terminal_flip"] else "RUNNING"
                return module.LIVE+"|"+status+"|0:0|32|hfacnormal01|iai806"
            return module.DONE+"|COMPLETED|0:0|32|hfacnormal01|iai806"
        if argv[0] == "squeue":
            return "" if state["terminal_flip"] else "28722810|hfo2-G2-flip-E067|RUNNING|hfacnormal01|32"
        if argv[:3] == ["scontrol","show","job"]:
            return body(argv[-1])
        state["writes"].append(argv)
        if argv[:2] == ["scontrol","update"]:
            state["updated"] = True
        elif argv[:2] == ["scontrol","release"] and state["ambiguous_release"]:
            raise TimeoutError("synthetic ambiguous accepted release")
        return ""
    monkeypatch.setattr(rules, "run", run)
    return module, state, module.sha(root/"audit_receipt.json")


def test_exact_existing_handle_changed_and_released_once(reviewed_workflow):
    module, state, digest = reviewed_workflow
    module.release_once(digest)
    assert state["writes"] == [["scontrol","update","JobId=28692775","Dependency=afterok:28722320"],
                               ["scontrol","release","28692775"]]
    assert all(argv[0] != "sbatch" for argv in state["all_calls"])
    with pytest.raises(FileExistsError):
        module.release_once(digest)
    assert len(state["writes"]) == 2


def test_stale_ack_never_repeats_update_or_releases(reviewed_workflow, monkeypatch):
    module, state, digest = reviewed_workflow
    state["stale"] = True
    times = iter([0.,100.])
    monkeypatch.setattr(module.time, "monotonic", lambda: next(times))
    with pytest.raises(TimeoutError, match="never repeat"):
        module.release_once(digest)
    assert state["writes"] == [["scontrol","update","JobId=28692775","Dependency=afterok:28722320"]]
    with pytest.raises(FileExistsError):
        module.release_once(digest)
    assert len(state["writes"]) == 1


def test_ambiguous_release_protected_against_retry(reviewed_workflow):
    module, state, digest = reviewed_workflow
    state["ambiguous_release"] = True
    with pytest.raises(TimeoutError):
        module.release_once(digest)
    with pytest.raises(FileExistsError):
        module.release_once(digest)
    assert len(state["writes"]) == 2


def test_both_parents_terminal_release_only_no_dependency_rewrite(reviewed_workflow):
    module, state, digest = reviewed_workflow
    state["terminal_flip"] = True
    module.release_once(digest)
    assert state["writes"] == [["scontrol","release","28692775"]]
    result = json.loads((module.ROOT/"rolling_handoff_receipt.json").read_text())
    assert result["dependency_writes"] == result["new_submissions"] == 0
    assert result["flip_accounting_state"] == "COMPLETED"


def test_actual_single_release_and_independent_startup():
    module = helper()
    case = Path(module.__file__).parent
    result = json.loads((case/"rolling_handoff_receipt.json").read_text())
    rows = [json.loads(line) for line in (case/"rolling_handoff_journal.jsonl").read_text().splitlines()]
    assert result["terminal_audit_sha256"] == module.sha(case/"audit_receipt.json")
    assert result["executed_helper_sha256"] == module.sha(Path(module.__file__))
    assert result["journal_sha256"] == module.sha(case/"rolling_handoff_journal.jsonl")
    assert result["new_submissions"] == result["dependency_writes"] == result["new_independent_chains"] == 0
    assert result["release_writes"] == 1 and result["target_job"] == "28692775"
    assert [r["job_id"] for r in rows if r["event"] == "existing_target_released"] == ["28692775"]
    assert not any(r["event"] == "dependency_update_intent" for r in rows)
    preflight = rows[0]["pending_seed_preflight"]
    assert preflight["geometry_images_checked"] == 9 and preflight["fresh_interiors"] == 7
    assert len(preflight["cached_endpoints"]) == 2 and preflight["new_DFT_calls"] == 0
    assert all(r["original_six_physical_bytes_native_DSIZE32_full_EFS_checked"] for r in preflight["cached_endpoints"])
    startup = json.loads((case/"startup_observation.json").read_text())
    assert startup["raw_accounting"] == "28692775|RUNNING|0:0|32|2026-10-11T10:17:58|hfacnormal01"
    assert "DSIZE = 32" in startup["native_first_fresh_SCF_DSIZE_line"]
    assert startup["current_running_study_allocations"] == 1 and startup["stderr_bytes_at_capture"] == 0
    assert not startup["first_complete_step_zero_or_fresh_SCF_completion_claimed"]
