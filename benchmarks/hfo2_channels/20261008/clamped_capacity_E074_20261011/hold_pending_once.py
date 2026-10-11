"""Single-use hold of four exact E060 waiters; never submit/release a job.

Run locally with a NEW receipt directory. If interrupted, reconcile by native
reads; never rerun this transaction. No production source or inputs are edited.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

CHANNELS = {
    "28709788": "PO_to_T", "28709789": "PO_to_M",
    "28709790": "PO_flip_T_pattern_preserving",
    "28709791": "PO_flip_T_pattern_reversing",
}
ROOT = "/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E060-r1"
PILOT = "/public/home/iai806/abacus/agent-runs/20261010-varneb-nested-controls-E055-r1/source/cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm"
DEPS = {"28709788": "28692775", "28709789": "28692775",
        "28709790": "28709788", "28709791": "28709788"}


def verify(body, job, *, held=False):
    row = dict(token.split("=", 1) for token in body.split() if "=" in token)
    channel = CHANNELS[job]
    first = DEPS[job]
    second = "28692776" if first == "28692775" else "28709789"
    expected = {
        "JobId": job, "JobName": "hfo2-G2-E060-"+channel,
        "JobState": "PENDING", "Partition": "hfacnormal01",
        "NumCPUs": "32", "NumTasks": "1", "CPUs/Task": "32",
        "Command": PILOT, "WorkDir": "/public/home/iai806",
        "StdOut": ROOT+"/"+channel+".slurm.out",
        "StdErr": ROOT+"/"+channel+".slurm.err",
        "Dependency": f"afterok:{first}(unfulfilled),afterok:{second}(unfulfilled)",
    }
    if any(row.get(key) != value for key, value in expected.items()):
        raise ValueError("exact pending identity/dependency changed; do not mutate")
    if not row.get("UserId", "").startswith("iai806("):
        raise ValueError("wrong owner")
    if held:
        if row.get("Reason") != "JobHeldUser" or row.get("Priority") != "0":
            raise ValueError("hold unconfirmed; reconcile by reads only")
    elif row.get("Reason") != "Dependency" or int(row.get("Priority", "0")) <= 0:
        raise ValueError("not an unheld dependency waiter")
    return row


def ssh(command):
    result = subprocess.run(["ssh", "hf", command], capture_output=True,
                            text=True, timeout=45, check=True)
    return result.stdout.strip()


def hold_once(output, *, run=ssh):
    output.mkdir()  # Exclusive: never overwrite/replay a transaction.
    events = []

    def record(event, **values):
        entry = dict(event=event, checked_UTC=datetime.now(timezone.utc).isoformat(), **values)
        events.append(entry)
        with (output/"journal.jsonl").open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(entry)+"\n")
            stream.flush()
            import os
            os.fsync(stream.fileno())

    record("intent", helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           targets=list(CHANNELS), maximum_study_allocations=2)
    try:
        for job in CHANNELS:  # All identities must pass before the first write.
            body = run("scontrol show job -o "+job)
            verify(body, job)
            record("preflight", job_id=job, raw=body)
        for job in CHANNELS:
            body = run("scontrol show job -o "+job)
            verify(body, job)
            record("hold_intent", job_id=job, raw=body)
            reply = run("scontrol hold "+job)
            body = run("scontrol show job -o "+job)
            verify(body, job, held=True)
            record("hold_verified", job_id=job, command_stdout=reply, raw=body)
    except Exception as exc:
        record("reconcile_do_not_rerun", error=type(exc).__name__, message=str(exc))
        raise
    receipt = dict(status="four_exact_pending_starters_user_held", held_jobs=list(CHANNELS),
                   running_jobs_changed=False, dependency_writes=0, new_submissions=0,
                   physical_input_changes=0, production_source_changes=0,
                   new_DFT_calls=0, holdout_generated_or_read=False,
                   policy="After terminal material audit, release at most one exact eligible waiter per verified free study slot.")
    with (output/"receipt.json").open("x", encoding="utf-8") as stream:
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    return receipt


if __name__ == "__main__":
    print(json.dumps(hold_once(Path(sys.argv[1]))))
