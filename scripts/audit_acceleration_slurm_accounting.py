"""Archive and audit Slurm ABACUS launches for the CPC acceleration table.

The capture subcommand is read-only on hf. The audit subcommand is fully
offline and counts completed Slurm steps through timestamps in the archived
optimizer logs. This measures process launches, not SCF-level convergence.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
import hashlib
import json
from pathlib import Path
import re
import subprocess


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "benchmarks" / "convergence" / "results"
JOBS = {
    "bto_global": (27675981, "bto_fire_baseline.log", "crossing"),
    "bto_block002": (27729004, "bto_blockfire_002.log", "crossing"),
    "bto_block001": (27729028, "bto_blockfire_001.log", "crossing"),
    "bto_scaled": (27729012, "bto_fire_scaled_002.log", "crossing"),
    "hfo2_global": (27678924, "hfo2_fire_baseline.log", "crossing"),
    "hfo2_staged_coarse": (27729066, "hfo2_fire_scaled.log", "stage_step_23"),
    "hfo2_staged_refine": (27729231, "hfo2_fire_staged_refine23.log", "crossing"),
}
FIELDS = ("job_id", "job_name", "state", "start", "end", "exit_code")
LOG_ROW = re.compile(
    r"^(?:FIRE|BlockFIRE):\s+(\d+)\s+(\d\d:\d\d:\d\d)\s+"
    r"[-+\d.eE]+\s+([-+\d.eE]+)\s*$"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def capture(raw: Path) -> None:
    if raw.exists():
        raise FileExistsError(raw)
    job_ids = ",".join(str(item[0]) for item in JOBS.values())
    command = [
        "ssh", "hf", "sacct", "-j", job_ids,
        "--format=JobID,JobName,State,Start,End,ExitCode", "-n", "-P",
    ]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    if not result.stdout.strip():
        raise ValueError("sacct returned no accounting rows")
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_text(result.stdout, encoding="utf-8", newline="\n")


def parse_accounting(raw_text: str) -> dict[str, dict]:
    jobs: dict[str, dict] = {}
    for line in raw_text.splitlines():
        fields = line.split("|")
        if len(fields) != len(FIELDS):
            raise ValueError(f"malformed sacct row: {line}")
        row = dict(zip(FIELDS, fields))
        label, separator, step = row["job_id"].partition(".")
        if not label.isdigit() or not label in {str(x[0]) for x in JOBS.values()}:
            raise ValueError(f"unexpected job ID: {row['job_id']}")
        bucket = jobs.setdefault(label, {"parent": None, "steps": []})
        if not separator:
            if bucket["parent"] is not None:
                raise ValueError(f"duplicate parent row: {label}")
            bucket["parent"] = row
        elif step.isdigit():
            bucket["steps"].append(row)
        elif step not in {"batch", "extern"}:
            raise ValueError(f"unexpected non-ABACUS step: {row['job_id']}")
    if set(jobs) != {str(x[0]) for x in JOBS.values()}:
        raise ValueError("missing job accounting")
    for label, bucket in jobs.items():
        parent = bucket["parent"]
        if parent is None or parent["state"] != "COMPLETED" or parent["exit_code"] != "0:0":
            raise ValueError(f"parent job not successfully completed: {label}")
        steps = sorted(bucket["steps"], key=lambda item: int(item["job_id"].split(".")[1]))
        indices = [int(row["job_id"].split(".")[1]) for row in steps]
        if indices != list(range(len(steps))) or not steps:
            raise ValueError(f"missing or repeated numbered Slurm step: {label}")
        for row in steps:
            if row["job_name"] != "abacus" or row["state"] != "COMPLETED" \
                    or row["exit_code"] != "0:0":
                raise ValueError(f"non-successful ABACUS launch: {row['job_id']}")
            start, end = datetime.fromisoformat(row["start"]), datetime.fromisoformat(row["end"])
            if end < start:
                raise ValueError(f"negative step time: {row['job_id']}")
        bucket["steps"] = steps
    return jobs


def cutoff_from_log(log: Path, policy: str, parent: dict) -> tuple[int, float, datetime]:
    rows = []
    for line in log.read_text(encoding="utf-8").splitlines():
        if match := LOG_ROW.match(line):
            rows.append((int(match[1]), match[2], float(match[3])))
    if not rows or [row[0] for row in rows] != list(range(len(rows))):
        raise ValueError(f"missing or noncontiguous optimizer rows: {log}")
    if policy == "crossing":
        candidates = [row for row in rows if row[2] <= 0.10]
        if not candidates:
            raise ValueError(f"no threshold crossing: {log}")
        step, clock, force = candidates[0]
    elif policy == "stage_step_23":
        step, clock, force = rows[23]
        if force <= 0.10:
            raise ValueError("coarse stage had already crossed the common threshold")
    else:
        raise ValueError(policy)
    start, end = datetime.fromisoformat(parent["start"]), datetime.fromisoformat(parent["end"])
    trial = datetime.fromisoformat(f"{start.date()}T{clock}")
    if trial < start:
        trial += timedelta(days=1)
    if not start <= trial <= end:
        raise ValueError(f"optimizer time outside Slurm parent interval: {log}")
    return step, force, trial


def audit(raw: Path) -> dict:
    accounting = parse_accounting(raw.read_text(encoding="utf-8"))
    routes = {}
    for name, (job_id, log_name, policy) in JOBS.items():
        bucket = accounting[str(job_id)]
        log = RESULTS / log_name
        step, force, cutoff = cutoff_from_log(log, policy, bucket["parent"])
        completed = [row for row in bucket["steps"]
                     if datetime.fromisoformat(row["end"]) <= cutoff]
        at_or_after = [row for row in bucket["steps"]
                       if datetime.fromisoformat(row["start"]) >= cutoff]
        if len(completed) + len(at_or_after) != len(bucket["steps"]):
            raise ValueError(f"ambiguous launch overlapping cutoff: {name}")
        routes[name] = {
            "slurm_job_id": job_id,
            "optimizer_log": log.relative_to(ROOT).as_posix(),
            "cutoff_kind": policy,
            "cutoff_step": step,
            "force_at_cutoff_eV_per_A": force,
            "cutoff_local_slurm_time": cutoff.isoformat(timespec="seconds"),
            "completed_launches_by_cutoff_second": len(completed),
            "launched_at_or_after_cutoff_second": len(at_or_after),
            "launched_in_same_cutoff_second": sum(
                datetime.fromisoformat(row["start"]) == cutoff for row in at_or_after
            ),
            "all_successful_abacus_launches": len(bucket["steps"]),
            "last_completed_step_id": completed[-1]["job_id"] if completed else None,
            "first_at_or_after_step_id": at_or_after[0]["job_id"] if at_or_after else None,
            "optimizer_log_sha256": sha256(log),
        }
    expected = {"bto_global": 205, "bto_block002": 79, "bto_block001": 44,
                "bto_scaled": 65, "hfo2_global": 217,
                "hfo2_staged_coarse": 122, "hfo2_staged_refine": 12}
    if any(routes[key]["completed_launches_by_cutoff_second"] != count
           for key, count in expected.items()):
        raise ValueError("Slurm launch totals disagree with the manuscript table")
    return {
        "kind": "seven_hf_slurm_jobs_acceleration_launch_accounting_not_scf_audit",
        "status": "all_numbered_abacus_steps_completed_and_cutoffs_timestamp_matched",
        "source_sha256": {"sacct_raw_psv": sha256(raw)},
        "sacct_fields": list(FIELDS),
        "routes": routes,
        "comparison": {
            "bto_block001_saved_percent": 100 * (205 - 44) / 205,
            "hfo2_staged_saved_percent": 100 * (217 - 122 - 12) / 217,
        },
        "limitations": [
            "Counts are successful ABACUS process launches, not individually SCF-audited evaluations.",
            "The threshold force is parsed from archived optimizer logs, not recomputed from raw image forces.",
            "Slurm and optimizer timestamps resolve only to seconds; starts in the cutoff second cannot be strictly ordered against the log line.",
            "A launch reduction is not a universal speedup or proof of identical minimum-energy paths.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    capture_parser = sub.add_parser("capture", help="read sacct on hf and archive its raw rows")
    capture_parser.add_argument("--raw", type=Path, required=True)
    audit_parser = sub.add_parser("audit", help="audit a saved raw sacct snapshot offline")
    audit_parser.add_argument("--raw", type=Path, required=True)
    audit_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.action == "capture":
        capture(args.raw)
        print(json.dumps({"raw": str(args.raw), "sha256": sha256(args.raw)}))
    else:
        if args.output.exists():
            raise FileExistsError(args.output)
        result = audit(args.raw)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result["comparison"]))


if __name__ == "__main__":
    main()
