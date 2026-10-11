"""Audit a completed ordinary G2 segment, without DFT or scheduler writes.

Only the terminal convergence gate is certified here. Discrete image peaks
are not continuous saddle points or measured barrier-error bounds.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import shlex
import subprocess
import tarfile

import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.export_hfo2_clamped_observation import export
from scripts.analyze_hfo2_clamped_residual import analyze


def completed_step(log):
    steps = []
    for line in log.read_text().splitlines():
        fields = line.split()
        if len(fields) == 5 and fields[0].endswith(":"):
            steps.append(int(fields[1]))
    if not steps or steps != list(range(steps[-1] + 1)):
        raise ValueError("one consecutive complete optimizer history required")
    return steps[-1]


def validate_terminal(job_id, accounting, summary, failure_exists, *, outcome="converged"):
    if not isinstance(job_id, str) or not job_id.isdecimal():
        raise ValueError("explicit numeric registered job ID required")
    if accounting != job_id + "|COMPLETED|0:0|32|hfacnormal01|iai806":
        raise ValueError("actual completed registered 32-rank HF segment required")
    force = summary.get("final_max_generalized_force_eV_per_A")
    finite = not isinstance(force,bool) and isinstance(force,(int,float)) and np.isfinite(force)
    if outcome == "converged":
        if (failure_exists or not finite or not 0 <= force <= .10
                or summary.get("climbing_image_active_final") is not False
                or summary.get("status") != "converged" or summary.get("converged") is not True
                or summary.get("termination") != "force_threshold"):
            raise ValueError("ordinary .10 convergence, not exit-zero or a step cap, required")
    elif outcome == "step-cap":
        if (failure_exists or not finite or force <= .10
                or summary.get("climbing_image_active_final") is not False
                or summary.get("status") != "max_steps_reached" or summary.get("converged") is not False
                or summary.get("termination") != "step_limit"):
            raise ValueError("explicit successful unconverged step cap required; no blind failed-run restart")
    else:
        raise ValueError("explicit converged or step-cap outcome required")


def validate_job_source(body, job_id, work, production_script):
    row = dict(token.split("=", 1) for token in body.split() if "=" in token)
    prefix = "flip" if "flip_continue" in production_script.name else work.parent.name
    root = work.parent.parent
    if (row.get("JobId") != job_id or row.get("JobState") != "COMPLETED"
            or not row.get("JobName", "").startswith("hfo2-")
            or not row.get("UserId", "").startswith("iai806(")
            or row.get("Partition") != "hfacnormal01" or row.get("NumCPUs") != "32"
            or row.get("Command") != str(production_script)
            or row.get("StdOut") != str(root/(prefix+".slurm.out"))
            or row.get("StdErr") != str(root/(prefix+".slurm.err"))):
        raise ValueError("registered job and native output/source namespaces differ")
    return row


def validate_accounting_source(body, job_id, work, production_script):
    """Use persistent native accounting after scontrol's terminal entry expires."""
    values = body.strip().split("|")
    if len(values) != 9:
        raise ValueError("complete native accounting source row required")
    job, state, user, cpus, partition, name, stdout, stderr, submit = values
    argv = shlex.split(submit)
    exports = [v[len("--export="):] for v in argv if v.startswith("--export=")]
    if len(exports) != 1:
        raise ValueError("one native submitted runtime environment required")
    environment = dict(v.split("=", 1) for v in exports[0].split(",") if "=" in v)
    if (not argv or argv[0] != "sbatch" or argv[-1] != str(production_script)
            or environment.get("WORKDIR") != str(work) or environment.get("RUN_DFT") != "1"):
        raise ValueError("native SubmitLine differs from audited work/production script")
    row = dict(JobId=job, JobState=state, UserId=user+"(accounting)", NumCPUs=cpus,
        Partition=partition, JobName=name, Command=argv[-1], StdOut=stdout, StdErr=stderr)
    return validate_job_source(" ".join(f"{k}={v}" for k,v in row.items()), job_id, work, production_script)


def verify_runtime(source, archive, archive_sha256):
    if sha256(archive) != archive_sha256:
        raise ValueError("registered immutable runtime archive changed")
    code = {}
    with tarfile.open(archive) as bundle:
        for member in bundle.getmembers():
            if (member.isfile() and member.name.endswith(".py")
                    and member.name.startswith(("vcneb/", "scripts/", "examples/"))):
                content = bundle.extractfile(member).read()
                if content != (source/member.name).read_bytes():
                    raise ValueError("immutable runtime source changed: " + member.name)
                code[member.name] = sha256(source/member.name)
    if not code or "vcneb/material_runner.py" not in code:
        raise ValueError("complete registered runtime, not an empty archive, required")
    return code


def audit(work, job_id, production_script, source, archive, archive_sha256, output, *, outcome="converged", seed_latest=False):
    if output.exists():
        raise FileExistsError("inspect an existing audit; do not overwrite or repeat it")
    accounting = subprocess.run(["sacct", "-X", "-n", "-P", "-j", job_id,
        "--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User"],
        check=True, capture_output=True, text=True, timeout=45).stdout.strip()
    summary_path = work/"vcneb_summary.json"
    summary_sha = sha256(summary_path)
    summary = json.loads(summary_path.read_text())
    validate_terminal(job_id, accounting, summary, (work/"vcneb_failure.json").exists(), outcome=outcome)
    if seed_latest and outcome != "step-cap":
        raise ValueError("a new geometry seed requires an explicit unconverged step cap")
    job_source = subprocess.run(["sacct", "-X", "-n", "-P", "-j", job_id,
        "--format=JobIDRaw,State,User,AllocCPUS,Partition,JobName%100,StdOut%1500,StdErr%1500,SubmitLine%4000"],
        check=True, capture_output=True, text=True, timeout=45).stdout.strip()
    validate_accounting_source(job_source, job_id, work, production_script)
    code = verify_runtime(source, archive, archive_sha256)
    step = completed_step(work/"vcneb.opt.log")
    calls = sorted(work.glob("image_*/scf_*/call_audit.json"))
    if len(calls) != 7*step:
        raise ValueError("fresh seven-interior evaluation count differs; inspect before proceeding")
    actual = subprocess.run
    def forbid(*args, **kwargs):
        raise RuntimeError("terminal audit forbids external launches and scheduler writes")
    subprocess.run = forbid
    try:
        call_reports = []
        for path in calls:
            record = json.loads(path.read_text())
            inputs = {n: sha256(path.parent/n) for n in (*CONTRACT, "STRU")}
            if (inputs != record["input_sha256"]
                    or any(inputs[n] != h for n, h in CONTRACT.items())
                    or sha256(path.parent/"OUT.ABACUS/running_scf.log") != record["raw_log_sha256"]):
                raise ValueError("native physical bytes or raw log differ from the call audit")
            raw = audited_results(path.parent)
            if any(not np.allclose(raw[k], record["results"][k], atol=1e-12, rtol=0)
                   for k in ("energy", "forces", "stress")):
                raise ValueError("native energy/forces/stress differ from pinned records")
            elapsed = record["elapsed_seconds"]
            if isinstance(elapsed, bool) or not np.isfinite(elapsed) or elapsed < 0:
                raise ValueError("recorded finite nonnegative SCF elapsed time required")
            call_reports.append(dict(path=str(path), audit_sha256=sha256(path),
                raw_log_sha256=record["raw_log_sha256"], input_sha256=inputs,
                elapsed_seconds=elapsed))
        terminal = output/"observations"/f"step_{step:04d}"
        result = export(work, step, terminal, job_id, production_script=production_script)
        if (result["ordinary_residual_pass"] != (outcome == "converged")
                or abs(result["replayed_fmax_eV_A"]-summary["final_max_generalized_force_eV_per_A"]) > 1e-10):
            raise ValueError("actual terminal force differs from same-boundary replay")
        residual = analyze(terminal)
        with (terminal/"residual.json").open("x") as stream:
            json.dump(residual, stream, indent=2, allow_nan=False)
        for image in result["raw_image_evaluations"]:
            raw = Path(image["raw_source"])
            dest = terminal/"raw"/f"image_{image['image_index']:04d}"
            (dest/"OUT.ABACUS").mkdir(parents=True)
            for name in ("INPUT", "KPT", "STRU", "OUT.ABACUS/running_scf.log"):
                shutil.copyfile(raw/name, dest/name)
            shutil.copyfile(image["audit_path"], dest/"source_audit.json")
        seed_manifest = None
        if seed_latest:
            from ase.io import read
            from scripts.prepare_hfo2_clamped_resume import prepare
            from examples.hfo2_fixed_input_factory import make_clamped_resume_cached_factory
            seed = output/"seed"
            prepare(terminal, seed)
            parameters = json.loads((seed/"factory_parameters.json").read_text())
            images = read(seed/"seed.traj", index=":")
            factory = make_clamped_resume_cached_factory(parameters=parameters,
                command="mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus")
            for i, image in enumerate(images):
                old = dict(image.calc.results)
                image.calc = factory(i, image, output/"cache_preflight"/f"image_{i:04d}")
                for key,value in (("energy",image.get_potential_energy()),
                        ("forces",image.get_forces()),("stress",image.get_stress())):
                    if image.calc.next_call or not np.allclose(value,old[key],atol=1e-12,rtol=0):
                        raise ValueError("all-nine exact latest-frame cache replay failed")
            seed_manifest = sha256(seed/"manifest.json")
    finally:
        subprocess.run = actual
    if sha256(summary_path) != summary_sha or (work/"vcneb_failure.json").exists():
        raise ValueError("terminal source changed during audit")
    shutil.copyfile(summary_path, output/"vcneb_summary.json")
    shutil.copyfile(work/"vcneb.opt.log", output/"vcneb.opt.log")
    receipt = dict(status=("audited_ordinary_converged_not_TS_or_sampling_certificate" if outcome == "converged"
                          else "audited_step_cap_not_converged_same_chain_seed"),
        check_UTC=datetime.now(timezone.utc).isoformat(), job_id=job_id,
        workdir=str(work), scheduler_row=accounting, scheduler_job_source_raw=job_source,
        scheduler_job_source_provider="persistent native sacct StdOut/StdErr/SubmitLine", terminal_step=step,
        terminal_summary_sha256=summary_sha, ordinary_converged=(outcome == "converged"),
        terminal_observation=result, fresh_interior_SCFs=len(calls),
        all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked=True,
        fresh_SCF_wall_seconds_sum=sum(r["elapsed_seconds"] for r in call_reports),
        fresh_SCF_core_hour_proxy_32=sum(r["elapsed_seconds"] for r in call_reports)*32/3600,
        runtime_archive_sha256=archive_sha256, runtime_code_sha256=code,
        executed_auditor_sha256=sha256(Path(__file__)),
        executed_exporter_sha256=result["export_script_sha256"],
        executed_residual_analyzer_sha256=sha256(Path(__file__).with_name("analyze_hfo2_clamped_residual.py")),
        calls=call_reports, new_DFT_calls=0, scheduler_mutations=0,
        physical_inputs_or_running_source_changed=False, holdout_generated_or_read=False,
        continuous_TS_or_barrier_error_bound_certified=False)
    if seed_latest:
        receipt.update(seed_manifest_sha256=seed_manifest, current_frame_exact_caches_checked=9,
                       new_independent_chains=0, FIRE_state_restored=False)
    with (output/"audit_receipt.json").open("x") as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("workdir", "production-script", "runtime-source", "runtime-archive", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--runtime-archive-sha256", required=True)
    parser.add_argument("--outcome", choices=("converged","step-cap"), default="converged")
    parser.add_argument("--seed-latest", action="store_true")
    args = parser.parse_args()
    receipt = audit(args.workdir, args.job_id, args.production_script, args.runtime_source,
                    args.runtime_archive, args.runtime_archive_sha256, args.output,
                    outcome=args.outcome, seed_latest=args.seed_latest)
    print(json.dumps({k: receipt[k] for k in ("status", "job_id", "terminal_step", "fresh_interior_SCFs")}))


if __name__ == "__main__":
    main()
