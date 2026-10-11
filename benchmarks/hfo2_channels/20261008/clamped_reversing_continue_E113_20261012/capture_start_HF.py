"""One immutable, read-only native startup snapshot; never launch DFT/jobs."""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import re

HELPER = Path('/public/home/iai806/abacus/agent-runs/20261012-varneb-clamped-flip-continue-E103-r1/continue_flips_HF.py')
HELPER_SHA = '6c64c0155800082e28cbbb280511ac761929bdb6ccc5757fa8d68d156b32eac7'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def capture(m):
    parent = '28709791'
    admitted = m.accepted(parent)
    job = admitted['job_id']
    p = m.paths(parent)
    accounting = m.run(['sacct', '-X', '-n', '-P', '-j', job,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User,Start,End'])
    fields = accounting.split('|')
    if fields[:6] != [job, 'RUNNING', '0:0', '32', 'hfacnormal01', 'iai806']:
        raise ValueError('actual admitted running handle required; inspect other states')
    queue = m.run(['squeue', '-h', '-u', 'iai806', '-o', '%i|%100j|%T|%P|%C|%R'])
    committed = m.slots(queue, parent, spare=False)
    if job not in committed:
        raise ValueError('accepted handle not independently in committed queue')
    preflight = p['work'] / 'vcneb_preflight.json'
    info = json.loads(preflight.read_text())
    caches = info['calculator_parameters']['seed_cache_records']
    if (info['factory'] != 'examples.hfo2_fixed_input_factory:make_clamped_resume_cached_factory'
            or len(caches) != 9 or any(not isinstance(c, dict) or not c for c in caches)):
        raise ValueError('exact nine-cache continuation startup required')
    log = (p['work'] / 'vcneb.opt.log').read_text()
    if not any(row.split()[:2] == ['CheckedFIRE:', '0'] for row in log.splitlines()):
        raise ValueError('complete initial chain required')
    raw = p['work'] / 'image_0001/scf_000000/OUT.ABACUS/running_scf.log'
    dsize = re.findall(r'DSIZE\s*=\s*(\d+)', raw.read_text())
    if not dsize or set(dsize) != {'32'}:
        raise ValueError('actual native first ABACUS DSIZE32 required')
    first = p['work'] / 'image_0001/scf_000000/call_audit.json'
    first_info = json.loads(first.read_text()) if first.exists() else None
    root = m.ROOT / parent
    return dict(check_UTC=datetime.now(timezone.utc).isoformat(), source_job_id=job,
        original_parent=parent, native_accounting=accounting, native_queue=queue,
        native_control=m.run(['scontrol', 'show', 'job', '-o', job]),
        committed_study_jobs=committed, workdir=str(p['work']),
        runtime_preflight_sha256=sha(preflight), factory=info['factory'],
        seed_cache_count=len(caches), none_seed_cache_count=0,
        optimizer_log_raw=log, native_first_DSIZE=dsize,
        fresh_call_audits_at_capture=len(list(p['work'].glob('image_*/scf_*/call_audit.json'))),
        first_fresh_call_audit_exists=first_info is not None,
        first_fresh_call_audit_sha256=sha(first) if first_info is not None else None,
        first_fresh_call_elapsed_seconds=first_info['elapsed_seconds'] if first_info else None,
        stderr_bytes=p['stderr'].stat().st_size,
        failure_exists=(p['work']/'vcneb_failure.json').exists(),
        terminal_summary_exists=(p['work']/'vcneb_summary.json').exists(),
        admission_receipt_sha256=sha(root/'dispatch_receipt.json'),
        admission_journal_sha256=sha(root/'dispatch_journal.jsonl'),
        parent_audit_sha256=sha(p['terminal']/'audit_receipt.json'),
        executed_capture_sha256=sha(Path(__file__)),
        startup_is_not_ordinary_convergence_or_full_fresh_call_audit=True,
        new_readonly_DFT_calls=0, scheduler_writes_for_capture=0,
        running_source_or_physical_inputs_changed=False, holdout_generated_or_read=False)


if __name__ == '__main__':
    if sha(HELPER) != HELPER_SHA:
        raise ValueError('pinned accepted-admission helper changed')
    spec = importlib.util.spec_from_file_location('E113_pinned_admission', HELPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    output = module.ROOT / '28709791/startup_observation_E113.json'
    if output.exists():
        raise FileExistsError('inspect existing snapshot; never overwrite/replay')
    report = capture(module)
    with output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(output=str(output), sha256=sha(output),
        job_id=report['source_job_id'], native_DSIZE=report['native_first_DSIZE'])))
