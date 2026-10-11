"""Recheck/release only the two existing +1% G2 flip handles.

Preflight never releases or runs DFT. Release requires audited zero-strain
flips and a free committed study slot. One exclusive journal per target
prevents retries after an ambiguous scheduler reply. No sbatch/update.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys

RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS/'20261011-varneb-clamped-held-flips-E100-r1'
WAIT = RUNS/'20261010-varneb-clamped-G2-E060-r1'
SOURCE = RUNS/'20261010-varneb-nested-controls-E055-r1/source'
PILOT = SOURCE/'cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm'
COMMON = RUNS/'20261011-varneb-clamped-M-terminal-E068-source-r2/release_existing_PO_to_T.py'
COMMON_SHA = '5ddca4a3bf73ed8ad1d8f348544d97bf170ef442354a63a1187ae49a5c56c0f0'
PREFLIGHT_SHA = 'b780585fbd5ccd05e85282f92dc746a5a79d48f6b25d923ff0bf3150098a2e8a'
CHANNELS = {'28709790': 'PO_flip_T_pattern_preserving',
            '28709791': 'PO_flip_T_pattern_reversing'}
LIVE = {'28727368': '28709789', '28727382': '28723655'}
PRESERVING_AUDIT = RUNS/'20261011-varneb-clamped-child-audit-E093-r1/28727382/terminal/audit_receipt.json'
REVERSING_AUDIT = RUNS/'20261011-varneb-clamped-reversing-terminal-E083-r1/audit_receipt.json'
REVERSING_SHA = '632e0e9b2b821b3bfc7d66b361c44c7db734988fb5d018ea2f31441fab14ed41'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(argv):
    return subprocess.run(argv, check=True, capture_output=True, text=True, timeout=45).stdout.strip()


def held_identity(body, job):
    if job not in CHANNELS:
        raise ValueError('only the two registered existing flip handles are admitted')
    row = dict(t.split('=', 1) for t in shlex.split(body) if '=' in t)
    channel = CHANNELS[job]
    expected = dict(JobId=job, JobName='hfo2-G2-E060-'+channel, JobState='PENDING',
        Priority='0', Reason='JobHeldUser', Partition='hfacnormal01', NumCPUs='32',
        NumTasks='1', Command=str(PILOT), StdOut=str(WAIT/(channel+'.slurm.out')),
        StdErr=str(WAIT/(channel+'.slurm.err')), TimeLimit='04:00:00',
        Dependency='(null)', **{'CPUs/Task': '32'})
    if (any(row.get(k) != v for k, v in expected.items())
            or not row.get('UserId', '').startswith('iai806(')
            or row.get('NumNodes') not in ('1', '1-1')
            or row.get('Restarts') != '0'):
        raise ValueError('existing held identity/resources/source/dependency changed')
    return row


def committed_slots(body, *, require_spare):
    """Runnable pending jobs reserve a slot too; unrelated pc-* stays untouched."""
    seen, committed = set(), set()
    for line in body.splitlines():
        parts = [p.strip() for p in line.split('|')]
        if len(parts) != 6:
            raise ValueError('malformed native queue')
        job, name, state, partition, cpus, reason = parts
        if not name.startswith('hfo2'):
            continue
        expected_name = ('hfo2-G2-E060-'+CHANNELS[job] if job in CHANNELS
                         else 'hfo2-G2-E091-'+LIVE[job] if job in LIVE else None)
        if (expected_name is None or name != expected_name or job in seen
                or partition != 'hfacnormal01' or cpus != '32'):
            raise ValueError('unexpected/duplicate study allocation')
        seen.add(job)
        if state in ('RUNNING', 'COMPLETING', 'CONFIGURING'):
            committed.add(job)
        elif state == 'PENDING':
            if reason not in ('JobHeldUser', '(JobHeldUser)'):
                committed.add(job)
        else:
            raise ValueError('unexpected study state')
    if len(committed) > (1 if require_spare else 2):
        raise ValueError('no free committed slot within the two-allocation ceiling')
    return sorted(committed)


def load_common():
    if sha(COMMON) != COMMON_SHA:
        raise ValueError('original endpoint/geometry preflight changed')
    spec = importlib.util.spec_from_file_location('E100_original_seed_preflight', COMMON)
    common = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(common)
    return common


def pending_readiness(job):
    if sha(WAIT/'preflight.json') != PREFLIGHT_SHA:
        raise ValueError('original submitted preflight changed')
    proof = json.loads((WAIT/'preflight.json').read_text())
    if (proof['source_root'] != str(SOURCE) or proof['condition'] != 'strain_p0100'
            or len(proof['source_code_sha256']) != 338
            or sha(PILOT) != proof['pilot_sha256']
            or any(sha(SOURCE/n) != h for n, h in proof['source_code_sha256'].items())):
        raise ValueError('original immutable production runtime changed')
    common = load_common()
    rules = common.load_rules()
    if sha(rules.BINARY) != proof['ABACUS_binary_sha256']:
        raise ValueError('original ABACUS binary changed')
    registered = [r for r in proof['reports'] if r['channel'] == CHANNELS[job]]
    if len(registered) != 1:
        raise ValueError('one original registered seed required')
    entry = registered[0]
    seed = Path(entry['seed_root'])
    if (entry['condition'] != 'strain_p0100'
            or sha(seed/'manifest.json') != entry['manifest_sha256']
            or any(sha(seed/n) != h for n, h in entry['seed_files_sha256'].items())
            or (WAIT/CHANNELS[job]/'band').exists()):
        raise ValueError('original pending seed changed or production already exists')
    run(['bash', '-n', str(PILOT)])
    sys.path.insert(0, str(SOURCE))
    original_run = subprocess.run
    try:
        def forbidden(*args, **kwargs):
            raise RuntimeError('DFT/external calls forbidden inside seed audit')
        subprocess.run = forbidden
        pending = common.pending_seed_preflight(seed)
    finally:
        subprocess.run = original_run
    for name in ('vcneb', 'examples.hfo2_fixed_input_factory', 'scripts.audit_hfo2_static_replica'):
        if not Path(sys.modules[name].__file__).resolve().is_relative_to(SOURCE.resolve()):
            raise ValueError('actual seed audit must use original runtime')
    return dict(original_runtime_files_checked=338, original_preflight_sha256=PREFLIGHT_SHA,
        seed_manifest_sha256=entry['manifest_sha256'], seed_audit=pending)


def preflight(job, *, require_spare=False):
    if job not in CHANNELS:
        raise ValueError('only the two registered existing flip handles are admitted')
    accounting = run(['sacct', '-X', '-n', '-P', '-j', job,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'])
    if accounting != job+'|PENDING|0:0|0|hfacnormal01|iai806':
        raise ValueError('existing unstarted pending native handle required')
    control = held_identity(run(['scontrol', 'show', 'job', '-o', job]), job)
    submitted = run(['sacct', '-X', '-n', '-P', '-j', job,
        '--format=JobIDRaw,SubmitLine%4000'])
    values = submitted.split('|')
    if len(values) != 2 or values[0] != job:
        raise ValueError('one original submitted native row required')
    argv = shlex.split(values[1])
    exports = [v[len('--export='):] for v in argv if v.startswith('--export=')]
    expected = dict(RUN_DFT='1', SOURCE_ROOT=str(SOURCE),
        SEED_ROOT=str(RUNS/'20261010-varneb-clamped-G2-remaining-E056-r1/seeds/strain_p0100'/CHANNELS[job]),
        WORKDIR=str(WAIT/CHANNELS[job]/'band'))
    if (not argv or argv[0] != 'sbatch' or argv[-1] != str(PILOT) or len(exports) != 1
            or dict(v.split('=', 1) for v in exports[0].split(',') if '=' in v) != expected):
        raise ValueError('original submitted source/seed/environment changed')
    queue = run(['squeue', '-h', '-u', 'iai806', '-o', '%i|%100j|%T|%P|%C|%R'])
    slots = committed_slots(queue, require_spare=require_spare)
    seed = pending_readiness(job)
    return dict(status='read_only_held_seed_rechecked_not_release_or_convergence',
        job_id=job, check_UTC=datetime.now(timezone.utc).isoformat(),
        native_accounting=accounting, native_control=control, native_queue=queue,
        native_submitted_row=submitted,
        committed_study_slots=slots, **seed, new_DFT_calls=0, scheduler_writes=0,
        new_submissions=0, physical_inputs_or_running_source_changed=False,
        holdout_generated_or_read=False, executed_helper_sha256=sha(Path(__file__)))


def verified_terminal(path, digest, job):
    if sha(path) != digest:
        raise ValueError('reviewed actual terminal audit changed')
    receipt = json.loads(path.read_text())
    observation = receipt.get('terminal_observation', {})
    force = observation.get('replayed_fmax_eV_A')
    expected = dict(job_id=job, ordinary_converged=True,
        status='audited_ordinary_converged_not_TS_or_sampling_certificate',
        new_DFT_calls=0, scheduler_mutations=0,
        physical_inputs_or_running_source_changed=False, holdout_generated_or_read=False,
        all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked=True)
    if (any(receipt.get(k) != v for k, v in expected.items())
            or observation.get('ordinary_residual_pass') is not True
            or type(force) not in (int, float) or not math.isfinite(force) or not 0 <= force <= .10
            or observation.get('pressure_GPa') != 0 or observation.get('climb') is not False
            or observation.get('n_total_images') != 9 or observation.get('n_active_images') != 7
            or sha(Path(receipt['workdir'])/'vcneb_summary.json') != receipt['terminal_summary_sha256']):
        raise ValueError('actual ordinary-converged original terminal material required')
    if run(['sacct', '-X', '-n', '-P', '-j', job,
            '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User']) != job+'|COMPLETED|0:0|32|hfacnormal01|iai806':
        raise ValueError('native successful terminal completion required')


def release_once(job, preserving_audit_sha256):
    if job not in CHANNELS:
        raise ValueError('only the two registered existing flip handles are admitted')
    journal_path = ROOT/(job+'.journal.jsonl')
    receipt_path = ROOT/(job+'.receipt.json')
    lock_path = ROOT/'admission.lock'
    if journal_path.exists() or receipt_path.exists() or lock_path.exists():
        raise FileExistsError('existing intent/receipt; reconcile without release retry')
    verified_terminal(PRESERVING_AUDIT, preserving_audit_sha256, '28727382')
    verified_terminal(REVERSING_AUDIT, REVERSING_SHA, '28692776')
    proof = preflight(job, require_spare=True)
    # Serialize the two handles as well as protecting each individual intent.
    # A failure leaves both lock and journal for read-only reconciliation.
    with lock_path.open('x') as lock:
        lock.write(json.dumps(dict(job_id=job, helper_sha256=sha(Path(__file__)))))
        lock.flush()
        os.fsync(lock.fileno())
    with journal_path.open('x') as journal:
        def record(event, **values):
            journal.write(json.dumps(dict(event=event,
                check_UTC=datetime.now(timezone.utc).isoformat(), **values))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        record('verified_original_held_target', proof=proof,
            preserving_audit_sha256=preserving_audit_sha256, reversing_audit_sha256=REVERSING_SHA)
        try:
            committed_slots(run(['squeue', '-h', '-u', 'iai806', '-o', '%i|%100j|%T|%P|%C|%R']), require_spare=True)
            held_identity(run(['scontrol', 'show', 'job', '-o', job]), job)
            record('one_existing_release_intent', job_id=job)
            reply = run(['scontrol', 'release', job])
            record('one_existing_release_accepted', job_id=job, reply=reply)
            # Keep the shared lock until an independent native read sees the
            # released handle reserving its slot, not merely an empty reply.
            row = dict(t.split('=', 1) for t in shlex.split(
                run(['scontrol', 'show', 'job', '-o', job])) if '=' in t)
            if (row.get('JobId') != job or row.get('Command') != str(PILOT)
                    or row.get('Partition') != 'hfacnormal01' or row.get('NumCPUs') != '32'
                    or row.get('JobState') not in ('PENDING', 'RUNNING', 'COMPLETING', 'CONFIGURING')
                    or row.get('Reason') == 'JobHeldUser' or row.get('Priority') == '0'):
                raise ValueError('accepted release not independently corroborated; retain lock')
            slots = committed_slots(run(['squeue', '-h', '-u', 'iai806', '-o', '%i|%100j|%T|%P|%C|%R']), require_spare=False)
            if job not in slots:
                raise ValueError('released target not yet a committed native slot; retain lock')
            record('released_slot_independently_corroborated', job_id=job,
                native_control=row, committed_slots=slots)
        except Exception as exc:
            record('reconcile_never_repeat_write', error=type(exc).__name__, message=str(exc))
            raise
    result = dict(status='one_existing_flip_handle_released_not_start_proof',
        job_id=job, channel=CHANNELS[job], new_submissions=0, dependency_writes=0,
        release_writes=1, new_independent_chains=0, admission_DFT_calls=0,
        physical_inputs_or_running_source_changed=False, holdout_generated_or_read=False,
        preserving_audit_sha256=preserving_audit_sha256, reversing_audit_sha256=REVERSING_SHA,
        max_committed_study_slots=2, executed_helper_sha256=sha(Path(__file__)),
        journal_sha256=sha(journal_path))
    with receipt_path.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    lock_path.unlink()  # Accepted intent/result remain in the immutable journal.
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job-id', choices=tuple(CHANNELS), required=True)
    parser.add_argument('--release', action='store_true')
    parser.add_argument('--preserving-audit-sha256')
    args = parser.parse_args()
    if args.release and not args.preserving_audit_sha256:
        parser.error('release requires a reviewed actual ordinary terminal audit SHA256')
    result = release_once(args.job_id, args.preserving_audit_sha256) if args.release else preflight(args.job_id)
    print(json.dumps(result, indent=2, allow_nan=False))
