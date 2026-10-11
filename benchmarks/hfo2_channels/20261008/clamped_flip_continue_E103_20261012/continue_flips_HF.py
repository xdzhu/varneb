"""Finite same-chain continuation of the two original +1% G2 flip starters.

Audit-only replay uses existing pinned E093 support and the original runtime.
Dispatch requires a successful native cap, all SCFs and nine latest caches.
It journals one held submission/release, reserves runnable pending jobs, and
shares E100's study admission lock. An ambiguous intent is never retried.
"""
import argparse
import contextlib
from datetime import datetime, timezone
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS/'20261012-varneb-clamped-flip-continue-E103-r1'
ORIGINAL = RUNS/'20261010-varneb-clamped-G2-E060-r1'
SOURCE = RUNS/'20261010-varneb-nested-controls-E055-r1/source'
PILOT = SOURCE/'cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm'
PILOT_SHA = '7bf764238c4cc857ca2b7e3e94642f0517e1f3b84f4e4c5b17833379ac1101ed'
SCRIPT = ROOT/'hf_hfo2_clamped_continue_E103_20261012.slurm'
SCRIPT_SHA = 'ab914119aaa7918bd0cb466ed8dd86daafecb20609a077a5d249d89a3228eed0'
LEGACY = RUNS/'20261011-varneb-clamped-child-audit-E093-r1/audit_child_HF.py'
LEGACY_SHA = '550bf525b6581a76c9f7634b5d3d8dcb75258193ad6189a57c2b6a66bf2e4ef3'
LOCK = RUNS/'20261011-varneb-clamped-held-flips-E100-r1/admission.lock'
BINARY = Path('/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus')
BINARY_SHA = '88ea7f91c9f3091ef5fb4f410ab88267b14ecfe7f35c86dfaee9bf9319cace7a'
CHANNELS = {'28709790': 'PO_flip_T_pattern_preserving',
            '28709791': 'PO_flip_T_pattern_reversing'}
EVENTS = ['one_held_submission_intent', 'one_accepted_handle',
          'one_release_intent', 'one_existing_handle_released',
          'released_slot_independently_corroborated']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(argv):
    return subprocess.run(argv, check=True, capture_output=True, text=True, timeout=45).stdout.strip()


def require_parent(parent):
    if parent not in CHANNELS:
        raise ValueError('only the two original +1% flip parents are admitted')


def paths(parent):
    require_parent(parent)
    channel = CHANNELS[parent]
    return dict(old_work=ORIGINAL/channel/'band', terminal=ROOT/parent/'terminal',
        work=ROOT/parent/channel/'band', stdout=ROOT/parent/(channel+'.slurm.out'),
        stderr=ROOT/parent/(channel+'.slurm.err'))


def environment(parent, receipt, digest):
    p = paths(parent)
    return dict(RUN_DFT='1', SOURCE_ROOT=str(SOURCE),
        SEED_ROOT=str(p['terminal']/'seed'), WORKDIR=str(p['work']),
        PARENT_JOB_ID=parent, PARENT_AUDIT_PATH=str(p['terminal']/'audit_receipt.json'),
        PARENT_AUDIT_SHA256=digest, SEED_MANIFEST_SHA256=receipt['seed_manifest_sha256'])


def identity(body, job, parent, *, child=False):
    """Native submitted ancestry is checked before terminal files/imports."""
    p = paths(parent)
    fields = body.strip().split('|')
    if len(fields) != 9:
        raise ValueError('one complete native accounting source row required')
    actual, state, user, cpus, partition, name, out, err, submit = fields
    expected_name = 'hfo2-G2-E103-'+parent if child else 'hfo2-G2-E060-'+CHANNELS[parent]
    root = ROOT/parent if child else ORIGINAL
    if (actual != job or state != 'COMPLETED' or user != 'iai806' or cpus != '32'
            or partition != 'hfacnormal01' or name != expected_name
            or out != str(root/(CHANNELS[parent]+'.slurm.out'))
            or err != str(root/(CHANNELS[parent]+'.slurm.err'))):
        raise ValueError('exact native successful terminal identity required')
    script = SCRIPT if child else PILOT
    argv = shlex.split(submit)
    exports = [v[len('--export='):] for v in argv if v.startswith('--export=')]
    if not argv or argv[0] != 'sbatch' or argv[-1] != str(script) or len(exports) != 1:
        raise ValueError('original native submitted wrapper required')
    terms = exports[0].split(',')
    pairs = [v.split('=', 1) for v in terms if '=' in v]
    if len({k for k, _ in pairs}) != len(pairs) or [v for v in terms if '=' not in v] != ['ALL']:
        raise ValueError('unambiguous original exported environment required')
    expected = (environment(parent, json.loads((p['terminal']/'audit_receipt.json').read_text()),
                            sha(p['terminal']/'audit_receipt.json')) if child else
        dict(RUN_DFT='1', SOURCE_ROOT=str(SOURCE),
            SEED_ROOT=str(RUNS/'20261010-varneb-clamped-G2-remaining-E056-r1/seeds/strain_p0100'/CHANNELS[parent]),
            WORKDIR=str(p['old_work'])))
    if dict(pairs) != expected:
        raise ValueError('original submitted ancestry/physical environment differs')


def load_auditor():
    if sha(LEGACY) != LEGACY_SHA or sha(PILOT) != PILOT_SHA:
        raise ValueError('pinned audit adapter or immutable original wrapper changed')
    spec = importlib.util.spec_from_file_location('E103_pinned_E093', LEGACY)
    legacy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)
    registration = legacy.REGISTRY['28727368']  # Same original 338-file E055 runtime.
    auditor = legacy.load_auditor(registration)
    return auditor, registration, legacy.SUPPORT_SHA


def audit(parent, *, child=False):
    p = paths(parent)
    admission = accepted(parent) if child else None
    job = admission['job_id'] if child else parent
    output = ROOT/parent/('child_terminal' if child else 'terminal')
    if output.exists():
        raise FileExistsError('inspect existing audit; never overwrite/replay it')
    body = run(['sacct', '-X', '-n', '-P', '-j', job,
        '--format=JobIDRaw,State,User,AllocCPUS,Partition,JobName%100,StdOut%1500,StdErr%1500,SubmitLine%4000'])
    identity(body, job, parent, child=child)
    auditor, registration, support = load_auditor()
    script, work, cap = (SCRIPT, p['work'], 40) if child else (PILOT, p['old_work'], 10)
    if child:
        if sha(SCRIPT) != SCRIPT_SHA:
            raise ValueError('accepted continuation wrapper changed')
        # Register only this pinned reviewed wrapper in this audit-only process.
        from scripts.export_hfo2_clamped_observation import REGISTERED_FACTORIES
        REGISTERED_FACTORIES[SCRIPT_SHA] = 'examples.hfo2_fixed_input_factory:make_clamped_resume_cached_factory'
    summary = json.loads((work/'vcneb_summary.json').read_text())
    step = auditor.completed_step(work/'vcneb.opt.log')
    capped = summary.get('status') == 'max_steps_reached'
    if not 0 <= step <= cap or (capped and step != cap):
        raise ValueError('unexpected finite history; inspect instead of silently continuing')
    receipt = auditor.audit(work, job, script, SOURCE, registration['archive'],
        registration['archive_sha'], output, outcome='step-cap' if capped else 'converged',
        seed_latest=capped)  # Dispatch admits only original starters, not a rolling child.
    if len(receipt['runtime_code_sha256']) != 338:
        raise ValueError('original immutable production runtime count differs')
    with (output/'driver_receipt.json').open('x') as stream:
        json.dump(dict(job_id=job, original_parent=parent, native_source_row=body,
            executed_driver_sha256=sha(Path(__file__)), support_sha256=support,
            terminal_audit_receipt_sha256=sha(output/'audit_receipt.json'),
            new_DFT_calls=0, scheduler_writes=0, new_independent_chains=0), stream, indent=2)
    return receipt


def accepted(parent):
    require_parent(parent)
    root = ROOT/parent
    result = json.loads((root/'dispatch_receipt.json').read_text())
    journal = root/'dispatch_journal.jsonl'
    records = [json.loads(line) for line in journal.read_text().splitlines()]
    if (result.get('status') != 'one_same_chain_continuation_submitted_released_and_slot_corroborated'
            or result.get('parent_job') != parent or not str(result.get('job_id', '')).isdecimal()
            or result.get('total_submissions') != 1 or result.get('release_writes') != 1
            or result.get('new_independent_chains') != 0 or result.get('script_sha256') != SCRIPT_SHA
            or result.get('executed_helper_sha256') != sha(Path(__file__))
            or result.get('journal_sha256') != sha(journal)
            or [r['event'] for r in records] != EVENTS
            or any(r.get('job_id') != result['job_id'] for r in records[1:])):
        raise ValueError('accepted admission evidence differs; reconcile before any write')
    return result


def slots(body, parent, *, spare, extra=None):
    require_parent(parent)
    names = {j: 'hfo2-G2-E060-'+c for j, c in CHANNELS.items()}
    names.update({'28727368': 'hfo2-G2-E091-28709789', '28727382': 'hfo2-G2-E091-28723655'})
    for original in CHANNELS:
        if (ROOT/original/'dispatch_receipt.json').exists():
            result = accepted(original)
            names[result['job_id']] = 'hfo2-G2-E103-'+original
    names.update(extra or {})
    seen, committed = set(), set()
    for line in body.splitlines():
        values = [v.strip() for v in line.split('|')]
        if len(values) != 6:
            raise ValueError('malformed native queue')
        job, name, state, partition, cpus, reason = values
        if not name.startswith('hfo2'):
            continue
        if (job in seen or names.get(job) != name or partition != 'hfacnormal01' or cpus != '32'):
            raise ValueError('unexpected/duplicate study allocation')
        seen.add(job)
        if state in ('RUNNING', 'COMPLETING', 'CONFIGURING'):
            committed.add(job)
        elif state == 'PENDING':
            if reason not in ('JobHeldUser', '(JobHeldUser)'):
                committed.add(job)
        else:
            raise ValueError('unexpected study state')
    if parent in committed or len(committed) > (1 if spare else 2):
        raise ValueError('terminal parent and a free committed study slot required')
    return sorted(committed)


def verify_audit(parent, reviewed_sha):
    p = paths(parent)
    path = p['terminal']/'audit_receipt.json'
    if sha(path) != reviewed_sha or sha(SCRIPT) != SCRIPT_SHA or sha(BINARY) != BINARY_SHA:
        raise ValueError('reviewed receipt, wrapper or binary differs')
    receipt = json.loads(path.read_text())
    if receipt['terminal_summary_sha256'] != sha(p['old_work']/'vcneb_summary.json') or p['work'].exists():
        raise ValueError('terminal summary or fresh work namespace differs')
    _, _, support = load_auditor()
    for key, name in (('executed_auditor_sha256', 'audit_hfo2_clamped_terminal.py'),
                      ('executed_exporter_sha256', 'export_hfo2_clamped_observation.py'),
                      ('executed_residual_analyzer_sha256', 'analyze_hfo2_clamped_residual.py')):
        if receipt[key] != support[name]:
            raise ValueError('actual terminal audit support differs')
    body = SCRIPT.read_text().split(" 'import hashlib,json,pathlib,sys\n", 1)[1]
    body = 'import hashlib,json,pathlib,sys\n'+body.split("' \"${seed_root}/manifest.json\"", 1)[0]
    old = sys.argv
    try:
        sys.argv = ['reviewed-startup', str(p['terminal']/'seed/manifest.json'),
            receipt['seed_manifest_sha256'], str(path), reviewed_sha, parent, str(SOURCE)]
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(body, str(SCRIPT), 'exec'), {})
    finally:
        sys.argv = old
    return receipt


def dispatch(parent, reviewed_sha):
    p = paths(parent)
    root = ROOT/parent
    journal_path = root/'dispatch_journal.jsonl'
    if journal_path.exists() or (root/'dispatch_receipt.json').exists() or LOCK.exists():
        raise FileExistsError('existing intent/lock; reconcile without write retries')
    receipt = verify_audit(parent, reviewed_sha)
    accounting = run(['sacct', '-X', '-n', '-P', '-j', parent,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'])
    if accounting != parent+'|COMPLETED|0:0|32|hfacnormal01|iai806':
        raise ValueError('audited parent must remain natively terminal')
    queue = ['squeue', '-h', '-u', 'iai806', '-o', '%i|%100j|%T|%P|%C|%R']
    slots(run(queue), parent, spare=True)
    run(['bash', '-n', str(SCRIPT)])
    name = 'hfo2-G2-E103-'+parent
    argv = ['sbatch', '--parsable', '--hold', '--job-name='+name,
        '--export=ALL,'+','.join(k+'='+v for k, v in environment(parent, receipt, reviewed_sha).items()),
        '--output='+str(p['stdout']), '--error='+str(p['stderr']), str(SCRIPT)]
    with LOCK.open('x') as stream:
        stream.write(json.dumps(dict(parent_job=parent, helper_sha256=sha(Path(__file__)))))
        stream.flush()
        os.fsync(stream.fileno())
    with journal_path.open('x') as journal:
        def record(event, **data):
            journal.write(json.dumps(dict(event=event, check_UTC=datetime.now(timezone.utc).isoformat(), **data))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        try:
            active = slots(run(queue), parent, spare=True)
            record(EVENTS[0], argv=argv, parent_accounting=accounting,
                audit_receipt_sha256=reviewed_sha, committed_slots=active)
            handle = run(argv)
            if not handle.isdecimal():
                raise ValueError('ambiguous submission reply; never resubmit')
            record(EVENTS[1], job_id=handle)
            row = dict(t.split('=', 1) for t in shlex.split(run(['scontrol', 'show', 'job', '-o', handle])) if '=' in t)
            expected = dict(JobId=handle, JobName=name, JobState='PENDING', Reason='JobHeldUser',
                Priority='0', Partition='hfacnormal01', NumCPUs='32', NumTasks='1',
                Command=str(SCRIPT), Dependency='(null)', StdOut=str(p['stdout']), StdErr=str(p['stderr']),
                TimeLimit='12:00:00', Restarts='0', **{'CPUs/Task': '32'})
            if (any(row.get(k) != v for k, v in expected.items())
                    or row.get('NumNodes') not in ('1', '1-1') or not row.get('UserId', '').startswith('iai806(')):
                raise ValueError('exact accepted held continuation not corroborated')
            slots(run(queue), parent, spare=True, extra={handle: name})
            record(EVENTS[2], job_id=handle, accepted_native_state=row)
            reply = run(['scontrol', 'release', handle])
            record(EVENTS[3], job_id=handle, reply=reply)
            after = dict(t.split('=', 1) for t in shlex.split(run(['scontrol', 'show', 'job', '-o', handle])) if '=' in t)
            if (any(after.get(k) != v for k, v in expected.items() if k not in ('JobState', 'Reason', 'Priority'))
                    or after.get('JobState') not in ('PENDING', 'RUNNING', 'COMPLETING', 'CONFIGURING')
                    or after.get('Reason') == 'JobHeldUser' or after.get('Priority') == '0'):
                raise ValueError('release not independently corroborated; retain shared lock')
            committed = slots(run(queue), parent, spare=False, extra={handle: name})
            if handle not in committed:
                raise ValueError('released target not a committed slot; retain shared lock')
            record(EVENTS[4], job_id=handle, native_control=after, committed_slots=committed)
        except Exception as exc:
            record('reconcile_never_repeat_write', error=type(exc).__name__, message=str(exc))
            raise
    result = dict(status='one_same_chain_continuation_submitted_released_and_slot_corroborated',
        job_id=handle, parent_job=parent, source_step=10, max_new_updates=40, wall_hours=12,
        total_submissions=1, release_writes=1, dependency_writes=0, new_independent_chains=0,
        audit_receipt_sha256=reviewed_sha, seed_manifest_sha256=receipt['seed_manifest_sha256'],
        script_sha256=SCRIPT_SHA, executed_helper_sha256=sha(Path(__file__)), journal_sha256=sha(journal_path),
        physical_inputs_or_running_source_changed=False, holdout_generated_or_read=False,
        FIRE_optimizer_state_restored=False, execution_started_claimed=False)
    with (root/'dispatch_receipt.json').open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    LOCK.unlink()  # Durable accepted intent/receipt are retained, never reset.
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('audit-parent', 'dispatch', 'audit-child'))
    parser.add_argument('--parent', choices=tuple(CHANNELS), required=True)
    parser.add_argument('--audit-sha256')
    args = parser.parse_args()
    if args.action == 'dispatch' and not args.audit_sha256:
        parser.error('dispatch requires the reviewed actual terminal receipt SHA256')
    result = dispatch(args.parent, args.audit_sha256) if args.action == 'dispatch' else audit(args.parent, child=args.action == 'audit-child')
    print(json.dumps(result, indent=2, allow_nan=False))
