"""One audited +1% T-to-M rolling handoff; read-only preflight is separate.

No sbatch, dependency update, retry, production-source or physical-input edit.
Only the registered M waiter may be released after an exclusive intent record.
"""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path

RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS/'20261011-varneb-clamped-plus1-T-terminal-E086-r1'
SHARED = RUNS/'20261011-varneb-clamped-reversing-terminal-E083-r1/release_existing_plus1_once.py'
SHARED_SHA = '8c6a6b24ffc1162ce3d860f4fc31708ea37b56f4b8b1d0895d0493cb6e4b955d'
WAIT = RUNS/'20261010-varneb-clamped-G2-E060-r1'
PREFLIGHT_SHA = 'b780585fbd5ccd05e85282f92dc746a5a79d48f6b25d923ff0bf3150098a2e8a'
SUBSTRATE_SHA = 'f1b88a625b851149201149faa91f5d4d5565b1e77022a0419974272d6c9cfa74'
DONE, LIVE, TARGET = '28709788', '28723655', '28709789'
WAITERS = {'28709789': 'PO_to_M', '28709790': 'PO_flip_T_pattern_preserving',
           '28709791': 'PO_flip_T_pattern_reversing'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_common():
    if sha(SHARED) != SHARED_SHA:
        raise ValueError('reviewed shared preflight changed')
    spec = importlib.util.spec_from_file_location('E083_read_only_shared', SHARED)
    shared = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(shared)
    if sha(shared.COMMON) != shared.COMMON_SHA:
        raise ValueError('reviewed endpoint preflight changed')
    return shared.load_common()


def verified_parent(receipt):
    observation = receipt.get('terminal_observation', {})
    force, step = observation.get('replayed_fmax_eV_A'), receipt.get('terminal_step')
    expected = dict(job_id=DONE, workdir=str(WAIT/'PO_to_T/band'), ordinary_converged=True,
        status='audited_ordinary_converged_not_TS_or_sampling_certificate',
        new_DFT_calls=0, scheduler_mutations=0, physical_inputs_or_running_source_changed=False,
        holdout_generated_or_read=False,
        all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked=True,
        initial_fresh_interior_SCFs=7)
    boundary = observation.get('mechanical_boundary', {})
    if (any(receipt.get(k) != v for k, v in expected.items())
            or observation.get('source_job_id') != DONE
            or observation.get('ordinary_residual_pass') is not True
            or isinstance(force, bool) or not isinstance(force, (int, float))
            or not math.isfinite(force) or not 0 <= force <= .10
            or type(step) is not int or step < 0
            or type(receipt.get('fresh_interior_SCFs')) is not int
            or receipt['fresh_interior_SCFs'] != 7*(step+1)
            or boundary.get('reference_file_sha256') != SUBSTRATE_SHA
            or boundary.get('kind') != 'clamped_plane'
            or boundary.get('allow_tilt') is not True
            or boundary.get('fixed_cell_rows') != [0, 1] or boundary.get('cell_dofs') != 3):
        raise ValueError('actual ordinary-converged registered +1% T audit required')


def held_waiter(body, job, rules):
    if job not in WAITERS:
        raise ValueError('unregistered remaining waiter')
    row = rules.fields(body)
    channel = WAITERS[job]
    dep = rules.dependencies(row.get('Dependency', ''))
    # Slurm can remove the now-completed T dependency. The still-held M
    # dependency must remain on both flipping starters. Never rewrite either.
    required = set() if job == TARGET else {TARGET}
    allowed = set() if job == TARGET else {DONE, TARGET}
    expected = dict(JobId=job, JobName='hfo2-G2-E060-'+channel, JobState='PENDING',
        Priority='0', Reason='JobHeldUser', Partition='hfacnormal01', NumCPUs='32',
        NumTasks='1', Command=rules.PILOT, TimeLimit='04:00:00',
        StdOut=str(WAIT/(channel+'.slurm.out')), StdErr=str(WAIT/(channel+'.slurm.err')))
    if (any(row.get(k) != v for k, v in expected.items())
            or not row.get('UserId', '').startswith('iai806(')
            or row.get('CPUs/Task') != '32' or not required.issubset(dep) or not dep.issubset(allowed)):
        raise ValueError('held identity/resources/original dependency differs')
    return row


def spare_slot(body, *, live_terminal):
    seen, active = set(), set()
    for line in body.splitlines():
        parts = [p.strip() for p in line.split('|')]
        if len(parts) != 5:
            raise ValueError('malformed native queue')
        job, name, state, partition, cpus = parts
        if not name.startswith('hfo2'):
            continue
        if (job not in {LIVE, *WAITERS} or job in seen
                or partition != 'hfacnormal01' or cpus != '32'):
            raise ValueError('unexpected registered-study allocation')
        seen.add(job)
        if state in {'RUNNING', 'COMPLETING', 'CONFIGURING'}:
            active.add(job)
        elif state != 'PENDING':
            raise ValueError('unexpected native study state')
    if active != (set() if live_terminal else {LIVE}) or not set(WAITERS).issubset(seen):
        raise ValueError('actual registered free slot and all three protected waiters required')
    return sorted(active)


def pending_inputs():
    common = load_common()
    rules = common.load_rules()
    if sha(WAIT/'preflight.json') != PREFLIGHT_SHA:
        raise ValueError('registered +1% preparation changed')
    proof = json.loads((WAIT/'preflight.json').read_text())
    source = Path(proof['source_root'])
    if (proof['condition'] != 'strain_p0100'
            or sha(Path(rules.PILOT)) != proof['pilot_sha256']
            or sha(rules.BINARY) != proof['ABACUS_binary_sha256']
            or not proof['source_code_sha256']
            or any(sha(source/n) != h for n, h in proof['source_code_sha256'].items())):
        raise ValueError('pending original runtime, binary or condition changed')
    entries = [r for r in proof['reports'] if r['channel'] == WAITERS[TARGET]]
    if len(entries) != 1:
        raise ValueError('one registered M seed required')
    entry = entries[0]
    seed = Path(entry['seed_root'])
    if (entry['condition'] != 'strain_p0100'
            or entry['seed_files_sha256'].get('substrate.vasp') != SUBSTRATE_SHA
            or sha(seed/'manifest.json') != entry['manifest_sha256']
            or any(sha(seed/n) != h for n, h in entry['seed_files_sha256'].items())
            or (WAIT/WAITERS[TARGET]/'band').exists()):
        raise ValueError('pending registered seed changed or production already exists')
    pending = common.pending_seed_preflight(seed)
    rows = {job: held_waiter(rules.run(['scontrol', 'show', 'job', '-o', job]), job, rules)
            for job in WAITERS}
    return rules, dict(target_job=TARGET, condition='strain_p0100',
        source_python_files_checked=len(proof['source_code_sha256']),
        pending_seed_preflight=pending, held_rows=rows,
        new_DFT_calls=0, scheduler_mutations=0, admission_or_parent_convergence_claimed=False)


def native_spare_slot(rules):
    fmt = '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'
    parent = rules.run(['sacct', '-X', '-n', '-P', '-j', DONE, fmt])
    if parent != DONE+'|COMPLETED|0:0|32|hfacnormal01|iai806':
        raise ValueError('native actual parent completion required')
    live = rules.run(['sacct', '-X', '-n', '-P', '-j', LIVE, fmt])
    if live not in {LIVE+'|'+s+'|0:0|32|hfacnormal01|iai806'
                   for s in ('RUNNING', 'COMPLETING', 'COMPLETED')}:
        raise ValueError('native registered companion accounting required')
    active = spare_slot(rules.run(['squeue', '-h', '-u', 'iai806', '-o', '%i|%100j|%T|%P|%C']),
                        live_terminal=live.split('|')[1] == 'COMPLETED')
    rows = {job: held_waiter(rules.run(['scontrol', 'show', 'job', '-o', job]), job, rules)
            for job in WAITERS}
    return dict(active_jobs=active, parent_accounting=parent, companion_accounting=live, held_rows=rows)


def release_once(audit_sha256):
    journal_path, audit_path = ROOT/'admission_journal.jsonl', ROOT/'audit_receipt.json'
    if journal_path.exists():
        raise FileExistsError('reconcile existing intent; never repeat release')
    if sha(audit_path) != audit_sha256:
        raise ValueError('reviewed actual terminal audit changed')
    receipt = json.loads(audit_path.read_text())
    verified_parent(receipt)
    if sha(Path(receipt['workdir'])/'vcneb_summary.json') != receipt['terminal_summary_sha256']:
        raise ValueError('terminal summary changed')
    rules, pending = pending_inputs()
    native = native_spare_slot(rules)
    with journal_path.open('x') as journal:
        def record(event, **values):
            journal.write(json.dumps(dict(event=event,
                checked_UTC=datetime.now(timezone.utc).isoformat(), **values))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        record('admission_verified', target=TARGET, parent_audit_sha256=audit_sha256,
               pending_preflight=pending, native_snapshot=native, executed_helper_sha256=sha(Path(__file__)))
        try:
            checked = native_spare_slot(rules)
            if sha(audit_path) != audit_sha256 or sha(Path(receipt['workdir'])/'vcneb_summary.json') != receipt['terminal_summary_sha256']:
                raise ValueError('terminal evidence changed before dispatch')
            record('release_intent', job_id=TARGET, rechecked_native=checked)
            rules.run(['scontrol', 'release', TARGET])
            record('release_accepted', job_id=TARGET)
        except Exception as exc:
            record('reconcile_do_not_repeat', error=type(exc).__name__, message=str(exc))
            raise
    result = dict(status='one_existing_plus1_M_starter_released', target_job=TARGET,
        audited_completed_parent=DONE, companion_job=LIVE, startup_claimed=False,
        release_writes=1, new_submissions=0, dependency_writes=0, new_independent_chains=0,
        other_two_starters_still_held=[j for j in WAITERS if j != TARGET],
        new_DFT_calls_in_admission=0, max_study_simultaneous_allocations=2,
        physical_inputs_or_running_source_changed=False, holdout_generated_or_read=False,
        journal_sha256=sha(journal_path), executed_helper_sha256=sha(Path(__file__)))
    with (ROOT/'admission_receipt.json').open('x') as stream:
        json.dump(result, stream, indent=2)
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--preflight-only', type=Path, metavar='NEW_REPORT')
    action.add_argument('--reviewed-audit-sha256')
    args = parser.parse_args()
    if args.preflight_only is not None:
        if args.preflight_only.exists():
            raise FileExistsError('refusing an existing preflight report')
        _, result = pending_inputs()
        with args.preflight_only.open('x') as stream:
            json.dump(result, stream, indent=2)
    else:
        result = release_once(args.reviewed_audit_sha256)
    print(json.dumps(result))
