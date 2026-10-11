"""Admit one existing +1% training starter after native terminal audit.

Exclusive journal before the sole scheduler write. Never submit a job,
rewrite a dependency, release the other three starters or replay an intent.
"""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path

RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS/'20261011-varneb-clamped-reversing-terminal-E083-r1'
COMMON = RUNS/'20261011-varneb-clamped-M-terminal-E068-source-r2/release_existing_PO_to_T.py'
COMMON_SHA = '5ddca4a3bf73ed8ad1d8f348544d97bf170ef442354a63a1187ae49a5c56c0f0'
DONE, PREVIOUS, LIVE, TARGET = '28692776', '28692775', '28723655', '28709788'
WAIT = RUNS/'20261010-varneb-clamped-G2-E060-r1'
PREFLIGHT_SHA = 'b780585fbd5ccd05e85282f92dc746a5a79d48f6b25d923ff0bf3150098a2e8a'
T_AUDIT = RUNS/'20261011-varneb-clamped-T-terminal-E075-r1/audit_receipt.json'
T_AUDIT_SHA = 'c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e'
CHANNELS = {'28709788': 'PO_to_T', '28709789': 'PO_to_M',
    '28709790': 'PO_flip_T_pattern_preserving', '28709791': 'PO_flip_T_pattern_reversing'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verified_audit(receipt):
    force = receipt.get('terminal_observation', {}).get('replayed_fmax_eV_A')
    step = receipt.get('terminal_step')
    expected = dict(job_id=DONE, ordinary_converged=True,
        status='audited_ordinary_converged_not_TS_or_sampling_certificate',
        new_DFT_calls=0, scheduler_mutations=0, physical_inputs_or_running_source_changed=False,
        holdout_generated_or_read=False,
        all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked=True,
        initial_fresh_interior_SCFs=7)
    if (any(receipt.get(k) != v for k, v in expected.items())
            or receipt.get('terminal_observation', {}).get('ordinary_residual_pass') is not True
            or isinstance(force, bool) or not isinstance(force, (int, float))
            or not math.isfinite(force) or not 0 <= force <= .10
            or type(step) is not int or step < 0
            or receipt.get('fresh_interior_SCFs') != 7*(step+1)):
        raise ValueError('actual completed ordinary fresh-starter audit required')


def one_used_slot(body):
    active, seen = set(), set()
    for line in body.splitlines():
        parts = [value.strip() for value in line.split('|')]
        if len(parts) != 5:
            raise ValueError('malformed native queue')
        job, name, state, partition, cpus = parts
        if not name.startswith('hfo2'):
            continue
        if (job not in {LIVE, *CHANNELS} or job in seen
                or partition != 'hfacnormal01' or cpus != '32'):
            raise ValueError('unexpected study allocation')
        seen.add(job)
        if state in {'RUNNING', 'COMPLETING', 'CONFIGURING'}:
            active.add(job)
        elif state != 'PENDING':
            raise ValueError('unexpected study state')
    if active != {LIVE} or not set(CHANNELS).issubset(seen):
        raise ValueError('one live continuation and four protected waiters required')
    return sorted(active)


def held_training(body, job, rules):
    if job not in CHANNELS:
        raise ValueError('unregistered training starter')
    row = rules.fields(body)
    channel = CHANNELS[job]
    dep = rules.dependencies(row.get('Dependency', ''))
    allowed = {DONE, PREVIOUS} if job in {'28709788', '28709789'} else {'28709788', '28709789'}
    required = set() if job in {'28709788', '28709789'} else allowed
    expected = dict(JobId=job, JobName='hfo2-G2-E060-'+channel,
        JobState='PENDING', Priority='0', Reason='JobHeldUser',
        Partition='hfacnormal01', NumCPUs='32', NumTasks='1',
        Command=rules.PILOT, StdOut=str(WAIT/(channel+'.slurm.out')),
        StdErr=str(WAIT/(channel+'.slurm.err')), TimeLimit='04:00:00')
    if (any(row.get(k) != v for k, v in expected.items())
            or not row.get('UserId', '').startswith('iai806(')
            or row.get('CPUs/Task') != '32'
            or not required.issubset(dep) or not dep.issubset(allowed)):
        raise ValueError('exact held starter identity/resources/dependencies differ')
    return row


def load_common():
    spec = importlib.util.spec_from_file_location('E068_shared_preflight', COMMON)
    common = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(common)
    return common


def release_once(audit_sha256):
    journal_path = ROOT/'admission_journal.jsonl'
    if journal_path.exists():
        raise FileExistsError('reconcile existing intent; never repeat release')
    if (sha(COMMON) != COMMON_SHA or sha(ROOT/'audit_receipt.json') != audit_sha256
            or sha(T_AUDIT) != T_AUDIT_SHA):
        raise ValueError('reviewed audit or pinned shared preflight changed')
    receipt = json.loads((ROOT/'audit_receipt.json').read_text())
    verified_audit(receipt)
    if sha(Path(receipt['workdir'])/'vcneb_summary.json') != receipt['terminal_summary_sha256']:
        raise ValueError('terminal summary changed')
    common = load_common()
    rules = common.load_rules()
    for job in (DONE, PREVIOUS):
        if rules.run(['sacct', '-X', '-n', '-P', '-j', job,
                '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User']) != job+'|COMPLETED|0:0|32|hfacnormal01|iai806':
            raise ValueError('both actual completed starters required')
    if rules.run(['sacct', '-X', '-n', '-P', '-j', LIVE,
            '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User']) != LIVE+'|RUNNING|0:0|32|hfacnormal01|iai806':
        raise ValueError('native registered live continuation required')
    if sha(WAIT/'preflight.json') != PREFLIGHT_SHA:
        raise ValueError('registered +1% preflight changed')
    proof = json.loads((WAIT/'preflight.json').read_text())
    source = Path(proof['source_root'])
    if (proof['condition'] != 'strain_p0100'
            or sha(Path(rules.PILOT)) != proof['pilot_sha256']
            or sha(rules.BINARY) != proof['ABACUS_binary_sha256']
            or any(sha(source/n) != h for n, h in proof['source_code_sha256'].items())):
        raise ValueError('pending immutable runtime, binary or training condition changed')
    entries = [r for r in proof['reports'] if r['channel'] == CHANNELS[TARGET]]
    if len(entries) != 1:
        raise ValueError('one registered target seed required')
    registered = entries[0]
    seed = Path(registered['seed_root'])
    if (registered['condition'] != 'strain_p0100'
            or sha(seed/'manifest.json') != registered['manifest_sha256']
            or any(sha(seed/n) != h for n, h in registered['seed_files_sha256'].items())
            or (WAIT/CHANNELS[TARGET]/'band').exists()):
        raise ValueError('pending seed changed or production already exists')
    pending = common.pending_seed_preflight(seed)
    rows = {job: held_training(rules.run(['scontrol', 'show', 'job', '-o', job]), job, rules)
            for job in CHANNELS}
    active = one_used_slot(rules.run(['squeue', '-h', '-u', 'iai806', '-o', '%i|%100j|%T|%P|%C']))
    with journal_path.open('x') as journal:
        def record(event, **values):
            journal.write(json.dumps(dict(event=event,
                checked_UTC=datetime.now(timezone.utc).isoformat(), **values))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        record('admission_verified', target=TARGET, parent_audit_sha256=audit_sha256,
            previous_T_audit_sha256=T_AUDIT_SHA, active_jobs=active,
            held_rows=rows, pending_seed_preflight=pending,
            executed_helper_sha256=sha(Path(__file__)))
        try:
            one_used_slot(rules.run(['squeue', '-h', '-u', 'iai806', '-o', '%i|%100j|%T|%P|%C']))
            for job in CHANNELS:
                held_training(rules.run(['scontrol', 'show', 'job', '-o', job]), job, rules)
            record('release_intent', job_id=TARGET)
            rules.run(['scontrol', 'release', TARGET])
            record('release_accepted', job_id=TARGET)
        except Exception as exc:
            record('reconcile_do_not_repeat', error=type(exc).__name__, message=str(exc))
            raise
    result = dict(status='one_existing_plus1_training_starter_released',
        target_job=TARGET, audited_completed_parent=DONE, live_continuation=LIVE,
        new_submissions=0, dependency_writes=0, release_writes=1,
        other_three_starters_still_held=[job for job in CHANNELS if job != TARGET],
        max_study_simultaneous_allocations=2, new_independent_chains=0,
        physical_inputs_or_running_source_changed=False, holdout_generated_or_read=False,
        startup_claimed=False, new_DFT_calls_in_admission=0,
        journal_sha256=sha(journal_path), executed_helper_sha256=sha(Path(__file__)))
    with (ROOT/'admission_receipt.json').open('x') as stream:
        json.dump(result, stream, indent=2)
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reviewed-audit-sha256', required=True)
    print(json.dumps(release_once(parser.parse_args().reviewed_audit_sha256)))
