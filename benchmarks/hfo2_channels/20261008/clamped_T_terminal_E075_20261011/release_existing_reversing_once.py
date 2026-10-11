"""One audited free-slot admission of existing reversing job; no submission."""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path

RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS/'20261011-varneb-clamped-T-terminal-E075-r1'
COMMON = RUNS/'20261011-varneb-clamped-M-terminal-E068-source-r2/release_existing_PO_to_T.py'
COMMON_SHA = '5ddca4a3bf73ed8ad1d8f348544d97bf170ef442354a63a1187ae49a5c56c0f0'
DONE, LIVE, TARGET = '28692775', '28723655', '28692776'
WAIT = RUNS/'20261010-varneb-clamped-G2-E058-r1'
CHANNEL = 'PO_flip_T_pattern_reversing'
DOWNSTREAM = {'28709788': 'PO_to_T', '28709789': 'PO_to_M',
    '28709790': 'PO_flip_T_pattern_preserving', '28709791': 'PO_flip_T_pattern_reversing'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verified_audit(receipt):
    force = receipt.get('terminal_observation', {}).get('replayed_fmax_eV_A')
    expected = dict(job_id=DONE, ordinary_converged=True,
        status='audited_ordinary_converged_not_TS_or_sampling_certificate',
        new_DFT_calls=0, scheduler_mutations=0, physical_inputs_or_running_source_changed=False,
        holdout_generated_or_read=False,
        all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked=True,
        initial_fresh_interior_SCFs=7)
    if (any(receipt.get(k) != v for k,v in expected.items())
            or receipt.get('terminal_observation', {}).get('ordinary_residual_pass') is not True
            or isinstance(force, bool) or not isinstance(force, (int,float))
            or not math.isfinite(force) or not 0 <= force <= .10
            or receipt.get('fresh_interior_SCFs') != 7*(receipt.get('terminal_step', -1)+1)):
        raise ValueError('actual converged fresh-starter audit required')


def one_used_slot(body):
    active = set()
    seen = set()
    for line in body.splitlines():
        parts = [v.strip() for v in line.split('|')]
        if len(parts) != 5:
            raise ValueError('malformed native queue')
        job, name, state, partition, cpus = parts
        if not name.startswith('hfo2'):
            continue
        if (job not in {LIVE, TARGET, *DOWNSTREAM} or job in seen
                or partition != 'hfacnormal01' or cpus != '32'):
            raise ValueError('unexpected study allocation')
        seen.add(job)
        if state in {'RUNNING','COMPLETING','CONFIGURING'}:
            active.add(job)
        elif state != 'PENDING':
            raise ValueError('unexpected study state')
    if active != {LIVE} or not {TARGET,*DOWNSTREAM}.issubset(seen):
        raise ValueError('one registered live continuation and five protected waiters required')
    return sorted(active)


def load_common():
    spec = importlib.util.spec_from_file_location('E068_shared_preflight', COMMON)
    common = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(common)
    return common


def release_once(audit_sha256):
    journal_path = ROOT/'admission_journal.jsonl'
    if journal_path.exists():
        raise FileExistsError('reconcile existing intent; never repeat release')
    if sha(COMMON) != COMMON_SHA or sha(ROOT/'audit_receipt.json') != audit_sha256:
        raise ValueError('reviewed audit or pinned shared preflight changed')
    common = load_common()
    receipt = json.loads((ROOT/'audit_receipt.json').read_text())
    verified_audit(receipt)
    if sha(Path(receipt['workdir'])/'vcneb_summary.json') != receipt['terminal_summary_sha256']:
        raise ValueError('terminal summary changed')
    rules = common.load_rules()
    if rules.run(['sacct','-X','-n','-P','-j',DONE,
                  '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User']) != DONE+'|COMPLETED|0:0|32|hfacnormal01|iai806':
        raise ValueError('native completed parent required')
    if rules.run(['sacct','-X','-n','-P','-j',LIVE,
                  '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User']) != LIVE+'|RUNNING|0:0|32|hfacnormal01|iai806':
        raise ValueError('native registered live continuation required')
    if sha(WAIT/'preflight.json') != common.PREFLIGHT_SHA:
        raise ValueError('registered pending preflight changed')
    proof = json.loads((WAIT/'preflight.json').read_text())
    source = Path(proof['source_root'])
    if (sha(Path(rules.PILOT)) != proof['pilot_sha256']
            or sha(rules.BINARY) != proof['ABACUS_binary_sha256']
            or any(sha(source/n) != h for n,h in proof['source_code_sha256'].items())):
        raise ValueError('pending immutable runtime or binary changed')
    seed_record = next(r for r in proof['reports'] if r['channel'] == CHANNEL)
    seed = Path(seed_record['seed_root'])
    if (sha(seed/'manifest.json') != seed_record['manifest_sha256']
            or any(sha(seed/n) != h for n,h in seed_record['seed_files_sha256'].items())
            or (WAIT/CHANNEL/'band').exists()):
        raise ValueError('pending seed changed or production already exists')
    pending = common.pending_seed_preflight(seed)
    row = rules.held_waiter(rules.run(['scontrol','show','job','-o',TARGET]), TARGET,
                            required=set(), allowed=set())
    for job, channel in DOWNSTREAM.items():
        values = rules.fields(rules.run(['scontrol','show','job','-o',job]))
        if (values.get('JobId') != job or values.get('JobName') != 'hfo2-G2-E060-'+channel
                or not values.get('UserId', '').startswith('iai806(')
                or values.get('JobState') != 'PENDING' or values.get('Priority') != '0'
                or values.get('Reason') != 'JobHeldUser' or values.get('NumCPUs') != '32'
                or values.get('Partition') != 'hfacnormal01' or values.get('Command') != rules.PILOT
                or values.get('StdOut') != str(RUNS/'20261010-varneb-clamped-G2-E060-r1'/(channel+'.slurm.out'))):
            raise ValueError('all four downstream user holds required')
    active = one_used_slot(rules.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']))
    with journal_path.open('x') as journal:
        def record(event, **values):
            journal.write(json.dumps(dict(event=event, checked_UTC=datetime.now(timezone.utc).isoformat(), **values))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        record('admission_verified', target=TARGET, parent_audit_sha256=audit_sha256,
               active_jobs=active, held_row=row, pending_seed_preflight=pending,
               executed_helper_sha256=sha(Path(__file__)))
        try:
            one_used_slot(rules.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']))
            rules.held_waiter(rules.run(['scontrol','show','job','-o',TARGET]), TARGET,
                              required=set(), allowed=set())
            record('release_intent', job_id=TARGET)
            rules.run(['scontrol','release',TARGET])
            record('release_accepted', job_id=TARGET)
        except Exception as exc:
            record('reconcile_do_not_repeat', error=type(exc).__name__, message=str(exc))
            raise
    report = dict(status='one_existing_reversing_candidate_released', target_job=TARGET,
        audited_completed_parent=DONE, live_continuation=LIVE, new_submissions=0,
        dependency_writes=0, release_writes=1, max_study_simultaneous_allocations=2,
        physical_inputs_or_running_source_changed=False, holdout_generated_or_read=False,
        new_independent_chains=0, startup_claimed=False, new_DFT_calls_in_admission=0,
        journal_sha256=sha(journal_path), executed_helper_sha256=sha(Path(__file__)))
    with (ROOT/'admission_receipt.json').open('x') as stream:
        json.dump(report, stream, indent=2)
    return report


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reviewed-audit-sha256', required=True)
    print(json.dumps(release_once(parser.parse_args().reviewed_audit_sha256)))
