"""Hold exactly two registered pending successors for terminal material audit.

Single use: ambiguous writes are journalled and must be reconciled by reads,
never by rerunning this helper. Running parents and downstream +1% jobs are
not mutated. This is an execution guard, not a physical convergence test.
"""
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261011-varneb-terminal-guard-E066-r1')
PARENTS = {'28661019', '28722320'}
WAITERS = {
    '28692775': 'PO_to_T',
    '28692776': 'PO_flip_T_pattern_reversing',
}
ALLOWED = PARENTS | WAITERS.keys() | {'28709788', '28709789', '28709790', '28709791'}
PILOT = '/public/home/iai806/abacus/agent-runs/20261010-varneb-nested-controls-E055-r1/source/cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm'
OUTPUT_ROOT = '/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E058-r1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(argv):
    return subprocess.run(argv, check=True, capture_output=True, text=True, timeout=45).stdout.strip()


def fields(body):
    return dict(token.split('=', 1) for token in body.split() if '=' in token)


def dependencies(value):
    parts = value.split(',')
    if any(not re.fullmatch(r'afterok:\d+(?:\((?:unfulfilled|fulfilled)\))?(?::\d+(?:\((?:unfulfilled|fulfilled)\))?)*', p)
           for p in parts):
        raise ValueError('registered afterok AND dependencies required')
    return set(re.findall(r'(?:afterok:|:)(\d+)', value))


def waiter(body, job, *, held=False):
    if job not in WAITERS:
        raise ValueError('not a registered successor')
    row = fields(body)
    channel = WAITERS[job]
    if (row.get('JobId') != job or row.get('JobName') != 'hfo2-G2-E058-'+channel
            or not row.get('UserId', '').startswith('iai806(')
            or row.get('JobState') != 'PENDING' or row.get('Partition') != 'hfacnormal01'
            or row.get('NumCPUs') != '32' or row.get('Command') != PILOT
            or row.get('StdOut') != OUTPUT_ROOT+'/'+channel+'.slurm.out'
            or row.get('StdErr') != OUTPUT_ROOT+'/'+channel+'.slurm.err'
            or dependencies(row.get('Dependency', '')) != PARENTS):
        raise ValueError('successor identity/state/dependency changed; do not mutate')
    if held:
        if row.get('Priority') != '0' or row.get('Reason') != 'JobHeldUser':
            raise ValueError('held successor not yet corroborated')
    elif (row.get('Reason') != 'Dependency'
            or not row.get('Priority', '').isdecimal() or int(row['Priority']) <= 0):
        raise ValueError('not an unheld dependency waiter; ownership unproved')
    return row


def queue(body):
    active = set()
    for line in body.splitlines():
        parts = [p.strip() for p in line.split('|')]
        if len(parts) != 5:
            raise ValueError('malformed scheduler queue')
        job, name, state, partition, cpus = parts
        if not name.startswith('hfo2'):
            continue
        if job not in ALLOWED or partition != 'hfacnormal01' or cpus != '32':
            raise ValueError('unexpected study allocation')
        if state in ('RUNNING', 'COMPLETING', 'CONFIGURING'):
            active.add(job)
        elif state != 'PENDING':
            raise ValueError('unexpected study state')
    if active != PARENTS:
        raise ValueError('both original parents must still be live before guard placement')
    return sorted(active)


def acknowledge(job, *, lookup=run, pause=time.sleep):
    for delay in (0, 1, 2, 4, 8):
        if delay:
            pause(delay)
        body = lookup(['scontrol', 'show', 'job', '-o', job])
        row = fields(body)
        if row.get('Priority') == '0' and row.get('Reason') == 'JobHeldUser':
            return waiter(body, job, held=True)
        # Only the identical pre-write state can be a visibility lag.
        waiter(body, job)
    raise TimeoutError('hold acknowledgement stale; reconcile existing journal, never repeat hold')


def hold_once(root=ROOT, *, lookup=run, pause=time.sleep):
    active = queue(lookup(['squeue', '-h', '-u', 'iai806', '-o', '%i|%100j|%T|%P|%C']))
    before = {job: waiter(lookup(['scontrol', 'show', 'job', '-o', job]), job) for job in WAITERS}
    held_jobs = []
    with (root/'hold_journal.jsonl').open('x') as journal:
        def record(event, **data):
            journal.write(json.dumps(dict(event=event, check_CST=datetime.now().isoformat(), **data))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        record('guard_intent', active_parents=active, original_waiters=before,
               helper_sha256=sha(Path(__file__)))
        try:
            for job in WAITERS:
                waiter(lookup(['scontrol', 'show', 'job', '-o', job]), job)
                record('hold_intent', job_id=job)
                lookup(['scontrol', 'hold', job])
                row = acknowledge(job, lookup=lookup, pause=pause)
                held_jobs.append(job)
                record('hold_verified', job_id=job, actual=row)
        except Exception as exc:
            record('reconcile_do_not_rerun', exception=type(exc).__name__,
                   message=str(exc), verified_held_jobs=held_jobs)
            raise
    receipt = dict(status='two_registered_successors_held_for_terminal_material_audit',
        check_CST=datetime.now().isoformat(), held_jobs=held_jobs, parent_jobs=sorted(PARENTS),
        dependencies_changed=False, running_jobs_changed=False, new_submissions=0,
        new_DFT_calls=0, physical_inputs_changed=False, holdout_generated_or_read=False,
        helper_sha256=sha(Path(__file__)), journal_sha256=sha(root/'hold_journal.jsonl'),
        release_policy='Only after terminal audit and verified continuation handoff, or proven no continuation needed.')
    with (root/'hold_receipt.json').open('x') as stream:
        json.dump(receipt, stream, indent=2)
        stream.write('\n')
    return receipt


if __name__ == '__main__':
    print(json.dumps(hold_once()))
