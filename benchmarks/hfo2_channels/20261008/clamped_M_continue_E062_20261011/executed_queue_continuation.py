"""One held continuation, dependency handoff, then release; no polling or retry."""
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261011-varneb-clamped-M-E062-r1')
WAITING = ('28692775', '28692776')
ALLOWED = {'28661019', *WAITING, '28709788', '28709789', '28709790', '28709791'}
SCRIPT = ROOT/'hf_hfo2_clamped_M_continue_E062_20261011.slurm'
SCRIPT_SHA = '2fa864368bc4e55f06e43e1bc0786e80d07f6d8113b82d485ab00d30751201d9'
BINARY = Path('/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus')
BINARY_SHA = '88ea7f91c9f3091ef5fb4f410ab88267b14ecfe7f35c86dfaee9bf9319cace7a'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(argv):
    return subprocess.run(argv, check=True, capture_output=True, text=True, timeout=45).stdout.strip()


def fields(body):
    return dict(field.split('=', 1) for field in body.split() if '=' in field)


def dependency_handles(value):
    # This exact handoff only accepts afterok AND dependencies, never OR or
    # afterany. Slurm removes fulfilled parents from the visible expression.
    if value == '(null)':
        return set()
    if not re.fullmatch(r'afterok:\d+(?:\((?:unfulfilled|fulfilled)\))?(?::\d+(?:\((?:unfulfilled|fulfilled)\))?)*', value):
        raise ValueError('unexpected dependency type/format')
    return set(re.findall(r'(?:afterok:|:)(\d+)', value))


def validate_waiter(body, job, expected):
    row = fields(body)
    if (row.get('JobId') != job or row.get('JobState') != 'PENDING'
            or row.get('Partition') != 'hfacnormal01' or row.get('NumCPUs') != '32'
            or not row.get('JobName', '').startswith('hfo2-G2-E058-')
            or not row.get('UserId', '').startswith('iai806(')
            or dependency_handles(row.get('Dependency', '')) != set(expected)):
        raise ValueError('waiting registered wave changed; reconcile, do not mutate')
    return row


def validate_active_queue(body, allowed):
    active = []
    for line in body.splitlines():
        parts = [p.strip() for p in line.split('|')]
        if len(parts) != 5:
            raise ValueError('malformed queue')
        if parts[1].startswith('hfo2'):
            if parts[0] not in allowed or parts[3:] != ['hfacnormal01', '32']:
                raise ValueError('unexpected study job')
            if parts[2] in ('RUNNING', 'COMPLETING', 'CONFIGURING'):
                active.append(parts[0])
            elif parts[2] != 'PENDING':
                raise ValueError('unrecognized study state')
    if len(active) > 1:
        raise ValueError('no spare slot for the single continuation')
    return active


def submit_once():
    receipt = json.loads((ROOT/'audit_receipt.json').read_text())
    seed = ROOT/'seed'
    manifest = json.loads((seed/'manifest.json').read_text())
    if (receipt['status'] != 'audited_same_chain_step20_continuation_ready'
            or not receipt['all_parent_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked']
            or receipt['current_frame_exact_caches_checked'] != 9 or receipt['fresh_parent_interior_SCFs'] != 140
            or receipt['ordinary_converged'] or receipt['new_DFT_calls']
            or receipt['physical_inputs_changed'] or receipt['holdout_generated_or_read']
            or receipt['seed_manifest_sha256'] != sha(seed/'manifest.json')
            or manifest['source_job_id'] != '28661020' or manifest['source_step'] != 20
            or any(sha(seed/n) != h for n,h in manifest['files_sha256'].items())
            or sha(SCRIPT) != SCRIPT_SHA or (ROOT/'PO_to_M/band').exists()):
        raise ValueError('audited seed/script changed')
    runtime = ROOT.parent/'20261010-varneb-clamped-G2-E054-r1/source-fixed'
    if sha(BINARY) != BINARY_SHA or any(sha(runtime/n) != h for n,h in receipt['runtime_code_sha256'].items()):
        raise ValueError('original binary or immutable runtime changed after audit')
    parent = run(['sacct','-X','-n','-P','-j','28661020',
                  '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'])
    if parent != '28661020|COMPLETED|0:0|32|hfacnormal01|iai806':
        raise ValueError('terminal parent proof changed')
    before = {job:validate_waiter(run(['scontrol','show','job','-o',job]), job, {'28661019'}) for job in WAITING}
    queue = run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C'])
    active = validate_active_queue(queue, ALLOWED)
    run(['bash','-n',str(SCRIPT)])
    jobs = []
    # An ambiguous or partial dispatch leaves a durable journal; never rerun.
    with (ROOT/'submission_journal.jsonl').open('x') as journal:
        def record(event, **data):
            journal.write(json.dumps(dict(event=event, check_CST=datetime.now().isoformat(), **data))+'\n')
            journal.flush(); os.fsync(journal.fileno())
        argv = ['sbatch','--parsable','--hold','--job-name=hfo2-G2-M-E062',
            '--export=ALL,RUN_DFT=1,SEED_ROOT='+str(seed)+',WORKDIR='+str(ROOT/'PO_to_M/band'),
            '--output='+str(ROOT/'PO_to_M.slurm.out'),'--error='+str(ROOT/'PO_to_M.slurm.err'),str(SCRIPT)]
        record('intent_held_submission', argv=argv, predecessor=parent, waiters=before,
               active_study_jobs=active, audit_receipt_sha256=sha(ROOT/'audit_receipt.json'))
        try:
            handle = run(argv)
            if not handle.isdecimal():
                raise ValueError('ambiguous accepted handle')
            record('accepted_held_handle', job_id=handle)
            held = fields(run(['scontrol','show','job','-o',handle]))
            if held.get('JobState') != 'PENDING' or held.get('Reason') != 'JobHeldUser':
                raise ValueError('held continuation not proved')
            dependency = 'afterok:28661019:'+handle
            for job in WAITING:
                validate_waiter(run(['scontrol','show','job','-o',job]), job, {'28661019'})
                record('dependency_update_intent', job_id=job, dependency=dependency)
                run(['scontrol','update','JobId='+job,'Dependency='+dependency])
                validate_waiter(run(['scontrol','show','job','-o',job]), job, {'28661019',handle})
                record('dependency_update_verified', job_id=job, dependency=dependency)
                jobs.append(job)
            # Successors cannot start while the inserted prerequisite is held.
            validate_active_queue(run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']), ALLOWED|{handle})
            record('release_intent', job_id=handle, both_successors_verified=True)
            run(['scontrol','release',handle])
            record('released', job_id=handle)
        except Exception as exc:
            record('reconcile_do_not_retry', exception=type(exc).__name__, message=str(exc),
                   waiters_already_updated=jobs)
            raise
    result = dict(status='same_chain_continuation_released_after_dependency_handoff',
        check_CST=datetime.now().isoformat(), job_id=handle, parent_job_id='28661020',
        dependency=dependency, waiting_jobs_updated=list(WAITING),
        max_study_simultaneous_allocations=2, new_independent_chains=0,
        physical_inputs_changed=False, running_source_overwritten=False,
        held_out_condition_generated_or_read=False, old_parent_not_modified=True,
        new_segment_steps=20, walltime_hours=6, MPI_ranks=32, ordinary_fmax=.10,
        audit_receipt_sha256=sha(ROOT/'audit_receipt.json'), script_sha256=sha(SCRIPT),
        queue_helper_sha256=sha(Path(__file__)), journal_sha256=sha(ROOT/'submission_journal.jsonl'))
    with (ROOT/'submission_receipt.json').open('x') as stream:
        json.dump(result, stream, indent=2); stream.write('\n')
    print(json.dumps(result))


if __name__ == '__main__':
    submit_once()
