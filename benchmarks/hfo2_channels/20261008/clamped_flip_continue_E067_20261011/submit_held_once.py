"""Submit one proved same-chain continuation held; no dependency/release writes."""
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261011-varneb-clamped-flip-E067-r1')
SOURCE = ROOT.parent/'20261010-varneb-clamped-G2-E054-r1/source-fixed'
WAITERS = {'28692775': 'PO_to_T', '28692776': 'PO_flip_T_pattern_reversing'}
OTHER = '28722320'
PARENT = '28661019'
ALLOWED = {OTHER, PARENT, *WAITERS, '28709788', '28709789', '28709790', '28709791'}
SCRIPT = ROOT/'hf_hfo2_clamped_flip_continue_E067_20261011.slurm'
SCRIPT_SHA = '91ee0f57dadbeb811d620f1287f346fae065ef496cead4b63e388d5110b92a09'
BINARY = Path('/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus')
BINARY_SHA = '88ea7f91c9f3091ef5fb4f410ab88267b14ecfe7f35c86dfaee9bf9319cace7a'
PILOT = '/public/home/iai806/abacus/agent-runs/20261010-varneb-nested-controls-E055-r1/source/cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm'
WAIT_ROOT = '/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E058-r1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(argv):
    return subprocess.run(argv, check=True, capture_output=True, text=True, timeout=45).stdout.strip()


def fields(body):
    return dict(token.split('=', 1) for token in body.split() if '=' in token)


def dependencies(value):
    if value == '(null)':
        return set()
    parts = value.split(',')
    if any(not re.fullmatch(r'afterok:\d+(?:\((?:unfulfilled|fulfilled)\))?(?::\d+(?:\((?:unfulfilled|fulfilled)\))?)*', p)
           for p in parts):
        raise ValueError('only registered afterok AND semantics allowed')
    return set(re.findall(r'(?:afterok:|:)(\d+)', value))


def held_waiter(body, job, *, required, allowed):
    if job not in WAITERS:
        raise ValueError('unregistered waiter')
    row = fields(body)
    channel = WAITERS[job]
    dep = dependencies(row.get('Dependency', ''))
    if (row.get('JobId') != job or row.get('JobName') != 'hfo2-G2-E058-'+channel
            or not row.get('UserId', '').startswith('iai806(')
            or row.get('JobState') != 'PENDING' or row.get('Partition') != 'hfacnormal01'
            or row.get('NumCPUs') != '32' or row.get('Priority') != '0'
            or row.get('Reason') != 'JobHeldUser' or row.get('Command') != PILOT
            or row.get('StdOut') != WAIT_ROOT+'/'+channel+'.slurm.out'
            or row.get('StdErr') != WAIT_ROOT+'/'+channel+'.slurm.err'
            or not set(required).issubset(dep) or not dep.issubset(set(allowed))):
        raise ValueError('held successor identity/dependencies changed')
    return row


def active_queue(body, *, allowed=ALLOWED):
    active = set()
    for line in body.splitlines():
        parts = [p.strip() for p in line.split('|')]
        if len(parts) != 5:
            raise ValueError('malformed queue')
        job, name, state, partition, cpus = parts
        if not name.startswith('hfo2'):
            continue
        if job not in allowed or partition != 'hfacnormal01' or cpus != '32':
            raise ValueError('unexpected study job')
        if state in ('RUNNING', 'COMPLETING', 'CONFIGURING'):
            active.add(job)
        elif state != 'PENDING':
            raise ValueError('unexpected study state')
    if active != {OTHER}:
        raise ValueError('exactly the registered M parent must remain live; reconcile other states')
    return sorted(active)


def verify_audit_seed():
    receipt = json.loads((ROOT/'audit_receipt.json').read_text())
    seed = ROOT/'seed'
    manifest = json.loads((seed/'manifest.json').read_text())
    names = {'initial.vasp', 'final.vasp', 'seed.traj', 'substrate.vasp', 'factory_parameters.json'}
    if (receipt['status'] != 'audited_same_flip_step20_continuation_ready'
            or receipt['parent_job'] != PARENT or receipt['ordinary_converged']
            or receipt['fresh_parent_interior_SCFs'] != 140
            or receipt['current_frame_exact_caches_checked'] != 9
            or not receipt['all_parent_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked']
            or receipt['new_DFT_calls'] or receipt['physical_inputs_changed'] or receipt['holdout_generated_or_read']
            or receipt['seed_manifest_sha256'] != sha(seed/'manifest.json')
            or manifest['source_job_id'] != PARENT or manifest['source_step'] != 20
            or set(manifest['files_sha256']) != names
            or any(sha(seed/n) != h for n,h in manifest['files_sha256'].items())
            or sha(SCRIPT) != SCRIPT_SHA or sha(BINARY) != BINARY_SHA
            or receipt['executed_audit_source_sha256'] != sha(ROOT/'audit_HF.py')
            or any(sha(ROOT/'audit_source/scripts'/n) != h for n,h in receipt['analysis_source_sha256'].items())
            or any(sha(SOURCE/n) != h for n,h in receipt['runtime_code_sha256'].items())
            or (ROOT/'PO_flip_T_pattern_preserving/band').exists()):
        raise ValueError('actual audited source/seed/script changed or work namespace already exists')
    return receipt


def validate_held_continuation(body, handle):
    row = fields(body)
    if (row.get('JobId') != handle or row.get('JobName') != 'hfo2-G2-flip-E067'
            or not row.get('UserId', '').startswith('iai806(')
            or row.get('JobState') != 'PENDING' or row.get('Reason') != 'JobHeldUser'
            or row.get('Priority') != '0' or row.get('Partition') != 'hfacnormal01'
            or row.get('NumCPUs') != '32' or row.get('Command') != str(SCRIPT)
            or row.get('StdOut') != str(ROOT/'flip.slurm.out')
            or row.get('StdErr') != str(ROOT/'flip.slurm.err')
            or dependencies(row.get('Dependency', ''))):
        raise ValueError('exact accepted held continuation not corroborated')
    return row


def submit_once():
    if (ROOT/'submission_journal.jsonl').exists() or (ROOT/'accepted_receipt.json').exists():
        raise FileExistsError('existing dispatch evidence; reconcile the accepted handle, never submit again')
    receipt = verify_audit_seed()
    parent = run(['sacct', '-X', '-n', '-P', '-j', PARENT,
                  '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'])
    if parent != PARENT+'|COMPLETED|0:0|32|hfacnormal01|iai806':
        raise ValueError('terminal audited parent no longer corroborated')
    before = {job: held_waiter(run(['scontrol','show','job','-o',job]), job,
                              required={OTHER}, allowed={OTHER,PARENT}) for job in WAITERS}
    active = active_queue(run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']))
    run(['bash','-n',str(SCRIPT)])
    with (ROOT/'submission_journal.jsonl').open('x') as journal:
        def record(event, **data):
            journal.write(json.dumps(dict(event=event, check_CST=datetime.now().isoformat(), **data))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        argv = ['sbatch', '--parsable', '--hold', '--job-name=hfo2-G2-flip-E067',
            '--export=ALL,RUN_DFT=1,SEED_ROOT='+str(ROOT/'seed')+',WORKDIR='+str(ROOT/'PO_flip_T_pattern_preserving/band')
            +',SEED_MANIFEST_SHA256='+receipt['seed_manifest_sha256'],
            '--output='+str(ROOT/'flip.slurm.out'), '--error='+str(ROOT/'flip.slurm.err'), str(SCRIPT)]
        record('single_held_submission_intent', argv=argv, parent_accounting=parent,
               active_jobs=active, held_waiters=before, audit_receipt_sha256=sha(ROOT/'audit_receipt.json'))
        try:
            handle = run(argv)
            if not handle.isdecimal():
                raise ValueError('ambiguous accepted handle; reconcile by identity, never resubmit')
            record('accepted_held_handle', job_id=handle)
            dispatch = dict(status='one_accepted_held_same_chain_continuation_not_released',
                check_CST=datetime.now().isoformat(), job_id=handle, parent_job_id=PARENT,
                total_submissions=1, new_independent_chains=0, new_DFT_calls_before_release=0,
                seed_manifest_sha256=receipt['seed_manifest_sha256'], audit_receipt_sha256=sha(ROOT/'audit_receipt.json'),
                executed_helper_sha256=sha(Path(__file__)), script_sha256=sha(SCRIPT),
                running_source_or_physical_inputs_changed=False, holdout_generated_or_read=False,
                dependencies_changed=False, execution_started_claimed=False)
            with (ROOT/'accepted_receipt.json').open('x') as stream:
                json.dump(dispatch, stream, indent=2)
                stream.write('\n')
        except Exception as exc:
            record('reconcile_existing_intent_do_not_resubmit', exception=type(exc).__name__, message=str(exc))
            raise
    print(json.dumps(dispatch))


if __name__ == '__main__':
    submit_once()
