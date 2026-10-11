"""One reviewed continuation per registered parent; never replay an intent."""
import argparse
import contextlib
from datetime import datetime, timezone
import io
import json
import os
from pathlib import Path
import subprocess
import sys

from audit_native import REGISTRY, ROOT, RUNS, SUPPORT, SUPPORT_SHA, sha

SCRIPT = ROOT/'hf_hfo2_clamped_continue_E090_20261011.slurm'
SCRIPT_SHA = 'bf8529bd7f4b182c90f5038aed13059eab157d7cdfc94a2e275abd58d4a45bb5'
BINARY = Path('/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus')
BINARY_SHA = '88ea7f91c9f3091ef5fb4f410ab88267b14ecfe7f35c86dfaee9bf9319cace7a'
WAIT_ROOT = RUNS/'20261010-varneb-clamped-G2-E060-r1'
WAITERS = {'28709790':'PO_flip_T_pattern_preserving','28709791':'PO_flip_T_pattern_reversing'}


def run(argv):
    return subprocess.run(argv,check=True,capture_output=True,text=True,timeout=45).stdout.strip()


def fields(body):
    return dict(t.split('=',1) for t in body.split() if '=' in t)


def held_waiters():
    result = {}
    for job, channel in WAITERS.items():
        row = fields(run(['scontrol','show','job','-o',job]))
        expected = dict(JobId=job,JobName='hfo2-G2-E060-'+channel,JobState='PENDING',
            Partition='hfacnormal01',NumCPUs='32',Priority='0',Reason='JobHeldUser',Dependency='(null)',
            Command=str(REGISTRY['28709789']['script']),StdOut=str(WAIT_ROOT/(channel+'.slurm.out')),
            StdErr=str(WAIT_ROOT/(channel+'.slurm.err')))
        if any(row.get(k)!=v for k,v in expected.items()) or not row.get('UserId','').startswith('iai806('):
            raise ValueError('original held flip identity or current dependency differs')
        # Compare identity/dependency invariants, not changing scheduler timestamps.
        result[job] = {k:row[k] for k in (*expected,'UserId')}
    return result


def admitted_other(parent):
    allowed = set(REGISTRY)|set(WAITERS)
    for other in set(REGISTRY)-{parent}:
        path = ROOT/other/'dispatch_receipt.json'
        if not path.exists():
            continue
        r = json.loads(path.read_text())
        journal = ROOT/other/'dispatch_journal.jsonl'
        records = [json.loads(line) for line in journal.read_text().splitlines()]
        if (r.get('status')!='one_same_chain_continuation_submitted_and_released_not_start_proof'
                or r.get('parent_job')!=other or not str(r.get('job_id','')).isdecimal()
                or r.get('total_submissions')!=1 or r.get('release_writes')!=1
                or r.get('new_independent_chains')!=0 or sha(journal)!=r.get('journal_sha256')
                or [e['event'] for e in records]!=['one_held_submission_intent','one_accepted_handle',
                    'one_release_intent','one_existing_handle_released']
                or any(e.get('job_id')!=r['job_id'] for e in records[1:])):
            raise ValueError('other admission is ambiguous; reconcile it before any write')
        allowed.add(r['job_id'])
    return allowed


def spare_slot(body, parent, allowed):
    active = set()
    for line in body.splitlines():
        parts = [p.strip() for p in line.split('|')]
        if len(parts)!=5:
            raise ValueError('malformed queue')
        job,name,state,partition,cpus = parts
        if not name.startswith('hfo2'):
            continue
        if job not in allowed or partition!='hfacnormal01' or cpus!='32':
            raise ValueError('unexpected study allocation')
        if state in ('RUNNING','COMPLETING','CONFIGURING'):
            active.add(job)
        elif state!='PENDING':
            raise ValueError('unexpected study state')
    if parent in active or len(active)>1:
        raise ValueError('terminated parent and one free study slot required')
    return sorted(active)


def verify_audit(parent, reviewed_sha):
    if parent not in REGISTRY:
        raise ValueError('only registered G2 parents are admitted')
    r = REGISTRY[parent]
    terminal = ROOT/parent/r['audit_subdir']
    path = terminal/'audit_receipt.json'
    if sha(path)!=reviewed_sha or sha(SCRIPT)!=SCRIPT_SHA or sha(BINARY)!=BINARY_SHA:
        raise ValueError('reviewed receipt, wrapper or binary differs')
    receipt = json.loads(path.read_text())
    if (receipt['terminal_summary_sha256']!=sha(r['work']/'vcneb_summary.json')
            or (ROOT/parent/r['channel']/'band').exists()):
        raise ValueError('original terminal summary or fresh work namespace differs')
    for key,name in (('executed_auditor_sha256','audit_hfo2_clamped_terminal.py'),
                     ('executed_exporter_sha256','export_hfo2_clamped_observation.py'),
                     ('executed_residual_analyzer_sha256','analyze_hfo2_clamped_residual.py')):
        if receipt[key]!=SUPPORT_SHA[name] or sha(SUPPORT/name)!=SUPPORT_SHA[name]:
            raise ValueError('actual terminal support differs')
    # Reuse the exact already-tested startup validator, rather than a weaker
    # second implementation of the nine-cache/runtime/physical contract.
    body = SCRIPT.read_text().split(" 'import hashlib,json,pathlib,sys\n",1)[1]
    body = 'import hashlib,json,pathlib,sys\n'+body.split("' \"${seed_root}/manifest.json\"",1)[0]
    old = sys.argv
    try:
        sys.argv = ['reviewed-startup',str(terminal/'seed/manifest.json'),receipt['seed_manifest_sha256'],
                    str(path),reviewed_sha,parent,str(r['source'])]
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(body,str(SCRIPT),'exec'),{})
    finally:
        sys.argv = old
    return receipt


def dispatch_once(parent, reviewed_sha):
    if parent not in REGISTRY:
        raise ValueError('only registered G2 parents are admitted')
    root = ROOT/parent
    journal_path = root/'dispatch_journal.jsonl'
    if journal_path.exists() or (root/'dispatch_receipt.json').exists():
        raise FileExistsError('existing intent/handle; reconcile without resubmission or release retry')
    receipt = verify_audit(parent,reviewed_sha)
    accounting = run(['sacct','-X','-n','-P','-j',parent,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'])
    if accounting!=parent+'|COMPLETED|0:0|32|hfacnormal01|iai806':
        raise ValueError('audited parent must remain natively terminal')
    before = held_waiters()
    allowed = admitted_other(parent)
    queue = ['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']
    active = spare_slot(run(queue),parent,allowed)
    run(['bash','-n',str(SCRIPT)])
    r = REGISTRY[parent]
    name = 'hfo2-G2-E091-'+parent
    terminal = root/r['audit_subdir']
    env = dict(RUN_DFT='1',SOURCE_ROOT=str(r['source']),SEED_ROOT=str(terminal/'seed'),
        WORKDIR=str(root/r['channel']/'band'),PARENT_JOB_ID=parent,
        PARENT_AUDIT_PATH=str(terminal/'audit_receipt.json'),PARENT_AUDIT_SHA256=reviewed_sha,
        SEED_MANIFEST_SHA256=receipt['seed_manifest_sha256'])
    argv = ['sbatch','--parsable','--hold','--job-name='+name,
        '--export=ALL,'+','.join(k+'='+v for k,v in env.items()),
        '--output='+str(root/(r['channel']+'.slurm.out')),
        '--error='+str(root/(r['channel']+'.slurm.err')),str(SCRIPT)]
    with journal_path.open('x') as journal:
        def record(event,**data):
            journal.write(json.dumps(dict(event=event,check_UTC=datetime.now(timezone.utc).isoformat(),**data))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        record('one_held_submission_intent',argv=argv,parent_accounting=accounting,
            audit_receipt_sha256=reviewed_sha,active_jobs=active,unchanged_held_flips=before)
        try:
            handle = run(argv)
            if not handle.isdecimal():
                record('ambiguous_submission_reply',reply=handle)
                raise ValueError('ambiguous accepted handle; never resubmit')
            record('one_accepted_handle',job_id=handle)
            row = fields(run(['scontrol','show','job','-o',handle]))
            expected = dict(JobId=handle,JobName=name,JobState='PENDING',Reason='JobHeldUser',Priority='0',
                Partition='hfacnormal01',NumCPUs='32',Command=str(SCRIPT),Dependency='(null)',
                StdOut=str(root/(r['channel']+'.slurm.out')),StdErr=str(root/(r['channel']+'.slurm.err')))
            if any(row.get(k)!=v for k,v in expected.items()) or not row.get('UserId','').startswith('iai806('):
                raise ValueError('exact accepted held continuation not corroborated')
            spare_slot(run(queue),parent,allowed|{handle})
            if held_waiters()!=before:
                raise ValueError('held flips changed before release')
            record('one_release_intent',job_id=handle,accepted_native_state=row)
            reply = run(['scontrol','release',handle])
            record('one_existing_handle_released',job_id=handle,reply=reply)
            result = dict(status='one_same_chain_continuation_submitted_and_released_not_start_proof',
                job_id=handle,parent_job=parent,source_step=r['cap'],max_new_updates=40,wall_hours=12,
                total_submissions=1,release_writes=1,dependency_writes=0,new_independent_chains=0,
                audit_receipt_sha256=reviewed_sha,seed_manifest_sha256=receipt['seed_manifest_sha256'],
                script_sha256=sha(SCRIPT),executed_helper_sha256=sha(Path(__file__)),
                journal_sha256=sha(journal_path),physical_inputs_running_source_or_holdout_changed=False,
                FIRE_optimizer_state_restored=False,execution_started_claimed=False)
            with (root/'dispatch_receipt.json').open('x') as stream:
                json.dump(result,stream,indent=2,allow_nan=False)
        except Exception as exc:
            def evidence(value):
                return value.decode('utf-8','replace') if isinstance(value,bytes) else value
            record('reconcile_existing_intent_never_repeat_writes',exception=type(exc).__name__,message=str(exc),
                   stdout=evidence(getattr(exc,'stdout',None)),stderr=evidence(getattr(exc,'stderr',None)))
            raise
    print(json.dumps(result))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent',choices=tuple(REGISTRY),required=True)
    parser.add_argument('--audit-sha256',required=True)
    args = parser.parse_args()
    dispatch_once(args.parent,args.audit_sha256)
