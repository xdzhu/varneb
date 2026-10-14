"""Reconcile the two accepted E150 handles; this file cannot submit a job."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path

ROOT=Path('/public/home/iai806/abacus/agent-runs/20261014-varneb-G3-local-curvature-E150-r1')
HANDLES={'T':'28855614','preserving':'28855615'}
ORIGINAL_SHA='1e34604ba6eefdc54aa57f758db5c4e4955ec6496be3fef66cfce17a403a7bfa'
PREFLIGHT_SHA='aebad0bc4931a892be4da9509cd87df6cbec2debe0ec0749e36371983391525e'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def original():
    path=ROOT/'dispatch_HF.py'
    if sha(path)!=ORIGINAL_SHA: raise ValueError('preserve the original executed dispatcher')
    spec=importlib.util.spec_from_file_location('executed_E150_dispatch',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def release():
    journal=ROOT/'release_existing_journal.jsonl'; target=ROOT/'reconciled_dispatch_receipt.json'
    if journal.exists() or target.exists():
        raise FileExistsError('existing release intent: inspect same handles, never replay')
    old=original(); reviewed,source,manifest,wrapper=old.verify(PREFLIGHT_SHA)
    path=ROOT/'dispatch_journal.jsonl'
    entries=[json.loads(line) for line in path.read_text().splitlines()]
    if ([e['event'] for e in entries]!=['one_submission_intent','one_handle_accepted']*2
            or {e['channel']:e['job_id'] for e in entries if e['event']=='one_handle_accepted'}!=HANDLES):
        raise ValueError('exactly two already accepted, unreleased handles required')
    accepted={j:'hfo2-G3-E150-'+c for c,j in HANDLES.items()}
    old.queue(accepted)
    native={}
    for channel,job in HANDLES.items():
        body=old.run(['scontrol','show','job','-o',job]);native[channel]=body
        tokens=dict(s.split('=',1) for s in body.split() if '=' in s)
        required=dict(JobId=job,JobName=accepted[job],Partition='hfacnormal01',JobState='PENDING',
            Reason='JobHeldUser',NumCPUs='32',NumTasks='1',Command=str(wrapper),WorkDir=str(ROOT),
            StdOut=str(ROOT/(channel+'.slurm.out')),StdErr=str(ROOT/(channel+'.slurm.err')))
        if (any(tokens.get(k)!=v for k,v in required.items()) or tokens.get('NumNodes') not in ('1','1-1')
                or tokens.get('CPUs/Task')!='32' or (ROOT/channel).exists()):
            raise ValueError('actual held source/resources/workspace differ')
    with journal.open('x') as stream:
        def record(event,**data):
            stream.write(json.dumps(dict(event=event,**data))+'\n');stream.flush();os.fsync(stream.fileno())
        record('reviewed_original_two_handle_state',initial_journal_sha256=sha(path),native_scontrol=native)
        for channel,job in HANDLES.items():
            record('one_existing_handle_release_intent',channel=channel,job_id=job)
            old.run(['scontrol','release',job])
            record('one_existing_handle_released',channel=channel,job_id=job)
        observed=old.queue(accepted)
        if any(line.rsplit('|',1)[-1].strip().strip('()')=='JobHeldUser' for line in observed.splitlines()):
            raise ValueError('release not observed: keep handles and reconcile without reissuing')
        record('release_states_independently_observed',queue=observed)
    result=dict(status='original_two_handles_reconciled_released_not_completed',handles=HANDLES,
        initial_dispatcher_sha256=ORIGINAL_SHA,executed_release_helper_sha256=sha(Path(__file__)),
        original_dispatch_journal_sha256=sha(path),release_journal_sha256=sha(journal),
        reviewed_preflight_sha256=PREFLIGHT_SHA,reviewed_runtime=reviewed,native_held_scontrol=native,
        post_release_queue=observed,total_submissions=2,additional_submissions=0,total_release_writes=2,
        root_cause='Slurm exact one-node pending range NumNodes=1-1 was refused by overly literal metadata gate',
        physical_inputs_wrapper_runtime_or_optimizer_changed=False,new_production_chains=0,
        physical_holdout_generated_or_read=False)
    with target.open('x') as stream:json.dump(result,stream,indent=2)
    print(json.dumps(dict(status=result['status'],handles=HANDLES)))
    return result


if __name__=='__main__': release()
