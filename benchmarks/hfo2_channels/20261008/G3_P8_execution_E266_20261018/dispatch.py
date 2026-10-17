"""One-shot held P8 submission and exact-handle release; no automatic retry."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import re
import subprocess

ROOT=Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-G3-P8-E266-r1')
command=lambda args:subprocess.check_output(args,text=True).strip()


def validate_held(raw,shard,job):
    fields=dict(re.findall(r'(\S+)=(\S+)',raw))
    required=dict(JobId=job,JobState='PENDING',Reason='JobHeldUser',Partition='hfacnormal01',
        NumCPUs='32',NumTasks='1',TimeLimit='04:00:00',JobName=f'hfo2-G3-E266-P{shard}')
    required['CPUs/Task']='32'
    if fields.get('NumNodes') not in ('1','1-1') or any(fields.get(k)!=v for k,v in required.items()):
        raise ValueError('held resources differ; keep held and inspect')
    return fields


def main(args):
    bundle=ROOT/'bundle';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    if sha(bundle/'config.json')!=args.config_sha256 or sha(bundle/'bundle_files.json')!=args.bundle_sha256:
        raise ValueError('exact frozen bundle/config required')
    files=json.loads((bundle/'bundle_files.json').read_bytes())
    if set(files)!={'worker.py','run.slurm','config.json','forecasts.json','starts.json','protocol.md'} or any(
            Path(n).name!=n or sha(bundle/n)!=h for n,h in files.items()):raise ValueError('immutable bundle drift')
    proof=json.loads((ROOT/'native_preflight.json').read_bytes())
    if proof['new_DFT']!=0 or len(proof['decks'])!=8 or proof['config_sha256']!=args.config_sha256:
        raise ValueError('complete eight-deck native preflight required')
    spec=importlib.util.spec_from_file_location('E266_dispatch_worker',bundle/'worker.py')
    w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w);old,_=w.old_worker()
    if args.mode=='submit-held':
        if w.queue()['committed']:raise ValueError('empty study queue required before two allocations')
        if any((ROOT/f'{name}_{s}.json').exists() for s in (0,1)
               for name in ('submission_intent','held_submission','accepted','released')):
            raise ValueError('prior intent/handle exists; inspect it, never resubmit')
        for shard in (0,1):
            if not w.queue()['ready']:raise ValueError('two-allocation cap')
            export=f'ALL,RUN_DFT=1,SHARD={shard},BUNDLE={bundle},BUNDLE_MANIFEST_SHA256={args.bundle_sha256},CONFIG_SHA256={args.config_sha256}'
            cmd=['sbatch','--hold','--parsable',f'--job-name=hfo2-G3-E266-P{shard}',f'--export={export}',
                 f'--output={ROOT}/slurm-%j.out',str(bundle/'run.slurm')]
            old.study.exclusive_json(ROOT/f'submission_intent_{shard}.json',dict(command=cmd,
                config_sha256=args.config_sha256,bundle_manifest_sha256=args.bundle_sha256))
            handle=command(cmd);job=handle.split(';')[0]
            if not job.isdigit():raise ValueError('unknown submission outcome; inspect queue, never resubmit')
            row=dict(shard=shard,job_id=job,config_sha256=args.config_sha256,
                bundle_manifest_sha256=args.bundle_sha256,submission_reply=handle,
                raw_scontrol=command(['scontrol','show','job','-o',job]))
            old.study.exclusive_json(ROOT/f'held_submission_{shard}.json',row);print(json.dumps(row),flush=True)
    else:
        rows=[]
        for shard in (0,1):
            row=json.loads((ROOT/f'held_submission_{shard}.json').read_bytes())
            if (ROOT/f'accepted_{shard}.json').exists() or (ROOT/f'released_{shard}.json').exists():
                raise ValueError('already accepted/released; inspect, never repeat')
            raw=command(['scontrol','show','job','-o',row['job_id']]);validate_held(raw,shard,row['job_id'])
            if row['config_sha256']!=args.config_sha256 or row['bundle_manifest_sha256']!=args.bundle_sha256:
                raise ValueError('held configuration differs; keep held')
            if not w.queue(row['job_id'],shard)['ready']:raise ValueError('committed allocation cap')
            rows.append((shard,row,raw))
        for shard,row,raw in rows:
            old.study.exclusive_json(ROOT/f'accepted_{shard}.json',dict(**row,validated_scontrol=raw,exact_one_node_range_accepted=True))
        for shard,row,_ in rows:
            old.study.exclusive_json(ROOT/f'release_intent_{shard}.json',dict(job_id=row['job_id']))
            reply=command(['scontrol','release',row['job_id']])
            result=dict(job_id=row['job_id'],reply=reply,raw_scontrol=command(['scontrol','show','job','-o',row['job_id']]))
            old.study.exclusive_json(ROOT/f'released_{shard}.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=('submit-held','release'))
    parser.add_argument('--config-sha256',required=True);parser.add_argument('--bundle-sha256',required=True)
    main(parser.parse_args())
