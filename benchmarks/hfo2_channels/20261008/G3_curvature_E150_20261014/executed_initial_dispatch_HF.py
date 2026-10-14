"""One two-handle admission, never replay an uncertain submission/release."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT=Path('/public/home/iai806/abacus/agent-runs/20261014-varneb-G3-local-curvature-E150-r1')
CHANNELS=('T','preserving')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(argv):
    return subprocess.run(argv,check=True,capture_output=True,text=True,timeout=45).stdout.strip()


def queue(accepted):
    body=run(['squeue','-h','-u','iai806','-p','hfacnormal01','-o','%i|%100j|%T|%C|%R'])
    rows=[line.split('|') for line in body.splitlines() if line.strip()]
    if (len(rows)>2 or {r[0] for r in rows} != set(accepted)
            or any(len(r)!=5 or r[1].strip()!=accepted[r[0]] for r in rows)):
        raise ValueError('unknown/missing allocation or two-slot conflict; reconcile, never resubmit')
    return body


def verify(preflight_sha):
    receipt=ROOT/'preflight_receipt.json'
    if sha(receipt)!=preflight_sha:
        raise ValueError('reviewed preflight receipt differs')
    record=json.loads(receipt.read_text()); source=ROOT/'source'
    if (record.get('status')!='native_two_center_preflight_passed_not_submitted'
            or record.get('new_DFT_calls')!=0 or record.get('scheduler_mutations')!=0
            or record.get('source_directory')!=str(source)
            or record.get('critical_source_text_hashes_checked')!=15
            or record.get('manifest_sha256')!='5f047f54ca83cb2dceb8cb71687eb9206a85c04205cdf26431565bc8a9b98a5a'):
        raise ValueError('actual reviewed finite native preflight required')
    manifest=source/'benchmarks/hfo2_channels/20261008/G3_curvature_E150_20261014/manifest.json'
    wrapper=source/'cluster/hf_hfo2_G3_local_curvature_E150_20261014.slurm'
    if (sha(manifest)!=record['manifest_sha256'] or sha(wrapper)!=record['wrapper_sha256']
            or sha(source/'runtime_files.json')!=record['runtime_manifest_sha256']):
        raise ValueError('immutable deployed inputs differ')
    files=json.loads((source/'runtime_files.json').read_text())
    if len(files)!=record['runtime_python_files'] or any(sha(source/n)!=h for n,h in files.items()):
        raise ValueError('runtime byte drift')
    return record,source,manifest,wrapper


def dispatch(preflight_sha):
    journal=ROOT/'dispatch_journal.jsonl'; receipt=ROOT/'dispatch_receipt.json'
    if journal.exists() or receipt.exists():
        raise FileExistsError('existing intent/handle: inspect journal+sacct, do not replay')
    reviewed,source,manifest,wrapper=verify(preflight_sha)
    queue({})
    run(['bash','-n',str(wrapper)])
    accepted={}; handles={}
    with journal.open('x') as stream:
        def record(event,**fields):
            stream.write(json.dumps(dict(event=event,**fields))+'\n');stream.flush();os.fsync(stream.fileno())
        for channel in CHANNELS:
            queue(accepted)
            name='hfo2-G3-E150-'+channel
            env=dict(RUN_DFT='1',SOURCE_ROOT=str(source),MANIFEST=str(manifest),
                MANIFEST_SHA256=reviewed['manifest_sha256'],CHANNEL=channel,
                WORKDIR=str(ROOT/channel),RUNTIME_MANIFEST_SHA256=reviewed['runtime_manifest_sha256'])
            if (ROOT/channel).exists():
                raise FileExistsError('fresh output namespace required')
            argv=['sbatch','--parsable','--hold','--job-name='+name,'--chdir='+str(ROOT),
                '--export=ALL,'+','.join(k+'='+v for k,v in env.items()),
                '--output='+str(ROOT/(channel+'.slurm.out')),
                '--error='+str(ROOT/(channel+'.slurm.err')),str(wrapper)]
            record('one_submission_intent',channel=channel,argv=argv)
            job=run(argv).split(';')[0]
            if not job.isdecimal():
                raise ValueError('uncertain sbatch outcome; inspect native accounting, never resubmit')
            record('one_handle_accepted',channel=channel,job_id=job)
            accepted[job]=name;handles[channel]=job
        queue(accepted)
        native={}
        for channel,job in handles.items():
            body=run(['scontrol','show','job','-o',job]);native[channel]=body
            required=('JobId='+job+' ','JobName='+accepted[job]+' ','Partition=hfacnormal01 ',
                'NumNodes=1 ','NumCPUs=32 ','NumTasks=1 ','CPUs/Task=32 ',
                'Command='+str(wrapper)+' ')
            if any(item not in body+' ' for item in required):
                raise ValueError('native held allocation/script differs; retain handles and inspect')
        for channel,job in handles.items():
            record('one_release_intent',channel=channel,job_id=job)
            run(['scontrol','release',job])
            record('one_existing_handle_released',channel=channel,job_id=job)
        observed=queue(accepted)
        if any(line.rsplit('|',1)[-1].strip()=='JobHeldUser' for line in observed.splitlines()):
            raise ValueError('release not yet independently observed; keep handles and reconcile')
        record('both_release_states_independently_observed',job_ids=list(handles.values()))
    result=dict(status='two_registered_curvature_batches_submitted_released_not_completed',
        handles=handles,preflight_receipt_sha256=preflight_sha,journal_sha256=sha(journal),
        helper_sha256=sha(Path(__file__)),reviewed_runtime=reviewed,
        native_held_scontrol=native,post_release_queue=observed,
        total_submissions=2,total_releases=2,new_production_chains=0,
        total_planned_static_attempts=296,physical_holdout_generated_or_read=False)
    with receipt.open('x') as output: json.dump(result,output,indent=2)
    print(json.dumps(dict(handles=handles,status=result['status'])))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preflight-sha256',required=True)
    dispatch(parser.parse_args().preflight_sha256)
