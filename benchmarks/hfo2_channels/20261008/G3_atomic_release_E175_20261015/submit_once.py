"""One manual held/validated/released admission; ambiguous outcomes never retry."""
from datetime import datetime,timezone
from pathlib import Path
import argparse,hashlib,importlib.util,json,os,re,subprocess,sys
ROOT=Path('/public/home/iai806/abacus/agent-runs/20261015-varneb-G3-atomic-release-E175-r2')
BUNDLE=ROOT/'bundle'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def capture(args):return subprocess.run(args,check=True,capture_output=True,text=True).stdout.strip()
def save(path,data):
    with path.open('x') as f:
        json.dump(data,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def now():return datetime.now(timezone.utc).isoformat()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('bundle-sha256','config-sha256','preflight-sha256'):p.add_argument('--'+n,required=True)
    a=p.parse_args()
    if (ROOT/'submission_intent.json').exists():
        if (ROOT/'accepted.json').exists():
            print(json.dumps(dict(status='existing_handle_no_resubmission',accepted=json.loads((ROOT/'accepted.json').read_bytes()))));return
        raise RuntimeError('ambiguous existing intent; reconcile accounting, never retry')
    assert sha(BUNDLE/'bundle_files.json')==a.bundle_sha256
    files=json.loads((BUNDLE/'bundle_files.json').read_bytes())
    assert len(files)==5 and all(sha(BUNDLE/n)==h for n,h in files.items())
    assert sha(BUNDLE/'config.json')==a.config_sha256
    assert sha(ROOT/'native_preflight.json')==a.preflight_sha256
    proof=json.loads((ROOT/'native_preflight.json').read_bytes())
    assert proof['native_frozen_SCFs_reparsed']==8 and proof['new_DFT']==0
    assert not (ROOT/'atomic_release').exists()
    sys.path.insert(0,'/public/home/iai806/abacus/agent-runs/20261015-varneb-G3-faces-E161-r1/source')
    spec=importlib.util.spec_from_file_location('frozen_atomic_release_worker',BUNDLE/'worker.py')
    w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
    w.validate(json.loads((BUNDLE/'config.json').read_bytes()),BUNDLE)
    admission=w.scheduler_snapshot()
    if not admission['ready']:
        print(json.dumps(dict(status='resource_gate_not_ready_no_submission',admission=admission)));return
    assert all(r['name']=='hfo2-G3-E161-preserving' for r in admission['committed'])
    variables=dict(RUN_DFT='1',BUNDLE=str(BUNDLE),BUNDLE_MANIFEST_SHA256=a.bundle_sha256,CONFIG_SHA256=a.config_sha256)
    command=['sbatch','--hold','--parsable','--chdir='+str(ROOT),
        '--output='+str(ROOT/'atomic.slurm.out'),'--error='+str(ROOT/'atomic.slurm.err'),
        '--export=ALL,'+','.join(k+'='+v for k,v in variables.items()),str(BUNDLE/'run.slurm')]
    save(ROOT/'submission_intent.json',dict(time_UTC=now(),command=command,admission=admission,
        bundle_manifest_sha256=a.bundle_sha256,config_sha256=a.config_sha256,
        complete_native_preflight_sha256=a.preflight_sha256,submission_driver_sha256=sha(Path(__file__))))
    response=capture(command)
    assert re.fullmatch('[0-9]+(?:;[A-Za-z0-9_.-]+)?',response),response
    job=response.split(';')[0]
    save(ROOT/'accepted.json',dict(time_UTC=now(),job_id=job,sbatch_response=response))
    details=capture(['scontrol','show','job','-o',job]);fields=dict(t.split('=',1) for t in details.split() if '=' in t)
    assert fields['JobState']=='PENDING' and fields['Reason']=='JobHeldUser',details
    assert fields['Partition']=='hfacnormal01' and fields['NumNodes'] in ('1','1-1'),details
    assert fields['NumCPUs']=='32' and fields['NumTasks']=='1' and fields['CPUs/Task']=='32',details
    assert fields['Command']==str(BUNDLE/'run.slurm') and fields['TimeLimit'] in ('05:00:00','5:00:00'),details
    own=w.scheduler_snapshot(job)
    assert own['ready'] and len(own['committed'])<=2
    save(ROOT/'held_validation.json',dict(time_UTC=now(),job_id=job,scontrol=details,admission=own))
    released=capture(['scontrol','release',job])
    save(ROOT/'released.json',dict(time_UTC=now(),job_id=job,scontrol_response=released))
    queue=capture(['squeue','-h','-j',job,'-o','%i|%100j|%T|%P|%C|%D|%R'])
    save(ROOT/'first_queue.json',dict(time_UTC=now(),job_id=job,queue=queue))
    print(json.dumps(dict(status='submitted_once_validated_released',job_id=job,queue=queue,maximum_allocations=2)))


if __name__=='__main__':main()
