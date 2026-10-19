"""One compact read-only snapshot of the two original handles; no polling loop."""
import json
import math
from pathlib import Path
import subprocess

ROOT=Path('/public/home/iai806/abacus/agent-runs/20261019-varneb-preserving-local-E276-r2')
JOBS={'minus':'28954771','plus':'28954772'}


def snapshot():
    cfg=json.loads((ROOT/'bundle/prepared/seed/seed.json').read_bytes())
    q=cfg['orthonormal_39_by_4'];e0=json.loads((ROOT/'preflight.json').read_bytes())['frozen_energy_eV_cell']
    rows=[]
    for sign,job in JOBS.items():
        folder=ROOT/'results'/sign;records=[];partial=False
        for path in sorted(folder.glob('evaluations/eval-*/result.json')):
            try:record=json.loads(path.read_bytes())
            except json.JSONDecodeError:partial=True;continue
            records.append(record)
        a=folder/'treatment_A.json';done=folder/'complete.json'
        row=dict(sign=sign,job_id=job,completed_native_calls=len(records),partially_written_result=partial,
            stage='complete' if done.is_file() else 'A' if not a.is_file() else 'checks_or_J')
        if records:
            last=records[-1];g=last['gradient'];coeff=[sum(q[i][k]*g[i] for i in range(39)) for k in range(4)]
            translations=[sum(g[j] for j in range(k,36,3))/12 for k in range(3)]
            rg=[g[i]-sum(q[i][k]*coeff[k] for k in range(4))-translations[i%3] for i in range(36)]
            row.update(latest_trial_norm_A=math.sqrt(sum(v*v for v in last['coordinates'])),
                latest_atomic_released_gradient_eV_A=math.sqrt(sum(v*v for v in rg)),
                latest_energy_minus_F_meV_fu=250*(last['enthalpy_eV']-e0),
                latest_call_seconds=last['validation']['elapsed_seconds'])
        if a.is_file():
            treatment=json.loads(a.read_bytes());row['A_status']=treatment['status']
            if treatment['candidate'] is not None:row['A_gradient_eV_A']=treatment['candidate']['released_gradient_norm_eV_A']
        if done.is_file():
            final=json.loads(done.read_bytes());row['new_calls']=final['new_calls']
            row['treatments']=[dict(treatment=r['treatment'],status=r['status'],
                dE_F_meV_fu=None if r['candidate'] is None else 250*(r['candidate']['energy_eV_cell']-e0),
                gradient_eV_A=None if r['candidate'] is None else r['candidate']['released_gradient_norm_eV_A'])
                for r in final['treatments']]
        rows.append(row)
    queue=subprocess.check_output(['squeue','-j',','.join(JOBS.values()),'-h','-o','%i|%T|%C|%N|%M'],text=True)
    result=dict(queue=queue.strip().splitlines(),rows=rows,new_DFT_calls=0,scheduler_mutations=0,
        latest_trial_is_not_automatically_a_candidate=True)
    if any(r['stage']=='complete' for r in rows) or not queue.strip():
        result['accounting']=subprocess.check_output(['sacct','-n','-X','-P','-j',','.join(JOBS.values()),
            '--format=JobIDRaw,State,ExitCode,ElapsedRaw,Start,End'],text=True).strip().splitlines()
    return result


if __name__=='__main__':print(json.dumps(snapshot()))
