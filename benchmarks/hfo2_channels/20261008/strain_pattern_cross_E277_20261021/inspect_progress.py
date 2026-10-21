"""Compact read-only progress for the two original handles, never resubmit."""
import json
from pathlib import Path
import subprocess
ROOT=Path('/public/home/iai806/abacus/agent-runs/20261021-varneb-strain-pattern-E277-r1')
rows=[]
for name in ('zero_q1','plus1_q0'):
    accepted=json.loads((ROOT/f'accepted_{name}.json').read_bytes());job=accepted['job_id']
    scheduler=subprocess.check_output(['sacct','-n','-X','-P','-j',job,
        '--format=JobIDRaw,State,ExitCode,Elapsed,Start,End'],text=True).strip()
    folder=ROOT/'results'/name;records=sorted(folder.glob('evaluations/eval-*/result.json'))
    row=dict(condition=name,job_id=job,scheduler=scheduler,completed_SCFs=len(records))
    if records:
        last=json.loads(records[-1].read_bytes());v=last['validation']
        row.update(last_energy_eV_cell=last['enthalpy_eV'],last_SCF_seconds=v['elapsed_seconds'],
            genuine_MPI_ranks=v['genuine_MPI_ranks'],physical_input_sha256=v['input_sha256'])
    if (folder/'complete.json').exists():
        terminal=json.loads((folder/'complete.json').read_bytes());candidate=terminal['candidate']
        row['terminal_result']={k:terminal[k] for k in ('status','new_calls','job_id','bundle_pin')}
        row['terminal_result']['candidate']=None if candidate is None else {
            k:candidate[k] for k in ('energy_eV_cell','released_gradient_eV_A','coordinate_norm_A')}
    rows.append(row)
print(json.dumps(rows,indent=2))
