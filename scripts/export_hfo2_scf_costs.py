"""Read-only cost export of the two completed E053 G2 pilot segments.

This deliberately does NOT read the running E054 continuations, held-out
strain, endpoint costs or Hessian probes. It cannot submit or run a DFT job.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import subprocess

import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.export_hfo2_clamped_observation import PILOT_SCRIPT_SHA256
from vcneb import ResponseEvaluationCost, recorded_response_dataset_cost


RUN_ROOT=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E053-r1')
JOBS={'PO_flip_T_pattern_preserving':'28574708','PO_to_M':'28574709'}
SACCT_COMMAND=['sacct','-j',','.join(JOBS.values()),'-X','--noheader','--parsable2',
               '--format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw,Partition']


def parse_scheduler(body):
    rows={}
    for line in body.splitlines():
        if not line.strip():
            continue
        values=line.strip().split('|')
        if len(values)!=6 or values[0] in rows:
            raise ValueError('unexpected or duplicate allocation-level scheduler row')
        job,state,exit_code,cores,elapsed,partition=values
        if job not in JOBS.values() or state!='COMPLETED' or exit_code!='0:0' or cores!='32' or partition!='hfacnormal01':
            raise ValueError('registered prior pilot allocation did not complete successfully on32HFCPUs')
        seconds=int(elapsed)
        if seconds<=0:
            raise ValueError('positive actual allocation duration required')
        rows[job]=dict(job_id=job,state=state,exit_code=exit_code,allocated_cpu_cores=int(cores),
                       allocation_seconds=seconds,partition=partition)
    if set(rows)!=set(JOBS.values()):
        raise ValueError('both actual completed prior allocations required')
    return rows


def export_costs(run_root, scheduler_body):
    root=Path(run_root).resolve()
    if root!=RUN_ROOT:
        raise ValueError('only the registered completed E053 pilot namespace is supported')
    production=root/'source/cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm'
    if sha256(production)!=PILOT_SCRIPT_SHA256:
        raise ValueError('registered MPI/thread1 physical driver changed')
    scheduling=parse_scheduler(scheduler_body)
    records=[]
    proofs=[]
    channels=[]
    for channel,job in JOBS.items():
        band=root/channel/'band'
        preflight=band/'vcneb_preflight.json'
        before=sha256(preflight)
        p=json.loads(preflight.read_text(encoding='utf-8'))
        if (p.get('n_images')!=9 or p.get('external_pressure_gpa')!=0.
                or p.get('climbing_image_requested') is not False
                or p.get('factory')!='examples.hfo2_fixed_input_factory:make_clamped_seed_cached_factory'
                or (p.get('mechanical_boundary') or {}).get('kind')!='clamped_plane'):
            raise ValueError('not the registered two-cached-endpoint seven-interior clamped pilot')
        channel_records=[]
        for i in range(1,8):
            image=band/f'image_{i:04d}'
            directories=sorted(image.glob('scf_*'))
            if len(directories)!=11 or any(not d.is_dir() for d in directories):
                raise ValueError('expected all11fresh evaluated frames per interior; cannot silently omit failures')
            for source in directories:
                audit_path=source/'call_audit.json'
                audit_sha=sha256(audit_path)
                audit=json.loads(audit_path.read_text(encoding='utf-8'))
                hashes=audit.get('input_sha256') or {}
                log=source/'OUT.ABACUS/running_scf.log'
                if (set(hashes)!={*CONTRACT,'STRU'} or any(hashes[k]!=v for k,v in CONTRACT.items())
                        or {k:sha256(source/k) for k in hashes}!=hashes
                        or sha256(log)!=audit.get('raw_log_sha256')):
                    raise ValueError('raw SCF input/log bytes differ from the recorded original physical contract')
                native=audited_results(source)  # complete E/F/stress, convergedSCF, oneDSIZE32
                if any(not np.allclose(native[k],audit.get('results',{}).get(k),atol=1e-12,rtol=0.)
                       for k in ('energy','forces','stress')):
                    raise ValueError('native raw results differ from the transport record')
                cost=ResponseEvaluationCost(evaluation_id='hf:'+str(source),raw_audit_sha256=audit_sha,
                    outcome='completed',elapsed_seconds=audit.get('elapsed_seconds'),
                    allocated_cpu_cores=scheduling[job]['allocated_cpu_cores'])
                if cost.elapsed_seconds is None or cost.elapsed_seconds<=0:
                    raise ValueError('this completed material cost export requires actual positive transport timing')
                if sha256(audit_path)!=audit_sha or sha256(log)!=audit['raw_log_sha256']:
                    raise ValueError('immutable completed record changed during inspection')
                records.append(cost);channel_records.append(cost)
                proofs.append(dict(evaluation_id=cost.evaluation_id,job_id=job,
                    raw_log_sha256=audit['raw_log_sha256'],STRU_sha256=hashes['STRU'],
                    six_original_physical_inputs_checked=True,native_MPI_ranks=32,native_complete_EFS_checked=True))
        if sha256(preflight)!=before:
            raise ValueError('completed pilot preflight changed during inspection')
        total=recorded_response_dataset_cost(channel_records,
            required_evaluation_ids=[r.evaluation_id for r in channel_records])
        allocation=scheduling[job]['allocation_seconds']*32/3600
        if total['recorded_transport_cpu_core_hours_sum']>allocation+1e-8:
            raise ValueError('serial SCF cost exceeds actual pilot allocation; review accounting')
        channels.append(dict(channel=channel,job_id=job,preflight_sha256=before,cost=total,
            scheduler_allocation_cpu_core_hours=allocation))
    total=recorded_response_dataset_cost(records,required_evaluation_ids=[r.evaluation_id for r in records])
    source_root=Path(__file__).resolve().parents[1]
    sources=['scripts/export_hfo2_scf_costs.py','vcneb/response_cost.py','scripts/audit_hfo2_static_replica.py']
    return dict(format_version=1,status='audited_completed_pilot_cost_not_prediction_cost',
        source_run_root=str(root),production_script_sha256=sha256(production),
        scheduler_command=SACCT_COMMAND,scheduler_raw_rows=scheduler_body,
        scheduler_records=list(scheduling.values()),records=[asdict(r) for r in records],
        per_call_proofs=proofs,channels=channels,cost=total,
        scheduler_allocation_cpu_core_hours=sum(r['allocation_seconds']*32/3600 for r in scheduling.values()),
        original_physical_contract_sha256=CONTRACT,
        export_source_sha256={s:sha256(source_root/s) for s in sources},
        new_DFT_calls=0,live_sources_read_or_modified=False,holdout_generated_or_read=False,
        scope='154fresh interior SCFs in two completed first pilots ONLY; no endpoint/preparation/continuation/probe costs',
        limitations=['Not the full study or complete G2 cost; earlier endpoint and preparation costs excluded explicitly',
                     'Transport wall time times allocation CPUs is a resource proxy, not measured CPU utilization',
                     'Job completion does not imply ordinary NEB convergence, TS identity or prediction accuracy'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('refusing an existing cost export artifact')
    if args.run_root.resolve()!=RUN_ROOT:
        raise ValueError('only the registered completed source namespace supported')
    scheduler=subprocess.run(SACCT_COMMAND,check=True,text=True,capture_output=True,timeout=30).stdout
    result=export_costs(args.run_root,scheduler)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8',newline='\n') as handle:
        json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(dict(status=result['status'],cost=result['cost'],new_DFT_calls=0)))


if __name__=='__main__':
    main()
