"""Collect original E276 handles only after terminal state; no scheduler writes."""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import numpy as np
from scripts.audit_hfo2_static_replica import sha256

ROOT=Path('/public/home/iai806/abacus/agent-runs/20261019-varneb-preserving-local-E276-r2')
PIN='16bd55d10c00d579026f6aff11e4b472db07d4526c177e0f0e8a6a0ae43422d1'
JOBS={'minus':'28954771','plus':'28954772'}


def save(path,data):
    with Path(path).open('x',encoding='utf-8',newline='\n') as f:
        json.dump(data,f,indent=2,allow_nan=False);f.write('\n')


def terminal_row(raw,sign):
    rows=[r.split('|') for r in raw.strip().splitlines() if r.split('|')[0]==JOBS[sign]]
    if (len(rows)!=1 or len(rows[0])<12
            or rows[0][1:7]!=['COMPLETED','0:0','32','hfacnormal01','iai806','1']
            or rows[0][10]!='hfo2-E276-'+sign):
        raise ValueError('exact original completed handle required; live, missing or timeout is not a result')
    return rows[0]


def collect(sign):
    accepted=json.loads((ROOT/f'accepted_{sign}.json').read_bytes())
    if accepted['job_id']!=JOBS[sign] or accepted['sign']!=sign or accepted['bundle_pin']!=PIN:
        raise ValueError('original admitted handle and runtime pin required')
    raw=subprocess.check_output(['sacct','-n','-X','-P','-j',JOBS[sign],
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User,AllocNodes,ElapsedRaw,Start,End,JobName%40,NodeList'],text=True)
    row=terminal_row(raw,sign)
    spec=importlib.util.spec_from_file_location('E276_terminal_worker',ROOT/'bundle/worker.py')
    w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
    study,helper,cfg=w.validate(PIN)
    seed=ROOT/'bundle/prepared/seed';chart=study.chart_for(seed);data=json.loads((seed/'seed.json').read_bytes())
    folder=ROOT/'results'/sign;complete=json.loads((folder/'complete.json').read_bytes())
    if complete['job_id']!=JOBS[sign] or complete['sign']!=sign or complete['bundle_pin']!=PIN or complete['policy']!=cfg['policy']:
        raise ValueError('actual finite terminal marker required')
    calls=[];seconds=[]
    for path in sorted(folder.glob('evaluations/eval-*/result.json')):
        native=json.loads(path.read_bytes());x,e,g=w.cached_start(chart,path.parent,data)
        calls.append(dict(relative_record=path.relative_to(folder).as_posix(),native_result_sha256=sha256(path),
            coordinate_norm_A=float(np.linalg.norm(x)),energy_eV_cell=e,
            seconds=native['validation']['elapsed_seconds']))
        seconds.append(native['validation']['elapsed_seconds'])
    if (len(calls)!=complete['new_calls'] or not 0<=len(calls)<=50
            or any(not np.isfinite(s) or s<=0 for s in seconds)
            or sum(r['new_calls'] for r in complete['treatments'])+complete['known_direction_check']['new_calls']!=len(calls)):
        raise ValueError('complete native call accounting differs')
    for arm,result in zip(('A','J'),complete['treatments']):
        candidate=result['candidate']
        if candidate is None:continue
        matches=[p for p in folder.glob('evaluations/eval-*/result.json')
            if json.loads(p.read_bytes())['coordinates']==candidate['coordinates']]
        if len(matches)!=1:raise ValueError('accepted candidate must identify one fresh native source')
        x,e,g=w.cached_start(chart,matches[0].parent,data)
        q=np.asarray(data['orthonormal_39_by_4']);t=np.zeros((39,3));t[:36]=np.tile(np.eye(3),(12,1))/np.sqrt(12)
        n=36 if arm=='A' else 39;released=g[:n]-q[:n]@(q[:n].T@g[:n])-t[:n]@(t[:n].T@g[:n])
        if (e!=candidate['energy_eV_cell'] or np.linalg.norm(released)>.03
                or not np.isclose(np.linalg.norm(released),candidate['released_gradient_norm_eV_A'],atol=2e-12,rtol=0)
                or np.max(np.abs(np.asarray(data['covectors_39_by_4']).T@x))>1e-10
                or np.linalg.norm(x)>1. or arm=='A' and np.any(x[36:])):
            raise ValueError('native candidate, retained observables or released-gradient criterion differs')
    dest=ROOT/'exports'/sign;dest.mkdir(parents=True,exist_ok=False);out=dest/'native_terminal';out.mkdir()
    for name in ('complete.json','treatment_A.json','treatment_J.json','known_direction.json'):
        shutil.copyfile(folder/name,out/name)
    for call in calls:
        source=folder/call['relative_record'];target=out/call['relative_record']
        target.parent.mkdir(parents=True)
        for name in ('result.json','scf_000000/INPUT','scf_000000/KPT','scf_000000/STRU',
                     'scf_000000/OUT.ABACUS/running_scf.log','scf_000000/call_audit.json',
                     'scf_000000/abacus.out','scf_000000/abacus.err'):
            src=source.parent/name;dst=target.parent/name
            if src.is_file():dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
    receipt=dict(job_id=JOBS[sign],sign=sign,bundle_pin=PIN,scheduler_row=row,new_SCF_calls=len(calls),
        native_call_proxy_core_hours=sum(seconds)*32/3600,allocation_core_hours=float(row[7])*32/3600,
        all_new_native_EFS_replayed_on_HF=True,original_assets_checked_by_fixed_runtime=True,
        old_cached_source_calls_not_charged_again=True,prior_costs_additional=True,
        collector_source_sha256=sha256(__file__),calls=calls,new_DFT_calls=0,scheduler_mutations=0)
    save(out/'receipt.json',receipt)
    files={p.relative_to(out).as_posix():sha256(p) for p in out.rglob('*') if p.is_file()};save(out/'files.json',files)
    archive=dest/'native_terminal.tar.gz'
    with tarfile.open(archive,'w:gz') as tar:tar.add(out,arcname='native_terminal')
    report={k:v for k,v in receipt.items() if k!='calls'}
    report.update(archive_sha256=sha256(archive),archive_bytes=archive.stat().st_size,exported_files=len(files))
    save(dest/'collection_receipt.json',report);return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sign',choices=tuple(JOBS),required=True)
    print(json.dumps(collect(p.parse_args().sign)))
