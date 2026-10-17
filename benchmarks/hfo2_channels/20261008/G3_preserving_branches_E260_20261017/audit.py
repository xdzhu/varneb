"""Read-only complete native audit of the fixed ten training candidates."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import argparse
import hashlib
import importlib.util
import json
import subprocess
import numpy as np

EXEC=Path('/public/home/iai806/abacus/agent-runs/20261017-varneb-G3-P-branches-E260-r2')
JOBS=['28934048','28934049']
CONFIG_SHA='63b1aa0eea88e7cd7c2f1a78ce9f4f3c7f39e16a9caf04d3f6d8c6ef5f84ab80'
MANIFEST_SHA='4422902d6474afe1056181b978f6098c9fc1284ff6cc66d11ff865d499d87eff'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def accounting(raw):
    rows={}
    for line in raw.splitlines():
        f=[v.strip() for v in line.split('|')]
        if (len(f)!=8 or f[0] not in JOBS or f[0] in rows
                or f[1:3]!=['COMPLETED','0:0'] or f[4]!='32'
                or f[6:]!=['hfacnormal01','1']):
            raise ValueError('both original jobs must complete normally; live/failed data remain pending audit')
        rows[f[0]]=dict(job_id=f[0],state=f[1],exit=f[2],elapsed_seconds=int(f[3]),
            CPUs=32,nodes=1,node=f[5],allocation_core_hours=int(f[3])*32/3600)
    if set(rows)!=set(JOBS):raise ValueError('both original job accounting rows required')
    return [rows[j] for j in JOBS]


def native():
    raw=subprocess.check_output(['sacct','-X','-n','-P','-j',','.join(JOBS),
        '-o','JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,NodeList,Partition,NNodes'],text=True)
    allocations=accounting(raw)
    bundle=EXEC/'bundle'
    if sha(bundle/'bundle_files.json')!=MANIFEST_SHA or sha(bundle/'config.json')!=CONFIG_SHA:
        raise ValueError('exact original bundle/config required')
    files=json.loads((bundle/'bundle_files.json').read_bytes())
    if len(files)!=8 or any(Path(n).name!=n or sha(bundle/n)!=h for n,h in files.items()):
        raise ValueError('immutable original source differs')
    spec=importlib.util.spec_from_file_location('E260_original_worker',bundle/'worker.py')
    worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
    config=json.loads((bundle/'config.json').read_bytes())
    reader,w,data,_,_,free,_,_=worker.native_context(config,bundle)
    direction=np.asarray(config['branch_direction_36'])
    all_native=[];rows=[];producers=[]
    for shard,job in enumerate(JOBS):
        accepted=json.loads((EXEC/f'accepted_{shard}.json').read_bytes())
        if accepted['job_id']!=job or accepted['config_sha256']!=CONFIG_SHA:
            raise ValueError('exact original accepted handle required')
        out=EXEC/f'shard_{shard}';summary=json.loads((out/'completed.json').read_bytes())
        targets=[p for p in config['points'] if p['id'] in config['policy']['shards'][shard]]
        if (summary['job_id']!=job or summary['shard']!=shard or summary['config_sha256']!=CONFIG_SHA
                or summary['new_independent_labels']!=0 or summary['physical_plus005_generated_or_read']
                or [r['id'] for r in summary['rows']]!=[p['id'] for p in targets]):
            raise ValueError('complete unchanged shard identity required')
        for point,stored in zip(targets,summary['rows']):
            point_path=out/f"point_{point['id']:03d}.json"
            row=json.loads(point_path.read_bytes())
            if row!=stored or row['q_A']!=point['q_A'] or row['start_sign']!=point['start_sign']:
                raise ValueError('literal original point/sign and complete producer row required')
            directory=out/f"point_{point['id']:03d}"
            contract=reader.physical_contract(data['chart'])
            budget=dict(schema_version=1,max_evaluation_attempts=12,physical_contract_sha256=contract['contract_sha256'])
            if (json.loads((directory/'contract.json').read_bytes())!=contract
                    or json.loads((directory/'evaluation-budget.json').read_bytes())!=budget
                    or list(directory.glob('eval-*/failure.json'))):
                raise ValueError('original native contract, budget and all failure records required')
            markers=sorted(directory.glob('attempt-*.json'));records=[];sources=[]
            if not 1<=len(markers)<=12 or len(markers)!=row['reserved_attempts']:
                raise ValueError('every charged attempt must be present, without rollover')
            for i,marker in enumerate(markers):
                if marker.name!=f'attempt-{i:06d}.json':raise ValueError('contiguous reservation order required')
                reservation=json.loads(marker.read_bytes());x=np.asarray(reservation['coordinates'])
                worker.geometry(data,x,point['q_A'])
                if i==0 and not np.array_equal(x,data['basis']@np.asarray(point['starting_internal_coordinates_A'])):
                    raise ValueError('the literal frozen initial start must be evaluated first')
                name=reader.reserve_identity(reservation,dict(coordinates_full_A=x.tolist()),i,budget)
                expected=directory/name/'result.json';literal=Path(row['native_sources'][i]['result_path'])
                if literal.resolve(strict=True)!=expected.resolve(strict=True):
                    raise ValueError('exact recorded native source required')
                record,source=reader.native_record(w,data,x,literal,directory)
                if source!=row['native_sources'][i]:raise ValueError('producer source differs from native audit')
                records.append(record);sources.append(source)
            if len(records)!=row['new_native_SCFs'] or len(list(directory.glob('eval-*')))!=len(records):
                raise ValueError('all evaluations and preparation costs required')
            final=row['result'];label=None
            if final is not None:
                if final['status']!=row['status']:raise ValueError('producer final status differs')
                x=np.asarray(final['coordinates']);full=data['basis']@x
                matched=[r for r in records if np.array_equal(r['coordinates'],full)]
                if len(matched)!=1 or matched[0]!=row['native_final']:
                    raise ValueError('unique exact final native evaluation required')
                g=data['basis'].T@np.asarray(matched[0]['gradient'])
                norm=float(np.linalg.norm(free.T@g));dg=float(direction@g)
                zeta=float(direction@(x-data['basis'].T@np.asarray(data['origin_record']['coordinates'])))
                if (final['energy']!=matched[0]['enthalpy_eV']
                        or not np.isclose(norm,final['released_gradient_norm'],atol=2e-12,rtol=0)
                        or not np.isclose(dg,final['branch_gradient'],atol=2e-12,rtol=0)
                        or not np.isclose(zeta,row['zeta_final_A'],atol=2e-12,rtol=0)):
                    raise ValueError('native final energy/dual-gradient/amplitude differs')
                converged=norm<=.03 and abs(dg)<=1e-4
                if (final['status']=='gradient_candidate_stability_unchecked')!=converged:
                    raise ValueError('both original stop tests required')
                label=dict(energy_eV_cell=final['energy'],zeta_A=zeta,
                    released_gradient_norm_eV_A=norm,zeta_gradient_eV_A=dg,
                    gradient_candidate=converged,coordinates_internal_A=x.tolist())
            elif row['status']!='attempt_budget_exhausted_no_retry':
                raise ValueError('missing training candidate outcome differs')
            rows.append(dict(id=point['id'],Q_id=point['Q_id'],q_A=point['q_A'],
                start_sign=point['start_sign'],status=row['status'],label=label,
                attempted_SCFs=len(markers),native_sources=sources,producer_sha256=sha(point_path)))
            all_native.extend(records)
        producers.append(dict(path=str(out/'completed.json'),sha256=sha(out/'completed.json')))
    rows.sort(key=lambda r:r['id'])
    if [r['id'] for r in rows]!=list(range(10)) or sum(r['attempted_SCFs'] for r in rows)>120:
        raise ValueError('all ten original starts and120 cap required')
    pairs=[]
    for i in range(5):
        plus,minus=rows[2*i:2*i+2];a,b=plus['label'],minus['label']
        valid=a is not None and b is not None and a['gradient_candidate'] and b['gradient_candidate']
        pairs.append(dict(Q_id=i,q_A=plus['q_A'],both_gradient_candidates=valid,
            plus_minus_energy_meV_per_fu=(250*(a['energy_eV_cell']-b['energy_eV_cell']) if valid else None),
            plus_minus_internal_distance_A=(float(np.linalg.norm(np.asarray(a['coordinates_internal_A'])-
                b['coordinates_internal_A'])) if valid else None),stable_distinct_branches_certified=False))
    return dict(status='complete_native_ten_training_branch_trial_not_stability_certificate',
        observed_CST=datetime.now(timezone(timedelta(hours=8))).isoformat(),raw_sacct=raw,
        config_sha256=CONFIG_SHA,bundle_manifest_sha256=MANIFEST_SHA,
        audit_LF_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        allocations=allocations,rows=rows,pairs=pairs,native_records=all_native,producer_sources=producers,
        new_native_SCFs_reparsed=len(all_native),training_candidates_processed=10,
        gradient_candidates=sum(r['label'] is not None and r['label']['gradient_candidate'] for r in rows),
        completed_call_proxy_core_hours=sum(r['validation']['elapsed_seconds'] for r in all_native)*32/3600,
        allocation_core_hours=sum(r['allocation_core_hours'] for r in allocations),all_prior_costs_additional=True,
        training_committed=dict(T=34,preserving=49),new_independent_labels=0,
        physical_plus005_or_private_arrays_accessed=False,new_DFT_or_scheduler_writes=0,
        original_2D_refusal_preserved=True,TS_stability_or_predictive_advantage_certified=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('fresh separate audit output required')
    result=native()
    with a.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:result[k] for k in ('status','new_native_SCFs_reparsed','gradient_candidates')}))
