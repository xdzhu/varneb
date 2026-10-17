"""Read-only native P8 audit and immutable all-arm scoring; never calls DFT."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import subprocess
import numpy as np

EXEC=Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-G3-P8-E266-r1')
JOBS=['28935677','28935679']
CONFIG_SHA='2ca8fd47f8a53208501c046a0afe38c3d42f8792db0ac989c1139839b6ab22ae'
BUNDLE_SHA='100f1a55d4c18f54a82412b346d26085aa33f1c57e3f43bbb956ac95fff548dc'
FORECAST_SHA='d1d368372537a3d57d6a0cf38716ac8ea67210dbcac463d0dd2296a13e7b2046'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def accounting(raw):
    result={}
    for line in raw.splitlines():
        f=[v.strip() for v in line.split('|')]
        if (len(f)!=8 or f[0] not in JOBS or f[0] in result or f[1:3]!=['COMPLETED','0:0']
                or f[4]!='32' or f[6:]!=['hfacnormal01','1']):
            raise ValueError('both original jobs must complete normally at exact resources; no partial audit')
        result[f[0]]=dict(job_id=f[0],state=f[1],exit=f[2],elapsed_seconds=int(f[3]),node=f[5],CPUs=32,nodes=1,
                          allocation_core_hours=int(f[3])*32/3600)
    if set(result)!=set(JOBS):raise ValueError('both original accounting rows required')
    return [result[j] for j in JOBS]


def context(bundle):
    if sha(bundle/'config.json')!=CONFIG_SHA or sha(bundle/'bundle_files.json')!=BUNDLE_SHA or sha(bundle/'forecasts.json')!=FORECAST_SHA:
        raise ValueError('exact pre-label deployment and forecasts required')
    files=json.loads((bundle/'bundle_files.json').read_bytes())
    if len(files)!=6 or any(Path(n).name!=n or sha(bundle/n)!=h for n,h in files.items()):raise ValueError('runtime bundle drift')
    spec=importlib.util.spec_from_file_location('E266_native_post_worker',bundle/'worker.py')
    worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
    config=json.loads((bundle/'config.json').read_bytes());old,values=worker.validate(config,bundle,native=True)
    return worker,config,old,values


def audit(startup=False):
    raw=subprocess.check_output(['sacct','-X','-n','-P','-j',','.join(JOBS),
        '-o','JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,NodeList,Partition,NNodes'],text=True)
    allocations=None if startup else accounting(raw)
    bundle=EXEC/'bundle';worker,config,old,values=context(bundle)
    reader,w,data,_,_,free,_,_=values
    direction=np.asarray(json.loads((worker.OLD/'config.json').read_bytes())['branch_direction_36'])
    rows=[];all_native=[];all_sources=[]
    for shard,job in enumerate(JOBS):
        accepted=json.loads((EXEC/f'accepted_{shard}.json').read_bytes())
        if accepted['job_id']!=job or accepted['config_sha256']!=CONFIG_SHA:raise ValueError('original accepted handle differs')
        out=EXEC/f'shard_{shard}';targets=[p for p in config['points'] if p['id'] in worker.POLICY['shards'][shard]]
        if startup:
            point=targets[0];directory=out/f"point_{point['id']:03d}"
            paths=sorted(directory.glob('eval-*/result.json'))[:1]
            for path in paths:
                x=data['basis']@np.asarray(point['starting_internal_coordinates_A'])
                old.geometry(data,x,point['q_A']);record,source=reader.native_record(w,data,x,path,directory)
                all_native.append(record);all_sources.append(source)
            continue
        summary=json.loads((out/'completed.json').read_bytes())
        if (summary['job_id']!=job or summary['shard']!=shard or summary['config_sha256']!=CONFIG_SHA
                or summary['physical_plus005_generated_or_read'] or summary['new_training_candidates']!=0
                or [r['id'] for r in summary['rows']]!=[p['id'] for p in targets]):
            raise ValueError('complete original four-target summary required')
        for point,saved in zip(targets,summary['rows']):
            row_path=out/f"point_{point['id']:03d}.json";row=json.loads(row_path.read_bytes())
            if row!=saved or row['q_A']!=point['q_A'] or row['training_parent_id']!=point['training_parent_id']:
                raise ValueError('original target/parent and full producer row required')
            directory=out/f"point_{point['id']:03d}";contract=reader.physical_contract(data['chart'])
            budget=dict(schema_version=1,max_evaluation_attempts=12,physical_contract_sha256=contract['contract_sha256'])
            if (json.loads((directory/'contract.json').read_bytes())!=contract or json.loads((directory/'evaluation-budget.json').read_bytes())!=budget
                    or list(directory.glob('eval-*/failure.json'))):raise ValueError('exact contract, budget and all failure records required')
            markers=sorted(directory.glob('attempt-*.json'));records=[];sources=[]
            if not 1<=len(markers)<=12 or len(markers)!=row['reserved_attempts'] or len(list(directory.glob('eval-*')))!=len(markers):
                raise ValueError('every charged attempt must be accounted, no rollover')
            for i,marker in enumerate(markers):
                if marker.name!=f'attempt-{i:06d}.json':raise ValueError('contiguous reservation order required')
                reservation=json.loads(marker.read_bytes());x=np.asarray(reservation['coordinates'])
                old.geometry(data,x,point['q_A'])
                if i==0 and not np.array_equal(x,data['basis']@np.asarray(point['starting_internal_coordinates_A'])):
                    raise ValueError('literal frozen start must be evaluated first')
                name=reader.reserve_identity(reservation,dict(coordinates_full_A=x.tolist()),i,budget)
                path=directory/name/'result.json';lexical=Path(row['native_sources'][i]['result_path'])
                if lexical.resolve(strict=True)!=path.resolve(strict=True):raise ValueError('exact original native file required')
                record,source=reader.native_record(w,data,x,lexical,directory)
                if source!=row['native_sources'][i]:raise ValueError('native source differs from producer')
                records.append(record);sources.append(source)
            if len(records)!=row['new_native_SCFs']:raise ValueError('all native SCFs must be charged')
            label=None;final=row['result']
            if final is None:
                if row['status']!='attempt_budget_exhausted_no_retry' or len(records)!=12:raise ValueError('capped status differs')
            else:
                full=data['basis']@np.asarray(final['coordinates']);matched=[r for r in records if np.array_equal(r['coordinates'],full)]
                if len(matched)!=1 or matched[0]!=row['native_final'] or row['status']!=final['status']:raise ValueError('unique complete native final state required')
                g=data['basis'].T@np.asarray(matched[0]['gradient']);norm=float(np.linalg.norm(free.T@g));dg=float(direction@g)
                if (final['energy']!=matched[0]['enthalpy_eV'] or not np.isclose(norm,final['released_gradient_norm'],atol=2e-12,rtol=0)
                        or not np.isclose(dg,final['branch_gradient'],atol=2e-12,rtol=0)):
                    raise ValueError('full native energy and physical dual gradients differ')
                passed=norm<=.03 and abs(dg)<=1e-4
                if (row['status']=='gradient_candidate_stability_unchecked')!=passed:raise ValueError('both frozen stopping tests required')
                if passed:label=dict(energy_change_eV_cell=final['energy']-data['center_energy'],energy_eV_cell=final['energy'],
                    released_gradient_norm_eV_A=norm,zeta_gradient_eV_A=dg,stability_certified=False)
            rows.append(dict(id=point['id'],q_A=point['q_A'],training_parent_id=point['training_parent_id'],status=row['status'],
                native_SCFs=len(records),label=label,producer_sha256=sha(row_path)))
            all_native.extend(records);all_sources.extend(sources)
    rows.sort(key=lambda r:r['id'])
    if not startup and ([r['id'] for r in rows]!=list(range(1,9)) or len(all_native)>96):raise ValueError('all eight targets and96cap required')
    return dict(status='native_startup_prefix_not_complete_P8' if startup else 'complete_native_P8_candidate_audit_not_stability_certificate',
        config_sha256=CONFIG_SHA,forecast_sha256=FORECAST_SHA,raw_sacct=raw,allocations=allocations,rows=rows,
        native_records=all_native,native_sources=all_sources,native_SCFs_reparsed=len(all_native),
        valid_candidate_labels=sum(r['label'] is not None for r in rows),complete_cohort_audited=not startup,
        completed_call_proxy_core_hours=sum(r['validation']['elapsed_seconds'] for r in all_native)*32/3600,
        allocation_core_hours=None if startup else sum(r['allocation_core_hours'] for r in allocations),
        source_LF_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        original_2D_refusal_preserved=True,physical_plus005_or_private_arrays_accessed=False,new_DFT_or_scheduler_operations=0)


def error_metrics(errors):
    if not errors:return dict(n=0,max_absolute_meV_fu=None,RMSE_meV_fu=None)
    values=np.asarray(errors)*250.
    return dict(n=len(errors),max_absolute_meV_fu=float(np.max(np.abs(values))),RMSE_meV_fu=float(np.sqrt(np.mean(values**2))))


def score(forecasts,native):
    if (not native['complete_cohort_audited'] or native['forecast_sha256']!=FORECAST_SHA
            or [r['id'] for r in native['rows']]!=list(range(1,9))):raise ValueError('complete native P8 cohort required')
    labels={r['id']:r['label']['energy_change_eV_cell'] for r in native['rows'] if r['label'] is not None};models=[]
    for model in forecasts['new_models']:
        models.append((model['model_id'],'new_same_algorithm_target',[
            dict(id=p['id'],raw=p['prediction_raw_eV_cell'],aligned=p['prediction_origin_aligned_eV_cell'],reason=p['abstention_reason']) for p in model['predictions']]))
    standard=forecasts['original_standard_preserving_record']
    for name in standard['targets'][0]['models']:
        models.append((name,'retained_old_control_new_target_not_old_stable_plane_validation',[
            dict(id=p['id'],raw=p['models'][name]['value'],aligned=p['models'][name]['seen_origin_calibrated_energy_change_eV_cell'],
                 reason=p['models'][name].get('reason')) for p in standard['targets']]))
    for model in forecasts['original_reference_preserving_records']:
        models.append((f"{model['reference']}@{model['step_A']:.2f}",'retained_old_reference_new_target_not_old_stable_plane_validation',[
            dict(id=p['id'],raw=p['prediction_raw_eV_cell'],aligned=p['prediction_calibrated_eV_cell'],reason=p['reason']) for p in model['points']]))
    if len(models)!=23:raise ValueError('all four new and nineteen original arms required')
    rows=[]
    for name,kind,predictions in models:
        if [p['id'] for p in predictions]!=list(range(1,9)):raise ValueError('all eight forecast records required')
        ids=[p['id'] for p in predictions if p['raw'] is not None and p['id'] in labels]
        rows.append(dict(model=name,kind=kind,denominator=8,prediction_support=sum(p['raw'] is not None for p in predictions),
            valid_labels=len(labels),scored_ids=ids,model_abstention_ids=[p['id'] for p in predictions if p['raw'] is None],
            missing_label_ids=sorted(set(range(1,9))-set(labels)),
            points=[dict(id=p['id'],prediction_raw_eV_cell=p['raw'],prediction_origin_aligned_eV_cell=p['aligned'],
                observed_energy_change_eV_cell=labels.get(p['id']),abstention_reason=p['reason'],
                signed_raw_error_eV_cell=(None if p['raw'] is None or p['id'] not in labels else p['raw']-labels[p['id']]),
                signed_origin_aligned_error_eV_cell=(None if p['aligned'] is None or p['id'] not in labels else p['aligned']-labels[p['id']]))
                for p in predictions],
            raw=error_metrics([p['raw']-labels[p['id']] for p in predictions if p['id'] in ids]),
            origin_aligned=error_metrics([p['aligned']-labels[p['id']] for p in predictions if p['id'] in ids])))
    common=set(labels)
    for row in rows[:4]:common.intersection_update(row['scored_ids'])
    all_common=set(labels)
    for row in rows:all_common.intersection_update(row['scored_ids'])
    for row in rows:
        row['raw_on_primary_common_ids']=(error_metrics([p['signed_raw_error_eV_cell'] for p in row['points'] if p['id'] in common])
            if common.issubset(row['scored_ids']) else None)
        row['raw_on_all23_common_ids']=error_metrics([p['signed_raw_error_eV_cell'] for p in row['points'] if p['id'] in all_common])
    return dict(status='all23_frozen_arms_on_revised_P8_candidate_target_not_conditional_sheet_validation',arms=rows,
        primary_four_arm_common_ids=sorted(common),all23_common_ids=sorted(all_common),valid_labels=len(labels),denominator=8,
        reporting_target_meV_fu=1.,reporting_target_is_not_physical_error_bound=True,
        old_2D_refusal_preserved=True,no_refit_step_selection_or_rebias=True,
        stable_sheet_TS_resolved_advantage_or_JCTC_readiness_claimed=False,new_DFT_or_scheduler_operations=0,
        new_native_test_SCFs=native['native_SCFs_reparsed'],allocation_core_hours=native['allocation_core_hours'],all_prior_costs_additional=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('startup','audit','score'));p.add_argument('--output',type=Path,required=True)
    p.add_argument('--audit',type=Path);p.add_argument('--audit-sha256');args=p.parse_args()
    if args.output.exists():raise FileExistsError('fresh read-only report required')
    if args.mode=='score':
        if sha(args.audit)!=args.audit_sha256:raise ValueError('external complete audit pin required')
        forecasts=EXEC/'bundle/forecasts.json'
        if sha(forecasts)!=FORECAST_SHA:raise ValueError('unchanged pre-label forecasts required')
        result=score(json.loads(forecasts.read_bytes()),json.loads(args.audit.read_bytes()))
        result.update(source_audit_sha256=args.audit_sha256,forecast_sha256=FORECAST_SHA,
            source_LF_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest())
    else:result=audit(startup=args.mode=='startup')
    with args.output.open('x',encoding='utf-8',newline='\n') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps(dict(status=result['status'],native_SCFs=result.get('native_SCFs_reparsed'),valid_labels=result.get('valid_candidate_labels',result.get('valid_labels')))))
