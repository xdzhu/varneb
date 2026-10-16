"""Read-only native T-test audit and fixed E185/E188 scoring; never run DFT."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import argparse
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
sys.dont_write_bytecode=True
import numpy as np

CONFIG_SHA='d78a39c60a6b7d2af52ab8a9cfabfe829635f4fcd94d6c682e0f63c98826506e'
MANIFEST_SHA='ef491d966740b3d82b970596e706d6701dc9673e197ee838b10dc0c0ca98c368'
ACCEPTED_SHA=('65d78497130e1cdadadd9925f113512399f7a10c1a4e0ea96738070327d30e56',
              '8b4f14ba647045a8ac5797d10750e270656c09e254db160f5c471124f5b18a9a')
JOBS=('28913247','28913248')
EXEC=Path('/public/home/iai806/abacus/agent-runs/20261017-varneb-G3-T8-E258-r1')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def pinned(path,digest):
    if sha(path)!=digest:raise ValueError('strict source SHA256 differs')
    return json.loads(Path(path).read_bytes())


def accounting(raw):
    rows=[]
    for line in raw.splitlines():
        f=[s.strip() for s in line.split('|')]
        if len(f)!=8 or f[0] not in JOBS:raise ValueError('unparsed original allocation')
        if f[1:3]!=['COMPLETED','0:0'] or f[4]!='32' or f[6:]!=['hfacnormal01','1']:
            raise ValueError('both original one-node32CPU jobs must complete normally')
        if not f[3].isdigit():raise ValueError('finite allocation seconds required')
        rows.append(dict(job_id=f[0],state=f[1],exit=f[2],elapsed_seconds=int(f[3]),CPUs=32,
            node=f[5],partition=f[6],nodes=1,allocation_core_hours=int(f[3])*32/3600))
    if sorted(r['job_id'] for r in rows)!=list(JOBS):raise ValueError('both original accounting rows required')
    return sorted(rows,key=lambda r:r['job_id'])


def native_audit(bundle):
    files=pinned(bundle/'bundle_files.json',MANIFEST_SHA)
    if len(files)!=9 or any(sha(bundle/n)!=h for n,h in files.items()):raise ValueError('original nine-file bundle drift')
    config=pinned(bundle/'config.json',CONFIG_SHA)
    raw=subprocess.check_output(['sacct','-X','-n','-P','-j',','.join(JOBS),
        '-o','JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,NodeList,Partition,NNodes'],text=True)
    allocations=accounting(raw)  # refuse partial/live data, never resubmit
    spec=importlib.util.spec_from_file_location('unchanged_E258_worker',bundle/'worker.py')
    worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
    reader,w,data=worker.native_context(config,bundle)
    rows=[];all_native=[];all_sources=[];producer_sources=[]
    for shard,job in enumerate(JOBS):
        accepted=pinned(EXEC/f'accepted_{shard}.json',ACCEPTED_SHA[shard])
        if accepted['job_id']!=job or accepted['config_sha256']!=CONFIG_SHA:raise ValueError('original accepted handle differs')
        out=EXEC/f'shard_{shard}';summary_path=out/'completed.json'
        summary=json.loads(summary_path.read_bytes())
        if (summary['job_id']!=job or summary['shard']!=shard or summary['config_sha256']!=CONFIG_SHA
            or summary['preserving_test_points_computed']!=0 or summary['physical_plus005_generated_or_read']):
            raise ValueError('original complete producer identity differs')
        targets=[p for p in config['points'] if p['id'] in config['policy']['shards'][shard]]
        if [r['id'] for r in summary['rows']]!=[p['id'] for p in targets]:raise ValueError('all four original shard targets required')
        for point,stored in zip(targets,summary['rows']):
            point_path=out/f"point_{point['id']:03d}.json";row=json.loads(point_path.read_bytes())
            if row!=stored or row['q_A']!=point['q_A'] or row['original_start_full_A']!=point['starting_full_coordinates_A']:
                raise ValueError('exact original target/start and complete producer row required')
            directory=out/f"point_{point['id']:03d}"
            contract=reader.physical_contract(data['chart'])
            budget=dict(schema_version=1,max_evaluation_attempts=12,physical_contract_sha256=contract['contract_sha256'])
            if (json.loads((directory/'contract.json').read_bytes())!=contract
                or json.loads((directory/'evaluation-budget.json').read_bytes())!=budget
                or list(directory.glob('eval-*/failure.json'))):raise ValueError('original native contract/budget or failure differs')
            markers=sorted(directory.glob('attempt-*.json'));paths=sorted(directory.glob('eval-*'))
            if not 1<=len(markers)<=12 or len(paths)!=len(markers) or row['reserved_attempts']!=len(markers):
                raise ValueError('all exact charged attempts required; no omitted failure or rollover')
            native=[];sources=[]
            for i,marker in enumerate(markers):
                if marker.name!=f'attempt-{i:06d}.json':raise ValueError('ordered original reservations required')
                reservation=json.loads(marker.read_bytes());x=np.asarray(reservation['coordinates'])
                name=reader.reserve_identity(reservation,dict(coordinates_full_A=x.tolist()),i,budget)
                path=directory/name/'result.json'
                old_source=row['native_sources'][i];lexical=Path(old_source['result_path'])
                if lexical.resolve(strict=True)!=path.resolve(strict=True):raise ValueError('strict same native file required')
                record,source=reader.native_record(w,data,x,lexical,directory)
                if source!=old_source:raise ValueError('producer native source differs from raw replay')
                native.append(record);sources.append(source)
            if row['new_native_SCFs']!=len(native):raise ValueError('all native SCF costs required')
            label=None;final=row['point']
            if final is not None:
                if row['status']!='single_start_gradient_candidate_stability_unchecked':raise ValueError('candidate status differs')
                x=np.asarray(final['coordinates_full_A']);internal=np.asarray(final['coordinates_internal_A'])
                matched=[r for r in native if np.array_equal(r['coordinates'],x)]
                if len(matched)!=1 or matched[0]!=row['native_candidate']:raise ValueError('unique complete final native candidate required')
                if (not np.allclose(data['basis']@internal,x,atol=2e-12,rtol=0)
                    or not np.allclose(data['axes'].T@internal,point['q_A'],atol=2e-12,rtol=0)
                    or not worker.faces.canary.geometry_guard(data['chart'],data['basis'],internal)):
                    raise ValueError('original Q/embedding/geometry guard differs')
                g=data['basis'].T@np.asarray(matched[0]['gradient'])
                orthogonal=g-data['axes']@(data['axes'].T@g);norm=float(np.linalg.norm(orthogonal))
                if not np.isclose(norm,final['orthogonal_gradient_norm_eV_A'],atol=2e-12,rtol=0) or norm>.03:
                    raise ValueError('original released-gradient candidate condition differs')
                if final['energy_eV_cell']!=matched[0]['enthalpy_eV']:raise ValueError('native candidate energy differs')
                label=dict(energy_eV_cell=final['energy_eV_cell'],released_gradient_norm_eV_A=norm,
                    coordinates_full_A=x.tolist(),native_source=sources[native.index(matched[0])])
            elif row['status']!='optimizer_not_converged_no_retry':raise ValueError('missing-label status differs')
            rows.append(dict(id=point['id'],q_A=point['q_A'],status=row['status'],label=label,
                reserved_attempts=len(markers),new_native_SCFs=len(native),
                producer_point_sha256=sha(point_path),native_sources=sources))
            all_native.extend(native);all_sources.extend(sources)
        producer_sources.append(dict(path=str(summary_path),sha256=sha(summary_path)))
    if [r['id'] for r in rows]!=list(range(1,9)) or sum(r['reserved_attempts'] for r in rows)>96:
        raise ValueError('all eight original targets and96 cap required')
    return dict(status='complete_native_original_T_test_cohort_not_global_stability',
        observed_CST=datetime.now(timezone(timedelta(hours=8))).isoformat(),
        audit_LF_code_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        config_sha256=CONFIG_SHA,bundle_manifest_sha256=MANIFEST_SHA,allocations=allocations,
        raw_sacct=raw,rows=rows,original_targets_processed=8,valid_gradient_candidate_labels=sum(r['label'] is not None for r in rows),
        new_native_SCFs_reparsed=len(all_native),native_records=all_native,native_sources=all_sources,producer_sources=producer_sources,
        completed_call_proxy_core_hours=sum(r['validation']['elapsed_seconds'] for r in all_native)*32/3600,
        allocation_core_hours=sum(r['allocation_core_hours'] for r in allocations),all_prior_costs_additional=True,
        preserving_independent_labels=0,total_original_independent_inventory=16,
        physical_plus005_or_private_author_arrays_accessed=False,new_DFT_or_scheduler_writes=0,
        full_stability_TS_or_whole_path_predictive_gain_certified=False)


def models(config):
    arms=[]
    for name in config['points'][0]['models']:
        arms.append(dict(name=name,kind='E185_standard',origin_bias_eV_cell=config['points'][0]['models'][name]['origin_bias_eV_cell'],points=[dict(id=p['id'],q_A=p['q_A'],
            raw_eV_cell=p['models'][name]['value'],calibrated_eV_cell=p['models'][name]['seen_origin_calibrated_energy_change_eV_cell'],
            formal_raw_eV_cell=p['models'][name].get('formal_energy_change_eV_cell',p['models'][name]['value']),
            reason=p['models'][name].get('reason')) for p in config['points']]))
    for r in config['reference_predictions']:
        arms.append(dict(name=f"{r['reference']}@{r['step_A']:.2f}",kind='E188_restricted_reference',origin_bias_eV_cell=r['origin_bias_eV_cell'],points=[dict(
            id=p['id'],q_A=p['q_A'],raw_eV_cell=p['prediction_raw_eV_cell'],
            calibrated_eV_cell=p['prediction_calibrated_eV_cell'],formal_raw_eV_cell=p['formal_energy_change_eV_cell'],reason=p['reason']) for p in r['points']]))
    return arms


def metrics(errors_eV):
    if not errors_eV:return None
    if not all(math.isfinite(x) for x in errors_eV):raise ValueError('finite signed errors required')
    value=np.asarray(errors_eV)*1000
    return dict(n=len(value),mean_signed_meV_cell=float(np.mean(value)),
        max_absolute_meV_cell=float(np.max(np.abs(value))),RMSE_meV_cell=float(np.sqrt(np.mean(value**2))),
        max_absolute_meV_per_fu=float(np.max(np.abs(value))/4),RMSE_meV_per_fu=float(np.sqrt(np.mean(value**2))/4))


def score_labels(config,labels):
    """Labels are energy changes relative to the original raw peak, not Q0."""
    if not set(labels)<=set(range(1,9)):raise ValueError('only original T target IDs allowed')
    if not all(math.isfinite(v) for v in labels.values()):raise ValueError('finite target labels required')
    arms=models(config);supported=[set(p['id'] for p in a['points'] if p['calibrated_eV_cell'] is not None) for a in arms]
    shared=sorted(set.intersection(*supported)&set(labels))
    full_support=[i for i,s in enumerate(supported) if len(s)==8]
    rows=[]
    for arm,support in zip(arms,supported):
        detail=[]
        for p in arm['points']:
            truth=labels.get(p['id']);pred=p['calibrated_eV_cell'];raw=p['raw_eV_cell']
            detail.append(dict(**p,observed_change_eV_cell=truth,
                signed_calibrated_error_eV_cell=None if truth is None or pred is None else pred-truth,
                signed_raw_error_eV_cell=None if truth is None or raw is None else raw-truth))
        available=[p for p in detail if p['signed_calibrated_error_eV_cell'] is not None]
        rows.append(dict(name=arm['name'],kind=arm['kind'],origin_bias_eV_cell=arm['origin_bias_eV_cell'],points=detail,denominator=8,
            prediction_support=len(support),valid_labels=len(labels),scored_ids=[p['id'] for p in available],
            model_abstention_ids=sorted(set(range(1,9))-support),missing_label_ids=sorted(set(range(1,9))-set(labels)),
            calibrated=metrics([p['signed_calibrated_error_eV_cell'] for p in available]),
            uncalibrated=metrics([p['signed_raw_error_eV_cell'] for p in detail if p['signed_raw_error_eV_cell'] is not None]),
            common_all_arms=metrics([p['signed_calibrated_error_eV_cell'] for p in available if p['id'] in shared])))
    return dict(status='original_T_local_independent_errors_not_whole_path_gain',arms=rows,
        all_19_arm_common_ids=shared,full_prediction_support_arm_names=[arms[i]['name'] for i in full_support],
        full_support_common_label_ids=sorted(labels),observed_origin_energy_change_eV_cell=config['standard_predictions']['observed_origin_energy_change_eV_cell'],
        T_denominator=8,preserving_valid_labels=0,preserving_denominator=8,
        total_original_independent_inventory=16,reporting_target_meV_per_fu=1.,reporting_target_is_not_physical_error_bound=True,
        no_refit_rebias_step_selection_or_branch_repair=True,no_statistical_replicates_from_reference_gauges=True,
        unmeasured_DFT_error_or_stability_resolution=True,resolved_advantage_TS_or_JCTC_readiness_claimed=False)


def score(config,audit):
    if (audit['status']!='complete_native_original_T_test_cohort_not_global_stability'
        or audit['config_sha256']!=CONFIG_SHA or audit['bundle_manifest_sha256']!=MANIFEST_SHA
        or [r['id'] for r in audit['rows']]!=list(range(1,9))):raise ValueError('complete native original T cohort required')
    peak=config['standard_predictions']['raw_peak_reference_energy_eV_cell']
    labels={r['id']:r['label']['energy_eV_cell']-peak for r in audit['rows'] if r['label'] is not None}
    result=score_labels(config,labels)
    result.update(new_native_SCFs=audit['new_native_SCFs_reparsed'],
        completed_call_proxy_core_hours=audit['completed_call_proxy_core_hours'],allocation_core_hours=audit['allocation_core_hours'],
        all_prior_costs_additional=True,source_config_sha256=CONFIG_SHA)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('native','score'))
    p.add_argument('--bundle',type=Path,default=Path(__file__).parent);p.add_argument('--audit',type=Path)
    p.add_argument('--audit-sha256');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('fresh report required')
    if a.mode=='native':
        if a.output.resolve().is_relative_to(EXEC/'bundle') or a.output.resolve().is_relative_to(EXEC/'shard_0') or a.output.resolve().is_relative_to(EXEC/'shard_1'):
            raise ValueError('read-only audit output must remain outside producer and immutable bundle')
        result=native_audit(a.bundle)
    else:
        audit=pinned(a.audit,a.audit_sha256);result=score(pinned(a.bundle/'config.json',CONFIG_SHA),audit)
        result['source_native_audit_sha256']=a.audit_sha256
    with a.output.open('x',encoding='utf-8',newline='\n') as out:json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps(dict(status=result['status'],new_DFT=0)))
