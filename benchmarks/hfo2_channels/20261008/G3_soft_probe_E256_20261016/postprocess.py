"""Native read-only four-point audit and offline signed energy/force analysis."""
from pathlib import Path
from datetime import datetime, timezone, timedelta
import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
sys.dont_write_bytecode = True
import numpy as np

EXEC = Path('/public/home/iai806/abacus/agent-runs/20261016-varneb-G3-soft4-E256-r1')
CONFIG_SHA = 'f7dc1053bff96bbd141d0a7e2e5c14d949f91f69832852d61bf7af9b932de7a7'
MANIFEST_SHA = 'ed53bc642ad6d7c1e3255a3d27a0648d8c0d93ae33b1e60687cf59b416dbc8b2'
SOFT_SHA = '1e5dcd01ec5c7731646c9e1f6b0f230fc8bcf45bbcad5a785125534a0d2ee065'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def pinned(path, digest):
    if sha(path) != digest: raise ValueError('frozen source bytes differ')
    return json.loads(path.read_bytes())


def checked_literal_path(expected, namespace, source):
    """Prove exact resolved identity, then keep the recorded provenance spelling."""
    path=Path(source['result_path'])
    resolved=path.resolve(strict=True)
    if resolved!=expected.resolve(strict=True) or not resolved.is_relative_to(namespace.resolve(strict=True)):
        raise ValueError('recorded native source is not the exact expected file')
    return path


def native_audit():
    if any(os.environ.get(k) != '1' for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS')):
        raise ValueError('single-thread native replay required')
    bundle=EXEC/'bundle'
    files=pinned(bundle/'bundle_files.json',MANIFEST_SHA)
    if set(files)!={'worker.py','run.slurm','config.json','protocol.md','soft_direction.json'} or any(
            sha(bundle/n)!=h for n,h in files.items()):
        raise ValueError('original immutable bundle differs')
    config=pinned(bundle/'config.json',CONFIG_SHA)
    account=subprocess.run(['sacct','-n','-P','-X','-j','28912687',
        '--format=JobIDRaw,State%30,ExitCode,ElapsedRaw,AllocCPUS,NodeList,Partition,NNodes'],
        check=True,capture_output=True,text=True).stdout.strip()
    f=[x.strip() for x in account.split('|')]
    if len(f)!=8 or f[:3]!=['28912687','COMPLETED','0:0'] or f[4]!='32' or f[6:]!=['hfacnormal01','1']:
        raise ValueError('actual normal one-node32CPU completion required')
    spec=importlib.util.spec_from_file_location('E257_original_four_probe_worker',bundle/'worker.py')
    w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
    reader,old,data,_,_,_=w.native_context(config,bundle)
    directory=EXEC/'probes/evaluations'
    contract=reader.physical_contract(data['chart'])
    budget=dict(schema_version=1,max_evaluation_attempts=4,
        physical_contract_sha256=contract['contract_sha256'])
    if (json.loads((directory/'contract.json').read_bytes())!=contract
            or json.loads((directory/'evaluation-budget.json').read_bytes())!=budget
            or len(list(directory.glob('attempt-*.json')))!=4
            or len(list(directory.glob('eval-*/result.json')))!=4
            or list(directory.glob('eval-*/failure.json'))):
        raise ValueError('original four-reservation contract/coverage differs')
    completion=json.loads((EXEC/'probes/completed.json').read_bytes())
    if len(completion['native_sources'])!=4:raise ValueError('four recorded native sources required')
    records=[];sources=[]
    for i,p in enumerate(config['points']):
        marker=directory/f'attempt-{i:06d}.json'
        name=reader.reserve_identity(json.loads(marker.read_bytes()),p,i,budget)
        literal=checked_literal_path(directory/name/'result.json',directory,completion['native_sources'][i])
        record,source=reader.native_record(old,data,np.asarray(p['coordinates_full_A']),literal,directory)
        records.append(record);sources.append(source)
    if (completion['status']!='four_native_training_statics_complete_not_branch_or_TS'
            or completion['job_id']!='28912687' or completion['config_sha256']!=CONFIG_SHA
            or completion['new_native_SCFs']!=4 or completion['new_center_SCFs']!=0
            or completion['records']!=records or completion['native_sources']!=sources
            or not completion['original_independent_admission_refusal_preserved']
            or completion['independent_points_computed']!=0
            or completion['physical_plus005_generated_or_read']):
        raise ValueError('producer summary differs from native records or declared scope')
    seconds=sum(r['validation']['elapsed_seconds'] for r in records)
    return dict(status='complete_native_four_training_probe_audit_not_stability_certificate',
        observed_CST=datetime.now(timezone(timedelta(hours=8))).isoformat(),
        job_id='28912687',accounting=dict(raw=account,elapsed_seconds=int(f[3]),CPUs=32,nodes=1,
            node_list=f[5],normal_completion=True,allocation_core_hours=int(f[3])*32/3600),
        config_sha256=CONFIG_SHA,bundle_manifest_sha256=MANIFEST_SHA,
        original_worker_sha256=files['worker.py'],completion_sha256=sha(EXEC/'probes/completed.json'),
        audit_text_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        old_center_SCFs_reparsed=1,new_probe_SCFs_reparsed=4,
        base_record=data['origin_record'],records=records,native_sources=sources,
        completed_call_proxy_core_hours=seconds*32/3600,
        earlier_training_derivatives_and_preparation_costs_additional=True,
        new_DFT_or_scheduler_writes=0,independent_points_computed=0,
        original_independent_admission_refusal_preserved=True,
        physical_plus005_generated_or_read=False,
        physical_error_bound_stability_TS_or_independent_accuracy_certified=False)


def analyze(audit, config, soft):
    if (audit['status']!='complete_native_four_training_probe_audit_not_stability_certificate'
            or audit['config_sha256']!=CONFIG_SHA or audit['new_probe_SCFs_reparsed']!=4
            or not audit['accounting']['normal_completion']
            or len(audit['records'])!=4 or len(config['points'])!=4
            or [p['t_A'] for p in config['points']]!=[.05,-.05,.10,-.10]):
        raise ValueError('actual complete original four-point cohort required')
    direction=np.asarray(soft['soft_direction_full40'])
    base=audit['base_record'];e0=base['enthalpy_eV'];g0=float(direction@base['gradient'])
    rows=[]
    for p,r in zip(config['points'],audit['records']):
        if not np.array_equal(p['coordinates_full_A'],r['coordinates']):
            raise ValueError('native coordinate order differs')
        delta=r['enthalpy_eV']-e0
        rows.append(dict(id=p['id'],t_A=p['t_A'],energy_eV_cell=r['enthalpy_eV'],
            delta_eV_cell=delta,delta_meV_cell=delta*1000,delta_meV_per_fu=delta*250,
            gradient_along_direction_eV_A=float(direction@r['gradient']),
            original_two_quadratic_delta_eV_cell=p['quadratic_delta_eV_cell'],
            original_two_quadratic_errors_eV_cell=[value-delta for value in p['quadratic_delta_eV_cell']]))
    pairs=[]
    for h in (.05,.10):
        plus=next(r for r in rows if r['t_A']==h);minus=next(r for r in rows if r['t_A']==-h)
        de=plus['energy_eV_cell']-minus['energy_eV_cell']
        gp=plus['gradient_along_direction_eV_A'];gm=minus['gradient_along_direction_eV_A']
        work=h*(gp+gm)
        pairs.append(dict(amplitude_A=h,
            energy_central_curvature_eV_A2=(plus['delta_eV_cell']+minus['delta_eV_cell'])/(h*h),
            gradient_central_curvature_eV_A2=(gp-gm)/(2*h),
            energy_odd_slope_eV_A=de/(2*h),
            paired_energy_difference_eV_cell=de,endpoint_trapezoid_work_eV_cell=work,
            energy_minus_endpoint_work_eV_cell=de-work))
    return dict(status='four_training_statics_signed_energy_force_diagnostic_not_independent_test',
        job_id=audit['job_id'],formula_units_per_cell=4,base_energy_eV_cell=e0,
        base_soft_gradient_eV_A=g0,rows=rows,pairs=pairs,
        predictions_refitted=False,amplitudes_or_signs_selected_after_labels=False,
        physical_error_bound_stable_branch_TS_or_independent_accuracy_certified=False,
        original_independent_admission_refusal_preserved=True,
        independent_points_computed=0,physical_plus005_generated_or_read=False,
        allocation_core_hours=audit['accounting']['allocation_core_hours'],
        completed_call_proxy_core_hours=audit['completed_call_proxy_core_hours'],
        earlier_training_derivatives_and_preparation_costs_additional=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('native','analyze'))
    p.add_argument('--output',type=Path,required=True);p.add_argument('--audit',type=Path)
    p.add_argument('--audit-sha256');p.add_argument('--bundle',type=Path,default=Path(__file__).parent)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('fresh report required')
    if a.mode=='native':
        if a.output.resolve().is_relative_to(EXEC/'bundle') or a.output.resolve().is_relative_to(EXEC/'probes'):
            raise ValueError('audit output must stay outside immutable bundle/native producer')
        result=native_audit()
    else:
        result=analyze(pinned(a.audit,a.audit_sha256),pinned(a.bundle/'config.json',CONFIG_SHA),
            pinned(a.bundle/'soft_direction.json',SOFT_SHA))
    with a.output.open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(status=result['status'],new_DFT=0)))
