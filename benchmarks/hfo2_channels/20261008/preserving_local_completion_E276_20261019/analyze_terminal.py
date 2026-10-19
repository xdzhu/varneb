"""Independent local replay of every new E276 call and matched A/J candidates."""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import numpy as np

REPO=Path(__file__).resolve().parents[4];sys.path.insert(0,str(REPO))
from scripts.audit_hfo2_static_replica import sha256
CASE=Path(__file__).resolve().parent;OLD=CASE.with_name('preserving_pattern_branches_E272_20261018')


def module(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def restore(sign,dest):
    collector=module(CASE/'collect_terminal.py','E276_local_admission')
    folder=CASE/'completed'/sign;outer=json.loads((folder/'collection_receipt.json').read_bytes())
    archive=folder/'native_terminal.tar.gz'
    if (outer['job_id']!=collector.JOBS[sign] or outer['sign']!=sign or outer['bundle_pin']!=collector.PIN
            or sha256(archive)!=outer['archive_sha256'] or outer['new_DFT_calls'] or outer['scheduler_mutations']
            or outer['collector_source_sha256']!=sha256(CASE/'collect_terminal.py')):
        raise ValueError('original exact terminal export and collector identity required')
    collector.terminal_row('|'.join(outer['scheduler_row']),sign)
    with tarfile.open(archive) as tar:
        for m in tar.getmembers():
            p=Path(m.name)
            if (not p.parts or p.parts[0]!='native_terminal' or p.is_absolute() or '..' in p.parts
                    or not(m.isfile() or m.isdir()) or not(dest/p).resolve().is_relative_to(dest.resolve())):
                raise ValueError('unsafe terminal native archive member')
        tar.extractall(dest)
    root=dest/'native_terminal';files=json.loads((root/'files.json').read_bytes())
    if len(files)!=outer['exported_files']:raise ValueError('complete exported native cohort required')
    for n,h in files.items():
        p=root/n
        if not p.resolve().is_relative_to(root.resolve()) or sha256(p)!=h:raise ValueError('native export bytes differ')
    return root,outer


def read_one(sign,destination):
    root,receipt=restore(sign,destination)
    reader=module(OLD/'analyze_prefix.py','E276_native_replay')
    motion=module(OLD/'analyze_terminal.py','E276_actual_motion')
    chart,data,e0,_=reader.setup(2);complete=json.loads((root/'complete.json').read_bytes())
    cfg=json.loads((CASE/'prepared/config.json').read_bytes())
    if (complete['job_id']!=receipt['job_id'] or complete['sign']!=sign or complete['bundle_pin']!=receipt['bundle_pin']
            or complete['policy']!=cfg['policy'] or complete['frozen_energy_eV_cell']!=e0):
        raise ValueError('same native F and finite numerical revision required')
    q=np.asarray(data['orthonormal_39_by_4']);c=np.asarray(data['covectors_39_by_4'])
    t=np.zeros((39,3));t[:36]=np.tile(np.eye(3),(12,1))/np.sqrt(12)
    native=[]
    for path in sorted(root.glob('evaluations/eval-*/result.json')):
        raw=json.loads(path.read_bytes());x=np.asarray(raw['coordinates'])
        _,e,g,origin=reader.replay_native(root,chart,raw['coordinates'])
        if np.linalg.norm(x)>1. or np.max(np.abs(c.T@x))>1e-10 or np.max(np.abs(t.T@x))>1e-10:
            raise ValueError('actual trial violated the revised window or original physical observables')
        native.append(dict(coordinates=x.tolist(),energy_eV_cell=e,gradient=g,
            seconds=raw['validation']['elapsed_seconds'],native_result_sha256=origin['native_result_sha256'],
            norm_A=float(np.linalg.norm(x))))
    if len(native)!=receipt['new_SCF_calls'] or len(native)!=complete['new_calls'] or len(native)>50:
        raise ValueError('all original fresh calls must be retained, including failed arms')
    if sum(r['new_calls'] for r in complete['treatments'])+complete['known_direction_check']['new_calls']!=len(native):
        raise ValueError('all stage counts must close to the complete fresh native cohort')
    arms=[]
    for arm,raw in zip(('A','J'),complete['treatments']):
        if not 0<=raw['new_calls']<=(32 if arm=='A' else 16):raise ValueError('per-arm new-call cap differs')
        candidate=raw['candidate'];detail=None
        if candidate is not None:
            matches=[r for r in native if r['coordinates']==candidate['coordinates']]
            if len(matches)!=1:raise ValueError('one exact native candidate required')
            r=matches[0];x=np.asarray(r['coordinates']);g=r['gradient'];n=36 if arm=='A' else 39
            rg=g[:n]-q[:n]@(q[:n].T@g[:n])-t[:n]@(t[:n].T@g[:n]);norm=float(np.linalg.norm(rg))
            if (norm>.03 or not np.isclose(norm,candidate['released_gradient_norm_eV_A'],atol=2e-12,rtol=0)
                    or r['energy_eV_cell']!=candidate['energy_eV_cell'] or arm=='A' and np.any(x[36:])):
                raise ValueError('native candidate criterion, A cell or energy differs')
            detail=dict(energy_eV_cell=r['energy_eV_cell'],dE_minus_F_meV_fu=250*(r['energy_eV_cell']-e0),
                released_gradient_eV_A=norm,coordinates=x.tolist(),norm_A=float(np.linalg.norm(x)),
                native_result_sha256=r['native_result_sha256'],actual_motion=motion.candidate_motion(chart,x))
        elif raw['status'] not in ('attempt_budget_exhausted','released_gradient_target_not_reached','not_attempted_A_unavailable'):
            raise ValueError('unknown candidate failure status')
        arms.append(dict(arm=arm,status=raw['status'],new_calls=raw['new_calls'],candidate=detail))
    if arms[0]['candidate'] is None and (arms[1]['candidate'] is not None or complete['known_direction_check']['new_calls']):
        raise ValueError('J and directional checks cannot precede qualified A')
    check=complete['known_direction_check'];direction=None
    if arms[0]['candidate'] is not None and check['new_calls']==2:
        prep=json.loads((CASE/'prepared/seed/branch_preparation.json').read_bytes())
        v=np.asarray(prep['signed_starts'][1])/.1;xa=np.asarray(arms[0]['candidate']['coordinates'])
        points=[]
        for s,x in zip((-1,1),check['coordinates']):
            if not np.allclose(x,xa+s*.01*v,atol=2e-12,rtol=0):raise ValueError('known-direction probes changed')
            points.append(next(r for r in native if r['coordinates']==x))
        kg=float((points[1]['gradient']-points[0]['gradient'])@v/.02)
        ke=float((points[1]['energy_eV_cell']-2*arms[0]['candidate']['energy_eV_cell']+points[0]['energy_eV_cell'])/.01**2)
        if not np.isclose(kg,check['gradient_secant_eV_A2'],atol=2e-12,rtol=0) or not np.isclose(ke,check['energy_secant_eV_A2'],atol=2e-12,rtol=0):
            raise ValueError('independent raw directional secants differ')
        direction=dict(gradient_secant_eV_A2=kg,energy_secant_eV_A2=ke,full_stability_certified=False)
    seconds=sum(r['seconds'] for r in native)
    if (not np.isclose(seconds*32/3600,receipt['native_call_proxy_core_hours'],atol=1e-10,rtol=0)
            or not np.isclose(float(receipt['scheduler_row'][7])*32/3600,receipt['allocation_core_hours'],atol=1e-10,rtol=0)):
        raise ValueError('new cost accounting differs')
    return dict(sign=sign,job_id=receipt['job_id'],arms=arms,directional_secants_at_A=direction,
        archive_sha256=receipt['archive_sha256'],original_F_energy_eV_cell=e0,
        new_SCF_calls=len(native),all_new_native_EFS_reparsed_locally=True,
        native_proxy_core_hours=receipt['native_call_proxy_core_hours'],allocated_core_hours=receipt['allocation_core_hours'],
        previous_cached_starts_and_prior_costs_not_charged_again=True,
        full_orbital_pseudo_bytes_rechecked_locally=False,old_protocol_success_claimed=False,
        original_search_window_A=.5,new_search_window_A=1.,old_path_labels_and_scores_unchanged=True,
        certified_TS_or_full_conditional_stability=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();rows=[]
    for sign in ('minus','plus'):
        with tempfile.TemporaryDirectory(prefix='E276-terminal-') as td:rows.append(read_one(sign,Path(td)))
    a.output.mkdir(parents=True,exist_ok=False)
    result=dict(status='both_original_signed_continuations_replayed',rows=rows,new_DFT_calls=0,
        new_calls=sum(r['new_SCF_calls'] for r in rows),physical_inputs_changed=False,
        source_reader_sha256=sha256(__file__),
        independent_prediction_or_new_optimizer_advantage=False)
    with (a.output/'analysis.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(status=result['status'],new_calls=result['new_calls'])))
