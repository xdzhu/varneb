"""Two signed, bounded +1% A/check/J continuations; unchanged ABACUS physics."""
from __future__ import annotations
import argparse
import importlib.util
import json
import os
from pathlib import Path
import numpy as np
from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results, sha256
from ase.calculators.singlepoint import SinglePointCalculator
from vcneb import CalculatorModeEvaluator, relax_orthogonal_at_q, validate_path_geometry

ROOT=Path('/public/home/iai806/abacus/agent-runs/20261019-varneb-preserving-local-E276-r2')
BASE=ROOT.with_name('20261018-varneb-matched-pattern-E270-r1')
HELPER=ROOT.with_name('20261018-varneb-matched-pattern-E271-r1')
BASE_PIN='f2aca3403e214fbad072ed53e69cd1afc496c560e52f2909aabe233dfb2b42ca'
HELPER_PIN='bcdb328562237d84522ddc5aa597b291ef5ba622716907e6692cd530cbda32f6'


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def save(path,value):
    with Path(path).open('x',encoding='utf-8',newline='\n') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def check_files(root,pin):
    bundle=root/'bundle';manifest=bundle/'bundle_files.json'
    if sha256(manifest)!=pin:raise ValueError('immutable bundle manifest differs')
    for name,digest in json.loads(manifest.read_bytes()).items():
        if Path(name).is_absolute() or '..' in Path(name).parts or sha256(bundle/name)!=digest:
            raise ValueError('immutable runtime source or prepared input changed')


def validate_config(config):
    p=config['policy']
    if (config['new_DFT_calls'] or config['execution_admitted'] or p['strain']!=.01
            or p['retained']!=['Qx','Qz','Qy','polar_proxy_z']
            or p['maximum_coordinate_norm_A']!=1. or p['atomic_attempts']!=32 or p['joint_attempts']!=16
            or p['maximum_new_calls_per_sign']!=50 or p['gradient_tolerance_eV_A']!=.03
            or p['initial_trust_radius_A']!=.05
            or any(config['physical_input_sha256'].get(n)!=h for n,h in CONTRACT.items())
            or p['physical_inputs_changed'] or p['new_physical_Hessian']):
        raise ValueError('finite numerical revision or fixed physical contract differs')


def validate(pin):
    for root,value in ((ROOT,pin),(BASE,BASE_PIN),(HELPER,HELPER_PIN)):check_files(root,value)
    study=load(BASE/'bundle/study.py','E276_original_physics')
    helper=load(HELPER/'bundle/repair.py','E276_original_controls')
    study.validate(BASE_PIN)
    config=json.loads((ROOT/'bundle/prepared/config.json').read_bytes());validate_config(config)
    return study,helper,config


def cached_start(chart,folder,data):
    result=json.loads((folder/'result.json').read_bytes());v=result['validation'];raw=folder/'scf_000000'
    # STRU changes with geometry, unlike the six immutable physics assets.
    if (result['status']!='complete' or v['genuine_MPI_ranks']!=32
            or any(v['input_sha256'][n]!=h for n,h in CONTRACT.items())):
        raise ValueError('original cached start physics/MPI differs')
    if (any(sha256(raw/n)!=v['input_sha256'][n] for n in ('INPUT','KPT','STRU'))
            or sha256(raw/'OUT.ABACUS/running_scf.log')!=v['raw_log_sha256']
            or sha256(raw/'call_audit.json')!=v['call_audit_sha256']):
        raise ValueError('original cached start bytes differ')
    x=np.asarray(result['coordinates']);atoms=chart.to_atoms(x)
    if not same_ordered_geometry(atoms,read_fixed_hfo2_stru(raw/'STRU')):raise ValueError('cached start geometry differs')
    parsed=audited_results(raw);audit=json.loads((raw/'call_audit.json').read_bytes())
    if any(not np.allclose(parsed[k],audit['results'][k],atol=1e-12,rtol=0) for k in parsed):raise ValueError('native EFS differs')
    atoms.calc=SinglePointCalculator(atoms,**parsed);energy,g=chart.energy_and_gradient(x,atoms)
    if energy!=result['enthalpy_eV'] or not np.allclose(g,result['gradient'],atol=2e-12,rtol=0):raise ValueError('native start gradient differs')
    return x,energy,g


class AttemptCap(RuntimeError):pass


class Cache:
    def __init__(self,delegate,known):
        self.delegate=delegate;self.values={np.asarray(x,dtype=float).tobytes():(float(e),np.asarray(g).copy()) for x,e,g in known}
        self.fresh=0;self.stage_start=0;self.stage_cap=0;self.stage='not_started'
    def begin(self,name,cap):self.stage=name;self.stage_start=self.fresh;self.stage_cap=cap
    def __call__(self,x):
        x=np.asarray(x,dtype=float);key=x.tobytes()
        if key not in self.values:
            if self.fresh-self.stage_start>=self.stage_cap or self.fresh>=50:raise AttemptCap('finite new-call budget exhausted')
            self.fresh+=1;e,g=self.delegate(x);self.values[key]=(float(e),np.asarray(g).copy())
            print(json.dumps(dict(stage=self.stage,new_calls=self.fresh,norm_A=float(np.linalg.norm(x)),energy_eV_cell=float(e))),flush=True)
        e,g=self.values[key];return e,g.copy()


def guard(study,chart,x):
    if not np.isfinite(x).all() or np.linalg.norm(x)>1.:return False
    try:
        report=validate_path_geometry([chart.reference_atoms,chart.to_atoms(x)],cell_scale=study.SCALE,
            minimum_distance=1.6,maximum_deformation=.25);return report['valid']
    except (ValueError,np.linalg.LinAlgError):return False


def arm(study,helper,chart,data,cache,hessian,start,name):
    if start is None:return dict(treatment=name,status='not_attempted_A_unavailable',candidate=None,new_calls=0)
    atomic=name=='A';cache.begin(name,32 if atomic else 16)
    plane,frozen,c,t=helper.release_space(data,atomic,np.asarray(start))
    def full(x):return np.r_[x,np.zeros(3)] if atomic else np.asarray(x)
    try:
        point=relax_orthogonal_at_q(plane,[0.,0.],lambda x:study.pullback(cache,full(x),atomic),
            frozen_directions=frozen,gradient_tolerance=.03,optimizer='safeguarded_bfgs',
            initial_trust_radius=.05,max_iterations=200,
            initial_hessian=hessian[:36,:36] if atomic else hessian,
            trial_validator=lambda x:guard(study,chart,full(x)))
        x=full(point.coordinates)
        if np.max(np.abs(c.T@x))>1e-10 or np.max(np.abs(t.T@x))>1e-10 or atomic and np.any(x[36:]):
            raise ValueError('retained physical coordinates, translations or A cell changed')
        candidate=dict(coordinates=x.tolist(),energy_eV_cell=point.energy,
            released_gradient_norm_eV_A=point.orthogonal_gradient_norm,physical_constraint_delta=(c.T@x).tolist(),
            coordinate_norm_A=float(np.linalg.norm(x)),stationary_TS_or_full_stability_certified=False)
        status='gradient_candidate_stability_unchecked'
    except AttemptCap:status,candidate='attempt_budget_exhausted',None
    except RuntimeError as exc:
        if str(exc)!='no starting point reached the orthogonal-gradient tolerance':raise
        status,candidate='released_gradient_target_not_reached',None
    return dict(treatment=name,status=status,candidate=candidate,new_calls=cache.fresh-cache.stage_start)


def check_direction(study,chart,cache,candidate,direction):
    if candidate is None:return dict(status='A_unavailable',new_calls=0,full_stability_certified=False)
    x=np.asarray(candidate['coordinates']);points=[x+s*.01*direction for s in (-1,1)]
    if not all(guard(study,chart,p) for p in points):return dict(status='guard_prevents_probes',new_calls=0,full_stability_certified=False)
    cache.begin('known_direction',2);values=[cache(p) for p in points]
    kg=float((values[1][1]-values[0][1])@direction/.02)
    ke=float((values[1][0]-2*candidate['energy_eV_cell']+values[0][0])/.01**2)
    return dict(status='known_direction_positive' if min(kg,ke)>0 else 'known_direction_not_positive_or_unresolved',
        gradient_secant_eV_A2=kg,energy_secant_eV_A2=ke,probe_step_A=.01,
        coordinates=[p.tolist() for p in points],new_calls=cache.fresh-cache.stage_start,full_stability_certified=False)


def preflight(pin):
    study,helper,cfg=validate(pin);seed=ROOT/'bundle/prepared/seed';chart=study.chart_for(seed)
    data=json.loads((seed/'seed.json').read_bytes());e0,g0=study.cached_frozen(seed,chart);rows=[]
    for sign in ('minus','plus'):
        folder=ROOT/'bundle/prepared/starts'/sign;x,e,g=cached_start(chart,folder,data)
        recorded=next(r for r in cfg['starts'] if r['sign']==sign)
        if sha256(folder/'result.json')!=recorded['source_result_sha256']:
            raise ValueError('original signed last-result identity differs')
        if not .49<np.linalg.norm(x)<.5 or not guard(study,chart,x):raise ValueError('original native last trial required')
        _,_,c,t=helper.release_space(data,True,x)
        if np.max(np.abs(c.T@x))>1e-10 or np.max(np.abs(t.T@x))>1e-10 or np.any(x[36:]):raise ValueError('original controls changed')
        rows.append(dict(sign=sign,coordinates=x.tolist(),energy_eV_cell=e,full_gradient_norm_eV_A=float(np.linalg.norm(g))))
    result=dict(bundle_pin=pin,new_DFT_calls=0,starts=rows,frozen_energy_eV_cell=e0,held_observables_unchanged=True)
    save(ROOT/'preflight.json',result);return result


def run(pin,sign):
    study,helper,cfg=validate(pin);accepted=json.loads((ROOT/f'accepted_{sign}.json').read_bytes())
    proof=json.loads((ROOT/'preflight.json').read_bytes())
    if (accepted['bundle_pin']!=pin or proof['bundle_pin']!=pin or proof['new_DFT_calls']
            or len(study.queue())>2 or os.environ.get('RUN_DFT')!='1'
            or os.environ.get('SLURM_JOB_ID')!=accepted['job_id']
            or os.environ.get('SLURM_JOB_NAME')!='hfo2-E276-'+sign
            or os.environ.get('SLURM_JOB_PARTITION')!='hfacnormal01'
            or os.environ.get('SLURM_JOB_NUM_NODES')!='1' or os.environ.get('SLURM_NTASKS')!='1'
            or os.environ.get('SLURM_CPUS_PER_TASK')!='32'
            or any(os.environ.get(n)!='1' for n in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'))):
        raise ValueError('accepted held allocation, preflight or resource cap differs')
    seed=ROOT/'bundle/prepared/seed';chart=study.chart_for(seed);data=json.loads((seed/'seed.json').read_bytes())
    prep=json.loads((seed/'branch_preparation.json').read_bytes());h=np.asarray(prep['numerical_hessian_39'])
    e0,g0=study.cached_frozen(seed,chart);x,e,g=cached_start(chart,ROOT/'bundle/prepared/starts'/sign,data)
    out=ROOT/'results'/sign;out.mkdir(parents=True,exist_ok=False)
    evaluator=CalculatorModeEvaluator(chart,study.make_factory(parameters={'source_directory':str(study.ASSETS)},command=study.COMMAND),
        workdir=out/'evaluations',calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
        minimum_distance_A=1.6,result_validator=study.result_validator,
        validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',max_evaluation_attempts=50)
    cache=Cache(evaluator,[(np.zeros(39),e0,g0),(x,e,g)])
    a=arm(study,helper,chart,data,cache,h,x,'A');save(out/'treatment_A.json',a)
    direction=np.asarray(prep['signed_starts'][1])/.1
    check=check_direction(study,chart,cache,a['candidate'],direction);save(out/'known_direction.json',check)
    j=arm(study,helper,chart,data,cache,h,None if a['candidate'] is None else np.asarray(a['candidate']['coordinates']),'J')
    save(out/'treatment_J.json',j)
    result=dict(job_id=accepted['job_id'],sign=sign,bundle_pin=pin,frozen_energy_eV_cell=e0,
        treatments=[a,j],known_direction_check=check,new_calls=cache.fresh,policy=cfg['policy'],
        source_start_energy_eV_cell=e,original_budget_failure_retained=True,old_labels_unchanged=True,
        status='finite_signed_continuation_finished_not_TS_or_global_minimum')
    save(out/'complete.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('preflight','run'))
    p.add_argument('--bundle-pin',required=True);p.add_argument('--sign',choices=('minus','plus'))
    a=p.parse_args()
    if a.mode=='preflight':print(json.dumps(preflight(a.bundle_pin)))
    elif a.sign is None:p.error('--sign is required for run')
    else:run(a.bundle_pin,a.sign)
