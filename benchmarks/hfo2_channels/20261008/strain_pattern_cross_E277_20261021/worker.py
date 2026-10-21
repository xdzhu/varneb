"""Exactly two finite joint releases; no physical input changes or retries."""
from __future__ import annotations
import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import numpy as np
from scripts.audit_hfo2_static_replica import sha256
from vcneb import CalculatorModeEvaluator, relax_orthogonal_at_q, validate_path_geometry
from vcneb.mode_surface import _orthogonal_directions

ROOT=Path('/public/home/iai806/abacus/agent-runs/20261021-varneb-strain-pattern-E277-r1')
BASE=ROOT.with_name('20261018-varneb-matched-pattern-E270-r1')
HELPER=ROOT.with_name('20261018-varneb-matched-pattern-E271-r1')
BASE_PIN='f2aca3403e214fbad072ed53e69cd1afc496c560e52f2909aabe233dfb2b42ca'
HELPER_PIN='bcdb328562237d84522ddc5aa597b291ef5ba622716907e6692cd530cbda32f6'
NAMES=('zero_q1','plus1_q0')


def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def save(path,value):
    with Path(path).open('x',encoding='utf-8',newline='\n') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def check_bundle(root,pin):
    manifest=root/'bundle/bundle_files.json'
    if sha256(manifest)!=pin:raise ValueError('runtime manifest differs')
    for n,h in json.loads(manifest.read_bytes()).items():
        if Path(n).is_absolute() or '..' in Path(n).parts or sha256(root/'bundle'/n)!=h:
            raise ValueError('immutable source/input differs')


def validate(pin):
    for root,h in ((ROOT,pin),(BASE,BASE_PIN),(HELPER,HELPER_PIN)):check_bundle(root,h)
    study=load(BASE/'bundle/study.py','E277_fixed_study');study.validate(BASE_PIN)
    helper=load(HELPER/'bundle/repair.py','E277_release_helper')
    cfg=json.loads((ROOT/'bundle/prepared/config.json').read_bytes())
    prep=load(ROOT/'bundle/policy.py','E277_fixed_policy')
    if cfg['policy']!=prep.POLICY or cfg['new_DFT_calls']!=0:
        raise ValueError('only the registered two cross conditions admitted')
    return study,helper,cfg


def allowed(chart,x):
    if not np.isfinite(x).all() or np.linalg.norm(x)>1.:return False
    try:
        return validate_path_geometry([chart.reference_atoms,chart.to_atoms(x)],cell_scale=5.12968067458423,
            minimum_distance=1.6,maximum_deformation=.25)['valid']
    except (ValueError,np.linalg.LinAlgError):return False


class AttemptCap(RuntimeError):pass


class BoundedCache:
    def __init__(self,delegate):self.delegate=delegate;self.values={};self.fresh=0
    def __call__(self,x):
        x=np.asarray(x,dtype=float);key=x.tobytes()
        if key not in self.values:
            if self.fresh>=40:raise AttemptCap('finite forty-evaluation budget reached')
            self.fresh+=1;e,g=self.delegate(x);self.values[key]=(float(e),np.asarray(g).copy())
            print(json.dumps(dict(new_calls=self.fresh,energy_eV_cell=e,norm_A=float(np.linalg.norm(x)))),flush=True)
        e,g=self.values[key];return e,g.copy()


def release(chart,data,condition,helper,evaluate):
    start=np.asarray(condition['start_coordinates']);h=np.asarray(condition['numerical_hessian_39'])
    plane,frozen,c,t=helper.release_space(data,False,start)
    if not np.allclose(c.T@start,condition['constraint_offset_A'],atol=1e-12,rtol=0):
        raise ValueError('start fails cross-observable constraint')
    u=_orthogonal_directions(plane,frozen);np.linalg.cholesky(u.T@h@u)
    cache=BoundedCache(evaluate)
    try:
        point=relax_orthogonal_at_q(plane,[0.,0.],cache,frozen_directions=frozen,
            gradient_tolerance=.03,optimizer='safeguarded_bfgs',initial_trust_radius=.05,
            max_iterations=200,initial_hessian=h,trial_validator=lambda x:allowed(chart,x))
        x=np.asarray(point.coordinates)
        if (not np.allclose(c.T@x,condition['constraint_offset_A'],atol=1e-10,rtol=0)
                or np.max(np.abs(t.T@x))>1e-10 or not allowed(chart,x)):
            raise ValueError('candidate violates the held physical variables or geometry')
        candidate=dict(coordinates=x.tolist(),energy_eV_cell=point.energy,
            released_gradient_eV_A=point.orthogonal_gradient_norm,
            constraint_offset_A=(c.T@x).tolist(),coordinate_norm_A=float(np.linalg.norm(x)))
        status='gradient_candidate_stability_unchecked'
    except AttemptCap:status,candidate='attempt_budget_exhausted',None
    except RuntimeError as exc:
        if str(exc)!='no starting point reached the orthogonal-gradient tolerance':raise
        status,candidate='released_gradient_target_not_reached',None
    return dict(condition=condition['name'],status=status,candidate=candidate,new_calls=cache.fresh,
        full_stability_or_TS_or_global_barrier_certified=False)


def preflight(pin):
    study,helper,cfg=validate(pin);rows=[]
    for name in NAMES:
        folder=ROOT/'bundle/prepared'/name;chart=study.chart_for(folder/'seed')
        data=json.loads((folder/'seed/seed.json').read_bytes())
        condition=json.loads((folder/'condition.json').read_bytes());x=np.asarray(condition['start_coordinates'])
        study.cached_frozen(folder/'seed',chart)
        donor=condition['target_pattern'];diagonal=cfg['diagonals'][donor]
        # The existing diagonal is at this target plane, not the cross start.
        raw=ROOT/'bundle/prepared/diagonals'/donor
        native=json.loads((raw/'native/result.json').read_bytes())
        validator=load(ROOT/'bundle/native_start.py','E277_native_start')
        _,e,g=validator.cached_start(chart,raw/'native',data)
        if e!=diagonal['energy_eV_cell'] or native['coordinates']!=diagonal['coordinates']:
            raise ValueError('native diagonal cache differs')
        plane,frozen,c,t=helper.release_space(data,False,x)
        if (not allowed(chart,x) or not np.allclose(c.T@x,condition['constraint_offset_A'],atol=1e-12,rtol=0)
                or np.max(np.abs(t.T@x))>1e-10):raise ValueError('invalid cross-start geometry or observables')
        u=_orthogonal_directions(plane,frozen)
        np.linalg.cholesky(u.T@np.asarray(condition['numerical_hessian_39'])@u)
        rows.append(dict(condition=name,diagonal_energy_eV_cell=e,release_rank=u.shape[1],new_DFT_calls=0))
    result=dict(bundle_pin=pin,conditions=rows,new_DFT_calls=0)
    save(ROOT/'preflight.json',result);return result


def dispatch(mode,pin,name):
    study,_,_=validate(pin)
    proof=json.loads((ROOT/'preflight.json').read_bytes())
    if proof['bundle_pin']!=pin or proof['new_DFT_calls'] or len(proof['conditions'])!=2:raise ValueError('preflight required')
    if mode=='submit-held':
        if len(study.queue())>=2 or (ROOT/f'intent_{name}.json').exists():raise ValueError('slot occupied or prior intent; no duplicate')
        cmd=['sbatch','--hold','--parsable',f'--job-name=hfo2-E277-{name}',
             f'--export=ALL,RUN_DFT=1,BUNDLE_SHA256={pin},CONDITION={name}',
             f'--output={ROOT}/job-%j.out',str(ROOT/'bundle/run.slurm')]
        save(ROOT/f'intent_{name}.json',dict(command=cmd))
        reply=subprocess.check_output(cmd,text=True).strip();job=reply.split(';')[0]
        if not job.isdigit():raise ValueError('ambiguous submission; inspect accounting, never resubmit')
        result=dict(job_id=job,bundle_pin=pin,reply=reply);save(ROOT/f'held_{name}.json',result);return result
    held=json.loads((ROOT/f'held_{name}.json').read_bytes());job=held['job_id']
    raw=subprocess.check_output(['scontrol','show','job','-o',job],text=True)
    fields=dict(re.findall(r'(\S+)=(\S+)',raw))
    expected=dict(JobId=job,JobName=f'hfo2-E277-{name}',JobState='PENDING',Reason='JobHeldUser',
        Partition='hfacnormal01',NumCPUs='32',NumTasks='1',TimeLimit='02:00:00')
    expected['CPUs/Task']='32'
    if (held['bundle_pin']!=pin or len(study.queue())>2 or fields.get('NumNodes') not in ('1','1-1')
            or any(fields.get(k)!=v for k,v in expected.items())):raise ValueError('held resources differ')
    save(ROOT/f'accepted_{name}.json',dict(**held,raw_scontrol=raw))
    save(ROOT/f'release_intent_{name}.json',dict(job_id=job))
    reply=subprocess.check_output(['scontrol','release',job],text=True).strip()
    result=dict(job_id=job,bundle_pin=pin,reply=reply);save(ROOT/f'released_{name}.json',result);return result


def run(pin,name):
    study,helper,cfg=validate(pin);accepted=json.loads((ROOT/f'accepted_{name}.json').read_bytes())
    if (os.environ.get('RUN_DFT')!='1' or accepted['bundle_pin']!=pin
            or os.environ.get('SLURM_JOB_ID')!=accepted['job_id']
            or os.environ.get('SLURM_JOB_NAME')!=f'hfo2-E277-{name}'
            or os.environ.get('SLURM_JOB_PARTITION')!='hfacnormal01'
            or os.environ.get('SLURM_JOB_NUM_NODES')!='1' or os.environ.get('SLURM_NTASKS')!='1'
            or os.environ.get('SLURM_CPUS_PER_TASK')!='32' or len(study.queue())>2
            or any(os.environ.get(k)!='1' for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'))):
        raise ValueError('accepted two-allocation/32-MPI execution contract required')
    folder=ROOT/'bundle/prepared'/name;chart=study.chart_for(folder/'seed')
    data=json.loads((folder/'seed/seed.json').read_bytes());condition=json.loads((folder/'condition.json').read_bytes())
    output=ROOT/'results'/name;output.mkdir(parents=True,exist_ok=False)
    evaluator=CalculatorModeEvaluator(chart,
        study.make_factory(parameters={'source_directory':str(study.ASSETS)},command=study.COMMAND),
        workdir=output/'evaluations',calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
        minimum_distance_A=1.6,result_validator=study.result_validator,
        validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',max_evaluation_attempts=40)
    result=release(chart,data,condition,helper,evaluator)
    result.update(job_id=accepted['job_id'],bundle_pin=pin,policy=cfg['policy'])
    save(output/'complete.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('preflight','submit-held','release','run'))
    p.add_argument('--bundle-pin',required=True);p.add_argument('--condition',choices=NAMES)
    a=p.parse_args()
    if a.mode=='preflight':r=preflight(a.bundle_pin)
    elif a.condition is None:p.error('--condition required')
    elif a.mode=='run':r=run(a.bundle_pin,a.condition)
    else:r=dispatch(a.mode,a.bundle_pin,a.condition)
    if r is not None:print(json.dumps(r),flush=True)
