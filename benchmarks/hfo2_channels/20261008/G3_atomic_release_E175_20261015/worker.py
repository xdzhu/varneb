"""Bounded B3 atomic-only release at the existing two centers/five Q coordinates.

Cell is fixed to each frozen-Q seed, not to one cell shared across Q. This
is a release-treatment ablation, not a new reference, surface or TS search.
"""
from __future__ import annotations
import argparse,hashlib,importlib.util,json,os
from copy import deepcopy
from pathlib import Path
import subprocess,sys
sys.dont_write_bytecode=True
if (Path(__file__).resolve().parents[4]/'vcneb').is_dir():
    sys.path.insert(0,str(Path(__file__).resolve().parents[4]))
import numpy as np
from scripts import hfo2_G3_conditional_faces as faces
from vcneb.mode_surface import _safeguarded_bfgs

study=faces.study
OLD=study.RUNS/'20261015-varneb-G3-faces-E161-r1/source'
OLD_EXEC=study.RUNS/'20261015-varneb-G3-feasible-execution-E171-r1/bundle'
RUNTIME_SHA='d07bb011c9e313afe23ca4e78207356ff98adfc43ca231bdd047e6d5d0f1be20'
CONTROLS_SHA='f49c662d705d4733dbcf1d78bb4e0260a0165df4c0a0df0af3e468f56fc8ffc8'
ADMISSION_SHA='8f0635eecad789677d0313c759924f2665a5f05186b334e9cb84ef8afde37d3f'
POLICY=dict(reservation_id='E175-five-atomic-release-training-nodes-per-existing-face',
    channels=['T','preserving'],points_per_face=5,maximum_new_attempts_per_node=12,
    maximum_new_attempts_per_face=60,maximum_new_attempts_total=120,
    gradient_tolerance_eV_A=.03,initial_trust_radius_A=.05,maximum_iterations=200,
    maximum_internal_radius_A=.5,minimum_distance_A=1.6,
    extra_training_nodes_per_face=5,training_committed_after=dict(T=34,preserving=35),
    training_remaining_after=dict(T=47,preserving=46),independent_slots_untouched_per_face=8,
    original_Hessian_centers_only=True,pressure_GPa=0.,strain=0.,ordinary_NEB_fmax_eV_A=.10,
    job_name='hfo2-G3-E175-atoms',partition='hfacnormal01',nodes=1,MPI_ranks=32,
    allocation_hours=5,maximum_committed_study_allocations=2,climb=False)

sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def atomic_release_basis(basis,axes):
    """Exact zero cell rows in the measured block-ordered internal chart."""
    b,a=np.asarray(basis),np.asarray(axes)
    if (b.shape!=(40,36) or a.shape!=(36,2) or not np.isfinite(b).all()
            or not np.isfinite(a).all() or np.any(b[36:,:33]) or np.any(b[:36,33:])
            or not np.allclose(b.T@b,np.eye(36),atol=2e-12,rtol=0)
            or not np.allclose(a.T@a,np.eye(2),atol=2e-12,rtol=0)):
        raise ValueError('original orthonormal block-ordered atomic/open-cell chart required')
    _,s,vh=np.linalg.svd(a[:33].T,full_matrices=True)
    if len(s)!=2 or min(s)<=1e-8:raise ValueError('two independent atomic Q constraints required')
    r=np.zeros((36,31));r[:33]=vh[2:].T
    return r


def load_controls(bundle):
    p=bundle/'controls_summary.json';n=bundle/'admission_complete_native_preflight.json'
    if sha(p)!=CONTROLS_SHA or sha(n)!=ADMISSION_SHA:
        raise ValueError('exact completed eight-control/native-admission evidence required')
    controls=json.loads(p.read_bytes());native=json.loads(n.read_bytes())
    if (controls['new_SCFs']!=8 or controls['job_id']!='28865656'
            or native['controls']['sources']!=controls['rows']
            or native['controls']['original_native_SCFs_reparsed']!=8):
        raise ValueError('all eight original native frozen controls required')
    return controls


def prepare(bundle):
    controls=load_controls(bundle);out={}
    for channel in POLICY['channels']:
        d=faces.load_inputs(channel);r=atomic_release_basis(d['basis'],d['axes'])
        for block in d['curvature']['blocks']:
            raw=np.asarray(block['raw_projected_eV_A2'])[:36,:36];h=(raw+raw.T)/2
            np.linalg.cholesky(r.T@h@r)  # seed gate, not post-relaxation stability
        points=[dict(id=0,q_A=[0.,0.],frozen_internal_A=np.zeros(36).tolist(),frozen_full_A=np.zeros(40).tolist())]
        for p in [v for v in controls['rows'] if v['channel']==channel]:
            expected=d['basis']@(d['axes']@p['q_A'])
            if not np.allclose(expected,p['coordinates_full_A'],atol=2e-12,rtol=0):raise ValueError('original frozen seed differs')
            points.append(dict(id=p['point_id'],q_A=p['q_A'],frozen_internal_A=(d['axes']@p['q_A']).tolist(),
                frozen_full_A=p['coordinates_full_A']))
        if [p['id'] for p in points]!=list(range(5)):raise ValueError('original axial point order required')
        out[channel]=dict(points=points,atomic_release_basis_36_by_31=r.tolist(),
            old_input_pins=list(faces.PINS[channel]))
    return dict(schema_version=1,purpose='same_Q_atomic_only_release_training_ablation',policy=deepcopy(POLICY),
        channels=out,controls_summary_sha256=CONTROLS_SHA,native_admission_sha256=ADMISSION_SHA,
        six_input_sha256=deepcopy(study.CONTRACT),binary_sha256=study.BINARY_SHA,old_runtime_sha256=RUNTIME_SHA,
        worker_text_sha256=study.text_sha(Path(__file__)),wrapper_text_sha256=study.text_sha(bundle/'run.slurm'),
        physical_reserved_geometry_generated_or_read=False,independent_prediction_or_released_stability_certified=False)


def validate(config,bundle):
    expected=prepare(bundle)
    for channel in POLICY['channels']:
        a=np.asarray(config['channels'][channel]['atomic_release_basis_36_by_31'])
        b=np.asarray(expected['channels'][channel]['atomic_release_basis_36_by_31'])
        if (a.shape!=(36,31) or not np.isfinite(a).all() or np.any(a[33:])
                or not np.allclose(a.T@a,np.eye(31),atol=2e-12,rtol=0)
                or not np.allclose(a@a.T,b@b.T,atol=2e-12,rtol=0)):
            raise ValueError('declared pure-atomic fixed-Q release span differs')
        expected['channels'][channel]['atomic_release_basis_36_by_31']=a.tolist()
        points=config['channels'][channel]['points']
        if len(points)!=5:raise ValueError('exact original five frozen seeds required')
        for point,ref in zip(points,expected['channels'][channel]['points']):
            base=np.asarray(point['frozen_internal_A'])
            if (base.shape!=(36,) or not np.isfinite(base).all() or
                    not np.allclose(base,ref['frozen_internal_A'],atol=2e-12,rtol=0)):
                raise ValueError('declared internal affine anchor differs')
            ref['frozen_internal_A']=base.tolist()
    if json.dumps(expected,sort_keys=True)!=json.dumps(config,sort_keys=True):
        raise ValueError('source/physics/reservation/five-point contract differs')
    return {c:faces.load_inputs(c) for c in POLICY['channels']}


def scheduler_snapshot(own=None):
    q=subprocess.run(['squeue','-h','-u','iai806','-p','hfacnormal01','-o','%i|%100j|%T|%C|%D'],
        check=True,capture_output=True,text=True).stdout
    committed=[]
    for line in q.splitlines():
        f=[v.strip() for v in line.split('|')]
        if len(f)!=5:raise ValueError('unparsed queue; no admission')
        if f[1].startswith('hfo2-G3-'):
            if f[3:]!=['32','1']:raise ValueError('study allocation resources differ')
            committed.append(dict(job_id=f[0],name=f[1],state=f[2],CPUs=32,nodes=1))
    ready=len(committed)<2 if own is None else (len(committed)<=2 and
        any(r['job_id']==own and r['name']==POLICY['job_name'] for r in committed))
    return dict(ready=ready,committed=committed,own_job=own,scheduler_writes=0)


def preflight(config,bundle):
    data=validate(config,bundle)
    if study.ROOT.resolve()!=OLD.resolve():raise ValueError('original immutable E161 source required')
    files=faces.canary.pinned(OLD/'runtime_files.json',RUNTIME_SHA)
    if len(files)!=380 or any(sha(OLD/n)!=h for n,h in files.items()):raise ValueError('old380runtime drift')
    if sha(OLD_EXEC/'bundle_files.json')!='d4f7a2fbf1950001b725eed0ff1d286e98825c23d62d9a070ac5bf6d4aaa00d2':
        raise ValueError('existing native-auditor bundle drift')
    oldfiles=json.loads((OLD_EXEC/'bundle_files.json').read_bytes())
    if any(sha(OLD_EXEC/n)!=h for n,h in oldfiles.items()):raise ValueError('existing native auditor drift')
    w=module('original_E171_control_auditor',OLD_EXEC/'worker.py')
    m=module('original_E170_cached_sources',OLD_EXEC/'matched.py')
    proof=w.audit_controls(m)
    if proof['sources']!=load_controls(bundle)['rows']:raise ValueError('reparsed eight native controls differ')
    for c in POLICY['channels']:faces.native_preflight(faces.prepare(c),c)
    return data,dict(status='native_atomic_release_preflight_no_DFT',native_frozen_SCFs_reparsed=8,
        original_runtime_files=380,native_controls=proof,centers=list(data),
        admission=scheduler_snapshot(),new_DFT=0,physical_reserved_geometry_generated_or_read=False)


class CachedPullback:
    def __init__(self,data,rows,delegate,point):
        self.data=data;self.delegate=delegate;self.hits=0;self.requests=0
        self.base=np.asarray(point['frozen_internal_A']);self.anchor=np.asarray(point['frozen_full_A'])
        self.records={hashlib.sha256(np.asarray(p['coordinates_full_A']).tobytes()).hexdigest():p for p in rows}
        self.last={};self.first=True
    def __call__(self,x):
        d=self.data;full=self.anchor+d['basis']@(x-self.base);self.requests+=1
        key=hashlib.sha256(full.tobytes()).hexdigest()
        row=self.records.get(key)
        if not np.any(full):
            e,g=d['center_energy'],d['center_gradient'];self.hits+=1
        elif row is not None and np.array_equal(full,row['coordinates_full_A']):
            e,g=row['frozen_energy_eV_cell'],np.asarray(row['physical_gradient_eV_A']);self.hits+=1
        else:
            if self.first:raise ValueError('first callback must reuse the exact frozen seed; no duplicate SCF')
            e,g=self.delegate(full)
        self.first=False;gi=d['basis'].T@g
        self.last=dict(coordinates_internal_A=x.tolist(),coordinates_full_A=full.tolist(),
            energy_eV_cell=float(e),physical_gradient_eV_A=np.asarray(g).tolist())
        return float(e),gi


def relax_point(data,q,r,evaluate,*,base=None):
    """Original safeguarded BFGS on pure atomic free coordinates; unchanged guards."""
    base=data['axes']@q if base is None else np.asarray(base);h=r.T@data['hessian']@r
    factor=np.linalg.cholesky((h+h.T)/2)
    inv=np.linalg.solve(factor.T,np.linalg.solve(factor,np.eye(31)))
    screens=dict(checks=0,rejected=0)
    def embed(z):return base+r@z
    def valid(z):
        x=embed(z);screens['checks']+=1
        ok=(np.array_equal(x[33:],base[33:]) and
            faces.canary.geometry_guard(data['chart'],data['basis'],x))
        screens['rejected']+=int(not ok);return bool(ok)
    def call(z):
        x=embed(z)
        if not np.allclose(data['axes'].T@x,q,atol=2e-12,rtol=0):raise ValueError('fixed Q changed')
        e,g=evaluate(x);return e,r.T@g
    z,calls=_safeguarded_bfgs(call,valid,np.zeros(31),gradient_tolerance=.03,
        max_iterations=200,initial_trust_radius=.05,initial_inverse_hessian=inv)
    e,g=call(z)
    if np.linalg.norm(g)>.03:raise RuntimeError('atomic released-gradient target not reached')
    return dict(energy_eV_cell=e,atomic_released_gradient_norm_eV_A=float(np.linalg.norm(g)),
        coordinates_internal_A=embed(z).tolist(),optimizer_callbacks=calls+1,geometry_screens=screens)


def compute(config,data,workdir,factory):
    workdir.mkdir();rows=[];controls=load_controls(Path(__file__).parent)['rows']
    for channel in POLICY['channels']:
        d=data[channel];r=np.asarray(config['channels'][channel]['atomic_release_basis_36_by_31'])
        for point in config['channels'][channel]['points']:
            directory=workdir/channel/f"node_{point['id']:03d}"
            delegate=faces.canary.CalculatorModeEvaluator(d['chart'],factory,workdir=directory,
                calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
                minimum_distance_A=1.6,result_validator=study.result_validator,
                validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',max_evaluation_attempts=12)
            evaluate=CachedPullback(d,[p for p in controls if p['channel']==channel],delegate,point)
            candidate=None;reason=None
            try:candidate=relax_point(d,np.asarray(point['q_A']),r,evaluate,base=point['frozen_internal_A'])
            except RuntimeError as exc:
                reason=str(exc)
                if (list(directory.glob('eval-*/failure.json')) or reason not in
                        ('mode evaluation attempt budget exhausted; calculator was not called',
                         'atomic released-gradient target not reached')):raise
            row=dict(channel=channel,id=point['id'],q_A=point['q_A'],candidate=candidate,reason=reason,
                status='atomic_gradient_candidate_stability_unchecked' if candidate else 'optimizer_not_converged',
                new_completed_SCFs=delegate.n_new_evaluations,reserved_attempts=delegate.n_evaluation_attempts,
                cached_callbacks=evaluate.hits,all_callbacks=evaluate.requests,
                candidate_native_identity=evaluate.last if candidate else None)
            study.exclusive_json(directory.parent/f"point_{point['id']:03d}.json",row);rows.append(row)
            print(json.dumps(dict(channel=channel,id=point['id'],status=row['status'],new_SCFs=row['new_completed_SCFs'])),flush=True)
    summary=dict(status='bounded_atomic_release_ablation_finished_not_stability_or_prediction',rows=rows,
        policy=POLICY,job_id=os.environ.get('SLURM_JOB_ID'),
        new_completed_SCFs=sum(p['new_completed_SCFs'] for p in rows),
        reserved_attempts=sum(p['reserved_attempts'] for p in rows),
        physical_reserved_geometry_generated_or_read=False)
    study.exclusive_json(workdir/'atomic_release.json',summary);return summary


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('prepare','preflight','run'))
    p.add_argument('--bundle',type=Path,default=Path(__file__).parent);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--config-sha256');p.add_argument('--output',type=Path);p.add_argument('--workdir',type=Path)
    a=p.parse_args()
    if a.mode=='prepare':
        if a.config.exists():raise FileExistsError('fresh atomic-release config required')
        study.exclusive_json(a.config,prepare(a.bundle));return
    if not a.config_sha256 or sha(a.config)!=a.config_sha256:raise ValueError('exact registered config required')
    config=json.loads(a.config.read_bytes())
    if a.mode=='preflight':
        _,proof=preflight(config,a.bundle)
        if a.output is None or a.output.exists():raise FileExistsError('fresh preflight output required')
        study.exclusive_json(a.output,proof);print(json.dumps(dict(status='preflight_passed_no_DFT')));return
    expected=dict(RUN_DFT='1',SLURM_JOB_PARTITION='hfacnormal01',SLURM_JOB_NUM_NODES='1',
        SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='32',SLURM_JOB_NAME=POLICY['job_name'],
        OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    if any(os.environ.get(k)!=v for k,v in expected.items()) or not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('authorized HF single-node32MPI allocation required')
    if (a.workdir is None or a.workdir.exists() or
            a.workdir.resolve()!=a.bundle.resolve().parent/'atomic_release' or
            not a.workdir.resolve().is_relative_to(study.RUNS)):
        raise ValueError('fresh dedicated atomic-release sibling required')
    if not scheduler_snapshot(os.environ['SLURM_JOB_ID'])['ready']:raise RuntimeError('two-allocation admission not ready')
    data,proof=preflight(config,a.bundle)
    study.exclusive_json(a.bundle.parent/'run_preflight.json',proof)
    compute(config,data,a.workdir,study.make_factory(parameters={'source_directory':str(study.ASSETS)},command=study.COMMAND))


if __name__=='__main__':main()
