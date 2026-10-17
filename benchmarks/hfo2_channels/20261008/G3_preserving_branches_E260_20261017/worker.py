"""Ten independently seeded, fixed-Q training candidates; no test labels."""
from pathlib import Path
from copy import deepcopy
import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import sys
sys.dont_write_bytecode = True
LOCAL = Path(__file__).resolve().parents[4]
if (LOCAL/'vcneb').is_dir():
    sys.path.insert(0, str(LOCAL))
import numpy as np


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


HERE = Path(__file__).resolve().parent
legacy = module('E260_legacy_soft', HERE/'legacy_soft_worker.py')
solver = module('E260_branch_solver', HERE/'branch_relaxation.py')
faces, study = legacy.faces, legacy.study
EXEC = study.RUNS/'20261017-varneb-G3-P-branches-E260-r2'
OLD_SOFT = study.RUNS/'20261016-varneb-G3-soft4-E256-r1/bundle'
PROBE_SHA = '761bb4bd2784fca366cfbb32abc4800647b8545e8281b5ccd03cf93028ceb318'
LEGACY_SHA = '10fa88e50d7f5b5caa983ff150dd1855a556b18c0159729e0022e75ffec17779'
POLICY = dict(q_A=[[0.,0.],[-.1,0.],[0.,-.1],[0.,.1],[.1,0.]],
    zeta_start_A=[.1,-.1], shards=[[0,2,4,6,8],[1,3,5,7,9]],
    maximum_attempts_per_candidate=12, maximum_attempts_cohort=120,
    training_before=dict(T=34,preserving=39),training_committed=dict(T=34,preserving=49),
    training_cap_each=81,independent_inventory=16,T_independent_seen=8,preserving_independent_seen=0,
    released_gradient_tolerance_eV_A=.03,zeta_gradient_tolerance_eV_A=1e-4,
    initial_trust_radius_A=.05,maximum_iterations=200,maximum_internal_radius_A=.5,
    minimum_distance_A=1.6,partition='hfacnormal01',nodes=1,MPI_ranks=32,
    allocation_hours=4,maximum_committed_allocations=2,automatic_retries=0,
    climb=False,physical_plus005_generated_or_read=False)
sha = legacy.sha


def inputs(bundle):
    if study.text_sha(bundle/'legacy_soft_worker.py') != LEGACY_SHA:
        raise ValueError('original soft worker source differs')
    data, x0, direction, soft = legacy.inputs(bundle)
    probes = legacy.pinned(bundle/'soft_probe_analysis.json', PROBE_SHA)
    if (probes['status']!='four_training_statics_signed_energy_force_diagnostic_not_independent_test'
            or [r['t_A'] for r in probes['rows']] != [.05,-.05,.10,-.10]):
        raise ValueError('all four already seen training probes required')
    if (LOCAL/'vcneb').is_dir():
        base = LOCAL/'benchmarks/hfo2_channels/20261008'
        cp = base/'G3_released_curvature_E181_20261015/config.json'
        ap = base/'G3_released_complete_E254_20261016/native_audit.json'
    else:
        cp = study.RUNS/'20261015-varneb-G3-release-curvature-E181-r1/bundle/config.json'
        ap = study.RUNS/'20261015-varneb-G3-released-audit-E190-r1/complete_20261016_2331.json'
    old = legacy.pinned(cp, legacy.CONFIG_SHA)
    audit = legacy.pinned(ap, legacy.AUDIT_SHA)
    free = np.asarray(old['channels']['preserving']['released_basis_36_by_34'])
    row = next(r for r in audit['channels'] if r['channel']=='preserving')
    raw = np.asarray(row['blocks'][0]['raw_projected_eV_A2'])
    d = data['basis'].T @ direction
    axes = data['axes']
    transverse = np.linalg.svd(np.column_stack([axes,d]).T, full_matrices=True)[2][3:].T
    measured = free @ (.5*(raw+raw.T)) @ free.T
    stable = transverse.T @ measured @ transverse
    if np.linalg.eigvalsh(stable)[0] <= 0:
        raise ValueError('positive measured R33 restriction required, no clipping')
    pairs = []
    for h in (.05,.10):
        plus=next(r for r in probes['rows'] if r['t_A']==h)
        minus=next(r for r in probes['rows'] if r['t_A']==-h)
        pairs.append((plus['gradient_along_direction_eV_A']-minus['gradient_along_direction_eV_A'])/(2*h))
    cubic=(pairs[1]-pairs[0])/(.10**2-.05**2)
    linear=pairs[0]-cubic*.05**2
    stiffness=linear+3*cubic*.10**2
    if not np.isfinite(stiffness) or stiffness<=0:
        raise ValueError('outer-probe numerical stiffness must be positive')
    preconditioner=transverse@stable@transverse.T+stiffness*np.outer(d,d)+axes@axes.T
    preconditioner=.5*(preconditioner+preconditioner.T)
    return data,x0,direction,free,preconditioner,dict(
        odd_linear_eV_A2=float(linear),odd_cubic_eV_A4=float(cubic),
        outer_numerical_stiffness_eV_A2=float(stiffness),
        raw_R33_min_eV_A2=float(np.linalg.eigvalsh(stable)[0]),
        not_a_physical_Hessian_or_stability_certificate=True)


def geometry(data, x, q):
    b,a=data['basis'],data['axes']
    if (x.shape!=(40,) or not np.isfinite(x).all() or x[-1]!=0.
            or not np.allclose(b@(b.T@x),x,atol=2e-12,rtol=0)
            or not np.allclose(a.T@(b.T@x),q,atol=2e-12,rtol=0)
            or not faces.canary.geometry_guard(data['chart'],b,b.T@x)):
        raise ValueError('fixed Q/control/span/radius/minimum-distance contract violated')
    atoms=data['chart'].to_atoms(x)
    if (np.linalg.det(atoms.cell.array)<=0 or np.linalg.cond(atoms.cell.array)>1e12
            or not np.array_equal(atoms.cell.array[:2],data['chart'].reference_atoms.cell.array[:2])):
        raise ValueError('positive finite unchanged substrate required')
    return atoms


def prepare(bundle):
    data,x0,direction,free,preconditioner,diagnostic=inputs(bundle)
    points=[]
    for qi,q in enumerate(POLICY['q_A']):
        for sign,zeta in enumerate(POLICY['zeta_start_A']):
            start=x0+data['basis']@data['axes']@np.asarray(q)+zeta*direction
            geometry(data,start,q)
            points.append(dict(id=len(points),Q_id=qi,q_A=q,start_sign=1 if zeta>0 else -1,
                zeta_start_A=zeta,starting_full_coordinates_A=start.tolist(),
                starting_internal_coordinates_A=(data['basis'].T@start).tolist()))
    return dict(schema_version=1,purpose='author_approved_ten_training_branch_candidates_not_holdout',
        policy=deepcopy(POLICY),points=points,preconditioner_36_by_36=preconditioner.tolist(),
        free_basis_36_by_34=free.tolist(),branch_direction_36=(data['basis'].T@direction).tolist(),
        preconditioner_diagnostic=diagnostic,protocol_sha256=sha(bundle/'protocol.md'),
        source_LF_sha256={n:study.text_sha(bundle/n) for n in
            ('worker.py','branch_relaxation.py','legacy_soft_worker.py','run.slurm')},
        soft_direction_sha256=legacy.SOFT_SHA,soft_probe_analysis_sha256=PROBE_SHA,
        complete_native_audit_sha256=legacy.AUDIT_SHA,physical_input_sha256=study.CONTRACT,
        binary_sha256=study.BINARY_SHA,original_2D_refusal_preserved=True,
        independent_or_stationary_or_global_stability_certified=False)


def validate(config,bundle):
    expected=prepare(bundle)
    numerical=('preconditioner_36_by_36','preconditioner_diagnostic','branch_direction_36','points')
    left={k:v for k,v in config.items() if k not in numerical}
    right={k:v for k,v in expected.items() if k not in numerical}
    if (set(config)!=set(expected) or
            json.dumps(left,sort_keys=True,allow_nan=False)!=json.dumps(right,sort_keys=True,allow_nan=False)):
        raise ValueError('frozen training sources/coordinates/policy differ')
    if (not np.allclose(config[numerical[0]],expected[numerical[0]],atol=2e-12,rtol=0)
            or set(config[numerical[1]])!=set(expected[numerical[1]])
            or any(not np.isclose(config[numerical[1]][k],v,atol=2e-12,rtol=0)
                   for k,v in expected[numerical[1]].items())):
        raise ValueError('frozen numerical preconditioner differs from source reconstruction')
    if not np.allclose(config['branch_direction_36'],expected['branch_direction_36'],atol=2e-12,rtol=0):
        raise ValueError('frozen branch direction differs from arithmetic reconstruction')
    arrays=('starting_full_coordinates_A','starting_internal_coordinates_A')
    if len(config['points'])!=10:raise ValueError('frozen training points differ')
    for p,e in zip(config['points'],expected['points']):
        if ({k:v for k,v in p.items() if k not in arrays}!={k:v for k,v in e.items() if k not in arrays}
                or any(not np.allclose(p[k],e[k],atol=2e-12,rtol=0) for k in arrays)):
            raise ValueError('frozen training point definition differs')
    values=list(inputs(bundle))
    values[4]=np.asarray(config['preconditioner_36_by_36'])
    return tuple(values)  # use frozen bytes, never a platform-regenerated optimizer seed


def queue(own=None,shard=None):
    proof=legacy.queue()
    committed=proof['committed']
    proof['ready']=(len(committed)<2 if own is None else len(committed)<=2 and any(
        r['job_id']==own and r['name']==f'hfo2-G3-E260-P{shard}' for r in committed))
    return proof


def native_context(config,bundle):
    local=validate(config,bundle)
    original=legacy.pinned(OLD_SOFT/'config.json','f7dc1053bff96bbd141d0a7e2e5c14d949f91f69832852d61bf7af9b932de7a7')
    reader,w,data,_,_,_=legacy.native_context(original,OLD_SOFT)
    if not np.array_equal(local[0]['origin_record']['coordinates'],data['origin_record']['coordinates']):
        raise ValueError('exact unchanged native origin required')
    return reader,w,*local


def preflight(config,bundle):
    if bundle.resolve()!=EXEC/'bundle' or (EXEC/'inert_decks').exists():
        raise ValueError('fresh exact inert namespace required')
    _,_,data,*_=native_context(config,bundle)
    decks=[]
    for p in config['points']:
        atoms=geometry(data,data['basis']@np.asarray(p['starting_internal_coordinates_A']),p['q_A'])
        d=EXEC/'inert_decks'/f"point_{p['id']:03d}";d.mkdir(parents=True)
        for name in study.CONTRACT:
            if name in ('INPUT','KPT'):shutil.copyfile(study.ASSETS/name,d/name)
            else:(d/name).symlink_to(study.ASSETS/name)
        legacy.write_structure(d/'STRU',atoms);legacy.geometry_roundtrip(d/'STRU',atoms)
        hashes={n:sha(d/n) for n in (*study.CONTRACT,'STRU')}
        if any(hashes[n]!=h for n,h in study.CONTRACT.items()):raise ValueError('inert input drift')
        decks.append(dict(id=p['id'],input_sha256=hashes))
    return dict(status='ten_native_inert_branch_starts_passed_no_DFT',decks=decks,
        config_sha256=sha(bundle/'config.json'),admission=queue(),new_DFT=0)


def run(config,bundle,shard):
    job=os.environ.get('SLURM_JOB_ID')
    required=dict(RUN_DFT='1',SLURM_JOB_PARTITION='hfacnormal01',SLURM_JOB_NUM_NODES='1',
        SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='32',SLURM_JOB_NAME=f'hfo2-G3-E260-P{shard}',
        OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    if shard not in (0,1) or not job or any(os.environ.get(k)!=v for k,v in required.items()):
        raise ValueError('authorized one-node32MPI allocation required')
    out=EXEC/f'shard_{shard}'
    if bundle.resolve()!=EXEC/'bundle' or out.exists():raise ValueError('fresh exact output required')
    accepted=json.loads((EXEC/f'accepted_{shard}.json').read_bytes())
    if accepted['job_id']!=job or accepted['config_sha256']!=sha(bundle/'config.json') or not queue(job,shard)['ready']:
        raise ValueError('accepted handle/config/allocation gate failed')
    reader,w,data,_,_,free,seed,_=native_context(config,bundle)
    out.mkdir();rows=[]
    for p in config['points']:
        if p['id'] not in POLICY['shards'][shard]:continue
        d=out/f"point_{p['id']:03d}"
        evaluator=faces.canary.CalculatorModeEvaluator(data['chart'],
            study.make_factory(parameters={'source_directory':str(study.ASSETS)},command=study.COMMAND),
            workdir=d,calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
            minimum_distance_A=1.6,result_validator=study.result_validator,
            validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',max_evaluation_attempts=12)
        evaluate=faces.canary.InternalEvaluator(evaluator,data['basis'],data['center_energy'],data['center_gradient'])
        def guard(x):
            try:geometry(data,data['basis']@x,p['q_A']);return True
            except ValueError:return False
        row=dict(id=p['id'],q_A=p['q_A'],start_sign=p['start_sign'],result=None)
        try:
            row['result']=solver.relax_branch_candidate(evaluate,np.asarray(p['starting_internal_coordinates_A']),
                free,np.asarray(config['branch_direction_36']),initial_hessian=seed,valid=guard,
                gradient_tolerance=.03,branch_gradient_tolerance=1e-4,initial_trust_radius=.05,max_iterations=200)
            row['status']=row['result']['status']
        except RuntimeError as exc:
            if str(exc)!='mode evaluation attempt budget exhausted; calculator was not called' or list(d.glob('eval-*/failure.json')):
                raise
            row.update(status='attempt_budget_exhausted_no_retry',reason=str(exc))
        native=[];sources=[]
        for path in sorted(d.glob('eval-*/result.json')):
            x=np.asarray(json.loads(path.read_bytes())['coordinates']);geometry(data,x,p['q_A'])
            record,source=reader.native_record(w,data,x,path,d);native.append(record);sources.append(source)
        if len(native)!=evaluator.n_new_evaluations or evaluator.n_evaluation_attempts>12:
            raise ValueError('charged native call accounting differs')
        row.update(reserved_attempts=evaluator.n_evaluation_attempts,new_native_SCFs=len(native),native_sources=sources)
        if row['result'] is not None:
            final=data['basis']@np.asarray(row['result']['coordinates'])
            matches=[r for r in native if np.array_equal(r['coordinates'],final)]
            if len(matches)!=1:raise ValueError('unique complete native final candidate required')
            record=matches[0];internal_g=data['basis'].T@np.asarray(record['gradient'])
            norm=float(np.linalg.norm(free.T@internal_g));dg=float(np.asarray(config['branch_direction_36'])@internal_g)
            if not np.isclose(norm,row['result']['released_gradient_norm'],atol=2e-12,rtol=0) or not np.isclose(dg,row['result']['branch_gradient'],atol=2e-12,rtol=0):
                raise ValueError('native final physical gradient differs')
            row.update(native_final=record,zeta_final_A=float(np.asarray(config['branch_direction_36'])@
                (np.asarray(row['result']['coordinates'])-data['basis'].T@np.asarray(data['origin_record']['coordinates']))))
        study.exclusive_json(out/f"point_{p['id']:03d}.json",row);rows.append(row)
        print(json.dumps(dict(id=p['id'],status=row['status'],SCFs=len(native),zeta=row.get('zeta_final_A'))),flush=True)
    if [r['id'] for r in rows]!=POLICY['shards'][shard] or sum(r['reserved_attempts'] for r in rows)>60:
        raise ValueError('five original starts and60-attempt shard cap required')
    study.exclusive_json(out/'completed.json',dict(status='five_training_branch_starts_processed_stability_unchecked',
        job_id=job,shard=shard,config_sha256=sha(bundle/'config.json'),rows=rows,
        new_independent_labels=0,physical_plus005_generated_or_read=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('prepare','preflight','run'))
    p.add_argument('--bundle',type=Path,required=True);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--config-sha256');p.add_argument('--shard',type=int,choices=(0,1));p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.mode=='prepare':study.exclusive_json(a.config,prepare(a.bundle))
    else:
        c=legacy.pinned(a.config,a.config_sha256)
        if a.mode=='preflight':study.exclusive_json(a.output,preflight(c,a.bundle))
        else:run(c,a.bundle,a.shard)
