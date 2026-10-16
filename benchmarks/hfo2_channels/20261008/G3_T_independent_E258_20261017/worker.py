"""Execute only the eight unchanged E185 T tests; preserving remains refused."""
from copy import deepcopy
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[4]
if (ROOT/'vcneb').is_dir(): sys.path.insert(0,str(ROOT))
import numpy as np
from scripts import hfo2_G3_conditional_faces as faces
from examples.hfo2_fixed_input_factory import write_structure, geometry_roundtrip

study=faces.study
EXEC=study.RUNS/'20261017-varneb-G3-T8-E258-r1'
RUNTIME=study.RUNS/'20261015-varneb-G3-faces-E161-r1/source'
READER=study.RUNS/'20261015-varneb-G3-released-audit-E190-r1/audit.py'
READER_SHA='51d0452968394ff03d26efface7705ee31c7679215cb2204f6c69e13a5391d4f'
PINS={
 'standard_forecast.json':'aa77cdd7460ad809a5fbe06e80c684f7f73ce118835067c380ec219b96abd2d8',
 'reference_forecast.json':'97692cb0989a4bafdff86349a732126cce5d7d93d9a7d47926ce844fb1f0b285',
 'complete_native.json':'a389154749a07109a3a91499343b86f06a2c185fd310611570637d71193e3381',
 'T_curvature.json':'66f7cdc9648afbfc3b61fc63a3835ec5868350c61aa921e56feea578a8498767',
 'B3_native.json':'8a2fdf2a2beec306bb2f39aed8b347019a4f049b2d8cbc42cd79de1fce922b92'}
POLICY=dict(channel='T',target_ids=list(range(1,9)),shards=[[1,2,3,4],[5,6,7,8]],
 maximum_new_attempts_per_point=12,maximum_new_attempts_this_cohort=96,
 original_maximum_attempts_two_faces=192,optimizer='safeguarded_bfgs',
 gradient_tolerance_eV_A=.03,initial_trust_radius_A=.05,maximum_iterations=200,
 maximum_internal_displacement_A=.5,minimum_distance_A=1.6,
 original_independent_reserved_per_face=8,physical_plus005_generated_or_read=False,
 preserving_admission_refused=True,training_committed=dict(T=34,preserving=39),
 automatic_retries=0,partition='hfacnormal01',nodes=1,MPI_ranks=32,
 allocation_hours_per_shard=4,maximum_committed_study_allocations=2,climb=False)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def pinned(path,digest):
    if sha(path)!=digest: raise ValueError('frozen source bytes differ')
    return json.loads(Path(path).read_bytes())


def prepare(bundle):
    f,ref,full,curv,b3=[pinned(bundle/n,h) for n,h in PINS.items()]
    policy=f['test_execution_policy']
    if (policy['maximum_new_attempts_per_point']!=12 or policy['maximum_new_attempts_two_faces']!=192
        or policy['optimizer']!='safeguarded_bfgs' or policy['gradient_tolerance_eV_A']!=.03
        or policy['initial_trust_radius_A']!=.05 or policy['maximum_iterations']!=200
        or policy['maximum_internal_displacement_A']!=.5 or policy['minimum_distance_A']!=1.6
        or f['local_test_labels_seen'] or not b3['all10points_covered'] or b3['pending'] is not None
        or b3['native_new_SCFs_reparsed']!=36 or not full['accounting']['normal_completion']
        or full['native_new_stencil_SCFs_reparsed']!=272 or curv['channel']!='T'
        or curv['new_native_SCFs']!=136):
        raise ValueError('original prediction policy and complete pre-label audits required')
    row=next(r for r in f['channels'] if r['channel']=='T')
    data=faces.load_inputs('T'); targets=deepcopy(row['targets'])
    if len(targets)!=8 or [p['id'] for p in targets]!=POLICY['target_ids']:
        raise ValueError('all original eight targets required; no subset selection')
    minima=[]
    for block in curv['blocks']:
        raw=np.asarray(block['raw_projected_eV_A2'])
        if raw.shape!=(34,34) or block['complete_pairs']!=34:
            raise ValueError('complete raw transverse block required')
        w=np.linalg.eigvalsh(.5*(raw+raw.T))
        saved=np.asarray(block['symmetric_projected_eigenvalues_eV_A2'])
        if not np.allclose(w,saved,atol=2e-12,rtol=0):raise ValueError('stored raw T spectrum differs')
        minima.append(float(saved[0]))  # use original serialized values, not platform-dependent eigensolver bytes
        if np.any(w<=0) or np.any(saved<=0):raise ValueError('T necessary raw sign admission refused')
    if [b['step_A'] for b in curv['blocks']]!=[.01,.02]:
        raise ValueError('both original step amplitudes required')
    refs=[r for r in ref['rows'] if r['channel']=='T']
    if len(refs)!=14:raise ValueError('all original references and both steps required')
    for point in targets:
        x=np.asarray(point['starting_internal_coordinates_A'])
        if (point['independent_DFT_label'] is not None or point['training_parent_id']!=0
            or not np.allclose(data['basis']@x,point['starting_full_coordinates_A'],atol=2e-12,rtol=0)
            or not np.allclose(data['axes'].T@x,point['q_A'],atol=2e-12,rtol=0)
            or not faces.canary.geometry_guard(data['chart'],data['basis'],x)):
            raise ValueError('original single-start branch/target/guard required')
    return dict(schema_version=1,purpose='original_T_eight_independent_tests_not_replacement_or_physical_strain_holdout',
        policy=deepcopy(POLICY),points=targets,standard_predictions=deepcopy(row),
        reference_predictions=deepcopy(refs),original_scoring=deepcopy(f['scoring']),
        original_source_sha256=deepcopy(PINS),raw_T_minima_eV_A2=minima,
        final_models_for_this_original_local_cohort='all five E185 arms and all fourteen E188 T reference/step rows unchanged',
        full_B0_B5_physical_strain_freeze_or_JCTC_readiness=False)


def validate(config,bundle):
    if config!=prepare(bundle):raise ValueError('frozen original targets/starts/models/policy differ')
    return faces.load_inputs('T')


def queue(own=None,shard=None):
    raw=subprocess.run(['squeue','-h','-u','iai806','-p','hfacnormal01','-o','%i|%100j|%T|%C|%D'],
        check=True,capture_output=True,text=True).stdout
    rows=[]
    for line in raw.splitlines():
        f=[s.strip() for s in line.split('|')]
        if len(f)!=5:raise ValueError('unparsed study queue')
        if f[1].startswith('hfo2-G3-'):
            if f[3:]!=['32','1']:raise ValueError('study resource layout differs')
            rows.append(dict(job_id=f[0],name=f[1],state=f[2]))
    ready=len(rows)<2 if own is None else len(rows)<=2 and any(
        r['job_id']==own and r['name']==f'hfo2-G3-E258-T{shard}' for r in rows)
    return dict(ready=ready,committed=rows,scheduler_writes=0)


def native_context(config,bundle):
    validate(config,bundle)
    if study.ROOT.resolve()!=RUNTIME or sha(READER)!=READER_SHA:
        raise ValueError('exact original runtime/reader required')
    spec=importlib.util.spec_from_file_location('E258_original_reader',READER)
    reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
    w,_,original=reader.load_worker()  # checks original380 files/binary/six physical assets
    data=original['T'];old=data['origin_record'];raw=Path(old['validation']['raw_directory'])
    checked,_=reader.native_record(w,data,np.asarray(old['coordinates']),raw.parent/'result.json',raw.parent.parent)
    if checked!=old:raise ValueError('original native cached origin differs')
    return reader,w,data


def preflight(config,bundle):
    if bundle.resolve()!=EXEC/'bundle' or (EXEC/'inert_decks').exists():
        raise ValueError('fresh exact inert namespace required')
    _,_,data=native_context(config,bundle);decks=[]
    for point in config['points']:
        atoms=data['chart'].to_atoms(np.asarray(point['starting_full_coordinates_A']))
        d=EXEC/'inert_decks'/f"point_{point['id']:03d}";d.mkdir(parents=True)
        for n in study.CONTRACT:
            if n in ('INPUT','KPT'):shutil.copyfile(study.ASSETS/n,d/n)
            else:(d/n).symlink_to(study.ASSETS/n)
        write_structure(d/'STRU',atoms);geometry_roundtrip(d/'STRU',atoms)
        hashes={n:sha(d/n) for n in (*study.CONTRACT,'STRU')}
        if any(hashes[n]!=h for n,h in study.CONTRACT.items()):raise ValueError('inert physical inputs drift')
        decks.append(dict(id=point['id'],input_sha256=hashes))
    return dict(status='eight_original_T_inert_decks_passed_no_DFT',decks=decks,admission=queue(),new_DFT=0)


def optimize_point(data,point,evaluate,guard):
    plane=faces.canary.ModePlane(np.zeros(36),data['axes'],np.ones(36),np.eye(2),
        ('original coupled internal1','original coupled internal2'),'local length','E185-original-test')
    return faces.canary.relax_orthogonal_at_q(plane,point['q_A'],evaluate,
        starts=[np.asarray(point['starting_internal_coordinates_A'])],include_frozen_start=False,
        optimizer='safeguarded_bfgs',trial_validator=guard,initial_hessian=data['hessian'],
        initial_trust_radius=.05,gradient_tolerance=.03,max_iterations=200)


def run(config,bundle,shard):
    job=os.environ.get('SLURM_JOB_ID'); name=f'hfo2-G3-E258-T{shard}'
    required=dict(RUN_DFT='1',SLURM_JOB_PARTITION='hfacnormal01',SLURM_JOB_NUM_NODES='1',
        SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='32',SLURM_JOB_NAME=name,
        OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    if shard not in (0,1) or not job or any(os.environ.get(k)!=v for k,v in required.items()):
        raise ValueError('authorized one-node32MPI allocation required')
    out=EXEC/f'shard_{shard}'
    if bundle.resolve()!=EXEC/'bundle' or out.exists():raise ValueError('fresh exact independent output required')
    accepted=json.loads((EXEC/f'accepted_{shard}.json').read_bytes())
    if accepted['job_id']!=job or accepted['config_sha256']!=sha(bundle/'config.json') or not queue(job,shard)['ready']:
        raise ValueError('accepted handle/config/committed-allocation gate failed')
    reader,w,data=native_context(config,bundle);out.mkdir();rows=[]
    for point in config['points']:
        if point['id'] not in POLICY['shards'][shard]:continue
        d=out/f"point_{point['id']:03d}"
        evaluator=faces.canary.CalculatorModeEvaluator(data['chart'],
            study.make_factory(parameters={'source_directory':str(study.ASSETS)},command=study.COMMAND),
            workdir=d,calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
            minimum_distance_A=1.6,result_validator=study.result_validator,
            validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',max_evaluation_attempts=12)
        evaluate=faces.canary.InternalEvaluator(evaluator,data['basis'],data['center_energy'],data['center_gradient'])
        guard=lambda x:faces.canary.geometry_guard(data['chart'],data['basis'],x)
        row=dict(id=point['id'],q_A=point['q_A'],original_start_full_A=point['starting_full_coordinates_A'],point=None)
        try:
            result=optimize_point(data,point,evaluate,guard)
            row.update(status='single_start_gradient_candidate_stability_unchecked',point=dict(
                coordinates_internal_A=result.coordinates.tolist(),coordinates_full_A=(data['basis']@result.coordinates).tolist(),
                energy_eV_cell=result.energy,orthogonal_gradient_norm_eV_A=result.orthogonal_gradient_norm))
        except RuntimeError as exc:
            if str(exc) not in ('mode evaluation attempt budget exhausted; calculator was not called',
                               'no starting point reached the orthogonal-gradient tolerance') or list(d.glob('eval-*/failure.json')):
                raise
            row.update(status='optimizer_not_converged_no_retry',reason=str(exc))
        paths=sorted(d.glob('eval-*/result.json'));native=[];sources=[]
        if len(paths)!=evaluator.n_new_evaluations or evaluator.n_evaluation_attempts>12:
            raise ValueError('original point accounting differs')
        for path in paths:
            x=np.asarray(json.loads(path.read_bytes())['coordinates'])
            record,source=reader.native_record(w,data,x,path,d);native.append(record);sources.append(source)
        row.update(reserved_attempts=evaluator.n_evaluation_attempts,new_native_SCFs=len(native),native_sources=sources)
        if row['point'] is not None:
            matches=[r for r in native if np.array_equal(r['coordinates'],row['point']['coordinates_full_A'])]
            if len(matches)!=1:raise ValueError('exact complete native candidate required')
            row['native_candidate']=matches[0]
        study.exclusive_json(out/f"point_{point['id']:03d}.json",row);rows.append(row)
        print(json.dumps(dict(id=point['id'],status=row['status'],SCFs=len(native))),flush=True)
    if [r['id'] for r in rows]!=POLICY['shards'][shard] or sum(r['reserved_attempts'] for r in rows)>48:
        raise ValueError('original shard budget differs')
    study.exclusive_json(out/'completed.json',dict(status='four_original_T_test_targets_processed_not_global_stability',
        job_id=job,shard=shard,config_sha256=sha(bundle/'config.json'),rows=rows,
        preserving_test_points_computed=0,physical_plus005_generated_or_read=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('prepare','preflight','run'))
    p.add_argument('--bundle',type=Path,required=True);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--config-sha256');p.add_argument('--shard',type=int,choices=(0,1));p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.mode=='prepare':study.exclusive_json(a.config,prepare(a.bundle))
    else:
        config=pinned(a.config,a.config_sha256)
        if a.mode=='preflight':study.exclusive_json(a.output,preflight(config,a.bundle))
        else:run(config,a.bundle,a.shard)
