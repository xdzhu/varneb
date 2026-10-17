"""Execute exactly the eight frozen P candidate tests, never physical strain."""
from pathlib import Path
from copy import deepcopy
import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
sys.dont_write_bytecode=True
import numpy as np

ROOT=Path(__file__).resolve().parents[4]
EXEC=Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-G3-P8-E266-r1')
OLD=(ROOT/'benchmarks/hfo2_channels/20261008/G3_preserving_branches_E260_20261017'
     if (ROOT/'vcneb').is_dir() else EXEC.parent/'20261017-varneb-G3-P-branches-E260-r2/bundle')
PINS=dict(forecasts='d1d368372537a3d57d6a0cf38716ac8ea67210dbcac463d0dd2296a13e7b2046',
          starts='90e63ea901c15abc19cb1e552bde55f85221fc928032947be6b30b401bb8b2fc',
          protocol='269c9385eec045f1de980ee013eb04e94a34556007085cb02a7b19d429b6ae82')
OLD_CONFIG='63b1aa0eea88e7cd7c2f1a78ce9f4f3c7f39e16a9caf04d3f6d8c6ef5f84ab80'
OLD_MANIFEST='4422902d6474afe1056181b978f6098c9fc1284ff6cc66d11ff865d499d87eff'
POLICY=dict(shards=[[1,2,3,4],[5,6,7,8]],maximum_new_attempts_per_point=12,
    maximum_new_attempts_cohort=96,gradient_tolerance_eV_A=.03,zeta_gradient_tolerance_eV_A=1e-4,
    initial_trust_radius_A=.05,maximum_iterations=200,maximum_internal_radius_A=.5,minimum_distance_A=1.6,
    partition='hfacnormal01',nodes=1,MPI_ranks=32,allocation_hours=4,maximum_committed_allocations=2,
    automatic_retries=0,training_committed=dict(T=34,preserving=49),T8_already_seen=True,
    physical_plus005_generated_or_read=False,climb=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);value=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value);return value


def old_worker():
    if sha(OLD/'config.json')!=OLD_CONFIG or sha(OLD/'bundle_files.json')!=OLD_MANIFEST:
        raise ValueError('exact original training configuration and runtime bundle required')
    files=json.loads((OLD/'bundle_files.json').read_bytes())
    if len(files)!=8 or any(Path(n).name!=n or sha(OLD/n)!=h for n,h in files.items()):
        raise ValueError('original immutable eight-file bundle differs')
    return module('E266_original_E260_worker',OLD/'worker.py'),json.loads((OLD/'config.json').read_bytes())


def prepare(bundle):
    if any(sha(bundle/(name+'.json' if name!='protocol' else 'protocol.md'))!=digest for name,digest in PINS.items()):
        raise ValueError('exact published forecasts, training-only starts and protocol required')
    forecast=json.loads((bundle/'forecasts.json').read_bytes());starts=json.loads((bundle/'starts.json').read_bytes())
    old,original=old_worker();data,*_=old.validate(original,OLD)
    if (forecast['independent_P8_labels_seen'] or forecast['physical_plus005_or_private_arrays_accessed']
            or forecast['source_sha256']['P8_protocol']!=PINS['protocol'] or not forecast['original_conditional_2D_refusal_preserved']
            or len(forecast['new_models'])!=4 or len(forecast['original_reference_preserving_records'])!=14
            or starts['complete_native_training_audit_sha256']!='30971f61ab7d12d8e9bfa729a517603566c20cefa6bf442b06a4100264bd8811'
            or forecast['physical_input_sha256']!=original['physical_input_sha256']):
        raise ValueError('complete pre-label panel and same original physics required')
    points=deepcopy(starts['points'])
    if [p['id'] for p in points]!=list(range(1,9)) or [p['training_parent_id'] for p in points]!=[2,4,7,8,2,2,4,7]:
        raise ValueError('all original eight targets and training-only parents required')
    for point in points:
        x=np.asarray(point['starting_internal_coordinates_A']);full=data['basis']@x
        if point['independent_DFT_label'] is not None or not np.allclose(full,point['starting_full_coordinates_A'],atol=2e-12,rtol=0):
            raise ValueError('frozen internal starts and existing arithmetic reconstruction rule required')
        old.geometry(data,full,point['q_A'])
    return dict(schema_version=1,purpose='remaining_P8_frozen_candidate_target_not_stable_sheet_or_physical_strain',
        policy=deepcopy(POLICY),points=points,source_input_sha256=deepcopy(PINS),
        original_config_sha256=OLD_CONFIG,original_bundle_manifest_sha256=OLD_MANIFEST,
        physical_input_sha256=deepcopy(original['physical_input_sha256']),binary_sha256=original['binary_sha256'],
        source_LF_sha256={n:hashlib.sha256((bundle/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                         for n in ('worker.py','run.slurm')},
        original_2D_refusal_preserved=True,physical_barrier_panel_frozen=False,stable_conditional_model_admitted=False)


def validate(config,bundle,native=False):
    if config!=prepare(bundle):raise ValueError('immutable local test configuration differs')
    old,original=old_worker()
    if native:return old,old.native_context(original,OLD)
    return old,old.validate(original,OLD)


def queue(own=None,shard=None):
    raw=subprocess.check_output(['squeue','-h','-u','iai806','-p','hfacnormal01','-o','%i|%100j|%T|%C|%D'],text=True)
    rows=[]
    for line in raw.splitlines():
        f=[s.strip() for s in line.split('|')]
        if len(f)!=5:raise ValueError('unparsed scheduler allocation')
        if f[1].startswith('hfo2-'):
            if f[3:]!=['32','1']:raise ValueError('study resource layout differs')
            rows.append(dict(job_id=f[0],name=f[1],state=f[2]))
    ready=(len(rows)<2 if own is None else len(rows)<=2 and any(
        r['job_id']==own and r['name']==f'hfo2-G3-E266-P{shard}' for r in rows))
    return dict(ready=ready,committed=rows,scheduler_writes=0)


def preflight(config,bundle):
    if bundle.resolve()!=EXEC/'bundle' or (EXEC/'inert_decks').exists():raise ValueError('fresh exact inert namespace required')
    old,context=validate(config,bundle,native=True);_,_,data,*_=context;decks=[]
    for point in config['points']:
        atoms=old.geometry(data,data['basis']@np.asarray(point['starting_internal_coordinates_A']),point['q_A'])
        d=EXEC/'inert_decks'/f"point_{point['id']:03d}";d.mkdir(parents=True)
        for name in old.study.CONTRACT:
            if name in ('INPUT','KPT'):shutil.copyfile(old.study.ASSETS/name,d/name)
            else:(d/name).symlink_to(old.study.ASSETS/name)
        old.legacy.write_structure(d/'STRU',atoms);old.legacy.geometry_roundtrip(d/'STRU',atoms)
        hashes={n:sha(d/n) for n in (*old.study.CONTRACT,'STRU')}
        if any(hashes[n]!=h for n,h in old.study.CONTRACT.items()):raise ValueError('physical input drift')
        decks.append(dict(id=point['id'],input_sha256=hashes))
    return dict(status='eight_native_P8_inert_starts_passed_no_DFT',decks=decks,new_DFT=0,
        config_sha256=sha(bundle/'config.json'),admission=queue(),original_runtime_and_assets_rechecked=True)


def run(config,bundle,shard):
    job=os.environ.get('SLURM_JOB_ID')
    required=dict(RUN_DFT='1',SLURM_JOB_PARTITION='hfacnormal01',SLURM_JOB_NUM_NODES='1',
        SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='32',SLURM_JOB_NAME=f'hfo2-G3-E266-P{shard}',
        OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    if shard not in (0,1) or not job or any(os.environ.get(k)!=v for k,v in required.items()):
        raise ValueError('authorized one-node32MPI allocation required')
    out=EXEC/f'shard_{shard}'
    if bundle.resolve()!=EXEC/'bundle' or out.exists():raise ValueError('fresh exact independent output required')
    accepted=json.loads((EXEC/f'accepted_{shard}.json').read_bytes())
    if accepted['job_id']!=job or accepted['config_sha256']!=sha(bundle/'config.json') or not queue(job,shard)['ready']:
        raise ValueError('accepted handle/config/allocation gate failed')
    old,context=validate(config,bundle,native=True);reader,w,data,_,_,free,seed,_=context
    original=json.loads((OLD/'config.json').read_bytes());direction=np.asarray(original['branch_direction_36'])
    out.mkdir();rows=[]
    for point in config['points']:
        if point['id'] not in POLICY['shards'][shard]:continue
        d=out/f"point_{point['id']:03d}"
        evaluator=old.faces.canary.CalculatorModeEvaluator(data['chart'],
            old.study.make_factory(parameters={'source_directory':str(old.study.ASSETS)},command=old.study.COMMAND),
            workdir=d,calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
            minimum_distance_A=1.6,result_validator=old.study.result_validator,
            validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',max_evaluation_attempts=12)
        evaluate=old.faces.canary.InternalEvaluator(evaluator,data['basis'],data['center_energy'],data['center_gradient'])
        def guard(x):
            try:old.geometry(data,data['basis']@x,point['q_A']);return True
            except ValueError:return False
        row=dict(id=point['id'],q_A=point['q_A'],training_parent_id=point['training_parent_id'],result=None)
        try:
            row['result']=old.solver.relax_branch_candidate(evaluate,np.asarray(point['starting_internal_coordinates_A']),
                free,direction,initial_hessian=seed,valid=guard,gradient_tolerance=.03,
                branch_gradient_tolerance=1e-4,initial_trust_radius=.05,max_iterations=200)
            row['status']=row['result']['status']
        except RuntimeError as exc:
            if str(exc)!='mode evaluation attempt budget exhausted; calculator was not called' or list(d.glob('eval-*/failure.json')):raise
            row.update(status='attempt_budget_exhausted_no_retry',reason=str(exc))
        native=[];sources=[]
        for path in sorted(d.glob('eval-*/result.json')):
            x=np.asarray(json.loads(path.read_bytes())['coordinates']);old.geometry(data,x,point['q_A'])
            record,source=reader.native_record(w,data,x,path,d)
            native.append(record);sources.append(source)
        if len(native)!=evaluator.n_new_evaluations or evaluator.n_evaluation_attempts>12:raise ValueError('charged native accounting differs')
        row.update(reserved_attempts=evaluator.n_evaluation_attempts,new_native_SCFs=len(native),native_sources=sources)
        if row['result'] is not None:
            matches=[r for r in native if np.array_equal(r['coordinates'],data['basis']@np.asarray(row['result']['coordinates']))]
            if len(matches)!=1:raise ValueError('unique native final candidate required')
            g=data['basis'].T@np.asarray(matches[0]['gradient']);norm=float(np.linalg.norm(free.T@g));dg=float(direction@g)
            if not np.isclose(norm,row['result']['released_gradient_norm'],atol=2e-12,rtol=0) or not np.isclose(dg,row['result']['branch_gradient'],atol=2e-12,rtol=0):
                raise ValueError('final full native dual gradient differs')
            row['native_final']=matches[0]
        old.study.exclusive_json(out/f"point_{point['id']:03d}.json",row);rows.append(row)
        print(json.dumps(dict(id=point['id'],status=row['status'],SCFs=len(native))),flush=True)
    if [r['id'] for r in rows]!=POLICY['shards'][shard] or sum(r['reserved_attempts'] for r in rows)>48:
        raise ValueError('four original targets and48-attempt shard cap required')
    old.study.exclusive_json(out/'completed.json',dict(status='four_P8_targets_processed_stability_unchecked',
        job_id=job,shard=shard,config_sha256=sha(bundle/'config.json'),rows=rows,
        new_independent_targets_processed=4,physical_plus005_generated_or_read=False,new_training_candidates=0))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=('prepare','preflight','run'))
    parser.add_argument('--bundle',type=Path,required=True);parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--config-sha256');parser.add_argument('--shard',type=int,choices=(0,1));parser.add_argument('--output',type=Path)
    args=parser.parse_args();old,_=old_worker()
    if args.mode=='prepare':old.study.exclusive_json(args.config,prepare(args.bundle))
    else:
        if sha(args.config)!=args.config_sha256:raise ValueError('exact frozen configuration required')
        config=json.loads(args.config.read_bytes())
        if args.mode=='preflight':old.study.exclusive_json(args.output,preflight(config,args.bundle))
        else:run(config,args.bundle,args.shard)
