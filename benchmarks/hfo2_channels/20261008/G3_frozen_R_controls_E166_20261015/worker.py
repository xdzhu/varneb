"""Eight matched frozen-R static controls, no optimizer or dispatcher.

Four axial Q points per existing face; R=0 about its original sampled center.
Prescribed substrate, old six physical inputs and cached Q0 are unchanged.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT));sys.dont_write_bytecode=True
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from scripts import hfo2_G3_conditional_faces as faces
from scripts.prepare_hfo2_channel_work_probes import write_structure
from scripts.prepare_hfo2_gamma_and_seeds import geometry_roundtrip
from vcneb.mode_evaluator import _fingerprint

study=faces.study
REL='benchmarks/hfo2_channels/20261008/G3_frozen_R_controls_E166_20261015'
POLICY=dict(reservation_id='E166-four-frozen-R-axial-training-controls-per-face',
    channel_order=['preserving','T'],points_per_channel=4,total_new_SCF_attempts=8,
    attempts_per_point=1,first_batch_nodes_per_face=25,new_training_slots_per_face=4,
    training_slots_committed_after_reservation_per_face=29,training_slots_remaining_per_face=52,
    independent_slots_remaining_per_face=8,total_training_allowance_per_face=81,
    pressure_GPa=0.,strain=0.,ordinary_NEB_fmax_eV_A=.10,minimum_distance_A=1.6,
    allocation_nodes=1,genuine_MPI_ranks=32,allocation_hours=1,
    maximum_committed_study_allocations=2,job_name='hfo2-G3-E166-frozen',climb=False)
PREFIX_PINS={'T':'1a7b3530bb8d389196b9ad747df4f642e405ca1cdbb2d8da38d0ad529fefeda6',
    'preserving':'3c2c90ab703c2e1cc7c9c699daee5e0ded8f0d54d1b5d038c33a862e1b89b8ee'}
DECISION_SHA='c2e51c37b185a00676a6ff21660ea132384c0a22cfac2cebd18e30a5e5a9b745'
PARENTS={'T':'28862319','preserving':'28862441'}
CODE=(*faces.canary.CODE,'scripts/hfo2_G3_conditional_faces.py',REL+'/worker.py')


def load_data():
    faces.canary.pinned(faces.BASE/'G3_quadratic_agreement_E165_20261015/analysis.json',DECISION_SHA)
    data={}
    for channel,pin in PREFIX_PINS.items():
        d=faces.load_inputs(channel)
        prefix=faces.canary.pinned(faces.BASE/'G3_face_prefix_audit_E163_20261015'/f'{channel}_native_prefix.json',pin)
        matched={r['id']:r for r in prefix['candidate_rows']}
        if any(matched[i]['candidate'] is None or matched[i]['q_A']!=faces.ordered_points()[i]['q_A']
               or matched[i]['replay']!=dict(original_order_cache_only_replay_passed=True,new_DFT=0)
               for i in range(1,5)):
            raise ValueError('four original matched axial gradient candidates required')
        d['matched']=matched;data[channel]=d
    return data


def prepare():
    data=load_data();channels={}
    for channel in POLICY['channel_order']:
        d=data[channel];points=[]
        for spec in faces.ordered_points()[1:5]:
            internal=d['axes']@spec['q_A'];full=d['basis']@internal
            if full[-1]!=0 or not faces.canary.geometry_guard(d['chart'],d['basis'],internal):
                raise ValueError('unsafe frozen-R point; no repair or DFT')
            points.append(dict(id=spec['id'],q_A=spec['q_A'],coordinates_full_A=full.tolist(),
                matched_released_energy_eV_cell=d['matched'][spec['id']]['candidate']['energy_eV_cell']))
        channels[channel]=points
    return dict(schema_version=1,purpose='matched_frozen_R_training_control_not_relaxation_or_prediction',
        policy=deepcopy(POLICY),points=channels,original_manifest_sha256=faces.canary.PARENT_SHA,
        original_six_input_sha256=deepcopy(study.CONTRACT),original_binary_sha256=study.BINARY_SHA,
        matched_native_prefix_sha256=deepcopy(PREFIX_PINS),decision_sha256=DECISION_SHA,
        normalized_source_sha256={n:study.text_sha(ROOT/n) for n in CODE},
        wrapper_text_sha256=study.text_sha(ROOT/REL/'run.slurm'),
        zero_reference_new_SCFs=0,physical_reserved_geometry_generated_or_read=False,
        independent_test_or_stable_conditional_minima_or_TS_certified=False)


def validate(config):
    expected=prepare();left=deepcopy(config);right=deepcopy(expected)
    for channel in POLICY['channel_order']:
        if len(left['points'][channel])!=4:raise ValueError('exactly four registered frozen controls required')
        for row,ref in zip(left['points'][channel],right['points'][channel]):
            x=np.asarray(row.pop('coordinates_full_A'),dtype=float)
            target=np.asarray(ref.pop('coordinates_full_A'))
            if (x.shape!=(40,) or not np.isfinite(x).all() or x[-1]!=0
                    or not np.allclose(x,target,atol=2e-12,rtol=0)):
                raise ValueError('frozen Q/R coordinates or exact control differ')
    if json.dumps(left,sort_keys=True,allow_nan=False)!=json.dumps(right,sort_keys=True,allow_nan=False):
        raise ValueError('frozen point/physical/source/budget contract differs')
    return load_data()


def admission_ready(parent_rows,queue,own_job=None):
    """No inference from a missing handle; failed/timed-out parents need review."""
    if set(parent_rows)!=set(PARENTS.values()):raise ValueError('both original accounting handles required')
    if any(not isinstance(v,(tuple,list)) or len(v)!=2 for v in parent_rows.values()):
        raise ValueError('exact parent state/exit observations required')
    observed=[tuple(v) for v in parent_rows.values()]
    good={'RUNNING','PENDING','COMPLETED'}
    if any(state not in good or (state=='COMPLETED' and code!='0:0') for state,code in observed):
        return False
    if not any(row==('COMPLETED','0:0') for row in observed):return False
    if any(r['CPUs']!=32 or r['nodes']!=1 for r in queue):raise ValueError('study resource identity differs')
    if own_job is None:return len(queue)<2
    return len(queue)<=2 and sum(r['job_id']==own_job and r['name']==POLICY['job_name'] for r in queue)==1


def scheduler_snapshot(own_job=None):
    body=subprocess.run(['sacct','-n','-P','-X','-j',','.join(PARENTS.values()),
        '--format=JobIDRaw,State,ExitCode'],check=True,capture_output=True,text=True).stdout
    parents={}
    for line in body.splitlines():
        f=line.strip().split('|')
        if f[0] in PARENTS.values():parents[f[0]]=tuple(f[1:3])
    bodyq=subprocess.run(['squeue','-h','-u','iai806','-p','hfacnormal01','-o','%i|%.100j|%T|%C|%D'],
        check=True,capture_output=True,text=True).stdout
    queue=[]
    for line in bodyq.splitlines():
        f=[s.strip() for s in line.split('|')]
        if len(f)!=5:raise ValueError('unparsed queue observation; no admission inferred')
        if f[1].startswith('hfo2-G3-'):
            queue.append(dict(job_id=f[0],name=f[1],state=f[2],CPUs=int(f[3]),nodes=int(f[4])))
    return dict(parent_states=parents,committed_G3_queue=queue,
                ready=admission_ready(parents,queue,own_job),own_job_id=own_job,scheduler_writes=0)


def native_preflight(config,decks=None):
    if decks is not None and (decks.exists() or any(decks.resolve().is_relative_to(p.resolve())
            for p in (ROOT,study.ASSETS,study.RUNS/'20261015-varneb-G3-faces-E161-r1'))
            or decks.resolve()==study.ROOT.parent/'frozen_R'):
        raise ValueError('fresh inert decks must stay outside protected source/assets/live study and calculation target')
    data=validate(config)
    for channel in POLICY['channel_order']:
        cfg=faces.prepare(channel)
        faces.native_preflight(cfg,channel)
    previews=[]
    if decks is not None:
        if decks.exists():raise FileExistsError('fresh inert native input-deck directory required')
        decks.mkdir(parents=True)
        for channel in POLICY['channel_order']:
            for point in config['points'][channel]:
                directory=decks/channel/f"point_{point['id']:03d}";directory.mkdir(parents=True)
                for name in study.CONTRACT:
                    if name in ('INPUT','KPT'):shutil.copyfile(study.ASSETS/name,directory/name)
                    else:(directory/name).symlink_to(study.ASSETS/name)
                atoms=data[channel]['chart'].to_atoms(point['coordinates_full_A'])
                write_structure(directory/'STRU',atoms);minimum=geometry_roundtrip(directory/'STRU',atoms)
                hashes={n:study.sha256(directory/n) for n in (*study.CONTRACT,'STRU')}
                if any(hashes[n]!=h for n,h in study.CONTRACT.items()):raise ValueError('inert physical deck drift')
                previews.append(dict(channel=channel,point_id=point['id'],input_sha256=hashes,
                    minimum_distance_A=minimum,original_asset_symlinks_not_copies=True))
    return data,dict(status='frozen_R_native_preflight_no_DFT',preview_decks=previews,
        admission=scheduler_snapshot(),new_DFT_calls=0,physical_reserved_geometry_generated_or_read=False)


def native_point(data,point,directory):
    """Reparse the one original native SCF and its conservative reservation."""
    x=np.asarray(point['coordinates_full_A']);digest=hashlib.sha256(x.tobytes()).hexdigest()
    path=directory/f'eval-000000-{digest[:12]}'/'result.json'
    if (len(list(directory.glob('eval-*')))!=1 or list(directory.glob('eval-*/failure.json'))
            or len(list(directory.glob('attempt-*.json')))!=1):
        raise ValueError('exactly one complete static attempt per frozen point required')
    row=json.loads(path.read_bytes());contract=json.loads((directory/'contract.json').read_bytes())
    reference=data['chart'].reference_atoms
    expected=dict(schema_version=1,coordinate_convention=data['chart'].coordinate_convention,
        reference=dict(numbers=reference.numbers.tolist(),cell_A=reference.cell.array.tolist(),
            scaled_positions=reference.get_scaled_positions(wrap=False).tolist(),pbc=reference.pbc.tolist()),
        calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
        validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',pressure_eV_per_A3=0.,
        minimum_distance_A=1.6,mechanical_boundary=data['chart'].mechanical_contract)
    fingerprint=_fingerprint(expected)
    budget=dict(schema_version=1,max_evaluation_attempts=1,physical_contract_sha256=fingerprint)
    marker=dict(schema_version=1,evaluation_index=0,coordinate_sha256=digest,
        coordinates=x.tolist(),budget_sha256=_fingerprint(budget))
    if (contract!=dict(expected,contract_sha256=fingerprint)
            or _fingerprint({k:v for k,v in contract.items() if k!='contract_sha256'})!=fingerprint
            or json.loads((directory/'evaluation-budget.json').read_bytes())!=budget
            or json.loads((directory/'attempt-000000.json').read_bytes())!=marker
            or row['status']!='complete' or row['coordinate_sha256']!=digest
            or row['contract_sha256']!=fingerprint or not np.array_equal(row['coordinates'],x)):
        raise ValueError('native static contract/reservation/coordinate identity differs')
    atoms=data['chart'].to_atoms(x)
    atoms.calc=SinglePointCalculator(atoms,**study.audited_results(path.parent/'scf_000000'))
    checked=study.result_validator(0,atoms,path.parent);e,g=data['chart'].energy_and_gradient(x,atoms)
    if (checked!=row['validation'] or e!=row['enthalpy_eV']
            or not np.allclose(g,row['gradient'],atol=2e-12,rtol=0)
            or not np.array_equal(atoms.get_stress(voigt=False),row['stress_eV_per_A3'])):
        raise ValueError('frozen static raw EFS/geometry/gradient differs')
    return dict(point_id=point['id'],q_A=point['q_A'],coordinates_full_A=row['coordinates'],
        frozen_energy_eV_cell=e,physical_gradient_eV_A=row['gradient'],native_validation=checked,
        result_sha256=study.sha256(path),matched_released_energy_eV_cell=point['matched_released_energy_eV_cell'],
        frozen_minus_released_meV_cell=(e-point['matched_released_energy_eV_cell'])*1000,
        new_SCFs=1,attempts=1,signed_difference_not_forced_nonnegative=True)


def compute(config,data,workdir,*,admission=None):
    if (set(config['points'])!=set(POLICY['channel_order']) or any(
            len(config['points'][c])!=4 or [p['id'] for p in config['points'][c]]!=[1,2,3,4]
            for c in POLICY['channel_order'])):
        raise ValueError('exactly eight original axial controls; no extra static calls')
    workdir.mkdir();rows=[]
    if admission is not None:study.exclusive_json(workdir/'admission.json',admission)
    for channel in POLICY['channel_order']:
        for point in config['points'][channel]:
            directory=workdir/channel/f"point_{point['id']:03d}"
            evaluator=faces.canary.CalculatorModeEvaluator(data[channel]['chart'],
                study.make_factory(parameters={'source_directory':str(study.ASSETS)},command=study.COMMAND),
                workdir=directory,calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
                minimum_distance_A=1.6,result_validator=study.result_validator,
                validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',max_evaluation_attempts=1)
            evaluator(np.asarray(point['coordinates_full_A']))
            if evaluator.n_new_evaluations!=1 or evaluator.n_evaluation_attempts!=1:
                raise ValueError('one and only one native call per frozen point required')
            row=dict(channel=channel,**native_point(data[channel],point,directory))
            study.exclusive_json(directory.parent/f"point_{point['id']:03d}.json",row);rows.append(row)
            print(json.dumps(dict(channel=channel,point_id=point['id'],new_SCFs=len(rows))),flush=True)
    result=dict(status='eight_frozen_R_statics_completed_not_PES_or_independent_forecast',rows=rows,
        policy=POLICY,new_SCFs=8,total_attempts=8,cached_original_zero_new_SCFs=0,
        job_id=os.environ.get('SLURM_JOB_ID'),source_original_centers_and_released_origins_reused=True,
        independent_material_prediction_or_released_stability_or_TS_certified=False,
        physical_reserved_geometry_generated_or_read=False)
    study.exclusive_json(workdir/'controls.json',result)
    return result


def run(config,workdir):
    if (os.environ.get('RUN_DFT')!='1' or os.environ.get('SLURM_JOB_PARTITION')!='hfacnormal01'
            or os.environ.get('SLURM_JOB_NUM_NODES')!='1' or os.environ.get('SLURM_NTASKS')!='1'
            or os.environ.get('SLURM_CPUS_PER_TASK')!='32' or not os.environ.get('SLURM_JOB_ID')
            or os.environ.get('SLURM_JOB_NAME')!=POLICY['job_name']
            or any(os.environ.get(k)!='1' for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'))):
        raise ValueError('authorized hfacnormal01/one-node32MPI/one-thread allocation required')
    if (workdir.resolve()!=study.ROOT.parent/'frozen_R' or not workdir.resolve().is_relative_to(study.RUNS)
            or workdir.exists() or ROOT.resolve()!=study.ROOT.resolve()):
        raise ValueError('fresh frozen_R sibling of immutable source required')
    admission=scheduler_snapshot(os.environ['SLURM_JOB_ID'])
    if not admission['ready']:raise RuntimeError('original allocations/resource gate not ready; no DFT')
    data,_=native_preflight(config)
    result=compute(config,data,workdir,admission=admission)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('prepare','preflight','run'));p.add_argument('--config',type=Path,required=True)
    p.add_argument('--config-sha256');p.add_argument('--decks',type=Path);p.add_argument('--workdir',type=Path)
    args=p.parse_args()
    if args.mode=='prepare':
        if args.config.exists():raise FileExistsError('fresh single frozen-control reservation required')
        config=prepare();args.config.parent.mkdir(parents=True,exist_ok=True);study.exclusive_json(args.config,config)
        print(json.dumps(dict(status='prepared_eight_frozen_R_controls_no_DFT',new_DFT=0)))
    else:
        if not args.config_sha256:raise ValueError('exact frozen-control config SHA required')
        config=faces.canary.pinned(args.config,args.config_sha256)
        if args.mode=='preflight':
            _,report=native_preflight(config,args.decks)
            print(json.dumps(report))
        else:
            if args.workdir is None:raise ValueError('fresh static output directory required')
            result=run(config,args.workdir)
            print(json.dumps({k:result[k] for k in ('status','new_SCFs')}))


if __name__=='__main__':main()
