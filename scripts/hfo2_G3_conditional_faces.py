"""First bounded5x5 conditional-candidate sheets at the two existing centers.

No dispatcher, physical-input repair, whole band, third Hessian center or
physical strain holdout. One continued candidate branch is not a global PES.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import zipfile

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator

from scripts import hfo2_G3_conditional_canary as canary

study=canary.study
BASE=study.ROOT/'benchmarks/hfo2_channels/20261008'
CURVATURE=BASE/'G3_complete_curvature_E159_20261014'
COMPLETED=BASE/'G3_conditional_canary_E160_20261014'
PINS=dict(T=('a182e4c386889880c77c2069e939cd65665724d98d3838e251a077204cad1c3b',
    'f77ebe3708f89b14dbbda053505f56fbf26672a2830bd0231680d6798d59a739',
    'b4097af45e0de951c32e550141d745a2d73b6e13eeb4aa8fb3a62be302afa685',
    '82fc05455912488051f52da2396606cbc060fb9be1c06ab42e07487768b2a44f'),
    preserving=('5d137fd0e676ad4fb7e87076d76b53e618ab9de25bc86228e02d4bcf0285a778',
    'e8ce5429c24fb443f4a651eef14aac8cf0190a7399cb19a0abef838658be12b5',
    '336071744885656b25515f28ca9a6d6b4ecd09f9b411a739827ec7d1e56e5403',
    'c1955f42e1d9142235971a83eae0fc52c18b4c65521ec5e731e3ed1d2d4c28f9'))
POLICY=dict(axis_A=[-.2,-.1,0.,.1,.2],modeling_nodes_total=25,reused_origin_nodes=1,
    maximum_new_attempts_per_node=12,maximum_new_attempts_per_face=288,
    independent_test_points_used=0,remaining_training_node_allowance=56,
    remaining_independent_point_allowance=8,allocation_hours=12,
    optimizer='safeguarded_bfgs',gradient_tolerance_eV_A=.03,initial_trust_radius_A=.05,
    maximum_internal_displacement_A=.5,minimum_distance_A=1.6,maximum_iterations=200,
    starting_policy='nearest-safe-accepted-Q; ties by fixed point order; one continuation start',
    pressure_GPa=0.,physical_strain=0.,climb=False)


def ordered_points():
    indices=sorted(((i,j) for i in range(-2,3) for j in range(-2,3)),
                   key=lambda p:(abs(p[0])+abs(p[1]),p))
    return [dict(id=k,index=list(ij),q_A=[POLICY['axis_A'][i+2] for i in ij]) for k,ij in enumerate(indices)]


def load_inputs(channel):
    cp,sp,ap,zp=PINS[channel]
    manifest=canary.pinned(BASE/'G3_curvature_E150_20261014/manifest.json',canary.PARENT_SHA)
    old=canary.pinned(CURVATURE/(channel+'_config.json'),sp)
    native=canary.pinned(COMPLETED/(channel+'_native_audit.json'),ap)
    with zipfile.ZipFile(CURVATURE/'native_reports.zip') as archive:
        raw=archive.read(channel+'_native.json')
        if hashlib.sha256(raw).hexdigest()!=cp:raise ValueError('original complete curvature bytes differ')
        curvature=json.loads(raw)
    chart,basis,gradient,matrices=canary.validate_config(old,manifest,channel,curvature,cp)
    path=COMPLETED/(channel+'_worker_outputs.zip')
    if study.sha256(path)!=zp:raise ValueError('original completed canary archive differs')
    with zipfile.ZipFile(path) as archive:
        summary=json.loads(archive.read(channel+'/canary.json'));outcome=summary['arms'][1]
        point=outcome['point']
        if (summary['settings']!=canary.SETTINGS or outcome['arm']!='curvature_seed' or point is None
                or outcome['status']!='gradient_converged_candidate_curvature_unchecked'
                or point['orthogonal_gradient_norm_eV_A']>POLICY['gradient_tolerance_eV_A']):
            raise ValueError('original seeded origin candidate required')
        matches=[]
        for name in archive.namelist():
            if name.startswith(channel+'/curvature_seed/eval-') and name.endswith('/result.json'):
                raw=archive.read(name);row=json.loads(raw)
                if np.array_equal(row['coordinates'],point['coordinates_full_A']):matches.append((raw,row))
    if len(matches)!=1:raise ValueError('origin must match one original complete native record')
    raw,row=matches[0]
    if hashlib.sha256(raw).hexdigest() not in {r['result_sha256'] for r in native['native_record_sources']}:
        raise ValueError('origin native audit linkage differs')
    axes=np.asarray(old['internal_axes_36_by_2'])
    if not np.allclose(basis@point['coordinates_internal_A'],row['coordinates'],atol=2e-12,rtol=0):
        raise ValueError('origin embedding differs')
    return dict(manifest=manifest,old_config=old,curvature=curvature,chart=chart,basis=basis,
        axes=axes,hessian=matrices[0],center_energy=manifest['centers'][channel]['raw_source']['energy_eV_cell'],
        center_gradient=gradient,origin_point=point,origin_record=row)


def choose_start(data,accepted,q):
    axes=data['axes']
    for row in sorted(accepted,key=lambda r:(float(np.linalg.norm(np.asarray(r['q_A'])-q)),r['id'])):
        previous=np.asarray(row['point']['coordinates_internal_A'])
        start=axes@q+previous-axes@(axes.T@previous)
        if canary.geometry_guard(data['chart'],data['basis'],start):return row['id'],start
    return None


def prepare(channel):
    data=load_inputs(channel)
    origin=dict(id=0,q_A=[0.,0.],point=data['origin_point'])
    if any(choose_start(data,[origin],p['q_A']) is None for p in ordered_points()):
        raise ValueError('prospective grid contains an unsafe origin-continued start; no domain repair')
    return dict(schema_version=1,purpose='first_training_conditional_candidate_sheet_not_certified_PES',
        channel=channel,policy=deepcopy(POLICY),points=ordered_points(),
        parent_manifest_sha256=canary.PARENT_SHA,input_artifact_sha256=list(PINS[channel]),
        source_text_sha256={n:study.text_sha(study.ROOT/n) for n in (*canary.CODE,'scripts/hfo2_G3_conditional_faces.py')},
        sampling_axes_and_rectangle_locked_for_this_batch=True,
        final_predictive_representative_or_domain_locked=False,physical_reserved_geometry_generated_or_read=False,
        released_stability_or_competing_branches_or_TS_certified=False)


def validate(config,channel):
    expected=prepare(channel)
    if json.dumps(config,sort_keys=True,allow_nan=False)!=json.dumps(expected,sort_keys=True,allow_nan=False):
        raise ValueError('registered grid/source/physics/budget differs')
    return load_inputs(channel)


def native_preflight(config,channel):
    data=validate(config,channel)
    canary.native_canary_preflight(data['old_config'],data['manifest'],channel,data['curvature'],PINS[channel][0])
    row=data['origin_record'];directory=Path(row['validation']['raw_directory']).parent
    atoms=data['chart'].to_atoms(row['coordinates'])
    atoms.calc=SinglePointCalculator(atoms,**study.audited_results(directory/'scf_000000'))
    checked=study.result_validator(0,atoms,directory)
    e,g=data['chart'].energy_and_gradient(row['coordinates'],atoms)
    if (checked!=row['validation'] or e!=row['enthalpy_eV']
            or not np.allclose(g,row['gradient'],atol=2e-12,rtol=0)):
        raise ValueError('cached conditional origin native EFS/gradient differs')
    return data


def compute_face(config,data,workdir):
    """Finite point orchestration; native entry point performs allocation checks."""
    workdir.mkdir()
    basis,axes=data['basis'],data['axes']
    plane=canary.ModePlane(np.zeros(36),axes,np.ones(36),np.eye(2),
        ('coupled internal1','coupled internal2'),'local length','G3-first-conditional-sheet')
    origin=dict(**config['points'][0],status='cached_origin_gradient_candidate_curvature_unchecked',
        point=data['origin_point'],new_completed_SCFs=0,reserved_attempts=0,
        physical_gradient_eV_A=data['origin_record']['gradient'],origin_native_source=data['origin_record']['validation'])
    study.exclusive_json(workdir/'point_000.json',origin)
    rows=[origin];accepted=[origin]
    for spec in config['points'][1:]:
        row=dict(**spec,status='no_safe_continuation_start_no_DFT',point=None,
                 new_completed_SCFs=0,reserved_attempts=0)
        selected=choose_start(data,accepted,spec['q_A'])
        if selected is not None:
            parent,start=selected;row['continuation_parent_id']=parent
            row['starting_internal_coordinates_A']=start.tolist()
            directory=workdir/f"node_{spec['id']:03d}"
            evaluator=canary.CalculatorModeEvaluator(data['chart'],
                study.make_factory(parameters={'source_directory':str(study.ASSETS)},command=study.COMMAND),
                workdir=directory,calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
                minimum_distance_A=POLICY['minimum_distance_A'],result_validator=study.result_validator,
                validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',
                max_evaluation_attempts=POLICY['maximum_new_attempts_per_node'])
            evaluate=canary.InternalEvaluator(evaluator,basis,data['center_energy'],data['center_gradient'])
            screens=dict(checks=0,rejected=0)
            def guard(x):
                screens['checks']+=1
                ok=canary.geometry_guard(data['chart'],basis,x);screens['rejected']+=int(not ok)
                return ok
            try:
                point=canary.relax_orthogonal_at_q(plane,spec['q_A'],evaluate,starts=[start],include_frozen_start=False,
                    optimizer='safeguarded_bfgs',trial_validator=guard,initial_hessian=data['hessian'],
                    initial_trust_radius=POLICY['initial_trust_radius_A'],
                    gradient_tolerance=POLICY['gradient_tolerance_eV_A'],max_iterations=POLICY['maximum_iterations'])
                row.update(status='gradient_candidate_release_stability_unchecked',point=dict(
                    coordinates_internal_A=point.coordinates.tolist(),coordinates_full_A=(basis@point.coordinates).tolist(),
                    energy_eV_cell=point.energy,orthogonal_gradient_norm_eV_A=point.orthogonal_gradient_norm,
                    callback_evaluations=point.total_n_evaluations))
            except RuntimeError as exc:
                known=str(exc) in ('mode evaluation attempt budget exhausted; calculator was not called',
                                  'no starting point reached the orthogonal-gradient tolerance')
                if (list(directory.glob('eval-*/failure.json')) or not known
                        or (str(exc).startswith('mode evaluation attempt budget exhausted')
                            and evaluator.n_evaluation_attempts!=POLICY['maximum_new_attempts_per_node'])):raise
                row.update(status='optimizer_not_converged',reason=str(exc))
            records=[json.loads(p.read_bytes()) for p in sorted(directory.glob('eval-*/result.json'))]
            if (len(records)!=evaluator.n_new_evaluations
                    or evaluator.n_evaluation_attempts>POLICY['maximum_new_attempts_per_node']):
                raise ValueError('native point accounting differs')
            row.update(new_completed_SCFs=evaluator.n_new_evaluations,reserved_attempts=evaluator.n_evaluation_attempts,
                callback_requests=evaluate.callback_requests,center_cache_hits=evaluate.center_cache_hits,
                geometry_screens=screens,sum_completed_native_call_seconds=sum(r['validation']['elapsed_seconds'] for r in records))
            if row['point'] is not None:
                final=[r for r in records if np.array_equal(r['coordinates'],row['point']['coordinates_full_A'])]
                if len(final)!=1:raise ValueError('candidate must match an exact complete native SCF')
                row['physical_gradient_eV_A']=final[0]['gradient']
                row['native_source']=final[0]['validation'];accepted.append(row)
        study.exclusive_json(workdir/f"point_{spec['id']:03d}.json",row);rows.append(row)
        print(json.dumps(dict(point_id=spec['id'],q_A=spec['q_A'],status=row['status'],
            new_SCFs=row['new_completed_SCFs'],completed_nodes=len(rows))),flush=True)
    attempts=sum(r['reserved_attempts'] for r in rows)
    if attempts>POLICY['maximum_new_attempts_per_face']:raise ValueError('face budget exceeded')
    result=dict(status='completed_first_5x5_candidate_sheet_not_certified_PES',channel=config['channel'],
        policy=POLICY,rows=rows,gradient_candidates=len(accepted),new_reserved_attempts=attempts,
        new_completed_SCFs=sum(r['new_completed_SCFs'] for r in rows),cached_origin_new_SCFs=0,
        allocation_core_hours=None,original_curvature_canary_and_G2_preparation_costs_additional=True,
        released_stability_global_branch_TS_contour_or_prediction_certified=False,
        physical_reserved_geometry_generated_or_read=False)
    study.exclusive_json(workdir/'face.json',result)
    return result


def run(config,channel,workdir):
    if (os.environ.get('RUN_DFT')!='1' or os.environ.get('SLURM_JOB_PARTITION')!='hfacnormal01'
            or os.environ.get('SLURM_JOB_NUM_NODES')!='1' or os.environ.get('SLURM_NTASKS')!='1'
            or os.environ.get('SLURM_CPUS_PER_TASK')!='32' or not os.environ.get('SLURM_JOB_ID')
            or any(os.environ.get(n)!='1' for n in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'))):
        raise ValueError('authorized one-node hfacnormal01/32MPI allocation required')
    if (workdir.resolve()!=study.ROOT.parent/channel or workdir.exists()
            or not workdir.resolve().is_relative_to(study.RUNS)):
        raise ValueError('fresh face directory beside pinned source required')
    return compute_face(config,native_preflight(config,channel),workdir)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('prepare','preflight','run'))
    p.add_argument('--channel',choices=tuple(PINS),required=True)
    p.add_argument('--config',type=Path,required=True);p.add_argument('--config-sha256')
    p.add_argument('--workdir',type=Path)
    args=p.parse_args()
    if args.mode=='prepare':
        if args.config.exists():raise FileExistsError('fresh face configuration required')
        result=prepare(args.channel);args.config.parent.mkdir(parents=True,exist_ok=True)
        study.exclusive_json(args.config,result)
        print(json.dumps(dict(status='prepared_first_5x5_training_candidate_sheet',new_DFT=0)))
    else:
        if not args.config_sha256:raise ValueError('exact face config SHA required')
        config=canary.pinned(args.config,args.config_sha256)
        if args.mode=='preflight':
            native_preflight(config,args.channel)
            print(json.dumps(dict(status='native_5x5_preflight_passed_no_DFT',channel=args.channel,new_DFT=0)))
        else:
            if args.workdir is None:raise ValueError('fresh output directory required')
            result=run(config,args.channel,args.workdir)
            print(json.dumps({k:result[k] for k in ('status','gradient_candidates','new_completed_SCFs')}))


if __name__=='__main__':main()
