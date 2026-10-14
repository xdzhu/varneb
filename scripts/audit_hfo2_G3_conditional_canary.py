"""Read-only completion audit of the two E159 same-input DFT canaries.

Completed trial SCFs are not accepted optimizer iterates. Only the worker's
explicit final point, matched to original native data, can satisfy its common
gradient target. No automatic submission, optimizer restart or TS claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator

from scripts import hfo2_G3_conditional_canary as worker
from vcneb.mode_evaluator import _fingerprint

study=worker.study
RUNTIME_SHA='76ad3da30ae0642b72b9cda2e78d91c8eb0196112c17137769f4389c2d9ebfc8'
PINS=dict(T='a182e4c386889880c77c2069e939cd65665724d98d3838e251a077204cad1c3b',
          preserving='5d137fd0e676ad4fb7e87076d76b53e618ab9de25bc86228e02d4bcf0285a778')
CONFIGS=dict(T='f77ebe3708f89b14dbbda053505f56fbf26672a2830bd0231680d6798d59a739',
             preserving='e8ce5429c24fb443f4a651eef14aac8cf0190a7399cb19a0abef838658be12b5')


class CacheOnly:
    """Strict original-order replay. A missing point never launches a calculator."""
    fixed_coordinate_covectors=np.eye(40)[:,-1:]

    def __init__(self,records):
        self.records=records;self.seen={}

    def __call__(self,x):
        key=hashlib.sha256(np.asarray(x).tobytes()).hexdigest()
        if key not in self.seen:
            index=len(self.seen)
            if index==worker.SETTINGS['maximum_attempts_per_arm']:
                raise RuntimeError('mode evaluation attempt budget exhausted; calculator was not called')
            if index>=len(self.records) or not np.array_equal(x,self.records[index]['coordinates']):
                raise ValueError('cache-only optimizer replay requires a missing or differently ordered trial; no DFT')
            self.seen[key]=self.records[index]
        row=self.seen[key]
        return row['enthalpy_eV'],np.array(row['gradient'],copy=True)


def replay_arm(chart,basis,config,energy,gradient,hessian,outcome,records):
    """Reexecute the pinned optimizer on cached E/g only, with no file writes."""
    oracle=CacheOnly(records)
    evaluate=worker.InternalEvaluator(oracle,basis,energy,gradient)
    plane=worker.ModePlane(np.zeros(36),np.asarray(config['internal_axes_36_by_2']),np.ones(36),np.eye(2),
        ('coupled internal1','coupled internal2'),'local length','G3-center-canary-only')
    screens=dict(checks=0,rejected=0)
    def guard(x):
        screens['checks']+=1
        ok=worker.geometry_guard(chart,basis,x);screens['rejected']+=int(not ok)
        return ok
    point=None;reason=None
    try:
        point=worker.relax_orthogonal_at_q(plane,[0.,0.],evaluate,optimizer='safeguarded_bfgs',
            trial_validator=guard,initial_trust_radius=worker.SETTINGS['initial_trust_radius_A'],
            initial_hessian=hessian if outcome['arm']=='curvature_seed' else None,
            gradient_tolerance=worker.SETTINGS['gradient_tolerance_eV_A'],
            max_iterations=worker.SETTINGS['maximum_iterations'])
    except RuntimeError as exc:
        reason=str(exc)
        if reason not in ('mode evaluation attempt budget exhausted; calculator was not called',
                          'no starting point reached the orthogonal-gradient tolerance'):raise
    if (len(oracle.seen)!=len(records) or evaluate.callback_requests!=outcome['callback_requests']
            or evaluate.center_cache_hits!=outcome['center_cache_hits']
            or screens!=outcome['trial_geometry_and_domain_screens']
            or (point is None)!=(outcome['point'] is None)):
        raise ValueError('cache-only replay differs from original optimizer decisions/counters')
    if point is None:
        if reason!=outcome['reason']:raise ValueError('cache-only replay stop reason differs')
    elif (not np.array_equal(point.coordinates,outcome['point']['coordinates_internal_A'])
            or point.energy!=outcome['point']['energy_eV_cell']
            or not np.isclose(point.orthogonal_gradient_norm,outcome['point']['orthogonal_gradient_norm_eV_A'],atol=2e-12,rtol=0)):
        raise ValueError('cache-only replay final candidate differs')
    return dict(original_order_and_cache_only_optimizer_replay_passed=True,new_DFT=0,
        replayed_callback_requests=evaluate.callback_requests,replayed_center_cache_hits=evaluate.center_cache_hits,
        replayed_geometry_screens=screens)


def arm_diagnostics(basis,axes,center_energy,center_gradient,outcome,records):
    """Arithmetic/bookkeeping only; caller must establish raw-data provenance."""
    cap=worker.SETTINGS['maximum_attempts_per_arm']
    n=len(records)
    if (not 0<=n<=cap or type(outcome['new_completed_SCFs']) is not int or outcome['new_completed_SCFs']!=n
            or type(outcome['reserved_attempts']) is not int or outcome['reserved_attempts']!=n
            or type(outcome['center_cache_hits']) is not int or outcome['center_cache_hits']<1
            or type(outcome['callback_requests']) is not int
            or outcome['callback_requests']<n+outcome['center_cache_hits']):
        raise ValueError('completed canary attempts/callback accounting differs')
    def norm(g):
        internal=basis.T@np.asarray(g)
        return float(np.linalg.norm(internal-axes@(axes.T@internal)))
    trials=[]
    for index,r in enumerate(records):
        x=np.asarray(r['coordinates']);g=np.asarray(r['gradient'])
        if (r['status']!='complete' or x.shape!=(40,) or g.shape!=(40,)
                or not np.isfinite(x).all() or not np.isfinite(g).all() or x[-1]!=0
                or not np.isfinite(r['enthalpy_eV'])
                or not np.allclose(basis@(basis.T@x),x,atol=2e-12,rtol=0)
                or not np.allclose(axes.T@(basis.T@x),0,atol=2e-12,rtol=0)):
            raise ValueError('finite same-Q internal trial with exact fixed control required')
        trials.append(dict(index=index,energy_eV_cell=r['enthalpy_eV'],
            released_gradient_norm_eV_A=norm(g),accepted_iterate_known=False))
    times=[r['validation']['elapsed_seconds'] for r in records]
    if any(not np.isfinite(t) or t<=0 for t in times):raise ValueError('finite positive native call times required')
    seconds=sum(times)
    if (not np.isfinite(seconds) or seconds<0
            or not np.isclose(seconds,outcome['sum_completed_native_call_seconds'],atol=1e-10,rtol=0)
            or not np.isclose(seconds*32/3600,outcome['completed_call_proxy_core_hours'],atol=1e-10,rtol=0)):
        raise ValueError('native completed-call cost differs')
    point=outcome['point'];final=None
    if point is None:
        if outcome['status']!='not_converged' or outcome.get('reason') not in (
                'mode evaluation attempt budget exhausted; calculator was not called',
                'no starting point reached the orthogonal-gradient tolerance'):
            raise ValueError('nonconvergence category differs; no last-trial promotion')
        if outcome['reason'].startswith('mode evaluation attempt budget exhausted') and n!=cap:
            raise ValueError('budget outcome requires actual cap exhaustion')
    else:
        if outcome['status']!='gradient_converged_candidate_curvature_unchecked':
            raise ValueError('candidate category differs')
        internal=np.asarray(point['coordinates_internal_A']);full=np.asarray(point['coordinates_full_A'])
        if (internal.shape!=(36,) or full.shape!=(40,) or not np.isfinite(internal).all()
                or not np.allclose(basis@internal,full,atol=2e-12,rtol=0)
                or not np.allclose(axes.T@internal,0,atol=2e-12,rtol=0)):
            raise ValueError('candidate coordinate or same-Q constraint differs')
        matches=[r for r in records if np.array_equal(r['coordinates'],full)]
        if not np.any(internal):
            e,g=center_energy,center_gradient
        elif len(matches)==1:
            e,g=matches[0]['enthalpy_eV'],matches[0]['gradient']
        else:raise ValueError('candidate must match a unique complete native record or the exact origin')
        residual=norm(g)
        if (point['energy_eV_cell']!=e or residual>worker.SETTINGS['gradient_tolerance_eV_A']
                or not np.isclose(residual,point['orthogonal_gradient_norm_eV_A'],atol=2e-12,rtol=0)
                or point['callback_evaluations']!=outcome['callback_requests']):
            raise ValueError('candidate energy/gradient/common target differs')
        final=dict(energy_eV_cell=e,released_gradient_norm_eV_A=residual,
                   coordinates_internal_A=internal.tolist(),curvature_stability_unchecked=True)
    return dict(arm=outcome['arm'],worker_status=outcome['status'],completed_SCFs=n,
        conservative_attempt_reservations=outcome['reserved_attempts'],final_candidate=final,
        native_trial_rows=trials,sum_audited_call_seconds=seconds,
        completed_call_proxy_core_hours=seconds*32/3600,
        callback_and_center_hit_counts_worker_declared_not_trace_replayed=True,
        last_complete_trial_is_not_assumed_accepted=True)


def audit(channel,curvature_report):
    root=study.ROOT.parent;workdir=root/channel
    config=worker.pinned(root/(channel+'_config.json'),CONFIGS[channel])
    runtime=study.ROOT/'runtime_files.json'
    if study.sha256(runtime)!=RUNTIME_SHA:raise ValueError('pinned runtime manifest differs')
    files=json.loads(runtime.read_bytes())
    if len(files)!=378 or any(study.sha256(study.ROOT/n)!=h for n,h in files.items()):
        raise ValueError('pinned runtime Python bytes differ')
    manifest=worker.pinned(study.ROOT/'benchmarks/hfo2_channels/20261008/G3_curvature_E150_20261014/manifest.json',worker.PARENT_SHA)
    report=worker.pinned(curvature_report,PINS[channel])
    chart,basis,energy,gradient,matrices=worker.native_canary_preflight(config,manifest,channel,report,PINS[channel])
    accepted=json.loads((root/(channel+'_accepted.json')).read_bytes())
    job=accepted['job_id']
    accounting=subprocess.run(['sacct','-n','-P','-X','-j',job,
        '--format=JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,NodeList,Partition'],
        check=True,capture_output=True,text=True).stdout.strip()
    fields=accounting.split('|')
    if len(fields)<7 or fields[:3]!=[job,'COMPLETED','0:0'] or fields[4]!='32' or fields[6]!='hfacnormal01':
        raise RuntimeError('same handle not yet normally complete; do not restart or call it failed: '+accounting)
    summary=json.loads((workdir/'canary.json').read_bytes())
    center=json.loads((workdir/'center.json').read_bytes())
    if (summary['status']!='bounded_conditional_canary_not_PES_TS_or_prediction'
            or summary['channel']!=channel or summary['job_id']!=job
            or summary['settings']!=worker.SETTINGS or center['source']!=manifest['centers'][channel]['raw_source']
            or center['new_center_SCFs']!=0 or center['energy_eV_cell']!=energy
            or not np.allclose(center['full_gradient_eV_A'],gradient,atol=2e-12,rtol=0)):
        raise ValueError('center/summary/policy differs')
    rows=[];sources=[]
    for arm in worker.SETTINGS['arms']:
        directory=workdir/arm
        contract=json.loads((directory/'contract.json').read_bytes())
        fingerprint=_fingerprint({k:v for k,v in contract.items() if k!='contract_sha256'})
        reference=chart.reference_atoms
        if (fingerprint!=contract['contract_sha256'] or contract['reference']!=dict(
                numbers=reference.numbers.tolist(),cell_A=reference.cell.array.tolist(),
                scaled_positions=reference.get_scaled_positions(wrap=False).tolist(),pbc=reference.pbc.tolist())
                or contract['coordinate_convention']!=chart.coordinate_convention
                or contract['mechanical_boundary']!=chart.mechanical_contract
                or contract['pressure_eV_per_A3']!=0
                or contract['calculator_id']!='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2'
                or contract['validation_id']!='E150-native-six-input-DSIZE32-EFS-one-SCF-v1'):
            raise ValueError('native evaluator contract differs')
        budget=json.loads((directory/'evaluation-budget.json').read_bytes())
        if budget!=dict(schema_version=1,max_evaluation_attempts=12,physical_contract_sha256=fingerprint):
            raise ValueError('original per-arm budget differs')
        paths=sorted(directory.glob('eval-*/result.json'))
        if list(directory.glob('eval-*/failure.json')) or len(list(directory.glob('eval-*')))!=len(paths):
            raise ValueError('failed/incomplete native attempt; not a complete benchmark')
        markers=sorted(directory.glob('attempt-*.json'))
        if len(markers)!=len(paths):raise ValueError('conservative reservations differ from complete records')
        records=[]
        for index,(path,marker) in enumerate(zip(paths,markers)):
            row=json.loads(path.read_bytes());x=np.asarray(row['coordinates'])
            digest=hashlib.sha256(x.tobytes()).hexdigest();reserved=json.loads(marker.read_bytes())
            if (path.parent.name!=f'eval-{index:06d}-{digest[:12]}'
                    or marker.name!=f'attempt-{index:06d}.json'
                    or reserved!=dict(schema_version=1,evaluation_index=index,coordinate_sha256=digest,
                        coordinates=x.tolist(),budget_sha256=_fingerprint(budget))
                    or row['contract_sha256']!=fingerprint or row['coordinate_sha256']!=digest):
                raise ValueError('native record/reservation identity differs')
            atoms=chart.to_atoms(x);atoms.calc=SinglePointCalculator(atoms,**study.audited_results(path.parent/'scf_000000'))
            validation=study.result_validator(index,atoms,path.parent)
            e,g=chart.energy_and_gradient(x,atoms)
            if (validation!=row['validation'] or e!=row['enthalpy_eV']
                    or not np.allclose(g,row['gradient'],atol=2e-12,rtol=0)
                    or not np.array_equal(atoms.get_stress(voigt=False),row['stress_eV_per_A3'])
                    or not worker.geometry_guard(chart,basis,basis.T@x)):
                raise ValueError('native EFS/gradient/geometry replay differs')
            records.append(row);sources.append(dict(arm=arm,result_sha256=study.sha256(path),validation=validation))
        outcome=json.loads((workdir/(arm+'_outcome.json')).read_bytes())
        if outcome!=summary['arms'][len(rows)] or outcome['arm']!=arm:raise ValueError('arm outcome differs')
        row=arm_diagnostics(basis,np.asarray(config['internal_axes_36_by_2']),energy,gradient,outcome,records)
        row['optimizer_replay']=replay_arm(chart,basis,config,energy,gradient,matrices[0],outcome,records)
        row['callback_and_center_hit_counts_worker_declared_not_trace_replayed']=False
        rows.append(row)
    if summary['new_probe_attempts']!=sum(r['completed_SCFs'] for r in rows):raise ValueError('total attempt count differs')
    if sum(r['sum_audited_call_seconds'] for r in rows)>int(fields[3])+1:
        raise ValueError('sequential native call time exceeds allocation time')
    both=all(r['final_candidate'] is not None for r in rows)
    comparison=None
    if both:
        a,b=[r['final_candidate'] for r in rows]
        comparison=dict(energy_difference_seed_minus_identity_eV_cell=b['energy_eV_cell']-a['energy_eV_cell'],
            internal_coordinate_difference_A=float(np.linalg.norm(np.asarray(b['coordinates_internal_A'])-a['coordinates_internal_A'])),
            common_gradient_target_reached=True,same_branch_or_material_speedup_certified=False)
    return dict(status='native_complete_conditional_canary_audit_not_PES_TS_or_prediction',channel=channel,
        job_id=job,sacct=accounting,allocation_core_hours=int(fields[3])*32/3600,
        arms=rows,comparison=comparison,raw_SCFs_reparsed=len(sources),native_record_sources=sources,
        center_source=manifest['centers'][channel]['raw_source'],native_audit_sha256=PINS[channel],
        configuration_sha256=CONFIGS[channel],runtime_manifest_sha256=RUNTIME_SHA,
        original_148_curvature_SCFs_per_center_and_prior_G2_costs_additional=True,
        physical_reserved_geometry_generated_or_read=False,new_DFT_or_scheduler_writes=0,
        conditional_stability_TS_full_face_or_material_speedup_certified=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--channel',choices=tuple(PINS),required=True)
    p.add_argument('--curvature-report',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise FileExistsError('fresh external audit required')
    if any(args.output.resolve().is_relative_to(path.resolve()) for path in (study.ROOT.parent,study.ASSETS)):
        raise ValueError('audit must not write inside the study/source/assets')
    result=audit(args.channel,args.curvature_report)
    args.output.parent.mkdir(parents=True,exist_ok=True);study.exclusive_json(args.output,result)
    print(json.dumps(dict(status=result['status'],channel=args.channel,raw_SCFs_reparsed=result['raw_SCFs_reparsed'],
        arms=[dict(arm=r['arm'],SCFs=r['completed_SCFs'],candidate=r['final_candidate'] is not None) for r in result['arms']])))


if __name__=='__main__':main()
