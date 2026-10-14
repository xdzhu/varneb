"""Bounded same-input center conditional canary, never an automatic dispatcher.

Prepare only after a complete two-amplitude native audit. One modeling node
per center, two separately capped optimizer arms, no face/holdout/TS claim.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import os
from pathlib import Path

import numpy as np

from scripts import hfo2_G3_local_curvature as study
from vcneb import CalculatorModeEvaluator, ModePlane, relax_orthogonal_at_q

PARENT_SHA = '5f047f54ca83cb2dceb8cb71687eb9206a85c04205cdf26431565bc8a9b98a5a'
SETTINGS = dict(gradient_tolerance_eV_A=.03, initial_trust_radius_A=.05,
    maximum_internal_displacement_A=.5, minimum_distance_A=1.6,
    maximum_attempts_per_arm=12, maximum_attempts_per_center=24,
    maximum_iterations=200, arms=['identity', 'curvature_seed'], q=[0., 0.],
    modeling_nodes_used=1, independent_test_points_used=0,
    optimizer='safeguarded_bfgs', pressure_GPa=0., strain=0., climb=False)
CODE = tuple(dict.fromkeys((*study.CODE, 'vcneb/mode_surface.py',
                           'scripts/hfo2_G3_conditional_canary.py')))


def pinned(path, expected):
    if study.sha256(path) != expected:
        raise ValueError('original source bytes differ')
    return json.loads(path.read_bytes())


def measured_pair(manifest, channel, report):
    """Raw two-step consistency screen, not an error/stationary certificate."""
    if (report.get('status') != 'complete_local_curvature_diagnostic'
            or report.get('channel') != channel
            or any(report.get(k) != 148 for k in ('completed_records', 'planned_records', 'raw_SCFs_reparsed'))
            or len(report.get('native_record_sources', [])) != 148
            or report.get('failures_present') != []
            or report.get('center_source') != manifest['centers'][channel]['raw_source']
            or report.get('physical_reserved_geometry_generated_or_read') is not False
            or report.get('new_DFT_or_scheduler_operations') != 0
            or report.get('physical_error_bound_or_TS_certificate') is not False):
        raise ValueError('complete original two-amplitude native audit required; no partial-plane selection')
    chart, _ = study.validate(manifest, channel)
    basis = np.asarray(manifest['centers'][channel]['directions'])[:, :36]
    gradient = np.asarray(report['raw_center_gradient_eV_A'], dtype=float)
    if gradient.shape != (40,) or not np.isfinite(gradient).all():
        raise ValueError('finite complete center gradient required')
    blocks = report['blocks']
    if [(b['step_A'], b['complete_pairs'], b['full_stencil_at_this_step']) for b in blocks] != [(.01,37,True),(.02,37,True)]:
        raise ValueError('both original37-column signed stencils required')
    matrices = []
    for block in blocks:
        raw = np.asarray(block['raw_projected_eV_A2'], dtype=float)
        action = np.asarray(block['full_hessian_action_eV_A2'], dtype=float)
        if (raw.shape != (37,37) or action.shape != (40,37)
                or not np.isfinite(raw).all() or not np.isfinite(action).all()
                or not np.allclose(np.asarray(manifest['centers'][channel]['directions']).T@action,
                                   raw, atol=2e-12, rtol=0)):
            raise ValueError('complete measured actions and projections differ')
        h = .5*raw[:36,:36]+.5*raw[:36,:36].T
        if not np.allclose(np.linalg.eigvalsh(h),block['fixed_strain_internal_eigenvalues_eV_A2'],atol=2e-12,rtol=0):
            raise ValueError('measured raw spectrum differs')
        matrices.append(h)
    return chart, basis, gradient, matrices


def prepare(manifest, channel, report, report_sha):
    _, _, _, matrices = measured_pair(manifest, channel, report)
    spectra = [np.linalg.eigh(h) for h in matrices]
    if any(np.count_nonzero(w<0) != 2 or np.any(w==0) for w,_ in spectra):
        raise ValueError('the two-coordinate canary policy requires two raw negative directions at both steps')
    axes = spectra[0][1][:,:2]
    # Freeze the supplied orientation in the manifest. Later validation checks
    # the subspace, not LAPACK eigenvector signs, and never rematches atoms.
    complement = np.linalg.svd(axes.T,full_matrices=True)[2][2:].T
    minima = [float(np.linalg.eigvalsh(complement.T@h@complement)[0]) for h in matrices]
    if min(minima) <= 0:
        raise ValueError('one measured complement remains nonpositive; no forced-positive seed')
    return dict(schema_version=1, purpose='single_origin_conditional_canary_not_face_or_TS',
        parent_manifest_sha256=PARENT_SHA, channel=channel, native_audit_sha256=report_sha,
        settings=deepcopy(SETTINGS), internal_axes_36_by_2=axes.tolist(),
        selection='first-step complete two-negative eigenspace; second-step positive-complement raw screen',
        raw_complement_minima_eV_A2=minima, raw_sign_screen_not_error_certificate=True,
        full_internal_index_or_TS_certified=False, final_surface_or_forecast_domain_locked=False,
        physical_reserved_geometry_generated_or_read=False,
        normalized_runtime_source_sha256={n:study.text_sha(study.ROOT/n) for n in CODE})


def validate_config(config, manifest, channel, report, report_sha):
    candidate = prepare(manifest,channel,report,report_sha)
    if set(config)!=set(candidate):
        raise ValueError('canary configuration fields differ')
    # Exact policy/settings/source identities; the prepared axes are immutable
    # data, not regenerated from hardware-dependent eigenvector signs.
    for key in candidate:
        if (key not in ('internal_axes_36_by_2','raw_complement_minima_eV_A2')
                and json.dumps(config.get(key),sort_keys=True,allow_nan=False)
                != json.dumps(candidate[key],sort_keys=True,allow_nan=False)):
            raise ValueError('canary policy/source/physical contract differs')
    axes = np.asarray(config['internal_axes_36_by_2'],dtype=float)
    expected = np.asarray(candidate['internal_axes_36_by_2'])
    if (axes.shape!=(36,2) or not np.isfinite(axes).all()
            or not np.allclose(axes.T@axes,np.eye(2),atol=2e-12,rtol=0)
            or not np.allclose(axes@axes.T,expected@expected.T,atol=1e-10,rtol=0)
            or not np.allclose(config['raw_complement_minima_eV_A2'],candidate['raw_complement_minima_eV_A2'],atol=2e-12,rtol=0)):
        raise ValueError('prepared axes or raw screen differs')
    return measured_pair(manifest,channel,report)


class InternalEvaluator:
    """Exact linear pullback and audited zero-coordinate reuse, not new physics."""
    def __init__(self, evaluator, basis, energy, gradient):
        if np.iscomplexobj(basis) or np.iscomplexobj(gradient):
            raise ValueError('real audited internal span and gradient required')
        self.evaluator, self.basis = evaluator, np.array(basis,copy=True)
        self.energy, self.gradient = energy, np.array(gradient,copy=True)
        self.source_fixed_coordinate_covectors=np.array(evaluator.fixed_coordinate_covectors,copy=True)
        if (self.basis.shape!=(40,36) or self.gradient.shape!=(40,)
                or not np.isfinite(self.energy) or not np.isfinite(self.gradient).all()
                or not np.isfinite(self.basis).all() or np.any(self.basis[-1])
                or not np.array_equal(self.source_fixed_coordinate_covectors,np.eye(40)[:,-1:])
                or not np.allclose(self.basis.T@self.basis,np.eye(36),atol=2e-12,rtol=0)
                or not np.allclose(evaluator.fixed_coordinate_covectors.T@self.basis,0.,atol=2e-12,rtol=0)):
            raise ValueError('audited internal span/center and exact fixed-control exclusion required')
        self.basis.flags.writeable=False
        self.gradient.flags.writeable=False
        self.source_fixed_coordinate_covectors.flags.writeable=False
        self.center_cache_hits = 0
        self.callback_requests = 0

    @property
    def fixed_coordinate_covectors(self):
        pulled=self.basis.T@self.source_fixed_coordinate_covectors
        # The full40D control is retained above. A constraint whose pullback is
        # identically zero is satisfied by every coordinate in this embedding;
        # it is not a zero-length constraint in the solver's36D space.
        return pulled[:,np.any(pulled!=0.,axis=0)]

    def __call__(self, coordinates):
        if np.iscomplexobj(coordinates):
            raise ValueError('real internal coordinates required')
        values = np.asarray(coordinates,dtype=float)
        if values.shape!=(36,) or not np.isfinite(values).all():
            raise ValueError('finite declared internal coordinates required')
        full = self.basis@values
        self.callback_requests += 1
        if not np.any(values):
            self.center_cache_hits += 1
            return self.energy, self.basis.T@self.gradient
        energy, gradient = self.evaluator(full)
        return energy, self.basis.T@gradient


def geometry_guard(chart,basis,coordinates):
    if not np.isfinite(coordinates).all() or np.linalg.norm(coordinates)>SETTINGS['maximum_internal_displacement_A']:
        return False
    try:
        atoms = chart.to_atoms(basis@coordinates)
        distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances,np.inf)
        return bool(np.min(distances)>=SETTINGS['minimum_distance_A'])
    except (ValueError,np.linalg.LinAlgError):
        return False


def native_canary_preflight(config,manifest,channel,report,report_sha):
    _,basis,old_gradient,matrices=validate_config(config,manifest,channel,report,report_sha)
    chart,_,energy,gradient=study.native_preflight(manifest,channel)
    if not np.allclose(gradient,old_gradient,atol=2e-12,rtol=0):
        raise ValueError('audited native center gradient differs beyond arithmetic replay')
    if not geometry_guard(chart,basis,np.zeros(36)):
        raise ValueError('audited center violates declared geometry/domain guard')
    return chart,basis,energy,gradient,matrices


def run(config,manifest,channel,report,report_sha,workdir):
    if (os.environ.get('RUN_DFT')!='1' or os.environ.get('SLURM_JOB_PARTITION')!='hfacnormal01'
            or os.environ.get('SLURM_JOB_NUM_NODES')!='1' or os.environ.get('SLURM_NTASKS')!='1'
            or os.environ.get('SLURM_CPUS_PER_TASK')!='32'
            or not os.environ.get('SLURM_JOB_ID')
            or any(os.environ.get(n)!='1' for n in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'))):
        raise ValueError('authorized one-node hfacnormal01/32MPI allocation required')
    # Source is in RUN_ROOT/source; calculations only in its fresh sibling.
    if (workdir.resolve()!=study.ROOT.parent/channel or workdir.exists()
            or not workdir.resolve().is_relative_to(study.RUNS)):
        raise ValueError('fresh channel output beside the pinned source required')
    chart,basis,energy,gradient,matrices=native_canary_preflight(config,manifest,channel,report,report_sha)
    plane = ModePlane(np.zeros(36),np.asarray(config['internal_axes_36_by_2']),np.ones(36),np.eye(2),
        ('coupled internal1','coupled internal2'),'local length','G3-center-canary-only')
    workdir.mkdir()
    center=dict(source=manifest['centers'][channel]['raw_source'],energy_eV_cell=energy,
        full_gradient_eV_A=gradient.tolist(),new_center_SCFs=0,
        prescribed_control_excluded_exactly_by_internal_embedding=bool(not np.any(basis[-1])))
    study.exclusive_json(workdir/'center.json',center)
    results=[]
    for arm in SETTINGS['arms']:
        evaluator = CalculatorModeEvaluator(chart,
            study.make_factory(parameters={'source_directory':str(study.ASSETS)},command=study.COMMAND),
            workdir=workdir/arm,calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
            minimum_distance_A=SETTINGS['minimum_distance_A'],result_validator=study.result_validator,
            validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',
            max_evaluation_attempts=SETTINGS['maximum_attempts_per_arm'])
        evaluate=InternalEvaluator(evaluator,basis,energy,gradient)
        row=dict(arm=arm,status='not_converged',point=None)
        screened=dict(checks=0,rejected=0)
        def guard(x):
            screened['checks']+=1
            ok=geometry_guard(chart,basis,x)
            screened['rejected']+=int(not ok)
            return ok
        try:
            point=relax_orthogonal_at_q(plane,[0.,0.],evaluate,optimizer='safeguarded_bfgs',
                trial_validator=guard,
                initial_trust_radius=SETTINGS['initial_trust_radius_A'],
                initial_hessian=matrices[0] if arm=='curvature_seed' else None,
                gradient_tolerance=SETTINGS['gradient_tolerance_eV_A'],
                max_iterations=SETTINGS['maximum_iterations'])
            row.update(status='gradient_converged_candidate_curvature_unchecked',point=dict(
                coordinates_internal_A=point.coordinates.tolist(),coordinates_full_A=(basis@point.coordinates).tolist(),
                energy_eV_cell=point.energy,orthogonal_gradient_norm_eV_A=point.orthogonal_gradient_norm,
                callback_evaluations=point.total_n_evaluations))
        except RuntimeError as exc:
            failures=list((workdir/arm).glob('eval-*/failure.json'))
            known = str(exc) in ('mode evaluation attempt budget exhausted; calculator was not called',
                                 'no starting point reached the orthogonal-gradient tolerance')
            exhausted = str(exc).startswith('mode evaluation attempt budget exhausted')
            if (failures or not known
                    or (exhausted and evaluator.n_evaluation_attempts!=SETTINGS['maximum_attempts_per_arm'])):
                raise
            row['reason']=str(exc)
        records=[json.loads(p.read_text()) for p in sorted((workdir/arm).glob('eval-*/result.json'))]
        if (len(records)!=evaluator.n_new_evaluations
                or evaluator.n_evaluation_attempts>SETTINGS['maximum_attempts_per_arm']):
            raise ValueError('native completion/cumulative-attempt accounting differs')
        completed_seconds=sum(r['validation']['elapsed_seconds'] for r in records)
        row.update(center_cache_hits=evaluate.center_cache_hits,callback_requests=evaluate.callback_requests,
            sum_completed_native_call_seconds=completed_seconds,completed_call_proxy_core_hours=completed_seconds*32/3600,
            allocation_core_hours=None,failed_or_interrupted_attempts_not_refunded=True,
            trial_geometry_and_domain_screens=screened,
            reserved_attempts=evaluator.n_evaluation_attempts,new_completed_SCFs=evaluator.n_new_evaluations)
        study.exclusive_json(workdir/(arm+'_outcome.json'),row)
        results.append(row)
    result=dict(status='bounded_conditional_canary_not_PES_TS_or_prediction',channel=channel,
        job_id=os.environ['SLURM_JOB_ID'],settings=SETTINGS,arms=results,
        new_probe_attempts=sum(r['reserved_attempts'] for r in results),new_center_SCFs=0,
        all_arms_gradient_converged=all(r['point'] is not None for r in results),
        same_branch_or_material_speedup_certified=False,
        parent_curvature_and_center_preparation_costs_additional=True,
        physical_reserved_geometry_generated_or_read=False,
        physical_release_stability_or_TS_certified=False,full_surface_or_forecast_completed=False)
    study.exclusive_json(workdir/'canary.json',result)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('prepare','preflight','run'))
    p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--native-audit',type=Path,required=True)
    p.add_argument('--native-audit-sha256',required=True)
    p.add_argument('--channel',choices=tuple(study.SOURCES),required=True)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--config-sha256')
    p.add_argument('--workdir',type=Path)
    args=p.parse_args()
    if args.mode=='prepare' and args.config.exists():
        raise FileExistsError('fresh configuration required')
    manifest=pinned(args.parent,PARENT_SHA)
    report=pinned(args.native_audit,args.native_audit_sha256)
    if args.mode=='prepare':
        config=prepare(manifest,args.channel,report,args.native_audit_sha256)
        args.config.parent.mkdir(parents=True,exist_ok=True);study.exclusive_json(args.config,config)
        print(json.dumps(dict(status='prepared_single_canary_not_forecast',channel=args.channel,new_DFT=0)))
    elif args.mode=='preflight':
        if not args.config_sha256:
            raise ValueError('exact configuration SHA required')
        config=pinned(args.config,args.config_sha256)
        chart,basis,energy,gradient,_=native_canary_preflight(config,manifest,args.channel,report,args.native_audit_sha256)
        axes=np.asarray(config['internal_axes_36_by_2']);g=basis.T@gradient
        print(json.dumps(dict(status='native_canary_preflight_passed_no_DFT',channel=args.channel,
            cached_center_energy_eV_cell=energy,starting_released_gradient_norm_eV_A=float(np.linalg.norm(g-axes@(axes.T@g))),
            checked_source_files=len(CODE),new_DFT=0,physical_reserved_geometry_generated_or_read=False,
            physical_release_stability_or_TS_certified=False)))
    else:
        if not args.config_sha256 or args.workdir is None:
            raise ValueError('exact configuration SHA and fresh workdir required')
        result=run(pinned(args.config,args.config_sha256),manifest,args.channel,report,
                   args.native_audit_sha256,args.workdir)
        print(json.dumps(result))


if __name__=='__main__':main()
