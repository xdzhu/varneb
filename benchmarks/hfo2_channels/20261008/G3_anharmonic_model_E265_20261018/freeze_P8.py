"""Freeze same-target surrogate rollouts and retain old controls; zero DFT."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
BASE=ROOT/'benchmarks/hfo2_channels/20261008'
TRAINING=BASE/'G3_preserving_branches_E260_20261017'
FIT_SHA='9b301edb4d96646a324293f90b5a8fa6b85d8b4ee4ef6ca9eed791a511f7c935'
START_SHA='90e63ea901c15abc19cb1e552bde55f85221fc928032947be6b30b401bb8b2fc'
MESSAGE='mode evaluation attempt budget exhausted; calculator was not called'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value


def evaluate(model,x,polynomial):
    """Physical internal energy and dual gradient of the frozen ansatz."""
    dx=np.asarray(x)-np.asarray(model['origin_internal_A'])
    transform=np.asarray(model['chart_transform_36_by_36']);h=np.asarray(model['symmetric_measured_extension_eV_A2'])
    origin_g=np.asarray(model['origin_gradient_internal_eV_A']);beta=np.asarray(model['fit']['coefficients'])
    values,jacobian=polynomial.monomials((dx@transform/.1)[None,:],np.asarray(model['exponents']))
    energy=model['origin_energy_eV_cell']+origin_g@dx+.5*dx@h@dx+float(values[0]@beta)
    gradient=origin_g+h@dx+transform@(jacobian[0]@beta)/.1
    return float(energy),gradient


def freeze():
    fit_path=HERE/'training_fit.json';start_path=BASE/'G3_P_branch_independent_E263_20261018/starts_E264.json'
    if sha(fit_path)!=FIT_SHA or sha(start_path)!=START_SHA:raise ValueError('exact selected fits and starts required')
    fitted=json.loads(fit_path.read_bytes());starts=json.loads(start_path.read_bytes())
    source=module('E265_frozen_polynomial',HERE/'fit.py')
    if (hashlib.sha256((HERE/'fit.py').read_bytes().replace(b'\r\n',b'\n')).hexdigest()!=fitted['source_LF_sha256']
            or sha(HERE/'protocol.md')!=fitted['protocol_sha256']):raise ValueError('registered fit source/protocol changed')
    worker=module('E265_original_rollout_worker',TRAINING/'worker.py')
    config=json.loads((TRAINING/'config.json').read_bytes())
    data,_,_,free,seed,_=worker.validate(config,TRAINING)
    basis=data['basis'];direction=np.asarray(config['branch_direction_36']);models=[]
    for model in fitted['models']:
        if not model['fit']['identified']:raise ValueError('rank-deficient selected model cannot silently forecast')
        predictions=[]
        for point in starts['points']:
            cache={}
            def oracle(x):
                key=(basis@x).tobytes()
                if key not in cache:
                    if len(cache)==12:raise RuntimeError(MESSAGE)
                    cache[key]=evaluate(model,x,source)
                return cache[key]
            def valid(x):
                try:worker.geometry(data,basis@x,point['q_A']);return True
                except ValueError:return False
            result=None;reason=None
            try:
                result=worker.solver.relax_branch_candidate(oracle,np.asarray(point['starting_internal_coordinates_A']),
                    free,direction,initial_hessian=seed,valid=valid,gradient_tolerance=.03,
                    branch_gradient_tolerance=1e-4,initial_trust_radius=.05,max_iterations=200)
                if result['status']!='gradient_candidate_stability_unchecked':reason=result['status']
            except RuntimeError as exc:
                if str(exc)!=MESSAGE:raise
                reason='surrogate_unique_oracle_budget_exhausted'
            value=(result['energy']-data['center_energy']) if reason is None else None
            predictions.append(dict(id=point['id'],q_A=point['q_A'],training_parent_id=point['training_parent_id'],
                prediction_raw_eV_cell=value,prediction_origin_aligned_eV_cell=value,
                prediction_eV_fu=(None if value is None else value/4),abstention_reason=reason,
                model_unique_oracle_calls=len(cache),formal_model_outcome=result,independent_DFT_label=None,
                model_or_DFT_conditional_stability_certified=False))
        models.append(dict(model_id=f"{model['family']}@{model['step_A']:.2f}",
            family=model['family'],step_A=model['step_A'],parameters=model['fit']['parameters'],
            predictions=predictions,seen_origin_bias_eV_cell=0.,
            definition='same frozen twelve-call single-start candidate algorithm with surrogate E/g, not certified conditional minimization'))
    legacy=BASE/'G3_T_independent_E258_20261017'
    standard_path=legacy/'standard_forecast.json';reference_path=legacy/'reference_forecast.json'
    if (sha(standard_path)!='aa77cdd7460ad809a5fbe06e80c684f7f73ce118835067c380ec219b96abd2d8'
            or sha(reference_path)!='97692cb0989a4bafdff86349a732126cce5d7d93d9a7d47926ce844fb1f0b285'):
        raise ValueError('original numerical controls must remain exact')
    original=next(r for r in json.loads(standard_path.read_bytes())['channels'] if r['channel']=='preserving')
    references=[r for r in json.loads(reference_path.read_bytes())['rows'] if r['channel']=='preserving']
    if len(original['targets'])!=8 or len(references)!=14 or len(models)!=4:
        raise ValueError('complete target/model/control panel required')
    return dict(status='prospective_remaining_P8_candidate_forecasts_not_stable_sheet_or_physical_strain_freeze',
        source_sha256=dict(training_fit=FIT_SHA,start_plan=START_SHA,
            original_standard_forecast=sha(standard_path),original_reference_forecast=sha(reference_path),
            P8_protocol=sha(HERE/'P8_PROTOCOL.md')),
        source_LF_sha256={p.name:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                          for p in (Path(__file__),HERE/'fit.py')},
        new_models=models,original_standard_preserving_record=original,original_reference_preserving_records=references,
        original_conditional_2D_refusal_preserved=True,legacy_values_do_not_revalidate_old_target=True,
        strongest_full_T_optical_cell_harmonic_reexpression='same full-space harmonic reference, not an independent gain',
        candidate_target_stability_unverified=True,stable_conditional_model_admission=False,
        physical_input_sha256=config['physical_input_sha256'],binary_sha256=config['binary_sha256'],
        pressure_GPa=0.,physical_substrate_strain=0.,ordinary_NEB_target_eV_A=.10,
        domain=dict(exact_Q_targets_A=[p['q_A'] for p in starts['points']],
            original_Q_training_hull_L1_radius_A=.1,maximum_total_internal_radius_A=.5,
            radius_center='original internal chart zero, not released origin',minimum_distance_A=1.6,
            free_basis_36_by_34=free.tolist(),zeta_direction_36=direction.tolist(),
            numerical_candidate_target_only=True,stable_conditional_predictive_domain=[],
            failed_training_candidates_preserved=[0,1,6],certified_interpolation_across_failed_center=False),
        common_preparation_cost=dict(new_training_SCFs=81,training_allocation_core_hours=85.6,
            original_derivative_SCFs=296,post_release_derivative_SCFs=272,T_tests_already_seen_SCFs=18,
            other_paths_faces_probes_additional=True,total_method_cost_not_yet_complete=True,
            shared_preparations_not_multiplied_by_model_count=True),
        original_denominator=8,maximum_new_P8_attempts=96,new_DFT_or_scheduler_operations=0,
        independent_P8_labels_seen=False,T8_labels_already_seen=True,
        physical_B0_B5_barrier_panel_frozen=False,physical_plus005_or_private_arrays_accessed=False,
        actual_DFT_deployment_preflight_and_hash_pinning_still_required=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('fresh prospective forecast artifact required')
    result=freeze()
    with args.output.open('x',encoding='utf-8',newline='\n') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps([dict(model=m['model_id'],predicted=sum(p['prediction_raw_eV_cell'] is not None for p in m['predictions']),
        refusals=[dict(id=p['id'],reason=p['abstention_reason']) for p in m['predictions'] if p['abstention_reason']])
        for m in result['new_models']],indent=2))
