"""Retrospective raw-quadratic agreement on already audited training nodes.

No fitted coefficients, new DFT or independent test labels. A numerical
zero floor here defines only a formal quadratic, not physical R stability.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT))
from scripts import hfo2_G3_conditional_faces as faces
from vcneb.quadratic_reduction import condition_quadratic_energy

PINS={'T':'1a7b3530bb8d389196b9ad747df4f642e405ca1cdbb2d8da38d0ad529fefeda6',
      'preserving':'3c2c90ab703c2e1cc7c9c699daee5e0ded8f0d54d1b5d038c33a862e1b89b8ee'}


def evaluate_points(model,points,center_energy):
    """Prediction minus saved DFT-candidate value; no last-trial promotion."""
    zero=next(p for p in points if p['id']==0)
    result=[]
    for point in points:
        q=np.asarray(point['q_A']);observed=point['energy_eV_cell']-center_energy
        change=model.energy_change(q)
        result.append(dict(id=point['id'],q_A=point['q_A'],
            model_energy_change_eV_cell=change,observed_energy_change_eV_cell=observed,
            energy_error_meV_cell=(change-observed)*1000,
            zero_aligned_shape_error_meV_cell=(change-model.energy_change(np.zeros(2))-
                (point['energy_eV_cell']-zero['energy_eV_cell']))*1000,
            model_retained_gradient_eV_A=model.retained_gradient(q).tolist(),
            observed_retained_partial_gradient_eV_A=point['retained_gradient_eV_A'],
            retained_partial_gradient_error_norm_eV_A=float(np.linalg.norm(
                model.retained_gradient(q)-point['retained_gradient_eV_A']))))
    return result


def statistics(rows,ids):
    selected=[r for r in rows if r['id'] in ids]
    if len(selected)!=len(ids):raise ValueError('complete declared common numeric-Q cohort required')
    return dict(n=len(selected),ids=sorted(ids),
        max_absolute_energy_error_meV_cell=max(abs(r['energy_error_meV_cell']) for r in selected),
        max_absolute_zero_aligned_shape_error_meV_cell=max(abs(r['zero_aligned_shape_error_meV_cell']) for r in selected),
        RMS_zero_aligned_shape_error_meV_cell=float(np.sqrt(np.mean([
            r['zero_aligned_shape_error_meV_cell']**2 for r in selected]))),
        max_retained_partial_gradient_error_norm_eV_A=max(r['retained_partial_gradient_error_norm_eV_A'] for r in selected))


def axial_differences(points,h=.1):
    """Finite candidate-energy differences and raw partial-gradient secants."""
    zero=next(p for p in points if p['id']==0)
    energies=[];columns=[]
    for axis in range(2):
        q=np.eye(2)[axis]*h
        plus=next(p for p in points if np.array_equal(p['q_A'],q))
        minus=next(p for p in points if np.array_equal(p['q_A'],-q))
        energies.append((plus['energy_eV_cell']+minus['energy_eV_cell']-2*zero['energy_eV_cell'])/h**2)
        columns.append((np.asarray(plus['retained_gradient_eV_A'])-minus['retained_gradient_eV_A'])/(2*h))
    secants=np.column_stack(columns)
    return dict(step_A=h,candidate_energy_axial_second_differences_eV_A2=energies,
        raw_retained_partial_gradient_secants_eV_A2=secants.tolist(),
        raw_secant_reciprocity_Frobenius_defect_eV_A2=float(np.linalg.norm(secants-secants.T)),
        infinitesimal_conditional_Hessian_or_physical_error_bound=False)


def analyze():
    data={};points={};reports={}
    for channel,pin in PINS.items():
        path=faces.BASE/'G3_face_prefix_audit_E163_20261015'/f'{channel}_native_prefix.json'
        raw=path.read_bytes()
        if hashlib.sha256(raw).hexdigest()!=pin:raise ValueError('audited native prefix bytes differ')
        reports[channel]=json.loads(raw);d=faces.load_inputs(channel);data[channel]=d
        report=reports[channel]
        if (report['audit_mode']!='prefix' or report['new_DFT_or_scheduler_writes']!=0
                or report['physical_reserved_geometry_generated_or_read'] is not False
                or report['gradient_candidates']!=report['completed_nodes_audited']
                or any(r['candidate'] is None or r['replay']!=dict(
                    original_order_cache_only_replay_passed=True,new_DFT=0) for r in report['candidate_rows'])):
            raise ValueError('complete audited gradient-candidate prefix required')
        b,q=d['basis'],d['axes'];origin=d['origin_point']
        points[channel]=[dict(id=0,q_A=[0.,0.],energy_eV_cell=origin['energy_eV_cell'],
            retained_gradient_eV_A=(q.T@(b.T@d['origin_record']['gradient'])).tolist())]
        points[channel]+=[dict(id=r['id'],q_A=r['q_A'],energy_eV_cell=r['candidate']['energy_eV_cell'],
            retained_gradient_eV_A=r['candidate']['retained_gradient_components_eV_A']) for r in report['candidate_rows']]
        if [p['id'] for p in points[channel]]!=list(range(report['completed_nodes_audited'])):
            raise ValueError('original consecutive native point identities required')
    common=set.intersection(*[{p['id'] for p in ps} for ps in points.values()])
    for i in common:
        if points['T'][i]['q_A']!=points['preserving'][i]['q_A']:
            raise ValueError('common numeric-Q policy differs; do not compare unmatched cohorts')
    rows=[]
    for channel,d in data.items():
        b,q=d['basis'],d['axes'];r=np.linalg.svd(q.T,full_matrices=True)[2][2:].T
        for block in d['curvature']['blocks']:
            raw=np.asarray(block['raw_projected_eV_A2'])[:36,:36];h=(raw+raw.T)/2
            # Algebraic gate only: no measured physical stability floor is available.
            model=condition_quadratic_energy(h,b.T@d['center_gradient'],d['center_energy'],q,r,stability_floor=0.)
            values=evaluate_points(model,points[channel],d['center_energy'])
            observed_zero=points[channel][0]['energy_eV_cell']-d['center_energy']
            rows.append(dict(channel=channel,curvature_probe_step_A=block['step_A'],
                raw_formal_quadratic_origin_release_meV_cell=model.reference_energy_change*1000,
                observed_DFT_candidate_origin_release_meV_cell=observed_zero*1000,
                retained_released_coupling_Frobenius_norm_eV_A2=float(np.linalg.norm(q.T@h@r)),
                formal_Schur_softening_Frobenius_norm_eV_A2=float(np.linalg.norm(model.curvature.softening)),
                model_retained_curvature_eV_A2=model.curvature.relaxed.tolist(),
                points=values,common_numeric_Q_statistics=statistics(values,common),
                all_available_point_statistics=statistics(values,{p['id'] for p in points[channel]})))
    return dict(status='retrospective_training_quadratic_disagreement_not_independent_prediction',rows=rows,
        common_numeric_Q_ids=sorted(common),available_nodes={k:len(p) for k,p in points.items()},
        common_numeric_Q_does_not_mean_same_physical_configurations=True,
        candidate_axial_differences={k:axial_differences(p) for k,p in points.items()},
        native_prefix_sha256=PINS,original_curvature_and_canary_pins=faces.PINS,
        all_model_parameters_from_original_center_no_coefficient_fit=True,
        zero_alignment_uses_seen_origin_and_is_only_a_shape_diagnostic=True,
        algebraic_stability_floor_zero_not_physical_resolution=True,
        new_DFT_or_scheduler_operations=0,new_raw_SCF_reparse=False,
        conditional_stability_full_faces_TS_or_material_prediction_gain_certified=False,
        physical_reserved_geometry_generated_or_read=False,
        limitations=[
            'All observations are seen training candidates, not independent validation.',
            'First-step axes are Hessian eigenvectors; vanishing Q-R coupling is a construction, not a method gain.',
            'A released-gradient norm pass does not certify stable minima or a smooth unique branch.',
            'Retained gradients are physical partials at finite-residual candidates, not certified envelope derivatives.',
            'Candidate axial differences at0.1Angstrom are not infinitesimal Hessians or stationary-index tests.',
            'Model disagreement can combine nonlinearity, finite-step effects, branch changes and residual relaxation.',
            'Two probe steps and training errors do not provide physical uncertainty bounds.',
            'No assertion that these pointwise energy discrepancies are whole-path barrier errors.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('fresh retrospective analysis output required')
    result=analyze();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(status=result['status'],common_Q_nodes=len(result['common_numeric_Q_ids']),new_DFT=0)))


if __name__=='__main__':main()
