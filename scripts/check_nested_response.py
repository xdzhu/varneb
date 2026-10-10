"""Bounded analytic checks of nested controls; no material forecast or DFT."""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np

from vcneb import (
    MeasuredStationaryCentre, StationaryResponseContract,
    nested_stationary_controls, nested_stationary_gap_responses,
)


def synthetic_centre(index, rotated):
    h=np.array([[2. if index==0 else -2.,.2,.3,.4,0.],
                [.2,4.,.1,.8,0.],[.3,.1,3.,1.1,0.],
                [.4,.8,1.1,6.,0.],[0.,0.,0.,0.,-7.]])
    g=np.array([.1,-.2,.05,.06,2.])
    u=np.linalg.qr(np.random.default_rng(192+index).normal(size=(5,5)))[0] if rotated else np.eye(5)
    basis=u.T
    known_h=u.T@h@u
    known_g=u.T@g
    centre=MeasuredStationaryCentre(
        measured_basis=basis[:,:4],hessian_columns=known_h@basis[:,:4],gradient=known_g,
        reference_potential_eV_cell=-10000.+.8*index,
        retained_internal_basis=basis[:,:1],admissible_internal_basis=basis[:,:3],
        atomic_internal_basis=basis[:,:2],control_direction=basis[:,3],expected_index=index,
        stability_floor_eV_A2=.01,internal_curvature_floor_eV_A2=.01,
        symmetry_tolerance_eV_A2=1e-8,source_audit_sha256=('a' if index==0 else 'b')*64,
        measurement_ids=(f'analytic-centre-{index}',))
    return centre,known_h,known_g


def direct_energy(centre,h,g,branch,t):
    internal=np.column_stack((branch.conditional.retained_basis[:,:-1],branch.conditional.eliminated_basis))
    c=centre.control_direction
    y=internal@(-np.linalg.solve(internal.T@h@internal,internal.T@(g+h@c*t)))+c*t
    return float(g@y+.5*y@h@y),y


def run_benchmark():
    contract=StationaryResponseContract({'synthetic_Hamiltonian':'a'*64},'analytic_clamped','b'*64,
        'synthetic_strain','1',0.,'synthetic_shared_zero',4,0.,'E')
    errors=dict(signed_gap=0.,gap_response=0.,first_derivative=0.,second_derivative=0.,
                full_displacement=0.,stable_B4_B5_equivalence=0.)
    rows=[]
    for rotated in (False,True):
        for distinct_scales in (False,True):
            ai,hi,gi=synthetic_centre(0,rotated)
            ats,hs,gs=synthetic_centre(1,rotated)
            domains=[dict(scale_A_per_parameter_unit=3.+(2.*j if distinct_scales else 0.),
                          parameter_interval=(-.03,.03),maximum_full_displacement_A=2.) for j in (0,1)]

            def responses(x):
                return nested_stationary_gap_responses(ai,ats,x,contract=contract,
                    initial_domain=domains[0],bottleneck_domain=domains[1],reference_gap_eV_cell=.8,
                    reference_gap_identity_tolerance_eV_cell=1e-8)

            for x in (-.012,0.,.019):
                current=responses(x)
                lower,upper=responses(x-1e-5),responses(x+1e-5)
                for result,lo,up in zip(current,lower,upper):
                    point=result.response
                    if point is None or lo.response is None or up.response is None:
                        raise ValueError('synthetic resolved quadratic unexpectedly unavailable')
                    ei,yi=direct_energy(ai,hi,gi,result.initial.branch,domains[0]['scale_A_per_parameter_unit']*x)
                    es,ys=direct_energy(ats,hs,gs,result.bottleneck.branch,domains[1]['scale_A_per_parameter_unit']*x)
                    ei0,_=direct_energy(ai,hi,gi,result.initial.branch,0.)
                    es0,_=direct_energy(ats,hs,gs,result.bottleneck.branch,0.)
                    errors['signed_gap']=max(errors['signed_gap'],abs(point.signed_stationary_gap_eV_fu-(.8+es-ei)/4))
                    errors['gap_response']=max(errors['gap_response'],abs(point.gap_response_eV_fu-((es-es0)-(ei-ei0))/4))
                    errors['first_derivative']=max(errors['first_derivative'],abs(
                        (up.response.signed_stationary_gap_eV_fu-lo.response.signed_stationary_gap_eV_fu)/2e-5-
                        point.gap_derivative_eV_fu_per_parameter_unit))
                    errors['second_derivative']=max(errors['second_derivative'],abs(
                        (up.response.gap_derivative_eV_fu_per_parameter_unit-lo.response.gap_derivative_eV_fu_per_parameter_unit)/2e-5-
                        point.gap_curvature_eV_fu_per_parameter_unit2))
                    for centre,h,g,r,y,t in ((ai,hi,gi,result.initial,yi,domains[0]['scale_A_per_parameter_unit']*x),
                                            (ats,hs,gs,result.bottleneck,ys,domains[1]['scale_A_per_parameter_unit']*x)):
                        errors['full_displacement']=max(errors['full_displacement'],float(np.max(np.abs(r.branch.evaluate(t).full_displacement-y))))
                    rows.append(dict(rotated=rotated,distinct_scales=distinct_scales,parameter=x,control=result.control,
                        signed_gap_eV_fu=point.signed_stationary_gap_eV_fu,gap_response_eV_fu=point.gap_response_eV_fu,
                        initial_response_eV_cell=point.initial_response_eV_cell,bottleneck_response_eV_cell=point.bottleneck_response_eV_cell,
                        anchor_correction_eV_cell=point.stationary_anchor_correction_eV_cell))
                errors['stable_B4_B5_equivalence']=max(errors['stable_B4_B5_equivalence'],
                    abs(current[-2].response.signed_stationary_gap_eV_fu-current[-1].response.signed_stationary_gap_eV_fu))

    c,h,g=synthetic_centre(1,False)
    partial=replace(c,measured_basis=c.measured_basis[:,[0,3]],hessian_columns=c.hessian_columns[:,[0,3]])
    partial_results=nested_stationary_controls(partial)
    if partial_results[0].branch is None or any(r.branch is not None for r in partial_results[1:]):
        raise ValueError('missing-curvature gate failed')
    # Instability belongs to an omitted cell coordinate, not an external control.
    h=h.copy();h[0,0]=2.;h[2,2]=-.6
    unstable=replace(c,hessian_columns=h@c.measured_basis)
    promoted=nested_stationary_controls(unstable)[-1]
    if promoted.branch is None or promoted.promoted_dimension!=1:
        raise ValueError('training-only instability promotion failed')
    energy,y=direct_energy(unstable,h,g,promoted.branch,.025)
    if abs(promoted.branch.evaluate(.025).energy_change-energy)>1e-12:
        raise ValueError('promoted branch disagrees with independent full solve')
    if max(errors.values())>1e-9:
        raise ValueError('independent analytic check failed; no material tuning')
    root=Path(__file__).resolve().parents[1]
    sources=['vcneb/nested_response.py','vcneb/response_cost.py','vcneb/stationary_branch.py',
             'vcneb/stationary_gap.py','vcneb/quadratic_reduction.py','scripts/check_nested_response.py']
    return dict(format_version=1,benchmark='analytic_nested_response_not_material_forecast',groups=4,
        paired_model_points=len(rows),maximum_errors=errors,rows=rows,
        partial_measurement_statuses=[r.status for r in partial_results],
        training_instability_promoted_dimensions=promoted.promoted_dimension,
        source_sha256={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sources},
        numpy_version=np.__version__,DFT_calls=0,calculator_calls=0,
        B2_to_B5_material_forecast_frozen=False,HfO2_advantage_proven=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('refusing an existing analytic benchmark artifact')
    result=run_benchmark()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8',newline='\n') as handle:
        json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:result[k] for k in ('benchmark','paired_model_points','maximum_errors','DFT_calls')}))


if __name__=='__main__':
    main()
