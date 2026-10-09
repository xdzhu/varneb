"""Bounded analytic paired-response check; not HfO2 material evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from vcneb.stationary_branch import stationary_quadratic_branch
from vcneb.stationary_gap import (
    ControlledStationaryModel, StationaryResponseContract, restricted_stationary_gap_response,
)


def run_benchmark():
    rows = []
    errors = {'signed_gap': 0., 'gap_response': 0., 'first_derivative': 0., 'second_derivative': 0.}
    context = StationaryResponseContract({'synthetic_Hamiltonian': 'a'*64}, 'analytic_clamped', 'b'*64,
                                        'strain', '1', 0., 'synthetic_shared_zero', 4, 0., 'E')
    for release in (False, True):
        for rotated in (False, True):
            for distinct_scales in (False, True):
                models = []
                for index in (0, 1):
                    h = np.array([[2. if index == 0 else -2., .2, .4, 0.],
                                  [.2, 4., .3, 0.], [.4, .3, 3., 0.], [0., 0., 0., -7.]])
                    g = np.array([.15, -.12, .09, 2.])
                    u = np.linalg.qr(np.random.default_rng(51+index).normal(size=(4,4)))[0] if rotated else np.eye(4)
                    basis = u.T
                    branch = stationary_quadratic_branch(u.T @ h @ u, u.T @ g, -10000.+.8*index,
                        basis[:,:1], basis[:,1:2] if release else np.empty((4,0)), basis[:,2],
                        admissible_internal_basis=basis[:,:2], expected_index=index,
                        stability_floor=.01, internal_curvature_floor=.01)
                    models.append(ControlledStationaryModel(branch, context, 3.+2.*index if distinct_scales else 3.,
                                                           (-.03,.03), 2., 'c'*64))

                def response(x):
                    return restricted_stationary_gap_response(*models, x, reference_gap_eV_cell=.8,
                                                             reference_gap_identity_tolerance_eV_cell=1e-8)

                def direct(model, x):
                    branch = model.branch
                    q = branch.conditional.retained_basis[:,:-1]
                    r = branch.conditional.eliminated_basis
                    internal = np.column_stack((q,r)); c = branch.conditional.retained_basis[:,-1]
                    t = model.scale_A_per_parameter_unit*x
                    y = internal @ (-np.linalg.solve(internal.T @ branch.hessian @ internal,
                                    internal.T @ (branch.gradient+branch.hessian @ c*t)))+c*t
                    return float(branch.gradient @ y+.5*y @ branch.hessian @ y)

                e0 = direct(models[1],0.)-direct(models[0],0.)
                for x in (-.012,0.,.019):
                    point = response(x)
                    gap = (.8+direct(models[1],x)-direct(models[0],x))/4
                    dx = 1e-5; lower, upper = response(x-dx), response(x+dx)
                    derivative = (upper.signed_stationary_gap_eV_fu-lower.signed_stationary_gap_eV_fu)/(2*dx)
                    curvature = (upper.gap_derivative_eV_fu_per_parameter_unit-lower.gap_derivative_eV_fu_per_parameter_unit)/(2*dx)
                    errors['signed_gap'] = max(errors['signed_gap'],abs(gap-point.signed_stationary_gap_eV_fu))
                    errors['gap_response'] = max(errors['gap_response'],abs(gap-(.8+e0)/4-point.gap_response_eV_fu))
                    errors['first_derivative'] = max(errors['first_derivative'],abs(derivative-point.gap_derivative_eV_fu_per_parameter_unit))
                    errors['second_derivative'] = max(errors['second_derivative'],abs(curvature-point.gap_curvature_eV_fu_per_parameter_unit2))
                    rows.append({'release':release,'rotated':rotated,'distinct_scales':distinct_scales,'parameter':x,
                                 'stationary_anchor_gap_eV_fu':point.stationary_anchor_gap_eV_fu,
                                 'signed_stationary_gap_eV_fu':point.signed_stationary_gap_eV_fu,
                                 'gap_response_eV_fu':point.gap_response_eV_fu,
                                 'initial_response_eV_cell':point.initial_response_eV_cell,
                                 'bottleneck_response_eV_cell':point.bottleneck_response_eV_cell,
                                 'first_derivative':point.gap_derivative_eV_fu_per_parameter_unit,
                                 'second_derivative':point.gap_curvature_eV_fu_per_parameter_unit2})
    if max(errors[k] for k in ('signed_gap','gap_response')) > 1e-12 or max(errors.values()) > 1e-9:
        raise ValueError('analytic paired-gap check failed; retain failure rather than tune on material labels')
    root = Path(__file__).resolve().parents[1]
    sources = ['vcneb/stationary_gap.py','vcneb/stationary_branch.py','vcneb/quadratic_reduction.py',
               'vcneb/relaxed_curvature.py','scripts/check_stationary_gap.py']
    return {'format_version':1,'benchmark':'analytic_restricted_stationary_gap_not_material_forecast',
            'groups':8,'points':len(rows),'maximum_errors':errors,'rows':rows,
            'source_sha256':{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sources},
            'numpy_version':np.__version__,'DFT_calls':0,'calculator_calls':0,
            'B2_to_B5_material_forecast_frozen':False,'HfO2_advantage_proven':False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('refusing an existing stationary-gap benchmark report')
    result = run_benchmark()
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'groups':result['groups'],'points':result['points'],'maximum_errors':result['maximum_errors']}))


if __name__ == '__main__':
    main()
