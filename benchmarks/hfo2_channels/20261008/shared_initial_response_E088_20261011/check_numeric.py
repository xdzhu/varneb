"""Bounded NumPy-only analytic replay; not material or unseen-label evidence."""

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from vcneb import (
    ControlledStationaryModel, SharedInitialEnergyResponse, StationaryResponseContract,
    stationary_gap_with_shared_initial_energy, stationary_quadratic_branch,
)
import vcneb.shared_initial_response
import vcneb.stationary_branch
import vcneb.stationary_gap


def make_model(index, release, rotate, scale, shift=0.):
    h = np.array([[2. if index == 0 else -2., .2, .4, 0.],
                  [.2, 4., .3, 0.], [.4, .3, 3., 0.], [0., 0., 0., -7.]])
    g = np.array([.15, -.12, .09, 2.])
    u = np.linalg.qr(np.random.default_rng(51+index).normal(size=(4, 4)))[0] if rotate else np.eye(4)
    b = u.T
    branch = stationary_quadratic_branch(u.T@h@u, u.T@g, -10000.+.8*index+shift,
        b[:, :1], b[:, 1:2] if release else np.empty((4, 0)), b[:, 2],
        admissible_internal_basis=b[:, :2], expected_index=index,
        stability_floor=.01, internal_curvature_floor=.01)
    contract = StationaryResponseContract({'INPUT': 'a'*64}, 'synthetic_substrate', 'c'*64,
        'biaxial_strain', '1', 0., 'synthetic_common_zero', 4, 0., 'E')
    return ControlledStationaryModel(branch, contract, scale, (-.03, .03), 2., 'd'*64)


def direct_energy(model, x):
    branch = model.branch
    b = np.column_stack((branch.conditional.retained_basis[:, :-1], branch.conditional.eliminated_basis))
    c = branch.conditional.retained_basis[:, -1]
    t = model.scale_A_per_parameter_unit*x
    y = b@(-np.linalg.solve(b.T@branch.hessian@b, b.T@(branch.gradient+branch.hessian@c*t)))+c*t
    return float(branch.gradient@y+.5*y@branch.hessian@y)


def run_check():
    rows = []
    for release in (False, True):
        for rotate in (False, True):
            initial = make_model(0, release, rotate, 3.)
            saddle = make_model(1, release, rotate, 5.)
            ei0 = direct_energy(initial, 0.)
            for x in (-.012, 0., .019):
                increment = direct_energy(initial, x)-ei0
                anchor = -10000.+ei0
                feature = SharedInitialEnergyResponse(initial.contract, x, anchor, anchor+increment,
                    increment, 1e-7, 'e'*64, 'f'*64, 7)
                result = stationary_gap_with_shared_initial_energy(feature, saddle,
                    reference_gap_eV_cell=.8-ei0, reference_gap_identity_tolerance_eV_cell=1e-7)
                truth = (.8+direct_energy(saddle, x)-direct_energy(initial, x))/4
                rows.append(dict(release=release, rotate=rotate, parameter=x,
                    signed_gap_eV_fu=result.signed_stationary_gap_eV_fu,
                    independent_gap_eV_fu=truth, absolute_residual_eV_fu=abs(result.signed_stationary_gap_eV_fu-truth),
                    shared_preparation_calls=feature.target_endpoint_DFT_calls))
    maximum = max(row['absolute_residual_eV_fu'] for row in rows)
    if maximum > 1e-14:
        raise AssertionError(f'analytic identity residual {maximum}')
    # Visible nonlinear input: not a manufactured derivative or an initial Hessian.
    altered = replace(feature, recorded_energy_change_eV_cell=increment+.004,
                      target_energy_eV_cell=anchor+increment+.004)
    nonlinear = stationary_gap_with_shared_initial_energy(altered, saddle,
        reference_gap_eV_cell=.8-ei0, reference_gap_identity_tolerance_eV_cell=1e-7)
    if abs(nonlinear.signed_stationary_gap_eV_fu-result.signed_stationary_gap_eV_fu+.001) > 1e-14:
        raise AssertionError('visible nonlinear initial increment not preserved')
    modules = (vcneb.shared_initial_response, vcneb.stationary_branch, vcneb.stationary_gap)
    return dict(kind='synthetic_implementation_replay_only', rows=rows, points=len(rows),
        maximum_identity_residual_eV_fu=maximum, nonlinear_initial_increment_preserved=True,
        DFT_calls=0, calculator_calls=0, physical_initial_Hessian_computations=0,
        known_synthetic_initial_Hessian_used_only_for_independent_truth=True,
        material_forecast_or_holdout_access=False, module_sha256={m.__name__: hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest() for m in modules},
        Python=sys.version, NumPy=np.__version__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh report filename required')
    result = run_check()
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ('points', 'maximum_identity_residual_eV_fu', 'DFT_calls', 'NumPy')}))


if __name__ == '__main__':
    main()
