"""Independent full-joint harmonic null and actual-chart cost lower bounds.

This is not a material Hessian, a G3 probe registration or a forecast. No
calculator is constructed, no geometry written and no scheduler called.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np

from scripts.check_nested_response import synthetic_centre
from scripts.relax_clamped_ase_endpoint import load_seed
from vcneb import (
    BiaxialClampedCurvatureCoordinates, StationaryResponseContract,
    nested_stationary_controls, nested_stationary_gap_responses,
)


ROOT = Path(__file__).resolve().parents[1]
VARIANTS = ('axis', 'tilted_atomic', 'omitted_negative', 'mixed_atomic_cell', 'complete_internal')
JOINT = ('B4_joint_release', 'B5_training_instability_promotion')


def partition(centre, name):
    """Change the retained subspace, not just the gauge within that space."""
    b = centre.admissible_internal_basis
    if name == 'complete_internal':
        q = b
    else:
        angle = {'axis': 0., 'tilted_atomic': .7, 'omitted_negative': 1.2,
                 'mixed_atomic_cell': .4}[name]
        v = b[:, 0] * np.cos(angle) + b[:, 1] * np.sin(angle)
        q = v[:, None]
        if name == 'mixed_atomic_cell':
            q = np.column_stack((v, b[:, 2]))
    return replace(centre, retained_internal_basis=q)


def independent_full_stationary(centre, h, g, t):
    """Solve in the original entire B, without any nested-result Q/R data."""
    b, c = centre.admissible_internal_basis, centre.control_direction
    k = b.T @ h @ b
    cross = b.T @ h @ c
    internal = -np.linalg.solve(k, b.T @ g + cross * t)
    y = b @ internal + c * t
    energy = float(g @ y + .5 * y @ h @ y)
    derivative = float(c @ (g + h @ y))
    curvature = float(c @ h @ c - cross @ np.linalg.solve(k, cross))
    return y, energy, derivative, curvature


def check_partitions():
    errors = dict(full_displacement_A=0., energy_eV_cell=0.,
                  control_gradient_eV_A=0., control_curvature_eV_A2=0.,
                  paired_signed_gap_eV_fu=0., paired_response_eV_fu=0.,
                  paired_derivative_eV_fu=0., paired_curvature_eV_fu=0.)
    rows, pair_rows, refusals = [], [], []
    frozen_difference = 0.
    contract = StationaryResponseContract(
        {'analytic_Hamiltonian': 'a' * 64}, 'synthetic_clamped', 'b' * 64,
        'synthetic_strain', '1', 0., 'shared_analytic_zero', 4, 0., 'E',
    )
    domains = [dict(scale_A_per_parameter_unit=s, parameter_interval=(-.03, .03),
                    maximum_full_displacement_A=3.) for s in (3., 5.)]
    for rotated in (False, True):
        original = [synthetic_centre(i, rotated) for i in (0, 1)]
        for name in VARIANTS:
            centres = [partition(c, name) for c, _, _ in original]
            for index, ((_, h, g), centre) in enumerate(zip(original, centres)):
                controls = {r.control: r for r in nested_stationary_controls(centre)}
                for label in JOINT:
                    item = controls[label]
                    if item.branch is None:
                        refusals.append(dict(rotated=rotated, index=index, partition=name,
                                             control=label, reason=item.reason))
                        if not (index == 1 and name == 'omitted_negative' and label == JOINT[0]):
                            raise AssertionError('unexpected unavailable fully measured control')
                        continue
                    for t in (-.04, 0., .04):
                        y, e, d, k = independent_full_stationary(centre, h, g, t)
                        point = item.branch.evaluate(t)
                        residuals = dict(
                            full_displacement_A=float(np.max(np.abs(point.full_displacement-y))),
                            energy_eV_cell=abs(point.energy_change-e),
                            control_gradient_eV_A=abs(point.control_gradient-d),
                            control_curvature_eV_A2=abs(item.branch.control_curvature-k),
                        )
                        for field, value in residuals.items():
                            errors[field] = max(errors[field], value)
                        rows.append(dict(rotated=rotated, index=index, partition=name,
                                         retained_dimension=centre.retained_internal_basis.shape[1],
                                         control=label, promoted_dimension=item.promoted_dimension,
                                         parameter_A=t, errors=residuals))
                if name == 'axis':
                    f = controls['B2_frozen'].branch
                    _, e, _, _ = independent_full_stationary(centre, h, g, .04)
                    frozen_difference = max(frozen_difference, abs(f.evaluate(.04).energy_change-e))
            for parameter in (-.01, 0., .01):
                results = nested_stationary_gap_responses(
                    *centres, parameter, contract=contract, initial_domain=domains[0],
                    bottleneck_domain=domains[1], reference_gap_eV_cell=.8,
                    reference_gap_identity_tolerance_eV_cell=1e-8,
                )
                direct = [independent_full_stationary(c, h, g, s * parameter)
                          for (c, h, g), s in zip(original, (3., 5.))]
                anchors = [independent_full_stationary(c, h, g, 0.)
                           for c, h, g in original]
                target = dict(
                    paired_signed_gap_eV_fu=(.8 + direct[1][1] - direct[0][1]) / 4,
                    paired_response_eV_fu=((direct[1][1]-anchors[1][1])-
                                          (direct[0][1]-anchors[0][1])) / 4,
                    paired_derivative_eV_fu=(5.*direct[1][2]-3.*direct[0][2]) / 4,
                    paired_curvature_eV_fu=(25.*direct[1][3]-9.*direct[0][3]) / 4,
                )
                for item in results:
                    if item.control not in JOINT or item.response is None:
                        continue
                    point = item.response
                    actual = dict(
                        paired_signed_gap_eV_fu=point.signed_stationary_gap_eV_fu,
                        paired_response_eV_fu=point.gap_response_eV_fu,
                        paired_derivative_eV_fu=point.gap_derivative_eV_fu_per_parameter_unit,
                        paired_curvature_eV_fu=point.gap_curvature_eV_fu_per_parameter_unit2,
                    )
                    residuals = {field: abs(actual[field]-value) for field, value in target.items()}
                    for field, value in residuals.items():
                        errors[field] = max(errors[field], value)
                    pair_rows.append(dict(rotated=rotated, partition=name, control=item.control,
                                          parameter=parameter, errors=residuals))
    if max(errors.values()) > 1e-9 or frozen_difference < 1e-5 or len(refusals) != 2:
        raise AssertionError('partition null, omission contrast or refusal boundary failed')
    return dict(maximum_errors=errors, stationary_points=len(rows), paired_points=len(pair_rows),
                rows=rows, paired_rows=pair_rows, correct_joint_release_refusals=refusals,
                frozen_vs_joint_max_energy_difference_eV_cell=frozen_difference)


def chart_cost_bounds():
    """Count full central-paired columns in existing starters, not TS centres.

    One or two step lengths mean 1+2*m or 1+4*m geometric evaluations when
    no valid centre/probe cache exists. These are design counts, not actual
    calls, certified minimum costs of all possible algorithms or permissions.
    The 81 surface-node cap is not silently reassigned as Hessian allowance.
    """
    base = ROOT/'benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds'
    manifest = json.loads((base/'clamped_seed_manifest.json').read_text())
    rows = []
    for item in manifest['seeds']:
        atoms, boundary, record = load_seed(base/item['manifest_file'])
        chart = BiaxialClampedCurvatureCoordinates(
            atoms, boundary=boundary, cell_scale_A=record['cell_scale_A'],
            anchor_strain=record['strain'], strain_scale_A=record['cell_scale_A'],
        )
        b = chart.internal_relaxation_basis()
        a = np.eye(chart.size)[:, :3*len(atoms)]
        atomic_rank = int(np.linalg.matrix_rank(b.T@a, tol=1e-10))
        m = np.column_stack((b, chart.controlled_direction()))
        count = int(np.linalg.matrix_rank(m, tol=1e-10))
        rows.append(dict(phase=record['phase_label'], strain=record['strain'],
                         starter_sha256=record['seed_sha256'], joint_chart_dimension=chart.size,
                         translation_free_internal_dimension=b.shape[1], atomic_dimension=atomic_rank,
                         open_cell_dimension=b.shape[1]-atomic_rank,
                         full_internal_plus_control_measured_columns=count,
                         full_central_pair_one_step_uncached_evaluations=1+2*count,
                         full_central_pair_two_steps_uncached_evaluations=1+4*count,
                         actual_new_DFT_calls=0, actual_bottleneck_or_TS=False))
    if len(rows) != 10 or any((r['translation_free_internal_dimension'],
                             r['full_internal_plus_control_measured_columns']) != (36,37) for r in rows):
        raise AssertionError('registered Hf4O8 clamped chart dimensions changed')
    return rows


def run_check():
    sources = ['scripts/check_response_partition_invariance.py', 'vcneb/nested_response.py',
               'vcneb/stationary_branch.py', 'vcneb/stationary_gap.py',
               'vcneb/active_curvature.py', 'vcneb/biaxial_curvature.py']
    return dict(format_version=1, check='full_joint_harmonic_partition_null_not_material_accuracy',
                analytic=check_partitions(), registered_geometry_chart_design_counts=chart_cost_bounds(),
                source_SHA256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                DFT_calls=0, calculator_calls=0, job_submissions=0, held_out_condition_accessed=False,
                full_Hessian_budget_registered=False, physical_parameters_changed=False,
                material_prediction_advantage_proven=False,
                next_action='Finish finite G2 labels; explicitly cost and register actual G3 probes. '
                            'Seek independent anharmonic/branch or cost evidence, not a coordinate-partition advantage.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh evidence namespace required')
    result = run_check()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps(dict(check=result['check'], stationary_points=result['analytic']['stationary_points'],
                          paired_points=result['analytic']['paired_points'],
                          maximum_errors=result['analytic']['maximum_errors'], DFT_calls=0)))


if __name__ == '__main__':
    main()
