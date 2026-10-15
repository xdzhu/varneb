"""Prospective elementary controls for eight uncomputed Q points per face.

Only the original two training centers are read. This is not a freeze of
B0--B5 strain/barrier predictions, a stability certificate, or a DFT runner.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import itertools
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from scripts import hfo2_G3_conditional_faces as faces
from vcneb.quadratic_reduction import condition_quadratic_energy

BASE = ROOT / 'benchmarks/hfo2_channels/20261008'
CONTROL = BASE / 'G3_coarsening_E184_20261015/analyze.py'
CONTROL_SHA = '9c8c8a68cd77b8a468f90df7a839afd249461f74bc84a1f3fdb723d1804209e3'
if hashlib.sha256(CONTROL.read_bytes()).hexdigest() != CONTROL_SHA:
    raise ValueError('original elementary control implementation differs')
spec = importlib.util.spec_from_file_location('E184_elementary', CONTROL)
elementary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(elementary)
MODELS = list(elementary.MODELS)
FIRST_SHA = '283f2b9de82f6162b58a214ab3f19722f835921eba016d942e04a1dd4164db68'


def lf_text_sha(path):
    """Git-normalized text identity, not a native DFT/binary file hash."""
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def targets():
    indices = sorted((k for k in itertools.product((-1, 0, 1), repeat=2) if k != (0, 0)),
                     key=lambda k: (abs(k[0]) + abs(k[1]), k))
    return [dict(id=i + 1, half_grid_index=list(k), q_A=[v * .05 for v in k])
            for i, k in enumerate(indices)]


def fixed_start(data, training, target):
    """Nearest central training Q; exact integer distances, then original ID.

    Test results never become continuation parents. Unsafe starts are not
    repaired, and the original single-start branch definition is preserved.
    """
    a = data['axes']
    key = np.asarray(target['half_grid_index'])
    ordered = sorted(training, key=lambda r: (
        sum((2 * i - j) ** 2 for i, j in zip(r['original_index'], key)), r['id']))
    for row in ordered:
        previous = np.asarray(row['coordinates_internal_A'])
        start = a @ np.asarray(target['q_A']) + previous - a @ (a.T @ previous)
        if faces.canary.geometry_guard(data['chart'], data['basis'], start):
            return row['id'], start
    raise ValueError('no registered safe training start; no repair or new branch')


def freeze():
    source = elementary.pinned(elementary.SOURCE, elementary.PIN)
    rows = []
    for channel in ('T', 'preserving'):
        data = faces.load_inputs(channel)
        proof_name, proof_sha = elementary.PROOFS[channel]
        proof = elementary.pinned(BASE / proof_name, proof_sha)
        original_specs = faces.ordered_points()
        selected = [p for p in original_specs if max(map(abs, p['index'])) <= 1]
        if len(selected) != 9:
            raise ValueError('fixed central nine training points required')
        observed_rows = next(r for r in source['rows'] if r['channel'] == channel
                             and r['reference'] == 'T_all_optical' and r['step_A'] == .01)
        training, labels = [], {}
        for p in selected:
            if p['id'] == 0:
                candidate = data['origin_point']
                calls = 0  # old canary preparation remains additional
            else:
                native = next(r for r in proof['candidate_rows'] if r['id'] == p['id'])
                candidate, calls = native['candidate'], native['completed_SCFs']
                if native['q_A'] != p['q_A']:
                    raise ValueError('training observable identity differs')
            if candidate is None:
                raise ValueError('missing central training label; no fill or replacement')
            energy = candidate['energy_eV_cell'] - data['center_energy']
            if abs(energy - observed_rows['points'][p['id']]['DFT_energy_change_eV_cell']) > 1e-12:
                raise ValueError('source scalar and audited candidate energy disagree')
            training.append(dict(id=p['id'], original_index=p['index'], q_A=p['q_A'],
                energy_change_eV_cell=energy, historical_face_SCFs=calls,
                coordinates_internal_A=candidate['coordinates_internal_A']))
            labels[tuple(2 * i for i in p['index'])] = energy
        beta, rank = elementary.fit_quadratic(labels)
        if rank != 6 or beta is None:
            raise ValueError('identified full quadratic required')
        b, a = data['basis'], data['axes']
        complement = np.linalg.svd(a.T, full_matrices=True)[2][2:].T
        models = {}
        for block in data['curvature']['blocks']:
            raw = np.asarray(block['raw_projected_eV_A2'])[:36, :36]
            model = condition_quadratic_energy((raw + raw.T) / 2,
                b.T @ data['center_gradient'], data['center_energy'], a, complement,
                stability_floor=0.)
            models[f"original_H_{block['step_A']:.2f}"] = model
        if set(models) != {'original_H_0.01', 'original_H_0.02'}:
            raise ValueError('both original measured steps required')
        origin_y = labels[(0, 0)]
        origin_models = {k: elementary.interpolate((0, 0), labels, k)['value']
                         for k in ('bilinear', 'tensor_quadratic')}
        origin_models['refitted_total_quadratic'] = float(beta[0])
        origin_models.update({k: m.energy_change([0., 0.]) for k, m in models.items()})
        predictions = []
        for target in targets():
            q, index = target['q_A'], target['half_grid_index']
            if any(q == p['q_A'] for p in original_specs):
                raise ValueError('test coordinate was already a training target')
            parent, start = fixed_start(data, training, target)
            values = {k: elementary.interpolate(index, labels, k)
                      for k in ('bilinear', 'tensor_quadratic')}
            values['refitted_total_quadratic'] = dict(
                value=float(elementary.features(index) @ beta), reason=None)
            for k, model in models.items():
                x = model.full_displacement(q)
                eligible = faces.canary.geometry_guard(data['chart'], b, x)
                values[k] = dict(value=model.energy_change(q) if eligible else None,
                    reason=None if eligible else 'outside_original_geometry_guard',
                    formal_energy_change_eV_cell=model.energy_change(q),
                    formal_internal_norm_A=float(np.linalg.norm(x)))
            for k, value in values.items():
                value['origin_model_energy_change_eV_cell'] = origin_models[k]
                value['origin_bias_eV_cell'] = origin_models[k] - origin_y
                value['seen_origin_calibrated_energy_change_eV_cell'] = (
                    None if value['value'] is None else value['value'] - origin_models[k] + origin_y)
                value['relative_to_seen_origin_eV_cell'] = (
                    None if value['value'] is None else value['value'] - origin_models[k])
            predictions.append(dict(**target, training_parent_id=parent,
                starting_internal_coordinates_A=start.tolist(),
                starting_full_coordinates_A=(b @ start).tolist(),
                starting_internal_norm_A=float(np.linalg.norm(start)),
                models=values, independent_DFT_label=None))
        rows.append(dict(channel=channel, training=training, targets=predictions,
            formula_units=4, raw_peak_reference_energy_eV_cell=data['center_energy'],
            observed_origin_energy_change_eV_cell=origin_y,
            refitted_total_quadratic_coefficients=beta.tolist(),
            energy_interpolants_training_labels=9,
            original_H_preparation_SCFs_per_center=148,
            historical_face_SCFs_at_nine_training_targets=sum(r['historical_face_SCFs'] for r in training),
            original_whole_training_face_SCFs=proof['new_SCFs_reparsed'],
            source_native_face_sha256=proof_sha, original_input_artifact_sha256=list(faces.PINS[channel])))
    code = sorted(set((*faces.canary.CODE, 'scripts/hfo2_G3_conditional_faces.py',
        'vcneb/quadratic_reduction.py', CONTROL.relative_to(ROOT).as_posix(),
        Path(__file__).resolve().relative_to(ROOT).as_posix())))
    first = elementary.pinned(Path(__file__).with_name('forecast.json'), FIRST_SHA)
    if rows != first['channels']:
        raise ValueError('metadata amendment must not change any original numerical prediction or start')
    return dict(status='prospective_local_elementary_predictions_locked_before_test_labels',
        revision=2, original_forecast_sha256=FIRST_SHA,
        revision_reason='canonical LF text-code hashes for Git archive portability; all original numerical rows and policies unchanged',
        channels=rows, models=MODELS, scalar_predictions=sum(len(p['models']) for r in rows for p in r['targets']),
        energy_unit='eV per 12-atom Hf4O8 cell; 4 formula units',
        training_scalar_source_sha256=elementary.PIN,
        code_lf_text_sha256={name: lf_text_sha(ROOT / name) for name in code},
        selection='all eight nonzero points of the central {-0.05,0,0.05} A square; no energy selection',
        branch_start_policy='nearest central training point by exact integer Q distance; ties original ID; one safe start; never use test labels',
        scoring=dict(primary='seen-origin-aligned signed error, max absolute error and RMSE; all five arms',
            secondary='uncalibrated signed errors and origin biases, never discard either step',
            reporting_target_meV_formula_unit=1., reporting_target_meV_12atom_cell=4.,
            target_is_engineering_resolution_not_DFT_error_bound_or_literature_standard=True,
            target_selected_after_seen_training_before_unseen_test_labels=True,
            coverage_denominator_per_face=8,
            cohort='same valid candidate labels for all arms; also report per-arm abstention and missing/failed labels out of eight',
            no_label_or_model_support_means_null_not_zero=True,
            uncertainty_overlap_means_no_resolved_advantage=True,
            no_replacement_refit_step_selection_or_branch_repair_after_test_labels=True),
        test_execution_policy=dict(independent_points_reserved_per_face=8,
            independent_points_computed_per_face=0, maximum_new_attempts_per_point=12,
            maximum_new_attempts_two_faces=192, optimizer='safeguarded_bfgs',
            initial_hessian='original step 0.01 symmetrized full internal Hessian',
            gradient_tolerance_eV_A=.03, initial_trust_radius_A=.05,
            maximum_internal_displacement_A=.5, minimum_distance_A=1.6,
            maximum_iterations=200, no_automatic_retry_or_rollover=True,
            prerequisites='complete native B3 and released-Q0 curvature audits; both raw R34 steps positive as necessary local evidence, not global certification; final candidate predictions locked before labels',
            labels_are_single_start_gradient_candidates_not_global_minimum_guarantees=True,
            no_new_allocation_while_two_existing_study_allocations_committed=True),
        modeling_training_commitments_unchanged={'T':34,'preserving':35},
        all_training_labels_seen_before_choice=True, local_test_labels_seen=False,
        standard_interpolation_and_Schur_algebra_are_not_innovations=True,
        no_equal_information_or_equal_preparation_cost_claim=True,
        origin_canary_Gamma_original296_G2_and_all_branch_preparation_additional=True,
        strong_T_Cmma_full_B0_B5_physical_strain_predictions_and_domain_freeze_complete=False,
        whole_path_error_gain_stability_TS_or_JCTC_validation_claimed=False,
        physical_reserved_plus005_geometry_generated_or_read=False,
        new_DFT_or_scheduler_writes=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('new lock file required; never replace a forecast')
    result = freeze()
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(scalar_predictions=result['scalar_predictions'],
        independent_points_reserved_per_face=8, new_DFT=0)))


if __name__ == '__main__':
    main()
