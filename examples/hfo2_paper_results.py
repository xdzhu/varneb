"""Export archived HfO2 paper tables offline; never starts DFT or SSH.

Run from a source checkout: python -m examples.hfo2_paper_results --output NEW_DIR
This replays audited report bytes, not raw ABACUS logs or unpublished arrays.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    'network': ('benchmarks/hfo2_channels/20261008/clamped_plus1_reversing_terminal_E148_20261014/analysis/network_final.json',
                '6b8a206d2199b6727601955cef32b86513e22805c7be5f309e66a4258602145b'),
    'local_errors': ('benchmarks/hfo2_channels/20261008/G3_T_independent_E258_20261017/score_E259.json',
                     '272bb11781be401e74ad333e4a211973410217ccce68a69fc9f2bd41f41a14cc'),
    'P8_native': ('benchmarks/hfo2_channels/20261008/G3_P8_execution_E266_20261018/native_complete_E267.json',
                  'dd245ac8f7335442ef5355b484b5c01dcba03ae3201ec1ad02131a3ae16ac654'),
    'P8_errors': ('benchmarks/hfo2_channels/20261008/G3_P8_execution_E266_20261018/score_E267.json',
                  'c51dec44d4e92e05ad824bd7952a9e1cee16d662476684675c9ae0c7489c2607'),
    'P8_forecasts': ('benchmarks/hfo2_channels/20261008/G3_P8_execution_E266_20261018/forecasts.json',
                     'd1d368372537a3d57d6a0cf38716ac8ea67210dbcac463d0dd2296a13e7b2046'),
}


def load_sources(root):
    result = {}
    for name, (relative, digest) in SOURCES.items():
        raw = (Path(root)/relative).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError(f'archived {name} source bytes changed; no silent replacement')
        result[name] = json.loads(raw)
    return result


def preserving_tables(sources):
    """Replay all revised P8 candidates and forecasts, not stable-sheet labels."""
    native, score, forecast = (sources[k] for k in ('P8_native', 'P8_errors', 'P8_forecasts'))
    if (not native['complete_cohort_audited'] or native['valid_candidate_labels'] != 8
            or native['native_SCFs_reparsed'] != 24 or len(native['native_records']) != 24
            or score['source_audit_sha256'] != SOURCES['P8_native'][1]
            or score['forecast_sha256'] != SOURCES['P8_forecasts'][1]
            or len(score['arms']) != 23 or score['valid_labels'] != 8):
        raise ValueError('complete pinned eight-target P8 cohort and all23 arms required')
    frozen = {}
    for model in forecast['new_models']:
        frozen[model['model_id']] = [(p['id'], p['prediction_raw_eV_cell'],
            p['prediction_origin_aligned_eV_cell']) for p in model['predictions']]
    for name in forecast['original_standard_preserving_record']['targets'][0]['models']:
        frozen[name] = [(p['id'], p['models'][name]['value'],
            p['models'][name]['seen_origin_calibrated_energy_change_eV_cell'])
            for p in forecast['original_standard_preserving_record']['targets']]
    for model in forecast['original_reference_preserving_records']:
        frozen[f"{model['reference']}@{model['step_A']:.2f}"] = [(p['id'],
            p['prediction_raw_eV_cell'], p['prediction_calibrated_eV_cell']) for p in model['points']]
    labels = {r['id']: r for r in native['rows']}
    if list(labels) != list(range(1, 9)) or sum(r['native_SCFs'] for r in labels.values()) != 24:
        raise ValueError('all original candidate IDs and charged SCFs required')
    candidates = []
    for row in native['rows']:
        label = row['label']
        if (row['status'] != 'gradient_candidate_stability_unchecked' or label['stability_certified']
                or label['released_gradient_norm_eV_A'] > .03 or abs(label['zeta_gradient_eV_A']) > 1e-4):
            raise ValueError('original dual-gradient candidate target, never stability certification')
        candidates.append(dict(point_id=row['id'], q1_A=row['q_A'][0], q2_A=row['q_A'][1],
            energy_change_eV_cell=label['energy_change_eV_cell'], native_SCFs=row['native_SCFs'],
            released_gradient_norm_eV_A=label['released_gradient_norm_eV_A'],
            zeta_gradient_eV_A=label['zeta_gradient_eV_A'], stability_certified=False))
    models, points = [], []
    if set(frozen) != {a['model'] for a in score['arms']}:
        raise ValueError('all four new and nineteen old model identities required')
    for arm in score['arms']:
        observed_forecasts = [(p['id'], p['prediction_raw_eV_cell'], p['prediction_origin_aligned_eV_cell']) for p in arm['points']]
        if observed_forecasts != frozen[arm['model']] or [p['id'] for p in arm['points']] != list(range(1, 9)):
            raise ValueError('predictions must remain identical to the pre-label panel')
        errors = {'raw': [], 'origin_aligned': []}
        for point in arm['points']:
            row = labels[point['id']]; observed = row['label']['energy_change_eV_cell']
            if point['observed_energy_change_eV_cell'] != observed:
                raise ValueError('every model must use the same native label')
            for mode in errors:
                predicted = point[f'prediction_{mode}_eV_cell']
                if predicted is not None: errors[mode].append((predicted-observed)*250.)
            points.append(dict(model=arm['model'], point_id=point['id'], q1_A=row['q_A'][0], q2_A=row['q_A'][1],
                observed_change_eV_cell=observed, raw_prediction_eV_cell=point['prediction_raw_eV_cell'],
                origin_aligned_prediction_eV_cell=point['prediction_origin_aligned_eV_cell'],
                abstention_reason=point['abstention_reason']))
        values = dict(model=arm['model'], kind=arm['kind'], denominator=8,
            supported_predictions=arm['prediction_support'], valid_labels=8)
        for mode, err in errors.items():
            maximum = max(map(abs, err)) if err else None
            rmse = math.sqrt(sum(v*v for v in err)/len(err)) if err else None
            expected = arm[mode]
            if (arm['denominator'] != 8 or len(err) != expected['n']
                    or (err and (not math.isclose(maximum, expected['max_absolute_meV_fu'], abs_tol=1e-10)
                        or not math.isclose(rmse, expected['RMSE_meV_fu'], abs_tol=1e-10)))):
                raise ValueError('P8 signed-error units, coverage or statistics do not replay')
            values[mode+'_max_error_meV_fu'] = maximum
            values[mode+'_RMSE_meV_fu'] = rmse
        models.append(values)
    return candidates, models, points


def replay(root=ROOT):
    """Recompute table values/coverage from every archived row, not winners."""
    sources = load_sources(root)
    network, scores = sources['network'], sources['local_errors']
    required = {(name, strain) for name in ('T', 'M', 'preserving', 'reversing') for strain in (0., .01)}
    if {(r['candidate'], r['substrate_strain']) for r in network['rows']} != required or len(network['rows']) != 8:
        raise ValueError('complete eight-path cohort required')
    paths, common = [], {}
    for row in network['rows']:
        profile = row['raw_relative_energy_meV_fu']
        forward, reverse, delta = row['forward_meV_fu'], row['reverse_meV_fu'], row['final_minus_initial_meV_fu']
        if (not 0 <= row['ordinary_fmax_eV_A'] <= .10 or len(profile) != 9
                or not math.isclose(max(profile), forward, abs_tol=1e-10)
                or not math.isclose(forward-reverse, delta, abs_tol=1e-10)
                or not math.isclose(profile[-1], delta, abs_tol=1e-10)):
            raise ValueError('ordinary residual, normalization or complete energy profile differs')
        strain = row['substrate_strain']
        initial = row['matched_initial_energy_eV_cell']
        if strain in common and common[strain] != initial:
            raise ValueError('one exact common initial energy per strain required')
        common[strain] = initial
        paths.append(dict(channel=row['candidate'], substrate_strain=strain,
            sampled_forward_meV_fu=forward, sampled_reverse_meV_fu=reverse,
            final_minus_initial_meV_fu=delta, ordinary_fmax_eV_A=row['ordinary_fmax_eV_A'],
            peak_image=row['sampled_peak_image'], barrier_error_bound_meV_fu=row['barrier_error_bound_meV_fu'],
            terminal_audit_sha256=row['terminal_audit_sha256']))
    paths.sort(key=lambda r: (r['substrate_strain'], r['channel']))
    if len(scores['arms']) != 19 or scores['T_denominator'] != 8:
        raise ValueError('all nineteen arms and eight original T tests required')
    models, points, labels = [], [], {}
    for arm in scores['arms']:
        if [p['id'] for p in arm['points']] != list(range(1, 9)):
            raise ValueError('all original test IDs required, including abstentions')
        errors = {'calibrated': [], 'uncalibrated': []}
        for point in arm['points']:
            observed = point['observed_change_eV_cell']
            if point['id'] in labels and labels[point['id']] != observed:
                raise ValueError('all models must see the same exact independent label')
            labels[point['id']] = observed
            for mode, key in (('calibrated', 'calibrated_eV_cell'), ('uncalibrated', 'raw_eV_cell')):
                value = point[key]
                if value is not None:
                    errors[mode].append((value-observed)*250.)  # 1000 meV/eV / 4 formula units
            points.append(dict(model=arm['name'], point_id=point['id'], q1_A=point['q_A'][0], q2_A=point['q_A'][1],
                observed_change_eV_cell=observed, raw_prediction_eV_cell=point['raw_eV_cell'],
                origin_aligned_prediction_eV_cell=point['calibrated_eV_cell'], abstention_reason=point['reason']))
        values = dict(model=arm['name'], kind=arm['kind'], denominator=8,
                      supported_predictions=arm['prediction_support'], valid_labels=arm['valid_labels'])
        for mode, err in errors.items():
            if len(err) != arm[mode]['n']:
                raise ValueError('coverage count differs; abstentions must not be dropped')
            maximum = max(map(abs, err)) if err else None
            rmse = math.sqrt(sum(v*v for v in err)/len(err)) if err else None
            if err and (not math.isclose(maximum, arm[mode]['max_absolute_meV_per_fu'], abs_tol=1e-10)
                        or not math.isclose(rmse, arm[mode]['RMSE_meV_per_fu'], abs_tol=1e-10)):
                raise ValueError('independent errors do not replay in meV/formula unit')
            values[mode+'_max_error_meV_fu'] = maximum
            values[mode+'_RMSE_meV_fu'] = rmse
        models.append(values)
    P_candidates, P_models, P_points = preserving_tables(sources)
    return dict(format_version=2, evidence_kind='offline_archived_report_replay_not_new_native_audit',
        sources={k: dict(path=p, sha256=h) for k, (p, h) in SOURCES.items()},
        paths=paths, local_models=models, local_points=points,
        preserving_candidates=P_candidates, preserving_models=P_models, preserving_points=P_points,
        preserving_native_SCFs=24, preserving_allocation_core_hours=sources['P8_native']['allocation_core_hours'],
        preserving_prior_preparation_costs_additional=True,
        preserving_primary_four_common_ids=sources['P8_errors']['primary_four_arm_common_ids'],
        preserving_all23_common_ids=sources['P8_errors']['all23_common_ids'],
        preserving_target='revised single-start dual-gradient candidate, not original stable-plane validation',
        preserving_original_2D_refusal_preserved=True,
        original_T_tests_seen=True, local_test_native_SCFs=scores['new_native_SCFs'],
        local_test_allocation_core_hours=scores['allocation_core_hours'],
        local_test_prior_preparation_costs_additional=True,
        physical_plus005_results_included=False, stable_TS_or_error_bound_or_method_advantage_certified=False,
        new_DFT_or_scheduler_operations=0)


def export(output, root=ROOT):
    """Write a fresh local result directory, refusing any existing target."""
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError('choose a new output directory; archived results are never overwritten')
    result = replay(root)  # validate every input before creating output
    output.mkdir(parents=True, exist_ok=False)
    for name, key in (('channel_comparison.csv', 'paths'), ('T_local_model_errors.csv', 'local_models'),
                      ('T_local_all_predictions.csv', 'local_points'),
                      ('P_local_candidates.csv', 'preserving_candidates'),
                      ('P_local_model_errors.csv', 'preserving_models'),
                      ('P_local_all_predictions.csv', 'preserving_points')):
        with (output/name).open('x', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(result[key][0]))
            writer.writeheader(); writer.writerows(result[key])
    with (output/'replay.json').open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = export(args.output)
    print(json.dumps(dict(output=str(args.output), paths=len(result['paths']),
                          T_model_arms=len(result['local_models']),
                          P_model_arms=len(result['preserving_models']), new_DFT=0)))


if __name__ == '__main__':
    main()
