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
}


def load_sources(root):
    result = {}
    for name, (relative, digest) in SOURCES.items():
        raw = (Path(root)/relative).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError(f'archived {name} source bytes changed; no silent replacement')
        result[name] = json.loads(raw)
    return result


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
    return dict(format_version=1, evidence_kind='offline_archived_report_replay_not_new_native_audit',
        sources={k: dict(path=p, sha256=h) for k, (p, h) in SOURCES.items()},
        paths=paths, local_models=models, local_points=points,
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
                      ('T_local_all_predictions.csv', 'local_points')):
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
                          T_model_arms=len(result['local_models']), new_DFT=0)))


if __name__ == '__main__':
    main()
