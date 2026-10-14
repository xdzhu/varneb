"""Eight native-audited training paths: sampled responses, not certified TSs.

Only the fixed0/+1% matrix is opened. Unknown barrier errors stay unknown;
neither an independent prediction nor a globally lowest mechanism is inferred.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from scripts.analyze_hfo2_clamped_reference import terminal_images
from scripts.audit_hfo2_manuscript_work_cohort import BASE, CASES
from scripts.audit_hfo2_static_replica import sha256

# Explicit final job/edge/condition identities, not an energy-based label fit.
IDENTITIES = {
    '28722320': ('M', 0.), '28692775': ('T', 0.),
    '28692776': ('reversing', 0.), '28727382': ('preserving', 0.),
    '28709788': ('T', .01), '28727368': ('M', .01),
    '28812952': ('preserving', .01), '28852667': ('reversing', .01),
}


def summarize(groups):
    if set(groups) != {0., .01} or any(set(g) != {'M', 'T', 'preserving', 'reversing'} for g in groups.values()):
        raise ValueError('the exact two-condition/four-candidate training matrix is required')
    if any(isinstance(row['forward_meV_fu'], bool) or not np.isfinite(row['forward_meV_fu'])
           or row['forward_meV_fu'] < 0 for group in groups.values() for row in group.values()):
        raise ValueError('finite nonnegative sampled barriers required')
    changes = {key: groups[.01][key]['forward_meV_fu']-groups[0.][key]['forward_meV_fu']
               for key in ('M', 'T', 'preserving', 'reversing')}
    minima = {}
    for strain, group in groups.items():
        switching = min(('preserving', 'reversing'), key=lambda k: group[k]['forward_meV_fu'])
        escape = min(('M', 'T'), key=lambda k: group[k]['forward_meV_fu'])
        bsw = group[switching]['forward_meV_fu']; besc = group[escape]['forward_meV_fu']
        minima[str(strain)] = dict(sampled_switching_candidate=switching, sampled_escape_candidate=escape,
            sampled_switching_meV_fu=bsw, sampled_escape_meV_fu=besc, sampled_S_escape_minus_switch_meV_fu=besc-bsw,
            physical_ordering_resolved=False)
    zero, plus = minima['0.0'], minima['0.01']
    ds = plus['sampled_switching_meV_fu']-zero['sampled_switching_meV_fu']
    de = plus['sampled_escape_meV_fu']-zero['sampled_escape_meV_fu']
    return dict(candidate_changes_meV_fu=changes, sampled_minima=minima,
        sampled_minimum_switching_change_meV_fu=ds,
        sampled_minimum_escape_change_meV_fu=de,
        sampled_S_change_meV_fu=plus['sampled_S_escape_minus_switch_meV_fu']-zero['sampled_S_escape_minus_switch_meV_fu'],
        raw_H1_favorable_signature=bool(ds < 0 and de >= 0),
        H1_primary_status='not_established; barrier_error_bounds_unmeasured',
        H1_error_certified_rejection=False, barrier_error_bounds_measured=False,
        global_candidate_coverage_or_kinetic_lifetime_certified=False,
        G3_modes_or_holdout_edges_selected=False, physical_holdout_generated_or_read=False)


def audit():
    groups = {0.: {}, .01: {}}; rows = []; common = {0.: [], .01: []}
    for label, case, pin, _, _ in CASES:
        images, observation, residual, receipt = terminal_images(BASE/case, pin)
        job = receipt['job_id']
        if job not in IDENTITIES:
            raise ValueError('unexpected training terminal')
        key, strain = IDENTITIES[job]
        if key in groups[strain] or len(images) != 9:
            raise ValueError('duplicate edge or wrong seven-moving/two-fixed image convention')
        if any(Counter(image.get_chemical_symbols()) != {'Hf': 4, 'O': 8} for image in images):
            raise ValueError('same12-atom/four-formula-unit composition required')
        b = observation['mechanical_boundary']
        if (observation['pressure_GPa'] != 0 or observation['climb'] is not False
                or b['kind'] != 'clamped_plane' or b['fixed_cell_rows'] != [0, 1]
                or b['allow_tilt'] is not True or b['cell_dofs'] != 3
                or b['cell_scale_A'] != 5.12968067458423):
            raise ValueError('same pressure and registered tilt-open clamped boundary required')
        energy = np.asarray([image.get_potential_energy() for image in images])
        relative = (energy-energy[0])*250
        forward = float(np.max(relative)); well = float(relative[-1]); reverse = forward-well
        if (not np.isclose(forward, observation['sampled_forward_barrier_meV_fu'], atol=1e-9, rtol=0)
                or not np.isclose(reverse, observation['sampled_reverse_barrier_meV_fu'], atol=1e-9, rtol=0)):
            raise ValueError('native energy/enthalpy or formula-unit convention differs')
        raw = observation['raw_image_evaluations'][0]
        common[strain].append((raw['raw_log_sha256'], raw['input_sha256']['STRU'], float(energy[0])))
        plane = np.asarray(b['reference_cell_A'])[:2]
        for image in images:
            np.testing.assert_allclose(image.cell.array[:2], plane, atol=1e-12, rtol=0)
        row = dict(candidate=key, substrate_strain=strain, source_job_id=job,
            terminal_audit_sha256=pin, forward_meV_fu=forward, reverse_meV_fu=reverse,
            final_minus_initial_meV_fu=well, sampled_peak_image=int(np.argmax(relative)),
            ordinary_fmax_eV_A=residual['fmax_eV_A'], raw_relative_energy_meV_fu=relative.tolist(),
            fixed_substrate_rows_A=plane.tolist(), matched_initial_energy_eV_cell=float(energy[0]),
            barrier_error_bound_meV_fu=None)
        groups[strain][key] = row; rows.append(row)
    if len(rows) != 8 or any(len(set(records)) != 1 for records in common.values()):
        raise ValueError('all four candidates must share one exact initial cache within each condition')
    zero_plane = np.asarray(groups[0.]['M']['fixed_substrate_rows_A'])
    for key in groups[0.]:
        np.testing.assert_allclose(groups[0.][key]['fixed_substrate_rows_A'], zero_plane, atol=1e-12, rtol=0)
        np.testing.assert_allclose(groups[.01][key]['fixed_substrate_rows_A'], zero_plane*1.01, atol=1e-12, rtol=0)
    return dict(status='complete_ordinary_G2_sampled_network_not_TS_or_error_certified',
        source_script_sha256=sha256(Path(__file__)), rows=rows, response=summarize(groups),
        ordinary_audited_training_paths=8, portable_raw_frame_records_reparsed=72,
        pressure_GPa=0., initial_cache_groups=[dict(strain=s, native_identity=list(v[0]), occurrences=len(v))
                                             for s, v in common.items()],
        new_DFT_calls=0, scheduler_mutations=0, physical_parameters_changed=False,
        new_mode_or_reference_selection=False, unpublished_author_arrays_required=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): raise FileExistsError('fresh network output required')
    result = audit(); args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(result['response']))


if __name__ == '__main__':
    main()
