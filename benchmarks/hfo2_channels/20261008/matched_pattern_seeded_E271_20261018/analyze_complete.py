"""Independently replay all six final T candidates from both portable archives."""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tarfile
import tempfile

REPOSITORY = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPOSITORY))
import numpy as np
from scripts.audit_hfo2_static_replica import sha256

BASE = REPOSITORY/'benchmarks/hfo2_channels/20261008'
ARCHIVES = dict(E270=('matched_pattern_release_E270_20261018', '8008b9c1c7f716324e9cfc28be866734739e635ec771a19b813d3720f434a9fc'),
                E271=('matched_pattern_seeded_E271_20261018', 'ff8077d1ded11438e709adf42193a0dc5f53b1126c37b13e1289bc6dd161536a'))


def restore(case, digest, destination):
    archive = case/'native_terminal.tar.gz'; receipt = json.loads((case/'native_terminal_receipt.json').read_bytes())
    if sha256(archive) != digest or receipt['archive_sha256'] != digest:
        raise ValueError('original terminal archive changed')
    destination.mkdir()
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            name = Path(member.name)
            if (name.is_absolute() or '..' in name.parts or name.parts[0] != 'native_terminal'
                    or not (member.isfile() or member.isdir())
                    or not (destination/name).resolve().is_relative_to(destination.resolve())):
                raise ValueError('unsafe native export member')
        tar.extractall(destination)
    output = destination/'native_terminal'
    files = json.loads((output/'files.json').read_bytes())
    if len(files) != receipt['exported_files']:
        raise ValueError('portable export count differs')
    for name, digest in files.items():
        path = output/name
        if not path.resolve().is_relative_to(output.resolve()) or sha256(path) != digest:
            raise ValueError('portable native bytes changed')
    return output, receipt


def analyze():
    original = BASE/ARCHIVES['E270'][0]
    spec = importlib.util.spec_from_file_location('original_four_pattern_reader', original/'review_arm.py')
    reader = importlib.util.module_from_spec(spec); spec.loader.exec_module(reader)
    rows = []; replays = []
    with tempfile.TemporaryDirectory(prefix='varneb_T_matched_') as temporary:
        exports = {tag: restore(BASE/name, digest, Path(temporary)/tag)
                   for tag, (name, digest) in ARCHIVES.items()}
        audit = json.loads((exports['E271'][0]/'audit.json').read_bytes())
        if audit['job_id'] != '28947723' or audit['new_native_SCF_calls'] != 14:
            raise ValueError('the actual finite terminal batch is required')
        for index, strain in enumerate((0., .005, .01)):
            values = {}
            for arm in ('A', 'J'):
                treatment = json.loads((exports['E271'][0]/str(index)/f'treatment_{arm}.json').read_bytes())
                source = 'E270' if treatment.get('origin') == 'E270_reused' else 'E271'
                r = reader.review(original, exports[source][0]/str(index), index, arm)
                if r['status'] != 'gradient_candidate_stability_unchecked':
                    raise ValueError('all six final local gradient candidates required')
                recorded = audit['rows'][index]['arms'][0 if arm == 'A' else 1]
                if (r['candidate_energy_eV_cell'] != recorded['candidate_energy_eV_cell']
                        or not np.isclose(r['released_gradient_norm_eV_A'], recorded['released_gradient_eV_A'], atol=2e-12, rtol=0)):
                    raise ValueError('independent portable native replay differs')
                r.update(candidate_source_batch=source, new_repair_attempts=treatment['new_attempts'])
                replays.append(r); values[arm] = r
            a, j = values['A'], values['J']
            rows.append(dict(strain_percent=100*strain, dE_A_minus_F_meV_fu=a['energy_change_meV_fu'],
                dE_J_minus_A_meV_fu=j['energy_change_meV_fu']-a['energy_change_meV_fu'],
                dE_J_minus_F_meV_fu=j['energy_change_meV_fu'],
                released_gradient_A_eV_A=a['released_gradient_norm_eV_A'],
                released_gradient_J_eV_A=j['released_gradient_norm_eV_A'],
                complete_allowed_gradient_J_eV_A=j['complete_allowed_gradient_norm_eV_A']))
        original_cost, repair_cost = exports['E270'][1], exports['E271'][1]
        costs = dict(original_attempts=original_cost['new_SCF_calls'], repair_attempts=repair_cost['new_SCF_calls'],
            total_new_SCF_calls=original_cost['new_SCF_calls']+repair_cost['new_SCF_calls'],
            allocated_core_hours=original_cost['allocation_core_hours']+repair_cost['allocation_core_hours'],
            native_proxy_core_hours=original_cost['native_call_proxy_core_hours']+repair_cost['native_call_proxy_core_hours'],
            all_earlier_curvature_and_path_preparation_costs_additional=True)
    return dict(status='all_six_T_same_pattern_candidates_replayed', rows=rows, native_candidate_replays=replays,
        costs=costs, archives={tag: digest for tag, (_, digest) in ARCHIVES.items()}, new_DFT_calls=0,
        same_coordinate_definitions_not_equal_native_values=True, old_labels_and_scores_unchanged=True,
        full_stationarity_stability_or_corrected_barrier_certified=False,
        local_energy_lowering_is_not_a_barrier_error_bound=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); result = analyze(); args.output.mkdir(parents=True, exist_ok=False)
    with (args.output/'analysis.json').open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    with (args.output/'energies.csv').open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result['rows'][0])); writer.writeheader(); writer.writerows(result['rows'])
    print(json.dumps({k: v for k, v in result.items() if k != 'native_candidate_replays'}))
