"""Export the final HfO2 mechanism tables offline, without DFT or SSH.

Run from a source checkout:
    python -m examples.hfo2_mechanism_results --output NEW_DIRECTORY
Native-log replay remains available in the linked benchmark case.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASE = Path('benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021')


def report(root=ROOT):
    """Check the archived report and energy arithmetic; do not run a calculator."""
    case = Path(root) / CASE
    raw = (case / 'analysis/analysis.json').read_bytes()
    data = json.loads(raw)
    receipt = json.loads((case / 'completed/receipt.json').read_bytes())
    archive_hash = hashlib.sha256((case / 'completed/native_terminal.tar.gz').read_bytes()).hexdigest()
    if (data['archive_sha256'] != archive_hash or receipt['archive_sha256'] != archive_hash
            or data['physical_inputs_changed'] or data['new_DFT_calls']
            or data['contrasts']['status'] != 'four_candidate_contrasts'):
        raise ValueError('complete, unchanged-input, native-replayed cross data required')
    rows = data['rows']
    by = {(r['strain'], r['held_pattern']): r for r in rows}
    if len(rows) != 4 or set(by) != {(s, q) for s in (0., .01) for q in ('q0', 'q1')}:
        raise ValueError('all four crossed conditions required')
    for row in rows:
        if (row['released_gradient_eV_A'] > .03 or not math.isclose(
                row['relative_energy_meV_fu'],
                250 * (row['energy_eV_cell'] - row['initial_energy_eV_cell']), abs_tol=1e-9)):
            raise ValueError('local target or four-formula-unit energy conversion differs')
    energies = {key: row['relative_energy_meV_fu'] for key, row in by.items()}
    interaction = (energies[.01, 'q1'] - energies[.01, 'q0']
                   - energies[0., 'q1'] + energies[0., 'q0'])
    if not math.isclose(interaction, data['contrasts']['interaction_meV_fu'], abs_tol=1e-9):
        raise ValueError('finite interaction arithmetic differs')
    if sum(r['new_SCF_calls'] for r in rows) != data['new_native_SCF_calls']:
        raise ValueError('new and reused SCF accounting differs')
    return dict(local_cross=rows, whole_paths=data['whole_path_context'],
                contrasts=data['contrasts'], new_native_SCF_calls=data['new_native_SCF_calls'],
                source_sha256=hashlib.sha256(raw).hexdigest(), archive_sha256=archive_hash,
                replay_new_DFT_calls=0, native_logs_reparsed_by_this_export=False,
                local_points_are_not_new_barriers=True,
                stationary_TS_or_global_ordering_certified=False)


def export(output, root=ROOT):
    data = report(root)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    for name, rows in (('local_cross', data['local_cross']), ('whole_paths', data['whole_paths'])):
        with (output / f'{name}.csv').open('x', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    with (output / 'mechanism.json').open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write('\n')
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = export(args.output)
    print(json.dumps(dict(output=str(args.output), local_conditions=4,
                          new_DFT_calls=0, interaction_meV_fu=result['contrasts']['interaction_meV_fu'])))
