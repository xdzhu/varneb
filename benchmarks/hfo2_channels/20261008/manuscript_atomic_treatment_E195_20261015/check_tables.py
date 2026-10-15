"""Bind the matched atomic-treatment paragraph and table to completed native evidence."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / 'benchmarks/hfo2_channels/20261008/G3_atomic_complete_E194_20261015'
PINS = {
    'analysis.json': '3cc150ec3a6d0dbfa5f5879f7215a3692524309add3e5bb7f5abed384a197444',
    'atomic_native_complete.json': '8a2fdf2a2beec306bb2f39aed8b347019a4f049b2d8cbc42cd79de1fce922b92',
}


def source():
    for name, pin in PINS.items():
        if hashlib.sha256((BASE / name).read_bytes()).hexdigest() != pin:
            raise ValueError('completed original atomic-treatment bytes differ')
    return (json.loads((BASE / 'analysis.json').read_bytes()),
            json.loads((BASE / 'atomic_native_complete.json').read_bytes()))


def expected_rows(data, native):
    if (native['audit_mode'] != 'complete' or not native['all10points_covered']
            or native['completed_nodes_audited'] != 10 or native['gradient_candidates'] != 10
            or native['native_new_SCFs_reparsed'] != 36
            or native['sacct'].split('|')[1:3] != ['COMPLETED', '0:0']
            or data['atomic_audit_mode'] != 'complete' or len(data['rows']) != 4):
        raise ValueError('complete original matched cohort required')
    rows = []
    for r in data['rows']:
        stats = r['non_origin_coverage_statistics']['model_minus_observed']
        values = []
        for treatment in ('frozen', 'atomic', 'joint'):
            s = stats[treatment]['shape_eV_cell']
            if s['ids'] != [1,2,3,4] or not s['complete_four_axial_points'] or s['denominator'] != 4:
                raise ValueError('same four measured axial coordinates required')
            values.append(s['max_abs_meV_cell'])
        label = 'T' if r['channel'] == 'T' else 'Preserving'
        rows.append(f"| {label} | {r['probe_step_A']:.2f} | " +
                    ' | '.join(f'{v:.6f}' for v in values) + ' |')
    return rows


def require_rows(text, data, native):
    rows = expected_rows(data, native)
    if any(text.splitlines().count(row) != 1 for row in rows):
        raise ValueError('source-derived complete atomic table row missing, changed or duplicated')
    return rows


def check():
    paper = ROOT / 'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md'
    data, native = source()
    rows = require_rows(paper.read_text(encoding='utf-8'), data, native)
    return dict(status='complete_atomic_treatment_manuscript_rows_checked_not_independent_gain',
        source_pins=PINS, table_rows=rows, new_quantitative_rows=4,
        current_manuscript_sha256=hashlib.sha256(paper.read_bytes()).hexdigest(),
        new_DFT_or_scheduler_writes=0, structural_reorganization_applied=True,
        stable_surface_TS_independent_accuracy_or_JCTC_readiness_established=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = check()
    with args.output.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
