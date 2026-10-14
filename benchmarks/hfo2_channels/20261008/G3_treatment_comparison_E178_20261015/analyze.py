"""Signed frozen/atomic/joint training comparisons from hash-bound audit reports.

No calculator, scheduler, coefficient fit, winning-step selection or raw SCF
reparse. Missing atomic candidates do not hide known frozen/joint differences.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
BASE = ROOT / 'benchmarks/hfo2_channels/20261008'
TREATMENTS = ('frozen', 'atomic', 'joint')
PAIRS = (('atomic', 'frozen'), ('joint', 'atomic'), ('joint', 'frozen'))
PINS = {
    'G3_atomic_audit_E176_20261015/forecast.json':
        '321ba7007c7f950166aab4eaa2dbe579773164968f04a94034a8a942bee1318d',
    'G3_frozen_release_E173_20261015/analysis.json':
        '7d8860690556197f8b38d4f67f80562159320f31ae26d61f95198ec3fb801dfa',
    'G3_atomic_release_E175_20261015/config.json':
        '5bd401c6ca4044ac654f2b19a22371b20c17d3736d207ff7044336d9ed463f38',
    'G3_atomic_release_E175_20261015/controls_summary.json':
        'f49c662d705d4733dbcf1d78bb4e0260a0165df4c0a0df0af3e468f56fc8ffc8',
    'G3_atomic_audit_E176_20261015/audit.py':
        'c201592b999037486de3db628935fcb86edcf8fcbc5ddb5dd969f2f905873793',
}
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def finite(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError('finite numeric energy/cost required; missing values use None')
    return float(value)


def changes(energies, origins):
    """Every available pair retained; shape requires only that pair's origins."""
    if set(energies) != set(TREATMENTS) or set(origins) != set(TREATMENTS):
        raise ValueError('exact frozen/atomic/joint treatments required')
    for value in (*energies.values(), *origins.values()):
        if value is not None:
            finite(value)
    out = {}
    for b, a in PAIRS:
        absolute = None if any(energies[t] is None for t in (a, b)) else energies[b] - energies[a]
        shape = None if absolute is None or any(origins[t] is None for t in (a, b)) else (
            (energies[b] - origins[b]) - (energies[a] - origins[a]))
        out[b + '_minus_' + a] = dict(absolute_eV_cell=absolute, shape_eV_cell=shape)
    for field in ('absolute_eV_cell', 'shape_eV_cell'):
        terms = [out[b + '_minus_' + a][field] for b, a in PAIRS]
        if all(v is not None for v in terms) and abs(terms[0] + terms[1] - terms[2]) > 2e-11:
            raise ValueError('signed treatment identity differs')
    return out


def compare(model, observed, model_origins, observed_origins):
    if any(model[t] is None or model_origins[t] is None for t in TREATMENTS):
        raise ValueError('all formal model energies/origins required')
    expected = changes(model, model_origins)
    actual = changes(observed, observed_origins)
    errors = {}
    for t in TREATMENTS:
        absolute = None if observed[t] is None else model[t] - observed[t]
        shape = None if observed[t] is None or observed_origins[t] is None else (
            (model[t] - model_origins[t]) - (observed[t] - observed_origins[t]))
        errors[t] = dict(absolute_eV_cell=absolute, shape_eV_cell=shape)
    mismatch = {}
    for b, a in PAIRS:
        name = b + '_minus_' + a
        mismatch[name] = {}
        for field in ('absolute_eV_cell', 'shape_eV_cell'):
            v = actual[name][field]
            difference = None if v is None else expected[name][field] - v
            if difference is not None and abs(difference - (errors[b][field] - errors[a][field])) > 2e-11:
                raise ValueError('pointwise treatment-error decomposition differs')
            mismatch[name][field] = difference
    return dict(observed_changes=actual, model_changes=expected,
                model_minus_observed=errors, missed_treatment_change=mismatch)


def summarize(points):
    """Non-origin statistics with visible coverage, never rank partial cohorts."""
    out = {}
    for group in ('observed_changes', 'model_minus_observed', 'missed_treatment_change'):
        out[group] = {}
        for name in points[0][group]:
            out[group][name] = {}
            for field in ('absolute_eV_cell', 'shape_eV_cell'):
                available = [(p['id'], p[group][name][field]) for p in points
                             if p['id'] != 0 and p[group][name][field] is not None]
                out[group][name][field] = dict(ids=[i for i, _ in available], denominator=4,
                    complete_four_axial_points=len(available) == 4,
                    max_abs_meV_cell=None if not available else max(abs(v) for _, v in available) * 1000)
    return out


def validate_atomic(report, config):
    """Validate a supplied native audit, not replace its raw-SCF/optimizer checks."""
    if (report['status'] != 'read_only_atomic_treatment_native_prefix_not_stability_or_independent_prediction'
            or report['job_id'] != '28867177'
            or report['audit_mode'] not in ('prefix', 'complete')
            or report['new_DFT_or_scheduler_writes'] != 0
            or any(report[k] is not False for k in ('original_runtime_or_physics_changed',
                'physical_reserved_geometry_generated_or_read',
                'conditional_stability_TS_acceleration_or_independent_prediction_certified'))):
        raise ValueError('original read-only atomic audit required')
    account = report['sacct'].split('|')
    if len(account) != 7 or account[0] != '28867177' or account[6] != 'hfacnormal01':
        raise ValueError('original accounting identity required')
    rows = report['candidate_rows']
    expected = [(c, p) for c in config['policy']['channels'] for p in config['channels'][c]['points']]
    if (len(rows) > 10 or type(report['completed_nodes_audited']) is not int
            or report['completed_nodes_audited'] != len(rows)
            or report['all10points_covered'] is not (len(rows) == 10)):
        raise ValueError('original atomic prefix coverage required')
    for row, (channel, point) in zip(rows, expected):
        if (row['channel'] != channel or type(row['id']) is not int
                or row['id'] != point['id'] or row['q_A'] != point['q_A']
                or type(row['new_completed_SCFs']) is not int
                or not 0 <= row['new_completed_SCFs'] <= 12
                or finite(row['sum_completed_call_seconds']) < 0
                or row['replay'].get('exact_request_order_optimizer_replay') is not True):
            raise ValueError('atomic original-Q/order/budget/replay differs')
        if len(row['q_A']) != 2:
            raise ValueError('two finite numeric Q coordinates required')
        for v in row['q_A']:
            finite(v)
        p = row['candidate']
        if p is not None:
            finite(p['energy_eV_cell'])
            if (not 0 <= finite(p['atomic_released_gradient_norm_eV_A']) <= .03
                    or finite(p['joint_released_gradient_norm_eV_A']) < 0
                    or p['released_stability_unchecked'] is not True):
                raise ValueError('atomic gradient candidate, not joint/stability claim, required')
    if (sum(r['candidate'] is not None for r in rows) != report['gradient_candidates']
            or sum(r['new_completed_SCFs'] for r in rows) != report['native_new_SCFs_reparsed']
            or len(report['native_sources']) != report['native_new_SCFs_reparsed']
            or (rows and report['old_native_cache_preflight'] is None)):
        raise ValueError('atomic audit candidate/native accounting differs')
    pending = None if len(rows) == 10 else dict(channel=expected[len(rows)][0], id=expected[len(rows)][1]['id'])
    if ((pending is None) != (report['pending'] is None) or pending is not None and
            any(report['pending'].get(k) != v for k, v in pending.items())):
        raise ValueError('atomic pending point/coverage differs')
    if report['audit_mode'] == 'complete' and (len(rows) != 10 or account[1:3] != ['COMPLETED', '0:0']):
        raise ValueError('complete atomic audit requires normal completion/all10points')
    return {(r['channel'], r['id']): r for r in rows}


def analyze(atomic_path, atomic_sha):
    if sha(atomic_path) != atomic_sha:
        raise ValueError('explicit atomic-audit byte hash differs')
    for name, pin in PINS.items():
        if sha(BASE / name) != pin:
            raise ValueError('original source bytes differ: ' + name)
    load = lambda name: json.loads((BASE / name).read_bytes())
    forecast = load('G3_atomic_audit_E176_20261015/forecast.json')
    old = load('G3_frozen_release_E173_20261015/analysis.json')
    config = load('G3_atomic_release_E175_20261015/config.json')
    controls = load('G3_atomic_release_E175_20261015/controls_summary.json')
    atomic = json.loads(atomic_path.read_bytes())
    index = validate_atomic(atomic, config)
    # Reuse the hash-bound original reference metadata; no calculator is created.
    from scripts.hfo2_G3_conditional_faces import load_inputs
    centers = {c: load_inputs(c)['center_energy'] for c in config['policy']['channels']}
    rows = []
    for f in forecast['rows']:
        channel = f['channel']
        seen = next(r for r in old['rows'] if r['channel'] == channel
                    and r['curvature_probe_step_A'] == f['probe_step_A'])
        points = {p['point_id']: p for p in seen['points']}
        ao = index.get((channel, 0), {}).get('candidate')
        origins = dict(frozen=0., atomic=None if ao is None else ao['energy_eV_cell'] - centers[channel],
                       joint=seen['raw_origin_DFT_release_meV_cell'] / 1000)
        mo = dict(frozen=0., atomic=f['formal_atomic_origin_release_meV_cell'] / 1000,
                  joint=f['formal_joint_origin_release_meV_cell'] / 1000)
        out = []
        for p in f['points']:
            point_id = p['id']
            oldp = points.get(point_id)
            a = index.get((channel, point_id))
            candidate = None if a is None else a['candidate']
            obs = dict(frozen=0. if point_id == 0 else oldp['DFT_frozen_energy_change_eV_cell'],
                       atomic=None if candidate is None else candidate['energy_eV_cell'] - centers[channel],
                       joint=origins['joint'] if point_id == 0 else oldp['DFT_released_energy_change_eV_cell'])
            model = {t: p['model_' + t + '_energy_change_eV_cell'] for t in TREATMENTS}
            out.append(dict(id=point_id, q_A=p['q_A'], model=model, observed=obs,
                atomic_status='not_processed' if a is None else ('gradient_candidate_stability_unchecked'
                    if candidate is not None else 'processed_no_gradient_candidate'),
                atomic_reason=None if a is None else a['reason'],
                **compare(model, obs, mo, origins)))
        rows.append(dict(channel=channel, probe_step_A=f['probe_step_A'],
                         observed_origins=origins, model_origins=mo, points=out,
                         non_origin_coverage_statistics=summarize(out)))
    if len(rows) != 4 or sum(len(r['points']) for r in rows) != 20:
        raise ValueError('all two centers/two original steps/five points required')
    return dict(status='signed_training_treatment_comparison_not_independent_prediction', rows=rows,
        energy_unit='eV per Hf4O8 cell (four HfO2 formula units)', formula_units=4,
        atomic_audit_observed_CST=atomic['observed_CST'], atomic_audit_sha256=atomic_sha,
        atomic_source_sacct=atomic['sacct'], atomic_audit_mode=atomic['audit_mode'],
        source_pins=dict(PINS), new_raw_SCF_reparse=False, new_DFT_or_scheduler_writes=0,
        no_step_or_point_selection=True, fitted_coefficients=0,
        costs=dict(frozen_static_attempts=8, frozen_completed_call_proxy_core_hours=sum(
            p['native_validation']['elapsed_seconds'] for p in controls['rows']) * 32 / 3600,
            atomic_completed_nodes=atomic['completed_nodes_audited'],
            atomic_new_SCFs=atomic['native_new_SCFs_reparsed'],
            atomic_cost_scope='completed nodes in supplied audit only; excludes in-flight/unaudited work',
            atomic_completed_call_proxy_core_hours=sum(
                r['sum_completed_call_seconds'] for r in atomic['candidate_rows']) * 32 / 3600,
            original_curvature_SCFs_additional=296,
            original_joint_runs_and_cached_preparation_not_recounted_here=True,
            not_total_campaign_cost_or_speedup=True),
        physical_reserved_geometry_generated_or_read=False, independent_slots_used=0,
        unique_mode_energy_TS_conditional_stability_or_predictive_advantage_certified=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--atomic-audit', type=Path, required=True)
    p.add_argument('--atomic-audit-sha256', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('fresh treatment comparison output required')
    result = analyze(a.atomic_audit, a.atomic_audit_sha256)
    with a.output.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps(result['costs']))


if __name__ == '__main__':
    main()
