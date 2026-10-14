"""Read-only native audit of completed E161 candidate nodes, not a PES certificate.

Run this standalone auditor with the immutable E161 source on PYTHONPATH.
Prefix mode never reads the current unfinished node as an accepted result;
complete mode additionally requires the original normally completed job and face.
No calculator, dispatcher, source repair or restart is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator

from scripts import hfo2_G3_conditional_faces as worker
from scripts.audit_hfo2_G3_conditional_canary import CacheOnly
from vcneb.mode_evaluator import _fingerprint

study = worker.study
ROOT = study.RUNS / '20261015-varneb-G3-faces-E161-r1'
RUNTIME_SHA = 'd07bb011c9e313afe23ca4e78207356ff98adfc43ca231bdd047e6d5d0f1be20'
FIRST_AUDIT_SHA = '4bfa088f920980a4582bb6e1dedf5d370b898ca81a4ac75a7fc57122cc2dfc58'
CONFIGS = {'T': '44f3b34e77a307b5d965893309a11b9ebfd1ad8961b3035904544843c6ec94d6',
           'preserving': 'b50611685e045023eb1bcd533232da144cad0ca0da67a1b02cb1e5fe0ac5b5cd'}
JOBS = {'T': '28862319', 'preserving': '28862441'}


def load_prefix(directory, specs):
    """Snapshot a consecutive completed prefix; partial JSON is pending, not failure."""
    rows, hashes = [], {}
    pending = None
    for spec in specs:
        path = directory / f"point_{spec['id']:03d}.json"
        if not path.exists():
            pending = dict(point_id=spec['id'], reason='point_file_not_yet_present')
            break
        raw = path.read_bytes()
        try:
            row = json.loads(raw)
        except (json.JSONDecodeError, UnicodeDecodeError):
            pending = dict(point_id=spec['id'], reason='point_JSON_not_yet_complete')
            break
        if any(json.dumps(row.get(k), allow_nan=False) != json.dumps(v) for k, v in spec.items()):
            raise ValueError('point identity or fixed order differs')
        rows.append(row)
        hashes[str(path)] = hashlib.sha256(raw).hexdigest()
    return rows, hashes, pending


def node_diagnostics(data, row, records):
    """Check fixed Q, exact native final identity and released residual; no stability claim."""
    basis, axes = data['basis'], data['axes']
    q = np.asarray(row['q_A'])
    n = len(records)
    for key in ('new_completed_SCFs', 'reserved_attempts'):
        if type(row[key]) is not int or row[key] != n or not 0 <= n <= 12:
            raise ValueError('completed-node attempt accounting differs')
    for record in records:
        x, g = np.asarray(record['coordinates']), np.asarray(record['gradient'])
        if (record['status'] != 'complete' or x.shape != (40,) or g.shape != (40,)
                or not np.isfinite(x).all() or not np.isfinite(g).all() or x[-1] != 0
                or not np.isfinite(record['enthalpy_eV'])
                or not np.allclose(basis @ (basis.T @ x), x, atol=2e-12, rtol=0)
                or not np.allclose(axes.T @ (basis.T @ x), q, atol=2e-12, rtol=0)):
            raise ValueError('trial geometry/gradient/fixed Q or prescribed control differs')
    times = [r['validation']['elapsed_seconds'] for r in records]
    if any(isinstance(t, bool) or not np.isfinite(t) or t <= 0 for t in times):
        raise ValueError('positive native call seconds required')
    if not np.isclose(sum(times), row['sum_completed_native_call_seconds'], atol=1e-10, rtol=0):
        raise ValueError('native call cost differs')
    point = row['point']
    if point is None:
        if row['status'] != 'optimizer_not_converged' or row['reason'] not in (
                'mode evaluation attempt budget exhausted; calculator was not called',
                'no starting point reached the orthogonal-gradient tolerance'):
            raise ValueError('missing candidate is not a registered nonconvergence hole')
        if row['reason'].startswith('mode evaluation attempt budget exhausted') and n != 12:
            raise ValueError('budget hole does not exhaust its actual budget')
        return dict(candidate=None, completed_SCFs=n, call_seconds=sum(times))
    if row['status'] != 'gradient_candidate_release_stability_unchecked':
        raise ValueError('candidate status differs; no last-trial promotion')
    internal, full = np.asarray(point['coordinates_internal_A']), np.asarray(point['coordinates_full_A'])
    if (internal.shape != (36,) or full.shape != (40,) or not np.isfinite(internal).all()
            or not np.allclose(basis @ internal, full, atol=2e-12, rtol=0)
            or not np.allclose(axes.T @ internal, q, atol=2e-12, rtol=0)):
        raise ValueError('final candidate embedding or retained coordinates differ')
    matches = [r for r in records if np.array_equal(r['coordinates'], full)]
    if len(matches) != 1:
        raise ValueError('final candidate must match exactly one native result')
    final = matches[0]
    g = basis.T @ np.asarray(final['gradient'])
    residual = float(np.linalg.norm(g - axes @ (axes.T @ g)))
    if (point['energy_eV_cell'] != final['enthalpy_eV'] or residual > .03
            or not np.isclose(residual, point['orthogonal_gradient_norm_eV_A'], atol=2e-12, rtol=0)
            or not np.array_equal(row['physical_gradient_eV_A'], final['gradient'])
            or row['native_source'] != final['validation']
            or type(row['callback_requests']) is not int or type(row['center_cache_hits']) is not int
            or row['center_cache_hits'] < 0
            or row['callback_requests'] < n + row['center_cache_hits']
            or point['callback_evaluations'] != row['callback_requests']):
        raise ValueError('final native energy/gradient/common target or counters differ')
    return dict(candidate=dict(energy_eV_cell=final['enthalpy_eV'],
        released_gradient_norm_eV_A=residual, coordinates_internal_A=internal.tolist(),
        coordinates_full_A=full.tolist(), full_physical_gradient_eV_A=final['gradient'],
        retained_gradient_components_eV_A=(axes.T @ g).tolist(),
        released_stability_unchecked=True), completed_SCFs=n, call_seconds=sum(times))


def replay_node(data, row, records, start):
    """Strict original-order optimizer replay using only existing native E/g records."""
    oracle = CacheOnly(records)
    evaluate = worker.canary.InternalEvaluator(oracle, data['basis'], data['center_energy'], data['center_gradient'])
    plane = worker.canary.ModePlane(np.zeros(36), data['axes'], np.ones(36), np.eye(2),
        ('coupled internal1', 'coupled internal2'), 'local length', 'G3-first-conditional-sheet')
    screens = dict(checks=0, rejected=0)
    def guard(x):
        screens['checks'] += 1
        ok = worker.canary.geometry_guard(data['chart'], data['basis'], x)
        screens['rejected'] += int(not ok)
        return ok
    point, reason = None, None
    try:
        point = worker.canary.relax_orthogonal_at_q(plane, row['q_A'], evaluate, starts=[start],
            include_frozen_start=False, optimizer='safeguarded_bfgs', trial_validator=guard,
            initial_hessian=data['hessian'], initial_trust_radius=.05, gradient_tolerance=.03, max_iterations=200)
    except RuntimeError as exc:
        reason = str(exc)
        if reason not in ('mode evaluation attempt budget exhausted; calculator was not called',
                          'no starting point reached the orthogonal-gradient tolerance'):
            raise
    if (len(oracle.seen) != len(records) or evaluate.callback_requests != row['callback_requests']
            or evaluate.center_cache_hits != row['center_cache_hits'] or screens != row['geometry_screens']
            or (point is None) != (row['point'] is None)):
        raise ValueError('node replay decisions/counters differ')
    if point is None:
        if reason != row['reason']:
            raise ValueError('node replay stop reason differs')
    elif (not np.array_equal(point.coordinates, row['point']['coordinates_internal_A'])
            or point.energy != row['point']['energy_eV_cell']
            or not np.isclose(point.orthogonal_gradient_norm,
                              row['point']['orthogonal_gradient_norm_eV_A'], atol=2e-12, rtol=0)):
        raise ValueError('node replay final candidate differs')
    return dict(original_order_cache_only_replay_passed=True, new_DFT=0)


def native_records(data, directory):
    """Re-read E/F/stress, original six assets, reservation identities and exact contract."""
    chart, basis = data['chart'], data['basis']
    ref = chart.reference_atoms
    expected = dict(schema_version=1, coordinate_convention=chart.coordinate_convention,
        reference=dict(numbers=ref.numbers.tolist(), cell_A=ref.cell.array.tolist(),
            scaled_positions=ref.get_scaled_positions(wrap=False).tolist(), pbc=ref.pbc.tolist()),
        calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
        validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1', pressure_eV_per_A3=0.,
        minimum_distance_A=1.6, mechanical_boundary=chart.mechanical_contract)
    fingerprint = _fingerprint(expected)
    contract = json.loads((directory / 'contract.json').read_bytes())
    if (contract != dict(expected, contract_sha256=fingerprint)
            or _fingerprint({k:v for k,v in contract.items() if k != 'contract_sha256'}) != fingerprint):
        raise ValueError('native node contract differs')
    budget = dict(schema_version=1, max_evaluation_attempts=12, physical_contract_sha256=fingerprint)
    if json.loads((directory / 'evaluation-budget.json').read_bytes()) != budget:
        raise ValueError('native per-node budget differs')
    paths = sorted(directory.glob('eval-*/result.json'))
    markers = sorted(directory.glob('attempt-*.json'))
    if (list(directory.glob('eval-*/failure.json')) or len(list(directory.glob('eval-*'))) != len(paths)
            or len(markers) != len(paths) or len(paths) > 12):
        raise ValueError('completed node contains failed/incomplete or unreserved trials')
    records, sources = [], []
    for index, (path, marker) in enumerate(zip(paths, markers)):
        row = json.loads(path.read_bytes())
        x = np.asarray(row['coordinates'])
        digest = hashlib.sha256(x.tobytes()).hexdigest()
        reservation = dict(schema_version=1, evaluation_index=index, coordinate_sha256=digest,
                           coordinates=x.tolist(), budget_sha256=_fingerprint(budget))
        if (path.parent.name != f'eval-{index:06d}-{digest[:12]}' or marker.name != f'attempt-{index:06d}.json'
                or json.loads(marker.read_bytes()) != reservation
                or row['contract_sha256'] != fingerprint or row['coordinate_sha256'] != digest):
            raise ValueError('native record/reservation identity differs')
        atoms = chart.to_atoms(x)
        atoms.calc = SinglePointCalculator(atoms, **study.audited_results(path.parent / 'scf_000000'))
        checked = study.result_validator(index, atoms, path.parent)
        e, g = chart.energy_and_gradient(x, atoms)
        if (checked != row['validation'] or e != row['enthalpy_eV']
                or not np.allclose(g, row['gradient'], atol=2e-12, rtol=0)
                or not np.array_equal(atoms.get_stress(voigt=False), row['stress_eV_per_A3'])
                or not np.array_equal(atoms.cell.array[:2], chart.reference_atoms.cell.array[:2])
                or not worker.canary.geometry_guard(chart, basis, basis.T @ x)):
            raise ValueError('native EFS/gradient/fixed substrate/geometry replay differs')
        records.append(row)
        sources.append(dict(result_sha256=study.sha256(path), validation=checked))
    return records, sources


def audit(channel, mode):
    if study.ROOT.resolve() != (ROOT / 'source').resolve():
        raise ValueError('import original immutable E161 source; never overwrite it with the auditor')
    files = worker.canary.pinned(study.ROOT / 'runtime_files.json', RUNTIME_SHA)
    if len(files) != 380 or any(study.sha256(study.ROOT / n) != h for n,h in files.items()):
        raise ValueError('original runtime Python bytes differ')
    config = worker.canary.pinned(worker.BASE / 'G3_conditional_faces_E161_20261015' / f'{channel}_config.json', CONFIGS[channel])
    data = worker.native_preflight(config, channel)
    first = worker.canary.pinned(ROOT / 'first_native_SCF_audit.json', FIRST_AUDIT_SHA)
    receipt = next(r for r in first['rows'] if r['channel'] == channel)['receipts']['accepted']
    if json.loads((ROOT / f'{channel}_accepted.json').read_bytes()) != receipt or receipt['job_id'] != JOBS[channel]:
        raise ValueError('original accepted job identity differs')
    accounting = subprocess.run(['sacct', '-n', '-P', '-X', '-j', JOBS[channel],
        '--format=JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,NodeList,Partition'],
        check=True, capture_output=True, text=True).stdout.strip()
    fields = accounting.split('|')
    if len(fields) < 7 or fields[0] != JOBS[channel] or fields[4] != '32' or fields[6] != 'hfacnormal01':
        raise RuntimeError('original accounting handle/resource identity unavailable: ' + accounting)
    directory = ROOT / channel
    rows, hashes, pending = load_prefix(directory, config['points'])
    if not rows:
        raise RuntimeError('no complete node prefix yet; no failure or restart inferred')
    origin = dict(**config['points'][0], status='cached_origin_gradient_candidate_curvature_unchecked',
        point=data['origin_point'], new_completed_SCFs=0, reserved_attempts=0,
        physical_gradient_eV_A=data['origin_record']['gradient'], origin_native_source=data['origin_record']['validation'])
    if json.dumps(rows[0], sort_keys=True, allow_nan=False) != json.dumps(origin, sort_keys=True, allow_nan=False):
        raise ValueError('cached origin differs; no repeated origin SCF is allowed')
    accepted, audited, sources = [rows[0]], [], []
    for row in rows[1:]:
        selected = worker.choose_start(data, accepted, row['q_A'])
        node = directory / f"node_{row['id']:03d}"
        if selected is None:
            if (row['status'] != 'no_safe_continuation_start_no_DFT' or row['point'] is not None
                    or row['new_completed_SCFs'] != 0 or row['reserved_attempts'] != 0 or node.exists()):
                raise ValueError('unsafe continuation must remain a zero-DFT hole')
            audited.append(dict(id=row['id'], q_A=row['q_A'], candidate=None, completed_SCFs=0, call_seconds=0.))
            continue
        parent, start = selected
        if row['continuation_parent_id'] != parent or not np.array_equal(start, row['starting_internal_coordinates_A']):
            raise ValueError('nearest-safe continuation parent/start differs')
        records, native = native_records(data, node)
        result = node_diagnostics(data, row, records)
        result.update(id=row['id'], q_A=row['q_A'], continuation_parent_id=parent,
                      replay=replay_node(data, row, records, start))
        audited.append(result)
        sources.extend(dict(point_id=row['id'], **s) for s in native)
        if row['point'] is not None:
            accepted.append(row)
    attempts = sum(r['reserved_attempts'] for r in rows)
    if attempts > 288 or sum(r['call_seconds'] for r in audited) > int(fields[3]) + 1:
        raise ValueError('prefix attempts or sequential native cost exceeds original allocation')
    complete_sha = None
    if mode == 'complete':
        if fields[1:3] != ['COMPLETED', '0:0'] or len(rows) != 25 or pending is not None:
            raise RuntimeError('original face not normally complete; do not restart: ' + accounting)
        path = directory / 'face.json'
        face = json.loads(path.read_bytes())
        if (face['status'] != 'completed_first_5x5_candidate_sheet_not_certified_PES'
                or face['channel'] != channel or face['policy'] != worker.POLICY or face['rows'] != rows
                or face['gradient_candidates'] != len(accepted) or face['new_reserved_attempts'] != attempts
                or face['new_completed_SCFs'] != len(sources) or face['cached_origin_new_SCFs'] != 0
                or face['released_stability_global_branch_TS_contour_or_prediction_certified'] is not False
                or face['physical_reserved_geometry_generated_or_read'] is not False):
            raise ValueError('completed face summary differs from actual nodes')
        complete_sha = study.sha256(path)
    if any(study.sha256(Path(p)) != h for p,h in hashes.items()):
        raise ValueError('node snapshot changed during audit; re-observe, never repair')
    return dict(status='native_completed_node_prefix_audited_not_PES_TS_or_prediction', audit_mode=mode,
        channel=channel, job_id=JOBS[channel], sacct=accounting, completed_nodes_audited=len(rows),
        cached_origin_new_SCFs=0, gradient_candidates=len(accepted), new_SCFs_reparsed=len(sources),
        reserved_attempts_in_audited_prefix=attempts, candidate_rows=audited, native_record_sources=sources,
        point_snapshot_sha256=hashes, pending_node=pending, completed_face_sha256=complete_sha,
        unfinished_node_failure_paths=[str(p) for p in directory.glob('node_*/eval-*/failure.json')],
        completed_call_proxy_core_hours=sum(r['call_seconds'] for r in audited)*32/3600,
        allocation_core_hours_at_observation=int(fields[3])*32/3600,
        preparation_costs_additional=True, runtime_manifest_sha256=RUNTIME_SHA,
        auditor_sha256=study.sha256(Path(__file__)),
        released_stability_branch_uniqueness_contour_TS_or_forecast_certified=False,
        physical_reserved_geometry_generated_or_read=False, new_DFT_or_scheduler_writes=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--channel', choices=tuple(CONFIGS), required=True)
    parser.add_argument('--mode', choices=('prefix', 'complete'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh external audit output required')
    if any(args.output.resolve().is_relative_to(p.resolve()) for p in (ROOT, study.ASSETS)):
        raise ValueError('audit output must remain outside live study/source/assets')
    result = audit(args.channel, args.mode)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    study.exclusive_json(args.output, result)
    print(json.dumps({k:result[k] for k in ('status', 'channel', 'completed_nodes_audited', 'new_SCFs_reparsed')}))


if __name__ == '__main__':
    main()
