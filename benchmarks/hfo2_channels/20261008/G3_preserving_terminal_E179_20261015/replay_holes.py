"""Diagnose every completed original preserving hole without repeating SCF audit.

Consumes the exact complete E163 native report. Original result/point bytes
are rechecked, then only cache-only optimizer/guard replay is performed.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import importlib.util
import json
import os
import sys

sys.dont_write_bytecode = True
RUNS = Path('/public/home/iai806/abacus/agent-runs')
EXEC = RUNS / '20261015-varneb-G3-faces-E161-r1'
PROOF_SHA = '982d2657a7a2c536509d5a8c460c1ccc0081028b367de398847720fe46937d5a'
AUDITOR_SHA = 'ab041335a93b8eaa98fa5cf8a812db5126931105cf972f586e4d15f8d06742da'
OBSERVER_SHA = 'cc745b14fe8286d61a596cd4669d2e6bae87ba8544ceb7a6dec4bd0740597cbf'
RUNTIME_SHA = 'd07bb011c9e313afe23ca4e78207356ff98adfc43ca231bdd047e6d5d0f1be20'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def module(name, path, digest):
    if sha(path) != digest:
        raise ValueError('original source bytes differ')
    spec = importlib.util.spec_from_file_location(name, path)
    out = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(out)
    return out


def validate_proof(p):
    fields = p['sacct'].split('|')
    if (p['audit_mode'] != 'complete' or p['job_id'] != '28862441'
            or p['channel'] != 'preserving' or p['completed_nodes_audited'] != 25
            or len(fields) != 7 or fields[:3] != ['28862441', 'COMPLETED', '0:0']
            or fields[4] != '32' or fields[6] != 'hfacnormal01'
            or p['runtime_manifest_sha256'] != RUNTIME_SHA or p['auditor_sha256'] != AUDITOR_SHA
            or p['pending_node'] is not None or p['completed_face_sha256'] is None
            or p['unfinished_node_failure_paths'] or p['new_DFT_or_scheduler_writes'] != 0
            or p['physical_reserved_geometry_generated_or_read'] is not False
            or p['released_stability_branch_uniqueness_contour_TS_or_forecast_certified'] is not False):
        raise ValueError('complete original native audit required')
    rows = p['candidate_rows']
    if ([r['id'] for r in rows] != list(range(1, 25))
            or len(p['point_snapshot_sha256']) != 25
            or sum(r['candidate'] is not None for r in rows) + 1 != p['gradient_candidates']
            or sum(r['completed_SCFs'] for r in rows) != p['new_SCFs_reparsed']
            or len(p['native_record_sources']) != p['new_SCFs_reparsed']):
        raise ValueError('complete source coverage/counts differ')
    return {r['id'] for r in rows if r['candidate'] is None}


def observe(proof_path):
    if any(os.environ.get(k) != '1' for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS')):
        raise ValueError('original single-thread numerical replay environment required')
    if sha(proof_path) != PROOF_SHA:
        raise ValueError('exact complete native proof required')
    proof = json.loads(proof_path.read_bytes())
    expected_holes = validate_proof(proof)
    audit = module('original_E163_complete_auditor', RUNS / '20261015-varneb-G3-face-audit-E163-r1/audit_hfo2_G3_conditional_faces.py', AUDITOR_SHA)
    diagnostic = module('unchanged_E168_observer', RUNS / '20261015-varneb-G3-bound-audit-E168-r1/diagnostics.py', OBSERVER_SHA)
    if audit.study.ROOT.resolve() != (EXEC / 'source').resolve():
        raise ValueError('original immutable runtime required')
    files = audit.worker.canary.pinned(EXEC / 'source/runtime_files.json', RUNTIME_SHA)
    if len(files) != 380 or any(sha(EXEC / 'source' / n) != h for n, h in files.items()):
        raise ValueError('original runtime drift')
    config = audit.worker.canary.pinned(audit.worker.BASE / 'G3_conditional_faces_E161_20261015/preserving_config.json', audit.CONFIGS['preserving'])
    data = audit.worker.native_preflight(config, 'preserving')
    rows = []
    for p, digest in proof['point_snapshot_sha256'].items():
        path = Path(p)
        if path.parent.resolve() != (EXEC / 'preserving').resolve() or sha(path) != digest:
            raise ValueError('original point source identity/bytes differ')
        rows.append(json.loads(path.read_bytes()))
    if [r['id'] for r in rows] != list(range(25)):
        raise ValueError('all original ordered point snapshots required')
    accepted, holes = [rows[0]], []
    proxy = SimpleNamespace(worker=audit.worker, replay_node=audit.replay_node)
    for row in rows[1:]:
        selected = audit.worker.choose_start(data, accepted, row['q_A'])
        if row['id'] not in expected_holes:
            accepted.append(row)
            continue
        if selected is None:
            if row['new_completed_SCFs'] != 0:
                raise ValueError('unsafe start cannot have native trials')
            holes.append(dict(point_id=row['id'], q_A=row['q_A'], completed_trials=0,
                              zero_DFT_unsafe_start=True, new_DFT_calls=0))
            continue
        parent, start = selected
        import numpy as np
        if parent != row['continuation_parent_id'] or not np.array_equal(start, row['starting_internal_coordinates_A']):
            raise ValueError('original continued start differs')
        records = []
        for source in [r for r in proof['native_record_sources'] if r['point_id'] == row['id']]:
            path = Path(source['validation']['raw_directory']).parent / 'result.json'
            expected = EXEC / 'preserving' / f"node_{row['id']:03d}"
            if not path.resolve().is_relative_to(expected.resolve()) or sha(path) != source['result_sha256']:
                raise ValueError('original audited native record bytes/path differ')
            records.append(json.loads(path.read_bytes()))
        if len(records) != row['new_completed_SCFs']:
            raise ValueError('all original hole trials required')
        holes.append(diagnostic.observe_replay(proxy, data, row, records, start))
    if {r['point_id'] for r in holes} != expected_holes:
        raise ValueError('all original holes required; no outcome selection')
    return dict(status='complete_original_holes_cache_only_guard_diagnosis_not_KKT_or_stability',
        original_complete_native_report_sha256=PROOF_SHA, original_job='28862441',
        processed_points=25, gradient_candidates=proof['gradient_candidates'], holes=holes,
        raw_SCFs_reparsed_by_original_complete_audit=proof['new_SCFs_reparsed'],
        new_raw_SCF_reparse_this_diagnostic=0, original_runtime_files_checked=380,
        observer_sha256=OBSERVER_SHA, original_auditor_sha256=AUDITOR_SHA,
        source_sha256=sha(Path(__file__)), new_DFT_or_scheduler_writes=0,
        optimizer_guard_force_target_or_physics_changed=False,
        physical_reserved_geometry_generated_or_read=False,
        conditional_stability_boundary_KKT_unreachable_minimum_TS_or_prediction_certified=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--proof', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('fresh external diagnostic output required')
    if not a.output.resolve().is_relative_to(RUNS / '20261015-varneb-G3-face-audit-E179-r1'):
        raise ValueError('dedicated external E179 output required')
    result = observe(a.proof)
    with a.output.open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps(dict(holes=len(result['holes']), new_DFT=0)))


if __name__ == '__main__':
    main()
