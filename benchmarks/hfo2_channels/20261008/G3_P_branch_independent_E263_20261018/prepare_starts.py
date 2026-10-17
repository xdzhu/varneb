"""Materialize training-only P8 start arrays after a pinned full native audit.

This prepares geometry only. It does not freeze predictions, authorize DFT,
read test labels, alter physical inputs or generate physical-strain holdouts.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
TRAINING = ROOT/'benchmarks/hfo2_channels/20261008/G3_preserving_branches_E260_20261017'
CONFIG_SHA = '63b1aa0eea88e7cd7c2f1a78ce9f4f3c7f39e16a9caf04d3f6d8c6ef5f84ab80'
MANIFEST_SHA = '4422902d6474afe1056181b978f6098c9fc1284ff6cc66d11ff865d499d87eff'
AUDITOR_SHA = '1f3f8fd4b880a22d2bc67743e4aee11ee2a252e5b694e1a58b408adaa47484fa'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec); spec.loader.exec_module(value)
    return value


def prepare(audit_path, expected_audit_sha256):
    raw = Path(audit_path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected_audit_sha256:
        raise ValueError('externally pinned complete native audit required')
    audit = json.loads(raw)
    if (audit['status'] != 'complete_native_ten_training_branch_trial_not_stability_certificate'
            or audit['config_sha256'] != CONFIG_SHA or audit['bundle_manifest_sha256'] != MANIFEST_SHA
            or audit['audit_LF_sha256'] != AUDITOR_SHA or audit['training_candidates_processed'] != 10
            or audit['new_independent_labels'] != 0 or audit['physical_plus005_or_private_arrays_accessed']
            or [row['id'] for row in audit['rows']] != list(range(10))
            or not audit['original_2D_refusal_preserved']):
        raise ValueError('original complete ten-training audit, not a prefix or label cohort, required')
    if sha(TRAINING/'config.json') != CONFIG_SHA or sha(TRAINING/'bundle_files.json') != MANIFEST_SHA:
        raise ValueError('original training configuration and bundle identities required')
    files = json.loads((TRAINING/'bundle_files.json').read_bytes())
    if any(sha(TRAINING/name) != digest for name, digest in files.items()):
        raise ValueError('unchanged eight-file training source required')
    worker = module('E264_original_training_geometry', TRAINING/'worker.py')
    policy = module('E264_training_only_target_policy', HERE/'target_policy.py')
    data = worker.faces.load_inputs('preserving')
    basis, axes = data['basis'], data['axes']
    def guard(x):
        try: worker.geometry(data, basis@x, axes.T@x); return True
        except ValueError: return False
    points = policy.training_only_starts(audit['rows'], axes, guard)
    for point in points:
        point['starting_full_coordinates_A'] = (
            None if point['starting_internal_coordinates_A'] is None else
            (basis@np.asarray(point['starting_internal_coordinates_A'])).tolist())
        point['independent_DFT_label'] = None
    return dict(status='training_only_P8_geometry_preparation_not_final_prediction_freeze',
        points=points, complete_native_training_audit_sha256=expected_audit_sha256,
        original_training_config_sha256=CONFIG_SHA, original_bundle_manifest_sha256=MANIFEST_SHA,
        source_LF_sha256={p.name: hashlib.sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                          for p in (Path(__file__), HERE/'target_policy.py')},
        original_forecasts_and_2D_refusal_preserved=True, revised_target_not_original_forecast_success=True,
        original_independent_inventory=16, T_tests_seen=8, remaining_P_tests=8,
        maximum_new_attempts_per_P_point=12, maximum_new_P_attempts=96, automatic_retries=0,
        new_DFT_or_scheduler_operations=0, physical_plus005_or_private_arrays_accessed=False,
        final_models_domains_and_predictions_frozen=False, authorizes_independent_DFT=False,
        stable_branches_or_TS_certified=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--audit-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError('fresh planning output required')
    result = prepare(args.audit, args.audit_sha256)
    with args.output.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(dict(status=result['status'], points=len(result['points']),
                          admissible_starts=sum(p['training_parent_id'] is not None for p in result['points']), new_DFT=0)))
