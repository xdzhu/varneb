"""Offline ordered replay of the complete natively audited training cohort.

Consumes stored physical energy/gradients only; no calculator or scheduler.
Native provenance belongs to the separately pinned original complete auditor.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import numpy as np

HERE = Path(__file__).resolve().parent
AUDIT_SHA = '30971f61ab7d12d8e9bfa729a517603566c20cefa6bf442b06a4100264bd8811'


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec); spec.loader.exec_module(value)
    return value


def analyze():
    audit_path = HERE/'native_complete_E264.json'
    if hashlib.sha256(audit_path.read_bytes()).hexdigest() != AUDIT_SHA:
        raise ValueError('exact complete native audit required')
    audit = json.loads(audit_path.read_bytes())
    worker = module('E264_original_worker_offline', HERE/'worker.py')
    cached = module('E264_ordered_replay', HERE/'prefix_review_E262.py')
    decomposition = module('E264_radial_decomposition', HERE/'diagnose_prefix_E262.py')
    config = json.loads((HERE/'config.json').read_bytes())
    data, _, _, free, seed, _ = worker.validate(config, HERE)
    direction = np.asarray(config['branch_direction_36'])
    basis = data['basis']; origin = basis.T@np.asarray(data['origin_record']['coordinates'])
    by_id = {r['id']: r for r in audit['rows']}
    if sorted(by_id) != list(range(10)) or audit['new_native_SCFs_reparsed'] != 81:
        raise ValueError('all ten and all81 audited native calls required')
    cursor = 0; groups = {}
    for index in (0, 2, 4, 6, 8, 1, 3, 5, 7, 9):
        row = by_id[index]; count = row['attempted_SCFs']
        records = audit['native_records'][cursor:cursor+count]; cursor += count
        if len(records) != count or len(row['native_sources']) != count:
            raise ValueError('every charged native record and source required')
        for record, source in zip(records, row['native_sources']):
            if Path(record['validation']['raw_directory']) != Path(source['result_path']).parent/'scf_000000':
                raise ValueError('ordered exact native source identity differs')
        groups[index] = records
    if cursor != len(audit['native_records']): raise ValueError('unconsumed native records')
    rows = []
    for point in config['points']:
        row = by_id[point['id']]; records = groups[point['id']]
        def guard(x):
            try: worker.geometry(data, basis@x, point['q_A']); return True
            except ValueError: return False
        recovered = cached.replay(worker.solver, point, basis, free, direction, seed, records, guard)
        if recovered['status'] == 'budget_exception_replayed':
            x = np.asarray(recovered['last_accepted_internal_A'])
            energy = recovered['last_accepted_energy_eV_cell']; trace = recovered['accepted_trace']
            if row['label'] is not None: raise ValueError('capped replay cannot replace a returned candidate')
        else:
            final = recovered['result']; x = np.asarray(final['coordinates'])
            energy = final['energy']; trace = final['trace']
            if row['label'] is None or energy != row['label']['energy_eV_cell']:
                raise ValueError('same returned native energy required')
        matches = [r for r in records if r['enthalpy_eV'] == energy
                   and np.allclose(r['coordinates'], basis@x, atol=2e-12, rtol=0)]
        if len(matches) != 1: raise ValueError('unique audited last accepted state required')
        native = matches[0]; gradient = basis.T@np.asarray(native['gradient'])
        rows.append(dict(id=row['id'], q_A=row['q_A'], start_sign=row['start_sign'],
            status=row['status'], native_SCFs=len(records), replay_status=recovered['status'],
            accepted_states=len(trace), energy_rejected_native_calls=len(records)-len(trace),
            maximum_replay_coordinate_arithmetic_difference_A=recovered['maximum_coordinate_arithmetic_difference_A'],
            total_radius_A=float(np.linalg.norm(x)), radius_margin_A=.5-float(np.linalg.norm(x)),
            minimum_distance_A=native['minimum_distance_A'], energy_eV_cell=energy,
            zeta_A=float(direction@(x-origin)), zeta_gradient_eV_A=float(direction@gradient),
            **decomposition.radial_gradient(x, gradient, free)))
    return dict(status='complete_training_ordered_offline_replay_not_second_native_audit',
        complete_native_audit_sha256=AUDIT_SHA,
        source_LF_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
        rows=rows, pairs=audit['pairs'], native_SCFs_consumed=sum(r['native_SCFs'] for r in rows),
        gradient_candidates=audit['gradient_candidates'], completed_call_proxy_core_hours=audit['completed_call_proxy_core_hours'],
        allocation_core_hours=audit['allocation_core_hours'], all_prior_preparation_costs_additional=True,
        stable_branches_or_TS_or_predictive_gain_certified=False, new_DFT_or_scheduler_operations=0,
        physical_plus005_or_private_arrays_accessed=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): raise FileExistsError('fresh complete-cohort analysis output required')
    result = analyze()
    with args.output.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(dict(status=result['status'], native_SCFs_consumed=result['native_SCFs_consumed'],
        rows=[{k: r[k] for k in ('id', 'status', 'accepted_states', 'energy_rejected_native_calls',
            'released_gradient_norm_eV_A', 'radius_margin_A')} for r in result['rows']]), indent=2))
