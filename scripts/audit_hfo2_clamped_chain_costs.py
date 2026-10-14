"""Cumulative costs of seven completed G2 chains, not an acceleration claim.

Reuse native-audited per-call receipts and original scheduler accounting.
The fixed ordered lineage includes capped segments; cached frames are not
new SCFs. Endpoint preparation, unfinished paths and G3 remain outside scope.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path, PurePosixPath

from examples.hfo2_fixed_input_factory import CONTRACT
from vcneb import ResponseEvaluationCost, recorded_response_dataset_cost

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'benchmarks/hfo2_channels/20261008'
CASE = BASE/'clamped_chain_costs_E147_20261014'
ACCOUNTING_SHA = 'd86642a62c6ddd7263e2c8402adacec6da0088156a2294a2b82348f896c8c931'
PILOT = ('nested_response_E055_20261010/prior_pilot_scf_costs.json',
         '20df6dcd59111049f45c2c03c305f27de4d0e80ad24a34e7efe48b22334738c3')
# Explicit original same-chain ancestry from the reviewed continuation cases.
# No sequence is chosen from its barrier, convergence speed, or unseen label.
BANDS = (
    ('PO-to-M, zero strain', ('28574709', '28661020', '28722320')),
    ('PO-to-T, zero strain', ('28692775',)),
    ('Reversing flip, zero strain', ('28692776',)),
    ('PO-to-T, +1% strain', ('28709788',)),
    ('Preserving flip, zero strain', ('28574708', '28661019', '28722810', '28723655', '28727382')),
    ('PO-to-M, +1% strain', ('28709789', '28727368')),
    ('Preserving flip, +1% strain', ('28709790', '28737465', '28812952')),
)
# job: (relative native receipt, byte SHA256, updates, fresh initial interiors)
RECEIPTS = {
    '28661020': ('clamped_M_continue_E062_20261011/audit_receipt.json',
        '479deda75b88f94e53acc5657c4c467ea7324581acbc08853f7138306f8d553f', 20, 0),
    '28661019': ('clamped_flip_continue_E067_20261011/audit_receipt.json',
        'fd0d25eef00090100291ebef6618ae3b6e2eea929020a4ff0d29e1cdc0e228a8', 20, 0),
    '28722320': ('clamped_M_terminal_E068_20261011/audit_receipt.json',
        '56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c', 13, 0),
    '28692775': ('clamped_T_terminal_E079_20261011/audit_receipt.json',
        'c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e', 10, 7),
    '28692776': ('clamped_reversing_terminal_E083_20261011/audit_receipt.json',
        '632e0e9b2b821b3bfc7d66b361c44c7db734988fb5d018ea2f31441fab14ed41', 5, 7),
    '28709788': ('clamped_plus1_T_terminal_E087_20261011/audit_receipt.json',
        '32a6e49a618e6e62632ae71a5e0ed093095254f6858221be2ade3d83cc194a28', 5, 7),
    '28722810': ('clamped_flip_continue_E069_20261011/audit_receipt.json',
        '5f078aca099f9004949eaee7dd28307adec9fdced6300b592243a7d1138802e7', 20, 0),
    '28723655': ('clamped_rolling_native_E091_20261011/28723655/terminal-r2/audit_receipt.json',
        '61edc4548d85984e78bbafc29d8b9d366157cfef48ec8175cabf3adfc7ce7e93', 40, 0),
    '28727382': ('clamped_held_flips_E100_20261011/zero_strain_preserving/audit_receipt.json',
        '46145da14f24b2b8bd46aa4656f0982650d6b2dc661282f58d809aab8a6d3202', 15, 0),
    '28709789': ('clamped_rolling_native_E091_20261011/28709789/terminal/audit_receipt.json',
        '1a05735749d8893b2c6f992f0ab90f347da9cd65b674a8c1d920227dd8565178', 10, 7),
    '28727368': ('clamped_reversing_release_E104_20261012/plus1_M_terminal/audit_receipt.json',
        'b86cdb72eab9e3d8ff2125e16e6fbb2d0954b34d54c0f05099a581c066108389', 22, 0),
    '28709790': ('clamped_flip_continue_E103_20261012/plus1_preserving_starter/audit_receipt.json',
        '8eaac3f12f513caaa24fc6dac72ba124e3b42f00646b2d31bad8c2c5acde1d15', 10, 7),
    '28737465': ('clamped_preserving_cap_E140_20261012/completed_HF/audit_receipt.json',
        'df34573384f76650e123691dbb508a52d7d8ee529e6e50dca0fce5e4e6835c02', 40, 0),
    '28812952': ('clamped_plus1_preserving_terminal_E145_20261012/completed_HF/audit_receipt.json',
        '0c8d91532bc9051decdfd600eae12ce409fd84a86c677158077948fa352c9d8a', 17, 0),
}
JOBS = {job for _, jobs in BANDS for job in jobs}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pinned(path, digest):
    if sha(path) != digest:
        raise ValueError('original receipt byte hash differs')
    return json.loads(path.read_text())


def parse_accounting(body, jobs=JOBS):
    rows = {}
    for line in body.splitlines():
        fields = line.strip().split('|')
        if len(fields) != 11:
            raise ValueError('complete allocation-level accounting row required')
        job, state, exitcode, cpus, nodes, partition, user, seconds, start, end, cputime = fields
        if (job not in jobs or job in rows or state != 'COMPLETED' or exitcode != '0:0'
                or cpus != '32' or nodes != '1' or partition != 'hfacnormal01' or user != 'iai806'):
            raise ValueError('exact successful registered single-node32 allocation required')
        elapsed = int(seconds)
        if (elapsed <= 0 or int(cputime) != elapsed*32
                or (datetime.fromisoformat(end)-datetime.fromisoformat(start)).total_seconds() != elapsed):
            raise ValueError('native duration, start/end and allocated core-seconds disagree')
        rows[job] = dict(start=start, end=end, elapsed_seconds=elapsed,
                         allocation_core_hours=elapsed*32/3600)
    if set(rows) != set(jobs):
        raise ValueError('missing registered production segment')
    return rows


def validate_call_grid(records, per_image):
    actual = {(PurePosixPath(r.evaluation_id.removeprefix('hf:')).parent.name,
               PurePosixPath(r.evaluation_id.removeprefix('hf:')).name) for r in records}
    expected = {(f'image_{i:04d}', f'scf_{j:06d}') for i in range(1, 8) for j in range(per_image)}
    if len(records) != len(actual) or actual != expected:
        raise ValueError('complete exact seven-interior fresh-call grid required')


def call_records(rows):
    """Exact source call identity, never geometry/hash-based deduplication."""
    result = []
    for row in rows:
        if {k: row['input_sha256'][k] for k in CONTRACT} != CONTRACT:
            raise ValueError('recorded original physical contract differs')
        path = PurePosixPath(row['path'])
        if path.name != 'call_audit.json' or not path.parent.name.startswith('scf_'):
            raise ValueError('original explicit SCF audit source required')
        result.append(ResponseEvaluationCost('hf:'+str(path.parent), row['audit_sha256'],
                      'completed', row['elapsed_seconds'], 32))
    if len({r.evaluation_id for r in result}) != len(result):
        raise ValueError('duplicate SCF source inside a production segment')
    return result


def segment_cost(records, allocated_core_hours):
    cost = recorded_response_dataset_cost(records, required_evaluation_ids=[r.evaluation_id for r in records])
    if (not cost['cost_complete'] or cost['failed_evaluations']
            or cost['recorded_transport_cpu_core_hours_sum'] > allocated_core_hours+1e-8):
        raise ValueError('complete successful serial SCF timings within the native allocation required')
    return cost


def audit():
    accounting = pinned(CASE/'native_accounting.json', ACCOUNTING_SHA)
    native = parse_accounting(accounting['stdout'])
    pilot = pinned(BASE/PILOT[0], PILOT[1])
    if pilot['original_physical_contract_sha256'] != CONTRACT:
        raise ValueError('pilot original physical contract differs')
    proofs = {p['evaluation_id']: p for p in pilot['per_call_proofs']}
    grouped = defaultdict(list)
    for row in pilot['records']:
        proof = proofs[row['evaluation_id']]
        if (proof['native_MPI_ranks'] != 32 or not proof['native_complete_EFS_checked']
                or not proof['six_original_physical_inputs_checked']):
            raise ValueError('pilot native EFS/physical/MPI audit missing')
        grouped[proof['job_id']].append(ResponseEvaluationCost(**row))
    if set(grouped) != {'28574708', '28574709'} or any(len(r) != 77 for r in grouped.values()):
        raise ValueError('both original eleven-frame seven-interior pilots required')
    for records in grouped.values():
        validate_call_grid(records, 11)
    sources = {job: dict(receipt=PILOT[0], sha256=PILOT[1], updates=10, fresh_initial_interiors=7)
               for job in grouped}
    for job, (path, pin, updates, initial) in RECEIPTS.items():
        r = pinned(BASE/path, pin)
        if (r.get('job_id', r.get('parent_job')) != job
                or r['scheduler_row'] != job+'|COMPLETED|0:0|32|hfacnormal01|iai806'
                or r['new_DFT_calls'] != 0 or r['holdout_generated_or_read'] is not False):
            raise ValueError('cost receipt identity or scope differs')
        full = 'fresh_interior_SCFs' in r
        check = ('all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked' if full else
                 'all_parent_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked')
        count = r['fresh_interior_SCFs' if full else 'fresh_parent_interior_SCFs']
        if (not r[check] or count != 7*updates+initial or len(r['calls']) != count
                or (full and r['terminal_step'] != updates)):
            raise ValueError('full registered fresh-SCF segment history required')
        records = call_records(r['calls'])
        validate_call_grid(records, updates+initial//7)
        cost = segment_cost(records, native[job]['allocation_core_hours'])
        saved = r['fresh_SCF_core_hour_proxy_32' if full else 'parent_SCF_core_hours_32']
        if not math.isclose(saved, cost['recorded_transport_cpu_core_hours_sum'], abs_tol=1e-9, rel_tol=0):
            raise ValueError('per-call timing sum differs from original native audit')
        grouped[job] = records
        sources[job] = dict(receipt=path, sha256=pin, updates=updates, fresh_initial_interiors=initial)
    bands, all_records = [], []
    for label, jobs in BANDS:
        final = pinned(BASE/RECEIPTS[jobs[-1]][0], RECEIPTS[jobs[-1]][1])
        if final['ordinary_converged'] is not True:
            raise ValueError('only actual ordinary-converged terminal bands admitted')
        if any(native[a]['end'] > native[b]['start'] for a, b in zip(jobs[:-1], jobs[1:])):
            raise ValueError('same-chain registered segments overlap in native time')
        records = [r for job in jobs for r in grouped[job]]
        allocation = math.fsum(native[j]['allocation_core_hours'] for j in jobs)
        cost = segment_cost(records, allocation)
        if cost['unique_source_DFT_evaluations'] != len(records):
            raise ValueError('registered fresh segment histories unexpectedly reuse one source call')
        bands.append(dict(label=label, jobs=list(jobs), segments=len(jobs),
            optimizer_updates=sum(sources[j]['updates'] for j in jobs),
            fresh_initial_interiors=sum(sources[j]['fresh_initial_interiors'] for j in jobs),
            fresh_interior_SCFs=len(records), SCF_core_hour_proxy=cost['recorded_transport_cpu_core_hours_sum'],
            allocation_core_hours=allocation, last_segment_fresh_SCFs=len(grouped[jobs[-1]]),
            last_segment_SCF_core_hour_proxy=segment_cost(grouped[jobs[-1]], native[jobs[-1]]['allocation_core_hours'])['recorded_transport_cpu_core_hours_sum']))
        all_records.extend(records)
    total = recorded_response_dataset_cost(all_records, required_evaluation_ids=[r.evaluation_id for r in all_records])
    if total['unique_source_DFT_evaluations'] != len(all_records):
        raise ValueError('different production chains share a supposedly fresh SCF source')
    return dict(status='seven_completed_production_chain_costs_not_full_study_or_acceleration',
        source_script_sha256=sha(Path(__file__)), accounting_sha256=ACCOUNTING_SHA, sources=sources,
        bands=bands, total=total, allocation_core_hours=math.fsum(r['allocation_core_hours'] for r in native.values()),
        cumulative_production_segments=len(native), native_accounting=native,
        pinned_prior_native_audits_reused=True, all_native_assets_or_raw_EFS_reaudited_now=False,
        all_preparation_or_full_study_cost_complete=False,
        excluded_costs=['endpoint relaxation/cache preparation', 'G1 and starting-structure preparation',
            'reference phonons/polarization and earlier diagnostics', 'unfinished +1% reversing chain',
            'future G3 surfaces, failed starts, probes and independent predictions'],
        physical_parameters_changed=False, G3_or_holdout_access=False, new_DFT_calls=0, scheduler_mutations=0,
        acceleration_TS_physical_error_or_prediction_accuracy_certified=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh cost-audit output required')
    result = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(dict(bands=len(result['bands']), fresh_SCFs=result['total']['unique_source_DFT_evaluations'],
        SCF_core_hours=result['total']['recorded_transport_cpu_core_hours_sum'],
        allocation_core_hours=result['allocation_core_hours'], new_DFT_calls=0)))


if __name__ == '__main__':
    main()
