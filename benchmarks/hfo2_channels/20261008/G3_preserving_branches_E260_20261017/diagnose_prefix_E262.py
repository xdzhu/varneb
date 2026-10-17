"""Read-only geometry/gradient diagnosis of the natively audited first four.

No calculator, scheduler, interpolation, changed domain or convergence test.
The original guard bounds ||x||, not displacement from the relaxed origin.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
PREFIX_SHA = 'c122d71eb2dd46f3a3e892e49389b5d38a24a38eb86794ee0dd35ea7bc100772'
CONFIG_SHA = '63b1aa0eea88e7cd7c2f1a78ce9f4f3c7f39e16a9caf04d3f6d8c6ef5f84ab80'


def radial_gradient(x, gradient, free):
    """Decompose the physical released gradient about the declared zero center."""
    y = free.T @ x
    g = free.T @ gradient
    radius = float(np.linalg.norm(y))
    radial = float(g @ y / radius) if radius else 0.
    tangent = g - radial * y / radius if radius else g
    return dict(released_radius_A=radius,
                frozen_radius_A=float(np.linalg.norm(x - free @ y)),
                released_gradient_norm_eV_A=float(np.linalg.norm(g)),
                radial_gradient_eV_A=radial,
                tangent_gradient_norm_eV_A=float(np.linalg.norm(tangent)))


def pinned(path, digest):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError('exact audited source required')
    return json.loads(raw)


def diagnose():
    report = pinned(HERE/'prefix_review_E262.json', PREFIX_SHA)
    config = pinned(HERE/'config.json', CONFIG_SHA)
    if (report['scope_ids'] != [0, 1, 2, 3]
            or report['native_SCFs_reparsed'] != 36
            or report['complete_ten_candidate_audit']):
        raise ValueError('first-four native audit scope required')
    basis = np.asarray(report['basis_40_by_36'])
    free = np.asarray(report['free_basis_36_by_34'])
    origin = np.asarray(report['origin_internal_A'])
    direction = np.asarray(config['branch_direction_36'])
    rows = []
    for row in report['rows']:
        replay = row['replay']
        if replay['status'] == 'budget_exception_replayed':
            x = np.asarray(replay['last_accepted_internal_A'])
            energy = replay['last_accepted_energy_eV_cell']
            trace = replay['accepted_trace']
        else:
            x = np.asarray(replay['result']['coordinates'])
            energy = replay['result']['energy']
            trace = replay['result']['trace']
        # Identity is the exact accepted energy/index in an already ordered
        # native replay, never a nearest-coordinate cache search. Matrix
        # reconstruction on a different BLAS uses the existing arithmetic check.
        matches = [i for i, r in enumerate(row['native_records']) if energy == r['enthalpy_eV']]
        if len(matches) != 1:
            raise ValueError('unique native last accepted state required')
        final = row['native_records'][matches[0]]
        arithmetic_error = float(np.max(np.abs(basis @ x - final['coordinates'])))
        if not np.allclose(basis @ x, final['coordinates'], atol=2e-12, rtol=0):
            raise ValueError('accepted-state matrix reconstruction differs')
        g = basis.T @ np.asarray(final['gradient'])
        accepted = []
        cursor = 0
        for state in trace:
            indices = [i for i, r in enumerate(row['native_records'])
                       if i >= cursor and state['energy'] == r['enthalpy_eV']]
            if len(indices) != 1:
                raise ValueError('unambiguous ordered accepted energy required')
            accepted.append(indices[0]); cursor = indices[0] + 1
        if accepted[-1] != matches[0] or any(
                b['energy'] > a['energy'] for a, b in zip(trace, trace[1:])):
            raise ValueError('last accepted state or descent sequence differs')
        total = float(np.linalg.norm(x))
        if total > .5 or min(r['minimum_distance_A'] for r in row['native_records']) < 1.6:
            raise ValueError('original native radius/distance guard differs')
        result = dict(id=row['id'], q_A=row['q_A'], start_sign=row['start_sign'],
            producer_status=row['producer_status'], attempted_native_SCFs=len(row['native_records']),
            accepted_native_indices=accepted,
            energy_rejected_native_SCFs=len(row['native_records'])-len(accepted),
            last_accepted_native_index=matches[0], total_radius_A=total,
            accepted_state_reconstruction_max_difference_A=arithmetic_error,
            total_radius_margin_A=.5-total,
            minimum_distance_at_last_accepted_A=final['minimum_distance_A'],
            energy_drop_from_start_eV_cell=energy-row['native_records'][0]['enthalpy_eV'],
            zeta_A=float(direction @ (x-origin)), branch_gradient_eV_A=float(direction @ g),
            **radial_gradient(x, g, free))
        if 'uncomputed_next_proposal_internal_A' in replay:
            proposal = np.asarray(replay['uncomputed_next_proposal_internal_A'])
            result['uncomputed_next_proposal_radius_A'] = float(np.linalg.norm(proposal))
            result['uncomputed_next_proposal_is_not_DFT'] = True
        rows.append(result)
    successful = [r['replay']['result'] for r in report['rows'][2:]]
    pair = dict(q_A=[-.1, 0.],
                energy_difference_meV_cell=abs(successful[0]['energy']-successful[1]['energy'])*1000,
                internal_distance_A=float(np.linalg.norm(np.subtract(
                    successful[0]['coordinates'], successful[1]['coordinates']))),
                distinct_stable_branches_established=False)
    seconds = [r['validation']['elapsed_seconds'] for row in report['rows'] for r in row['native_records']]
    return dict(status='audited_prefix_radius_and_accepted_state_diagnosis',
                prefix_sha256=PREFIX_SHA, config_sha256=CONFIG_SHA,
                source_LF_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
                source_observed_CST=report['observed_CST'], scope_ids=[0, 1, 2, 3], rows=rows,
                signed_start_pair=pair, native_calls=36,
                native_call_elapsed_seconds=dict(min=min(seconds), median=float(np.median(seconds)), max=max(seconds)),
                completed_call_proxy_core_hours=sum(seconds)*32/3600,
                original_radius_center='internal coordinate zero, not relaxed origin',
                new_DFT_or_scheduler_writes=0, physical_inputs_or_domain_changed=False,
                complete_cohort_or_stability_or_predictive_claim=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh diagnostic output required')
    result = diagnose()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(result, indent=2, allow_nan=False))
