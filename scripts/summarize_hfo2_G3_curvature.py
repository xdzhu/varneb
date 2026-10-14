"""Summarize measured fixed-strain complements, not certify a TS or select axes.

The native auditor supplies raw provenance. This offline arithmetic replay
does not reread native SCFs, fill missing columns, modify forces or run DFT.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

MANIFEST_SHA256 = '5f047f54ca83cb2dceb8cb71687eb9206a85c04205cdf26431565bc8a9b98a5a'


def block_diagnostics(raw_projected, basis, gradient, tangent):
    """Raw signs in the declared metric; all tolerances here are algebraic only."""
    a, b, g, t = map(lambda x: np.asarray(x, dtype=float),
                     (raw_projected, basis, gradient, tangent))
    if (a.shape != (37, 37) or b.shape != (40, 37)
            or g.shape != (40,) or t.shape != (40,)
            or not all(np.isfinite(x).all() for x in (a, b, g, t))):
        raise ValueError('complete finite37-column measured block required')
    if (not np.allclose(b.T @ b, np.eye(37), atol=2e-12, rtol=0)
            or not np.array_equal(b[-1], np.eye(37)[-1])
            or np.any(b[:-1, -1]) or t[-1] != 0):
        raise ValueError('orthonormal36-internal plus separate control required')
    internal = (a[:36, :36] + a[:36, :36].T) / 2
    values, vectors = np.linalg.eigh(internal)
    gi = b[:, :36].T @ g
    ti = b[:, :36].T @ t
    ti_norm = np.linalg.norm(ti)
    if ti_norm == 0:
        raise ValueError('native tangent has no measured internal component')
    excluded_tangent = t-b[:, :36] @ ti
    ti = ti / ti_norm
    # SVD chooses only an orthogonal complement; its individual axes have no
    # mode identity and are never exported as phonon eigenvectors.
    complement = np.linalg.svd(ti[None, :], full_matrices=True)[2][1:].T
    transverse_values = np.linalg.eigvalsh(complement.T @ internal @ complement)
    negative = vectors[:, values < 0]
    positive = vectors[:, values > 0]
    return dict(
        fixed_strain_internal_dimension=36,
        external_control_excluded=True,
        symmetrized_internal_eigenvalues_eV_A2=values.tolist(),
        raw_negative_count=int(np.count_nonzero(values < 0)),
        counts_are_raw_not_resolution_certified=True,
        full_gradient_norm_eV_A=float(np.linalg.norm(g)),
        measured_internal_gradient_norm_eV_A=float(np.linalg.norm(gi)),
        prescribed_control_gradient_eV_A=float(g[-1]),
        unmeasured_translation_gradient_norm_eV_A=float(np.linalg.norm(g-b @ (b.T @ g))),
        native_tangent_norm=float(np.linalg.norm(t)),
        native_tangent_internal_norm=float(ti_norm),
        native_tangent_excluded_translation_norm=float(np.linalg.norm(excluded_tangent)),
        normalized_internal_tangent_lowest_direction_absolute_overlap=float(abs(vectors[:, 0] @ ti)),
        tangent_negative_subspace_squared_overlap=float(np.linalg.norm(negative.T @ ti)**2),
        tangent_complement_eigenvalues_eV_A2=transverse_values.tolist(),
        tangent_complement_raw_negative_count=int(np.count_nonzero(transverse_values < 0)),
        positive_spectral_complement_dimension=positive.shape[1],
        positive_spectral_complement_gradient_norm_eV_A=float(np.linalg.norm(positive.T @ gi)),
        positive_spectral_complement_minimum_eigenvalue_eV_A2=(
            float(values[values > 0][0]) if positive.shape[1] else None),
        minimum_retained_dimension_for_raw_positive_complement=int(np.count_nonzero(values <= 0)),
        retained_axes_or_conditional_plane_selected=False,
        physical_error_bound_or_stationary_index_or_TS_certificate=False)


def summarize(manifest, reports):
    rows = []
    for channel in ('T', 'preserving'):
        r = reports[channel]
        n = r['completed_records']
        expected_status = ('complete_local_curvature_diagnostic' if n == 148
                           else 'partial_local_curvature_diagnostic')
        if (r['channel'] != channel or r['planned_records'] != 148
                or not 0 <= n <= 148 or r['raw_SCFs_reparsed'] != n
                or r['status'] != expected_status
                or len(r['native_record_sources']) != n or r['failures_present']
                or r['physical_error_bound_or_TS_certificate']
                or r['physical_reserved_geometry_generated_or_read']
                or r['new_DFT_or_scheduler_operations'] != 0):
            raise ValueError('failed, remapped or non-native audited prefix')
        center = manifest['centers'][channel]
        if r['center_source'] != center['raw_source']:
            raise ValueError('native report belongs to a different physical center')
        expected_blocks = [(h, max(0, min(74, n-74*i))//2)
                           for i, h in enumerate((.01, .02))]
        if [(b['step_A'], b['complete_pairs']) for b in r['blocks']] != [
                (h, pairs) for h, pairs in expected_blocks if pairs]:
            raise ValueError('measured blocks differ from the ordered native prefix')
        blocks = []
        for block in r['blocks']:
            h = block['step_A']
            if h not in (.01, .02):
                raise ValueError('unregistered amplitude')
            if block['full_stencil_at_this_step'] != (block['complete_pairs'] == 37):
                raise ValueError('full-stencil flag differs from measured pair count')
            if not block['full_stencil_at_this_step']:
                if block['fixed_strain_internal_eigenvalues_eV_A2'] is not None:
                    raise ValueError('partial block cannot assign a full internal spectrum')
                continue
            if block['complete_pairs'] != 37:
                raise ValueError('complete37signed pairs required')
            diagnostic = block_diagnostics(block['raw_projected_eV_A2'], center['directions'],
                                           r['raw_center_gradient_eV_A'], center['ordered_local_tangent'])
            if not np.allclose(diagnostic['symmetrized_internal_eigenvalues_eV_A2'],
                               block['fixed_strain_internal_eigenvalues_eV_A2'], atol=2e-12, rtol=0):
                raise ValueError('cached internal spectrum differs from raw arithmetic replay')
            blocks.append(dict(step_A=h, **diagnostic,
                raw_reciprocity_relative_defect=block['reciprocity_relative_defect'],
                symmetrization_operator_change_eV_A2=block['symmetrization_operator_change_eV_A2'],
                maximum_energy_gradient_diagonal_defect_eV_A2=float(np.max(np.abs(
                    block['energy_gradient_diagonal_defect_eV_A2'])))))
        rows.append(dict(channel=channel, audited_prefix_records=n, planned_records=148,
            complete_amplitudes=[b['step_A'] for b in blocks], blocks=blocks,
            mean_audited_native_call_seconds=r['sum_audited_call_seconds']/n if n else None,
            two_step_comparison_available=r['two_step_comparison'] is not None))
    return dict(status='offline_raw_complement_diagnostics_not_stability_or_TS_certificate', rows=rows,
        new_DFT_or_scheduler_operations=0, new_raw_SCF_reparse=False,
        conditional_axes_or_forecasts_selected=False, physical_reserved_geometry_generated_or_read=False,
        limitations=[
            'Signs refer to a finite-step symmetrized Hessian in the declared nonstationary chart.',
            'One retained internal direction cannot leave a positive complement when there are two negative directions.',
            'The displayed positive spectral complement is an algebraic diagnostic, not a chosen conditional branch.',
            'The native tangent is projected only for this internal-space diagnostic; excluded translation norm is reported, not erased from the source.',
            'Step spread, reciprocity and energy/force disagreement do not bound total DFT or barrier error.',
            'No source force/ASR modification, phonon-frequency identification, release or surface minimization.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--T-report', type=Path, required=True)
    parser.add_argument('--T-sha256', required=True)
    parser.add_argument('--preserving-report', type=Path, required=True)
    parser.add_argument('--preserving-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh arithmetic report required')
    paths = dict(manifest=(args.manifest, MANIFEST_SHA256),
        T=(args.T_report, args.T_sha256), preserving=(args.preserving_report, args.preserving_sha256))
    loaded = {}
    for key, (path, expected) in paths.items():
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != expected:
            raise ValueError(f'{key} original input bytes differ')
        loaded[key] = json.loads(data)
    result = summarize(loaded['manifest'], loaded)
    result['input_sha256'] = {key: expected for key, (_, expected) in paths.items()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=result['status'], full_amplitudes=[r['complete_amplitudes'] for r in result['rows']])))


if __name__ == '__main__':
    main()
