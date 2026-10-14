"""Compare all registered strong-reference release blocks on one measured H.

Offline, raw finite-step diagnostics only. No author arrays are exported, no
new DFT is called and no stationary response, axes or forecast is selected.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from scripts import hfo2_G3_local_curvature as study
from scripts import summarize_hfo2_G3_curvature as measured
from scripts.verify_hfo2_reference_transport import (
    BASE, ROOT, MASS, MANIFEST_SHA, REGISTRATION_SHA256, load_registered_reference_bank)
from vcneb import transport_atomic_reference_span


def _real(values):
    if np.iscomplexobj(values):
        raise ValueError('declared real directions required; no imaginary-part removal')
    return np.array(values, dtype=float, copy=True)


def _complement(container, retained):
    # Retained columns already have audited orthonormality and full rank.
    if (retained.shape[1] > container.shape[1]
            or not np.allclose(container @ (container.T @ retained), retained, atol=2e-12, rtol=0)):
        raise ValueError('retained space must lie wholly in the release container')
    _, s, vh = np.linalg.svd(retained.T @ container, full_matrices=True)
    if not np.allclose(s, np.ones(retained.shape[1]), atol=2e-12, rtol=0):
        raise ValueError('retained rank cannot be dropped')
    return container @ vh[retained.shape[1]:].T


def release_diagnostics(internal_hessian, internal_basis, local_retained, full_gradient):
    """Atomic and joint complement restrictions of a complete measured block.

    Zero is used only for explicitly raw sign counts, never as a registered
    physical stability floor. A positive spectral block is not certified here.
    """
    h, b, q_local, g = map(_real,
        (internal_hessian, internal_basis, local_retained, full_gradient))
    if (h.shape != (36, 36) or b.shape != (40, 36) or g.shape != (40,)
            or q_local.ndim != 2 or q_local.shape[0] != 40
            or not 1 <= q_local.shape[1] <= 33
            or not all(np.isfinite(x).all() for x in (h, b, q_local, g))):
        raise ValueError('complete finite measured internal block and atomic reference required')
    k = q_local.shape[1]
    if (not np.allclose(h, h.T, atol=2e-12, rtol=0)
            or not np.allclose(b.T @ b, np.eye(36), atol=2e-12, rtol=0)
            or np.any(b[-1]) or np.any(q_local[36:])
            or not np.allclose(q_local.T @ q_local, np.eye(k), atol=2e-12, rtol=0)):
        raise ValueError('symmetric internal H, orthonormal spans and zero atomic-reference cell/control entries required')
    q = b.T @ q_local
    if not np.allclose(b @ q, q_local, atol=2e-12, rtol=0):
        raise ValueError('reference not wholly measured; do not fill missing directions')
    # The chart has33 translation-free atomic and3open-row coordinates.
    # Right singular vectors express that atomic subspace inside measured B.
    _, s, vh = np.linalg.svd(b[:36], full_matrices=False)
    if not np.allclose(s, np.r_[np.ones(33), np.zeros(3)], atol=2e-12, rtol=0):
        raise ValueError('registered atomic/open-cell decomposition differs')
    atomic = vh[:33].T
    gi = b.T @ g
    hi_values, hi_vectors = np.linalg.eigh(h)
    negative = hi_vectors[:, hi_values < 0]
    result = dict(retained_rank=k, complete_measured_internal_dimension=36,
        full_internal_raw_negative_count=int(np.count_nonzero(hi_values < 0)),
        mean_raw_negative_subspace_squared_coverage=(
            float(np.linalg.norm(q.T @ negative, ord='fro')**2/negative.shape[1]) if negative.shape[1] else None),
        retained_projected_raw_negative_count=int(np.count_nonzero(np.linalg.eigvalsh(q.T @ h @ q) < 0)),
        retained_gradient_norm_eV_A=float(np.linalg.norm(q.T @ gi)))
    for name, container in (('atomic_release', atomic), ('joint_release', np.eye(36))):
        r = _complement(container, q)
        hr = r.T @ h @ r
        values = np.linalg.eigvalsh(hr)
        nonpositive = int(np.count_nonzero(values <= 0))
        # Raw conditional matrix algebra is reported only for positive R.
        # This is not a stationary-index test or a physical stability decision.
        schur = None
        if not nonpositive:
            coupled = q.T @ h @ r
            reduced = q.T @ h @ q
            if r.shape[1]:
                reduced -= coupled @ np.linalg.solve(hr, coupled.T)
            schur = np.linalg.eigvalsh((reduced+reduced.T)/2)
        result[name] = dict(released_dimension=r.shape[1],
            raw_minimum_eigenvalue_eV_A2=float(values[0]) if values.size else None,
            raw_negative_count=int(np.count_nonzero(values < 0)),
            raw_nonpositive_count=nonpositive,
            released_gradient_norm_eV_A=float(np.linalg.norm(r.T @ gi)),
            raw_positive_block_Schur_retained_negative_count=(int(np.count_nonzero(schur < 0)) if schur is not None else None),
            empty_release_is_not_a_physical_stability_test=not bool(values.size))
    result.update(signs_are_raw_not_resolution_certified=True,
        released_branch_or_stationary_index_or_TS_certified=False)
    return result


def analyze(source_root, manifest, reports):
    measured.summarize(manifest, reports)  # original-center/prefix/completeness guards
    t, bank, registered, pins = load_registered_reference_bank(source_root)
    expected = [('T_geometric_triplet', 3), ('T_two_complete_optical_doublets', 4), ('T_all_optical', 33),
                *((f'Cmma_frame_{i}', 4) for i in (103, 106, 114, 127))]
    if [(name, columns.shape[1]) for name, columns in bank] != expected:
        raise ValueError('complete registered seven-reference bank required')
    p = np.eye(3)[[2, 0, 1]]
    rows = []
    for channel in ('T', 'preserving'):
        chart, _ = study.validate(manifest, channel)
        center = manifest['centers'][channel]
        if not np.allclose(np.asarray(center['substrate_reference_cell_A'])[:2],
                           (p @ t.cell.array @ p.T)[:2], atol=1e-10, rtol=0):
            raise ValueError('original two prescribed substrate rows differ')
        b = np.asarray(center['directions'])[:, :36]
        for block in reports[channel]['blocks']:
            if not block['full_stencil_at_this_step']:
                continue
            a = np.asarray(block['raw_projected_eV_A2'])[:36, :36]
            h = (a+a.T)/2
            for label, source in bank:
                transported = transport_atomic_reference_span(t, source, MASS, chart,
                    source_to_target_cell_basis=p)
                diagnostic = release_diagnostics(h, b, transported.local_basis,
                    reports[channel]['raw_center_gradient_eV_A'])
                rows.append(dict(channel=channel, step_A=block['step_A'], reference=label, **diagnostic))
    return dict(status='registered_strong_reference_raw_complement_diagnostics_not_prediction', rows=rows,
        source_author_sha256=pins, registration_sha256=REGISTRATION_SHA256,
        curvature_manifest_sha256=MANIFEST_SHA,
        production_T_SHA256=registered['production_T_SHA256'],
        production_Gamma_SHA256=registered['production_T_Gamma_SHA256'],
        new_DFT_or_scheduler_operations=0, new_raw_SCF_reaudit=False,
        raw_or_transformed_author_arrays_exported=False, source_modes_forces_or_ASR_changed=False,
        source_frequencies_used_as_PBE_curvature=False, conditional_axes_or_forecasts_selected=False,
        physical_reserved_geometry_generated_or_read=False,
        physical_sign_floor_error_bound_or_stationary_response_certified=False,
        code_sha256={n:study.sha256(ROOT/n) for n in ('scripts/analyze_hfo2_G3_reference_complements.py',
            'scripts/verify_hfo2_reference_transport.py','vcneb/reference_transport.py')})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    for channel in ('T', 'preserving'):
        parser.add_argument('--'+channel+'-report', type=Path, required=True)
        parser.add_argument('--'+channel+'-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh diagnostic report required')
    if args.output.resolve().is_relative_to(args.source_root.resolve()):
        raise ValueError('private author source remains read-only')
    manifest_path = BASE/'G3_curvature_E150_20261014/manifest.json'
    if study.sha256(manifest_path) != MANIFEST_SHA:
        raise ValueError('original frozen physical manifest differs')
    reports, hashes = {}, {}
    for channel in ('T', 'preserving'):
        path, expected = getattr(args, channel+'_report'), getattr(args, channel+'_sha256')
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != expected:
            raise ValueError('native audit report bytes differ')
        reports[channel] = json.loads(data)
        hashes[channel] = expected
    result = analyze(args.source_root, json.loads(manifest_path.read_bytes()), reports)
    result['native_report_sha256'] = hashes
    args.output.parent.mkdir(parents=True, exist_ok=True)
    study.exclusive_json(args.output, result)
    print(json.dumps(dict(status=result['status'], reference_center_amplitude_rows=len(result['rows']), new_DFT=0)))


if __name__ == '__main__':
    main()
