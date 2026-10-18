"""Reuse the measured zero-strain T matrix and actual first trial, without DFT.

This is a local numerical/mechanism diagnostic, not a stationary stability
certificate, a new Hessian measurement or permission to change the live job.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

REPOSITORY = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPOSITORY))
import numpy as np
from scripts.audit_hfo2_static_replica import sha256
from scripts.hfo2_G3_conditional_faces import load_inputs, PINS
from vcneb import ModePlane
from vcneb.mode_surface import _orthogonal_directions


def block_screen(matrix, atomic_rank):
    """Measured symmetric blocks only; no clipping, fitting or stability claim."""
    matrix = np.asarray(matrix, dtype=float)
    if (matrix.shape != (atomic_rank+3, atomic_rank+3) or atomic_rank < 1
            or not np.isfinite(matrix).all()
            or not np.allclose(matrix, matrix.T, atol=1e-12, rtol=0)):
        raise ValueError('finite symmetric atomic plus three-cell matrix required')
    a, c, d = matrix[:atomic_rank, :atomic_rank], matrix[atomic_rank:, atomic_rank:], matrix[:atomic_rank, atomic_rank:]
    wa, wc, wj = np.linalg.eigvalsh(a), np.linalg.eigvalsh(c), np.linalg.eigvalsh(matrix)
    schur = None if wa[0] <= 0 else c-d.T@np.linalg.solve(a, d)
    return dict(atomic_eigenvalue_range_eV_A2=[float(wa[0]), float(wa[-1])],
        atomic_condition_number=None if wa[0] <= 0 else float(wa[-1]/wa[0]),
        cell_block_eigenvalues_eV_A2=wc.tolist(),
        coupled_eigenvalue_range_eV_A2=[float(wj[0]), float(wj[-1])],
        coupled_negative_count=int(np.count_nonzero(wj < 0)),
        cell_Schur_eigenvalues_eV_A2=None if schur is None else np.linalg.eigvalsh(schur).tolist(),
        separate_blocks_positive=bool(wa[0] > 0 and wc[0] > 0),
        removing_mixed_block_removes_raw_negative=bool(wa[0] > 0 and wc[0] > 0 and wj[0] < 0),
        stationary_physical_stability_or_error_bound_certified=False)


def diagnose(case):
    case = Path(case)
    manifest = case/'executed_bundle_files.json'
    receipt = json.loads((case/'package_receipt.json').read_bytes())
    if sha256(manifest) != receipt['bundle_pin']:
        raise ValueError('executed bundle identity changed')
    files = json.loads(manifest.read_bytes())
    seed_path = case/'prepared/seeds/0/seed.json'
    if sha256(seed_path) != files['seeds/0/seed.json']:
        raise ValueError('the actual fixed-pattern source changed')
    spec = importlib.util.spec_from_file_location('E270_read_only_data', case/'study.py')
    study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)
    seed = json.loads(seed_path.read_bytes()); q = np.asarray(seed['orthonormal_39_by_4'])
    data = load_inputs('T')  # Reuses the complete original two-amplitude native audit.
    chart = study.chart_for(seed_path.parent)
    delta = chart.reference_atoms.get_scaled_positions(wrap=False)-data['chart'].reference_atoms.get_scaled_positions(wrap=False)
    reference_offset = np.rint(delta).astype(int)
    e0, g0 = study.cached_frozen(seed_path.parent, chart)
    if (not np.allclose(chart.reference_cell, data['chart'].reference_cell, atol=1e-12, rtol=0)
            or not np.allclose(delta, reference_offset, atol=1e-12, rtol=0)
            or seed['native_files_sha256']['STRU'] != data['manifest']['centers']['T']['raw_source']['input_sha256']['STRU']
            or e0 != data['center_energy']
            or not np.allclose(g0, data['center_gradient'][:39], atol=2e-12, rtol=0)):
        raise ValueError('the reused matrix must refer to this zero-strain native center')
    plane = ModePlane(np.zeros(36), q[:36, :2], np.ones(36), np.eye(2), ('1','2'), 'A', 'diagnostic-only')
    translations = np.tile(np.eye(3), (12, 1))
    ra = _orthogonal_directions(plane, np.column_stack([q[:36, 2:], translations]))
    basis = np.zeros((39, 32)); basis[:36, :29] = ra; basis[36:, 29:] = np.eye(3)
    original_basis = data['basis'][:39, :]
    rows, measured_matrices = [], []
    for block in data['curvature']['blocks']:
        raw = np.asarray(block['raw_projected_eV_A2'])[:36, :36]
        h = original_basis@(.5*(raw+raw.T))@original_basis.T
        measured_matrices.append(h)
        reduced = basis.T@h@basis
        _, vectors = np.linalg.eigh(reduced)
        weak = basis@vectors[:, 0]
        rows.append(dict(probe_step_A=block['step_A'], **block_screen(reduced, 29),
            weak_direction_Hf_squared_norm=float(weak[:12]@weak[:12]),
            weak_direction_oxygen_squared_norm=float(weak[12:36]@weak[12:36]),
            weak_direction_scaled_cell_squared_norm=float(weak[36:]@weak[36:]),
            weak_direction_full_gradient_eV_A=float(data['center_gradient'][:39]@weak)))
    first_path = case/'first_native_result.json'
    first = json.loads(first_path.read_bytes()); x = np.asarray(first['coordinates'])
    norm2 = float(x@x)
    if first['status'] != 'complete' or norm2 <= 0 or np.any(x[36:]):
        raise ValueError('actual completed first atomic trial required')
    trace = dict(step_norm_A=norm2**.5,
        energy_change_meV_fu=250*(first['enthalpy_eV']-e0),
        kappa_energy_eV_A2=2*(first['enthalpy_eV']-e0-float(g0@x))/norm2,
        kappa_force_eV_A2=float(x@(np.asarray(first['gradient'])-g0))/norm2,
        original_measured_Rayleigh_eV_A2=[float(x@h@x/norm2) for h in measured_matrices],
        not_an_accepted_step_or_converged_candidate=True)
    return dict(status='same_source_raw_block_and_first_trial_diagnostic', rows=rows, first_trial=trace,
        source_job='28692775', native_curvature_member='T_native.json', native_curvature_sha256=PINS['T'][0],
        constant_reference_integer_offset=reference_offset.tolist(), native_STRU_identity_checked=True,
        source_bundle_pin=receipt['bundle_pin'], first_native_result_sha256=sha256(first_path),
        new_DFT_calls=0, new_curvature_probes=0, live_source_or_parameters_changed=False,
        conditional_stability_or_TS_or_prediction_advantage_certified=False,
        numerical_fallback_execution_admitted=False,
        limitations=['nonstationary original sampled center, not a saddle-index test',
            'raw finite-step symmetrized matrices, not physical uncertainty bounds',
            'scaled-cell squared norm is neither substrate strain nor energy fraction',
            'zero-strain numerical block positivity does not prove it at strained or released points'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = diagnose(Path(__file__).resolve().parent)
    with args.output.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(dict(blocks=result['rows'], first_trial=result['first_trial'], new_DFT_calls=0)))
