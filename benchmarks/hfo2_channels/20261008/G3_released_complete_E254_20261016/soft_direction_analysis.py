"""Offline diagnosis of the measured preserving R34 sign refusal; no DFT."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import numpy as np
from scripts.hfo2_G3_conditional_faces import load_inputs

BASE = ROOT / 'benchmarks/hfo2_channels/20261008'
PINS = {
    'G3_released_curvature_E181_20261015/config.json':
        '71a2e1b8a46b2ef903e4023d329109c5e74975aaf5bab54f4d4638d15c1add39',
    'G3_released_complete_E254_20261016/native_audit.json':
        'a389154749a07109a3a91499343b86f06a2c185fd310611570637d71193e3381',
}


def analyze():
    sources = {}
    for name, digest in PINS.items():
        raw = (BASE / name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('original audited input bytes differ')
        sources[name] = json.loads(raw)
    config, audit = (sources[n] for n in PINS)
    if (audit['mode'] != 'complete' or audit['native_new_stencil_SCFs_reparsed'] != 272
            or not audit['accounting']['normal_completion']):
        raise ValueError('complete original native evidence required')
    spec = config['channels']['preserving']
    row = next(r for r in audit['channels'] if r['channel'] == 'preserving')
    directions = np.asarray(spec['literal_directions_40_by_34'])
    released = np.asarray(spec['released_basis_36_by_34'])
    anchor = np.asarray(spec['literal_native_anchor_40_A'])
    gradient = np.asarray(spec['native_anchor_result']['gradient'])
    data = load_inputs('preserving')  # audited public chart/reference, no evaluator
    species = data['chart'].reference_atoms.get_chemical_symbols()
    if species != ['Hf']*4 + ['O']*8 or directions.shape != (40, 34):
        raise ValueError('original ordered Hf4O8 length chart required')
    if not np.allclose(data['basis'] @ released, directions, atol=2e-12, rtol=0):
        raise ValueError('literal native embedding differs')
    raw_blocks = [np.asarray(b['raw_projected_eV_A2']) for b in row['blocks']]
    symmetric = [(h+h.T)/2 for h in raw_blocks]
    pairs = [np.linalg.eigh(h) for h in symmetric]
    if ([b['step_A'] for b in row['blocks']] != [.01, .02]
            or any(np.count_nonzero(w < 0) != 1 for w, _ in pairs)):
        raise ValueError('diagnosis requires the actual two-step one-negative cohort')
    z = pairs[0][1][:, 0].copy()  # first-step direction, never select a step by outcome
    full = directions @ z
    pivot = int(np.flatnonzero(np.abs(full) >= np.max(np.abs(full))-1e-12)[0])
    if full[pivot] < 0:
        z *= -1
        full *= -1  # deterministic sign, not a downhill/label-dependent choice
    if (not np.allclose(full @ full, 1., atol=2e-12, rtol=0)
            or full[-1] != 0.):
        raise ValueError('unit fixed-control internal direction required')
    old_axes = data['axes']
    internal = released @ z
    promoted = np.column_stack((old_axes, internal))
    if not np.allclose(promoted.T @ promoted, np.eye(3), atol=2e-12, rtol=0):
        raise ValueError('third direction must be orthogonal to unchanged Q1/Q2')
    # Restriction at the joint candidate, NOT another atomic-only relaxation.
    _, cell_singular, vt = np.linalg.svd(directions[36:39], full_matrices=True)
    rank = int(np.count_nonzero(cell_singular > 1e-12))
    if rank != 3:
        raise ValueError('original three open-cell constraints required')
    atomic = vt[rank:].T
    remaining = pairs[0][1][:, 1:]  # same first-step soft direction removed at both steps
    atomic_weights = full[:36].reshape(12, 3)**2
    scalar = float(anchor @ full)
    radius = config['policy']['maximum_internal_radius_A']
    discriminant = scalar*scalar + radius*radius - float(anchor @ anchor)
    if discriminant < 0:
        raise ValueError('original anchor is outside its original radius')
    return dict(
        schema_version=1,
        status='existing_data_soft_direction_diagnosis_not_physical_instability_certificate',
        sources=PINS, channel='preserving', step_A=[.01, .02],
        coordinate_metric='original reference-atomic/open-cell length chart; not mass weighted',
        sign_rule='first full40 component within1e-12 of maximum absolute value positive; first-step eigenvector',
        soft_direction_full40=full.tolist(), soft_direction_released34=z.tolist(),
        original_Q_overlap=(old_axes.T @ internal).tolist(),
        minimum_eigenvalues_eV_A2=[float(w[0]) for w, _ in pairs],
        next_eigenvalues_eV_A2=[float(w[1]) for w, _ in pairs],
        two_step_soft_squared_overlap=float((z @ pairs[1][1][:, 0])**2),
        first_direction_Rayleigh_eV_A2=[float(z @ h @ z) for h in symmetric],
        coordinate_weights=dict(Hf=float(atomic_weights[:4].sum()),
            O=float(atomic_weights[4:].sum()), open_cell=float(full[36:39] @ full[36:39]),
            atomic_reference_axes=atomic_weights.sum(axis=0).tolist()),
        coordinate_weights_are_energy_fractions_or_phonon_labels=False,
        soft_gradient_eV_A=float(gradient @ full),
        full_released_gradient_norm_eV_A=spec['base_released_gradient_norm_eV_A'],
        atomic_restriction_dimension=atomic.shape[1],
        cell_constraint_singular_values=cell_singular.tolist(),
        atomic_restriction_minima_eV_A2=[float(np.linalg.eigvalsh(atomic.T @ h @ atomic)[0]) for h in symmetric],
        after_formal_third_coordinate_remaining_dimension=remaining.shape[1],
        remaining_raw_minima_eV_A2=[float(np.linalg.eigvalsh(remaining.T @ h @ remaining)[0]) for h in symmetric],
        anchor_radius_A=float(np.linalg.norm(anchor)),
        radius_only_signed_interval_A=[-scalar-float(np.sqrt(discriminant)),
                                      -scalar+float(np.sqrt(discriminant))],
        radius_interval_is_complete_geometry_or_DFT_admission=False,
        full_new_retained_Hessian_stationarity_or_TS_certified=False,
        independent_accuracy_or_stable_nonlinear_branch_certified=False,
        original_independent_sign_admission_passed=False,
        new_DFT_inputs_or_scheduler_writes=0,
        new_training_reservations=0, original_independent_slots_computed=0,
        physical_plus005_generated_or_read=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh analysis output required')
    result = analyze()
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(result['status'])
