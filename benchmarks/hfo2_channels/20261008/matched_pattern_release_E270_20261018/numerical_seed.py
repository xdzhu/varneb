"""Inert same-source block-Jacobi numerical seed; no new physical curvature.

The off-diagonal block is removed only from this numerical initializer.
Actual DFT energies/forces/stress and their coupling remain unchanged.
No run is authorized by constructing these matrices.
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
from scripts.hfo2_G3_conditional_faces import load_inputs
from vcneb import ModePlane
from vcneb.mode_surface import _orthogonal_directions


def chart_transport(source, target):
    """Target reference-fractional/open-cell increments -> source increments.

    This is an affine chart map, not an assertion that strained structures
    have the same physical Hessian. Atomic patterns retain their meaning.
    """
    h0, h1 = source.reference_cell, target.reference_cell
    d0 = np.asarray(source.mechanical_contract['deformation_basis'])
    d1 = np.asarray(target.mechanical_contract['deformation_basis'])
    if d0.shape != (3, 3, 3) or d1.shape != d0.shape:
        raise ValueError('same three-open-cell chart family required')
    atomic = np.linalg.solve(h1, h0)
    transform = np.zeros((39, 39))
    for i in range(12):
        transform[3*i:3*i+3, 3*i:3*i+3] = atomic.T
    for j in range(3):
        deformation = np.linalg.solve(h0, h1@d1[j].T).T
        transform[36:, 36+j] = np.einsum('aij,ij->a', d0, deformation)
    if np.linalg.matrix_rank(transform) != 39:
        raise ValueError('numerical chart transport must be nonsingular')
    return transform


def seed_for(case, index, *, probe_index=0):
    """Both old probe steps may be diagnosed; execution must freeze one in advance."""
    case = Path(case)
    manifest_path = case/'executed_bundle_files.json'
    if sha256(manifest_path) != json.loads((case/'package_receipt.json').read_bytes())['bundle_pin']:
        raise ValueError('executed source identity changed')
    executed = json.loads(manifest_path.read_bytes())
    if sha256(case/'study.py') != executed['study.py']:
        raise ValueError('original numerical experiment source changed')
    for i in (0, index):
        if sha256(case/'prepared/seeds'/str(i)/'seed.json') != executed[f'seeds/{i}/seed.json']:
            raise ValueError('original source or target pattern constraints changed')
    spec = importlib.util.spec_from_file_location('E270_seed_data', case/'study.py')
    study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)
    source_path = case/'prepared/seeds/0'
    target_path = case/'prepared/seeds'/str(index)
    source_chart, target_chart = study.chart_for(source_path), study.chart_for(target_path)
    study.cached_frozen(source_path, source_chart)
    study.cached_frozen(target_path, target_chart)
    source = json.loads((source_path/'seed.json').read_bytes())
    target = json.loads((target_path/'seed.json').read_bytes())
    data = load_inputs('T')
    offset = source_chart.reference_atoms.get_scaled_positions(wrap=False)-data['chart'].reference_atoms.get_scaled_positions(wrap=False)
    if (source['native_files_sha256']['STRU'] != data['manifest']['centers']['T']['raw_source']['input_sha256']['STRU']
            or not np.allclose(offset, np.rint(offset), atol=1e-12, rtol=0)
            or not np.allclose(source_chart.reference_cell, data['chart'].reference_cell, atol=1e-12, rtol=0)):
        raise ValueError('numerical source matrix must match the original native input, up to fixed integer gauge')
    block = data['curvature']['blocks'][probe_index]
    raw = np.asarray(block['raw_projected_eV_A2'])[:36, :36]
    b = data['basis'][:39, :]
    h = b@(.5*(raw+raw.T))@b.T
    h[:36, 36:] = 0.; h[36:, :36] = 0.  # Numerical block-Jacobi, not physical decoupling.
    transform = chart_transport(source_chart, target_chart)
    cs, ct = (np.asarray(row['covectors_39_by_4']) for row in (source, target))
    if not np.allclose(transform.T@cs, ct, atol=1e-12, rtol=0):
        raise ValueError('numerical transport must preserve the four physical covectors')
    seed = transform.T@h@transform
    q = np.asarray(target['orthonormal_39_by_4']); translation = np.tile(np.eye(3), (12, 1))
    restrictions = []
    for atomic in (True, False):
        axes = q[:36, :2] if atomic else q[:, :2]
        gauge = translation if atomic else np.vstack([translation, np.zeros((3, 3))])
        frozen = np.column_stack([q[:len(axes), 2:], gauge])
        plane = ModePlane(np.zeros(len(axes)), axes, np.ones(len(axes)), np.eye(2), ('1','2'), 'A', 'numerical-seed-only')
        r = _orthogonal_directions(plane, frozen)
        restriction = r.T@(seed[:36, :36] if atomic else seed)@r
        np.linalg.cholesky(.5*(restriction+restriction.T))  # No clipping or added damping.
        eig = np.linalg.eigvalsh(restriction)
        restrictions.append(dict(treatment='A' if atomic else 'J', rank=len(eig),
            minimum_seed_eigenvalue_eV_A2=float(eig[0]), maximum_seed_eigenvalue_eV_A2=float(eig[-1])))
    return seed, dict(source_strain=0., target_strain=target['substrate_strain'],
        measured_probe_step_A=block['step_A'], restrictions=restrictions,
        numerical_block_jacobi_only=True, source_DFT_mixed_block_still_present=True,
        physical_strain_Hessian_or_stability_prediction=False, eigenvalue_clipping=False,
        new_DFT_calls=0, execution_admitted=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    case = Path(__file__).resolve().parent
    rows = [seed_for(case, i, probe_index=j)[1] for j in (0, 1) for i in range(3)]
    with args.output.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(dict(status='inert_numerical_seed_feasibility_not_a_new_run', rows=rows),
                  handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(rows))
