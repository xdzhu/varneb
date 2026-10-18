"""Inert two-sign controls for two measured preserving neighborhoods."""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import sys

REPOSITORY = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPOSITORY))
import numpy as np
from ase.io import read
from scripts.audit_hfo2_static_replica import sha256
from scripts.hfo2_G3_conditional_faces import load_inputs, PINS
from scripts.plot_hfo2_mechanism import rotated_t_triplet
from vcneb import ModePlane, validate_path_geometry
from vcneb.mode_surface import _orthogonal_directions

BASE = REPOSITORY/'benchmarks/hfo2_channels/20261008'
POLICY = dict(channel='preserving', prepared_strains=[0., .01], midpoint_admitted=False,
    retained=['Qx', 'Qz', 'Qy', 'polar_proxy_z'], source_probe_step_A=.01,
    signed_atomic_starts_A=[-.1, .1], maximum_attempts_per_arm=12,
    maximum_prospective_optimization_attempts=96, signed_direction_probe_step_A=.01,
    maximum_directional_check_attempts=8, maximum_prospective_attempts=104,
    gradient_tolerance_eV_A=.03,
    initial_trust_radius_A=.05, maximum_coordinate_norm_A=.5,
    numerical_seed='unit branch coordinate plus on-site Jacobi complement',
    numerical_branch_curvature_eV_A2=1.,
    eigenvalue_clipping=False, new_full_Hessian=False, execution_admitted=False)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec); spec.loader.exec_module(result)
    return result


def save(path, data):
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(data, handle, indent=2, allow_nan=False); handle.write('\n')


def release(data, n):
    q = np.asarray(data['orthonormal_39_by_4']); t = np.zeros((39, 3))
    t[:36] = np.tile(np.eye(3), (12, 1))
    plane = ModePlane(np.zeros(n), q[:n, :2], np.ones(n), np.eye(2), ('1', '2'), 'A', 'same-four-patterns')
    u = _orthogonal_directions(plane, np.column_stack([q[:n, 2:], t[:n]]))
    if u.shape[1] != (29 if n == 36 else 32):
        raise ValueError('physical release space changed')
    return u


def numerical_branch_initializer(onsite, direction):
    """Default BFGS unit scale along the branch, on-site metric elsewhere.

    This replaces no eigenvalue of the physical Hessian. It separates one
    explicitly chosen numerical coordinate from its complement.
    """
    direction = np.asarray(direction)
    if not np.isclose(direction@direction, 1., atol=1e-12, rtol=0):
        raise ValueError('unit physical branch direction required')
    np.linalg.cholesky(onsite)
    vv = np.outer(direction, direction); p = np.eye(len(direction))-vv
    result = p@onsite@p+POLICY['numerical_branch_curvature_eV_A2']*vv
    np.linalg.cholesky(.5*(result+result.T))
    return result


def prepare(output):
    study = module(BASE/'matched_pattern_release_E270_20261018/study.py', 'E272_frozen_reader')
    numeric = module(BASE/'matched_pattern_release_E270_20261018/numerical_seed.py', 'E272_transport')
    register_path = BASE/'midpoint_reverse_E269_20261018/six_neighborhoods.json'
    selected = [r for r in json.loads(register_path.read_bytes())['rows']
                if r['channel'] == 'preserving' and r['substrate_strain'] in POLICY['prepared_strains']]
    if [r['substrate_strain'] for r in selected] != POLICY['prepared_strains']:
        raise ValueError('only the original zero and plus1 preserving peaks required')
    output = Path(output); output.mkdir(parents=True, exist_ok=False)
    t_path = BASE/'reference_variants/T.vasp'; t = read(t_path)
    orientation_path = BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json'
    rotation = np.asarray(json.loads(orientation_path.read_bytes())['orientation']['cartesian_rotation_rows'])
    if not np.array_equal(rotation, np.eye(3)[[2, 0, 1]]):
        raise ValueError('original fixed crystal frame required')
    _, patterns, _ = rotated_t_triplet(t)
    rows = []
    for index, row in zip((0, 2), selected):
        obs_path = REPOSITORY/row['observation']
        if sha256(obs_path) != row['observation_sha256']:
            raise ValueError('registered original peak changed')
        obs = json.loads(obs_path.read_bytes()); peak = row['sampled_peak_image']
        if (peak != int(np.argmax(obs['relative_enthalpy_meV_fu'])) or obs['climb']
                or obs['pressure_GPa'] != 0 or not obs['ordinary_residual_pass']):
            raise ValueError('original ordinary zero-pressure sampled peak required')
        raw = obs['raw_image_evaluations'][peak]; seed = output/'seeds'/str(index)
        (seed/'raw/OUT.ABACUS').mkdir(parents=True)
        source = obs_path.parent/'raw'/f'image_{peak:04d}'
        hashes = {n: raw['input_sha256'][n] for n in ('INPUT', 'KPT', 'STRU')}
        hashes['OUT.ABACUS/running_scf.log'] = raw['raw_log_sha256']
        for name, h in hashes.items():
            if sha256(source/name) != h:
                raise ValueError('original native frozen source changed')
            shutil.copyfile(source/name, seed/'raw'/name)
        atoms = study.read_fixed_hfo2_stru(seed/'raw/STRU')
        transform = np.linalg.inv(atoms.cell.array)@rotation@t.cell.array
        c = np.zeros((39, 4))
        for j, m in enumerate((0, 2, 1)): c[:36, j] = (patterns[m]@transform.T).ravel()
        weights = np.r_[np.full(4, -.25), np.full(8, .125)]
        c[:36, 3] = (weights[:, None]*transform[:, 2]).ravel()
        q = np.vstack([np.linalg.qr(c[:36], mode='reduced')[0], np.zeros((3, 4))])
        data = dict(index=index, channel='preserving', substrate_strain=row['substrate_strain'],
            source_job=row['source_job'], observation=row['observation'], observation_sha256=row['observation_sha256'],
            native_files_sha256=hashes, physical_input_sha256=raw['input_sha256'],
            results=dict(energy=raw['energy_eV_cell'], forces=raw['forces_eV_A'], stress=raw['stress_ASE_voigt_eV_A3']),
            covectors_39_by_4=c.tolist(), orthonormal_39_by_4=q.tolist())
        save(seed/'seed.json', data)
        chart = study.chart_for(seed); energy, gradient = study.cached_frozen(seed, chart)
        rows.append(dict(index=index, seed=seed, data=data, chart=chart, energy=energy, gradient=gradient))
    measured = load_inputs('preserving'); origin = rows[0]; chart0 = origin['chart']
    offset = chart0.reference_atoms.get_scaled_positions(wrap=False)-measured['chart'].reference_atoms.get_scaled_positions(wrap=False)
    if (origin['data']['native_files_sha256']['STRU'] != measured['manifest']['centers']['preserving']['raw_source']['input_sha256']['STRU']
            or not np.allclose(chart0.reference_cell, measured['chart'].reference_cell, atol=1e-12, rtol=0)
            or not np.allclose(offset, np.rint(offset), atol=1e-12, rtol=0)):
        raise ValueError('raw measured source and fixed periodic gauge must match native zero-strain P')
    b = measured['basis'][:39, :]; u = release(origin['data'], 36)
    diagnostics = []; directions = []
    for block in measured['curvature']['blocks']:
        raw = np.asarray(block['raw_projected_eV_A2'])[:36, :36]
        h = b@(.5*(raw+raw.T))@b.T
        eig, vec = np.linalg.eigh(u.T@h[:36, :36]@u)
        v = np.r_[u@vec[:, 0], np.zeros(3)]
        anchor = 12+int(np.argmax(np.abs(v[12:36])))
        if v[anchor] < 0: v = -v
        directions.append(v)
        joint = release(origin['data'], 39)
        je = np.linalg.eigvalsh(joint.T@h@joint)
        diagnostics.append(dict(probe_step_A=block['step_A'], atomic_min_eigenvalue_eV_A2=float(eig[0]),
            joint_min_eigenvalue_eV_A2=float(je[0]), atomic_negative_count=int(np.sum(eig < 0)),
            gradient_along_signed_direction_eV_A=float(origin['gradient']@v),
            released_atomic_gradient_norm_eV_A=float(np.linalg.norm(u.T@origin['gradient'][:36])),
            oxygen_squared_norm_fraction=float(v[12:36]@v[12:36]),
            physical_stationarity_or_stability_certified=False))
        if block['step_A'] == .01: source_hessian = h
    if diagnostics[0]['probe_step_A'] != .01 or not all(d['atomic_negative_count'] == 1 for d in diagnostics):
        raise ValueError('preselected measured source and one raw negative atomic direction required')
    onsite = np.zeros((39, 39))
    for i in range(12): onsite[3*i:3*i+3, 3*i:3*i+3] = source_hessian[3*i:3*i+3, 3*i:3*i+3]
    onsite[36:, 36:] = source_hessian[36:, 36:]
    np.linalg.cholesky(onsite)
    numerical = numerical_branch_initializer(onsite, directions[0])
    reports = []
    for row in rows:
        transform = numeric.chart_transport(chart0, row['chart'])
        c0, c = (np.asarray(r['data']['covectors_39_by_4']) for r in (origin, row))
        if not np.allclose(transform.T@c0, c, atol=1e-12, rtol=0):
            raise ValueError('physical pattern meanings changed in chart transport')
        direction = np.linalg.solve(transform, directions[0]); direction /= np.linalg.norm(direction)
        matrix = transform.T@numerical@transform
        for n in (36, 39):
            ur = release(row['data'], n); np.linalg.cholesky(ur.T@matrix[:n, :n]@ur)
        starts = [amplitude*direction for amplitude in POLICY['signed_atomic_starts_A']]
        for start in starts:
            if np.any(start[36:]) or np.max(np.abs(c.T@start)) > 1e-12:
                raise ValueError('signed atomic seed must keep exact cell and four held physical quantities')
            validate_path_geometry([row['chart'].reference_atoms, row['chart'].to_atoms(start)],
                cell_scale=study.SCALE, minimum_distance=1.6, maximum_deformation=.25)
        report = dict(index=row['index'], strain=row['data']['substrate_strain'],
            frozen_energy_eV_cell=row['energy'], signed_starts=[s.tolist() for s in starts],
            numerical_hessian_39=matrix.tolist(), execution_admitted=False, new_DFT_calls=0)
        save(row['seed']/'branch_preparation.json', report)
        reports.append({k: v for k, v in report.items() if k not in ('signed_starts', 'numerical_hessian_39')})
    config = dict(policy=POLICY, rows=reports, source_curvature_native_sha256=PINS['preserving'][0],
        source_registration_sha256=sha256(register_path), T_reference_sha256=sha256(t_path),
        orientation_sha256=sha256(orientation_path), diagnostics=diagnostics,
        signed_direction_probe_overlap=float(directions[0]@directions[1]),
        new_DFT_calls=0, new_submissions=0, execution_admitted=False)
    save(output/'config.json', config)
    return config


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.output)))
