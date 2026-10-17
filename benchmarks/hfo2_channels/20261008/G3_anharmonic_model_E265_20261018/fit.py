"""One registered gradient-enhanced polynomial experiment; no DFT or labels."""
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
BASE = ROOT/'benchmarks/hfo2_channels/20261008'
PINS = {
    'G3_preserving_branches_E260_20261017/native_complete_E264.json': '30971f61ab7d12d8e9bfa729a517603566c20cefa6bf442b06a4100264bd8811',
    'G3_released_complete_E254_20261016/native_audit.json': 'a389154749a07109a3a91499343b86f06a2c185fd310611570637d71193e3381',
    'G3_released_curvature_E181_20261015/config.json': '71a2e1b8a46b2ef903e4023d329109c5e74975aaf5bab54f4d4638d15c1add39',
}
LENGTH = .1
RCOND = 1e-10


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec); spec.loader.exec_module(value)
    return value


def powers(anharmonic):
    """Exact registered monomials in Q1,Q2,zeta,r1,...,r33."""
    triples = [(2, 0, 0), (1, 1, 0), (0, 2, 0)]
    if anharmonic:
        triples += [(0, 0, 3), (0, 0, 4), (1, 0, 2), (0, 1, 2),
                    (2, 0, 2), (1, 1, 2), (0, 2, 2)]
    rows = [list(p)+[0]*33 for p in triples]
    if anharmonic:
        for j in range(33):
            for q in ((0, 0, 2), (1, 0, 2), (0, 1, 2)):
                p = list(q)+[0]*33; p[j+3] = 1; rows.append(p)
    return np.asarray(rows, dtype=int)


def monomials(x, exponents):
    """Values and analytic coordinate derivatives, including at zero."""
    x = np.asarray(x, dtype=float); exponents = np.asarray(exponents)
    if (x.ndim != 2 or exponents.ndim != 2 or x.shape[1] != exponents.shape[1]
            or not np.isfinite(x).all() or not np.isfinite(exponents).all()
            or np.any(exponents < 0) or np.any(exponents != np.floor(exponents))):
        raise ValueError('finite coordinates and nonnegative integer monomials required')
    values = np.ones((len(x), len(exponents)))
    jacobian = np.zeros((len(x), x.shape[1], len(exponents)))
    for k, p in enumerate(exponents):
        active = np.flatnonzero(p)
        for j in active: values[:, k] *= x[:, j]**p[j]
        for j in active:
            derivative = np.full(len(x), float(p[j]))
            for m in active: derivative *= x[:, m]**(p[m]-(m == j))
            jacobian[:, j, k] = derivative
    return values, jacobian


def solve(values, jacobian, energy, gradient, mask):
    """Column-scaled ordinary least squares; deficient fits cannot forecast."""
    design = np.concatenate((values[mask], jacobian[mask].reshape(-1, values.shape[1])))
    target = np.concatenate((energy[mask], gradient[mask].reshape(-1)))
    scale = np.linalg.norm(design, axis=0); scale[scale == 0] = 1.
    coefficients, _, rank, singular = np.linalg.lstsq(design/scale, target, rcond=RCOND)
    identified = int(rank) == design.shape[1]
    return dict(rank=int(rank), parameters=design.shape[1], rows=design.shape[0],
                identified=identified, coefficients=(coefficients/scale if identified else None),
                scaled_singular_values=singular,
                scaled_condition=(float(singular[0]/singular[-1]) if identified else None))


def metrics(values, jacobian, energy, gradient, coefficients, mask):
    if coefficients is None: return None
    e = (values[mask]@coefficients-energy[mask])*250.  # meV/HfO2
    g = (np.einsum('ndp,p->nd', jacobian[mask], coefficients)-gradient[mask])/LENGTH
    return dict(n=int(np.count_nonzero(mask)), energy_max_error_meV_fu=float(np.max(np.abs(e))),
                energy_RMSE_meV_fu=float(np.sqrt(np.mean(e*e))),
                full_gradient_max_norm_error_eV_A=float(np.max(np.linalg.norm(g, axis=1))),
                gradient_component_RMSE_eV_A=float(np.sqrt(np.mean(g*g))))


def run():
    source = {}
    for name, digest in PINS.items():
        raw = (BASE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest: raise ValueError('pinned training source required')
        source[name] = json.loads(raw)
    cohort, curvature, old = [source[n] for n in PINS]
    worker = module('E265_original_training_inputs', BASE/'G3_preserving_branches_E260_20261017/worker.py')
    config = json.loads((BASE/'G3_preserving_branches_E260_20261017/config.json').read_bytes())
    data, _, _, free, _, _ = worker.validate(config, BASE/'G3_preserving_branches_E260_20261017')
    basis, axes = data['basis'], data['axes']; direction = np.asarray(config['branch_direction_36'])
    stable = np.linalg.svd(np.column_stack([axes, direction]).T, full_matrices=True)[2][3:].T
    transform = np.column_stack([axes, direction, stable])
    if not np.allclose(transform.T@transform, np.eye(36), atol=2e-12, rtol=0):
        raise ValueError('original orthonormal chart required')
    original = data['origin_record']; origin = basis.T@np.asarray(original['coordinates'])
    specification = old['channels']['preserving']
    if (not np.array_equal(specification['literal_native_anchor_40_A'], original['coordinates'])
            or not np.allclose(basis@free, specification['literal_directions_40_by_34'], atol=2e-12, rtol=0)):
        raise ValueError('all curvature columns must belong to the same exact native origin')
    records = cohort['native_records']; rows = {r['id']: r for r in cohort['rows']}
    groups = np.asarray([rows[i]['Q_id'] for i in (0,2,4,6,8,1,3,5,7,9) for _ in range(rows[i]['attempted_SCFs'])])
    if len(records) != len(groups) or len(records) != 81: raise ValueError('all81 charged observations required')
    x = np.asarray([r['coordinates'] for r in records])@basis
    gradient = np.asarray([r['gradient'] for r in records])@basis
    energy = np.asarray([r['enthalpy_eV'] for r in records])-original['enthalpy_eV']
    displacement = x-origin; y = displacement@transform/LENGTH
    g0 = basis.T@np.asarray(original['gradient'])
    blocks = next(r for r in curvature['channels'] if r['channel']=='preserving')['blocks']
    results = []
    for block in blocks:
        action = basis.T@np.asarray(block['full_hessian_action_eV_A2'])
        raw = free.T@action; symmetric = .5*(raw+raw.T)
        corrected = action+free@(symmetric-raw)
        extension = corrected@free.T+free@corrected.T-free@symmetric@free.T
        minimum = float(np.linalg.eigvalsh(stable.T@extension@stable)[0])
        if minimum <= 0: raise ValueError('raw positive R33 source block required; no clipping')
        residual_energy = energy-displacement@g0-.5*np.einsum('ni,ij,nj->n', displacement, extension, displacement)
        residual_gradient = LENGTH*(gradient-g0-displacement@extension)@transform
        for anharmonic in (False, True):
            exponents = powers(anharmonic); values, jacobian = monomials(y, exponents)
            mask = np.ones(81, dtype=bool)
            fit = solve(values, jacobian, residual_energy, residual_gradient, mask)
            folds = []
            for group in range(5):
                train = groups != group; check = ~train
                folded = solve(values, jacobian, residual_energy, residual_gradient, train)
                folds.append(dict(held_out_training_Q_id=group, rank=folded['rank'],
                    parameters=folded['parameters'], identified=folded['identified'],
                    metrics=metrics(values, jacobian, residual_energy, residual_gradient, folded['coefficients'], check)))
            results.append(dict(family='anharmonic_retained_residual' if anharmonic else 'full_space_harmonic_reference',
                step_A=block['step_A'], raw_R33_minimum_eV_A2=minimum,
                raw_reciprocity_relative_defect=block['reciprocity_relative_defect'],
                exponents=exponents.tolist(), origin_energy_eV_cell=original['enthalpy_eV'],
                origin_internal_A=origin.tolist(), origin_gradient_internal_eV_A=g0.tolist(),
                chart_transform_36_by_36=transform.tolist(), symmetric_measured_extension_eV_A2=extension.tolist(),
                fit={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in fit.items()},
                in_sample=metrics(values,jacobian,residual_energy,residual_gradient,fit['coefficients'],mask),
                leave_one_training_Q_family_out=folds))
    return dict(status='one_finite_training_only_anharmonic_identifiability_experiment', sources=PINS,
        protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),
        source_LF_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        normalization_length_A=LENGTH, gradient_weight_length_A=LENGTH, scaled_SVD_relative_cutoff=RCOND,
        native_training_observations=81, total_new_fitted_models=4, models=results,
        origin_and_curvature_same_native_point_verified=True,
        training_validation_is_not_independent=True, independent_labels_seen=False,
        actual_P8_forecasts_or_final_domain_frozen=False, authorizes_independent_DFT=False,
        stable_DFT_sheet_TS_or_predictive_gain_certified=False, new_DFT_or_scheduler_operations=0,
        physical_plus005_or_private_arrays_accessed=False, all_prior_costs_additional=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists(): raise FileExistsError('fresh model report required')
    result = run()
    with args.output.open('x',encoding='utf-8',newline='\n') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps([dict(family=m['family'],step_A=m['step_A'],rank=m['fit']['rank'],
        parameters=m['fit']['parameters'],in_sample=m['in_sample'],folds=m['leave_one_training_Q_family_out'])
        for m in result['models']],indent=2))
