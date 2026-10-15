"""Freeze all seven registered quadratic reference controls on E185 targets.

Private mode arrays stay in memory. A restricted reference clamps omitted
directions; it is not a full native-chart conditional model or a TS proof.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
BASE = ROOT / 'benchmarks/hfo2_channels/20261008'
CONTROL = BASE / 'G3_shared_Q_controls_E180_20261015/analyze.py'
CONTROL_SHA = '033ce310190bf8498a4d4ad366877a73f9e5e706fa3baef4bd1eb1c24b282b8c'
FORECAST = BASE / 'G3_local_forecast_E185_20261015/forecast_v2.json'
FORECAST_SHA = 'aa77cdd7460ad809a5fbe06e80c684f7f73ce118835067c380ec219b96abd2d8'


def text_sha(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


if text_sha(CONTROL) != CONTROL_SHA:
    raise ValueError('original shared-reference implementation changed')
spec = importlib.util.spec_from_file_location('E180_reference_controls', CONTROL)
control = importlib.util.module_from_spec(spec)
spec.loader.exec_module(control)


def locked_forecast():
    if control.sha(FORECAST) != FORECAST_SHA:
        raise ValueError('original independent coordinates/elementary forecasts changed')
    return json.loads(FORECAST.read_bytes())


def summarize_coverage(rows):
    result = []
    for channel in ('T','preserving'):
        for step in (.01,.02):
            selected = [r for r in rows if r['channel'] == channel and r['step_A'] == step]
            supported = [{p['id'] for p in r['points'] if p['prediction_calibrated_eV_cell'] is not None}
                         for r in selected]
            common = sorted(set.intersection(*supported))
            result.append(dict(channel=channel, step_A=step, denominator=8,
                supported_per_reference={r['reference']:len(ids) for r,ids in zip(selected,supported)},
                common_all_reference_ids=common, independent_labels=0))
    return result


def freeze(source_root):
    locked = locked_forecast()
    t, bank, registered, private_pins = control.load_registered_reference_bank(source_root)
    if [(name, q.shape[1]) for name,q in bank] != control.NAMES:
        raise ValueError('all seven complete original reference spans required')
    rotation = np.eye(3)[[2,0,1]]
    rows, equalities = [], []
    for training in locked['channels']:
        channel = training['channel']
        data = control.faces.load_inputs(channel)
        b, a = data['basis'], data['axes']
        g = b.T @ data['center_gradient']
        origin_y = training['observed_origin_energy_change_eV_cell']
        if origin_y != data['origin_point']['energy_eV_cell'] - data['center_energy']:
            raise ValueError('same visible Q0 energy calibration required')
        allowed = []
        for name, source in bank:
            moved = control.transport_atomic_reference_span(t, source, control.MASS, data['chart'],
                source_to_target_cell_basis=rotation)
            u = np.column_stack((b.T @ moved.local_basis, b[36:39].T))
            if name == 'T_all_optical' and not np.allclose(u @ u.T, np.eye(36), atol=2e-12, rtol=0):
                raise ValueError('complete fixed-reference optical plus cell span required')
            allowed.append((name, moved.rank, u))
        for block in data['curvature']['blocks']:
            raw = np.asarray(block['raw_projected_eV_A2'])[:36,:36]
            h = (raw + raw.T) / 2
            for name, rank, u in allowed:
                model = control.constrain_quadratic(h, g, u, a, stability_floor=0.)
                origin_model = model.energy_change([0.,0.])
                origin_ok = bool(control.faces.canary.geometry_guard(data['chart'], b, model.coordinates([0.,0.])))
                points = []
                for p in training['targets']:
                    if p['independent_DFT_label'] is not None:
                        raise ValueError('test labels must remain unseen')
                    x, energy = model.coordinates(p['q_A']), model.energy_change(p['q_A'])
                    geometry = bool(control.faces.canary.geometry_guard(data['chart'], b, x))
                    points.append(dict(id=p['id'], q_A=p['q_A'],
                        formal_energy_change_eV_cell=energy, model_internal_norm_A=float(np.linalg.norm(x)),
                        passes_original_geometry_guard=geometry,
                        prediction_raw_eV_cell=energy if geometry else None,
                        prediction_calibrated_eV_cell=energy-origin_model+origin_y if geometry and origin_ok else None,
                        reason=None if geometry and origin_ok else 'target_or_calibration_origin_outside_original_guard',
                        algebraic_residuals=model.residuals(p['q_A']), independent_DFT_label=None))
                rows.append(dict(channel=channel, step_A=block['step_A'], reference=name,
                    nominal_atomic_rank=rank, allowed_dimension=u.shape[1], free_dimension=u.shape[1]-2,
                    omitted_directions_clamped=36-u.shape[1],
                    minimum_raw_constraint_null_curvature_eV_A2=model.minimum_free_curvature,
                    formal_origin_energy_change_eV_cell=origin_model,
                    formal_origin_internal_norm_A=float(np.linalg.norm(model.coordinates([0.,0.]))),
                    origin_geometry_eligible=origin_ok, observed_origin_energy_change_eV_cell=origin_y,
                    origin_bias_eV_cell=origin_model-origin_y, points=points))
                if name == 'T_all_optical':
                    key = f"original_H_{block['step_A']:.2f}"
                    differences = [abs(point['formal_energy_change_eV_cell']-p['models'][key]['value'])
                                   for point,p in zip(points, training['targets'])]
                    if max(differences) > 1e-10:
                        raise ValueError('complete reference differs from already frozen same-space Hessian')
                    equalities.append(dict(channel=channel, step_A=block['step_A'],
                        max_energy_difference_from_E185_complete_H_eV_cell=max(differences),
                        same_problem_identity_not_predictive_gain=True))
    if any(control.sha(source_root/name) != digest for name,digest in private_pins.items()):
        raise ValueError('private source bytes changed while reading')
    code = set((*control.faces.canary.CODE, 'scripts/hfo2_G3_conditional_faces.py',
        'scripts/verify_hfo2_reference_transport.py', 'scripts/analyze_hfo2_cmma_path_reference.py',
        'scripts/analyze_hfo2_parent_patterns.py', 'vcneb/reference_transport.py',
        'vcneb/mode_subspaces.py', 'vcneb/qe_modes.py', 'vcneb/phonons.py', 'vcneb/affine_quadratic.py',
        CONTROL.relative_to(ROOT).as_posix(), Path(__file__).resolve().relative_to(ROOT).as_posix()))
    return dict(status='prospective_shared_Q_reference_controls_not_independent_validation', rows=rows,
        scalar_model_points=sum(len(r['points']) for r in rows), coverage=summarize_coverage(rows),
        complete_space_equalities=equalities, formula_units=4, energy_unit='eV per 12-atom Hf4O8 cell',
        original_E185_forecast_sha256=FORECAST_SHA, inherited_scoring_policy=locked['scoring'],
        reference_order_and_ranks=control.NAMES, registration_sha256=control.REGISTRATION_SHA256,
        private_source_file_sha256=private_pins,
        production_T_sha256=registered['production_T_SHA256'], production_Gamma_sha256=registered['production_T_Gamma_SHA256'],
        code_lf_text_sha256={name:text_sha(ROOT/name) for name in sorted(code)},
        private_raw_or_transformed_arrays_exported=False, private_sources_written=False,
        LDA_source_frequencies_used_as_PBE_curvature=False,
        all_four_Cmma_frames_retained_no_postlabel_selection=True,
        restricted_shared_Q_controls_not_complete_native_Cmma_chart_models=True,
        same_H_g_Q_and_mechanical_boundary=True, stability_floor_zero_is_formal_not_physical_resolution=True,
        original296_curvature_and_118_172_face_SCFs_Gamma_G2_canary_preparation_additional=True,
        no_equal_preparation_cost_or_DFT_saving_or_new_harmonic_theory_claim=True,
        original16_test_targets_starts_scoring_and_192attempt_budget_unchanged=True,
        independent_slots_reserved_computed_each=[8,0], training_commitments_unchanged={'T':34,'preserving':35},
        full_B0_B5_final_physical_prediction_and_domain_freeze_complete=False,
        physical_plus005_geometry_or_label_access=False, new_DFT_or_scheduler_writes=0,
        independent_gain_stability_TS_or_JCTC_readiness_claimed=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh forecast file required; never replace a prediction')
    if args.output.resolve().is_relative_to(args.source_root.resolve()):
        raise ValueError('private sources remain read-only')
    report = freeze(args.source_root)
    with args.output.open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(report,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(model_points=report['scalar_model_points'],coverage=report['coverage'],new_DFT=0)))


if __name__ == '__main__':
    main()
