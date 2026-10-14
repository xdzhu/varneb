"""All registered references at the same two training observables, without DFT.

Private source arrays stay in memory. Scalar diagnostics distinguish the
complete-space harmonic control from restricted spans with omitted modes
clamped. Neither representation nor raw positivity certifies a DFT branch.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from scripts import hfo2_G3_conditional_faces as faces
from scripts.verify_hfo2_reference_transport import load_registered_reference_bank, MASS, REGISTRATION_SHA256
from vcneb import transport_atomic_reference_span
from vcneb.affine_quadratic import constrain_quadratic
from vcneb.quadratic_reduction import condition_quadratic_energy

BASE = ROOT / 'benchmarks/hfo2_channels/20261008'
FACE_PINS = {
    'T': ('G3_feasible_execution_E171_20261015/T_complete_native_audit.json',
          'b8e7836722e3c7cadfaed9771eeb90e1d53097a03fc5d19272109eca0a47c166'),
    'preserving': ('G3_preserving_terminal_E179_20261015/preserving_complete_native_audit.json',
          '982d2657a7a2c536509d5a8c460c1ccc0081028b367de398847720fe46937d5a'),
}
NAMES = [('T_geometric_triplet', 3), ('T_two_complete_optical_doublets', 4),
         ('T_all_optical', 33), *((f'Cmma_frame_{n}', 4) for n in (103, 106, 114, 127))]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def observed_points(channel, data):
    name, digest = FACE_PINS[channel]
    path = BASE / name
    if sha(path) != digest:
        raise ValueError('original complete native face bytes differ')
    proof = json.loads(path.read_bytes())
    if (proof['audit_mode'] != 'complete' or proof['channel'] != channel
            or proof['completed_nodes_audited'] != 25 or proof['completed_face_sha256'] is None
            or proof['new_DFT_or_scheduler_writes'] != 0
            or proof['physical_reserved_geometry_generated_or_read'] is not False):
        raise ValueError('complete original same-chart training face required')
    observed = {0: data['origin_point']['energy_eV_cell'] - data['center_energy']}
    specs = faces.ordered_points()
    for row in proof['candidate_rows']:
        if row['id'] != specs[row['id']]['id'] or row['q_A'] != specs[row['id']]['q_A']:
            raise ValueError('original fixed Q identity differs')
        if row['candidate'] is not None:
            if row['candidate']['released_gradient_norm_eV_A'] > .03:
                raise ValueError('native gradient candidate required')
            observed[row['id']] = row['candidate']['energy_eV_cell'] - data['center_energy']
    if len(observed) != proof['gradient_candidates']:
        raise ValueError('observed candidate coverage differs')
    return observed


def statistics(points, ids):
    selected = [p for p in points if p['id'] in ids and p['id'] != 0]
    if any(p['DFT_energy_change_eV_cell'] is None for p in selected):
        raise ValueError('statistics require all declared DFT candidates')
    return dict(ids=[p['id'] for p in selected], count=len(selected), denominator_non_origin=24,
        max_absolute_energy_discrepancy_meV_cell=None if not selected else max(
            abs(p['model_minus_DFT_energy_eV_cell']) for p in selected) * 1000,
        max_absolute_seen_origin_aligned_shape_discrepancy_meV_cell=None if not selected else max(
            abs(p['seen_origin_aligned_shape_discrepancy_eV_cell']) for p in selected) * 1000)


def analyze(source_root):
    t, bank, registered, private_pins = load_registered_reference_bank(source_root)
    if [(n, q.shape[1]) for n, q in bank] != NAMES:
        raise ValueError('all seven original complete reference spans required')
    rotation = np.eye(3)[[2, 0, 1]]
    specs = faces.ordered_points()
    rows, equalities = [], []
    for channel in ('T', 'preserving'):
        data = faces.load_inputs(channel)
        b, a = data['basis'], data['axes']
        g = b.T @ data['center_gradient']
        observed = observed_points(channel, data)
        cell = b[36:39].T
        allowed = []
        for label, source in bank:
            moved = transport_atomic_reference_span(t, source, MASS, data['chart'],
                source_to_target_cell_basis=rotation)
            u = np.column_stack((b.T @ moved.local_basis, cell))
            if label == 'T_all_optical' and not np.allclose(u @ u.T, np.eye(36), atol=2e-12, rtol=0):
                raise ValueError('complete T optical plus open-cell span must be complete')
            allowed.append((label, moved.local_basis.shape[1], u))
        allowed.append(('adaptive_two_directions_plus_complete_complement', None, np.eye(36)))
        for block in data['curvature']['blocks']:
            raw = np.asarray(block['raw_projected_eV_A2'])[:36, :36]
            h = (raw + raw.T) / 2
            r = np.linalg.svd(a.T, full_matrices=True)[2][2:].T
            schur = condition_quadratic_energy(h, g, data['center_energy'], a, r, stability_floor=0.)
            current, models = [], {}
            for label, rank, u in allowed:
                # Outside-reference directions stay clamped here, unlike the
                # E159 test that releases the complete outside-reference space.
                model = constrain_quadratic(h, g, u, a, stability_floor=0.)
                models[label] = model
                origin = model.energy_change([0., 0.])
                points = []
                for spec in specs:
                    q = np.array(spec['q_A'])
                    x = model.coordinates(q)
                    energy = model.energy_change(q)
                    dft = observed.get(spec['id'])
                    inside = bool(np.linalg.norm(x) <= .5)
                    geometry = bool(faces.canary.geometry_guard(data['chart'], b, x))
                    points.append(dict(id=spec['id'], q_A=spec['q_A'],
                        model_energy_change_eV_cell=energy, DFT_energy_change_eV_cell=dft,
                        model_minus_DFT_energy_eV_cell=None if dft is None else energy - dft,
                        seen_origin_aligned_shape_discrepancy_eV_cell=None if dft is None else
                            (energy - origin) - (dft - observed[0]),
                        model_total_internal_norm_A=float(np.linalg.norm(x)),
                        inside_original_radius=inside, passes_original_geometry_guard=geometry,
                        algebraic_residuals=model.residuals(q),
                        outside_guard_values_diagnostic_only_not_material_predictions=not geometry))
                current.append(dict(channel=channel, step_A=block['step_A'], reference=label,
                    nominal_atomic_reference_rank=rank, shared_observable_count=2,
                    allowed_dimension=u.shape[1], free_dimension=u.shape[1]-2,
                    unmodeled_dimensions_clamped=36-u.shape[1],
                    minimum_raw_constraint_null_curvature_eV_A2=model.minimum_free_curvature,
                    constraint_singular_values=model.constraint_singular_values.tolist(),
                    observed_DFT_candidate_ids=sorted(observed),
                    model_geometry_eligible_ids=[p['id'] for p in points if p['passes_original_geometry_guard']],
                    formal_origin_release_meV_cell=1000*origin,
                    observed_origin_release_meV_cell=1000*observed[0], points=points,
                    all_seen_diagnostic_statistics=statistics(points, set(observed)),
                    own_geometry_eligible_statistics=statistics(points, {p['id'] for p in points
                        if p['id'] in observed and p['passes_original_geometry_guard']})))
            common = {p['id'] for p in specs if p['id'] in observed and all(
                row['points'][p['id']]['passes_original_geometry_guard'] for row in current)}
            for row in current:
                row['common_all_reference_geometry_eligible_statistics'] = statistics(row['points'], common)
            full = models['T_all_optical']
            errors, coordinates = [], []
            for spec in specs:
                q = spec['q_A']
                errors.append(abs(full.energy_change(q) - schur.energy_change(q)))
                coordinates.append(float(np.linalg.norm(full.coordinates(q) - schur.full_displacement(q))))
            if max(errors) > 1e-10 or max(coordinates) > 1e-10:
                raise ValueError('complete fixed-reference and adaptive full-space quadratic disagree')
            equalities.append(dict(channel=channel, step_A=block['step_A'], numeric_Q_points=25,
                max_complete_T_vs_adaptive_Schur_energy_difference_eV_cell=max(errors),
                max_complete_T_vs_adaptive_Schur_coordinate_difference_A=max(coordinates),
                common_all_reference_non_origin_ids=sorted(common-{0}),
                no_non_origin_common_comparison=len(common-{0}) == 0,
                equality_is_standard_same_problem_algebra_not_an_accuracy_improvement=True))
            rows.extend(current)
    if any(sha(source_root / name) != digest for name, digest in private_pins.items()):
        raise ValueError('private reference sources changed during analysis')
    return dict(status='post_hoc_shared_Q_harmonic_reference_controls_not_independent_material_prediction',
        rows=rows, complete_space_equivalence=equalities, model_point_rows=sum(len(r['points']) for r in rows),
        energy_unit='eV per twelve-atom Hf4O8 cell; shape summaries in meV per cell', formula_units=4,
        same_source_Hessian_gradient_Q_and_mechanical_boundary=True,
        original_two_probe_steps_and_all_four_Cmma_frames_retained=True,
        restricted_reference_controls_clamp_unmodeled_directions=True,
        restricted_controls_are_not_full_outside_reference_release_B4=True,
        raw_zero_stability_floor_not_physical_resolution=True, fitted_coefficients=0,
        original_complete_face_sha256={c: pin for c, (_, pin) in FACE_PINS.items()},
        source_author_sha256=private_pins, source_registration_sha256=REGISTRATION_SHA256,
        production_T_sha256=registered['production_T_SHA256'],
        production_Gamma_sha256=registered['production_T_Gamma_SHA256'],
        original_curvature_SCFs_additional=296, original_conditional_face_SCFs_additional={'T':118, 'preserving':172},
        Gamma_canary_G2_reference_and_other_preparation_costs_additional=True,
        no_equal_preparation_cost_or_speedup_claim=True,
        new_DFT_or_scheduler_writes=0, new_raw_SCF_reparse=False,
        raw_transformed_author_basis_or_response_arrays_exported=False,
        private_source_root_not_written=True, physical_reserved_geometry_generated_or_read=False,
        independent_material_slots_used=0,
        stable_conditional_branch_TS_physical_error_or_independent_advantage_certified=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('fresh reference comparison output required')
    if a.output.resolve().is_relative_to(a.source_root.resolve()):
        raise ValueError('private sources remain read-only')
    result = analyze(a.source_root)
    with a.output.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps(dict(model_point_rows=result['model_point_rows'], equivalence=result['complete_space_equivalence'])))


if __name__ == '__main__':
    main()
