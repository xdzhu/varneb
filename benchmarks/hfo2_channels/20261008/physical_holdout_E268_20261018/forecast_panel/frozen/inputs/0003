"""Finite E268 sampled-peak protocol; never weakens the strict selector.

Registration reads only already-seen training evidence. Endpoint preparation
is an explicit later operation; it does not compute DFT or target interiors.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.prepare_hfo2_clamped_chains import endpoint_evidence
from scripts.prepare_hfo2_clamped_endpoints import geometry_guard
from scripts.prepare_hfo2_reference_variants import write_clean_poscar
from scripts.relax_clamped_ase_endpoint import load_seed
from vcneb import clamped_plane_vcneb_boundary, endpoint_structure_record


BASE = Path('benchmarks/hfo2_channels/20261008')
CASE = BASE / 'physical_holdout_E268_20261018'
NETWORK = BASE / 'clamped_plus1_reversing_terminal_E148_20261014/analysis/network_final.json'
NETWORK_SHA = '6b8a206d2199b6727601955cef32b86513e22805c7be5f309e66a4258602145b'
CHANNELS = {'preserving': 'PO_flip_T_pattern_preserving', 'T': 'PO_to_T'}
PHASES = ('PO_plus', 'PO_minus_T_preserving', 'T')
REPRESENTATIONS = ('path_adapted', 'T_triplet', 'T_two_optical_doublets',
                   'T_all_optical', 'Cmma_103', 'Cmma_106', 'Cmma_114', 'Cmma_127')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    with Path(path).open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(obj, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write('\n')


def training_rows(network):
    if network['status'] != 'complete_ordinary_G2_sampled_network_not_TS_or_error_certified':
        raise ValueError('only the explicitly sampled eight-band target is supported')
    rows = network['rows']
    keys = {(r['substrate_strain'], r['candidate']) for r in rows}
    if len(rows) != 8 or keys != {(s, c) for s in (0., .01) for c in ('T', 'M', 'preserving', 'reversing')}:
        raise ValueError('complete eight-band training required')
    for r in rows:
        values = np.asarray(r['raw_relative_energy_meV_fu'], dtype=float)
        if (values.shape != (9,) or not np.isfinite(values).all()
                or abs(values[0]) > 1e-8 or r['ordinary_fmax_eV_A'] > .10
                or r['ordinary_fmax_eV_A'] < 0
                or not np.isfinite(r['ordinary_fmax_eV_A'])
                or abs(values.max() - r['forward_meV_fu']) > 1e-8
                or abs(values[-1] - r['final_minus_initial_meV_fu']) > 1e-8
                or r['barrier_error_bound_meV_fu'] is not None):
            raise ValueError('sampled profile, residual or unresolved-error semantics changed')
    return {(r['substrate_strain'], r['candidate']): r for r in rows}


def sampled_selection(network):
    rows = training_rows(network)
    roles = {'switching': ('preserving', 'reversing'), 'decay': ('T', 'M')}
    selected = {}
    for role, names in roles.items():
        name = min(names, key=lambda n: (rows[0., n]['forward_meV_fu'], n))
        selected[role] = name
    return selected


def simple_values(anchor_peak, upper_peak, initial_shift=None, final_shift=None,
                  reaction=None):
    """All inputs eV/f.u.; retain signed invalid predictions as abstentions."""
    required = [anchor_peak, upper_peak]
    optional = [v for v in (initial_shift, final_shift, reaction) if v is not None]
    if not np.isfinite(required + optional).all():
        raise ValueError('finite energy inputs required')
    values = {'B0': .5 * (anchor_peak + upper_peak),
              'B1_fixed': None if initial_shift is None else anchor_peak - initial_shift,
              'B1_follow': None if initial_shift is None or final_shift is None
              else anchor_peak + final_shift - initial_shift}
    result = {}
    for name, raw in values.items():
        invalid = raw is not None and raw < max(0., reaction or 0.) - 1e-12
        reason = ('visible_endpoint_features_pending' if raw is None else
                  'below_endpoint_inclusive_peak_limit' if invalid else None)
        result[name] = {'raw_signed_prediction_eV_fu': raw,
                        'prediction_eV_fu': None if invalid else raw,
                        'abstention_reason': reason, 'physical_error_bound_eV_fu': None}
    return result


def register(repository, output):
    if output.exists():
        raise FileExistsError(output)
    source = repository / NETWORK
    if sha(source) != NETWORK_SHA:
        raise ValueError('pinned E148 training report changed')
    network = json.loads(source.read_bytes())
    rows, selected = training_rows(network), sampled_selection(network)
    if selected != {'switching': 'preserving', 'decay': 'T'}:
        raise ValueError('prospective selected pair changed')
    predictions = {}
    for name in selected.values():
        predictions[CHANNELS[name]] = simple_values(rows[0., name]['forward_meV_fu']/1000,
                                                    rows[.01, name]['forward_meV_fu']/1000)
    report = dict(format_version=1, status='pre_endpoint_registration_not_final_numeric_panel',
        target='endpoint_inclusive_nine_image_ordinary_sampled_peak', formula_units=4,
        target_strain=.005, training_strains=[0., .01], pressure_GPa=0., electric_field=0.,
        training_network_sha256=NETWORK_SHA, protocol_sha256=sha(repository/CASE/'protocol.md'),
        source_sha256=sha(Path(__file__)), selected=selected, predictions=predictions,
        physical_input_sha256=CONTRACT, endpoint_phases=list(PHASES),
        endpoint_shards=[list(PHASES[:2]), [PHASES[2]]],
        endpoint_recipe=dict(optimizer='BFGS', fmax=.03, open_traction_kbar=2.,
                             maxstep=.02, steps=80, calls_per_endpoint_max=81),
        path_recipe=dict(total_images=9, interior_images=7, optimizer='FIRE',
                         fmax=.10, k=.2, maxstep=.02, steps=80, maximum_cell_step=.02,
                         minimum_distance=1.6, maximum_deformation=.25,
                         candidate_step_retries=4, climb=False),
        unprovided_physical_transfers={
            name: dict(representations=list(REPRESENTATIONS), curvature_steps_A=[.01, .02],
                       prediction_eV_fu=None,
                       reason='no validated physical-epsilon continuation of the complete sampled path; '
                              'local centre is not a certified stationary branch')
            for name in ('B2_frozen', 'B3_atomic', 'B4_joint', 'B5_branch')},
        retained_nonlinear_model=dict(parameters=109, steps_A=[.01, .02],
            prediction_eV_fu=None, reason='zero-strain internal-Q fit has no validated physical-epsilon transfer'),
        original_2D_preserving_refusal_retained=True,
        strict_error_certified_selector_unchanged=True,
        physical_ordering_certified=False, physical_advantage_established=False,
        predeclared_visible_features='three audited target endpoints equally available to all arms; count every SCF',
        final_freeze_gate='substitute fixed B1 formulas and commit all arms before target-interior SCF',
        target_structures_generated=False, target_DFT_calls=0, target_path_labels_read=False)
    save(output, report)
    return report


def prepare_endpoints(repository, registration, expected_sha, output):
    if output.exists():
        raise FileExistsError(output)
    if sha(registration) != expected_sha:
        raise ValueError('external pre-endpoint registration pin changed')
    registered = json.loads(registration.read_bytes())
    if (registered['status'] != 'pre_endpoint_registration_not_final_numeric_panel'
            or registered['protocol_sha256'] != sha(repository/CASE/'protocol.md')
            or registered['physical_input_sha256'] != CONTRACT
            or registered['selected'] != {'switching': 'preserving', 'decay': 'T'}):
        raise ValueError('fixed registered physics/selected pair required')
    prepared = []
    for phase in PHASES:
        _, atoms, _, original, evidence = endpoint_evidence(repository/BASE, 'strain_0000', phase)
        reference = np.asarray(original['reference_cell_A']).copy()
        reference[:2] *= 1.005
        cell = atoms.cell.array.copy()
        cell[:2] = reference[:2]
        atoms.set_cell(cell, scale_atoms=True)
        atoms.calc = None
        boundary = clamped_plane_vcneb_boundary(12, reference, allow_tilt=True)
        boundary.validate_images([atoms]); geometry_guard([atoms])
        prepared.append((phase, atoms, original, evidence, reference))
    output.mkdir(parents=True, exist_ok=False)
    entries = []
    for phase, atoms, original, evidence, reference in prepared:
        folder = output/phase; folder.mkdir()
        write_clean_poscar(folder/'POSCAR.seed', atoms)
        manifest = dict(format_version=1, status='geometry_seed_not_relaxed', phase_label=phase,
            strain=.005, seed_file='POSCAR.seed', seed_sha256=sha(folder/'POSCAR.seed'),
            ordered_species=atoms.get_chemical_symbols(), reference_cell_A=reference.tolist(),
            allow_tilt=True, cell_scale_A=original['cell_scale_A'], pressure_gpa=0.,
            minimum_distance_guard_A=1.6, physical_contract_sha256=CONTRACT,
            source_zero_strain_evidence=evidence, ordered_geometry=endpoint_structure_record(atoms),
            registration_sha256=expected_sha, new_DFT_calls=0,
            old_endpoint_results_reused=False, target_interior_generated=False)
        save(folder/'endpoint_seed.json', manifest)
        load_seed(folder/'endpoint_seed.json')
        entries.append(dict(phase=phase, manifest_sha256=sha(folder/'endpoint_seed.json'),
                            seed_sha256=manifest['seed_sha256']))
    result = dict(status='three_visible_target_endpoint_seeds_only', registration_sha256=expected_sha,
                  entries=entries, target_strain=.005, new_DFT_calls=0,
                  target_interiors_generated=False, complete_forecast_panel_frozen=False)
    save(output/'seeds.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('register', 'endpoints'))
    p.add_argument('--repository', type=Path, default=Path.cwd())
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--registration', type=Path)
    p.add_argument('--registration-sha256')
    a = p.parse_args()
    if a.mode == 'register':
        result = register(a.repository, a.output)
    else:
        result = prepare_endpoints(a.repository, a.registration, a.registration_sha256, a.output)
    print(json.dumps({'status': result['status'], 'new_DFT_calls': 0}))
