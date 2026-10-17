"""Audited endpoint features, complete forecast freeze, then geometric seeds.

No calculator, SSH or scheduler call. The three operations have distinct fresh
namespaces; only the last operation may generate target interior structures.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np
from ase.io import read, write

from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.audit_hfo2_clamped_canary import replay_directories
from scripts.audit_hfo2_clamped_endpoint_variants import geometric_patterns
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.analyze_hfo2_clamped_reference import terminal_images
from scripts.prepare_hfo2_clamped_endpoints import AXES, geometry_guard, orient_long_axis_x
from scripts.prepare_hfo2_clamped_chains import HF_PHYSICS
from scripts.prepare_hfo2_reference_variants import write_clean_poscar
from scripts.relax_clamped_ase_endpoint import load_seed
from vcneb.forecast_freeze import freeze_forecast_panel, load_frozen_forecast_panel
from vcneb.periodic_path import minimum_image_path_lift, validate_periodic_path_lift

spec = importlib.util.spec_from_file_location('E268_registration', Path(__file__).with_name('study.py'))
study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)
sha, save = study.sha, study.save
REGISTRATION_PIN = 'faebe7588579f5dca9253e672755d4db5c85f8a33165685dc5cdce3cef2899ff'
REMOTE = '/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1'
CHANNEL_PHASES = {'PO_flip_T_pattern_preserving': 'PO_minus_T_preserving', 'PO_to_T': 'T'}
SOURCES = {
    'PO_flip_T_pattern_preserving': ('clamped_held_flips_E100_20261011/zero_strain_preserving',
        '46145da14f24b2b8bd46aa4656f0982650d6b2dc661282f58d809aab8a6d3202'),
    'PO_to_T': ('clamped_T_terminal_E079_20261011',
        'c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e'),
}


def registration(repository):
    case = repository/study.CASE
    if sha(case/'registration.json') != REGISTRATION_PIN:
        raise ValueError('original external registration pin changed')
    r = json.loads((case/'registration.json').read_bytes())
    if r['protocol_sha256'] != sha(case/'protocol.md') or sha(repository/study.NETWORK) != study.NETWORK_SHA:
        raise ValueError('original protocol or training changed')
    return r


def features(repository, exported, export_report, output):
    registration(repository)
    if output.exists(): raise FileExistsError(output)
    report = json.loads(export_report.read_bytes())
    if report['target_interiors_read'] or set(report['scheduler']) != {'28938150', '28938151'}:
        raise ValueError('only original three endpoint feature export allowed')
    for name, digest in report['exported_file_sha256'].items():
        path = Path(name)
        if path.is_absolute() or '..' in path.parts or sha(exported/path) != digest:
            raise ValueError('exported native file drift')
    for row in report['scheduler'].values():
        if row[1:6] != ['COMPLETED', '0:0', 'hfacnormal01', '32', '1']:
            raise ValueError('original one-node32MPI endpoint jobs not completed')
    network = study.training_rows(json.loads((repository/study.NETWORK).read_bytes()))
    t = read(repository/study.BASE/'reference_variants/T.vasp', format='vasp')
    old_parent, old_basis, _ = rotated_t_triplet(t)
    parent, basis = orient_long_axis_x(old_parent), old_basis[..., list(AXES)]
    entries, structures = {}, {}
    for phase in study.PHASES:
        seed_path = repository/study.CASE/'endpoint_seeds'/phase/'endpoint_seed.json'
        seed, boundary, m = load_seed(seed_path)
        folder = exported/'endpoints'/phase
        native = json.loads((folder/'native_audit.json').read_bytes())
        summary_path = folder/'endpoint_relax_summary.json'
        summary = json.loads(summary_path.read_bytes())
        calls = sorted((folder/'calculator/image_0000').glob('scf_*'))
        if (native['summary_sha256'] != sha(summary_path) or native['seed_manifest_sha256'] != sha(seed_path)
                or not native['endpoint_converged'] or not native['raw_full_physical_bytes_checked_here']
                or native['new_SCF_calls'] != len(calls) or len(calls) != summary['optimizer_steps']+1
                or len(calls) > 81 or summary['strain'] != .005 or not summary['converged']
                or summary['physical_contract_sha256'] != CONTRACT
                or summary['fmax_target_eV_per_A'] != .03 or summary['open_stress_target_kbar'] != 2.
                or summary['maxstep'] != .02 or summary['steps_requested'] != 80):
            raise ValueError('endpoint source/recipe/convergence changed')
        rows = replay_directories(calls, boundary, full_physical_bytes=False)
        last = rows[-1]
        if last['max_atomic_force_eV_A'] >= .03 or last['open_traction_norm_kbar'] >= 2.:
            raise ValueError('physical endpoint thresholds failed')
        if len(rows) != len(native['rows']) or any(a['input_sha256'] != b['input_sha256']
                or a['raw_log_sha256'] != b['raw_log_sha256']
                or abs(a['energy_eV_cell']-b['energy_eV_cell']) > 1e-9
                for a, b in zip(rows, native['rows'])):
            raise ValueError('native E/F/stress audit does not replay')
        lifted, lift = minimum_image_path_lift([seed] + [read_fixed_hfo2_stru(p/'STRU') for p in calls])
        atoms = lifted[-1]; validate_periodic_path_lift(lifted)
        if not same_ordered_geometry(atoms, read(folder/'CONTCAR', format='vasp')):
            raise ValueError('continuous endpoint lift changed geometry')
        patterns = geometric_patterns(atoms, parent, basis)
        labels = [x['number'] for x in last['structure_audit']['symmetry_sweep']]
        if not patterns['nearest_reference_chart_accepted']:
            raise ValueError('endpoint left the registered geometric reference chart')
        polar = np.asarray(patterns['oxygen_minus_hafnium_mean_displacement_A'])
        identity = (set(labels) <= {68, 137} and np.linalg.norm(polar) < .03) if phase == 'T' else (
            labels == [29, 29, 29] and patterns['Q_A'][0] < -.1
            and (polar[0] > .05 if phase == 'PO_plus' else polar[0] < -.05))
        if not identity:
            raise ValueError(f'endpoint phase/registered variant changed: {phase}; report, do not force symmetry')
        anchor = network[0., 'T']['matched_initial_energy_eV_cell']
        if phase != 'PO_plus':
            channel = 'T' if phase == 'T' else 'preserving'
            anchor += network[0., channel]['final_minus_initial_meV_fu']*4/1000
        entries[phase] = dict(anchor_energy_eV_cell=anchor, target_energy_eV_cell=last['energy_eV_cell'],
            native_audit_sha256=sha(folder/'native_audit.json'), summary_sha256=sha(summary_path),
            new_SCF_calls=len(rows), SCF_core_hours=sum(r['elapsed_seconds'] for r in rows)*32/3600,
            physical_input_sha256=CONTRACT, patterns=patterns,
            three_tolerance_phase_numbers=labels, phase_variant_screen_passed=True,
            full_Hessian_stability_certified=False, raw_full_physical_bytes_rechecked_offline=False,
            original_HF_full_physical_audit_preserved=True, continuous_lift=lift,
            cache=native['cache'], max_atomic_force_eV_A=last['max_atomic_force_eV_A'],
            open_traction_norm_kbar=last['open_traction_norm_kbar'])
        structures[phase] = atoms
    output.mkdir()
    for phase, atoms in structures.items(): write_clean_poscar(output/f'{phase}.vasp', atoms)
    result = dict(status='three_native_visible_endpoint_features_audited', registration_sha256=REGISTRATION_PIN,
        export_report_sha256=sha(export_report), endpoints=entries, target_strain=.005, formula_units=4,
        physical_input_sha256=CONTRACT, new_SCF_calls=sum(e['new_SCF_calls'] for e in entries.values()),
        SCF_core_hours=sum(e['SCF_core_hours'] for e in entries.values()),
        allocation_core_hours=report['allocation_core_hours'],
        geometry_sha256={phase: sha(output/f'{phase}.vasp') for phase in structures},
        target_interiors_read=False, all_models_equal_access=True)
    save(output/'features.json', result)
    return result


def model_definitions():
    definitions = {
        'B0': dict(formula='0.5*(sampled_peak_0+sampled_peak_1percent)', new_target_feature_calls_required=0),
        'B1_fixed': dict(formula='sampled_peak_0-delta_E_initial_per_fu', target_endpoint_inputs_required=True),
        'B1_follow': dict(formula='sampled_peak_0+delta_E_final_per_fu-delta_E_initial_per_fu', target_endpoint_inputs_required=True),
    }
    reason = ('No validated continuation of the complete sampled path in physical epsilon; '
              'local chart is not a certified stationary branch. Not a measured inaccurate forecast.')
    for arm in ('B2_frozen', 'B3_atomic', 'B4_joint', 'B5_branch'):
        for rep in study.REPRESENTATIONS:
            for step in (.01, .02):
                definitions[f'{arm}/{rep}/h{step}'] = dict(handling=arm, representation=rep,
                    curvature_step_A=step, physical_prediction_unavailable_reason=reason,
                    local_preserving_original_2D_negative_block_refusal_retained=True)
    for step in (.01, .02):
        definitions[f'B5_retained_anharmonic/h{step}'] = dict(parameters=109, curvature_step_A=step,
            physical_prediction_unavailable_reason='P8 model fits internal Q at zero physical strain; '
                'no physical-epsilon transfer was fitted or validated; P8 labels are not reused to add one')
    return definitions


def freeze(repository, feature_path, feature_pin, output):
    registration(repository)
    if output.exists(): raise FileExistsError(output)
    if sha(feature_path) != feature_pin: raise ValueError('externally audited endpoint feature pin changed')
    f = json.loads(feature_path.read_bytes())
    if (f['status'] != 'three_native_visible_endpoint_features_audited'
            or f['registration_sha256'] != REGISTRATION_PIN or f['target_interiors_read']
            or f['physical_input_sha256'] != CONTRACT or set(f['endpoints']) != set(study.PHASES)):
        raise ValueError('only three audited, predeclared endpoint inputs')
    rows = study.training_rows(json.loads((repository/study.NETWORK).read_bytes()))
    initial = f['endpoints']['PO_plus']
    di = (initial['target_energy_eV_cell']-initial['anchor_energy_eV_cell'])/4
    values = {}
    for key, channel in study.CHANNELS.items():
        endpoint = f['endpoints'][CHANNEL_PHASES[channel]]
        df = (endpoint['target_energy_eV_cell']-endpoint['anchor_energy_eV_cell'])/4
        reaction = (endpoint['target_energy_eV_cell']-initial['target_energy_eV_cell'])/4
        values[channel] = study.simple_values(rows[0., key]['forward_meV_fu']/1000,
            rows[.01, key]['forward_meV_fu']/1000, di, df, reaction)
    output.mkdir()
    cost = dict(visible_endpoint_SCFs=f['new_SCF_calls'], endpoint_SCF_core_hours=f['SCF_core_hours'],
        endpoint_allocation_core_hours=f['allocation_core_hours'],
        all_prior_training_curvature_faces_and_P8_costs_additional=True,
        prior_costs_not_zero_or_subtracted_for_reuse=True,
        prior_complete_cost_total_not_recomputed_by_this_freezer=True,
        B0_requires_no_new_endpoint_feature_calls=True,
        B1_requires_audited_endpoints_shared_with_ground_truth=True)
    save(output/'cost.json', cost)
    contract = dict(target='nine_image_ordinary_sampled_peak_not_TS', formula_units=4, energy_unit='eV/formula_unit',
        physical_input_sha256=CONTRACT, pressure_GPa=0., strain=.005, allow_tilt=True,
        protocol_sha256=sha(repository/study.CASE/'protocol.md'), registration_sha256=REGISTRATION_PIN,
        physical_error_bounds=None, reference_resolution_certified=False,
        all_target_endpoints_visible_to_every_arm=True, target_path_labels_read=False,
        endpoint_features_sha256=feature_pin, complete_candidate_network_claimed=False)
    save(output/'contract.json', contract)
    inputs = {
        'training_eight_bands': dict(path=repository/study.NETWORK, kind='training'),
        'visible_endpoints': dict(path=feature_path, kind='visible_endpoint'),
        'protocol': dict(path=repository/study.CASE/'protocol.md', kind='protocol'),
        'pre_endpoint_registration': dict(path=repository/study.CASE/'registration.json', kind='protocol'),
        'source': dict(path=Path(__file__), kind='source'),
        'registered_arithmetic_source': dict(path=Path(__file__).with_name('study.py'), kind='source'),
        'preparation_cost': dict(path=output/'cost.json', kind='preparation_cost'),
    }
    hashes = {k: sha(v['path']) for k, v in inputs.items()}
    definitions = model_definitions(); artifacts = {}
    for i, (name, definition) in enumerate(definitions.items()):
        predictions = []
        for channel in CHANNEL_PHASES:
            row = values[channel][name] if name in ('B0', 'B1_fixed', 'B1_follow') else dict(
                prediction_eV_fu=None, abstention_reason=definition['physical_prediction_unavailable_reason'])
            predictions.append(dict(case_id=channel, prediction_eV_fu=row['prediction_eV_fu'],
                                    abstention_reason=row['abstention_reason']))
        if name in ('B0', 'B1_fixed', 'B1_follow'):
            definition['raw_signed_predictions'] = {c: v[name]['raw_signed_prediction_eV_fu'] for c, v in values.items()}
        artifacts[name] = output/f'model_{i:03d}.json'
        save(artifacts[name], dict(format_version=1, energy_unit='eV/formula_unit', model_id=name,
            comparison_contract_sha256=sha(output/'contract.json'), model_definition=definition,
            input_artifact_sha256=hashes, forecasts=predictions))
    candidate = 'B5_branch/path_adapted/h0.01'
    save(output/'panel_registration.json', dict(format_version=1, energy_unit='eV/formula_unit',
        comparison_contract_sha256=sha(output/'contract.json'), candidate_model=candidate,
        baseline_models=[name for name in definitions if name != candidate], case_ids=list(CHANNEL_PHASES)))
    receipt = freeze_forecast_panel(output/'panel_registration.json', output/'contract.json',
        model_artifact_paths=artifacts, input_artifacts=inputs, output=output/'frozen', holdout_path_labels_unread=True)
    actual = load_frozen_forecast_panel(output/'frozen', expected_manifest_sha256=receipt['manifest_sha256'])
    result = dict(status='complete_physical_sampled_peak_panel_frozen', manifest_sha256=receipt['manifest_sha256'],
        models=len(definitions), forecast_records=len(actual), numeric_records=sum(r.prediction_eV_fu is not None for r in actual),
        baselines=values, new_DFT_calls=0, target_interiors_read=False,
        unavailable_transfer_records_are_not_inaccurate_models=True)
    save(output/'panel_receipt.json', result)
    return result


def deform_zero_band(images, endpoints, boundary):
    if len(images) != 9: raise ValueError('nine training images required')
    validate_periodic_path_lift(images)
    result = [a.copy() for a in images]
    for atoms in result:
        cell = atoms.cell.array.copy(); cell[:2] = boundary.reference_cell[:2]
        atoms.set_cell(cell, scale_atoms=True); atoms.calc = None
    corrections = [end.get_scaled_positions(wrap=False)-a.get_scaled_positions(wrap=False)
                   for a, end in zip((result[0], result[-1]), endpoints)]
    cells = [end.cell.array-a.cell.array for a, end in zip((result[0], result[-1]), endpoints)]
    for i, atoms in enumerate(result):
        t = i/8
        frac = atoms.get_scaled_positions(wrap=False)+(1-t)*corrections[0]+t*corrections[1]
        cell = atoms.cell.array+(1-t)*cells[0]+t*cells[1]
        cell[:2] = boundary.reference_cell[:2]
        atoms.set_cell(cell, scale_atoms=False); atoms.set_scaled_positions(frac)
    boundary.validate_images(result); geometry_guard(result); validate_periodic_path_lift(result)
    if not all(same_ordered_geometry(a, b) for a, b in zip((result[0], result[-1]), endpoints)):
        raise ValueError('target endpoint geometry drift')
    return result


def chains(repository, feature_path, panel, panel_pin, output):
    registration(repository)
    if output.exists(): raise FileExistsError(output)
    frozen = load_frozen_forecast_panel(panel, expected_manifest_sha256=panel_pin)
    expected = set(model_definitions())
    if len(frozen) != 138 or {r.model_id for r in frozen} != expected:
        raise ValueError('all69 declared model/variant records across two cases required')
    contract = json.loads((panel/'comparison_contract.json').read_bytes())
    if contract['endpoint_features_sha256'] != sha(feature_path):
        raise ValueError('frozen endpoint features changed')
    f = json.loads(feature_path.read_bytes()); output.mkdir()
    all_seeds = []
    for channel, phase in CHANNEL_PHASES.items():
        source, pin = SOURCES[channel]
        images, observation, residual, receipt = terminal_images(repository/study.BASE/source, pin)
        if residual['fmax_eV_A'] > .10 or observation['pressure_GPa'] != 0.:
            raise ValueError('audited ordinary zero-strain source required')
        endpoints = []
        for p, original in zip(('PO_plus', phase), (images[0], images[-1])):
            geometry = feature_path.parent/f'{p}.vasp'
            if sha(geometry) != f['geometry_sha256'][p]: raise ValueError('endpoint feature geometry changed')
            atoms = read(geometry, format='vasp')
            # Preserve the original G2 channel's integer lift. The feature's
            # history was lifted from the same zero-strain endpoint gauge.
            anchor = study.endpoint_evidence(repository/study.BASE, 'strain_0000', p)[1]
            delta = original.get_scaled_positions(wrap=False)-anchor.get_scaled_positions(wrap=False)
            lift = np.rint(delta).astype(int)
            if np.max(np.abs(delta-lift)) > 1e-9: raise ValueError('source endpoint gauge not registered')
            atoms.set_scaled_positions(atoms.get_scaled_positions(wrap=False)+lift)
            endpoints.append(atoms)
        _, boundary, m = load_seed(repository/study.CASE/'endpoint_seeds/PO_plus/endpoint_seed.json')
        new = deform_zero_band(images, endpoints, boundary)
        folder = output/channel; folder.mkdir()
        for name, atoms in [('initial.vasp', new[0]), ('final.vasp', new[-1]), ('substrate.vasp', new[0])]:
            write_clean_poscar(folder/name, atoms)
        write(folder/'seed.traj', new)
        caches = [f['endpoints']['PO_plus']['cache']]+[None]*7+[f['endpoints'][phase]['cache']]
        save(folder/'factory_parameters.json', dict(source_directory=HF_PHYSICS, seed_cache_records=caches))
        manifest = dict(status='frozen_prediction_target_seed_no_interior_SCF', channel=channel, strain=.005,
            forecast_manifest_sha256=panel_pin, endpoint_features_sha256=sha(feature_path),
            source_terminal_audit_sha256=pin, source_job_id=receipt['job_id'],
            source_step=receipt['terminal_step'], recipe=registration(repository)['path_recipe'],
            physical_input_sha256=CONTRACT, cell_scale_A=m['cell_scale_A'],
            reference_cell_A=m['reference_cell_A'], new_DFT_calls=0,
            target_interiors_generated=True, target_interiors_evaluated=False,
            files_sha256={name: sha(folder/name) for name in ('initial.vasp', 'final.vasp', 'substrate.vasp', 'seed.traj', 'factory_parameters.json')})
        save(folder/'manifest.json', manifest)
        all_seeds.append(dict(channel=channel, manifest_sha256=sha(folder/'manifest.json')))
    result = dict(status='two_registered_target_seed_chains_prepared', seeds=all_seeds,
                  forecast_manifest_sha256=panel_pin, new_DFT_calls=0)
    save(output/'chains.json', result); return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('features', 'freeze', 'chains'))
    p.add_argument('--repository', type=Path, default=Path.cwd())
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--exported', type=Path); p.add_argument('--export-report', type=Path)
    p.add_argument('--features', type=Path); p.add_argument('--feature-sha256')
    p.add_argument('--panel', type=Path); p.add_argument('--panel-sha256')
    a = p.parse_args()
    if a.mode == 'features': r = features(a.repository, a.exported, a.export_report, a.output)
    elif a.mode == 'freeze': r = freeze(a.repository, a.features, a.feature_sha256, a.output)
    else: r = chains(a.repository, a.features, a.panel, a.panel_sha256, a.output)
    print(json.dumps({k: v for k, v in r.items() if k not in ('endpoints', 'baselines')}))
