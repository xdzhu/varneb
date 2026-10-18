"""Prepare two upper-to-midpoint strain continuations; no DFT or submission.

Reuse the ordinary-converged +1% polygons and the very same E268 target
endpoints. This is an initial-path sensitivity check, not a new holdout.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from ase.io import read, write

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.analyze_hfo2_clamped_reference import terminal_images
from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_clamped_chains import endpoint_evidence
from scripts.prepare_hfo2_reference_variants import write_clean_poscar
from scripts.relax_clamped_ase_endpoint import load_seed

BASE = Path('benchmarks/hfo2_channels/20261008')
OLD = BASE/'physical_holdout_E268_20261018'
CASE = BASE/'midpoint_reverse_E269_20261018'
FEATURE_PIN = '0441a5e7e63ec86a68807068ea21db6901f4a762d740592473d91bf15bbd8f23'
SOURCES = {
    'PO_flip_T_pattern_preserving': (
        'clamped_plus1_preserving_terminal_E145_20261012/completed_HF',
        '0c8d91532bc9051decdfd600eae12ce409fd84a86c677158077948fa352c9d8a',
        'PO_minus_T_preserving'),
    'PO_to_T': ('clamped_plus1_T_terminal_E087_20261011',
        '32a6e49a618e6e62632ae71a5e0ed093095254f6858221be2ade3d83cc194a28', 'T'),
}


def save(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def historical_preparer(repository):
    spec = importlib.util.spec_from_file_location('E268_geometry_recipe',
        repository/OLD/'forecast_and_chain.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def transfer_endpoint_lift(source, target, zero_anchor):
    """Register small endpoint motion to the original parent lift, atom by atom.

    No species reassignment or image-wise gauge fitting. The target endpoint
    is already in the zero-anchor gauge. Large/ambiguous motion is rejected.
    """
    if not (source.get_chemical_symbols() == target.get_chemical_symbols()
            == zero_anchor.get_chemical_symbols()):
        raise ValueError('ordered endpoint species changed')
    delta = source.get_scaled_positions(wrap=False)-zero_anchor.get_scaled_positions(wrap=False)
    lift = np.rint(delta).astype(int)
    residual = delta-lift
    if not np.isfinite(residual).all() or np.max(np.abs(residual)) >= .2:
        raise ValueError('endpoint motion left the declared small-strain parent registration')
    result = target.copy()
    result.calc = None
    result.set_scaled_positions(target.get_scaled_positions(wrap=False)+lift)
    return result, dict(integer_lift=lift.tolist(),
        maximum_fractional_parent_displacement=float(np.max(np.abs(residual))))


def prepare(repository, output):
    repository, output = Path(repository), Path(output)
    if output.exists():
        raise FileExistsError(output)
    features_path = repository/OLD/'endpoint_features/features.json'
    if sha256(features_path) != FEATURE_PIN:
        raise ValueError('the same completed E268 endpoint features are required')
    features = json.loads(features_path.read_bytes())
    if features['physical_input_sha256'] != CONTRACT or features['target_strain'] != .005:
        raise ValueError('unchanged input and +0.5% boundary required')
    recipe = historical_preparer(repository)
    _, boundary, endpoint_manifest = load_seed(repository/OLD/'endpoint_seeds/PO_plus/endpoint_seed.json')
    records, prepared = [], []
    # Complete every input check before creating the output namespace.
    for channel, (case, pin, final_phase) in SOURCES.items():
        images, observation, residual, receipt = terminal_images(repository/BASE/case, pin)
        upper_plane = boundary.reference_cell[:2]*1.01/1.005
        if (residual['fmax_eV_A'] > .10 or observation['pressure_GPa'] != 0
                or not all(np.allclose(a.cell.array[:2], upper_plane, atol=1e-10, rtol=0)
                           for a in images)):
            raise ValueError('ordinary-converged +1% polygon required')
        endpoints, lifts = [], []
        for source, phase in zip((images[0], images[-1]), ('PO_plus', final_phase)):
            path = features_path.parent/f'{phase}.vasp'
            if sha256(path) != features['geometry_sha256'][phase]:
                raise ValueError('original target endpoint changed')
            target = read(path, format='vasp')
            zero_anchor = endpoint_evidence(repository/BASE, 'strain_0000', phase)[1]
            endpoint, lift = transfer_endpoint_lift(source, target, zero_anchor)
            endpoints.append(endpoint); lifts.append(lift)
        new = recipe.deform_zero_band(images, endpoints, boundary)
        caches = [features['endpoints']['PO_plus']['cache']]+[None]*7+[features['endpoints'][final_phase]['cache']]
        records.append(dict(channel=channel, source_case=case, source_job=receipt['job_id'],
            source_step=receipt['terminal_step'], source_terminal_audit_sha256=pin,
            endpoint_lift=lifts, source_strain=.01, target_strain=.005,
            physical_input_sha256=CONTRACT, endpoint_features_sha256=FEATURE_PIN,
            cell_scale_A=endpoint_manifest['cell_scale_A'], endpoint_cache_records=caches,
            new_DFT_calls=0, independent_forecast_test=False))
        prepared.append(new)
    output.mkdir(parents=True)
    for record, images in zip(records, prepared):
        folder = output/record['channel']; folder.mkdir()
        for name, atoms in (('initial.vasp', images[0]), ('final.vasp', images[-1]),
                            ('substrate.vasp', images[0])):
            write_clean_poscar(folder/name, atoms)
        write(folder/'seed.traj', images)
        old_parameters = repository/OLD/'path_seeds'/record['channel']/'factory_parameters.json'
        parameters = json.loads(old_parameters.read_bytes())
        if parameters['seed_cache_records'] != record['endpoint_cache_records']:
            raise ValueError('new and old midpoint chains must reuse identical endpoint records')
        save(folder/'factory_parameters.json', parameters)
        record['files_sha256'] = {n: sha256(folder/n) for n in
            ('initial.vasp', 'final.vasp', 'substrate.vasp', 'seed.traj', 'factory_parameters.json')}
        save(folder/'manifest.json', record)
    result = dict(status='upper_to_midpoint_two_seeds_no_DFT',
        recipe_source_sha256=sha256(repository/OLD/'forecast_and_chain.py'),
        preparation_source_sha256=sha256(Path(__file__)),
        endpoint_features_sha256=FEATURE_PIN, seeds=records,
        new_DFT_calls=0, original_E268_labels_unchanged=True, independent_forecast_test=False)
    save(output/'preparation.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, default=Path.cwd())
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.repository, args.output)
    print(json.dumps({k: v for k, v in report.items() if k != 'seeds'}))
