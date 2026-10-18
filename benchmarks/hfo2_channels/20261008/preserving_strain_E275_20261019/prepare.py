"""Transfer one actual lower band to two existing strain boundaries; zero DFT."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import tempfile

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
import numpy as np
from ase.io import read, write
from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts.audit_hfo2_static_replica import sha256
from scripts.plot_hfo2_preserving_connection import restore, native_images
from scripts.prepare_hfo2_midpoint_reverse import historical_preparer, transfer_endpoint_lift, FEATURE_PIN
from scripts.prepare_hfo2_clamped_chains import endpoint_evidence
from scripts.analyze_hfo2_clamped_reference import terminal_images
from scripts.prepare_hfo2_reference_variants import write_clean_poscar
from vcneb import clamped_plane_vcneb_boundary

CASE = Path(__file__).resolve().parent
BASE = CASE.parent
CHANNEL = 'PO_flip_T_pattern_preserving'
UPPER = BASE/'clamped_plus1_preserving_terminal_E145_20261012/completed_HF'
UPPER_PIN = '0c8d91532bc9051decdfd600eae12ce409fd84a86c677158077948fa352c9d8a'
SOURCE = BASE/'preserving_connection_E273_20261019'
SOURCE_ARCHIVE_PIN = '54bda58fa6ea2e2bf1bbc6ae1ef193379ebb90aa5a386f699e6546b0e6056d10'
ASSETS = ('/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/validation/hfo2_t_to_po/'
          'vcneb_n7_fire_distributed_resume27678004_v4/00')


def save(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as h:
        json.dump(value, h, indent=2, allow_nan=False); h.write('\n')


def prepare(output):
    output = Path(output)
    if output.exists(): raise FileExistsError(output)
    if sha256(SOURCE/'native_terminal.tar.gz') != SOURCE_ARCHIVE_PIN:
        raise ValueError('the complete actual E273 terminal archive is required')
    feature_path = BASE/'physical_holdout_E268_20261018/endpoint_features/features.json'
    if sha256(feature_path) != FEATURE_PIN: raise ValueError('original midpoint endpoints changed')
    features = json.loads(feature_path.read_bytes())
    if features['physical_input_sha256'] != CONTRACT: raise ValueError('physics changed')
    upper, _, upper_residual, upper_receipt = terminal_images(UPPER, UPPER_PIN)
    if upper_residual['fmax_eV_A'] > .10: raise ValueError('original upper band must be terminal')
    upper_parameters = json.loads((UPPER/'vcneb_summary.json').read_bytes())['calculator_parameters']
    recipe = historical_preparer(REPO)
    records = []
    with tempfile.TemporaryDirectory(prefix='E275-source-') as td:
        folder, _ = restore(9, Path(td), terminal=True)
        source, observation, receipt = native_images(folder, 9)
        if receipt['job_id'] != '28949435' or len(source) != 9:
            raise ValueError('only the actual lower zero-strain source is admissible')
        output.mkdir(parents=True)
        for label, strain, tag in (('strain_0050', .005, 'strain_0000'),
                                   ('strain_0100', .01, 'strain_p0100')):
            endpoints, caches, lifts = [], [], []
            for k, phase in enumerate(('PO_plus', 'PO_minus_T_preserving')):
                if strain == .005:
                    path = feature_path.parent/f'{phase}.vasp'
                    if sha256(path) != features['geometry_sha256'][phase]:
                        raise ValueError('midpoint endpoint geometry changed')
                    target = read(path)
                    cache = features['endpoints'][phase]['cache']
                else:
                    target = endpoint_evidence(BASE, tag, phase)[1]
                    if not same_ordered_geometry(target, upper[0 if k == 0 else -1]):
                        raise ValueError('upper endpoint source is not the existing path endpoint')
                    cache = upper_parameters['seed_cache_records'][0 if k == 0 else -1]
                anchor = endpoint_evidence(BASE, 'strain_0000', phase)[1]
                endpoint, lift = transfer_endpoint_lift(source[0 if k == 0 else -1], target, anchor)
                endpoints.append(endpoint); caches.append(cache); lifts.append(lift)
            boundary = clamped_plane_vcneb_boundary(12, endpoints[0].cell.array, allow_tilt=True)
            if not np.allclose(boundary.reference_cell[:2], source[0].cell.array[:2]*(1+strain),
                               atol=1e-10, rtol=0):
                raise ValueError('target must be the same T-referenced substrate boundary')
            images = recipe.deform_zero_band(source, endpoints, boundary)
            for image in images: image.calc = None
            seed = output/label; seed.mkdir()
            for name, atoms in (('initial.vasp', images[0]), ('final.vasp', images[-1]),
                                ('substrate.vasp', images[0])):
                write_clean_poscar(seed/name, atoms)
            write(seed/'seed.traj', images)
            save(seed/'factory_parameters.json', dict(source_directory=ASSETS,
                seed_cache_records=[caches[0]]+[None]*7+[caches[1]]))
            row = dict(label=label, target_strain=strain, source_strain=0., source_job='28949435',
                source_step=9, source_archive_sha256=SOURCE_ARCHIVE_PIN,
                source_observation_sha256=sha256(folder/'observation/observation.json'),
                source_summary_sha256=sha256(folder/'vcneb_summary.json'),
                physical_input_sha256=CONTRACT, endpoint_lift=lifts,
                target_endpoint_source='E268' if strain == .005 else 'E145',
                target_original_job='28938772' if strain == .005 else upper_receipt['job_id'],
                n_total_images=9, free_interiors=list(range(1, 8)), new_DFT_calls=0,
                original_labels_and_scores_unchanged=True, independent_forecast_test=False,
                cell_scale_A=5.12968067458423, retained_pattern_constraints_in_NEB=[],
                files_sha256={n: sha256(seed/n) for n in ('initial.vasp', 'final.vasp',
                    'substrate.vasp', 'seed.traj', 'factory_parameters.json')})
            save(seed/'manifest.json', row); records.append(row)
    result = dict(status='two_lower_candidate_strain_seeds_no_DFT', seeds=records,
        new_DFT_calls=0, scheduler_mutations=0, maximum_updates=60,
        maximum_new_interior_SCFs_per_chain=427, walltime='18:00:00',
        protocol_sha256=sha256(CASE/'protocol.md'), preparation_source_sha256=sha256(Path(__file__)),
        deformation_recipe_sha256=sha256(BASE/'physical_holdout_E268_20261018/forecast_and_chain.py'))
    save(output/'preparation.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    print(json.dumps(prepare(p.parse_args().output)))
