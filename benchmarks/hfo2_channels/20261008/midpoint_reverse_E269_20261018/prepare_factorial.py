"""Prepare four frozen atomic/cell crossover SCFs without launching any.

Each channel has two observed diagonal states (s_A,H_A) and (s_B,H_B).
Only (s_A,H_B) and (s_B,H_A) are new. Fractions are internal coordinates in
one fixed parent chart; holding them under a cell change is an affine strain,
not holding Cartesian positions fixed. No unique mode-energy split follows.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

import numpy as np
from ase.io import read

from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_channel_work_probes import write_structure
from scripts.prepare_hfo2_clamped_chains import HF_PHYSICS
from vcneb import clamped_plane_vcneb_boundary, validate_path_geometry

BASE = Path('benchmarks/hfo2_channels/20261008')
OLD_RESULTS = BASE/'physical_holdout_E268_20261018/completed_paths_20261018'


def save(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def crossover(atom_source, cell_source):
    """Change only H, preserving the declared unwrapped fractional pattern."""
    if atom_source.get_chemical_symbols() != cell_source.get_chemical_symbols():
        raise ValueError('ordered species differ')
    if not np.allclose(atom_source.cell.array[:2], cell_source.cell.array[:2], atol=1e-12, rtol=0):
        raise ValueError('the two states must have the same prescribed substrate plane')
    result = atom_source.copy(); result.calc = None
    result.set_cell(cell_source.cell.array, scale_atoms=True)
    if not np.allclose(result.get_scaled_positions(wrap=False),
                       atom_source.get_scaled_positions(wrap=False), atol=1e-12, rtol=0):
        raise ValueError('cell exchange changed internal fractional coordinates')
    return result


def contrasts(energies):
    """Finite four-state contrasts, eV/cell inputs -> meV/HfO2 outputs."""
    if set(energies) != {'AA', 'AB', 'BA', 'BB'}:
        raise ValueError('all four measured state energies are required')
    if any(isinstance(v, bool) or not isinstance(v, (float, int)) or not np.isfinite(v)
           for v in energies.values()):
        raise ValueError('missing/nonfinite energies cannot be replaced by zero or predictions')
    aa, ab, ba, bb = (energies[k] for k in ('AA', 'AB', 'BA', 'BB'))
    return dict(diagonal_change_meV_fu=(bb-aa)*250,
        internal_at_A_cell_meV_fu=(ba-aa)*250,
        cell_at_A_internal_meV_fu=(ab-aa)*250,
        internal_at_B_cell_meV_fu=(bb-ab)*250,
        cell_at_B_internal_meV_fu=(bb-ba)*250,
        finite_interaction_meV_fu=(bb-ba-ab+aa)*250,
        unique_elastic_or_phonon_energy_partition=False)


def source_peak(folder, *, initial=False):
    receipt = json.loads((folder/('initial_receipt.json' if initial else 'audit.json')).read_bytes())
    path = folder/f"step_{receipt['step']:04d}"/'observation.json'
    observation = json.loads(path.read_bytes())
    index = int(np.argmax(observation['relative_enthalpy_meV_fu']))
    record = observation['raw_image_evaluations'][index]
    poscar = path.parent/f'POSCAR_{index:02d}'
    if sha256(path) != receipt['observation_sha256'] or sha256(poscar) != record['snapshot_POSCAR_sha256']:
        raise ValueError('source sampled peak changed')
    return read(poscar, format='vasp'), record, path.parent/'raw'/f'image_{index:04d}', dict(
        observation_sha256=sha256(path), source_job=observation['source_job_id'], step=receipt['step'],
        sampled_peak_image=index, cached_energy_eV_cell=record['energy_eV_cell'],
        source_static_directory=record['raw_source'], source_input_sha256=record['input_sha256'],
        source_log_sha256=record['raw_log_sha256'])


def prepare(repository, exported, output, *, initial=False):
    if output.exists():
        raise FileExistsError(output)
    spec = importlib.util.spec_from_file_location('E269_fixed_comparison', Path(__file__).with_name('compare_paths.py'))
    comparison = importlib.util.module_from_spec(spec); spec.loader.exec_module(comparison)
    # Reuse the existing complete native comparison, including equal endpoints,
    # potential, target plane, source identities and explicit initial/final status.
    data = comparison.compare(repository, exported, initial=initial)
    prepared, channels = [], []
    for channel in ('preserving', 'T'):
        a, ar, raw_a, ap = source_peak(repository/OLD_RESULTS/channel)
        b, br, _, bp = source_peak(exported/channel, initial=initial)
        delta = b.get_scaled_positions(wrap=False)-a.get_scaled_positions(wrap=False)
        if np.max(np.abs(delta)) >= .25:
            raise ValueError('source fractions left the same parent lift; never infer new atom matching')
        ab, ba = crossover(a, b), crossover(b, a)
        boundary = clamped_plane_vcneb_boundary(12, a.cell.array, allow_tilt=True)
        boundary.validate_images([a, ab, ba, b])
        validate_path_geometry([a, ab, ba, b], cell_scale=5.12968067458423,
                               minimum_distance=1.6, maximum_deformation=.25)
        for name in ('INPUT', 'KPT'):
            if sha256(raw_a/name) != CONTRACT[name]:
                raise ValueError('source electronic input changed')
        for name, geometry in (('AB', ab), ('BA', ba)):
            prepared.append((channel, name, geometry, raw_a))
        channels.append(dict(channel=channel, sources={'A': ap, 'B': bp},
            corner_energy_eV_cell={'AA': ar['energy_eV_cell'], 'AB': None, 'BA': None, 'BB': br['energy_eV_cell']},
            maximum_source_fractional_difference=float(np.max(np.abs(delta))),
            additional_static_evaluations_required=2, measured_crossover_energies_available=False))
    output.mkdir(parents=True)
    decks = []
    for channel, corner, atoms, raw in prepared:
        folder = output/channel/corner; folder.mkdir(parents=True)
        for name in ('INPUT', 'KPT'):
            shutil.copyfile(raw/name, folder/name)
        write_structure(folder/'STRU', atoms)
        restored = read_fixed_hfo2_stru(folder/'STRU')
        if not same_ordered_geometry(restored, atoms):
            raise ValueError('case-native STRU roundtrip changed crossover geometry')
        decks.append(dict(channel=channel, corner=corner,
            path=folder.relative_to(output).as_posix(), input_sha256={n: sha256(folder/n) for n in ('INPUT', 'KPT', 'STRU')},
            geometry_roundtrip_with_case_reader=True, actual_HF_ASE_preflight_done=False))
    result = dict(status='four_inert_crossover_decks_no_DFT',
        upper_stage=data['upper_stage'], initializer_strains=[0., .01], target_strain=.005,
        formula_units=4, pressure_GPa=0., physical_input_sha256=CONTRACT,
        internal_coordinate_convention='fixed_ordered_unwrapped_fractions_in_common_parent_chart',
        affine_cell_change_not_fixed_Cartesian_atoms=True,
        channels=channels, decks=decks, additional_SCFs_required=4, new_DFT_calls=0,
        jobs_submitted=0, diagonals_recalculated=False, execution_admitted=False,
        original_endpoint_and_holdout_labels_unchanged=True,
        source_directory_on_HF=HF_PHYSICS, pseudo_orbital_bytes_not_distributed=True,
        unique_elastic_or_phonon_energy_partition=False,
        note='Initial-snapshot decks are diagnostic only. Do not execute automatically; '
             'review terminal E269 outcomes and regenerate from terminal sources if appropriate.')
    save(output/'preparation.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, default=Path(__file__).resolve().parents[4])
    parser.add_argument('--exported', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--initial', action='store_true', help='explicitly prepare diagnostic, unsubmitted initial-snapshot decks')
    args = parser.parse_args()
    result = prepare(args.repository, args.exported, args.output, initial=args.initial)
    print(json.dumps({k: v for k, v in result.items() if k not in ('channels', 'decks', 'physical_input_sha256')}))
