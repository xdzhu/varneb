"""Register the six physical bottleneck neighborhoods from existing results.

No new forces, derivatives, fit or electronic calculation. This table keeps
ordinary convergence separate from allowed-space physical stationarity.
"""
import argparse
import json
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

from scripts.analyze_hfo2_physical_holdout import analyze, BASE, RESULTS, OBSERVATION_PINS
from scripts.audit_hfo2_manuscript_work_cohort import CASES
from scripts.audit_hfo2_clamped_network import IDENTITIES
from scripts.audit_hfo2_static_replica import sha256
import numpy as np
from ase.io import read
from vcneb import VCNEB, clamped_plane_vcneb_boundary


def peak_physical_force(source, observation, index):
    """Replay old and new snapshots in the same explicit allowed-force API.

    Earlier archived schemas lack the explicit allowed-force fields. Never
    replace them with a raw force, which would include substrate reactions.
    """
    trajectory = source.parent/'evaluated_chain.traj'
    if sha256(trajectory) != observation['evaluated_chain_sha256']:
        raise ValueError('original evaluated trajectory changed')
    images = read(trajectory, index=':')
    for image, raw in zip(images, observation['raw_image_evaluations']):
        for name, expected in (('energy', raw['energy_eV_cell']), ('forces', raw['forces_eV_A']),
                               ('stress', raw['stress_ASE_voigt_eV_A3'])):
            if not np.allclose(image.calc.results[name], expected, atol=1e-12, rtol=0):
                raise ValueError('cached E/F/stress changed')
    boundary = clamped_plane_vcneb_boundary(12,
        np.asarray(observation['mechanical_boundary']['reference_cell_A']), allow_tilt=True)
    chain = VCNEB(images, k=.2, cell_scale=5.12968067458423, climb=False, pressure=0,
                  **boundary.vcneb_kwargs(images))
    saddle = chain.saddle_diagnostics()
    diagnostics = chain.path_diagnostics()
    if diagnostics['highest_image_index'] != index:
        raise ValueError('replayed peak differs from registered sampled maximum')
    return saddle, diagnostics['images'][index]


def register(repository):
    data = analyze(repository)
    sources = {}
    for _, case, pin, _, _ in CASES:
        receipt_path = repository/BASE/case/'audit_receipt.json'
        if sha256(receipt_path) != pin:
            raise ValueError('original completed source changed')
        receipt = json.loads(receipt_path.read_bytes())
        channel, strain = IDENTITIES[receipt['job_id']]
        if channel in ('preserving', 'T'):
            sources[channel, strain] = receipt_path.parent/'observations'/f"step_{receipt['terminal_step']:04d}"/'observation.json'
    for channel in ('preserving', 'T'):
        folder = repository/RESULTS/channel
        audit = json.loads((folder/'audit.json').read_bytes())
        path = folder/f"step_{audit['step']:04d}"/'observation.json'
        if sha256(path) != OBSERVATION_PINS[channel]:
            raise ValueError('original midpoint observation changed')
        sources[channel, .005] = path
    rows = []
    for channel in ('preserving', 'T'):
        for strain in (0., .005, .01):
            source = sources[channel, strain]
            observation = json.loads(source.read_bytes())
            band = next(b for b in data['bands'] if b['channel'] == channel and b['strain'] == strain)
            index = band['peak_image']; geometry = band['rows'][index]
            saddle, physical = peak_physical_force(source, observation, index)
            if physical['image_index'] != index:
                raise ValueError('physical-gradient row differs from sampled peak')
            rows.append(dict(channel=channel, substrate_strain=strain,
                source_job=observation['source_job_id'], sampled_peak_image=index,
                sampled_peak_meV_fu=geometry['relative_energy_meV_fu'],
                ordinary_chain_fmax_eV_A=observation['replayed_fmax_eV_A'],
                peak_allowed_physical_max_vector_eV_A=physical['allowed_true_generalized_force_max_vector_eV_per_A'],
                peak_allowed_physical_norm_eV_A=physical['allowed_true_generalized_force_norm_eV_per_A'],
                peak_allowed_tangential_force_eV_A=saddle['allowed_true_tangential_force_eV_per_A'],
                peak_raw_atom_force_eV_A=physical['max_atom_force_eV_per_A'],
                Qx_A=geometry['Q_x_A'], Qz_A=geometry['Q_z_A'],
                open_y_stretch_percent=100*(geometry['F_yy']-1),
                observation=source.relative_to(repository).as_posix(), observation_sha256=sha256(source),
                stationary_bottleneck_certified=False, Hessian_index_certified=False))
    return dict(status='six_existing_path_neighborhoods_for_mechanism_followup',
        rows=rows, new_DFT_calls=0, new_curvature_probes=0,
        ordinary_force_target_eV_A=.10, energy_unit='meV/HfO2',
        constrained_reaction_stress_excluded_from_allowed_force=True,
        force_residual_is_not_a_barrier_error_bound=True,
        original_barriers_and_holdout_scores_unchanged=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, default=Path(__file__).resolve().parents[4])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = register(args.repository)
    with args.output.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(result, allow_nan=False))
