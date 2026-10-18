"""Replay the frozen +0.5% physical test and connect it to the paper mechanism.

Uses completed native exports and the original forecast panel. No DFT, fitting,
model selection, interpolation of unseen labels or scheduler mutation occurs.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
from io import StringIO
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.plot_hfo2_mechanism import (
    BASE, COLORS, build_data, finish, geometry_figure, load, original_chart_images,
    path_coordinates, rotated_t_triplet, sha, style,
)
from vcneb.forecast_freeze import load_frozen_forecast_panel

CASE = BASE/'physical_holdout_E268_20261018'
RESULTS = CASE/'completed_paths_20261018'
FORECAST_PIN = '1dd8dc4cacb22b073a6402acb38ffb956b3e0dd3738a85db18fab2afb279a88b'
OBSERVATION_PINS = {
    'preserving': '31e4c8458e9e1facce7f7531ceb901a4b940f9a4b715219f54e5e810871c1810',
    'T': 'a3e654b8a2cfe5ffaeb4dc39249c6c04df20820fe45ce9586bbfc4792cacfe41',
}


def analyze(repository):
    repository = Path(repository)
    data = build_data(repository)
    forecasts = load_frozen_forecast_panel(repository/CASE/'forecast_panel/frozen',
        expected_manifest_sha256=FORECAST_PIN)
    spec = importlib.util.spec_from_file_location('E268_score', repository/CASE/'audit_paths.py')
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    t = read(repository/BASE/'reference_variants/T.vasp')
    parent, basis, _ = rotated_t_triplet(t)
    orientation = load(repository/BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json')['orientation']
    references, outcomes, initial_energies = {}, [], []
    for short, channel in zip(('preserving', 'T'), reader.CHANNELS):
        folder = repository/RESULTS/short
        audit = load(folder/'audit.json')
        observation_path = folder/f"step_{audit['step']:04d}"/'observation.json'
        if (sha(observation_path) != OBSERVATION_PINS[short]
                or audit['observation_sha256'] != OBSERVATION_PINS[short]
                or audit['forecast_manifest_sha256'] != FORECAST_PIN
                or not audit['ordinary_converged'] or not audit['native_all_calls_reparsed']
                or audit['channel'] != channel):
            raise ValueError('complete original two-path result required')
        observation = load(observation_path)
        if observation['formula_units'] != 4 or observation['climb'] or observation['pressure_GPa'] != 0:
            raise ValueError('registered physical target differs')
        images = []
        for i, r in enumerate(observation['raw_image_evaluations']):
            path = observation_path.parent/f'POSCAR_{i:02d}'
            if sha(path) != r['snapshot_POSCAR_sha256']:
                raise ValueError('ordered path geometry changed')
            images.append(read(path, format='vasp'))
        original = original_chart_images(images, orientation, t)
        q, polar, chart, residual = path_coordinates(original, parent, basis)
        energies = np.array([r['energy_eV_cell'] for r in observation['raw_image_evaluations']])
        initial_energies.append(float(energies[0]))
        relative = (energies-energies[0])*250
        if len(images) != 9 or not np.isclose(relative.max(), audit['sampled_peak_eV_fu']*1000, atol=1e-9):
            raise ValueError('nine-image sampled peak does not replay')
        arc = np.array(observation['extended_reaction_coordinate_A'])
        records = []
        for i in range(9):
            records.append(dict(channel=short, substrate_strain=.005, image=i, s=float(arc[i]/arc[-1]),
                relative_energy_meV_fu=float(relative[i]), Q_x_A=float(q[i, 0]), Q_y_A=float(q[i, 1]),
                Q_z_A=float(q[i, 2]), polar_dx_A=float(polar[i, 0]), polar_dy_A=float(polar[i, 1]),
                polar_dz_A=float(polar[i, 2]), F_xx=float(chart['deformation_gradients'][i, 0, 0]),
                F_yy=float(chart['deformation_gradients'][i, 1, 1]), F_zz=float(chart['deformation_gradients'][i, 2, 2]),
                Green_yy=float(chart['green_strains'][i, 1, 1]), residual_atomic_A=float(residual[i]),
                source_POSCAR=(observation_path.parent/f'POSCAR_{i:02d}').relative_to(repository).as_posix(),
                source_POSCAR_sha256=observation['raw_image_evaluations'][i]['snapshot_POSCAR_sha256']))
        band = dict(channel=short, strain=.005, peak_image=int(relative.argmax()), rows=records)
        data['bands'].append(band); data['geometry'].extend(records)
        anchor = next(b for b in data['bands'] if b['channel'] == short and b['strain'] == 0)
        peak = records[band['peak_image']]
        outcomes.append(dict(channel=short, sampled_peak_meV_fu=float(relative.max()),
            response_from_zero_meV_fu=float(relative.max()-max(r['relative_energy_meV_fu'] for r in anchor['rows'])),
            peak_Qx_A=peak['Q_x_A'], peak_Qz_A=peak['Q_z_A'], peak_open_stretch_percent=100*(peak['F_yy']-1),
            force_eV_A=audit['final_fmax_eV_A'], updates=audit['step'], native_SCFs=audit['new_SCF_calls'],
            SCF_core_hours=audit['SCF_core_hours'], allocation_core_hours=audit['allocation_core_hours']))
        references[channel] = audit
        data['source_sha256'][observation_path.relative_to(repository).as_posix()] = sha(observation_path)
        data['source_sha256'][(folder/'audit.json').relative_to(repository).as_posix()] = sha(folder/'audit.json')
    if np.ptp(initial_energies) > 1e-10:
        raise ValueError('both channels must start from the same cached PO+ well')
    scores = reader.point_scores(forecasts, references)
    gap = outcomes[1]['sampled_peak_meV_fu']-outcomes[0]['sampled_peak_meV_fu']
    gaps = []
    for strain in (0., .005, .01):
        values = {b['channel']: max(r['relative_energy_meV_fu'] for r in b['rows'])
                  for b in data['bands'] if b['strain'] == strain and b['channel'] in ('preserving', 'T')}
        gaps.append(dict(strain=strain, T_minus_preserving_meV_fu=values['T']-values['preserving']))
    endpoint_cost = load(repository/CASE/'forecast_panel/cost.json')
    cost = dict(endpoint_native_SCFs=endpoint_cost['visible_endpoint_SCFs'],
                path_native_SCFs=sum(r['native_SCFs'] for r in outcomes),
                endpoint_and_path_SCF_core_hours=endpoint_cost['endpoint_SCF_core_hours']+sum(r['SCF_core_hours'] for r in outcomes),
                endpoint_and_path_allocation_core_hours=endpoint_cost['endpoint_allocation_core_hours']+sum(r['allocation_core_hours'] for r in outcomes),
                all_earlier_training_and_model_preparation_additional=True)
    data.update(holdout_outcomes=outcomes, holdout_scores=scores, selected_pair_gaps=gaps, cost=cost,
        holdout_gap_meV_fu=gap, holdout_gap_change_from_zero_meV_fu=gap-gaps[0]['T_minus_preserving_meV_fu'],
        forecast_manifest_sha256=FORECAST_PIN, numeric_forecasts=sum(r.prediction_eV_fu is not None for r in forecasts),
        unavailable_forecasts=sum(r.prediction_eV_fu is None for r in forecasts),
        no_new_model_advantage=True, full_four_channel_network_at_half_strain=False)
    return data


def independent_figure(data):
    plt = style()
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 165/25.4))
    fig.subplots_adjust(left=.105, right=.985, bottom=.10, top=.90, wspace=.36, hspace=.62)
    conditions = ((0., COLORS[0], 'o'), (.005, '#3C8C78', 'D'), (.01, COLORS[1], 's'))
    for ax, channel in zip(axes[0], ('preserving', 'T')):
        for strain, color, marker in conditions:
            band = next(b for b in data['bands'] if b['channel'] == channel and b['strain'] == strain)
            ax.plot([r['s'] for r in band['rows']], [r['relative_energy_meV_fu'] for r in band['rows']],
                    color=color, marker=marker, markersize=3.4)
        ax.set_xlim(-.03, 1.03); ax.set_ylim(-5, 140); ax.set_yticks([0, 50, 100])
        ax.set_xlabel('Normalized path arc\n'+('Preserving flip' if channel == 'preserving' else 'PO to T'))
        ax.set_ylabel('Energy (meV/HfO$_2$)')
    fig.legend(handles=[Line2D([], [], color=c, marker=m, label=f'{s*100:g}% strain')
                       for s, c, m in conditions], loc='upper center', bbox_to_anchor=(.54, .995), ncol=3)
    scores = {r['model']: r for r in data['holdout_scores']}
    for j, (channel, case, label) in enumerate((('preserving', 'PO_flip_T_pattern_preserving', 'Preserving flip'),
                                             ('T', 'PO_to_T', 'PO to T'))):
        peaks = [max(r['relative_energy_meV_fu'] for r in next(b for b in data['bands']
                    if b['channel'] == channel and b['strain'] == s)['rows']) for s in (0., .005, .01)]
        axes[1, 0].plot([0, .5, 1], peaks, color=COLORS[j], marker=('o', 's')[j], label=label)
        forecast = scores['B0']['points'][case]['prediction_eV_fu']*1000
        axes[1, 0].plot(.5, forecast, color=COLORS[j], marker='D', mfc='white', ms=6, zorder=5)
        axes[1, 1].bar(np.arange(3)+(j-.5)*.34,
            [scores[m]['points'][case]['signed_error_meV_fu'] for m in ('B0', 'B1_fixed', 'B1_follow')],
            width=.34, color=COLORS[j], edgecolor='#333333', linewidth=.5)
    axes[1, 0].set_xlabel('Substrate strain (%)'); axes[1, 0].set_xticks([0, .5, 1])
    axes[1, 0].set_ylabel('Sampled peak (meV/HfO$_2$)'); axes[1, 0].set_ylim(78, 133)
    axes[1, 1].axhline(0, color='#555555', lw=.7)
    axes[1, 1].set_ylabel('Prediction − DFT (meV/HfO$_2$)'); axes[1, 1].set_ylim(-10, 15)
    axes[1, 1].set_xticks([0, 1, 2], ['B0\nlinear', 'B1\nfixed peak', 'B1\nfinal-follow'])
    handles, labels = axes[1, 0].get_legend_handles_labels()
    handles.append(Line2D([], [], color='#555555', marker='D', mfc='white', linestyle='none'))
    labels.append('Frozen B0 at +0.5%')
    fig.legend(handles, labels, loc='center', bbox_to_anchor=(.55, .485), ncol=3,
               columnspacing=1.1, handlelength=1.5, fontsize=9)
    return finish(fig, axes)


def export(data, output):
    import matplotlib.pyplot as plt
    output.mkdir(parents=True, exist_ok=True)
    for name, fig in [('hfo2_independent_response', independent_figure(data)),
                      ('hfo2_path_geometry', geometry_figure(data, conditions=(0., .005, .01)))]:
        fig.savefig(output/f'{name}.pdf', facecolor='white')
        fig.savefig(output/f'{name}.png', dpi=300, facecolor='white')
        buffer = StringIO(); fig.savefig(buffer, format='svg', facecolor='white')
        (output/f'{name}.svg').write_text('\n'.join(s.rstrip() for s in buffer.getvalue().splitlines())+'\n', encoding='utf-8')
        plt.close(fig)
    errors = [dict(model=r['model'], channel=case, **point) for r in data['holdout_scores']
              for case, point in r['points'].items()]
    for name, rows in [('outcomes', data['holdout_outcomes']), ('geometry', data['geometry']),
                       ('pair_gaps', data['selected_pair_gaps']), ('forecast_errors', errors)]:
        with (output/f'{name}.csv').open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, list(rows[0])); writer.writeheader(); writer.writerows(rows)
    with (output/'analysis.json').open('w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, allow_nan=False); f.write('\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repository', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); data = analyze(args.repository); export(data, args.output)
    print(json.dumps({k: data[k] for k in ('holdout_outcomes', 'selected_pair_gaps', 'numeric_forecasts', 'unavailable_forecasts')}))
