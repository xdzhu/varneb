"""Connect archived HfO2 channel energetics, ordered shuffles and strain.

Read only public training records. No calculator or scheduler is invoked.
The original native audits are reused rather than repeated as a new campaign.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from io import StringIO
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.analyze_hfo2_clamped_reference import original_chart_images
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.audit_hfo2_manuscript_work_cohort import CASES
from scripts.audit_hfo2_clamped_network import IDENTITIES
from vcneb.continuous_projection import continuous_reference_coordinates

BASE = Path('benchmarks/hfo2_channels/20261008')
NETWORK = BASE/'clamped_plus1_reversing_terminal_E148_20261014/analysis/network_final.json'
NETWORK_PIN = '6b8a206d2199b6727601955cef32b86513e22805c7be5f309e66a4258602145b'
ATOMIC = BASE/'G3_atomic_complete_E194_20261015/analysis.json'
ATOMIC_PIN = '3cc150ec3a6d0dbfa5f5879f7215a3692524309add3e5bb7f5abed384a197444'
CHANNELS = ('preserving', 'reversing', 'T', 'M')
LABELS = ('Preserving\nflip', 'Reversing\nflip', 'PO to T', 'PO to M')
COLORS = ('#306A94', '#D48345')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_bytes())


def path_coordinates(images, parent, basis):
    """Signed geometric coordinates in one fixed metric, not irrep labels."""
    if parent.get_chemical_symbols() != ['Hf']*4+['O']*8:
        raise ValueError('this case requires the original Hf4O8 ordered species')
    chart = continuous_reference_coordinates(images, parent)
    u = chart['displacements_A']
    centered = u-u.mean(axis=1, keepdims=True)
    q = np.einsum('aij,nij->na', basis, centered)
    polar = u[:, 4:].mean(axis=1)-u[:, :4].mean(axis=1)
    residual = centered-np.einsum('na,aij->nij', q, basis)
    return q, polar, chart, np.linalg.norm(residual.reshape(len(u), -1), axis=1)


def response_rows(rows):
    pairs = {(r['substrate_strain'], r['candidate']): r for r in rows}
    if set(pairs) != {(s, c) for s in (0., .01) for c in CHANNELS} or len(rows) != 8:
        raise ValueError('all eight training candidates required')
    answer = []
    for c in CHANNELS:
        lo, hi = pairs[0., c], pairs[.01, c]
        initial = (hi['matched_initial_energy_eV_cell']-lo['matched_initial_energy_eV_cell'])*250
        response = hi['forward_meV_fu']-lo['forward_meV_fu']
        endpoint = hi['final_minus_initial_meV_fu']-lo['final_minus_initial_meV_fu']
        answer.append(dict(channel=c, zero_peak_meV_fu=lo['forward_meV_fu'],
            plus1_peak_meV_fu=hi['forward_meV_fu'], delta_barrier_meV_fu=response,
            delta_initial_meV_fu=initial, delta_absolute_sampled_peak_meV_fu=response+initial,
            endpoint_only_delta_barrier_meV_fu=endpoint))
    if np.ptp([r['delta_initial_meV_fu'] for r in answer]) > 1e-8:
        raise ValueError('common initial well required')
    return answer


def build_data(repository):
    repository = Path(repository)
    if sha(repository/NETWORK) != NETWORK_PIN:
        raise ValueError('archived eight-band summary changed')
    network = load(repository/NETWORK)
    by_id = {r['source_job_id']: r for r in network['rows']}
    seeds_path = repository/BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json'
    seeds = load(seeds_path)
    t_path = repository/BASE/'reference_variants/T.vasp'
    t = read(t_path, format='vasp')
    parent, basis, _ = rotated_t_triplet(t)
    inputs = {NETWORK.as_posix(): NETWORK_PIN, seeds_path.relative_to(repository).as_posix(): sha(seeds_path),
              t_path.relative_to(repository).as_posix(): sha(t_path)}
    bands, geometry = [], []
    for _, case, pin, _, _ in CASES:
        receipt_path = repository/BASE/case/'audit_receipt.json'
        if sha(receipt_path) != pin:
            raise ValueError('archived terminal receipt differs')
        receipt = load(receipt_path)
        job = receipt['job_id']; c, strain = IDENTITIES[job]
        summary = by_id[job]
        folder = receipt_path.parent/'observations'/f"step_{receipt['terminal_step']:04d}"
        observation_path = folder/'observation.json'
        observation = load(observation_path)
        if observation != receipt['terminal_observation']:
            raise ValueError('geometry source no longer matches terminal receipt')
        images = []
        for i, record in enumerate(observation['raw_image_evaluations']):
            p = folder/f'POSCAR_{i:02d}'
            if sha(p) != record['snapshot_POSCAR_sha256']:
                raise ValueError('ordered source geometry changed')
            images.append(read(p, format='vasp'))
        if len(images) != 9:
            raise ValueError('seven interior plus two endpoint convention required')
        original = original_chart_images(images, seeds['orientation'], t)
        q, polar, chart, residual = path_coordinates(original, parent, basis)
        energy = np.array([r['energy_eV_cell'] for r in observation['raw_image_evaluations']])
        relative = (energy-energy[0])*250
        np.testing.assert_allclose(relative, summary['raw_relative_energy_meV_fu'], atol=1e-8, rtol=0)
        arc = np.asarray(observation['extended_reaction_coordinate_A'])
        if arc.shape != (9,) or arc[0] != 0 or np.any(np.diff(arc) <= 0):
            raise ValueError('monotonic native generalized arc required')
        records = []
        for i in range(9):
            record = dict(channel=c, substrate_strain=strain, image=i, s=float(arc[i]/arc[-1]),
                relative_energy_meV_fu=float(relative[i]), Q_x_A=float(q[i, 0]),
                Q_y_A=float(q[i, 1]), Q_z_A=float(q[i, 2]),
                polar_dx_A=float(polar[i, 0]), polar_dy_A=float(polar[i, 1]), polar_dz_A=float(polar[i, 2]),
                F_xx=float(chart['deformation_gradients'][i, 0, 0]),
                F_yy=float(chart['deformation_gradients'][i, 1, 1]),
                F_zz=float(chart['deformation_gradients'][i, 2, 2]),
                Green_yy=float(chart['green_strains'][i, 1, 1]), residual_atomic_A=float(residual[i]),
                source_POSCAR=(folder/f'POSCAR_{i:02d}').relative_to(repository).as_posix(),
                source_POSCAR_sha256=observation['raw_image_evaluations'][i]['snapshot_POSCAR_sha256'])
            records.append(record); geometry.append(record)
        bands.append(dict(channel=c, strain=strain, peak_image=int(np.argmax(relative)), rows=records))
        inputs[receipt_path.relative_to(repository).as_posix()] = pin
        inputs[observation_path.relative_to(repository).as_posix()] = sha(observation_path)
    if sha(repository/ATOMIC) != ATOMIC_PIN:
        raise ValueError('archived matched atomic-treatment result changed')
    atomic = load(repository/ATOMIC)
    treatments = []
    for c in ('preserving', 'T'):
        twins = [r for r in atomic['rows'] if r['channel'] == c]
        if len(twins) != 2 or any(len(r['points']) != 5 for r in twins):
            raise ValueError('all five observed local controls required for each center')
        for a, b in zip(twins[0]['points'], twins[1]['points']):
            if a['observed'] != b['observed'] or a['q_A'] != b['q_A']:
                raise ValueError('physical observations must not depend on curvature probe step')
            e = a['observed']
            treatments.append(dict(channel=c, id=a['id'], Q1_A=a['q_A'][0], Q2_A=a['q_A'][1],
                atomic_minus_frozen_meV_fu=(e['atomic']-e['frozen'])*250,
                joint_minus_atomic_meV_fu=(e['joint']-e['atomic'])*250,
                joint_minus_frozen_meV_fu=(e['joint']-e['frozen'])*250))
    inputs[ATOMIC.as_posix()] = ATOMIC_PIN
    return dict(responses=response_rows(network['rows']), bands=bands, geometry=geometry,
        treatments=treatments, basis_cartesian_A_unit_vectors=basis.tolist(),
        source_sha256=inputs, new_DFT_calls=0, full_native_audit_repeated=False,
        coordinate_kind='fixed-metric signed geometric patterns, not assigned phonon irreps',
        polar_proxy='oxygen minus hafnium mean displacement in fixed parent metric; not polarization')


def style():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 10.5, 'axes.labelsize': 11, 'xtick.labelsize': 9.5, 'ytick.labelsize': 9.5,
        'legend.fontsize': 9.5, 'axes.linewidth': .8, 'lines.linewidth': 1.6,
        'svg.fonttype': 'none', 'pdf.fonttype': 42, 'axes.spines.top': True,
        'axes.spines.right': True, 'axes.unicode_minus': True})
    return plt


def finish(fig, axes):
    axes = np.asarray(axes).ravel()
    for i, ax in enumerate(axes):
        ax.text(-.16, 1.055, f'({chr(97+i)})', transform=ax.transAxes,
                fontsize=13, fontweight='normal', va='bottom')
        ax.tick_params(direction='in', top=True, right=True, length=3.5)
        for spine in ax.spines.values(): spine.set_visible(True)
    for legend in fig.legends:
        legend.get_frame().set_facecolor('white'); legend.get_frame().set_alpha(.85)
        legend.get_frame().set_edgecolor('#888888')
    fig.canvas.draw()
    return fig


def competition_figure(data):
    plt = style()
    fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 96/25.4))
    fig.subplots_adjust(left=.105, right=.985, bottom=.23, top=.78, wspace=.40)
    x = np.arange(4); width = .34
    for j, (key, label) in enumerate((('zero_peak_meV_fu', '0% strain'), ('plus1_peak_meV_fu', '+1% strain'))):
        axes[0].bar(x+(j-.5)*width, [r[key] for r in data['responses']], width=width,
                    label=label, color=COLORS[j], edgecolor='#333333', linewidth=.5)
    axes[0].set_ylim(0, 380); axes[0].set_yticks([0, 100, 200, 300])
    axes[0].set_ylabel('Sampled peak (meV/HfO$_2$)')
    for j, (key, label) in enumerate((('delta_barrier_meV_fu', 'Observed peak change'),
                                     ('endpoint_only_delta_barrier_meV_fu', 'Final-well-only change'))):
        axes[1].bar(x+(j-.5)*width, [r[key] for r in data['responses']], width=width,
                    label=label, color=('#506D82', '#C9CED1')[j], edgecolor='#333333', linewidth=.5,
                    hatch=None if j == 0 else '//')
    axes[1].axhline(0, color='#555555', lw=.7)
    axes[1].set_ylim(-53, 40); axes[1].set_yticks([-40, -20, 0, 20, 40])
    axes[1].set_ylabel('Change at +1% (meV/HfO$_2$)')
    for ax in axes:
        ax.set_xticks(x, LABELS, fontsize=9.5)
    for ax in axes:
        pos = ax.get_position()
        fig.legend(*ax.get_legend_handles_labels(), loc='lower center',
                   bbox_to_anchor=(pos.x0+pos.width/2, .84), frameon=True, borderpad=.5)
    return finish(fig, axes)


def geometry_figure(data):
    plt = style()
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(3, 3, figsize=(183/25.4, 196/25.4))
    fig.subplots_adjust(left=.105, right=.99, bottom=.09, top=.89, wspace=.33, hspace=.48)
    for col, c in enumerate(('preserving', 'reversing', 'T')):
        for j, s in enumerate((0., .01)):
            band = next(b for b in data['bands'] if b['channel'] == c and b['strain'] == s)
            rows = band['rows']; arc = [r['s'] for r in rows]
            kw = dict(color=COLORS[j], marker=('o', 's')[j], markersize=3.5)
            axes[0, col].plot(arc, [r['relative_energy_meV_fu'] for r in rows], **kw)
            axes[1, col].plot(arc, [r['Q_x_A'] for r in rows], **kw)
            axes[1, col].plot(arc, [r['Q_z_A'] for r in rows], linestyle='--', **kw)
            axes[2, col].plot(arc, [100*(r['F_yy']-1) for r in rows], **kw)
        axes[0, col].set_ylim(-5, 370 if c == 'reversing' else 145)
        axes[0, col].set_yticks([0, 100, 200, 300] if c == 'reversing' else [0, 50, 100])
        axes[1, col].set_ylim(-1.17, 1.05); axes[1, col].set_yticks([-1, 0, 1])
        axes[1, col].axhline(0, color='#AAAAAA', lw=.5, zorder=0)
        axes[2, col].set_ylim(-1.6, 2.05); axes[2, col].set_yticks([-1, 0, 1, 2])
        axes[2, col].axhline(0, color='#AAAAAA', lw=.5, zorder=0)
        for ax in axes[:, col]: ax.set_xlim(-.035, 1.035); ax.set_xticks([0, .5, 1])
        axes[2, col].set_xlabel('Normalized path arc\n'+('Preserving flip', 'Reversing flip', 'PO to T')[col])
    axes[0, 0].set_ylabel('Energy (meV/HfO$_2$)')
    axes[1, 0].set_ylabel('Shuffle amplitude (Å)')
    axes[2, 0].set_ylabel('Open-axis stretch (%)')
    handles = [Line2D([], [], color=COLORS[j], marker=('o', 's')[j], label=label)
               for j, label in enumerate(('0% strain', '+1% strain'))]
    handles += [Line2D([], [], color='#555555', linestyle=ls, label=label)
                for ls, label in (('-', '$Q_x$ (shuffle row)'), ('--', '$Q_z$ (shuffle row)'))]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.54, .994),
               ncol=2, frameon=True, columnspacing=1.8)
    return finish(fig, axes)


def treatment_figure(data):
    plt = style()
    fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 94/25.4))
    fig.subplots_adjust(left=.10, right=.985, bottom=.26, top=.78, wspace=.38)
    for ax, c in zip(axes, ('preserving', 'T')):
        points = [p for p in data['treatments'] if p['channel'] == c]
        if [(p['Q1_A'], p['Q2_A']) for p in points] != [(0., 0.), (-.1, 0.), (0., -.1), (0., .1), (.1, 0.)]:
            raise ValueError('plotted labels require the original ordered five-coordinate cohort')
        x = np.arange(5)
        atomic = [-r['atomic_minus_frozen_meV_fu'] for r in points]
        extra = [-r['joint_minus_atomic_meV_fu'] for r in points]
        ax.bar(x, atomic, width=.64, color=COLORS[0], label='Atomic release', edgecolor='#333333', linewidth=.5)
        ax.bar(x, extra, bottom=atomic, width=.64, color=COLORS[1], label='Additional joint release', edgecolor='#333333', linewidth=.5)
        ax.set_xticks(x, ['(0, 0)', '(−.1, 0)', '(0, −.1)', '(0, .1)', '(.1, 0)'], rotation=40, ha='right')
        ax.set_xlabel('Fixed $(Q_1,Q_2)$ (Å): '+('preserving' if c == 'preserving' else 'T escape'))
        ax.set_ylabel('Energy lowering (meV/HfO$_2$)')
        ax.set_ylim(0, 13.5 if c == 'preserving' else .53)
    fig.legend(*axes[0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(.55, .98),
               ncol=2, frameon=True)
    return finish(fig, axes)


def export_figures(data, output):
    import matplotlib.pyplot as plt
    report = {}
    for name, builder in (('hfo2_competition', competition_figure), ('hfo2_path_geometry', geometry_figure),
                          ('hfo2_relaxation', treatment_figure)):
        fig = builder(data)
        for ext in ('pdf', 'svg', 'png'):
            if ext == 'svg':
                # Matplotlib emits insignificant trailing spaces in SVG paths.
                # Normalize generated text without altering coordinates or labels.
                buffer = StringIO()
                fig.savefig(buffer, format='svg', dpi=300, facecolor='white')
                (output/f'{name}.svg').write_text(
                    '\n'.join(line.rstrip() for line in buffer.getvalue().splitlines())+'\n', encoding='utf-8')
            else:
                fig.savefig(output/f'{name}.{ext}', dpi=300, facecolor='white')
        report[name] = dict(dimensions_mm=(fig.get_size_inches()*25.4).tolist(),
            axes_positions=[list(ax.get_position().bounds) for ax in fig.axes],
            normal_weight_panel_labels=True, all_spines_visible=all(s.get_visible() for ax in fig.axes for s in ax.spines.values()),
            individual_titles_empty=all(not ax.get_title() for ax in fig.axes),
            visual_review='not performed by exporter; see dated manual review in README')
        plt.close(fig)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repository', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    data = build_data(a.repository)
    a.output.mkdir(parents=True, exist_ok=True)
    for key in ('responses', 'geometry', 'treatments'):
        with (a.output/f'{key}.csv').open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(data[key][0]))
            writer.writeheader(); writer.writerows(data[key])
    with (a.output/'analysis.json').open('w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, allow_nan=False); f.write('\n')
    report = export_figures(data, a.output)
    with (a.output/'figure_checks.json').open('w', encoding='utf-8') as f:
        json.dump(report, f, indent=2); f.write('\n')
    print(json.dumps(dict(figures=list(report), archived_bands=len(data['bands']), new_DFT_calls=0)))


if __name__ == '__main__':
    main()
