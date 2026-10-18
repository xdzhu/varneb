"""Connect existing F/A/J energies to the paper's fixed oxygen-pattern chart.

Thirty archived geometries, ten physical points, zero new DFT. The three
geometric T patterns are not the locally retained joint Q1/Q2 coordinates.
Displacement norms and their squared fractions are never energy fractions.
"""
from __future__ import annotations

import argparse
import csv
from io import StringIO
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.hfo2_G3_conditional_faces import load_inputs
from scripts.plot_hfo2_mechanism import (
    ATOMIC, ATOMIC_PIN, BASE, COLORS, finish, load, original_chart_images,
    path_coordinates, rotated_t_triplet, sha, style,
)

PINS = {
    BASE/'G3_curvature_E150_20261014/manifest.json':
        '5f047f54ca83cb2dceb8cb71687eb9206a85c04205cdf26431565bc8a9b98a5a',
    ATOMIC: ATOMIC_PIN,
    BASE/'G3_atomic_complete_E194_20261015/atomic_native_complete.json':
        '8a2fdf2a2beec306bb2f39aed8b347019a4f049b2d8cbc42cd79de1fce922b92',
    BASE/'G3_atomic_release_E175_20261015/config.json':
        '5bd401c6ca4044ac654f2b19a22371b20c17d3736d207ff7044336d9ed463f38',
    BASE/'G3_face_prefix_audit_E163_20261015/preserving_native_prefix.json':
        '3c2c90ab703c2e1cc7c9c699daee5e0ded8f0d54d1b5d038c33a862e1b89b8ee',
    BASE/'G3_face_prefix_audit_E163_20261015/T_native_prefix.json':
        '1a7b3530bb8d389196b9ad747df4f642e405ca1cdbb2d8da38d0ad529fefeda6',
}
POINTS = [[0., 0.], [-.1, 0.], [0., -.1], [0., .1], [.1, 0.]]


def displacement_parts(first, second, patterns):
    """Translation-free geometry in one reference metric, not mass weighting."""
    first, second, patterns = (np.asarray(x, dtype=float) for x in (first, second, patterns))
    if (first.shape != (12, 3) or second.shape != first.shape or patterns.shape != (3, 12, 3)
            or not all(np.isfinite(x).all() for x in (first, second, patterns))):
        raise ValueError('finite Hf4O8 displacements and three registered patterns required')
    flat = patterns.reshape(3, 36)
    if (not np.allclose(flat@flat.T, np.eye(3), atol=1e-12, rtol=0)
            or not np.allclose(patterns.sum(axis=1), 0, atol=1e-12, rtol=0)):
        raise ValueError('orthonormal translation-free reference patterns required')
    delta = second-first
    delta -= delta.mean(axis=0)
    q = np.einsum('aij,ij->a', patterns, delta)
    residual = delta-np.einsum('a,aij->ij', q, patterns)
    total2 = float(np.sum(delta**2))
    hf2, oxygen2 = float(np.sum(delta[:4]**2)), float(np.sum(delta[4:]**2))
    triplet2, residual2 = float(q@q), float(np.sum(residual**2))
    if not np.isclose(total2, triplet2+residual2, atol=1e-13, rtol=0):
        raise ValueError('geometric orthogonal decomposition differs')
    return dict(total_norm_A=total2**.5, Hf_norm_A=hf2**.5, oxygen_norm_A=oxygen2**.5,
        delta_Qx_A=float(q[0]), delta_Qy_A=float(q[1]), delta_Qz_A=float(q[2]),
        triplet_norm_A=triplet2**.5, complement_norm_A=residual2**.5,
        oxygen_squared_norm_fraction=None if total2 == 0 else oxygen2/total2,
        complement_squared_norm_fraction=None if total2 == 0 else residual2/total2,
        geometric_fractions_are_not_energy_fractions=True)


def analyze(repository):
    repository = Path(repository)
    pinned = {}
    for path, digest in PINS.items():
        if sha(repository/path) != digest:
            raise ValueError(f'archived treatment source changed: {path}')
        pinned[path] = load(repository/path)
    energy = pinned[ATOMIC]
    atomic = pinned[BASE/'G3_atomic_complete_E194_20261015/atomic_native_complete.json']
    config = pinned[BASE/'G3_atomic_release_E175_20261015/config.json']
    t_path = repository/BASE/'reference_variants/T.vasp'
    orientation_path = repository/BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json'
    t = read(t_path, format='vasp')
    parent, patterns, _ = rotated_t_triplet(t)
    orientation = load(orientation_path)['orientation']
    states, changes, centers = [], [], []
    for channel in ('preserving', 'T'):
        data = load_inputs(channel)  # Existing chart, reference and native-audit checks; no calculator.
        prefix_path = BASE/'G3_face_prefix_audit_E163_20261015'/f'{channel}_native_prefix.json'
        joint = {p['id']: p['candidate'] for p in pinned[prefix_path]['candidate_rows']}
        candidates = {r['id']: r['candidate'] for r in atomic['candidate_rows'] if r['channel'] == channel}
        twins = [r for r in energy['rows'] if r['channel'] == channel]
        # Model predictions at the two probe steps differ, but physical observations cannot.
        if (len(twins) != 2 or any(len(r['points']) != 5 for r in twins)
                or any(a['observed'] != b['observed'] or a['q_A'] != b['q_A']
                       for a, b in zip(twins[0]['points'], twins[1]['points']))):
            raise ValueError('two curvature steps cannot create new physical observations')
        points = config['channels'][channel]['points']
        if [p['id'] for p in points] != list(range(5)) or [p['q_A'] for p in points] != POINTS:
            raise ValueError('all five original matched coordinates required')
        centers.append(dict(channel=channel, source_job=data['manifest']['centers'][channel]['source_job_id'],
                            substrate_strain=0., sampled_peak_image=4))
        for point, observed in zip(points, twins[0]['points']):
            i, q = point['id'], point['q_A']
            a = candidates[i]
            j = data['origin_point'] if i == 0 else joint[i]
            full = [np.asarray(point['frozen_full_A']), np.asarray(a['full_coordinate']),
                    np.asarray(j['coordinates_full_A'])]
            if observed['id'] != i or observed['q_A'] != q:
                raise ValueError('physical energy and geometry identities differ')
            expected = [data['center_energy']+observed['observed'][key] for key in ('frozen', 'atomic', 'joint')]
            if not np.allclose([a['energy_eV_cell'], j['energy_eV_cell']], expected[1:], atol=2e-12, rtol=0):
                raise ValueError('candidate energies differ from the original treatment table')
            for x in full:
                if (x.shape != (40,) or not np.isfinite(x).all() or x[-1] != 0
                        or not np.allclose(data['axes'].T@(data['basis'].T@x), q, atol=2e-12, rtol=0)):
                    raise ValueError('retained local Q or prescribed strain changed')
            atoms = [data['chart'].to_atoms(x) for x in full]
            if not np.array_equal(atoms[0].cell.array, atoms[1].cell.array):
                raise ValueError('atomic-only treatment must retain the exact frozen-Q seed cell')
            original = original_chart_images(atoms, orientation, t)
            q_xyz, polar, chart, _ = path_coordinates(original, parent, patterns)
            for k, treatment in enumerate(('F', 'A', 'J')):
                states.append(dict(channel=channel, point_id=i, Q1_A=q[0], Q2_A=q[1], treatment=treatment,
                    energy_eV_cell=expected[k], Qx_A=float(q_xyz[k, 0]), Qy_A=float(q_xyz[k, 1]),
                    Qz_A=float(q_xyz[k, 2]), polar_proxy_z_A=float(polar[k, 2]),
                    open_y_stretch_percent=100*float(chart['deformation_gradients'][k, 1, 1]-1)))
            for name, first, second in (('A_minus_F', 0, 1), ('J_minus_A', 1, 2)):
                parts = displacement_parts(chart['displacements_A'][first], chart['displacements_A'][second], patterns)
                if not np.allclose([parts[f'delta_Q{d}_A'] for d in 'xyz'],
                                   q_xyz[second]-q_xyz[first], atol=1e-12, rtol=0):
                    raise ValueError('path and treatment projection conventions differ')
                changes.append(dict(channel=channel, point_id=i, Q1_A=q[0], Q2_A=q[1],
                    contrast=name, energy_change_meV_fu=250*(expected[second]-expected[first]), **parts))
    source_pins = {p.as_posix(): h for p, h in PINS.items()}
    source_pins.update({p.relative_to(repository).as_posix(): sha(p) for p in (t_path, orientation_path)})
    return dict(status='existing_zero_strain_treatment_geometry_not_new_conditional_surface',
        centers=centers, states=states, changes=changes, source_sha256=source_pins,
        pattern_basis=patterns.tolist(), physical_points=10, geometries=30,
        energy_unit='meV/HfO2 (four formula units per twelve-atom cell)',
        displacement_metric='original T Cartesian metric; common translation removed',
        new_DFT_calls=0, new_native_audit=False, independent_prediction_test=False,
        energy_partition_or_TS_or_stable_surface_certified=False)


def figure(data):
    plt = style()
    fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 175/25.4))
    fig.subplots_adjust(left=.105, right=.985, bottom=.13, top=.90, wspace=.40, hspace=1.08)
    labels = ['(0, 0)', '(-.1, 0)', '(0, -.1)', '(0, .1)', '(.1, 0)']
    for col, channel in enumerate(('preserving', 'T')):
        rows = {name: [r for r in data['changes'] if r['channel'] == channel and r['contrast'] == name]
                for name in ('A_minus_F', 'J_minus_A')}
        if any([p['point_id'] for p in r] != list(range(5)) for r in rows.values()):
            raise ValueError('all ten original physical points must remain visible')
        x = np.arange(5)
        atomic = [-r['energy_change_meV_fu'] for r in rows['A_minus_F']]
        extra = [-r['energy_change_meV_fu'] for r in rows['J_minus_A']]
        axes[0, col].bar(x, atomic, width=.64, color=COLORS[0], edgecolor='#333333', lw=.5,
                         label='Atomic release: F to A')
        axes[0, col].bar(x, extra, bottom=atomic, width=.64, color=COLORS[1], edgecolor='#333333', lw=.5,
                         label='Additional release: A to J')
        axes[0, col].set_ylabel('Energy lowering (meV/HfO$_2$)')
        axes[0, col].set_ylim(0, 13.5 if channel == 'preserving' else .53)
        for offset, key, label, color in ((-.18, 'triplet_norm_A', 'T-pattern triplet', '#568778'),
                                          (.18, 'complement_norm_A', 'Orthogonal complement', '#8B7B9B')):
            axes[1, col].bar(x+offset, [r[key] for r in rows['A_minus_F']], width=.34,
                label=label, color=color, edgecolor='#333333', lw=.5)
        axes[1, col].set_ylabel('F to A displacement norm (Å)')
        axes[1, col].set_ylim(0, .34 if channel == 'preserving' else .08)
        for ax in axes[:, col]:
            ax.set_xticks(x, labels, rotation=35, ha='right')
            ax.set_xlabel('Fixed $(Q_1,Q_2)$ (Å): '+('preserving' if channel == 'preserving' else 'T escape'))
    for row, top in ((0, .995), (1, .478)):
        fig.legend(*axes[row, 0].get_legend_handles_labels(), loc='upper center',
                   bbox_to_anchor=(.54, top), ncol=2, frameon=True, columnspacing=1.3)
    finish(fig, axes)
    # Align labels with the plot edge; long rotated y labels occupy the left gutter.
    for ax in axes.ravel():
        for text in ax.texts:
            if text.get_text() in ('(a)', '(b)', '(c)', '(d)'):
                text.set_x(0.)
    return fig


def export(data, output, *, draw=True):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    for name in ('states', 'changes'):
        with (output/f'{name}.csv').open('x', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[name][0]))
            writer.writeheader(); writer.writerows(data[name])
    with (output/'analysis.json').open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(data, handle, indent=2, allow_nan=False); handle.write('\n')
    if draw:
        fig = figure(data)
        try:
            for ext in ('png', 'pdf', 'svg'):
                if ext == 'svg':
                    stream = StringIO(); fig.savefig(stream, format='svg', facecolor='white')
                    (output/f'hfo2_accommodation.{ext}').write_text(
                        '\n'.join(s.rstrip() for s in stream.getvalue().splitlines())+'\n', encoding='utf-8')
                else:
                    fig.savefig(output/f'hfo2_accommodation.{ext}', dpi=300, facecolor='white')
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.repository)
    export(result, args.output)
    print(json.dumps(dict(physical_points=10, geometries=30, new_DFT_calls=0,
        origins=[r for r in result['changes'] if r['point_id'] == 0 and r['contrast'] == 'A_minus_F'])))
