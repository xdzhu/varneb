"""Completed same-pattern T controls: local energies and actual atomic motion.

All six native candidates are replayed. No interpolation, statistical
replicates, corrected barrier, or cross-definition P/T energy subtraction.
"""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from ase.io import read
from scripts.analyze_hfo2_local_accommodation import displacement_parts
from scripts.audit_hfo2_static_replica import sha256
from scripts.plot_hfo2_mechanism import COLORS, original_chart_images, path_coordinates, rotated_t_triplet

BASE = ROOT/'benchmarks/hfo2_channels/20261008'
CASE = BASE/'matched_pattern_seeded_E271_20261018'


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec); spec.loader.exec_module(result)
    return result


def build_data():
    replay = module(CASE/'analyze_complete.py', 'matched_T_native_analysis').analyze()
    original = BASE/'matched_pattern_release_E270_20261018'
    study = module(original/'study.py', 'matched_T_coordinates')
    config = json.loads((original/'prepared/config.json').read_bytes())
    t_path = BASE/'reference_variants/T.vasp'
    orientation_path = BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json'
    if sha256(t_path) != config['T_reference_sha256'] or sha256(orientation_path) != config['orientation_sha256']:
        raise ValueError('original physical coordinate frame changed')
    t = read(t_path); parent, patterns, _ = rotated_t_triplet(t)
    orientation = json.loads(orientation_path.read_bytes())['orientation']
    rows = []
    with tarfile.open(CASE/'native_terminal.tar.gz') as archive:
        for index, energy in enumerate(replay['rows']):
            chart = study.chart_for(original/'prepared/seeds'/str(index))
            candidates = [json.load(archive.extractfile(f'native_terminal/{index}/treatment_{arm}.json'))['candidate']
                          for arm in ('A', 'J')]
            atoms = [chart.reference_atoms]+[chart.to_atoms(c['coordinates']) for c in candidates]
            if not np.array_equal(atoms[0].cell.array, atoms[1].cell.array):
                raise ValueError('atomic-only cell changed')
            images = original_chart_images(atoms, orientation, t)
            q, polar, physical, _ = path_coordinates(images, parent, patterns)
            defect = max(float(np.max(np.abs(q[1:]-q[0]))), float(np.max(np.abs(polar[1:,2]-polar[0,2]))))
            if defect > 1e-12:
                raise ValueError('same-pattern candidate changed the actual displayed physical variables')
            u = physical['displacements_A']
            a = displacement_parts(u[0], u[1], patterns)
            j = displacement_parts(u[0], u[2], patterns)
            aj = displacement_parts(u[1], u[2], patterns)
            # Same original-T open-axis convention used by the path figure.
            f = [np.linalg.solve(t.cell.array, image.cell.array).T for image in images]
            rows.append(dict(energy, atomic_displacement_F_A_norm_A=a['total_norm_A'],
                atomic_displacement_F_J_norm_A=j['total_norm_A'],
                atomic_displacement_A_J_norm_A=aj['total_norm_A'],
                joint_minus_frozen_open_y_stretch_percent=100*(f[2][1,1]-f[0][1,1]),
                measured_four_pattern_defect_A=defect))
    return dict(rows=rows, archives=replay['archives'], costs=replay['costs'],
        native_candidates_independently_replayed=6, new_DFT_calls=0,
        held_definitions_match_but_native_values_differ_across_strain=True,
        older_two_axis_controls_not_nested_with_this_four_pattern_control=True,
        geometric_norms_are_not_RMS_or_energy_fractions=True,
        local_energy_lowering_is_not_barrier_error_bound_or_stability_certificate=True)


def draw(data, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    plt.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':12, 'font.weight':'normal', 'axes.labelsize':12, 'axes.labelweight':'normal',
        'axes.linewidth':.9, 'svg.fonttype':'none', 'pdf.fonttype':42,
        'xtick.labelsize':11.5, 'ytick.labelsize':11.5, 'legend.fontsize':11})
    fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 94/25.4))
    fig.subplots_adjust(left=.10, right=.985, bottom=.21, top=.76, wspace=.34)
    rows = data['rows']; x = np.arange(3)
    atomic = -np.array([r['dE_A_minus_F_meV_fu'] for r in rows])
    joint_increment = -np.array([r['dE_J_minus_A_meV_fu'] for r in rows])
    axes[0].bar(x, atomic, width=.48, color=COLORS[0], edgecolor='#25323D', linewidth=.6, label=r'$F-A$')
    axes[0].bar(x, joint_increment, bottom=atomic, width=.48, color=COLORS[1],
                edgecolor='#25323D', linewidth=.6, label=r'$A-J$')
    axes[0].set_ylabel('Energy lowering\n(meV/HfO$_2$)')
    axes[0].set_ylim(0, .30); axes[0].yaxis.set_major_locator(MultipleLocator(.1))
    axes[1].bar(x-.14, [r['atomic_displacement_F_A_norm_A'] for r in rows], width=.26,
        color=COLORS[0], edgecolor='#25323D', linewidth=.6, label=r'$F\to A$')
    axes[1].bar(x+.14, [r['atomic_displacement_F_J_norm_A'] for r in rows], width=.26,
        color=COLORS[1], edgecolor='#25323D', linewidth=.6, label=r'$F\to J$')
    axes[1].set_ylabel('Atomic displacement\nnorm (Å)')
    axes[1].set_ylim(0, .12); axes[1].yaxis.set_major_locator(MultipleLocator(.04))
    labels = []
    for ax, letter in zip(axes, ('a', 'b')):
        ax.set_xticks(x, ['0', '0.5', '1'])
        ax.set_xlim(-.65, 2.65); ax.set_xlabel('Substrate strain (%)', labelpad=6)
        ax.tick_params(which='both', top=True, right=True, direction='in', length=4)
        for spine in ax.spines.values(): spine.set_visible(True)
        legend = ax.legend(loc='lower center', bbox_to_anchor=(.54, 1.10), ncol=2,
            frameon=True, facecolor='white', edgecolor='#A8ADB2', framealpha=.85,
            borderpad=.35, handlelength=1.3, handletextpad=.45, columnspacing=.9)
        labels.append(ax.text(0, 1.055, f'({letter})', transform=ax.transAxes,
                              fontsize=13, fontweight='normal', va='bottom'))
    fig.canvas.draw()
    height = fig.canvas.get_width_height()[1]
    bounds = [ax.get_position().bounds for ax in axes]
    if not np.allclose(np.array(bounds)[:,[1,3]], np.array(bounds)[0,[1,3]], atol=1e-12, rtol=0):
        raise ValueError('panel axes are not vertically aligned')
    for text in labels:
        box = text.get_window_extent(fig.canvas.get_renderer())
        if box.y1 > height: raise ValueError('panel letter clipped')
    for extension in ('pdf', 'svg', 'png'):
        fig.savefig(output/f'matched_T_release.{extension}', dpi=300, facecolor='white')
    plt.close(fig)
    return dict(size_mm=[183,94], axes_bounds=bounds, panel_letter_weight='normal',
        panel_letter_size_pt=13, top_and_right_spines=True, no_subpanel_titles=True,
        framed_white_legend_alpha=.85, smoothing_or_error_bars_or_synthetic_data=False,
        backend='Python/matplotlib', native_candidates_independently_replayed=6)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); data = build_data(); args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/'source_data.csv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=list(data['rows'][0])); writer.writeheader(); writer.writerows(data['rows'])
    qa = draw(data,args.output)
    qa.update(plot_source_sha256=sha256(Path(__file__)), source_data_sha256=sha256(args.output/'source_data.csv'))
    for name, value in [('analysis.json',data),('figure_qa.json',qa)]:
        (args.output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(data))
