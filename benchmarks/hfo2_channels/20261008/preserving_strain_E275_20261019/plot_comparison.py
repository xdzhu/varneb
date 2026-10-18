"""Plot one actually collected strain in the existing Python/matplotlib style."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np


def draw(source, strain, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    data = json.loads(Path(source).read_bytes())
    available = {r['strain'] for r in data['outcomes']}
    if strain not in available: raise ValueError('only actual completed strain outcomes can be plotted')
    bands = {b['kind']: b for b in data['bands'] if b['strain'] == strain}
    kinds = ('registered_preserving', 'released_seed', 'T_escape', 'upper_seed_midpoint')
    if not set(kinds[:3]) <= set(bands): raise ValueError('all three same-strain candidates required')
    style = {
        'registered_preserving': ('#306A94', 'o', '-', 'Registered preserving'),
        'released_seed': ('#398C76', 'D', '-', 'Released-seed preserving'),
        'T_escape': ('#D17A40', 's', '-', 'PO to T'),
        'upper_seed_midpoint': ('#82728D', '^', '--', 'Upper-seed preserving'),
    }
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 9.5, 'axes.labelsize': 10.5, 'xtick.labelsize': 9.5, 'ytick.labelsize': 9.5,
        'svg.fonttype': 'none', 'pdf.fonttype': 42, 'axes.linewidth': .8})
    fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 155/25.4))
    fig.subplots_adjust(left=.125, right=.985, bottom=.105, top=.84, wspace=.34, hspace=.43)
    for kind in kinds:
        if kind not in bands: continue
        color, marker, ls, _ = style[kind]; rows = bands[kind]['rows']
        for ax, field in zip(axes.flat, ('energy_meV_fu', 'polar_proxy_y_A',
                                        'open_y_stretch_percent', 'NEB_residual_eV_A')):
            selected = [r for r in rows if r[field] is not None]
            ax.plot([r['s'] for r in selected], [r[field] for r in selected], color=color,
                marker=marker, ms=3.8, lw=1.05, ls=ls, mfc='white' if kind == 'upper_seed_midpoint' else color,
                mew=.8)
    axes[0, 0].set_ylabel('Energy (meV/HfO$_2$)')
    axes[0, 1].set_ylabel('Transverse proxy ($\\AA$)')
    axes[1, 0].set_ylabel('Open-axis response (%)')
    axes[1, 1].set_ylabel('NEB residual (eV/$\\AA$)')
    axes[1, 1].axhline(.10, color='#777777', lw=.8, ls=':', zorder=0)
    axes[1, 1].set_ylim(bottom=0)
    for label, ax in zip('abcd', axes.flat):
        ax.set_xlim(-.025, 1.025); ax.set_xticks([0, .5, 1]); ax.set_xlabel('Normalized generalized arc')
        ax.tick_params(direction='in', top=True, right=True)
        for spine in ax.spines.values(): spine.set_visible(True)
        ax.text(-.12, 1.045, f'({label})', transform=ax.transAxes, fontsize=13, fontweight='normal')
    handles = [Line2D([], [], color=style[k][0], marker=style[k][1], ls=style[k][2],
               label=style[k][3], ms=3.8) for k in kinds if k in bands]
    legend = fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.54, .998), ncol=2,
        frameon=True, facecolor='white', edgecolor='#888888', framealpha=.85, fontsize=9.2)
    fig.align_labels(); fig.canvas.draw()
    # Check complete axes/text bounds in the selected backend before exporting.
    renderer = fig.canvas.get_renderer(); canvas = fig.bbox
    text_artists = [legend]+[text for ax in axes.flat for text in [ax.xaxis.label, ax.yaxis.label, *ax.texts]]
    if any(not canvas.contains(*artist.get_window_extent(renderer).get_points()[0])
           or not canvas.contains(*artist.get_window_extent(renderer).get_points()[1]) for artist in text_artists):
        plt.close(fig); raise ValueError('labels/legend extend outside the declared figure dimensions')
    output = Path(output); output.mkdir(parents=True, exist_ok=False)
    for extension in ('pdf', 'svg', 'png'):
        path = output/f'same_strain_comparison.{extension}'
        fig.savefig(path, dpi=300)
        if extension == 'svg':
            path.write_bytes(('\n'.join(line.rstrip() for line in path.read_text(encoding='utf-8').splitlines())+'\n').encode('utf-8'))
    plt.close(fig)
    rows = [r for kind in kinds if kind in bands for r in bands[kind]['rows']]
    with (output/'source_data.csv').open('x', encoding='utf-8', newline='') as h:
        writer = csv.DictWriter(h, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    qa = dict(strain=strain, source_analysis_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest(),
        backend='Python/matplotlib', size_mm=[183, 155], panel_labels='normal13pt(a)--(d)',
        all_axes_and_legend_bounds_inside=True, source_rows=len(rows), same_native_initial_well=True,
        smoothing=False, fitted_peak=False, statistical_replicates=False, inferred_energy_error_bars=False,
        structural_proxy_not_electronic_polarization=True, pending=data['pending_terminal_archives'],
        old_labels_and_scores_unchanged=data['old_labels_and_scores_unchanged'],
        physical_barrier_error_bound=None, stationary_TS_certified=False)
    with (output/'qa.json').open('x', encoding='utf-8', newline='\n') as h:
        json.dump(qa, h, indent=2, allow_nan=False); h.write('\n')
    return qa


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--source', type=Path, required=True)
    p.add_argument('--strain', type=float, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); print(json.dumps(draw(args.source, args.strain, args.output)))
