"""Three-strain material summary; refuses a missing original target outcome."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

KINDS = ('registered_preserving', 'released_seed', 'T_escape')
CONDITIONS = (0., .005, .01)


def summary_rows(data):
    if (not data['complete_selected_target_pair_available'] or data['pending_terminal_archives']
            or {r['strain'] for r in data['outcomes']} != set(CONDITIONS)):
        raise ValueError('complete original three-strain target comparison required; never fill a missing result')
    rows = []
    for strain in CONDITIONS:
        group = {b['kind']: b for b in data['bands'] if b['strain'] == strain and b['kind'] in KINDS}
        if set(group) != set(KINDS): raise ValueError('all three actual candidates required')
        for kind in KINDS:
            b = group[kind]; peak = b['rows'][b['sampled_peak_image']]
            rows.append(dict(strain=strain, kind=kind, source_job=b['source_job'],
                peak_meV_fu=b['sampled_polygon_peak_meV_fu'],
                T_minus_preserving_meV_fu=None if kind == 'T_escape' else
                    group['T_escape']['sampled_polygon_peak_meV_fu']-b['sampled_polygon_peak_meV_fu'],
                peak_transverse_proxy_A=peak['polar_proxy_y_A'],
                peak_open_y_percent=peak['open_y_stretch_percent'],
                ordinary_residual_pass=b['ordinary_residual_pass'], fmax_eV_A=b['fmax_eV_A']))
    values = [[r[k] for k in ('peak_meV_fu', 'peak_transverse_proxy_A', 'peak_open_y_percent', 'fmax_eV_A')]
              for r in rows]
    if not np.isfinite(values).all(): raise ValueError('finite actual energies/geometry/residuals required')
    return rows


def draw(source, output):
    data = json.loads(Path(source).read_bytes()); rows = summary_rows(data)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    colors = ('#306A94', '#398C76', '#D17A40'); markers = ('o', 'D', 's')
    labels = ('Registered preserving', 'Released-source preserving', 'PO to T')
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 9.5, 'axes.labelsize': 10.5, 'xtick.labelsize': 9.5, 'ytick.labelsize': 9.5,
        'svg.fonttype': 'none', 'pdf.fonttype': 42, 'axes.linewidth': .8})
    fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 155/25.4))
    fig.subplots_adjust(left=.13, right=.98, bottom=.105, top=.85, wspace=.37, hspace=.43)
    x = np.arange(3)
    for j, kind in enumerate(KINDS):
        selected = [next(r for r in rows if r['kind'] == kind and r['strain'] == s) for s in CONDITIONS]
        bars = axes[0, 0].bar(x+(j-1)*.23, [r['peak_meV_fu'] for r in selected], width=.21,
            color=colors[j], edgecolor='#333333', linewidth=.55)
        for bar, r in zip(bars, selected):
            if not r['ordinary_residual_pass']: bar.set_hatch('//'); bar.set_alpha(.45)
        if kind != 'T_escape':
            bars = axes[0, 1].bar(x+(-.14 if j == 0 else .14),
                [r['T_minus_preserving_meV_fu'] for r in selected], width=.25,
                color=colors[j], edgecolor='#333333', linewidth=.55)
            for bar, r in zip(bars, selected):
                if not r['ordinary_residual_pass']: bar.set_hatch('//'); bar.set_alpha(.45)
        for ax, field in zip(axes[1], ('peak_transverse_proxy_A', 'peak_open_y_percent')):
            for i, r in enumerate(selected):
                ax.plot(i+(j-1)*.10, r[field], marker=markers[j], ls='none', color=colors[j],
                    ms=5.5, mfc=colors[j] if r['ordinary_residual_pass'] else 'white', mew=1.)
    axes[0, 0].set_ylabel('Sampled peak (meV/HfO$_2$)'); axes[0, 0].set_ylim(bottom=0)
    axes[0, 1].set_ylabel('T $-$ preserving\n(meV/HfO$_2$)')
    axes[0, 1].axhline(0, color='#777777', lw=.8, zorder=0)
    axes[1, 0].set_ylabel('Transverse proxy\nat peak ($\\AA$)')
    axes[1, 1].set_ylabel('Open-axis response\nat peak (%)')
    for letter, ax in zip('abcd', axes.flat):
        ax.set_xlim(-.5, 2.5); ax.set_xticks(x, ['0', '0.5', '1'])
        ax.set_xlabel('Substrate strain (%)'); ax.tick_params(direction='in', top=True, right=True)
        for spine in ax.spines.values(): spine.set_visible(True)
        ax.text(-.13, 1.045, f'({letter})', transform=ax.transAxes, fontsize=13, fontweight='normal')
    handles = [Line2D([], [], color=c, marker=m, ls='none', label=label, ms=5.)
               for c, m, label in zip(colors, markers, labels)]
    legend = fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.54, .995),
        ncol=2, frameon=True, facecolor='white', edgecolor='#888888', framealpha=.85, fontsize=9.2)
    fig.align_labels(); fig.canvas.draw(); renderer = fig.canvas.get_renderer(); canvas = fig.bbox
    artists = [legend]+[t for ax in axes.flat for t in [ax.xaxis.label, ax.yaxis.label, *ax.texts]]
    if any(not canvas.contains(*a.get_window_extent(renderer).get_points()[0]) or
           not canvas.contains(*a.get_window_extent(renderer).get_points()[1]) for a in artists):
        plt.close(fig); raise ValueError('axes/legend labels extend outside the declared size')
    output = Path(output); output.mkdir(parents=True, exist_ok=False)
    for ext in ('pdf', 'svg', 'png'):
        path = output/f'strain_mechanism_summary.{ext}'; fig.savefig(path, dpi=300)
        if ext == 'svg':
            path.write_bytes(('\n'.join(line.rstrip() for line in path.read_text(encoding='utf-8').splitlines())+'\n').encode())
    plt.close(fig)
    with (output/'source_data.csv').open('x', encoding='utf-8', newline='') as h:
        writer = csv.DictWriter(h, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    qa = dict(source_analysis_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest(),
        source_rows=len(rows), complete_three_strain_comparison=True, backend='Python/matplotlib',
        size_mm=[183, 155], all_axes_legend_bounds_inside=True, panel_labels='normal13pt(a)--(d)',
        interpolation_or_fit=False, statistical_replicates=False, energy_error_bars_inferred=False,
        marker_offsets_are_categorical_not_different_physical_strains=True,
        ordinary_pass_by_source=[r['ordinary_residual_pass'] for r in rows],
        old_labels_and_scores_unchanged=data['old_labels_and_scores_unchanged'],
        no_global_MEP_TS_or_ordering_certificate=True)
    with (output/'qa.json').open('x', encoding='utf-8', newline='\n') as h:
        json.dump(qa, h, indent=2, allow_nan=False); h.write('\n')
    return qa


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--source', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path); args = p.parse_args()
    print(json.dumps(draw(args.source, args.output)))
