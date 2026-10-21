"""Connect archived atom motion, complete paths and the final crossed controls.

No calculator or network is invoked. Reference-metric atom projections retain
ordered identities and the existing periodic lift; arrows are not magnified.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from ase.io import read

BASE = ROOT / 'benchmarks/hfo2_channels/20261008'
INPUTS = {
    'motion': BASE / 'preserving_connection_E273_20261019/atom_motion/analysis.json',
    'paths': BASE / 'preserving_strain_E275_20261019/completed_analysis_full/analysis.json',
    'cross': BASE / 'strain_pattern_cross_E277_20261021/analysis/analysis.json',
    'reference': BASE / 'reference_variants/T.vasp',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def evidence():
    motion, paths, cross = (json.loads(INPUTS[k].read_bytes()) for k in ('motion', 'paths', 'cross'))
    if any(d['new_DFT_calls'] for d in (motion, paths, cross)):
        raise ValueError('archived zero-new-DFT analyses required')
    if cross['whole_path_context_sha256'] != sha(INPUTS['paths']):
        raise ValueError('cross and complete-path references differ')
    if motion['reference_sha256'] != sha(INPUTS['reference']):
        raise ValueError('original-T metric changed')
    old = np.asarray(motion['bands']['original']['fractional_positions'])
    low = np.asarray(motion['bands']['lower']['fractional_positions'])
    cell = read(INPUTS['reference']).cell.array
    delta = (low[4] - old[4]) @ cell
    delta -= delta.mean(axis=0)
    np.testing.assert_allclose(delta, motion['center_lower_minus_original_translation_free_A'], atol=1e-12)
    np.testing.assert_allclose(delta.mean(axis=0), 0., atol=1e-14)
    # One registration inherited from the common initial state, never image-wise wrapping.
    gauge = np.floor(old[0] + 1e-10)
    start = (old[4] - gauge) @ cell
    end = start + delta
    selected = [{k: b[k] for k in ('strain', 'kind', 'sampled_polygon_peak_meV_fu', 'fmax_eV_A')}
                for b in paths['bands'] if b['kind'] in ('registered_preserving', 'released_seed', 'T_escape')]
    expected = {(s, k) for s in (0., .005, .01)
                for k in ('registered_preserving', 'released_seed', 'T_escape')}
    if {(r['strain'], r['kind']) for r in selected} != expected or len(selected) != 9:
        raise ValueError('three complete candidate categories at all three strains required')
    if any(r['fmax_eV_A'] > .10 for r in selected):
        raise ValueError('ordinary completed candidates required')
    rows = cross['rows']
    if len(rows) != 4 or cross['contrasts']['status'] != 'four_candidate_contrasts':
        raise ValueError('actual four-condition local controls required')
    margins = {k: motion['bands'][k]['minimum_layer_band_margin_fractional_z'] for k in ('original', 'lower')}
    return dict(sources={k: {'path': p.relative_to(ROOT).as_posix(), 'sha256': sha(p)} for k,p in INPUTS.items()},
                ordered_symbols=motion['atom_order'], old_center_reference_A=start.tolist(),
                lower_center_aligned_reference_A=end.tolist(), center_translation_free_difference_A=delta.tolist(),
                inherited_integer_registration=gauge.astype(int).tolist(), path_rows=selected, cross_rows=rows,
                cross_contrasts=cross['contrasts'], minimum_ordered_layer_margins=margins,
                Hf_max_abs_dx_A=float(abs(delta[:4, 0]).max()),
                O_max_abs_dz_A=float(abs(delta[4:, 2]).max()),
                new_DFT_calls=0, atom_order_changed=False, imagewise_remapping=False,
                atom_arrows_magnified=False, affine_cell_motion_in_atom_arrows=False,
                literature_path_identity_certified=False, TS_or_global_barrier_ordering_certified=False)


def draw(data, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':12, 'axes.labelsize':13, 'xtick.labelsize':11, 'ytick.labelsize':11,
        'legend.fontsize':10.5, 'pdf.fonttype':42, 'svg.fonttype':'none', 'axes.linewidth':1.0})
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.8))
    fig.subplots_adjust(left=.10, right=.97, bottom=.09, top=.89, wspace=.34, hspace=.40)
    old, new, delta = (np.asarray(data[k]) for k in ('old_center_reference_A',
        'lower_center_aligned_reference_A', 'center_translation_free_difference_A'))
    ax = axes[0,0]
    for atom in range(12):
        color = '#78818C' if atom < 4 else '#B95D54'
        marker = 's' if atom < 4 else 'o'
        ax.scatter(old[atom,1], old[atom,2], s=72, marker=marker, facecolors='white', edgecolors=color, linewidth=1.2, zorder=3)
        ax.annotate('', xy=new[atom,[1,2]], xytext=old[atom,[1,2]],
                    arrowprops={'arrowstyle':'->', 'color':color, 'lw':1.25}, zorder=4)
        ax.scatter(new[atom,1], new[atom,2], s=32, marker=marker, color=color, zorder=5)
    ax.set(xlabel=r'$y$ in reference metric ($\mathrm{\AA}$)', ylabel=r'$z$ in reference metric ($\mathrm{\AA}$)',
           ylim=(-.4,4.4))
    # Equal physical units with unchanged panel boxes; limits adjust at draw time.
    ax.set_aspect('equal', adjustable='datalim')

    ax = axes[0,1]
    x = np.arange(12)
    for j,(label,color) in enumerate(zip(('x','y','z'),('#326D95','#D18D48','#69917B'))):
        ax.bar(x+(j-1)*.24, delta[:,j], width=.22, color=color, label=r'$\Delta '+label+'$')
    ax.axhline(0, color='#687079', lw=.65)
    ax.axvline(3.5, color='#B7BCC2', lw=.8, linestyle=':')
    ax.set(xlabel='Ordered atom index (1–12)', ylabel=r'Lower minus original ($\mathrm{\AA}$)', ylim=(-.26,.31),
           xticks=x, xticklabels=[str(i+1) for i in x], xlim=(-.7,11.7))
    ax.text(1.5,.275,'Hf',ha='center',va='center',fontsize=11)
    ax.text(7.5,.275,'O',ha='center',va='center',fontsize=11)
    ax.legend(loc='lower center', ncol=3, framealpha=.85, facecolor='white', edgecolor='#B0B5BB')

    ax = axes[1,0]
    specs = [('registered_preserving','Original preserving','#326D95','o'),
             ('released_seed','Lower preserving','#C87F3B','s'),('T_escape','T escape','#727A84','D')]
    for kind,label,color,marker in specs:
        rows=sorted((r for r in data['path_rows'] if r['kind']==kind),key=lambda r:r['strain'])
        ax.plot([r['strain']*100 for r in rows],[r['sampled_polygon_peak_meV_fu'] for r in rows],
                marker=marker,color=color,lw=1.7,markersize=6,label=label)
    ax.set(xlabel='Imposed plane strain (%)',ylabel=r'Sampled peak (meV/HfO$_2$)',
           xlim=(-.12,1.12),ylim=(76,143),xticks=[0,.5,1.])
    ax.legend(loc='upper left',framealpha=.85,facecolor='white',edgecolor='#B0B5BB')

    ax = axes[1,1]
    offsets={(0.,'q0'):(12,-10),(.01,'q0'):(12,4),(0.,'q1'):(-12,8),(.01,'q1'):(12,-4)}
    for row in data['cross_rows']:
        key=(row['strain'],row['held_pattern']); color='#326D95' if key[1]=='q0' else '#C87F3B'
        marker='o' if key[1]=='q0' else 's'
        ax.scatter(row['open_y_stretch_percent'],row['transverse_proxy_A'],s=70,marker=marker,
                   facecolors='white' if key[0]==0 else color,edgecolors=color,linewidth=1.5,zorder=4)
        shift=offsets[key]
        label=f"{'0' if key[0]==0 else '+1'}%, {key[1]}\n{row['relative_energy_meV_fu']:.2f}"
        ax.annotate(label,(row['open_y_stretch_percent'],row['transverse_proxy_A']),xytext=shift,
                    textcoords='offset points',ha='right' if shift[0]<0 else 'left',
                    va='bottom' if shift[1]>0 else 'top',fontsize=10.5)
    ax.set(xlabel='Open-$y$ stretch (%)',ylabel=r'Transverse O–Hf proxy ($\mathrm{\AA}$)',
           xlim=(-.72,.83),ylim=(.085,.255))
    for tag,ax in zip('abcd',axes.ravel()):
        ax.tick_params(direction='in', top=True, right=True, length=4)
        ax.text(-.14,1.04,f'({tag})',transform=ax.transAxes,fontsize=16,fontweight='normal',va='bottom')
        for spine in ax.spines.values():spine.set_visible(True)
    handles=[Line2D([],[],marker='s',ls='',color='#78818C',label='Hf'),
             Line2D([],[],marker='o',ls='',color='#B95D54',label='O'),
             Line2D([],[],marker='o',ls='',color='#59616B',markerfacecolor='white',label='Original center'),
             Line2D([],[],marker='o',ls='',color='#59616B',label='Lower center')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,.99),ncol=4,
               framealpha=.85,facecolor='white',edgecolor='#B0B5BB',fontsize=11)
    fig.align_ylabels(axes[:,0]);fig.align_ylabels(axes[:,1])
    fig.canvas.draw()
    positions=[list(ax.get_position().bounds) for ax in axes.ravel()]
    assert np.isclose(positions[0][0],positions[2][0]) and np.isclose(positions[1][0],positions[3][0])
    assert np.isclose(positions[0][1],positions[1][1]) and np.isclose(positions[2][1],positions[3][1])
    for ext in ('pdf','svg','png'):fig.savefig(output/f'hfo2_mechanism_synthesis.{ext}',dpi=250,facecolor='white')
    plt.close(fig)
    return dict(source_script_sha256=sha(Path(__file__)),panel_positions=positions,
                panel_labels_normal_weight=True,top_and_right_axes=True,
                plotted_atom_arrows_scale=1.,reference_projection_equal_unit_aspect=True,
                physical_distances_from_native_cells=False,
                source_data='analysis.json',statistical_replicates_or_error_bars=False,
                visual_inspection_required=True)


def main(output):
    data=evidence();output=Path(output);output.mkdir(parents=True,exist_ok=False)
    with (output/'analysis.json').open('x',encoding='utf-8',newline='\n') as handle:
        json.dump(data,handle,indent=2,allow_nan=False);handle.write('\n')
    with (output/'center_motion.csv').open('x',encoding='utf-8',newline='') as handle:
        writer=csv.writer(handle);writer.writerow(['atom','species','dx_A','dy_A','dz_A'])
        for i,(s,v) in enumerate(zip(data['ordered_symbols'],data['center_translation_free_difference_A'])):
            writer.writerow([i+1,s]+v)
    qa=draw(data,output)
    with (output/'qa.json').open('x',encoding='utf-8',newline='\n') as handle:
        json.dump(qa,handle,indent=2);handle.write('\n')
    print(json.dumps(dict(output=str(output),new_DFT_calls=0,Hf_max_dx=data['Hf_max_abs_dx_A'],O_max_dz=data['O_max_abs_dz_A'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
