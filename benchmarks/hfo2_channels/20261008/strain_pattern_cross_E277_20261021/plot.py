"""Four measured cross conditions; no interpolation or statistical replicas."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
REPO=Path(__file__).resolve().parents[4];sys.path.insert(0,str(REPO))
from scripts.audit_hfo2_static_replica import sha256


def plot(source,output):
    source=Path(source);d=json.loads(source.read_bytes());rows=d['rows']
    if (d['contrasts']['status']!='four_candidate_contrasts' or len(rows)!=4
            or any(r.get('relative_energy_meV_fu') is None for r in rows)
            or d['physical_inputs_changed'] or d['new_DFT_calls']):
        raise ValueError('all four actual native-replayed local candidates required')
    by={(r['strain'],r['held_pattern']):r for r in rows}
    if set(by)!={(s,q) for s in (0.,.01) for q in ('q0','q1')}:raise ValueError('complete crossed design required')
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':13,'axes.labelsize':14,'xtick.labelsize':12,'ytick.labelsize':12,
        'legend.fontsize':12,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':1.0})
    fig,axes=plt.subplots(2,2,figsize=(10.4,8.1))
    fig.subplots_adjust(left=.12,right=.97,bottom=.10,top=.89,wspace=.35,hspace=.37)
    specs=[(axes[0,0],'relative_energy_meV_fu',r'$E-E_{\rm IS}$ (meV/HfO$_2$)'),
           (axes[1,0],'transverse_proxy_A',r'Transverse displacement ($\mathrm{\AA}$)'),
           (axes[1,1],'open_y_stretch_percent',r'Open-$y$ stretch (%)')]
    handles=[]
    for ax,key,label in specs:
        for q,color,marker in [('q0','#306A94','o'),('q1','#C77C3C','s')]:
            line,=ax.plot([0,1],[by[s,q][key] for s in (0.,.01)],color=color,marker=marker,
                markersize=7,linewidth=1.7,label=r'$\mathbf{q}_{'+q[1]+r'}$ retained')
            if ax is axes[0,0]:handles.append(line)
        ax.set_ylabel(label);ax.margins(y=.18)
    ax=axes[0,1];contrast=d['contrasts']['fixed_strain_pattern_contrast_meV_fu']
    values=[contrast['0.0'],contrast['0.01']]
    ax.bar([0,1],values,width=.38,color='#74818F',edgecolor='#3D4854',linewidth=.7)
    ax.axhline(0,color='#6D7379',lw=.7)
    for x,v in zip((0,1),values):
        ax.annotate(f'{v:.2f}',(x,v),xytext=(0,5 if v>=0 else -6),textcoords='offset points',
            ha='center',va='bottom' if v>=0 else 'top',fontsize=12)
    ax.set_ylabel(r'$E(\mathbf{q}_1)-E(\mathbf{q}_0)$ (meV/HfO$_2$)')
    low=min(0.,*values);high=max(0.,*values);span=max(high-low,1.)
    ax.set_ylim(low-.22*span,high+.35*span)
    ax.text(.5,.96,f"Difference of contrasts: {d['contrasts']['interaction_meV_fu']:+.2f} meV/HfO$_2$",
        transform=ax.transAxes,ha='center',va='top',fontsize=11,
        bbox=dict(facecolor='white',edgecolor='#C3C7CA',alpha=.85,pad=4))
    for tag,ax in zip('abcd',axes.ravel()):
        ax.set_xlim(-.25,1.25);ax.set_xticks([0,.5,1.],['0','0.5','1.0'])
        ax.set_xlabel('Imposed plane strain (%)')
        ax.tick_params(direction='in',top=True,right=True,length=4)
        for spine in ax.spines.values():spine.set_visible(True)
        ax.text(-.16,1.06,f'({tag})',transform=ax.transAxes,fontsize=16,fontweight='normal',va='bottom')
    legend=fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.985),ncol=2,
        frameon=True,facecolor='white',edgecolor='#AAB0B5',framealpha=.85)
    fig.align_ylabels(axes[:,0]);fig.align_ylabels(axes[:,1]);fig.align_xlabels(axes[0]);fig.align_xlabels(axes[1])
    for ext in ('pdf','svg','png'):fig.savefig(output/f'strain_pattern_cross.{ext}',dpi=240,facecolor='white')
    fig.canvas.draw();positions=[list(ax.get_position().bounds) for ax in axes.ravel()]
    # Box alignment is checked numerically in addition to human visual review.
    assert np.isclose(positions[0][0],positions[2][0]) and np.isclose(positions[1][0],positions[3][0])
    assert np.isclose(positions[0][1],positions[1][1]) and np.isclose(positions[2][1],positions[3][1])
    qa=dict(source_sha256=sha256(source),plotter_sha256=sha256(Path(__file__)),
        measured_conditions=4,new_cross_conditions=2,statistical_replicates=False,
        connecting_lines_are_guides=True,interpolated_or_fitted_energies=False,
        panel_positions=positions,normal_weight_panel_labels=True,top_and_right_axes=True,
        source_rows=rows,visual_inspection_required=True)
    with (output/'qa.json').open('x',encoding='utf-8',newline='\n') as f:json.dump(qa,f,indent=2);f.write('\n')
    plt.close(fig);return qa


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();plot(a.source,a.output)
