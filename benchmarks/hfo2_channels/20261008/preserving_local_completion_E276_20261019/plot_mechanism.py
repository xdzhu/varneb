"""Measured local energy, mixed motion and open-cell response; no imputation."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import numpy as np

REPO=Path(__file__).resolve().parents[4];sys.path.insert(0,str(REPO))
from scripts.audit_hfo2_static_replica import sha256


def figure_rows(data):
    rows=data['rows']
    if not data.get('all_four_triplets_A_J_qualified') or len(rows)!=4:
        raise ValueError('four actual qualified triplets required, no missing-point substitution')
    if [(r['strain_percent'],r['sign']) for r in rows]!=[(0.,'minus'),(0.,'plus'),(1.,'minus'),(1.,'plus')]:
        raise ValueError('both original signs and both original mechanical conditions required')
    fields=('atomic_lowering_meV_fu','additional_joint_lowering_meV_fu','atomic_norm_F_A_A',
            'atomic_norm_F_J_A','transverse_proxy_change_F_A_A','transverse_proxy_change_F_J_A',
            'open_y_change_F_A_percentage_points','open_y_change_F_J_percentage_points')
    for r in rows:
        if not np.isfinite([r[k] for k in fields]).all():raise ValueError('finite native energy/geometry required')
        if not 0<=r['A_released_gradient_eV_A']<=.03 or not 0<=r['J_released_gradient_eV_A']<=.03:
            raise ValueError('ordinary NEB convergence cannot replace the local gradient criterion')
    return rows


def draw(source,output):
    source,output=Path(source),Path(output)
    data=json.loads(source.read_bytes());rows=figure_rows(data)
    if output.exists():raise FileExistsError(output)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':11.5,'axes.labelsize':12,'axes.labelweight':'normal',
        'xtick.labelsize':11,'ytick.labelsize':11,'legend.fontsize':10.5,
        'axes.linewidth':.85,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(2,2,figsize=(183/25.4,160/25.4))
    fig.subplots_adjust(left=.145,right=.98,bottom=.13,top=.82,wspace=.46,hspace=.53)
    x=np.arange(4);blue,orange='#306A94','#D17A40'
    atomic=np.array([r['atomic_lowering_meV_fu'] for r in rows])
    joint=np.array([r['additional_joint_lowering_meV_fu'] for r in rows])
    axes[0,0].bar(x,atomic,.54,color=blue,edgecolor='#333333',lw=.6,label=r'$F-A$')
    axes[0,0].bar(x,joint,.54,bottom=atomic,color=orange,edgecolor='#333333',lw=.6,label=r'$A-J$')
    axes[0,0].set_ylabel('Energy lowering\n(meV/HfO$_2$)')
    axes[0,1].bar(x-.15,[r['atomic_norm_F_A_A'] for r in rows],.27,color=blue,edgecolor='#333333',lw=.6,label=r'$F\to A$')
    axes[0,1].bar(x+.15,[r['atomic_norm_F_J_A'] for r in rows],.27,color=orange,edgecolor='#333333',lw=.6,label=r'$F\to J$')
    axes[0,1].set_ylabel('Atomic displacement\nnorm (Å)')
    for ax,prefix in ((axes[1,0],'transverse_proxy_change'),(axes[1,1],'open_y_change')):
        suffix='_A' if prefix.startswith('transverse') else '_percentage_points'
        ax.plot(x-.08,[r[prefix+'_F_A'+suffix] for r in rows],'o',ls='none',color=blue,ms=5.5)
        ax.plot(x+.08,[r[prefix+'_F_J'+suffix] for r in rows],'D',ls='none',color=orange,ms=5.5)
        ax.axhline(0,color='#888888',lw=.7,zorder=0)
    axes[1,0].set_ylabel('Transverse proxy change\nrelative to F (Å)')
    axes[1,1].set_ylabel('Open-$y$ response\nrelative to F (pp)')
    legends=[];letters=[]
    for ax in axes[0]:
        legends.append(ax.legend(loc='lower center',bbox_to_anchor=(.54,1.13),ncol=2,
            frameon=True,facecolor='white',edgecolor='#A8ADB2',framealpha=.85,
            handlelength=1.1,handletextpad=.4,columnspacing=.8,borderpad=.35))
    for ax,letter in zip(axes.ravel(),'abcd'):
        ax.set_xticks(x,['0, −','0, +','1, −','1, +']);ax.set_xlim(-.6,3.6)
        ax.set_xlabel('Strain (%), start sign',labelpad=6)
        ax.tick_params(which='both',top=True,right=True,direction='in',length=4)
        for spine in ax.spines.values():spine.set_visible(True)
        ax.yaxis.set_label_coords(-.20,.5)
        letters.append(ax.text(0,1.055,'('+letter+')',transform=ax.transAxes,fontsize=13,fontweight='normal',va='bottom'))
    fig.canvas.draw();renderer=fig.canvas.get_renderer();width,height=fig.canvas.get_width_height()
    for artist in legends+letters+[a.xaxis.label for a in axes.ravel()]+[a.yaxis.label for a in axes.ravel()]:
        box=artist.get_window_extent(renderer)
        if box.x0<0 or box.y0<0 or box.x1>width or box.y1>height:raise ValueError('figure label/legend clipped')
    bounds=np.array([a.get_position().bounds for a in axes.ravel()])
    if not np.allclose(bounds[[0,2],0],bounds[0,0]) or not np.allclose(bounds[[1,3],0],bounds[1,0]):
        raise ValueError('column axes not aligned')
    if not np.allclose(bounds[[0,1],1],bounds[0,1]) or not np.allclose(bounds[[2,3],1],bounds[2,1]):
        raise ValueError('row axes not aligned')
    output.mkdir(parents=True)
    for ext in ('pdf','svg','png'):fig.savefig(output/('matched_preserving_release.'+ext),dpi=300,facecolor='white')
    plt.close(fig)
    # Canonical SVG whitespace does not alter appearance or scientific data.
    svg=output/'matched_preserving_release.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    qa=dict(source_analysis_sha256=sha256(source),source_rows=4,backend='Python/matplotlib',
        size_mm=[183,160],axes_bounds=bounds.tolist(),all_axes_labels_and_legends_inside=True,
        panel_labels='normal13pt(a)--(d)',top_and_right_axes=True,subpanel_titles=False,
        old_radius_A=.5,new_radius_A=1.,interpolation_or_fit=False,
        statistical_replicates=False,energy_error_bars_from_force=False,
        local_points_not_corrected_barriers_or_full_stability_certificates=True)
    with (output/'qa.json').open('x',encoding='utf-8',newline='\n') as f:json.dump(qa,f,indent=2);f.write('\n')
    return qa


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();print(json.dumps(draw(a.source,a.output)))
