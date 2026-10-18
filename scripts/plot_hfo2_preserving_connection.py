"""Same-parameter free-chain snapshots: energy, structure and NEB residual.

Figure contract: a released central configuration is not a relaxed whole
path. Four quantitative panels link its local gain and transverse motion
to neighboring force imbalance. One chain per initialization, no statistical
replicates, fitted smoothing, barrier error bars or acceleration claim.
"""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tarfile
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read
from examples.hfo2_fixed_input_factory import CONTRACT,read_fixed_hfo2_stru,same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results,sha256
from scripts.analyze_hfo2_clamped_residual import decompose
from scripts.plot_hfo2_mechanism import original_chart_images,path_coordinates,rotated_t_triplet
from vcneb import VCNEB,clamped_plane_vcneb_boundary,validate_path_geometry,validate_periodic_path_lift

BASE=ROOT/'benchmarks/hfo2_channels/20261008'
CASE=BASE/'preserving_connection_E273_20261019'
PARENT=BASE/'preserving_pattern_branches_E272_20261018'
SCALE=5.12968067458423
PIN='a452d5b125d6d95c729ea60607524ce5d35370d5aab57ac6f2b1b7869d61991a'


def restore(step,destination):
    archive=CASE/f'snapshot_step_{step:04d}.tar.gz'
    receipt=json.loads((CASE/f'snapshot_step_{step:04d}_receipt.json').read_bytes())
    if (sha256(archive)!=receipt['archive_sha256'] or receipt['job_id']!='28949435'
            or receipt['step']!=step or receipt['new_DFT_calls']!=0 or receipt['scheduler_mutations']!=0):
        raise ValueError('exact selected native-frame archive required')
    expected=f'step_{step:04d}'
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            p=Path(member.name)
            if (p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0]!=expected
                    or not(member.isfile() or member.isdir())
                    or not(destination/p).resolve().is_relative_to(destination.resolve())):
                raise ValueError('unsafe native snapshot member')
        tar.extractall(destination)
    folder=destination/expected;files=json.loads((folder/'files.json').read_bytes())
    if len(files)!=receipt['exported_files']:raise ValueError('selected native file count differs')
    for name,digest in files.items():
        p=folder/name
        if not p.resolve().is_relative_to(folder.resolve()) or sha256(p)!=digest:
            raise ValueError('selected native bytes changed')
    return folder,receipt


def native_images(folder,step):
    receipt=json.loads((folder/'receipt.json').read_bytes());obsdir=folder/'observation'
    obs=json.loads((obsdir/'observation.json').read_bytes())
    if (receipt['job_id']!='28949435' or receipt['bundle_pin']!=PIN or receipt['step']!=step
            or obs['source_job_id']!='28949435' or obs['snapshot_step']!=step
            or obs['n_total_images']!=9 or obs['formula_units']!=4 or obs['pressure_GPa']!=0 or obs['climb']
            or obs['cell_scale_A']!=SCALE or obs['k_eV_A2']!=.2 or obs['fmax_target_eV_A']!=.10
            or sha256(obsdir/'runtime_preflight.json')!=obs['runtime_preflight_sha256']):
        raise ValueError('same ordinary nine-image zero-pressure input contract required')
    images=[]
    for i,row in enumerate(obs['raw_image_evaluations']):
        poscar=obsdir/f'POSCAR_{i:02d}';raw=obsdir/'raw'/f'image_{i:04d}'
        if (row['image_index']!=i or sha256(poscar)!=row['snapshot_POSCAR_sha256']
                or any(row['input_sha256'][n]!=h for n,h in CONTRACT.items())
                or any(sha256(raw/n)!=row['input_sha256'][n] for n in ('INPUT','KPT','STRU'))
                or sha256(raw/'OUT.ABACUS/running_scf.log')!=row['raw_log_sha256']):
            raise ValueError('native physical input/geometry/output changed')
        image=read(poscar,format='vasp');parsed=audited_results(raw)
        values=dict(energy=row['energy_eV_cell'],forces=row['forces_eV_A'],stress=row['stress_ASE_voigt_eV_A3'])
        if (not same_ordered_geometry(image,read_fixed_hfo2_stru(raw/'STRU'))
                or any(not np.allclose(parsed[k],values[k],atol=1e-12,rtol=0) for k in parsed)):
            raise ValueError('native ordered energy/forces/stress replay differs')
        image.calc=SinglePointCalculator(image,**parsed);images.append(image)
    if len(images)!=9:raise ValueError('all nine evaluated geometries required')
    return images,obs,receipt


def band(images,reference,t,parent,patterns,orientation,label):
    boundary=clamped_plane_vcneb_boundary(12,reference,allow_tilt=True)
    boundary.validate_images(images);validate_periodic_path_lift(images)
    geometry=validate_path_geometry(images,cell_scale=SCALE,minimum_distance=1.6,maximum_deformation=.25)
    chain=VCNEB(images,cell_scale=SCALE,k=.2,climb=False,pressure=0.,**boundary.vcneb_kwargs(images))
    force=decompose(chain);arc=chain.reaction_coordinate();relative=(chain.enthalpies-chain.enthalpies[0])*250
    # Full vectors are checked inside decompose and remain reproducible from
    # the native archive; retain compact physical diagnostics in this plot.
    force['images']=[{k:v for k,v in r.items() if not k.endswith('_generalized_vector_eV_A')}
                     for r in force['images']]
    q,polar,_,_=path_coordinates(original_chart_images(images,orientation,t),parent,patterns)
    rows=[]
    for i in range(9):
        f=None if i in (0,8) else force['images'][i-1]['residual_max_vector_eV_A']
        rows.append(dict(series=label,image=i,s=float(arc[i]/arc[-1]),energy_meV_fu=float(relative[i]),
            Qx_A=float(q[i,0]),Qy_A=float(q[i,1]),Qz_A=float(q[i,2]),
            polar_proxy_x_A=float(polar[i,0]),polar_proxy_y_A=float(polar[i,1]),polar_proxy_z_A=float(polar[i,2]),
            NEB_residual_eV_A=f))
    return dict(label=label,rows=rows,fmax_eV_A=force['fmax_eV_A'],ordinary_residual_pass=force['fmax_eV_A']<=.10,
        sampled_polygon_peak_meV_fu=float(relative.max()),sampled_peak_image=int(relative.argmax()),
        residual_decomposition=force,geometry_gate=geometry)


def build_data(step=1):
    if not 1<=step<=24:raise ValueError('one actual post-initialization complete update required')
    spec=importlib.util.spec_from_file_location('E273_original_native_frames',PARENT/'analyze_prefix.py')
    original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
    _,seed=original.analyze();old,old_obs=original.original_images()
    config=json.loads((BASE/'matched_pattern_release_E270_20261018/prepared/config.json').read_bytes())
    tpath=BASE/'reference_variants/T.vasp';frame=BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json'
    if sha256(tpath)!=config['T_reference_sha256'] or sha256(frame)!=config['orientation_sha256']:
        raise ValueError('the common original-T metric changed')
    t=read(tpath);parent,patterns,_=rotated_t_triplet(t);orientation=json.loads(frame.read_bytes())['orientation']
    reference=old_obs['mechanical_boundary']['reference_cell_A']
    with tempfile.TemporaryDirectory(prefix='varneb_connection_frame_') as tmp:
        folder,outer=restore(step,Path(tmp));current,obs,receipt=native_images(folder,step)
        for i in (0,8):
            if (not np.allclose(current[i].get_scaled_positions(wrap=False),old[i].get_scaled_positions(wrap=False),atol=1e-12,rtol=0)
                    or not np.array_equal(current[i].calc.results['forces'],old[i].calc.results['forces'])
                    or current[i].get_potential_energy()!=old[i].get_potential_energy()):
                raise ValueError('fixed endpoint ordered lift and native cache changed')
        bands=[band(images,reference,t,parent,patterns,orientation,label) for images,label in
               ((old,'Original band'),(seed,'Released-center seed'),(current,f'Free-chain step {step}'))]
        if abs(bands[-1]['fmax_eV_A']-obs['replayed_fmax_eV_A'])>1e-10:
            raise ValueError('local full-chain force replay differs from native snapshot')
    return dict(status='actual_free_chain_progress_not_a_final_barrier',step=step,bands=bands,
        source_archive_sha256=outer['archive_sha256'],source_observation_sha256=original.OBSERVATION_SHA,
        same_physical_parameters_and_endpoint_caches=True,full_orbital_bytes_rechecked_locally=False,
        new_DFT_calls=0,old_labels_and_forecast_scores_unchanged=True,stationary_TS_or_sampling_bound_certified=False,
        current_scheduler_row_at_export=receipt['scheduler_row_at_export'],selected_new_SCF_calls=len(receipt['selected_fresh_calls']),
        selected_frame_native_proxy_core_hours=receipt['selected_frame_native_proxy_core_hours'],earlier_costs_additional=True,
        structural_proxies_are_not_electronic_polarization=True,one_chain_not_statistical_replicates=True)


def draw(data,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':12,'axes.labelsize':12,'xtick.labelsize':11.5,'ytick.labelsize':11.5,
        'legend.fontsize':11,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.9})
    fig,axes=plt.subplots(2,2,figsize=(183/25.4,155/25.4))
    fig.subplots_adjust(left=.13,right=.98,bottom=.13,top=.87,wspace=.39,hspace=.48)
    styles=[dict(color='#306A94',marker='o',linestyle='-'),dict(color='#D48345',marker='^',linestyle='--',markerfacecolor='white'),
            dict(color='#3C8C78',marker='D',linestyle='-')]
    for series,style in zip(data['bands'],styles):
        rows=series['rows'];kw=dict(style,markersize=4,linewidth=1.3)
        for ax,key in zip(axes.flat,('energy_meV_fu','polar_proxy_y_A','polar_proxy_z_A','NEB_residual_eV_A')):
            selected=rows[1:-1] if key=='NEB_residual_eV_A' else rows
            ax.plot([r['s'] for r in selected],[r[key] for r in selected],**kw)
    labels=['Energy (meV/HfO$_2$)','Polar y proxy (Å)','Polar z proxy (Å)','NEB residual (eV/Å)']
    for ax,label,letter in zip(axes.flat,labels,'abcd'):
        ax.set_xlim(-.035,1.035);ax.set_xticks([0,.5,1]);ax.set_xlabel('Normalized generalized arc',labelpad=6)
        ax.set_ylabel(label,labelpad=6);ax.tick_params(top=True,right=True,direction='in')
        for spine in ax.spines.values():spine.set_visible(True)
        ax.text(0,1.055,f'({letter})',transform=ax.transAxes,fontsize=13,fontweight='normal',va='bottom')
    axes[0,0].set_ylim(-5,105);axes[0,0].set_yticks([0,50,100])
    axes[0,1].set_ylim(-.01,.14);axes[0,1].set_yticks([0,.05,.10])
    axes[1,1].set_ylim(0,.155);axes[1,1].set_yticks([0,.05,.10,.15])
    target=axes[1,1].axhline(.10,color='#686C70',linewidth=.9,linestyle=':')
    axes[1,1].legend([target],['0.10 target'],loc='upper right',frameon=True,framealpha=.85,
        facecolor='white',edgecolor='#A8ADB2',borderpad=.3,handlelength=1.3)
    handles=[Line2D([],[],**style,markersize=4,linewidth=1.3) for style in styles]
    fig.legend(handles,[b['label'] for b in data['bands']],loc='upper center',bbox_to_anchor=(.55,.992),ncol=3,
        frameon=True,framealpha=.85,facecolor='white',edgecolor='#A8ADB2',borderpad=.4,columnspacing=.9,handlelength=1.5)
    fig.canvas.draw();bounds=[ax.get_position().bounds for ax in axes.flat]
    for row in (bounds[:2],bounds[2:]):
        if not np.allclose(np.asarray(row)[:,[1,3]],np.asarray(row)[0,[1,3]],atol=1e-12,rtol=0):
            raise ValueError('same-row panel axes are not aligned')
    for column in ((bounds[0],bounds[2]),(bounds[1],bounds[3])):
        if not np.allclose(np.asarray(column)[:,[0,2]],np.asarray(column)[0,[0,2]],atol=1e-12,rtol=0):
            raise ValueError('same-column panel axes are not aligned')
    renderer=fig.canvas.get_renderer();width,height=fig.canvas.get_width_height()
    for artist in list(fig.legends)+[a.get_legend() for a in axes.flat if a.get_legend()]:
        box=artist.get_window_extent(renderer)
        if min(box.x0,box.y0)<0 or box.x1>width or box.y1>height:
            raise ValueError('legend clipped outside the figure')
    for extension in ('pdf','svg','png'):fig.savefig(output/f'connection_step_{data["step"]:04d}.{extension}',dpi=300,facecolor='white')
    plt.close(fig)
    return dict(archetype='quantitative grid',size_mm=[183,155],axes_bounds=bounds,backend='Python/matplotlib',
        panel_letter_weight='normal',panel_letter_size_pt=13,top_and_right_spines=True,subpanel_titles=False,
        legend_white_background_alpha=.85,interpolation_or_error_bars_or_synthetic_points=False,
        intended_role='mechanistic progress only; not a final or acceleration benchmark')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--step',type=int,default=1)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    data=build_data(args.step);args.output.mkdir(parents=True,exist_ok=True)
    with (args.output/'source_data.csv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(data['bands'][0]['rows'][0]));writer.writeheader()
        writer.writerows(row for b in data['bands'] for row in b['rows'])
    qa=draw(data,args.output);qa.update(plot_source_sha256=sha256(Path(__file__)),source_data_sha256=sha256(args.output/'source_data.csv'))
    for name,value in (('analysis.json',data),('figure_qa.json',qa)):
        (args.output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:data[k] for k in ('status','step','new_DFT_calls')}))
