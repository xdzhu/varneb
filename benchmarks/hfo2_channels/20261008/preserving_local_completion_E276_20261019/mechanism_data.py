"""Link the completed matched local controls to retained whole-path context."""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
from ase.io import read

REPO=Path(__file__).resolve().parents[4];sys.path.insert(0,str(REPO))
from examples.hfo2_fixed_input_factory import read_fixed_hfo2_stru
from scripts.audit_hfo2_static_replica import sha256
CASE=Path(__file__).resolve().parent;BASE=CASE.parent


def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def transverse_independence(seed):
    """Qy is a reference shuffle, not the transverse mean-O-minus-Hf proxy."""
    data=json.loads((seed/'seed.json').read_bytes());atoms=read_fixed_hfo2_stru(seed/'raw/STRU')
    t=read(BASE/'reference_variants/T.vasp')
    rotation=np.asarray(json.loads((BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json').read_bytes())
        ['orientation']['cartesian_rotation_rows'])
    transform=np.linalg.inv(atoms.cell.array)@rotation@t.cell.array
    weights=np.r_[np.full(4,-.25),np.full(8,.125)]
    py=np.r_[(weights[:,None]*transform[:,1]).ravel(),np.zeros(3)]
    c=np.asarray(data['covectors_39_by_4']);q=np.asarray(data['orthonormal_39_by_4'])
    rank=int(np.linalg.matrix_rank(np.column_stack([c,py])))
    if np.linalg.matrix_rank(c)!=4 or rank!=5:raise ValueError('transverse proxy must be independent of retained four observables')
    return dict(strain=data['substrate_strain'],retained_rank=4,with_transverse_proxy_rank=rank,
        transverse_covector=py.tolist(),released_transverse_covector_norm=float(np.linalg.norm(py-q@(q.T@py))),
        geometric_Qy_not_mean_polar_proxy_y=True,not_electronic_polarization_or_a_new_eigenmode=True)


def measured_row(strain,sign,e0,a,j,ma,mj,ga,gj,guard,job):
    if ma['measured_four_pattern_defect_A']>1e-12 or mj['measured_four_pattern_defect_A']>1e-12:
        raise ValueError('the original four physical observables must be independently preserved')
    if not 0<=ga<=.03 or not 0<=gj<=.03:raise ValueError('both local raw released-gradient targets required')
    if abs(ma['open_y_stretch_change_percentage_points'])>1e-10:raise ValueError('A must retain its entire original F cell')
    r=dict(strain_percent=strain,sign=sign,source_job=job,numerical_radius_A=guard,
        F_energy_eV_cell=e0,A_energy_eV_cell=a,J_energy_eV_cell=j,
        atomic_lowering_meV_fu=250*(e0-a),additional_joint_lowering_meV_fu=250*(a-j),
        total_lowering_meV_fu=250*(e0-j),A_released_gradient_eV_A=ga,J_released_gradient_eV_A=gj,
        atomic_norm_F_A_A=ma['atomic_motion_F_to_candidate']['total_norm_A'],
        atomic_norm_F_J_A=mj['atomic_motion_F_to_candidate']['total_norm_A'],
        oxygen_squared_displacement_fraction_A=ma['atomic_motion_F_to_candidate']['oxygen_squared_norm_fraction'],
        transverse_proxy_change_F_A_A=ma['structural_polar_proxy_change_A'][1],
        transverse_proxy_change_F_J_A=mj['structural_polar_proxy_change_A'][1],
        open_y_change_F_A_percentage_points=ma['open_y_stretch_change_percentage_points'],
        open_y_change_F_J_percentage_points=mj['open_y_stretch_change_percentage_points'])
    if not np.isfinite([v for v in r.values() if isinstance(v,(int,float))]).all():raise ValueError('finite native values required')
    return r


def complete_new_rows(source):
    data=json.loads(source.read_bytes())
    if (data['status']!='both_original_signed_continuations_replayed' or data['new_DFT_calls']
            or data['physical_inputs_changed'] or data['source_reader_sha256']!=sha256(CASE/'analyze_terminal.py')):
        raise ValueError('complete original native replay required, never fabricated target rows')
    if [r['sign'] for r in data['rows']]!=['minus','plus']:raise ValueError('both original signs required')
    rows=[]
    for row,job in zip(data['rows'],('28954771','28954772')):
        archive=CASE/'completed'/row['sign']/'native_terminal.tar.gz'
        if row['job_id']!=job or not row['all_new_native_EFS_reparsed_locally'] or sha256(archive)!=row['archive_sha256']:
            raise ValueError('actual original terminal native archive required')
        if [r['arm'] for r in row['arms']]!=['A','J'] or any(r['candidate'] is None for r in row['arms']):
            raise ValueError('qualified A and J required; never substitute a budget-failed last trial')
        a,j=[r['candidate'] for r in row['arms']]
        rows.append(measured_row(1.,row['sign'],row['original_F_energy_eV_cell'],a['energy_eV_cell'],j['energy_eV_cell'],
            a['actual_motion'],j['actual_motion'],a['released_gradient_eV_A'],j['released_gradient_eV_A'],1.,job))
    return rows,data


def build():
    source=CASE/'complete_analysis/analysis.json'
    # Missing/live targets fail before re-reading old cohorts, and are never filled.
    new,data=complete_new_rows(source)
    oldcase=BASE/'preserving_pattern_branches_E272_20261018'
    old=load(oldcase/'analyze_terminal.py','E276_zero_controls').analyze()
    reader=load(oldcase/'analyze_prefix.py','E276_zero_F')
    _,_,e0,_=reader.setup(0);rows=[]
    for r in old['rows']:
        if r['index']!=0:continue
        a,j=r['arms']
        if a['energy_eV_cell'] is None or j['energy_eV_cell'] is None:raise ValueError('original zero controls required')
        rows.append(measured_row(0.,r['sign'],e0,a['energy_eV_cell'],j['energy_eV_cell'],
            a['actual_geometry'],j['actual_geometry'],a['released_gradient_eV_A'],j['released_gradient_eV_A'],.5,'28948425'))
    rows+=new
    if [(r['strain_percent'],r['sign']) for r in rows]!=[(0.,'minus'),(0.,'plus'),(1.,'minus'),(1.,'plus')]:
        raise ValueError('exact zero/+1% signed comparison required')
    context=BASE/'preserving_strain_E275_20261019/completed_analysis_full/analysis.json'
    if sha256(context)!='829e524e083b9622ddb41da9a0029b0246deeaa05e97b88733f59aef1738ed7f':
        raise ValueError('retained whole-path response source changed')
    path=json.loads(context.read_bytes())
    for row in rows:
        strain=row['strain_percent']/100
        original=next(b for b in path['bands'] if b['strain']==strain and b['kind']=='registered_preserving')
        connected=next(b for b in path['bands'] if b['strain']==strain and b['kind']=='released_seed')
        old_peak=original['rows'][original['sampled_peak_image']]
        new_peak=connected['rows'][connected['sampled_peak_image']]
        outcome=next(r for r in path['outcomes'] if r['strain']==strain)
        initial=original['endpoint_native_records'][0]['energy_eV_cell']
        if not np.isclose(250*(row['F_energy_eV_cell']-initial),original['sampled_polygon_peak_meV_fu'],atol=1e-7,rtol=0):
            raise ValueError('local F must be the same actual original sampled bottleneck, with the same initial well')
        row.update(J_point_energy_relative_IS_meV_fu=250*(row['J_energy_eV_cell']-initial),
            connected_free_chain_sampled_peak_meV_fu=outcome['released_seed_meV_fu'],
            J_point_minus_connected_free_peak_meV_fu=250*(row['J_energy_eV_cell']-initial)-outcome['released_seed_meV_fu'])
        row.update(retained_Qx_at_original_F_A=old_peak['Qx_A'],connected_free_peak_Qx_A=new_peak['Qx_A'],
            Qx_original_F_to_free_peak_change_A=new_peak['Qx_A']-old_peak['Qx_A'],
            connected_free_peak_transverse_proxy_A=new_peak['polar_proxy_y_A'],
            connected_free_peak_open_y_percent=new_peak['open_y_stretch_percent'])
    return dict(rows=rows,new_source_sha256=sha256(source),zero_native_archive_sha256=old['archive_sha256'],
        descriptor_independence=[transverse_independence(oldcase/'prepared/seeds/0'),transverse_independence(CASE/'prepared/seed')],
        whole_path_context_sha256=sha256(context),whole_path_outcomes=path['outcomes'],
        new_native_calls=data['new_calls'],new_DFT_calls=0,all_four_triplets_A_J_qualified=True,
        original_plus1_failure_preserved=True,old_and_new_numerical_radii_differ=True,
        held_observable_definitions_same_native_values_different_across_strain=True,
        local_energy_is_not_corrected_path_barrier_or_unique_partition=True,
        displacement_norms_and_fractions_are_not_energy_fractions=True,
        J_point_to_free_peak_difference_is_not_an_additive_causal_mode_energy=True,
        statistical_replicates_or_inferred_error_bars=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();d=build();a.output.mkdir(parents=True,exist_ok=False)
    with (a.output/'analysis.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
    with (a.output/'source_data.csv').open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(d['rows'][0]));w.writeheader();w.writerows(d['rows'])
    print(json.dumps(d['rows']))
