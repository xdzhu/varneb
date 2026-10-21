"""Native replay and finite two-factor contrasts, never a new barrier estimate."""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import numpy as np
REPO=Path(__file__).resolve().parents[4];sys.path.insert(0,str(REPO))
from scripts.audit_hfo2_static_replica import sha256
CASE=Path(__file__).resolve().parent;BASE=CASE.parent


def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def contrasts(rows):
    values={(r['strain'],r['held_pattern']):r['relative_energy_meV_fu'] for r in rows}
    if len(rows)!=4 or set(values)!={(s,q) for s in (0.,.01) for q in ('q0','q1')}:
        raise ValueError('exactly the two by two measured conditions required')
    if any(v is None for v in values.values()):return dict(status='incomplete_no_factor_contrasts')
    strain={q:values[.01,q]-values[0.,q] for q in ('q0','q1')}
    pattern={str(s):values[s,'q1']-values[s,'q0'] for s in (0.,.01)}
    interaction=strain['q1']-strain['q0']
    if not np.isclose(interaction,pattern['0.01']-pattern['0.0'],atol=1e-12,rtol=0):
        raise ValueError('two finite contrast orders must agree')
    return dict(status='four_candidate_contrasts',fixed_pattern_strain_response_meV_fu=strain,
        fixed_strain_pattern_contrast_meV_fu=pattern,interaction_meV_fu=interaction,
        unique_mode_energy_partition=False,statistical_or_certified_barrier_error=False)


def analyze():
    prep=load(CASE/'prepare.py','E277_preparation')
    study=load(BASE/'matched_pattern_release_E270_20261018/study.py','E277_analysis_chart')
    reader=load(BASE/'preserving_pattern_branches_E272_20261018/analyze_prefix.py','E277_analysis_native')
    cached=load(BASE/'preserving_local_completion_E276_20261019/worker.py','E277_cache_native')
    collector=load(CASE/'collect_terminal.py','E277_collection')
    cfg=json.loads((CASE/'prepared/config.json').read_bytes())
    receipt=json.loads((CASE/'completed/receipt.json').read_bytes());archive=CASE/'completed/native_terminal.tar.gz'
    if (sha256(archive)!=receipt['archive_sha256'] or receipt['bundle_pin']!=collector.PIN
            or receipt['collector_sha256']!=sha256(CASE/'collect_terminal.py')
            or receipt['new_DFT_calls'] or receipt['scheduler_mutations']):raise ValueError('terminal archive identity differs')
    context_path=BASE/'preserving_strain_E275_20261019/completed_analysis_full/analysis.json'
    if sha256(context_path)!='829e524e083b9622ddb41da9a0029b0246deeaa05e97b88733f59aef1738ed7f':
        raise ValueError('same-strain initial wells must retain their original source')
    context=json.loads(context_path.read_bytes())
    wells={s:next(b for b in context['bands'] if b['strain']==s and b['kind']=='registered_preserving')
        ['endpoint_native_records'][0]['energy_eV_cell'] for s in (0.,.01)}
    t=prep.read(BASE/'reference_variants/T.vasp')
    def result_row(label,strain,pattern,chart,x,e,g,source,new_calls):
        seed=CASE/'prepared'/('zero_q1' if strain==0. else 'plus1_q0')/'seed'
        data=json.loads((seed/'seed.json').read_bytes());q=np.asarray(data['orthonormal_39_by_4'])
        tr=chart.rigid_translation_directions()/np.sqrt(12)
        residual=float(np.linalg.norm(g-q@(q.T@g)-tr@(tr.T@g)))
        if residual>.03:raise ValueError('unqualified local candidate')
        held,polar,cell=prep.observables(chart.to_atoms(x))
        if not np.allclose(held,cfg['diagonals'][pattern]['held_values_A'],atol=1e-10,rtol=0):
            raise ValueError('crossed four-vector must equal the declared donor values')
        deformation=np.linalg.solve(t.cell.array,cell).T
        return dict(condition=label,strain=strain,held_pattern=pattern,origin=source,
            energy_eV_cell=e,initial_energy_eV_cell=wells[strain],relative_energy_meV_fu=250*(e-wells[strain]),
            released_gradient_eV_A=residual,Qx_A=float(held[0]),Qz_A=float(held[1]),Qy_A=float(held[2]),
            polar_proxy_z_A=float(held[3]),transverse_proxy_A=float(polar[1]),
            open_y_stretch_percent=100*(deformation[1,1]-1),coordinate_norm_A=float(np.linalg.norm(x)),
            new_SCF_calls=new_calls)
    rows=[];calls=[]
    for pattern,name in (('q0','zero_q1'),('q1','plus1_q0')):
        chart=study.chart_for(CASE/'prepared'/name/'seed');d=cfg['diagonals'][pattern]
        x,e,g=cached.cached_start(chart,CASE/'prepared/diagonals'/pattern/'native',{})
        rows.append(result_row('diagonal_'+pattern,d['strain'],pattern,chart,x,e,g,'reused_native_J',0))
    with tempfile.TemporaryDirectory(prefix='E277-results-') as td:
        dest=Path(td)
        with tarfile.open(archive) as tar:
            for m in tar.getmembers():
                p=Path(m.name)
                if (p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0]!='native_terminal'
                        or not(m.isfile() or m.isdir()) or not(dest/p).resolve().is_relative_to(dest.resolve())):
                    raise ValueError('unsafe native member')
            tar.extractall(dest)
        native=dest/'native_terminal';files=json.loads((native/'files.json').read_bytes())
        if len(files)!=receipt['files'] or any(sha256(native/n)!=h for n,h in files.items()):
            raise ValueError('native files changed')
        for record in receipt['rows']:
            name=record['condition'];folder=native/name;condition=json.loads((CASE/'prepared'/name/'condition.json').read_bytes())
            complete=json.loads((folder/'complete.json').read_bytes());chart=study.chart_for(CASE/'prepared'/name/'seed')
            if (complete['job_id']!=record['job_id'] or complete['bundle_pin']!=collector.PIN
                    or complete['policy']!=prep.POLICY or record['scheduler'][:6]!=
                    [record['job_id'],'COMPLETED','0:0','hfacnormal01','1','32']):raise ValueError('original terminal outcome required')
            c=np.asarray(json.loads((CASE/'prepared'/name/'seed/seed.json').read_bytes())['covectors_39_by_4'])
            values={};seconds=0.
            for p in sorted(folder.glob('evaluations/eval-*/result.json')):
                raw=json.loads(p.read_bytes());x=np.asarray(raw['coordinates'])
                _,e,g,evidence=reader.replay_native(folder,chart,raw['coordinates'])
                if not np.allclose(c.T@x,condition['constraint_offset_A'],atol=1e-10,rtol=0) or np.linalg.norm(x)>1.:
                    raise ValueError('a native evaluation violated the predeclared physical constraints/guard')
                values[x.tobytes()]=(e,g);seconds+=raw['validation']['elapsed_seconds']
                calls.append(dict(condition=name,result_sha256=evidence['native_result_sha256'],energy_eV_cell=e))
            if len(values)!=complete['new_calls'] or not 0<len(values)<=40:raise ValueError('finite native call accounting differs')
            candidate=complete['candidate']
            if candidate is None:
                rows.append(dict(condition=name,strain=condition['substrate_strain'],held_pattern=condition['held_pattern'],
                    relative_energy_meV_fu=None,status=complete['status'],new_SCF_calls=len(values)))
            else:
                x=np.asarray(candidate['coordinates']);e,g=values[x.tobytes()]
                if e!=candidate['energy_eV_cell']:raise ValueError('candidate energy differs from native')
                row=result_row(name,condition['substrate_strain'],condition['held_pattern'],chart,x,e,g,'new_cross_J',len(values))
                if not np.isclose(row['released_gradient_eV_A'],candidate['released_gradient_eV_A'],atol=2e-12,rtol=0):
                    raise ValueError('released gradient differs from native')
                rows.append(row)
            record['native_proxy_core_hours']=seconds*32/3600
    rows.sort(key=lambda r:(r['strain'],r['held_pattern']))
    path_context=[{k:b[k] for k in ('strain','kind','sampled_polygon_peak_meV_fu','fmax_eV_A')}
                  for b in context['bands']]
    for row in rows:
        peak=next(b['sampled_polygon_peak_meV_fu'] for b in path_context
                  if b['strain']==row['strain'] and b['kind']=='released_seed')
        row['connected_free_sampled_peak_meV_fu']=peak
        row['local_candidate_minus_free_sampled_peak_meV_fu']=(None
            if row['relative_energy_meV_fu'] is None else row['relative_energy_meV_fu']-peak)
    return dict(rows=rows,contrasts=contrasts(rows),native_calls=calls,new_native_SCF_calls=len(calls),
        whole_path_context=path_context,local_minus_path_is_not_a_unique_mode_energy=True,
        costs=receipt['rows'],archive_sha256=receipt['archive_sha256'],
        source_code_sha256=sha256(Path(__file__)),prepared_config_sha256=sha256(CASE/'prepared/config.json'),
        whole_path_context_sha256=sha256(context_path),new_DFT_calls=0,physical_inputs_changed=False,
        full_stability_TS_global_barrier_or_new_prediction_claimed=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();d=analyze();a.output.mkdir(parents=True,exist_ok=False)
    with (a.output/'analysis.json').open('x',encoding='utf-8',newline='\n') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
    keys=list(dict.fromkeys(k for r in d['rows'] for k in r))
    with (a.output/'source_data.csv').open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(d['rows'])
    print(json.dumps(dict(rows=d['rows'],contrasts=d['contrasts'],new_calls=d['new_native_SCF_calls'])))
