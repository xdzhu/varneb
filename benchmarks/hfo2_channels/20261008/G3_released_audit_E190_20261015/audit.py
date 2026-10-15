# E190: preserve E183; fix lexical provenance replay after exact resolved-path checks.
"""Read-only E181 native stencils: prefix, complete, and one old-record smoke.

Use E161 PYTHONPATH and single-thread BLAS. Never construct a calculator,
submit, retry, repair live files, or turn an incomplete block into a result.
"""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import argparse,hashlib,importlib.util,json,os,subprocess,sys
sys.dont_write_bytecode=True
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from vcneb.mode_evaluator import _fingerprint

EXEC=Path('/public/home/iai806/abacus/agent-runs/20261015-varneb-G3-release-curvature-E181-r1')
MANIFEST_SHA='2e61cfc1ee403e00b91363453fdb4b485dcdeb2eb492ea184eba6f4dd100283d'
CONFIG_SHA='71a2e1b8a46b2ef903e4023d329109c5e74975aaf5bab54f4d4638d15c1add39'
ACCEPTED_SHA='4e49130bbb897b38062fc3d9122e6a5200972195b87d549465a8e5eb1ec73643'
HELD_SHA='4eaac24691c8794f0e936c948a068eb23d5a26512ca1d8f407960c14643e7312'
PREFLIGHT_SHA='9630d04159a7213453d48bb6831235a38ec15d48e9b4948884bc102cb1ec3d16'
JOB='28868827'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def pinned(path,digest):
    if sha(path)!=digest:raise ValueError('immutable source bytes differ')
    return json.loads(path.read_bytes())


def load_worker():
    bundle=EXEC/'bundle';files=pinned(bundle/'bundle_files.json',MANIFEST_SHA)
    if set(files)!={'worker.py','run.slurm','config.json','protocol.md'} or any(sha(bundle/n)!=h for n,h in files.items()):
        raise ValueError('original four-file bundle drift')
    config=pinned(bundle/'config.json',CONFIG_SHA)
    if pinned(EXEC/'accepted.json',ACCEPTED_SHA)['job_id']!=JOB:raise ValueError('accepted handle differs')
    pinned(EXEC/'held_validation.json',HELD_SHA);pinned(EXEC/'native_preflight.json',PREFLIGHT_SHA)
    s=importlib.util.spec_from_file_location('original_E181_frozen_worker',bundle/'worker.py')
    w=importlib.util.module_from_spec(s);s.loader.exec_module(w)
    data=w.validate(config,bundle)
    if w.study.ROOT.resolve()!=w.OLD.resolve():raise ValueError('original E161 runtime required')
    sources=pinned(w.OLD/'runtime_files.json',w.RUNTIME_SHA)
    if len(sources)!=380 or any(sha(w.OLD/n)!=h for n,h in sources.items()):raise ValueError('original380runtime drift')
    if sha(w.study.BINARY)!=w.study.BINARY_SHA or any(sha(w.study.ASSETS/n)!=h for n,h in w.study.CONTRACT.items()):
        raise ValueError('original executable or six assets drift')
    return w,config,data


def physical_contract(chart):
    ref=chart.reference_atoms
    c=dict(schema_version=1,coordinate_convention=chart.coordinate_convention,
        reference=dict(numbers=ref.numbers.tolist(),cell_A=ref.cell.array.tolist(),
            scaled_positions=ref.get_scaled_positions(wrap=False).tolist(),pbc=ref.pbc.tolist()),
        calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
        validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',pressure_eV_per_A3=0.,
        minimum_distance_A=1.6,mechanical_boundary=chart.mechanical_contract)
    return dict(c,contract_sha256=_fingerprint(c))


def accounting(raw):
    f=[v.strip() for v in raw.strip().split('|')]
    if len(f)!=8 or f[0]!=JOB or f[6]!='hfacnormal01':raise ValueError('original accounting identity unavailable')
    if not f[3].isdigit() or not f[4].isdigit() or not f[7].isdigit():raise ValueError('invalid allocation accounting')
    cpus,nodes=int(f[4]),int(f[7])
    if cpus not in (0,32) or nodes not in (0,1):raise ValueError('unexpected actual allocation layout')
    normal=f[1:3]==['COMPLETED','0:0'] and cpus==32 and nodes==1
    return dict(raw=raw.strip(),state=f[1],exit=f[2],elapsed_seconds=int(f[3]),CPUs=cpus,
        node_list=f[5],nodes=nodes,normal_completion=normal,
        allocation_core_hours=int(f[3])*cpus/3600)


def native_record(w,data,x,path,namespace):
    """Reparse one complete native E/F/stress record, never evaluate it."""
    row=json.loads(path.read_bytes());chart=data['chart'];x=np.asarray(x)
    digest=hashlib.sha256(x.tobytes()).hexdigest();contract=physical_contract(chart)
    if (row.get('status')!='complete' or row['coordinate_sha256']!=digest
            or row['contract_sha256']!=contract['contract_sha256']
            or not np.array_equal(row['coordinates'],x)
            or not path.resolve(strict=True).is_relative_to(namespace.resolve(strict=True))):
        raise ValueError('complete record coordinate/contract/namespace identity differs')
    expected=path.parent/'scf_000000';raw=Path(row['validation']['raw_directory'])
    if raw.name!='scf_000000' or raw.resolve(strict=True)!=expected.resolve(strict=True):
        raise ValueError('recorded raw directory differs from its native calculation')
    atoms=chart.to_atoms(x);atoms.calc=SinglePointCalculator(atoms,**w.study.audited_results(expected))
    # Resolved identity was checked above; replay the stored lexical source.
    # Do not mistake bundle/../ provenance spelling for a different SCF.
    checked=w.study.result_validator(0,atoms,raw.parent);e,g=chart.energy_and_gradient(x,atoms)
    if json.loads((expected/'input_sha256.json').read_bytes())!=checked['input_sha256']:
        raise ValueError('native initial input manifest differs')
    distances=atoms.get_all_distances(mic=True);np.fill_diagonal(distances,np.inf)
    if (checked!=row['validation'] or e!=row['enthalpy_eV']
            or not np.allclose(g,row['gradient'],atol=2e-12,rtol=0)
            or not np.array_equal(atoms.get_stress(voigt=False),row['stress_eV_per_A3'])
            or not np.isclose(atoms.get_volume(),row['volume_A3'],atol=2e-12,rtol=0)
            or not np.isclose(distances.min(),row['minimum_distance_A'],atol=2e-12,rtol=0)
            or not np.isclose(np.max(np.linalg.norm(atoms.get_forces(),axis=1)),row['maximum_atomic_force_eV_per_A'],atol=2e-12,rtol=0)):
        raise ValueError('native energy/gradient/stress/volume/force/distance differs')
    return row,dict(result_path=str(path),result_sha256=sha(path),validation=checked)


def reserve_identity(row,probe,index,budget):
    x=np.asarray(probe['coordinates_full_A']);digest=hashlib.sha256(x.tobytes()).hexdigest()
    expected=dict(schema_version=1,evaluation_index=index,coordinate_sha256=digest,
        coordinates=x.tolist(),budget_sha256=_fingerprint(budget))
    if json.dumps(row,sort_keys=True,allow_nan=False)!=json.dumps(expected,sort_keys=True,allow_nan=False):
        raise ValueError('exact ordered stencil reservation differs')
    return f'eval-{index:06d}-{digest[:12]}'


def reconstruct(spec,records):
    if len(records)!=136:raise ValueError('all136records and both amplitudes required')
    directions=np.asarray(spec['literal_directions_40_by_34']);g=np.asarray([r['gradient'] for r in records])
    if g.shape!=(136,40) or not np.isfinite(g).all():raise ValueError('complete finite full gradients required')
    blocks=[]
    for j,h in enumerate((.01,.02)):
        part=g[j*68:(j+1)*68];action=(part[::2]-part[1::2]).T/(2*h)
        raw=directions.T@action;symmetric=(raw+raw.T)/2
        eig=np.linalg.eigvalsh(symmetric)
        blocks.append(dict(step_A=h,complete_pairs=34,full_hessian_action_eV_A2=action.tolist(),
            raw_projected_eV_A2=raw.tolist(),symmetric_projected_eigenvalues_eV_A2=eig.tolist(),
            raw_negative_count=int(np.count_nonzero(eig<0)),
            reciprocity_relative_defect=float(np.linalg.norm(raw-raw.T)/max(np.linalg.norm(raw),1e-30)),
            symmetrization_operator_change_eV_A2=float(np.linalg.norm(raw-symmetric,ord=2))))
    return blocks


def check_summary(summary,channel,spec,records,blocks):
    expected=dict(status='complete_post_release_raw_transverse_curvature_not_stability_certificate',
        channel=channel,new_native_SCFs=136,new_center_SCFs=0,
        base_released_gradient_norm_eV_A=spec['base_released_gradient_norm_eV_A'],
        physical_stability_TS_or_independent_prediction_certified=False)
    if any(type(summary.get(k)) is not type(v) or summary.get(k)!=v for k,v in expected.items()):
        raise ValueError('channel summary scope/count/source differs')
    arrays=dict(energies_eV_cell=[r['enthalpy_eV'] for r in records],
        full_gradients_eV_A=[r['gradient'] for r in records],
        sum_native_call_seconds=sum(r['validation']['elapsed_seconds'] for r in records))
    arrays['call_time_proxy_core_hours']=arrays['sum_native_call_seconds']*32/3600
    arrays['two_step_raw_block_operator_difference_eV_A2']=float(np.linalg.norm(
        np.asarray(blocks[0]['raw_projected_eV_A2'])-blocks[1]['raw_projected_eV_A2'],ord=2))
    for k,v in arrays.items():
        if not np.allclose(summary[k],v,atol=2e-12,rtol=0):raise ValueError('summary differs from native records')
    if len(summary['blocks'])!=2:raise ValueError('both original amplitudes required')
    for old,new in zip(summary['blocks'],blocks):
        for k in old:
            if k not in new or not np.allclose(old[k],new[k],atol=2e-12,rtol=0):
                raise ValueError('stored block differs from independently reconstructed native gradients')


def channel_prefix(w,data,spec,directory,*,complete):
    if not directory.exists():
        if complete:raise ValueError('complete channel has no evaluation directory')
        return dict(records=[],sources=[],reserved=0,pending='not_started',failures=[],blocks=None)
    contract=physical_contract(data['chart']);budget=dict(schema_version=1,max_evaluation_attempts=136,
        physical_contract_sha256=contract['contract_sha256'])
    if not (directory/'contract.json').exists() or not (directory/'evaluation-budget.json').exists():
        if complete:raise ValueError('complete channel missing evaluator contract')
        return dict(records=[],sources=[],reserved=0,pending='contract_not_complete',failures=[],blocks=None)
    if json.loads((directory/'contract.json').read_bytes())!=contract or json.loads((directory/'evaluation-budget.json').read_bytes())!=budget:
        raise ValueError('physical evaluator contract or136attempt budget differs')
    # Freeze one inventory; later append-only outputs belong to a later audit.
    markers=sorted(directory.glob('attempt-*.json'));paths=sorted(directory.glob('eval-*'))
    if len(markers)>136 or len(paths)>len(markers):raise ValueError('unreserved or over-budget calculation')
    records=[];sources=[];pending=None;failures=[]
    for i,marker in enumerate(markers):
        if marker.name!=f'attempt-{i:06d}.json':raise ValueError('missing or differently ordered reservation')
        try:reservation=json.loads(marker.read_bytes())
        except (json.JSONDecodeError,UnicodeDecodeError):
            if complete:raise ValueError('incomplete terminal reservation')
            pending='reservation_JSON_inflight';break
        name=reserve_identity(reservation,spec['probes'][i],i,budget);parent=directory/name
        if i<len(paths) and paths[i].name!=name:raise ValueError('directory coordinate identity differs')
        if (parent/'failure.json').exists():
            failures.append(dict(index=i,path=str(parent/'failure.json'),sha256=sha(parent/'failure.json')));break
        result=parent/'result.json'
        if not result.exists():pending=f'reserved_probe_{i}_output_pending';break
        w.safe_probe(data,np.asarray(spec['probes'][i]['coordinates_full_A']))
        row,source=native_record(w,data,spec['probes'][i]['coordinates_full_A'],result,directory)
        records.append(row);sources.append(source)
    if complete and (pending is not None or failures or len(records)!=136 or len(paths)!=136 or len(markers)!=136):
        raise ValueError('complete channel has incomplete, failed or missing native probes')
    blocks=None if len(records)!=136 else reconstruct(spec,records)
    summary=directory.parent/'curvature.json'
    if complete and not summary.exists():raise ValueError('terminal channel summary missing')
    if blocks is not None and summary.exists():check_summary(json.loads(summary.read_bytes()),directory.parent.name,spec,records,blocks)
    return dict(records=records,sources=sources,reserved=len(markers),pending=pending,failures=failures,blocks=blocks)


def audit(mode):
    if any(os.environ.get(n)!='1' for n in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS')):
        raise ValueError('original single-thread audit environment required')
    w,config,data=load_worker()
    if mode=='smoke':
        d=data['T'];record=d['origin_record'];raw=Path(record['validation']['raw_directory'])
        row,source=native_record(w,d,record['coordinates'],raw.parent/'result.json',raw.parent.parent)
        if row!=record:raise ValueError('old native smoke record differs from pinned released base')
        return dict(status='one_existing_native_record_smoke_not_new_stencil_data',native_old_SCFs_reparsed=1,
            native_new_stencil_SCFs_reparsed=0,source=source,new_DFT_or_scheduler_writes=0,
            stable_branch_TS_or_independent_prediction_certified=False)
    raw=subprocess.run(['sacct','-n','-P','-X','-j',JOB,
        '--format=JobIDRaw,State%30,ExitCode,ElapsedRaw,AllocCPUS,NodeList,Partition,NNodes'],
        check=True,capture_output=True,text=True).stdout
    account=accounting(raw)
    if mode=='complete' and not account['normal_completion']:raise RuntimeError('complete audit requires normal Slurm completion')
    rows=[];sources=[]
    for c in config['policy']['channels']:
        r=channel_prefix(w,data[c],config['channels'][c],EXEC/'released_curvature'/c/'evaluations',complete=mode=='complete')
        if mode=='complete' and (len(r['records'])!=136 or r['reserved']!=136 or r['blocks'] is None):
            raise ValueError('both channels must supply all136 native records before completion')
        rows.append(dict(channel=c,completed_native_records=len(r['records']),reserved_attempts=r['reserved'],
            pending=r['pending'],failures=r['failures'],blocks=r['blocks'],
            completed_call_proxy_core_hours=sum(v['validation']['elapsed_seconds'] for v in r['records'])*32/3600))
        sources.extend(r['sources'])
    if mode=='complete':
        terminal=json.loads((EXEC/'released_curvature/completed.json').read_bytes())
        expected=dict(status='complete_both_released_stencils',new_native_SCFs=272,new_center_SCFs=0,
            job_id=JOB,allocation_core_hours=None,all_other_training_and_preparation_costs_additional=True,
            channels=[dict(channel=r['channel'],new_native_SCFs=136,call_time_proxy_core_hours=r['completed_call_proxy_core_hours']) for r in rows])
        if terminal!=expected:raise ValueError('whole-job completion/cost/source summary differs')
    return dict(status='complete_native_released_stencils_not_stability_certificate' if mode=='complete' else 'native_released_stencil_prefix_not_completion',
        mode=mode,observed_CST=datetime.now(timezone(timedelta(hours=8))).isoformat(),job_id=JOB,accounting=account,
        config_sha256=CONFIG_SHA,bundle_manifest_sha256=MANIFEST_SHA,native_new_stencil_SCFs_reparsed=len(sources),
        channels=rows,native_sources=sources,inventory_is_start_snapshot_new_appends_require_later_audit=True,
        all_other_curvature_training_preparation_costs_additional=True,new_DFT_or_scheduler_writes=0,
        physical_reserved_geometry_generated_or_read=False,stable_branch_TS_or_independent_prediction_certified=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mode',choices=('prefix','complete','smoke'),required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('fresh external audit output required')
    if a.output.resolve().is_relative_to(EXEC.resolve()):raise ValueError('report must remain outside the live producer namespace')
    result=audit(a.mode)
    with a.output.open('x',encoding='utf-8',newline='\n') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k in ('status','native_old_SCFs_reparsed','native_new_stencil_SCFs_reparsed','new_DFT_or_scheduler_writes')}))
