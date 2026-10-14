"""Read-only fixed-Q/fixed-cell B3 native prefix and exact optimizer replay.

Use original E161 PYTHONPATH and OMP/MKL/OpenBLAS=1. Never submit, resume,
repair a running source, launch a calculator or promote an incomplete trial.
"""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import argparse,hashlib,importlib.util,json,os,subprocess,sys
sys.dont_write_bytecode=True
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator

EXEC=Path('/public/home/iai806/abacus/agent-runs/20261015-varneb-G3-atomic-release-E175-r2')
MANIFEST_SHA='99e462af901e05e3665f3d034973fae478e391cd2f89740bb3a8851358a50096'
CONFIG_SHA='5bd401c6ca4044ac654f2b19a22371b20c17d3736d207ff7044336d9ed463f38'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def load_prefix(directory,config):
    rows=[];pending=None
    for channel in config['policy']['channels']:
        for p in config['channels'][channel]['points']:
            path=directory/channel/f"point_{p['id']:03d}.json"
            if not path.exists():return rows,dict(channel=channel,id=p['id'],reason='point_not_present')
            try:row=json.loads(path.read_bytes())
            except (json.JSONDecodeError,UnicodeDecodeError):
                return rows,dict(channel=channel,id=p['id'],reason='point_JSON_incomplete')
            if any(json.dumps(row.get(k),allow_nan=False)!=json.dumps(v,allow_nan=False)
                   for k,v in dict(channel=channel,id=p['id'],q_A=p['q_A']).items()):
                raise ValueError('completed prefix channel/id/Q differs')
            rows.append(row)
    return rows,pending


class CacheOnly:
    def __init__(self,records):self.records=records;self.seen={}
    def __call__(self,x):
        key=hashlib.sha256(np.asarray(x).tobytes()).hexdigest()
        if key not in self.seen:
            i=len(self.seen)
            if i==12:raise RuntimeError('mode evaluation attempt budget exhausted; calculator was not called')
            if i>=len(self.records) or not np.array_equal(x,self.records[i]['coordinates']):
                raise ValueError('missing or differently ordered replay trial; no DFT')
            self.seen[key]=self.records[i]
        p=self.seen[key];return p['enthalpy_eV'],np.array(p['gradient'],copy=True)


def replay_node(w,data,point,row,records,r,controls):
    if (type(row['new_completed_SCFs']) is not int or type(row['reserved_attempts']) is not int
            or row['new_completed_SCFs']!=len(records) or row['reserved_attempts']!=len(records)
            or len(records)>12):raise ValueError('completed-node native attempt accounting differs')
    delegate=CacheOnly(records)
    evaluate=w.CachedPullback(data,controls,delegate,point)
    candidate=None;reason=None
    try:candidate=w.relax_point(data,np.asarray(point['q_A']),r,evaluate,base=point['frozen_internal_A'])
    except RuntimeError as exc:
        reason=str(exc)
        if reason not in ('mode evaluation attempt budget exhausted; calculator was not called',
                          'atomic released-gradient target not reached'):raise
    status='atomic_gradient_candidate_stability_unchecked' if candidate else 'optimizer_not_converged'
    if (candidate!=row['candidate'] or reason!=row['reason'] or row['status']!=status
            or len(delegate.seen)!=len(records) or evaluate.hits!=row['cached_callbacks']
            or evaluate.requests!=row['all_callbacks']
            or row['candidate_native_identity']!=(evaluate.last if candidate else None)):
        raise ValueError('strict cached replay candidate/counters/source identity differs')
    return dict(exact_request_order_optimizer_replay=True,new_DFT=0,
        cached_callbacks=evaluate.hits,all_callbacks=evaluate.requests)


def native_records(w,data,point,directory):
    """All reservations/contract plus raw EFS and exact literal seed cell."""
    from vcneb.mode_evaluator import _fingerprint
    chart=data['chart'];ref=chart.reference_atoms
    expected=dict(schema_version=1,coordinate_convention=chart.coordinate_convention,
        reference=dict(numbers=ref.numbers.tolist(),cell_A=ref.cell.array.tolist(),
            scaled_positions=ref.get_scaled_positions(wrap=False).tolist(),pbc=ref.pbc.tolist()),
        calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
        validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',pressure_eV_per_A3=0.,
        minimum_distance_A=1.6,mechanical_boundary=chart.mechanical_contract)
    fingerprint=_fingerprint(expected)
    budget=dict(schema_version=1,max_evaluation_attempts=12,physical_contract_sha256=fingerprint)
    if (json.loads((directory/'contract.json').read_bytes())!=dict(expected,contract_sha256=fingerprint)
            or json.loads((directory/'evaluation-budget.json').read_bytes())!=budget):
        raise ValueError('exact original physical evaluator contract/budget differs')
    paths=sorted(directory.glob('eval-*/result.json'));markers=sorted(directory.glob('attempt-*.json'))
    if (list(directory.glob('eval-*/failure.json')) or len(list(directory.glob('eval-*')))!=len(paths)
            or len(markers)!=len(paths) or len(paths)>12):raise ValueError('completed node has failed/incomplete/unreserved trials')
    records=[];sources=[];seed=chart.to_atoms(point['frozen_full_A'])
    for i,(p,marker) in enumerate(zip(paths,markers)):
        row=json.loads(p.read_bytes());x=np.asarray(row['coordinates']);digest=hashlib.sha256(x.tobytes()).hexdigest()
        reservation=dict(schema_version=1,evaluation_index=i,coordinate_sha256=digest,
            coordinates=x.tolist(),budget_sha256=_fingerprint(budget))
        if (p.parent.name!=f'eval-{i:06d}-{digest[:12]}' or marker.name!=f'attempt-{i:06d}.json'
                or json.loads(marker.read_bytes())!=reservation or row['status']!='complete'
                or row['coordinate_sha256']!=digest or row['contract_sha256']!=fingerprint or x[-1]!=0):
            raise ValueError('new native coordinate/contract/reservation identity differs')
        raw=Path(row['validation']['raw_directory']);expected_raw=p.parent/'scf_000000'
        if (raw.name!='scf_000000' or raw.resolve(strict=True)!=expected_raw.resolve(strict=True)
                or not raw.resolve(strict=True).is_relative_to((EXEC/'atomic_release').resolve(strict=True))):
            raise ValueError('recorded raw directory not exactly expected native calculation')
        atoms=chart.to_atoms(x);atoms.calc=SinglePointCalculator(atoms,**w.study.audited_results(expected_raw))
        checked=w.study.result_validator(i,atoms,raw.parent);e,g=chart.energy_and_gradient(x,atoms)
        if (checked!=row['validation'] or e!=row['enthalpy_eV']
                or not np.allclose(g,row['gradient'],atol=2e-12,rtol=0)
                or not np.array_equal(atoms.get_stress(voigt=False),row['stress_eV_per_A3'])
                or not np.array_equal(atoms.cell.array,seed.cell.array)
                or not np.allclose(data['axes'].T@(data['basis'].T@x),point['q_A'],atol=2e-12,rtol=0)):
            raise ValueError('native EFS/full gradient/exact seed cell/fixed Q differs')
        records.append(row);sources.append(dict(result_sha256=sha(p),validation=checked))
    return records,sources


def audit(mode):
    if any(os.environ.get(k)!='1' for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS')):
        raise ValueError('original single-thread numerical replay environment required')
    bundle=EXEC/'bundle'
    if sha(bundle/'bundle_files.json')!=MANIFEST_SHA:raise ValueError('immutable five-file manifest differs')
    files=json.loads((bundle/'bundle_files.json').read_bytes())
    if len(files)!=5 or any(sha(bundle/n)!=h for n,h in files.items()):raise ValueError('immutable bundle drift')
    if sha(bundle/'config.json')!=CONFIG_SHA:raise ValueError('original E175config required')
    if json.loads((EXEC/'accepted.json').read_bytes())['job_id']!='28867177':raise ValueError('original accepted handle differs')
    config=json.loads((bundle/'config.json').read_bytes())
    account=subprocess.run(['sacct','-n','-P','-X','-j','28867177',
        '--format=JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,NodeList,Partition'],
        check=True,capture_output=True,text=True).stdout.strip()
    f=account.split('|')
    if len(f)<7 or f[0]!='28867177' or f[6]!='hfacnormal01':raise RuntimeError('original accounting identity unavailable')
    rows,pending=load_prefix(EXEC/'atomic_release',config)
    if mode=='complete' and (f[1:3]!=['COMPLETED','0:0'] or len(rows)!=10):
        raise RuntimeError('complete audit requires original normal completion/all10points')
    audited=[];sources=[];old_proof=None
    if rows:
        w=module('unchanged_E175_atomic_worker',bundle/'worker.py');data,old_proof=w.preflight(config,bundle)
        controls=w.load_controls(bundle)['rows']
        for row in rows:
            channel=row['channel'];d=data[channel];p=config['channels'][channel]['points'][row['id']]
            r=np.asarray(config['channels'][channel]['atomic_release_basis_36_by_31'])
            records,native=native_records(w,d,p,EXEC/'atomic_release'/channel/f"node_{row['id']:03d}")
            replay=replay_node(w,d,p,row,records,r,[p for p in controls if p['channel']==channel])
            point=None
            if row['candidate'] is not None:
                final=row['candidate_native_identity'];g=d['basis'].T@final['physical_gradient_eV_A']
                norm=float(np.linalg.norm(r.T@g))
                if norm>.03 or not np.isclose(norm,row['candidate']['atomic_released_gradient_norm_eV_A'],atol=2e-12,rtol=0):
                    raise ValueError('candidate full atomic-released gradient differs')
                point=dict(energy_eV_cell=final['energy_eV_cell'],atomic_released_gradient_norm_eV_A=norm,
                    joint_released_gradient_norm_eV_A=float(np.linalg.norm(g-d['axes']@(d['axes'].T@g))),
                    full_physical_gradient_eV_A=final['physical_gradient_eV_A'],
                    full_coordinate=final['coordinates_full_A'],released_stability_unchecked=True)
            audited.append(dict(channel=channel,id=row['id'],q_A=row['q_A'],candidate=point,
                reason=row['reason'],new_completed_SCFs=len(records),
                sum_completed_call_seconds=sum(r['validation']['elapsed_seconds'] for r in records),replay=replay))
            sources.extend(native)
    if mode=='complete':
        summary=json.loads((EXEC/'atomic_release/atomic_release.json').read_bytes())
        if (summary['rows']!=rows or summary['new_completed_SCFs']!=sum(r['new_completed_SCFs'] for r in rows)
                or summary['reserved_attempts']!=summary['new_completed_SCFs'] or summary['job_id']!='28867177'):
            raise ValueError('terminal summary/source accounting differs')
    return dict(status='read_only_atomic_treatment_native_prefix_not_stability_or_independent_prediction',
        observed_CST=datetime.now(timezone(timedelta(hours=8))).isoformat(),audit_mode=mode,job_id='28867177',sacct=account,
        completed_nodes_audited=len(audited),gradient_candidates=sum(p['candidate'] is not None for p in audited),
        native_new_SCFs_reparsed=len(sources),candidate_rows=audited,native_sources=sources,pending=pending,
        old_native_cache_preflight=old_proof,all10points_covered=len(rows)==10,
        new_DFT_or_scheduler_writes=0,original_runtime_or_physics_changed=False,
        physical_reserved_geometry_generated_or_read=False,
        conditional_stability_TS_acceleration_or_independent_prediction_certified=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mode',choices=('prefix','complete'),required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('fresh external audit report required')
    result=audit(a.mode)
    with a.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(audited_nodes=result['completed_nodes_audited'],new_native_SCFs=result['native_new_SCFs_reparsed'],
        pending=result['pending'],new_DFT=0)))
