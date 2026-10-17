"""Native first-four training prefix and cached optimizer replay; zero DFT."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
import numpy as np

EXEC=Path('/public/home/iai806/abacus/agent-runs/20261017-varneb-G3-P-branches-E260-r2')
CONFIG_SHA='63b1aa0eea88e7cd7c2f1a78ce9f4f3c7f39e16a9caf04d3f6d8c6ef5f84ab80'
MANIFEST_SHA='4422902d6474afe1056181b978f6098c9fc1284ff6cc66d11ff865d499d87eff'
MESSAGE='mode evaluation attempt budget exhausted; calculator was not called'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def replay(solver,point,basis,free,direction,seed,records,guard):
    """Consume native records in order; never select a nearby cached point."""
    index=0;maximum_difference=0.;captured={}
    start=np.asarray(point['starting_internal_coordinates_A'])
    def evaluate(x):
        nonlocal index,maximum_difference
        if index==len(records):raise RuntimeError(MESSAGE)
        record=records[index];expected=np.asarray(record['coordinates']);actual=basis@x
        error=float(np.max(np.abs(actual-expected)))
        if not np.allclose(actual,expected,atol=2e-12,rtol=0):
            raise ValueError('ordered recorded proposal differs; no nearest-cache substitution')
        maximum_difference=max(maximum_difference,error);index+=1
        return record['enthalpy_eV'],basis.T@np.asarray(record['gradient'])
    code=solver.relax_branch_candidate.__code__
    def trace(frame,event,arg):
        if frame.f_code is code and event=='exception' and str(arg[1])==MESSAGE:
            state=frame.f_locals
            captured.update(accepted_trace=[dict(r) for r in state['trace']],
                last_accepted_internal_A=(state['anchor']+state['basis']@state['z']).tolist(),
                last_accepted_energy_eV_cell=float(state['energy']),
                last_accepted_gradient_internal_eV_A=np.asarray(state['full_gradient']).tolist(),
                trust_radius_A=float(state['radius']),
                uncomputed_next_proposal_internal_A=(state['anchor']+state['basis']@state['proposal']).tolist())
        return trace
    previous=sys.gettrace();sys.settrace(trace)
    try:
        result=solver.relax_branch_candidate(evaluate,start,free,direction,
            initial_hessian=seed,valid=guard,gradient_tolerance=.03,
            branch_gradient_tolerance=1e-4,initial_trust_radius=.05,max_iterations=200)
        outcome=dict(status='normal_return_replayed',result=result)
    except RuntimeError as exc:
        if str(exc)!=MESSAGE:raise
        outcome=dict(status='budget_exception_replayed',**captured)
    finally:sys.settrace(previous)
    if index!=len(records):raise ValueError('every recorded SCF must be consumed in replay')
    outcome.update(recorded_SCFs_consumed=index,maximum_coordinate_arithmetic_difference_A=maximum_difference,
        new_DFT=0,uncomputed_proposal_is_not_a_DFT_point=True)
    return outcome


def native_prefix():
    bundle=EXEC/'bundle'
    if sha(bundle/'config.json')!=CONFIG_SHA or sha(bundle/'bundle_files.json')!=MANIFEST_SHA:
        raise ValueError('original immutable configuration required')
    files=json.loads((bundle/'bundle_files.json').read_bytes())
    if any(sha(bundle/n)!=h for n,h in files.items()):raise ValueError('original bundle bytes differ')
    spec=importlib.util.spec_from_file_location('E262_original_training_worker',bundle/'worker.py')
    worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
    config=json.loads((bundle/'config.json').read_bytes())
    reader,w,data,_,_,free,seed,_=worker.native_context(config,bundle)
    direction=np.asarray(config['branch_direction_36']);basis=data['basis'];rows=[]
    for point in config['points'][:4]:
        shard=point['id']%2;directory=EXEC/f"shard_{shard}/point_{point['id']:03d}"
        producer_path=directory.parent/f"point_{point['id']:03d}.json"
        row=json.loads(producer_path.read_bytes())
        if row['id']!=point['id'] or row['q_A']!=point['q_A'] or row['start_sign']!=point['start_sign']:
            raise ValueError('complete first-four literal producer identity required')
        contract=reader.physical_contract(data['chart'])
        budget=dict(schema_version=1,max_evaluation_attempts=12,physical_contract_sha256=contract['contract_sha256'])
        if (json.loads((directory/'contract.json').read_bytes())!=contract
                or json.loads((directory/'evaluation-budget.json').read_bytes())!=budget
                or list(directory.glob('eval-*/failure.json'))):raise ValueError('original contract/budget or failure differs')
        markers=sorted(directory.glob('attempt-*.json'));records=[];sources=[]
        if not 1<=len(markers)<=12 or len(markers)!=row['reserved_attempts']:
            raise ValueError('all ordered attempted SCFs required')
        for i,marker in enumerate(markers):
            if marker.name!=f'attempt-{i:06d}.json':raise ValueError('reservation order differs')
            reservation=json.loads(marker.read_bytes());x=np.asarray(reservation['coordinates'])
            worker.geometry(data,x,point['q_A'])
            name=reader.reserve_identity(reservation,dict(coordinates_full_A=x.tolist()),i,budget)
            expected=directory/name/'result.json';literal=Path(row['native_sources'][i]['result_path'])
            if literal.resolve(strict=True)!=expected.resolve(strict=True):raise ValueError('exact native file required')
            record,source=reader.native_record(w,data,x,literal,directory)
            if source!=row['native_sources'][i]:raise ValueError('native source differs from producer')
            records.append(record);sources.append(source)
        if len(records)!=row['new_native_SCFs']:raise ValueError('all native calls required')
        def guard(x):
            try:worker.geometry(data,basis@x,point['q_A']);return True
            except ValueError:return False
        recovered=replay(worker.solver,point,basis,free,direction,seed,records,guard)
        if row['result'] is None:
            if row['status']!='attempt_budget_exhausted_no_retry' or recovered['status']!='budget_exception_replayed' or len(records)!=12:
                raise ValueError('original capped outcome did not replay')
        else:
            a,b=row['result'],recovered.get('result')
            if (b is None or a['status']!=b['status']
                    or not np.allclose(a['coordinates'],b['coordinates'],atol=2e-12,rtol=0)
                    or a['energy']!=b['energy']
                    or not np.isclose(a['released_gradient_norm'],b['released_gradient_norm'],atol=2e-12,rtol=0)
                    or not np.isclose(a['branch_gradient'],b['branch_gradient'],atol=2e-12,rtol=0)):
                raise ValueError('original returned candidate did not replay')
        rows.append(dict(id=point['id'],q_A=point['q_A'],start_sign=point['start_sign'],
            producer_status=row['status'],producer_sha256=sha(producer_path),native_sources=sources,
            native_records=records,replay=recovered,
            completed_call_proxy_core_hours=sum(r['validation']['elapsed_seconds'] for r in records)*32/3600))
    raw=subprocess.check_output(['sacct','-X','-n','-P','-j','28934048,28934049',
        '-o','JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,NodeList,Partition'],text=True)
    return dict(status='native_first_four_training_prefix_and_ordered_optimizer_replay',
        observed_CST=datetime.now(timezone(timedelta(hours=8))).isoformat(),raw_sacct=raw,rows=rows,
        native_SCFs_reparsed=sum(len(r['native_records']) for r in rows),scope_ids=[0,1,2,3],
        basis_40_by_36=basis.tolist(),free_basis_36_by_34=free.tolist(),
        origin_internal_A=(basis.T@np.asarray(data['origin_record']['coordinates'])).tolist(),
        source_LF_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        complete_ten_candidate_audit=False,stable_branch_or_TS_or_predictive_gain_certified=False,
        new_DFT_or_scheduler_writes=0,physical_plus005_or_private_arrays_accessed=False,
        original_2D_refusal_preserved=True,all_prior_and_other_training_costs_additional=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('fresh independent audit output required')
    result=native_prefix()
    with a.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(status=result['status'],native_SCFs_reparsed=result['native_SCFs_reparsed'],
        rows=[dict(id=r['id'],status=r['producer_status'],replay_status=r['replay']['status'],
              accepted_steps=(len(r['replay']['accepted_trace'])-1 if 'accepted_trace' in r['replay'] else None)) for r in result['rows']])))
