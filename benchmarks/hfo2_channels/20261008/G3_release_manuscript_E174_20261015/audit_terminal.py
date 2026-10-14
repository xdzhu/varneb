"""Read-only terminal native/replay audit; no calculator, scheduler writes or resume."""
from datetime import datetime, timezone, timedelta
from pathlib import Path
import hashlib, importlib.util, json, os, subprocess, sys
sys.dont_write_bytecode = True
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator

RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS / '20261015-varneb-G3-feasible-execution-E171-r1'
OUT = RUNS / '20261015-varneb-G3-terminal-audit-E174-r1'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def require_capped_outcome(row, n):
    if (row['job_id'] != '28866132' or row['status'] != 'optimizer_not_converged'
            or row['candidate'] is not None
            or row['reason'] != 'mode evaluation attempt budget exhausted; calculator was not called'
            or n != 12 or row['new_completed_SCFs'] != n or row['new_reserved_attempts'] != n
            or row['net_speedup_ratio'] is not None
            or row['old_baseline_recomputed_SCFs'] != 0
            or row['physical_reserved_geometry_generated_or_read']
            or row['stable_conditional_minimum_or_TS_or_independent_prediction_certified']):
        raise ValueError('exact capped nonconvergence and no false promotion required')


def trial_diagnostics(data, rows):
    trials = []
    for i, row in enumerate(rows):
        x = np.asarray(row['coordinates']); g = np.asarray(row['gradient'])
        internal = data['basis'].T @ x; gi = data['basis'].T @ g
        q = data['axes'].T @ internal
        if (not np.allclose(q, [-.2,-.1], atol=2e-12, rtol=0) or x[-1] != 0
                or not np.isfinite(gi).all()):
            raise ValueError('trial fixed Q/control/finite gradient differs')
        r = internal-data['axes']@q
        gr = gi-data['axes']@(data['axes'].T@gi)
        direction = r/np.linalg.norm(r)
        radial = float(direction@gr)
        trials.append(dict(index=i, coordinate_sha256=row['coordinate_sha256'],
            energy_eV_cell=row['enthalpy_eV'], internal_radius_A=float(np.linalg.norm(internal)),
            released_gradient_norm_eV_A=float(np.linalg.norm(gr)),
            radial_released_gradient_eV_A=radial,
            tangential_released_gradient_eV_A=float(np.linalg.norm(gr-direction*radial)),
            accepted_iterate_known=False))
    return trials


def require_original_embedding(basis, internal, full, guard):
    """Use exact replayed optimizer coordinates, not a second basis round trip.

    No domain epsilon, inward move or raw-result edit is permitted. A full
    coordinate must be bitwise the embedding of a strictly admissible original
    internal proposal. The reprojected norm remains a separate diagnostic.
    """
    if not np.array_equal(basis @ internal, full) or not guard(internal):
        raise ValueError('exact original embedding and unchanged strict guard required')


def native_with_original_embedding(audit, data, directory, internals):
    """Existing E163 native contract/reservation/EFS checks, original-coordinate guard."""
    chart,basis=data['chart'],data['basis'];ref=chart.reference_atoms
    expected=dict(schema_version=1,coordinate_convention=chart.coordinate_convention,
        reference=dict(numbers=ref.numbers.tolist(),cell_A=ref.cell.array.tolist(),
            scaled_positions=ref.get_scaled_positions(wrap=False).tolist(),pbc=ref.pbc.tolist()),
        calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
        validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',pressure_eV_per_A3=0.,
        minimum_distance_A=1.6,mechanical_boundary=chart.mechanical_contract)
    fingerprint=audit._fingerprint(expected)
    contract=json.loads((directory/'contract.json').read_bytes())
    budget=dict(schema_version=1,max_evaluation_attempts=12,physical_contract_sha256=fingerprint)
    if (contract!=dict(expected,contract_sha256=fingerprint)
            or json.loads((directory/'evaluation-budget.json').read_bytes())!=budget):
        raise ValueError('exact original physical contract/budget required')
    paths=sorted(directory.glob('eval-*/result.json'));markers=sorted(directory.glob('attempt-*.json'))
    if (list(directory.glob('eval-*/failure.json')) or len(list(directory.glob('eval-*')))!=len(paths)
            or len(markers)!=len(paths) or len(paths)!=12):
        raise ValueError('all12 complete reserved new native records required')
    records=[];sources=[]
    for index,(p,marker) in enumerate(zip(paths,markers)):
        row=json.loads(p.read_bytes());x=np.asarray(row['coordinates']);digest=hashlib.sha256(x.tobytes()).hexdigest()
        reservation=dict(schema_version=1,evaluation_index=index,coordinate_sha256=digest,
            coordinates=x.tolist(),budget_sha256=audit._fingerprint(budget))
        if (p.parent.name!=f'eval-{index:06d}-{digest[:12]}' or marker.name!=f'attempt-{index:06d}.json'
                or json.loads(marker.read_bytes())!=reservation or row['contract_sha256']!=fingerprint
                or row['coordinate_sha256']!=digest or row['status']!='complete' or digest not in internals):
            raise ValueError('native identity/reservation/original internal proposal missing')
        atoms=chart.to_atoms(x)
        atoms.calc=SinglePointCalculator(atoms,**audit.study.audited_results(p.parent/'scf_000000'))
        checked=audit.study.result_validator(index,atoms,p.parent);e,g=chart.energy_and_gradient(x,atoms)
        if (checked!=row['validation'] or e!=row['enthalpy_eV']
                or not np.allclose(g,row['gradient'],atol=2e-12,rtol=0)
                or not np.array_equal(atoms.get_stress(voigt=False),row['stress_eV_per_A3'])
                or not np.array_equal(atoms.cell.array[:2],ref.cell.array[:2])):
            raise ValueError('native energy/force/stress/full gradient/substrate differs')
        require_original_embedding(basis,internals[digest],x,
            lambda v:audit.worker.canary.geometry_guard(chart,basis,v))
        records.append(row);sources.append(dict(result_sha256=sha(p),validation=checked,
            original_internal_A=internals[digest].tolist(),
            original_internal_norm_A=float(np.linalg.norm(internals[digest])),
            roundtrip_internal_norm_A=float(np.linalg.norm(basis.T@x)),
            roundtrip_guard=bool(audit.worker.canary.geometry_guard(chart,basis,basis.T@x))))
    return records,sources


def main():
    if any(os.environ.get(k)!='1' for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS')):
        raise ValueError('exact original single-thread numerical replay environment required')
    if OUT.exists(): raise FileExistsError('fresh audit namespace required')
    bundle = ROOT/'bundle'
    if sha(bundle/'bundle_files.json') != 'd4f7a2fbf1950001b725eed0ff1d286e98825c23d62d9a070ac5bf6d4aaa00d2':
        raise ValueError('immutable bundle manifest differs')
    files = json.loads((bundle/'bundle_files.json').read_bytes())
    if len(files)!=9 or any(sha(bundle/n)!=h for n,h in files.items()): raise ValueError('bundle drift')
    w=module('terminal_original_execution_worker',bundle/'worker.py')
    m,c,point,old,data=w.validate(json.loads((bundle/'execution.json').read_bytes()),bundle)
    w.verify_runtime(m)
    account=subprocess.run(['sacct','-n','-P','-X','-j','28866132',
        '--format=JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,NodeList,Partition'],
        check=True,capture_output=True,text=True).stdout.strip()
    fields=account.split('|')
    if (len(fields)<7 or fields[:3]!=['28866132','COMPLETED','0:0']
            or fields[4]!='32' or fields[6]!='hfacnormal01'):
        raise RuntimeError('exact normally completed original allocation required: '+account)
    if json.loads((ROOT/'accepted.json').read_bytes())['job_id']!='28866132':
        raise ValueError('accepted handle differs')
    summary=json.loads((ROOT/'arm/arm.json').read_bytes())
    canonical=ROOT/'arm/new_evaluations';paths=sorted(canonical.glob('eval-*/result.json'))
    require_capped_outcome(summary,len(paths))
    alias=module('terminal_verified_recorded_alias',Path(__file__).with_name('audit_first.py'))
    recorded=[]
    for p in paths:
        row=json.loads(p.read_bytes())
        recorded.append(alias.recorded_evaluation_directory(row['validation']['raw_directory'],p.parent,ROOT/'arm'))
    directory=recorded[0].parent
    if (directory.resolve(strict=True)!=canonical.resolve(strict=True)
            or any(p.parent!=directory for p in recorded)):
        raise ValueError('recorded native root is not one exact expected directory')
    audit=module('terminal_existing_E163_auditor',bundle/'audit_faces.py')
    fresh=[json.loads(p.read_bytes()) for p in paths]
    # Every original cached point is also re-read with the immutable auditor.
    old_native,_=audit.native_records(data,Path(m.POINT_PATH).parent/'node_013')
    if old_native!=old: raise ValueError('old cache differs from native baseline')
    baseline=audit.replay_node(data,point,old,np.array(point['starting_internal_coordinates_A']))
    cache=audit.CacheOnly(fresh)
    cache.contract=json.loads((directory/'contract.json').read_bytes())
    type(cache).n_new_evaluations=property(lambda self:len(self.seen))
    counted=w.CountedDelegate(cache);oracle=m.SharedExactCache(old,counted)
    evaluate=m.faces.canary.InternalEvaluator(w.ExactStartingOracle(oracle,old[0]['coordinates']),
        data['basis'],data['center_energy'],data['center_gradient'])
    plane=m.faces.canary.ModePlane(np.zeros(36),data['axes'],np.ones(36),np.eye(2),
        ('coupled internal1','coupled internal2'),'local length','G3-first-conditional-sheet')
    screens=dict(checks=0,rejected=0);internals={}
    def guard(x):
        screens['checks']+=1
        ok=m.faces.canary.geometry_guard(data['chart'],data['basis'],x)
        screens['rejected']+=int(not ok)
        if ok:
            digest=hashlib.sha256((data['basis']@x).tobytes()).hexdigest()
            if digest in internals and not np.array_equal(internals[digest],x):
                raise ValueError('one full coordinate maps to different replayed internal proposals')
            internals[digest]=np.array(x,copy=True)
        return ok
    solver=module('vcneb.feasible_relaxation',bundle/'feasible_relaxation.py')
    candidate=None;reason=None
    try:
        candidate=solver.relax_orthogonal_in_ball(plane,m.POLICY['q_A'],evaluate,
            start=np.array(point['starting_internal_coordinates_A']),
            orthogonal_radius_bound=c['released_radius_A'],trial_validator=guard,
            initial_hessian=data['hessian'],initial_trust_radius=.05,
            gradient_tolerance=.03,max_iterations=200)
    except RuntimeError as exc: reason=str(exc)
    checks=dict(capped_reason=reason==summary['reason'],no_candidate=candidate is None,
        all12_new_used=len(cache.seen)==12,old_cache_hits=oracle.old_cache_hits==summary['old_cache_hits'],
        physical_request_sequence=oracle.requests==summary['physical_requests'],
        delegate_requests=counted.requests==summary['new_delegate_requests'],
        internal_callbacks=evaluate.callback_requests==summary['internal_callback_requests'],
        center_hits=evaluate.center_cache_hits==summary['center_cache_hits'],
        screens=screens==summary['geometry_screens'],
        native_call_seconds=np.isclose(sum(r['validation']['elapsed_seconds'] for r in fresh),
            summary['sum_new_completed_call_seconds'],atol=1e-10,rtol=0))
    checks={k:bool(v) for k,v in checks.items()}
    if not all(checks.values()): raise ValueError('cache-only exact decisions/counters differ: '+str(checks))
    native,sources=native_with_original_embedding(audit,data,directory,internals)
    if native!=fresh: raise ValueError('native records differ from cache-only replay source')
    w.verify_runtime(m)
    trials=trial_diagnostics(data,fresh);oldtrials=trial_diagnostics(data,old)
    result=dict(status='normally_completed_capped_nonconvergence_native_and_replay_audited',
        observed_CST=datetime.now(timezone(timedelta(hours=8))).isoformat(),sacct=account,
        job_id='28866132',native_new_SCFs_reparsed=12,native_old_SCFs_reparsed=12,
        original_runtime_files_rechecked=380,immutable_bundle_files_rechecked=9,
        arm_summary_sha256=sha(ROOT/'arm/arm.json'),arm_summary=summary,native_sources=sources,
        old_baseline_replay=baseline,new_cache_only_replay_checks=checks,
        exact_original_internal_coordinates_and_unchanged_guard_used=True,
        roundtrip_guard_failures_kept_as_diagnostic=sum(not r['roundtrip_guard'] for r in sources),
        domain_epsilon_or_coordinate_repair_used=False,
        numerical_replay_thread_environment={k:os.environ[k] for k in
            ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS')},
        new_trials=trials,old_trials=oldtrials,
        minimum_new_trial_released_gradient_eV_A=min(r['released_gradient_norm_eV_A'] for r in trials),
        minimum_old_trial_released_gradient_eV_A=min(r['released_gradient_norm_eV_A'] for r in oldtrials),
        new_completed_call_proxy_core_hours=summary['sum_new_completed_call_seconds']*32/3600,
        new_allocation_core_hours=int(fields[3])*32/3600,
        additional_old_curvature_and_all_prior_preparation_costs=True,
        material_convergence_or_net_speedup_certified=False,
        original_censored_baseline_not_speedup_denominator=True,
        physical_reserved_geometry_generated_or_read=False,new_DFT_or_scheduler_writes=0)
    OUT.mkdir()
    with (OUT/'terminal_audit.json').open('x') as f: json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:result[k] for k in ('status','observed_CST','minimum_new_trial_released_gradient_eV_A',
        'minimum_old_trial_released_gradient_eV_A','new_completed_call_proxy_core_hours','new_allocation_core_hours')}))


if __name__=='__main__':main()
