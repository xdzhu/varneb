"""Released-curvature audit gates; synthetic records below are not DFT."""
from pathlib import Path
from copy import deepcopy
from types import SimpleNamespace
import importlib.util,json
import numpy as np
import pytest
from scripts import hfo2_G3_conditional_faces as faces

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/G3_released_audit_E183_20261015'
spec=importlib.util.spec_from_file_location('E183_audit',CASE/'audit.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
CONFIG=json.loads((ROOT/'benchmarks/hfo2_channels/20261008/G3_released_curvature_E181_20261015/config.json').read_bytes())


def account(state='COMPLETED',exit='0:0',cpus=32,nodes=1):
    return f'28868827|{state}|{exit}|100|{cpus}|node321|hfacnormal01|{nodes}'


def test_pending_is_not_completed_and_actual_allocation_cost_is_not_budget():
    p=a.accounting(account('PENDING',cpus=0,nodes=0));assert not p['normal_completion']
    assert p['allocation_core_hours']==0
    p=a.accounting(account());assert p['normal_completion']
    assert p['allocation_core_hours']==pytest.approx(100*32/3600)
    assert not a.accounting(account(exit='1:0'))['normal_completion']


@pytest.mark.parametrize('ranks',['1','0','32\nDSIZE = 32'])
def test_native_log_singleton_or_duplicate_MPI_groups_rejected_before_EFS_parser(tmp_path,ranks):
    from scripts.audit_hfo2_static_replica import audited_results
    out=tmp_path/'OUT.ABACUS';out.mkdir()
    (out/'running_scf.log').write_text('charge density convergence is achieved\nDSIZE = '+ranks+'\n')
    with pytest.raises(ValueError,match='one genuine DSIZE=32'):audited_results(tmp_path)


@pytest.mark.parametrize('broken',[account(cpus=64),account(nodes=2),account().replace('28868827','1'),account().replace('hfacnormal01','other'),''])
def test_wrong_or_missing_actual_accounting_refuses(broken):
    with pytest.raises(ValueError):a.accounting(broken)


def records():
    s=CONFIG['channels']['T'];d=np.asarray(s['literal_directions_40_by_34'])
    eigen=np.arange(1,35,dtype=float);eigen[0]=-2.;h=d@np.diag(eigen)@d.T
    base=np.asarray(s['literal_native_anchor_40_A']);g=np.linspace(-.1,.1,40)
    rows=[dict(gradient=(g+h@(np.asarray(p['coordinates_full_A'])-base)).tolist(),
        enthalpy_eV=float(i),validation=dict(elapsed_seconds=1.)) for i,p in enumerate(s['probes'])]
    return s,rows,eigen


def test_both_signed_amplitudes_full_actions_and_negative_sign_retained():
    s,r,eigen=records();blocks=a.reconstruct(s,r)
    assert [b['step_A'] for b in blocks]==[.01,.02]
    for b in blocks:
        assert b['raw_negative_count']==1
        np.testing.assert_allclose(b['raw_projected_eV_A2'],np.diag(eigen),atol=2e-12,rtol=0)
    with pytest.raises(ValueError):a.reconstruct(s,r[:68])
    r[2]['gradient'][0]=np.nan
    with pytest.raises(ValueError):a.reconstruct(s,r)


@pytest.mark.parametrize('fault',['coordinate','index','budget','type'])
def test_reservation_is_exact_literal_stencil_not_nearest_cache(fault):
    p=CONFIG['channels']['T']['probes'][0];x=np.asarray(p['coordinates_full_A'])
    b=dict(schema_version=1,max_evaluation_attempts=136,physical_contract_sha256='exact')
    r=dict(schema_version=1,evaluation_index=0,coordinate_sha256=a.hashlib.sha256(x.tobytes()).hexdigest(),
        coordinates=x.tolist(),budget_sha256=a._fingerprint(b))
    assert a.reserve_identity(r,p,0,b).startswith('eval-000000-')
    if fault=='coordinate':r['coordinates'][0]=np.nextafter(r['coordinates'][0],1.)
    elif fault=='index':r['evaluation_index']=1
    elif fault=='budget':r['budget_sha256']='changed'
    else:r['evaluation_index']=False
    with pytest.raises(ValueError):a.reserve_identity(r,p,0,b)


def test_pending_reservation_failure_and_terminal_refusal_are_separate(tmp_path,monkeypatch):
    d=faces.load_inputs('T');s=CONFIG['channels']['T'];p=tmp_path/'T/evaluations'
    w=SimpleNamespace(safe_probe=lambda *args:None)
    assert a.channel_prefix(w,d,s,p,complete=False)['pending']=='not_started'
    with pytest.raises(ValueError):a.channel_prefix(w,d,s,p,complete=True)


@pytest.fixture
def synthetic_directory(tmp_path):
    d=faces.load_inputs('T');s=CONFIG['channels']['T'];p=tmp_path/'T/evaluations';p.mkdir(parents=True)
    c=a.physical_contract(d['chart']);b=dict(schema_version=1,max_evaluation_attempts=136,physical_contract_sha256=c['contract_sha256'])
    (p/'contract.json').write_text(json.dumps(c));(p/'evaluation-budget.json').write_text(json.dumps(b))
    return d,s,p,b


def test_charged_pending_is_not_failure_and_bad_JSON_is_not_terminal_success(synthetic_directory):
    d,s,p,b=synthetic_directory;x=np.asarray(s['probes'][0]['coordinates_full_A']);dg=a.hashlib.sha256(x.tobytes()).hexdigest()
    marker=dict(schema_version=1,evaluation_index=0,coordinate_sha256=dg,coordinates=x.tolist(),budget_sha256=a._fingerprint(b))
    (p/'attempt-000000.json').write_text(json.dumps(marker))
    w=SimpleNamespace(safe_probe=lambda *args:None)
    r=a.channel_prefix(w,d,s,p,complete=False)
    assert r['reserved']==1 and r['pending']=='reserved_probe_0_output_pending' and r['failures']==[]
    directory=p/f'eval-000000-{dg[:12]}';directory.mkdir();(directory/'failure.json').write_text('{}')
    r=a.channel_prefix(w,d,s,p,complete=False)
    assert len(r['failures'])==1 and not r['records']
    with pytest.raises(ValueError):a.channel_prefix(w,d,s,p,complete=True)
    (p/'attempt-000000.json').write_text('{')
    assert a.channel_prefix(w,d,s,p,complete=False)['pending']=='reservation_JSON_inflight'
    with pytest.raises(ValueError):a.channel_prefix(w,d,s,p,complete=True)


def test_complete_synthetic_inventory_is_checked_not_only_empty_prefix(synthetic_directory,monkeypatch):
    d,s,p,b=synthetic_directory;_,rows,_=records()
    for i,probe in enumerate(s['probes']):
        x=np.asarray(probe['coordinates_full_A']);dg=a.hashlib.sha256(x.tobytes()).hexdigest()
        marker=dict(schema_version=1,evaluation_index=i,coordinate_sha256=dg,coordinates=x.tolist(),budget_sha256=a._fingerprint(b))
        (p/f'attempt-{i:06d}.json').write_text(json.dumps(marker))
        directory=p/f'eval-{i:06d}-{dg[:12]}';directory.mkdir();(directory/'result.json').write_text(json.dumps(rows[i]))
    def native(w,data,x,path,namespace):
        i=int(path.parent.name.split('-')[1]);return rows[i],dict(result_sha256=a.sha(path))
    monkeypatch.setattr(a,'native_record',native)
    w=SimpleNamespace(safe_probe=lambda *args:None)
    # Complete inventory without a terminal channel summary is still refused.
    with pytest.raises(ValueError,match='summary missing'):a.channel_prefix(w,d,s,p,complete=True)
    r=a.channel_prefix(w,d,s,p,complete=False)
    assert len(r['records'])==r['reserved']==136 and len(r['blocks'])==2
    summary=dict(status='complete_post_release_raw_transverse_curvature_not_stability_certificate',channel='T',
        new_native_SCFs=136,new_center_SCFs=0,base_released_gradient_norm_eV_A=s['base_released_gradient_norm_eV_A'],
        physical_stability_TS_or_independent_prediction_certified=False,
        energies_eV_cell=[q['enthalpy_eV'] for q in rows],full_gradients_eV_A=[q['gradient'] for q in rows],
        sum_native_call_seconds=136.,call_time_proxy_core_hours=136.*32/3600,blocks=r['blocks'],
        two_step_raw_block_operator_difference_eV_A2=float(np.linalg.norm(np.asarray(r['blocks'][0]['raw_projected_eV_A2'])-r['blocks'][1]['raw_projected_eV_A2'],ord=2)))
    (p.parent/'curvature.json').write_text(json.dumps(summary))
    assert len(a.channel_prefix(w,d,s,p,complete=True)['records'])==136
    for key,value in [('new_native_SCFs',135),('new_center_SCFs',1),('physical_stability_TS_or_independent_prediction_certified',True)]:
        changed=deepcopy(summary);changed[key]=value
        with pytest.raises(ValueError):a.check_summary(changed,'T',s,rows,r['blocks'])
    changed=deepcopy(summary);changed['blocks'][0]['raw_projected_eV_A2'][0][0]*=-1
    with pytest.raises(ValueError):a.check_summary(changed,'T',s,rows,r['blocks'])


def test_complete_audit_requires_normal_job_before_reading_any_record(monkeypatch):
    config=dict(policy=dict(channels=['T','preserving']))
    monkeypatch.setattr(a,'load_worker',lambda:(None,config,{}))
    monkeypatch.setattr(a.subprocess,'run',lambda *args,**kwargs:SimpleNamespace(stdout=account('PENDING',cpus=0,nodes=0)))
    monkeypatch.setattr(a,'channel_prefix',lambda *args,**kwargs:pytest.fail('must refuse before data access'))
    for n in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):monkeypatch.setenv(n,'1')
    with pytest.raises(RuntimeError,match='normal Slurm completion'):a.audit('complete')


def test_normal_job_with_missing_channel_never_becomes_a_complete_report(monkeypatch):
    config=dict(policy=dict(channels=['T','preserving']),channels={'T':{},'preserving':{}})
    monkeypatch.setattr(a,'load_worker',lambda:(None,config,{'T':{},'preserving':{}}))
    monkeypatch.setattr(a.subprocess,'run',lambda *args,**kwargs:SimpleNamespace(stdout=account()))
    monkeypatch.setattr(a,'channel_prefix',lambda *args,**kwargs:dict(records=[],reserved=0,blocks=None))
    for n in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):monkeypatch.setenv(n,'1')
    with pytest.raises(ValueError,match='all136'):a.audit('complete')
