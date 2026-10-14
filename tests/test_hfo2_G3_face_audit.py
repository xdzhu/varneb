"""Read-only audit tests: injected arithmetic/optimizer records are not DFT."""
from copy import deepcopy
import hashlib
import json

import numpy as np
import pytest

from scripts import audit_hfo2_G3_conditional_faces as audit


@pytest.fixture
def example():
    basis = np.eye(40)[:, :36]
    axes = np.eye(36)[:, :2]
    internal = np.zeros(36)
    internal[0], internal[2] = -.1, .01
    full = basis @ internal
    gradient = np.zeros(40)
    gradient[2] = .01
    record = dict(status='complete', coordinates=full.tolist(), gradient=gradient.tolist(),
                  enthalpy_eV=-4., validation=dict(elapsed_seconds=2., explicit_toy_not_native=True))
    row = dict(id=1, index=[-1, 0], q_A=[-.1, 0.], status='gradient_candidate_release_stability_unchecked',
        new_completed_SCFs=1, reserved_attempts=1, callback_requests=2, center_cache_hits=0,
        sum_completed_native_call_seconds=2., physical_gradient_eV_A=gradient.tolist(),
        native_source=record['validation'], geometry_screens=dict(checks=1, rejected=0),
        point=dict(coordinates_internal_A=internal.tolist(), coordinates_full_A=full.tolist(),
            energy_eV_cell=-4., orthogonal_gradient_norm_eV_A=.01, callback_evaluations=2))
    data=dict(basis=basis, axes=axes, chart=None, hessian=np.eye(36),
              center_energy=0., center_gradient=np.zeros(40))
    return data, row, [record]


def test_nonzero_Q_candidate_uses_exact_native_final_not_last_trial(example):
    data,row,records=deepcopy(example)
    last=deepcopy(records[0]);last['coordinates'][2]=.04
    last['gradient'][2]=2.;last['enthalpy_eV']=99.
    records.append(last)
    row.update(new_completed_SCFs=2,reserved_attempts=2,callback_requests=3,sum_completed_native_call_seconds=4.)
    row['point']['callback_evaluations']=3
    result=audit.node_diagnostics(data,row,records)
    assert result['candidate']['energy_eV_cell']==-4.
    assert result['candidate']['released_gradient_norm_eV_A']==pytest.approx(.01)
    assert result['candidate']['released_stability_unchecked']
    assert result['completed_SCFs']==2


@pytest.mark.parametrize('fault', ['budget','counter_type','control','different_Q','nonfinite',
    'wrong_final','wrong_embedding','wrong_energy','gradient_target','wrong_gradient',
    'wrong_source','negative_seconds','wrong_cost','callback_type','callback_count','false_status'])
def test_corrupt_node_is_refused_without_promoting_a_trial(example,fault):
    data,row,records=deepcopy(example)
    if fault=='budget':row['reserved_attempts']=2
    elif fault=='counter_type':row['new_completed_SCFs']=True
    elif fault=='control':records[0]['coordinates'][-1]=.005
    elif fault=='different_Q':records[0]['coordinates'][0]=.1
    elif fault=='nonfinite':records[0]['gradient'][3]=float('nan')
    elif fault=='wrong_final':row['point']['coordinates_full_A'][2]=.02
    elif fault=='wrong_embedding':row['point']['coordinates_internal_A'][3]=.1
    elif fault=='wrong_energy':row['point']['energy_eV_cell']+=1
    elif fault=='gradient_target':records[0]['gradient'][2]=.04
    elif fault=='wrong_gradient':row['physical_gradient_eV_A'][2]=.02
    elif fault=='wrong_source':row['native_source']={}
    elif fault=='negative_seconds':records[0]['validation']['elapsed_seconds']=-1
    elif fault=='wrong_cost':row['sum_completed_native_call_seconds']=3
    elif fault=='callback_type':row['callback_requests']=True
    elif fault=='callback_count':row['point']['callback_evaluations']=1
    elif fault=='false_status':row['status']='certified_TS'
    with pytest.raises(ValueError):audit.node_diagnostics(data,row,records)


def test_budget_exhausted_hole_stays_hole_even_when_a_trial_has_small_gradient(example):
    data,row,records=deepcopy(example)
    records*=12
    row.update(point=None,status='optimizer_not_converged',new_completed_SCFs=12,reserved_attempts=12,
        sum_completed_native_call_seconds=24.,reason='mode evaluation attempt budget exhausted; calculator was not called')
    result=audit.node_diagnostics(data,row,records)
    assert result['candidate'] is None and result['completed_SCFs']==12
    row['reserved_attempts']=row['new_completed_SCFs']=1
    row['sum_completed_native_call_seconds']=2.
    with pytest.raises(ValueError,match='exhaust'):audit.node_diagnostics(data,row,records[:1])


def test_nonconvergence_reason_cannot_hide_a_backend_error(example):
    data,row,records=deepcopy(example)
    row.update(point=None,status='optimizer_not_converged',reason='ABACUS failed')
    with pytest.raises(ValueError,match='registered nonconvergence'):audit.node_diagnostics(data,row,records)


def test_complete_JSON_prefix_skips_inflight_file_without_a_failure_claim(tmp_path):
    specs=audit.worker.ordered_points()
    for spec in specs[:3]:
        (tmp_path/f"point_{spec['id']:03d}.json").write_text(json.dumps(dict(spec,status='explicit_toy')))
    (tmp_path/'point_003.json').write_text('{"point_id":')
    rows,hashes,pending=audit.load_prefix(tmp_path,specs)
    assert len(rows)==len(hashes)==3
    assert pending==dict(point_id=3,reason='point_JSON_not_yet_complete')
    assert all(hashlib.sha256(open(path,'rb').read()).hexdigest()==sha for path,sha in hashes.items())


@pytest.mark.parametrize('fault',['wrong_id','wrong_Q','boolean_id'])
def test_complete_but_wrong_point_identity_is_not_treated_as_pending(tmp_path,fault):
    spec=deepcopy(audit.worker.ordered_points()[0]);row=deepcopy(spec)
    if fault=='wrong_id':row['id']=1
    elif fault=='wrong_Q':row['q_A']=[.005,0.]
    else:row['id']=False
    (tmp_path/'point_000.json').write_text(json.dumps(row))
    with pytest.raises(ValueError,match='identity'):audit.load_prefix(tmp_path,[spec])


def test_missing_point_stops_prefix_not_skip_to_later_nodes(tmp_path):
    specs=audit.worker.ordered_points()
    (tmp_path/'point_000.json').write_text(json.dumps(specs[0]))
    (tmp_path/'point_002.json').write_text(json.dumps(specs[2]))
    rows,_,pending=audit.load_prefix(tmp_path,specs)
    assert len(rows)==1 and pending['point_id']==1 and pending['reason']=='point_file_not_yet_present'


def test_cache_only_nonzero_Q_replay_and_changed_trial_refusal(example,monkeypatch):
    data,row,records=deepcopy(example)
    monkeypatch.setattr(audit.worker.canary,'geometry_guard',lambda *a:True)
    monkeypatch.setattr(audit.study,'make_factory',lambda *a,**kw:pytest.fail('auditor must never create calculator'))
    start=np.asarray(row['point']['coordinates_internal_A'])
    result=audit.replay_node(data,row,records,start)
    assert result==dict(original_order_cache_only_replay_passed=True,new_DFT=0)
    wrong=deepcopy(records);wrong[0]['coordinates'][2]=.02
    with pytest.raises(ValueError,match='missing or differently ordered'):audit.replay_node(data,row,wrong,start)


def test_existing_output_refuses_before_native_or_scheduler_access(tmp_path,monkeypatch):
    path=tmp_path/'old.json';path.write_text('preserve')
    monkeypatch.setattr('sys.argv',['audit','--channel','T','--mode','prefix','--output',str(path)])
    monkeypatch.setattr(audit,'audit',lambda *a:pytest.fail('must refuse before native access'))
    with pytest.raises(FileExistsError):audit.main()
    assert path.read_text()=='preserve'


def test_audit_cannot_import_working_checkout_as_pinned_native_runtime():
    with pytest.raises(ValueError,match='immutable E161'):audit.audit('T','prefix')


@pytest.mark.parametrize('channel,pin,nodes,SCFs',[
    ('T','1a7b3530bb8d389196b9ad747df4f642e405ca1cdbb2d8da38d0ad529fefeda6',9,29),
    ('preserving','3c2c90ab703c2e1cc7c9c699daee5e0ded8f0d54d1b5d038c33a862e1b89b8ee',6,19)])
def test_actual_native_prefix_preserves_candidates_full_gradients_and_scope(channel,pin,nodes,SCFs):
    case=audit.worker.BASE/'G3_face_prefix_audit_E163_20261015'
    raw=(case/f'{channel}_native_prefix.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest()==pin
    report=json.loads(raw)
    assert report['audit_mode']=='prefix' and report['job_id']==audit.JOBS[channel]
    assert report['completed_nodes_audited']==report['gradient_candidates']==nodes
    assert report['new_SCFs_reparsed']==report['reserved_attempts_in_audited_prefix']==SCFs
    assert len(report['native_record_sources'])==SCFs
    assert report['cached_origin_new_SCFs']==report['new_DFT_or_scheduler_writes']==0
    assert report['sacct'].split('|')[1:3]==['RUNNING','0:0']
    assert report['sacct'].split('|')[4]=='32' and report['sacct'].split('|')[6]=='hfacnormal01'
    data=audit.worker.load_inputs(channel);basis,axes=data['basis'],data['axes']
    specs=audit.worker.ordered_points()
    assert [r['id'] for r in report['candidate_rows']]==list(range(1,nodes))
    for row in report['candidate_rows']:
        point=row['candidate'];u=np.asarray(point['coordinates_internal_A'])
        np.testing.assert_allclose(basis@u,point['coordinates_full_A'],atol=2e-12,rtol=0)
        np.testing.assert_allclose(axes.T@u,specs[row['id']]['q_A'],atol=2e-12,rtol=0)
        assert point['coordinates_full_A'][-1]==0.
        g=basis.T@np.asarray(point['full_physical_gradient_eV_A'])
        residual=np.linalg.norm(g-axes@(axes.T@g))
        assert residual==pytest.approx(point['released_gradient_norm_eV_A'],abs=2e-12,rel=0)
        assert residual<=.03 and point['released_stability_unchecked']
        np.testing.assert_allclose(axes.T@g,point['retained_gradient_components_eV_A'],atol=2e-12,rtol=0)
        assert row['replay']==dict(original_order_cache_only_replay_passed=True,new_DFT=0)
    for row in report['native_record_sources']:
        native=row['validation']
        assert native['genuine_MPI_ranks']==32 and native['raw_SCFs']==1
        assert all(native['input_sha256'][name]==sha for name,sha in audit.study.CONTRACT.items())
    assert report['completed_face_sha256'] is None and report['pending_node']['point_id']==nodes
    assert not report['unfinished_node_failure_paths']
    assert not report['physical_reserved_geometry_generated_or_read']
    assert not report['released_stability_branch_uniqueness_contour_TS_or_forecast_certified']
