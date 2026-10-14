"""Synthetic bookkeeping tests; completed trials never prove accepted minima."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import zipfile

from scripts import audit_hfo2_G3_conditional_canary as audit

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'benchmarks/hfo2_channels/20261008'


@pytest.fixture
def example():
    manifest=json.loads((BASE/'G3_curvature_E150_20261014/manifest.json').read_bytes())
    config=json.loads((BASE/'G3_complete_curvature_E159_20261014/T_config.json').read_bytes())
    basis=np.asarray(manifest['centers']['T']['directions'])[:,:36]
    axes=np.asarray(config['internal_axes_36_by_2'])
    stable=np.linalg.svd(axes.T,full_matrices=True)[2][2:].T
    internal=stable[:,0]*.02
    row=dict(status='complete',coordinates=(basis@internal).tolist(),
        gradient=(basis@stable[:,1]*.01).tolist(),enthalpy_eV=-4.,validation=dict(elapsed_seconds=2.))
    outcome=dict(arm='identity',status='gradient_converged_candidate_curvature_unchecked',
        new_completed_SCFs=1,reserved_attempts=1,center_cache_hits=1,callback_requests=2,
        sum_completed_native_call_seconds=2.,completed_call_proxy_core_hours=2.*32/3600,
        point=dict(coordinates_internal_A=internal.tolist(),coordinates_full_A=row['coordinates'],
            energy_eV_cell=-4.,orthogonal_gradient_norm_eV_A=.01,callback_evaluations=2))
    # All energies/gradients above are synthetic; only the coordinate span comes
    # from the measured manifest. No native SCF provenance is claimed.
    return basis,axes,-3.,np.zeros(40),outcome,[row]


def test_explicit_final_candidate_matches_native_coordinate_not_the_last_trial(example):
    *args,outcome,records=deepcopy(example)
    extra=deepcopy(records[0]);extra['coordinates']=(np.asarray(extra['coordinates'])*2).tolist()
    extra['gradient']=(np.asarray(extra['gradient'])*100).tolist();extra['enthalpy_eV']=100.
    records.append(extra)
    outcome.update(new_completed_SCFs=2,reserved_attempts=2,callback_requests=3,
        sum_completed_native_call_seconds=4.,completed_call_proxy_core_hours=4.*32/3600)
    outcome['point']['callback_evaluations']=3
    result=audit.arm_diagnostics(*args,outcome,records)
    assert result['final_candidate']['energy_eV_cell']==-4.
    assert result['final_candidate']['released_gradient_norm_eV_A']==pytest.approx(.01)
    assert result['native_trial_rows'][-1]['released_gradient_norm_eV_A']>audit.worker.SETTINGS['gradient_tolerance_eV_A']
    assert all(not row['accepted_iterate_known'] for row in result['native_trial_rows'])


def test_low_gradient_trial_is_not_promoted_after_real_budget_exhaustion(example):
    *args,outcome,records=deepcopy(example)
    records=records*12
    outcome.update(status='not_converged',point=None,new_completed_SCFs=12,reserved_attempts=12,
        callback_requests=13,sum_completed_native_call_seconds=24.,completed_call_proxy_core_hours=24.*32/3600,
        reason='mode evaluation attempt budget exhausted; calculator was not called')
    result=audit.arm_diagnostics(*args,outcome,records)
    assert result['final_candidate'] is None and result['conservative_attempt_reservations']==12
    assert all(row['released_gradient_norm_eV_A']<.03 for row in result['native_trial_rows'])


@pytest.mark.parametrize('fault',['false_budget','wrong_count','false_convergence','wrong_final_record',
    'fixed_control','different_Q','wrong_energy','wrong_gradient','negative_time','wrong_cost','boolean_counter'])
def test_corrupted_or_misclassified_outcomes_are_refused(example,fault):
    *args,outcome,records=deepcopy(example)
    if fault=='false_budget':
        outcome.update(status='not_converged',point=None,
            reason='mode evaluation attempt budget exhausted; calculator was not called')
    elif fault=='wrong_count':outcome['reserved_attempts']=0
    elif fault=='false_convergence':records[0]['gradient']=(np.asarray(records[0]['gradient'])*100).tolist()
    elif fault=='wrong_final_record':outcome['point']['coordinates_full_A'][0]+=.001
    elif fault=='fixed_control':records[0]['coordinates'][-1]=.001
    elif fault=='different_Q':
        records[0]['coordinates']=(np.asarray(records[0]['coordinates'])+args[0]@args[1][:,0]*.01).tolist()
    elif fault=='wrong_energy':outcome['point']['energy_eV_cell']+=1.
    elif fault=='wrong_gradient':outcome['point']['orthogonal_gradient_norm_eV_A']+=.01
    elif fault=='negative_time':records[0]['validation']['elapsed_seconds']=-1.
    elif fault=='wrong_cost':outcome['completed_call_proxy_core_hours']+=1.
    elif fault=='boolean_counter':outcome['new_completed_SCFs']=True
    with pytest.raises(ValueError):audit.arm_diagnostics(*args,outcome,records)


def test_existing_output_refuses_before_native_read_or_scheduler_query(tmp_path,monkeypatch):
    output=tmp_path/'existing.json';output.write_text('preserve')
    monkeypatch.setattr('sys.argv',['audit','--channel','T','--curvature-report','unused','--output',str(output)])
    calls=[];monkeypatch.setattr(audit,'audit',lambda *a:calls.append(True))
    with pytest.raises(FileExistsError):audit.main()
    assert not calls and output.read_text()=='preserve'


def test_cache_only_oracle_preserves_order_cache_and_hard_attempt_cap(example):
    records=example[-1];oracle=audit.CacheOnly(records)
    x=np.asarray(records[0]['coordinates'])
    first=oracle(x);again=oracle(x)
    assert first[0]==again[0] and len(oracle.seen)==1
    first[1][:]=99
    assert not np.any(oracle(x)[1]==99)
    with pytest.raises(ValueError,match='missing or differently ordered'):oracle(x*2)
    assert len(oracle.seen)==1
    capped=audit.CacheOnly([dict(records[0],coordinates=(x*(i+1)).tolist()) for i in range(12)])
    for i in range(12):capped(x*(i+1))
    with pytest.raises(RuntimeError,match='attempt budget exhausted'):capped(x*13)
    assert len(capped.seen)==12


def test_exact_origin_cache_replay_never_calls_an_oracle_or_changes_geometry(example,monkeypatch):
    basis,axes,energy,_,_,_=example
    config=dict(internal_axes_36_by_2=axes.tolist())
    outcome=dict(arm='identity',point=dict(coordinates_internal_A=np.zeros(36).tolist(),
        energy_eV_cell=energy,orthogonal_gradient_norm_eV_A=0.),
        callback_requests=2,center_cache_hits=2,trial_geometry_and_domain_screens=dict(checks=1,rejected=0))
    monkeypatch.setattr(audit.worker,'geometry_guard',lambda *a:True)
    result=audit.replay_arm(None,basis,config,energy,np.zeros(40),np.eye(36),outcome,[])
    assert result['original_order_and_cache_only_optimizer_replay_passed'] and result['new_DFT']==0


@pytest.mark.parametrize('channel,zip_sha,audit_sha,seed_count,residual',[
    ('preserving','c1955f42e1d9142235971a83eae0fc52c18b4c65521ec5e731e3ed1d2d4c28f9',
     '336071744885656b25515f28ca9a6d6b4ecd09f9b411a739827ec7d1e56e5403',5,.01959108478103186),
    ('T','82fc05455912488051f52da2396606cbc060fb9be1c06ab42e07487768b2a44f',
     'b4097af45e0de951c32e550141d745a2d73b6e13eeb4aa8fb3a62be302afa685',2,.02629995445855942)])
def test_actual_complete_native_rows_and_bounded_outcome(channel,zip_sha,audit_sha,seed_count,residual):
    case=BASE/'G3_conditional_canary_E160_20261014'
    data=(case/(channel+'_worker_outputs.zip')).read_bytes()
    assert hashlib.sha256(data).hexdigest()==zip_sha
    raw=(case/(channel+'_native_audit.json')).read_bytes()
    assert hashlib.sha256(raw).hexdigest()==audit_sha
    native=json.loads(raw)
    manifest=json.loads((BASE/'G3_curvature_E150_20261014/manifest.json').read_bytes())
    config=json.loads((BASE/'G3_complete_curvature_E159_20261014'/(channel+'_config.json')).read_bytes())
    basis=np.asarray(manifest['centers'][channel]['directions'])[:,:36]
    axes=np.asarray(config['internal_axes_36_by_2'])
    with zipfile.ZipFile(BASE/'G3_complete_curvature_E159_20261014/native_reports.zip') as original:
        center=json.loads(original.read(channel+'_native.json'))
    replayed=[]
    with zipfile.ZipFile(case/(channel+'_worker_outputs.zip')) as archive:
        summary=json.loads(archive.read(channel+'/canary.json'))
        for index,arm in enumerate(audit.worker.SETTINGS['arms']):
            names=sorted(n for n in archive.namelist() if n.startswith(channel+'/'+arm+'/eval-') and n.endswith('/result.json'))
            records=[json.loads(archive.read(n)) for n in names]
            row=audit.arm_diagnostics(basis,axes,center['center_source']['energy_eV_cell'],
                center['raw_center_gradient_eV_A'],summary['arms'][index],records)
            assert row['completed_SCFs']==native['arms'][index]['completed_SCFs']
            np.testing.assert_allclose([r['released_gradient_norm_eV_A'] for r in row['native_trial_rows']],
                [r['released_gradient_norm_eV_A'] for r in native['arms'][index]['native_trial_rows']],atol=2e-12,rtol=0)
            assert native['arms'][index]['optimizer_replay']['original_order_and_cache_only_optimizer_replay_passed']
            replayed.append(row)
    assert [r['completed_SCFs'] for r in replayed]==[12,seed_count]
    assert replayed[0]['final_candidate'] is None
    assert replayed[1]['final_candidate']['released_gradient_norm_eV_A']==pytest.approx(residual,abs=2e-12)
    assert native['comparison'] is None and not native['conditional_stability_TS_full_face_or_material_speedup_certified']


def test_actual_ordered_release_diagnostic_CLI_uses_both_candidates_without_DFT():
    case=BASE/'G3_conditional_canary_E160_20261014'
    result=subprocess.run([sys.executable,str(case/'analyze_preserving_release.py')],
        cwd=ROOT,check=True,capture_output=True,text=True)
    actual=json.loads(result.stdout);expected=json.loads((case/'preserving_release_diagnostic.json').read_bytes())
    for key in ('energy_change_eV_cell','energy_difference_from_T_conditional_candidate_eV_cell',
                'internal_displacement_norm_A','retained_gradient_norm_eV_A','released_gradient_norm_eV_A'):
        assert actual[key]==pytest.approx(expected[key],abs=2e-12,rel=0)
    assert len(actual['comparisons'])==4
    for a,e in zip(actual['comparisons'],expected['comparisons']):
        assert a.keys()==e.keys()
        for key in a:
            if isinstance(e[key],float):assert a[key]==pytest.approx(e[key],abs=2e-12,rel=0)
            else:assert a[key]==e[key]
    assert not actual['same_phase_branch_TS_or_barrier_reduction_certified'] and actual['new_DFT']==0
