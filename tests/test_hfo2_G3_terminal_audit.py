"""Terminal audit gates and trial-only diagnostics, not a new material benchmark."""
from pathlib import Path
from copy import deepcopy
import importlib.util
import hashlib,json,zipfile
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('terminal_audit_tests',ROOT/
    'benchmarks/hfo2_channels/20261008/G3_release_manuscript_E174_20261015/audit_terminal.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def capped():
    return dict(job_id='28866132',status='optimizer_not_converged',candidate=None,
        reason='mode evaluation attempt budget exhausted; calculator was not called',
        new_completed_SCFs=12,new_reserved_attempts=12,net_speedup_ratio=None,
        old_baseline_recomputed_SCFs=0,physical_reserved_geometry_generated_or_read=False,
        stable_conditional_minimum_or_TS_or_independent_prediction_certified=False)


def test_standalone_auditor_loads_solver_with_its_real_package_relative_imports():
    solver=m.module('vcneb.feasible_relaxation',ROOT/'vcneb/feasible_relaxation.py')
    assert callable(solver.relax_orthogonal_in_ball)


@pytest.mark.parametrize('key,value',[
    ('job_id','other'),('status','converged'),('candidate',{}),('reason','backend failed'),
    ('new_completed_SCFs',11),('new_reserved_attempts',13),('net_speedup_ratio',1.2),
    ('old_baseline_recomputed_SCFs',12),('physical_reserved_geometry_generated_or_read',True),
    ('stable_conditional_minimum_or_TS_or_independent_prediction_certified',True)])
def test_no_partial_wrong_handle_failure_or_candidate_promotion(key,value):
    row=capped();m.require_capped_outcome(row,12)
    bad=deepcopy(row);bad[key]=value
    with pytest.raises(ValueError,match='exact capped'):
        m.require_capped_outcome(bad,12)
    with pytest.raises(ValueError,match='exact capped'):
        m.require_capped_outcome(row,11)


def test_full_gradient_not_only_tangent_and_trials_are_not_accepted_iterates():
    basis=np.eye(40)[:,:36];axes=np.eye(36)[:,:2]
    x=np.zeros(40);x[:3]=[-.2,-.1,.4]
    g=np.zeros(40);g[2:4]=[-.04,.03]
    row=dict(coordinates=x.tolist(),gradient=g.tolist(),coordinate_sha256='synthetic',enthalpy_eV=1.)
    out=m.trial_diagnostics(dict(basis=basis,axes=axes),[row])[0]
    assert out['released_gradient_norm_eV_A']==pytest.approx(.05)
    assert out['radial_released_gradient_eV_A']==pytest.approx(-.04)
    assert out['tangential_released_gradient_eV_A']==pytest.approx(.03)
    assert not out['accepted_iterate_known']
    row['coordinates'][-1]=1e-4
    with pytest.raises(ValueError,match='fixed Q/control'):
        m.trial_diagnostics(dict(basis=basis,axes=axes),[row])


def test_original_embedding_guard_is_exact_not_an_epsilon_or_inward_repair():
    # A finite-precision basis round trip need not preserve a boundary norm.
    basis=np.eye(2)*np.nextafter(1.,2.)
    internal=np.array([.5,0.]);full=basis@internal
    guard=lambda x:np.linalg.norm(x)<=.5
    assert guard(internal) and not guard(basis.T@full)
    m.require_original_embedding(basis,internal,full,guard)
    changed=full.copy();changed[0]=np.nextafter(changed[0],0.)
    with pytest.raises(ValueError,match='exact original embedding'):
        m.require_original_embedding(basis,internal,changed,guard)
    outside=np.array([np.nextafter(.5,1.),0.])
    with pytest.raises(ValueError,match='unchanged strict guard'):
        m.require_original_embedding(basis,outside,basis@outside,guard)


def test_actual_terminal_native_proof_keeps_all_trials_costs_and_roundtrip_failure():
    case=ROOT/'benchmarks/hfo2_channels/20261008/G3_release_manuscript_E174_20261015'
    proof=json.loads((case/'terminal_audit.json').read_bytes())
    m.require_capped_outcome(proof['arm_summary'],12)
    assert proof['sacct'].split('|')[:3]==['28866132','COMPLETED','0:0']
    assert proof['native_new_SCFs_reparsed']==proof['native_old_SCFs_reparsed']==12
    assert all(proof['new_cache_only_replay_checks'].values())
    assert proof['old_baseline_replay']['original_order_cache_only_replay_passed']
    assert proof['exact_original_internal_coordinates_and_unchanged_guard_used']
    assert proof['roundtrip_guard_failures_kept_as_diagnostic']==1
    assert not proof['domain_epsilon_or_coordinate_repair_used']
    assert len(proof['native_sources'])==len(proof['new_trials'])==len(proof['old_trials'])==12
    assert all(s['original_internal_norm_A']<=.5 for s in proof['native_sources'])
    assert all(s['validation']['genuine_MPI_ranks']==32 for s in proof['native_sources'])
    assert not any(r['accepted_iterate_known'] for r in proof['new_trials']+proof['old_trials'])
    assert proof['minimum_new_trial_released_gradient_eV_A']==pytest.approx(.04585524009718605)
    assert proof['minimum_new_trial_released_gradient_eV_A']>.03
    assert proof['new_allocation_core_hours']==pytest.approx(1417*32/3600)
    assert proof['new_completed_call_proxy_core_hours']==pytest.approx(
        proof['arm_summary']['sum_new_completed_call_seconds']*32/3600)
    assert proof['new_DFT_or_scheduler_writes']==0
    assert not proof['material_convergence_or_net_speedup_certified']
    assert hashlib.sha256((case/'new_cache.zip').read_bytes()).hexdigest()==(
        'a5ee615cd46e93fa42a56e050c99274dbe46aeedb770ddad0b19106ffcaff67d')
    with zipfile.ZipFile(case/'new_cache.zip') as z:
        results=sorted(n for n in z.namelist() if n.endswith('/result.json'))
        assert len(results)==12
        rows=[json.loads(z.read(n)) for n in results]
        assert [r['coordinate_sha256'] for r in rows]==[
            r['coordinate_sha256'] for r in proof['new_trials']]
