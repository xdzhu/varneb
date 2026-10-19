"""Actual E272 source replay and bounded numerical continuation safeguards."""
import importlib.util
import json
from pathlib import Path
import numpy as np
import pytest

BASE=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008'
CASE=BASE/'preserving_local_completion_E276_20261019'


def load(path,name):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


w=load(CASE/'worker.py','E276_test_worker')
study=load(BASE/'matched_pattern_release_E270_20261018/study.py','E276_test_native')
helper=load(BASE/'matched_pattern_seeded_E271_20261018/repair.py','E276_test_space')


@pytest.fixture(scope='module')
def prepared():
    seed=CASE/'prepared/seed';data=json.loads((seed/'seed.json').read_bytes())
    h=json.loads((seed/'branch_preparation.json').read_bytes())
    return study.chart_for(seed),data,np.asarray(h['numerical_hessian_39']),np.asarray(h['signed_starts'][1])/.1


@pytest.mark.parametrize('sign',['minus','plus'])
def test_actual_failed_start_replays_unchanged_native_EFS_and_original_controls(prepared,sign):
    chart,data,_,_=prepared;x,e,g=w.cached_start(chart,CASE/'prepared/starts'/sign,data)
    assert .4999<np.linalg.norm(x)<.5 and not np.any(x[36:])
    q=np.asarray(data['orthonormal_39_by_4']);t=np.zeros((39,3));t[:36]=np.tile(np.eye(3),(12,1))/np.sqrt(12)
    rg=g-q@(q.T@g)-t@(t.T@g);rg[36:]=0
    assert .37<np.linalg.norm(rg)<.39
    assert np.max(np.abs(np.asarray(data['covectors_39_by_4']).T@x))<1e-12
    assert e<json.loads((CASE/'prepared/config.json').read_bytes())['original_F_energy_eV_cell']


def test_diagnosis_does_not_convert_boundary_crawling_to_SCF_failure_or_convergence():
    d=json.loads((CASE/'original_failure_diagnosis/diagnosis.json').read_bytes())
    assert d['local_native_calls_reparsed']==24 and d['new_DFT_calls']==0
    for a in d['arms']:
        last=a['trace'][-1]
        assert len(a['trace'])==12 and last['tangent_to_sphere_gradient_eV_A']>.37
        assert last['step_from_previous_trial_A']<5e-5
        assert last['released_gradient_eV_A']>.03 and a['final_is_not_a_qualified_candidate']


def test_existing_band_scale_requires_one_registered_lift_not_peakwise_MIC():
    d=json.loads((CASE/'registered_connection_scale.json').read_bytes())
    bad=json.loads((CASE/'connection_scale.json').read_bytes())
    assert not bad['same_atom_order_frame_and_continuous_lift']
    assert d['registration_fixed_by_initial_and_verified_at_final_endpoint']
    assert d['projected_same_observable_atomic_norm_A']==pytest.approx(.81006375896916)
    assert d['projected_energy_or_gradient'] is None and d['new_DFT_calls']==0


def test_cache_reuses_exact_points_and_enforces_stage_and_total_new_call_caps():
    calls=[]
    def evaluate(x):calls.append(x.copy());return float(x@x),2*x
    c=w.Cache(evaluate,[(np.zeros(2),0.,np.zeros(2))]);c.begin('A',2)
    assert c(np.zeros(2))[0]==0 and c.fresh==0
    c(np.ones(2));c(2*np.ones(2))
    with pytest.raises(w.AttemptCap):c(3*np.ones(2))
    assert c.fresh==2 and len(calls)==2
    c.begin('J',100)
    for i in range(3,51):c(i*np.ones(2))
    with pytest.raises(w.AttemptCap):c(51*np.ones(2))
    assert c.fresh==50 and len(calls)==50


@pytest.mark.parametrize('sign',[-1.,1.])
def test_same_BFGS_can_reach_a_model_minimum_beyond_old_guard_without_changing_controls(prepared,sign):
    # Synthetic numerical regression only, not a hafnia energy or speed result.
    chart,data,h,v=prepared;target=sign*.81*v;start=sign*.49997*v
    def model(x):
        dx=x-target;return float(.5*dx@h@dx),h@dx
    e,g=model(start);cache=w.Cache(model,[(start,e,g)])
    a=w.arm(study,helper,chart,data,cache,h,start,'A')
    assert a['candidate'] is not None and a['new_calls']<=32
    x=np.asarray(a['candidate']['coordinates'])
    assert .75<np.linalg.norm(x)<1. and not np.any(x[36:])
    assert a['candidate']['released_gradient_norm_eV_A']<=.03
    assert np.max(np.abs(np.asarray(data['covectors_39_by_4']).T@x))<1e-10
    assert not a['candidate']['stationary_TS_or_full_stability_certified']


def test_missing_A_blocks_checks_and_J_without_any_call(prepared):
    chart,data,h,v=prepared
    c=w.Cache(lambda x:pytest.fail('no calculation without A'),[])
    assert w.arm(study,helper,chart,data,c,h,None,'J')['new_calls']==0
    assert w.check_direction(study,chart,c,None,v)['new_calls']==0


def test_policy_and_submission_script_keep_physics_resources_and_explicit_revision():
    d=json.loads((CASE/'prepared/config.json').read_bytes());p=d['policy']
    assert p['original_norm_guard_A']==.5 and p['maximum_coordinate_norm_A']==1.
    assert p['atomic_attempts']+p['joint_attempts']+p['direction_probes']==50
    assert not p['physical_inputs_changed'] and not p['original_protocol_success_claimed']
    assert not d['execution_admitted'] and d['new_DFT_calls']==0
    script=(CASE/'run.slurm').read_text()
    for value in ('--cpus-per-task=32','--ntasks=1','--time=04:00:00','I_MPI_HYDRA_BOOTSTRAP=fork','RUN_DFT','ARM_SIGN'):
        assert value in script
    assert '235' not in script and 'srun ' not in script and 'sbatch ' not in script


def test_native_config_contains_geometry_hash_but_only_six_assets_are_fixed():
    import copy
    d=json.loads((CASE/'prepared/config.json').read_bytes())
    assert 'STRU' in d['physical_input_sha256']
    w.validate_config(d)
    for name in w.CONTRACT:
        bad=copy.deepcopy(d);bad['physical_input_sha256'][name]='changed'
        with pytest.raises(ValueError,match='fixed physical contract'):w.validate_config(bad)
    changed_geometry=copy.deepcopy(d);changed_geometry['physical_input_sha256']['STRU']='different geometry'
    # Geometry is checked against its native source separately, not as one of six assets.
    w.validate_config(changed_geometry)


@pytest.mark.parametrize('state',['RUNNING','PENDING','TIMEOUT','FAILED'])
def test_terminal_collector_refuses_nonterminal_or_failed_original_handle(state):
    collector=load(CASE/'collect_terminal.py','E276_collector_test_'+state)
    row='28954771|'+state+'|0:0|32|hfacnormal01|iai806|1|60|start|end|hfo2-E276-minus|node117'
    with pytest.raises(ValueError,match='exact original completed'):collector.terminal_row(row,'minus')


def test_terminal_reader_does_not_impute_a_still_missing_original_archive(tmp_path):
    import shutil
    # Keep the real collector code but remove the actual export namespace.
    # This checks missing data, not a missing Python module, after completion too.
    reader=load(CASE/'analyze_terminal.py','E276_offline_missing_test')
    missing=tmp_path/'missing_original_exports';missing.mkdir()
    shutil.copyfile(CASE/'collect_terminal.py',missing/'collect_terminal.py')
    reader.CASE=missing
    with pytest.raises(FileNotFoundError,match='collection_receipt.json'):reader.read_one('minus',tmp_path)
