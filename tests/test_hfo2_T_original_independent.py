"""Original test identity and per-center admission, not synthetic DFT evidence."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/G3_T_independent_E258_20261017'


@pytest.fixture(scope='module')
def setup():
    spec=importlib.util.spec_from_file_location('original_T_test_worker',CASE/'worker.py')
    w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
    c=w.prepare(CASE)
    return w,c


def test_every_original_target_start_and_prediction_is_retained(setup):
    w,c=setup
    f=json.loads((CASE/'standard_forecast.json').read_bytes())
    row=next(r for r in f['channels'] if r['channel']=='T')
    assert c['points']==row['targets'] and c['standard_predictions']==row
    assert c['original_scoring']==f['scoring']
    assert len(c['points'])==8 and all(p['independent_DFT_label'] is None for p in c['points'])
    assert {p['training_parent_id'] for p in c['points']}=={0}
    assert c['reference_predictions']==[r for r in json.loads((CASE/'reference_forecast.json').read_bytes())['rows'] if r['channel']=='T']


def test_both_T_steps_pass_without_repairing_preserving(setup):
    _,c=setup
    assert c['raw_T_minima_eV_A2']==pytest.approx([.1640645284775998,.1575959364082007],abs=2e-12)
    assert c['policy']['preserving_admission_refused']
    full=json.loads((CASE/'complete_native.json').read_bytes())
    p=next(r for r in full['channels'] if r['channel']=='preserving')
    assert all(np.linalg.eigvalsh(.5*(np.asarray(b['raw_projected_eV_A2'])+np.asarray(b['raw_projected_eV_A2']).T))[0]<0 for b in p['blocks'])


@pytest.mark.parametrize('change',['target','start','budget','preserving','prediction'])
def test_changed_frozen_test_contract_is_refused(setup,change):
    w,base=setup;c=deepcopy(base)
    if change=='target':c['points'][0]['q_A'][0]+=.001
    if change=='start':c['points'][0]['starting_internal_coordinates_A'][0]+=.001
    if change=='budget':c['policy']['maximum_new_attempts_per_point']=13
    if change=='preserving':c['policy']['channel']='preserving'
    if change=='prediction':c['reference_predictions']=c['reference_predictions'][:-1]
    with pytest.raises(ValueError,match='frozen original'):w.validate(c,CASE)


def test_shards_are_disjoint_and_do_not_expand_original_attempt_or_training_budget(setup):
    _,c=setup;p=c['policy']
    assert p['shards']==[[1,2,3,4],[5,6,7,8]]
    assert p['maximum_new_attempts_this_cohort']==8*p['maximum_new_attempts_per_point']==96
    assert p['original_maximum_attempts_two_faces']==192
    assert p['training_committed']==dict(T=34,preserving=39)
    assert p['automatic_retries']==0 and not p['physical_plus005_generated_or_read']


def test_optimizer_gets_literal_start_without_extra_start_or_new_label_parent(setup,monkeypatch):
    w,c=setup;data=w.faces.load_inputs('T');seen=[]
    def optimizer(plane,q,evaluate,**kw):
        seen.append((q,kw));return SimpleNamespace()
    monkeypatch.setattr(w.faces.canary,'relax_orthogonal_at_q',optimizer)
    for p in c['points'][:2]:w.optimize_point(data,p,None,lambda x:True)
    for point,(q,kw) in zip(c['points'],seen):
        assert q==point['q_A'] and kw['include_frozen_start'] is False
        np.testing.assert_array_equal(kw['starts'][0],point['starting_internal_coordinates_A'])
        np.testing.assert_array_equal(kw['initial_hessian'],data['hessian'])
        assert len(kw['starts'])==1 and kw['initial_trust_radius']==.05
        assert kw['gradient_tolerance']==.03 and kw['max_iterations']==200


def test_run_refuses_local_shell_before_any_calculator_or_output(setup,monkeypatch):
    w,c=setup;monkeypatch.delenv('SLURM_JOB_ID',raising=False)
    with pytest.raises(ValueError,match='authorized'):w.run(c,CASE,0)
