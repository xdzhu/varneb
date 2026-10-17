"""Original sources, fixed training starts and unseen-data separation."""
from copy import deepcopy
import importlib.util
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/G3_preserving_branches_E260_20261017'


@pytest.fixture(scope='module')
def setup():
    spec=importlib.util.spec_from_file_location('E260_test_worker',CASE/'worker.py')
    w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
    return w,w.prepare(CASE)


def test_ten_literal_two_sign_starts_keep_old_Q_and_substrate(setup):
    w,c=setup;data,anchor,direction,*_=w.inputs(CASE)
    assert len(c['points'])==10
    for point in c['points']:
        x=np.asarray(point['starting_full_coordinates_A'])
        np.testing.assert_allclose(x,anchor+data['basis']@data['axes']@point['q_A']+
            point['zeta_start_A']*direction,atol=0,rtol=0)
        w.geometry(data,x,point['q_A'])
    assert [p['start_sign'] for p in c['points']]==[1,-1]*5


def test_budget_and_test_inventory_are_explicit(setup):
    _,c=setup;p=c['policy']
    assert p['maximum_attempts_cohort']==10*p['maximum_attempts_per_candidate']==120
    assert p['training_committed']==dict(T=34,preserving=49)
    assert p['training_before']==dict(T=34,preserving=39)
    assert p['T_independent_seen']==8 and p['preserving_independent_seen']==0
    assert p['independent_inventory']==16 and not p['physical_plus005_generated_or_read']
    assert sorted(sum(p['shards'],[]))==list(range(10))


def test_seed_is_positive_numerical_preconditioner_not_clipped_physics(setup):
    w,c=setup;data,_,_,free,h,d=w.inputs(CASE)
    assert np.linalg.eigvalsh(free.T@h@free)[0]>0
    assert d['outer_numerical_stiffness_eV_A2']==pytest.approx(.04986124402846963,abs=2e-12)
    assert d['raw_R33_min_eV_A2']>.38
    assert d['not_a_physical_Hessian_or_stability_certificate']
    assert c['original_2D_refusal_preserved']
    assert c['policy']['released_gradient_tolerance_eV_A']==.03
    assert c['policy']['zeta_gradient_tolerance_eV_A']==1e-4


@pytest.mark.parametrize('defect',['start','tolerance','budget','inventory'])
def test_frozen_protocol_mutations_refuse(setup,defect):
    w,base=setup;c=deepcopy(base)
    if defect=='start':c['points'][0]['zeta_start_A']=.2
    if defect=='tolerance':c['policy']['zeta_gradient_tolerance_eV_A']=.01
    if defect=='budget':c['policy']['maximum_attempts_cohort']=121
    if defect=='inventory':c['policy']['preserving_independent_seen']=8
    with pytest.raises(ValueError,match='frozen training'):w.validate(c,CASE)


def test_local_shell_cannot_launch_DFT(setup,monkeypatch):
    w,c=setup;monkeypatch.delenv('SLURM_JOB_ID',raising=False)
    with pytest.raises(ValueError,match='authorized'):w.run(c,CASE,0)


def test_arithmetic_replay_keeps_frozen_optimizer_seed(setup):
    w,base=setup;c=deepcopy(base)
    c['preconditioner_36_by_36'][0][0]+=1e-14
    values=w.validate(c,CASE)
    np.testing.assert_array_equal(values[4],c['preconditioner_36_by_36'])
    changed=deepcopy(base);changed['points'][0]['starting_internal_coordinates_A'][0]+=.001
    with pytest.raises(ValueError,match='frozen training'):w.validate(changed,CASE)
