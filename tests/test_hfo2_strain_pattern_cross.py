"""The final cross design preserves coordinates, inputs and its finite budget."""
import importlib.util
import json
from pathlib import Path
import numpy as np
import pytest

REPO=Path(__file__).resolve().parents[1]
BASE=REPO/'benchmarks/hfo2_channels/20261008'
CASE=BASE/'strain_pattern_cross_E277_20261021'


def module(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


@pytest.fixture(scope='module')
def modules():
    return (module(CASE/'prepare.py','cross_prepare'),module(CASE/'worker.py','cross_worker'),
        module(BASE/'matched_pattern_release_E270_20261018/study.py','cross_chart'),
        module(BASE/'matched_pattern_seeded_E271_20261018/repair.py','cross_helper'))


@pytest.mark.parametrize('name',['zero_q1','plus1_q0'])
def test_cross_geometry_and_unchanged_inputs(modules,name):
    prep,worker,study,helper=modules
    folder=CASE/'prepared'/name;chart=study.chart_for(folder/'seed')
    condition=json.loads((folder/'condition.json').read_bytes())
    data=json.loads((folder/'seed/seed.json').read_bytes())
    study.cached_frozen(folder/'seed',chart)
    x=np.asarray(condition['start_coordinates']);atoms=chart.to_atoms(x)
    q,_,_=prep.observables(atoms)
    np.testing.assert_allclose(q,condition['held_values_A'],atol=1e-12,rtol=0)
    assert worker.allowed(chart,x)
    np.testing.assert_allclose(atoms.cell.array[:2],chart.reference_cell[:2],atol=1e-12,rtol=0)
    plane,frozen,c,t=helper.release_space(data,False,x)
    np.testing.assert_allclose(c.T@x,condition['constraint_offset_A'],atol=1e-12,rtol=0)
    np.testing.assert_allclose(t.T@x,0,atol=1e-12,rtol=0)
    assert worker._orthogonal_directions(plane,frozen).shape==(39,32)
    for key,h in study.CONTRACT.items():assert data['physical_input_sha256'][key]==h


def test_exact_cache_and_failed_attempt_budget(modules):
    worker=modules[1];calls=[]
    def f(x):calls.append(x.copy());return float(x@x),2*x
    cache=worker.BoundedCache(f)
    for i in range(40):cache(np.array([float(i)]))
    cache(np.array([39.]));assert len(calls)==40
    with pytest.raises(worker.AttemptCap):cache(np.array([40.]))
    broken=worker.BoundedCache(lambda x:(_ for _ in ()).throw(ValueError('SCF failed')))
    with pytest.raises(ValueError):broken(np.zeros(1))
    assert broken.fresh==1


@pytest.mark.parametrize('name',['zero_q1','plus1_q0'])
def test_joint_release_keeps_nonzero_cross_constraints(modules,name):
    _,worker,study,helper=modules
    folder=CASE/'prepared'/name;chart=study.chart_for(folder/'seed')
    data=json.loads((folder/'seed/seed.json').read_bytes())
    condition=json.loads((folder/'condition.json').read_bytes());start=np.asarray(condition['start_coordinates'])
    plane,frozen,_,_=helper.release_space(data,False,start)
    u=worker._orthogonal_directions(plane,frozen);target=start+.015*u[:,0]
    # Identity is only an analytic test initializer, never a production override.
    condition['numerical_hessian_39']=np.eye(39).tolist()
    result=worker.release(chart,data,condition,helper,lambda x:(float((x-target)@(x-target)),2*(x-target)))
    assert result['status']=='gradient_candidate_stability_unchecked'
    np.testing.assert_allclose(result['candidate']['constraint_offset_A'],condition['constraint_offset_A'],atol=1e-10,rtol=0)
    assert result['new_calls']<=40


def test_only_two_conditions_and_no_automatic_extension(modules):
    prep=modules[0];cfg=json.loads((CASE/'prepared/config.json').read_bytes())
    assert cfg['policy']==prep.POLICY
    assert prep.POLICY['conditions']==['zero_q1','plus1_q0']
    assert prep.POLICY['maximum_new_calls']==80
    assert prep.POLICY['sign']=='minus'
    assert not prep.POLICY['physical_inputs_changed']
    assert not prep.POLICY['automatic_extension']
