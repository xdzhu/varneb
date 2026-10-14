"""Fixed-Q atomic-only release: exact cell constraints, cache and finite policy."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import importlib.util,json
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/G3_atomic_release_E175_20261015'
spec=importlib.util.spec_from_file_location('atomic_release_case',CASE/'worker.py')
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)


@pytest.fixture(scope='module')
def prepared():
    config=w.prepare(CASE)
    return config,w.validate(config,CASE)


def test_frozen_config_file_is_exact_source_bound_not_only_a_regenerated_fixture():
    frozen=json.loads((CASE/'config.json').read_bytes())
    w.validate(frozen,CASE)


def test_existing_two_centers_atomic_constraint_and_all_five_coordinates(prepared):
    c,data=prepared
    assert c['policy']['training_committed_after']==dict(T=34,preserving=35)
    assert c['policy']['training_remaining_after']==dict(T=47,preserving=46)
    assert c['policy']['independent_slots_untouched_per_face']==8
    for channel,d in data.items():
        r=np.asarray(c['channels'][channel]['atomic_release_basis_36_by_31'])
        assert r.shape==(36,31) and not np.any(r[33:])
        np.testing.assert_allclose(d['axes'].T@r,0,atol=2e-12,rtol=0)
        np.testing.assert_allclose(r.T@r,np.eye(31),atol=2e-12,rtol=0)
        for p in c['channels'][channel]['points']:
            x=np.asarray(p['frozen_internal_A'])
            np.testing.assert_allclose(d['basis']@x,p['frozen_full_A'],atol=2e-12,rtol=0)
            y=x+r@np.linspace(-.01,.01,31)
            np.testing.assert_array_equal(y[33:],x[33:])
            np.testing.assert_array_equal(d['chart'].to_atoms(d['basis']@y).cell.array,
                                          d['chart'].to_atoms(d['basis']@x).cell.array)


@pytest.mark.parametrize('fault',['cutoff','budget','holdout','source','extra_point','cell_motion','nonfinite'])
def test_mutated_policy_sources_coordinates_or_release_space_fail_before_DFT(prepared,fault):
    c=deepcopy(prepared[0])
    if fault=='cutoff':c['six_input_sha256']['INPUT']='0'*64
    elif fault=='budget':c['policy']['maximum_new_attempts_per_node']=13
    elif fault=='holdout':c['physical_reserved_geometry_generated_or_read']=True
    elif fault=='source':c['worker_text_sha256']='0'*64
    elif fault=='extra_point':c['channels']['T']['points'].append(c['channels']['T']['points'][0])
    elif fault=='cell_motion':c['channels']['T']['atomic_release_basis_36_by_31'][33][0]=1e-16
    else:c['channels']['T']['atomic_release_basis_36_by_31'][0][0]=float('nan')
    with pytest.raises(ValueError):w.validate(c,CASE)


def test_config_does_not_alias_physical_contract_or_policy(prepared):
    c=w.prepare(CASE);c['policy']['channels'].reverse();c['six_input_sha256']['INPUT']='bad'
    assert w.POLICY['channels']==['T','preserving']
    assert w.study.CONTRACT['INPUT'].startswith('dc6684ff')


@pytest.mark.parametrize('channel',['T','preserving'])
def test_quadratic_atomic_release_solves_the_actual31D_gradient_with_exact_cell(prepared,channel):
    c,data=prepared;d=deepcopy(data[channel])
    r=np.asarray(c['channels'][channel]['atomic_release_basis_36_by_31'])
    h=np.eye(36)*2;g=r@np.linspace(-.015,.015,31);q=np.array([.1,0.])
    d['hessian']=h
    def evaluate(x):return float(g@x+.5*x@h@x),g+h@x
    out=w.relax_point(d,q,r,evaluate)
    x=np.asarray(out['coordinates_internal_A']);base=d['axes']@q
    np.testing.assert_array_equal(x[33:],base[33:])
    assert out['atomic_released_gradient_norm_eV_A']<=.03
    np.testing.assert_allclose(d['axes'].T@x,q,atol=2e-12,rtol=0)
    assert out['geometry_screens']['rejected']==0


def test_frozen_seed_and_zero_are_reused_and_wrong_first_callback_never_calls_backend(prepared):
    c,data=prepared;d=data['T'];controls=w.load_controls(CASE)['rows'];calls=[]
    def backend(x):calls.append(x);raise AssertionError('no DFT expected')
    oracle=w.CachedPullback(d,[p for p in controls if p['channel']=='T'],backend,c['channels']['T']['points'][0])
    e,g=oracle(np.zeros(36));assert e==d['center_energy'] and oracle.hits==1
    point=c['channels']['T']['points'][1]
    oracle=w.CachedPullback(d,[p for p in controls if p['channel']=='T'],backend,point)
    x=np.asarray(point['frozen_internal_A'])
    e,g=oracle(x);assert oracle.hits==1 and not calls
    oracle=w.CachedPullback(d,[],backend,point)
    with pytest.raises(ValueError,match='exact frozen seed'):oracle(x)
    assert not calls


def test_literal_native_anchor_avoids_cross_machine_rebuilding_of_frozen_geometry(prepared):
    c,data=prepared;d=data['T'];p=c['channels']['T']['points'][1]
    rows=[r for r in w.load_controls(CASE)['rows'] if r['channel']=='T']
    oracle=w.CachedPullback(d,rows,lambda x:pytest.fail('cached seed must not invoke DFT'),p)
    # A literal cached seed stays exact even if one basis element changes by1ULP.
    d=deepcopy(d);d['basis'][0,0]=np.nextafter(d['basis'][0,0],1.)
    oracle.data=d
    oracle(np.asarray(p['frozen_internal_A']))
    np.testing.assert_array_equal(oracle.last['coordinates_full_A'],p['frozen_full_A'])
    assert oracle.hits==1


def test_all_pending_and_held_study_allocations_count_toward_limit(monkeypatch):
    def queue(text):
        monkeypatch.setattr(w.subprocess,'run',lambda *a,**k:SimpleNamespace(stdout=text))
    queue('1|hfo2-G3-E161-preserving|RUNNING|32|1\n')
    assert w.scheduler_snapshot()['ready']
    queue('1|hfo2-G3-E161-preserving|RUNNING|32|1\n2|hfo2-G3-E175-atoms|PENDING|32|1\n')
    assert not w.scheduler_snapshot()['ready'] and w.scheduler_snapshot('2')['ready']
    queue('1|hfo2-G3-E161-preserving|RUNNING|32|1\n2|hfo2-G3-E175-atoms|PENDING|32|1\n3|hfo2-G3-other|PENDING|32|1\n')
    assert not w.scheduler_snapshot('2')['ready']
    queue('1|hfo2-G3-E175-atoms|RUNNING|64|1\n')
    with pytest.raises(ValueError,match='resources'):w.scheduler_snapshot('1')


def test_backend_failure_is_not_reclassified_as_an_optimizer_hole(prepared,tmp_path,monkeypatch):
    class Fake:
        def __init__(self,*args,workdir,**kwargs):
            self.workdir=workdir;workdir.mkdir(parents=True)
            self.n_new_evaluations=0;self.n_evaluation_attempts=0
        def __call__(self,x):
            self.n_evaluation_attempts+=1
            (self.workdir/'eval-test').mkdir()
            (self.workdir/'eval-test/failure.json').write_text('{}')
            raise RuntimeError('atomic released-gradient target not reached')
    monkeypatch.setattr(w.faces.canary,'CalculatorModeEvaluator',Fake)
    with pytest.raises(RuntimeError,match='atomic released-gradient'):
        w.compute(*prepared,tmp_path/'fresh',lambda *a:None)
    assert not list((tmp_path/'fresh').glob('*/point_*.json'))
