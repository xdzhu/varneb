"""Pure arithmetic/replay gates; synthetic records below are not DFT evidence."""
from pathlib import Path
from copy import deepcopy
import hashlib,importlib.util,json
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/G3_atomic_audit_E176_20261015'
OLD=ROOT/'benchmarks/hfo2_channels/20261008/G3_atomic_release_E175_20261015'
def load(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
f=load('atomic_expectation_tests',CASE/'forecast.py')
a=load('atomic_auditor_tests',CASE/'audit.py')
w=load('original_atomic_release_for_audit_tests',OLD/'worker.py')


def assert_formal_model_replay(actual, expected):
    """Only formal-model floating arithmetic; raw inputs and metadata stay exact."""
    assert type(actual) is type(expected)
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys()
        for key in expected:
            assert_formal_model_replay(actual[key], expected[key])
    elif isinstance(expected, list):
        assert len(actual) == len(expected)
        for a, e in zip(actual, expected):
            assert_formal_model_replay(a, e)
    elif isinstance(expected, float):
        assert actual == pytest.approx(expected, abs=1e-12, rel=0)
    else:
        assert actual == expected


@pytest.mark.parametrize('seed',range(6))
def test_signed_three_treatment_identity_and_common_energy_zero(seed):
    values=np.random.default_rng(seed).normal(size=6)
    out=f.signed_treatment_changes(*values[:3],origins=values[3:])
    again=f.signed_treatment_changes(*(values[:3]+1000),origins=values[3:]+1000)
    for key in ('absolute','shape'):
        assert out[key]['atomic_minus_frozen_eV']+out[key]['joint_minus_atomic_eV']==pytest.approx(out[key]['joint_minus_frozen_eV'])
        for label,v in out[key].items():assert again[key][label]==pytest.approx(v,abs=1e-11)


def test_positive_release_and_missing_origin_are_preserved_not_hidden_or_clipped():
    out=f.signed_treatment_changes(0.,.2,.3,origins=[0.,None,0.])
    assert out['absolute']['atomic_minus_frozen_eV']==pytest.approx(.2) and out['shape'] is None
    missing=f.signed_treatment_changes(0.,None,.3,origins=[0.,None,0.])
    assert missing['absolute'] is None
    with pytest.raises(ValueError,match='finite'):
        f.signed_treatment_changes(0.,float('nan'),.3,origins=[0.,0.,0.])


def test_all_twenty_unfitted_training_expectations_are_reproducible_and_not_holdout_results():
    out=f.forecast();saved=json.loads((CASE/'forecast.json').read_bytes())
    assert hashlib.sha256((CASE/'forecast.json').read_bytes()).hexdigest() == '321ba7007c7f950166aab4eaa2dbe579773164968f04a94034a8a942bee1318d'
    assert_formal_model_replay(out, saved)
    assert [(r['channel'],r['probe_step_A']) for r in out['rows']]==[
        ('T',.01),('T',.02),('preserving',.01),('preserving',.02)]
    assert sum(len(r['points']) for r in out['rows'])==20
    assert out['fitted_coefficients']==out['independent_slots_used']==out['new_DFT']==0
    assert not out['atomic_treatment_result_labels_read']
    assert not out['physical_reserved_geometry_generated_or_read']
    assert not out['post_relaxation_stability_or_physical_error_or_prediction_advantage_certified']
    assert [r['formal_atomic_origin_release_meV_cell'] for r in out['rows']]==pytest.approx([
        -.3893452356121704,-.3974376947445344,-50.38743850676271,-50.515678344441795])


@pytest.mark.parametrize('changed', [dict(energy=1e-9, source='same'), dict(energy=0., source='changed')])
def test_formal_replay_still_rejects_numeric_or_exact_metadata_change(changed):
    with pytest.raises(AssertionError):
        assert_formal_model_replay(changed, dict(energy=0., source='same'))


def test_partial_json_or_missing_point_is_pending_not_a_false_failure(tmp_path):
    c=json.loads((OLD/'config.json').read_bytes())
    rows,pending=a.load_prefix(tmp_path,c)
    assert not rows and pending==dict(channel='T',id=0,reason='point_not_present')
    p=tmp_path/'T/point_000.json';p.parent.mkdir();p.write_text('{')
    rows,pending=a.load_prefix(tmp_path,c)
    assert not rows and pending['reason']=='point_JSON_incomplete'
    p.write_text(json.dumps(dict(channel='T',id=0,q_A=[0.,0.])))
    rows,pending=a.load_prefix(tmp_path,c)
    assert len(rows)==1 and pending['id']==1
    p.write_text(json.dumps(dict(channel='T',id=False,q_A=[0.,0.])))
    with pytest.raises(ValueError,match='channel/id/Q'):a.load_prefix(tmp_path,c)


def test_exact_cache_order_and_twelve_attempt_denial_never_call_a_calculator():
    records=[dict(coordinates=[float(i)],enthalpy_eV=float(i),gradient=[1.]) for i in range(12)]
    oracle=a.CacheOnly(records)
    for i in range(12):assert oracle(np.array([float(i)]))[0]==float(i)
    assert oracle(np.array([0.]))[0]==0.
    with pytest.raises(RuntimeError,match='budget exhausted'):oracle(np.array([12.]))
    with pytest.raises(ValueError,match='differently ordered'):a.CacheOnly(records)(np.array([1.]))


@pytest.fixture(scope='module')
def synthetic_node():
    c=json.loads((OLD/'config.json').read_bytes());d=deepcopy(w.validate(c,OLD)['T'])
    p=c['channels']['T']['points'][0];r=np.array(c['channels']['T']['atomic_release_basis_36_by_31'])
    h=np.eye(36)*2;g=r@np.linspace(-.02,.02,31)
    d['hessian']=h;d['center_energy']=0.;d['center_gradient']=d['basis']@g
    records=[]
    def backend(x):
        for old in records:
            if np.array_equal(x,old['coordinates']):
                return old['enthalpy_eV'],np.array(old['gradient'])
        u=d['basis'].T@x;energy=float(g@u+.5*u@h@u);fullg=d['basis']@(g+h@u)
        records.append(dict(coordinates=x.tolist(),enthalpy_eV=energy,gradient=fullg.tolist()))
        return energy,fullg
    evaluate=w.CachedPullback(d,[],backend,p)
    result=w.relax_point(d,np.array(p['q_A']),r,evaluate,base=p['frozen_internal_A'])
    row=dict(candidate=result,reason=None,status='atomic_gradient_candidate_stability_unchecked',
        new_completed_SCFs=len(records),reserved_attempts=len(records),cached_callbacks=evaluate.hits,
        all_callbacks=evaluate.requests,candidate_native_identity=evaluate.last)
    return d,p,row,records,r


def test_actual_optimizer_replay_on_explicit_synthetic_quadratic(synthetic_node):
    d,p,row,records,r=synthetic_node
    out=a.replay_node(w,d,p,row,records,r,[])
    assert out['exact_request_order_optimizer_replay'] and out['new_DFT']==0


@pytest.mark.parametrize('key', ['status','candidate','cached_callbacks','all_callbacks','candidate_native_identity',
                                'new_completed_SCFs','reserved_attempts'])
def test_wrong_candidate_counters_or_fake_final_trial_cannot_pass(synthetic_node,key):
    d,p,original,records,r=synthetic_node;row=deepcopy(original)
    if key=='status':row[key]='converged_TS'
    elif key in ('candidate','candidate_native_identity'):row[key]=None
    else:row[key]+=1
    with pytest.raises(ValueError):a.replay_node(w,d,p,row,records,r,[])
