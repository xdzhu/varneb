"""Synthetic cached replays; no calculator or material result generated."""
from pathlib import Path
import importlib.util
from types import SimpleNamespace
import sys
import numpy as np
import pytest
from vcneb.branch_relaxation import relax_branch_candidate

P=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/G3_preserving_branches_E260_20261017/prefix_review_E262.py'
spec=importlib.util.spec_from_file_location('E262_cached_replay_tests',P)
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)


def fixture(capped):
    records=[];basis=np.eye(3);free=basis[:,1:];direction=basis[:,2];seed=np.eye(3)
    point=dict(starting_internal_coordinates_A=[.17,.02,.12])
    def evaluate(x):
        if capped and len(records)==3:raise RuntimeError(w.MESSAGE)
        if capped:energy=x[1];gradient=np.array([100.,1.,0.])
        else:energy=.5*(x[1]**2+x[2]**2);gradient=np.array([100.,x[1],x[2]])
        records.append(dict(coordinates=x.tolist(),enthalpy_eV=float(energy),gradient=gradient.tolist()))
        return energy,gradient
    try:
        result=relax_branch_candidate(evaluate,np.asarray(point['starting_internal_coordinates_A']),
            free,direction,initial_hessian=seed,valid=lambda x:True)
    except RuntimeError as exc:
        assert str(exc)==w.MESSAGE;result=None
    return records,point,basis,free,direction,seed,result


def run(values):
    records,point,basis,free,direction,seed,_=values
    return w.replay(SimpleNamespace(relax_branch_candidate=relax_branch_candidate),
        point,basis,free,direction,seed,records,lambda x:True)


def test_cap_replay_recovers_last_accepted_state_without_new_call():
    values=fixture(True);result=run(values)
    assert result['status']=='budget_exception_replayed'
    assert result['recorded_SCFs_consumed']==3 and result['new_DFT']==0
    assert len(result['accepted_trace'])==3
    np.testing.assert_array_equal(result['last_accepted_internal_A'],values[0][-1]['coordinates'])
    assert result['last_accepted_energy_eV_cell']==values[0][-1]['enthalpy_eV']
    assert result['uncomputed_proposal_is_not_a_DFT_point']


def test_normal_candidate_replays_same_energy_and_stop_tests():
    values=fixture(False);result=run(values)
    assert result['status']=='normal_return_replayed'
    assert result['result']==values[-1]


def test_changed_ordered_coordinate_refuses_instead_of_searching_nearest():
    values=fixture(True);values[0][1]['coordinates'][1]+=.001
    old=sys.gettrace()
    with pytest.raises(ValueError,match='ordered recorded proposal'):run(values)
    assert sys.gettrace() is old


def test_extra_unconsumed_records_are_refused():
    values=fixture(False);values[0].append(dict(values[0][-1]))
    with pytest.raises(ValueError,match='every recorded SCF'):run(values)
