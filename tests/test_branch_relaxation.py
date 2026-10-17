"""Analytic optimizer checks, not material evidence."""
import numpy as np
import pytest
from vcneb.branch_relaxation import relax_branch_candidate


def case(sign=1, rotation=None, guard=lambda x: True, cap=None):
    rotation = np.eye(3) if rotation is None else rotation
    start = rotation @ np.array([.17,.02,sign*.12])
    basis = rotation[:,1:]
    direction = rotation[:,2]
    hessian = rotation @ np.diag([1.,2.,.05]) @ rotation.T
    seen=[]
    def evaluate(x):
        if cap is not None and len(seen)>=cap:
            raise RuntimeError('native attempt cap')
        assert guard(x)
        seen.append(x.copy())
        q,r,z=rotation.T@x
        return r*r-.01*z*z+.5*z**4,rotation@np.array([100.,2*r,-.02*z+2*z**3])
    result=relax_branch_candidate(evaluate,start,basis,direction,initial_hessian=hessian,
        valid=guard,gradient_tolerance=.03,branch_gradient_tolerance=1e-4)
    return result,seen


@pytest.mark.parametrize('sign',[-1,1])
def test_soft_direction_must_converge_even_when_total_norm_would_pass(sign):
    result,seen=case(sign)
    assert result['status']=='gradient_candidate_stability_unchecked'
    assert len(seen)>1
    assert abs(result['branch_gradient'])<=1e-4
    assert result['released_gradient_norm']<=.03
    assert sign*result['coordinates'][2] == pytest.approx(.1,abs=.003)
    assert all(x[0]==.17 for x in seen)
    assert result['stability_or_TS_certified'] is False


def test_rotation_keeps_physical_solution_and_energy():
    rotation=np.linalg.qr(np.array([[1.,2.,3.],[-2.,3.,1.],[3.,1.,-2.]]))[0]
    base,_=case();rotated,seen=case(rotation=rotation)
    np.testing.assert_allclose(rotation.T@rotated['coordinates'],base['coordinates'],atol=1e-10,rtol=0)
    assert rotated['energy']==pytest.approx(base['energy'],abs=1e-13)
    assert all(abs(rotation[:,0]@x-.17)<1e-12 for x in seen)


def test_geometry_rejection_precedes_every_calculator_call():
    result,seen=case(guard=lambda x: abs(x[1])<=.02000001 and abs(x[2])<=.121)
    assert seen and all(abs(x[2])<=.121 for x in seen)
    assert result['released_gradient_norm']<=.03


def test_calculator_budget_exception_is_not_retried_or_swallowed():
    with pytest.raises(RuntimeError,match='native attempt cap'):
        case(cap=1)


@pytest.mark.parametrize('defect',['direction','negative_seed','complex','bad_start'])
def test_bad_contract_refuses_before_calculation(defect):
    start=np.array([.17,.02,.12]);basis=np.eye(3)[:,1:];direction=np.array([0.,0.,1.]);hessian=np.eye(3)
    if defect=='direction': direction=np.array([1.,0.,0.])
    if defect=='negative_seed': hessian[2,2]=-1
    if defect=='complex': start=start.astype(complex)+1j
    calls=[]
    with pytest.raises((ValueError,np.linalg.LinAlgError)):
        relax_branch_candidate(lambda x:calls.append(x),start,basis,direction,
            initial_hessian=hessian,valid=lambda x:defect!='bad_start')
    assert not calls
