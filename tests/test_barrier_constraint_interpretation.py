"""Analytic interpretation fixture, never material DFT or a predictor."""
import numpy as np
import pytest


def potential(x,y):
    return (x*x-1)**2+(1+x*x)*y*y


def gradient(x,y):
    return np.array([4*x*(x*x-1)+2*x*y*y,2*(1+x*x)*y])


def hessian(x,y):
    return np.array([[12*x*x-4+2*y*y,4*x*y],[4*x*y,2*(1+x*x)]])


def test_added_constraint_raises_pointwise_energy_not_necessarily_barrier():
    x=np.linspace(-2,2,301)
    assert np.all(potential(x,1)>=potential(x,0))
    free_well=potential(1,0); free_peak=potential(0,0)
    constrained_well=potential(np.sqrt(.5),1); constrained_peak=potential(0,1)
    assert (free_well,free_peak)==(0,1)
    assert constrained_well==pytest.approx(7/4)
    assert constrained_peak==2
    assert constrained_peak-constrained_well==pytest.approx(1/4)
    assert constrained_peak>free_peak
    assert constrained_peak-constrained_well<free_peak-free_well


def test_stationary_points_and_prescribed_reactions_are_distinct():
    for x in (-1,1):
        np.testing.assert_allclose(gradient(x,0),0,atol=1e-15)
        assert np.linalg.eigvalsh(hessian(x,0)).min()>0
    assert np.count_nonzero(np.linalg.eigvalsh(hessian(0,0))<0)==1
    for x in (-np.sqrt(.5),np.sqrt(.5)):
        assert gradient(x,1)[0]==pytest.approx(0,abs=1e-15)
        assert hessian(x,1)[0,0]>0
        assert gradient(x,1)[1]==pytest.approx(3)
    assert gradient(0,1)[0]==0 and gradient(0,1)[1]==2
    assert hessian(0,1)[0,0]==-2
    assert not np.array_equal(np.array([1,0]),np.array([np.sqrt(.5),1]))


@pytest.mark.parametrize('fixed_y',[0.,.25,.5,1.,1.25])
def test_fixed_y_family_and_stationary_envelope_derivative(fixed_y):
    rho=fixed_y
    x=np.sqrt(1-rho*rho/2)
    well=potential(x,rho); peak=potential(0,rho)
    assert gradient(x,rho)[0]==pytest.approx(0,abs=2e-15)
    assert well==pytest.approx(2*rho*rho-rho**4/4)
    assert peak==pytest.approx(1+rho*rho)
    assert peak-well==pytest.approx((1-rho*rho/2)**2)
    assert hessian(x,rho)[0,0]>0 and hessian(0,rho)[0,0]<0
    difference=gradient(0,rho)[1]-gradient(x,rho)[1]
    assert difference==pytest.approx(-2*rho+rho**3)
    q=np.linspace(-x,x,501)
    assert np.max(potential(q,rho))==pytest.approx(peak)


@pytest.mark.parametrize('point',[(-.7,.6),(.3,1.),(1.2,.1)])
def test_analytic_gradient_and_hessian_match_numeric_checks(point):
    x,y=point; delta=1e-5
    numeric=np.array([(potential(x+delta,y)-potential(x-delta,y))/(2*delta),
                      (potential(x,y+delta)-potential(x,y-delta))/(2*delta)])
    np.testing.assert_allclose(numeric,gradient(x,y),atol=2e-9,rtol=0)
    numeric_h=np.column_stack([(gradient(x+delta,y)-gradient(x-delta,y))/(2*delta),
                              (gradient(x,y+delta)-gradient(x,y-delta))/(2*delta)])
    np.testing.assert_allclose(numeric_h,hessian(x,y),atol=2e-9,rtol=0)
