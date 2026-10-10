"""Independent scalar-energy finite differences and rotation/convention checks."""
import numpy as np
import pytest
from vcneb.strain_work import configuration_work


@pytest.mark.parametrize('pressure',[0.,.037])
@pytest.mark.parametrize('seed',[1,2,3,4])
def test_configuration_work_matches_independent_energy_finite_difference(seed,pressure):
    rng=np.random.default_rng(seed)
    h=np.diag([3.,4.,5.])+rng.normal(scale=.1,size=(3,3))
    s=rng.normal(size=(4,3))*.2
    dh=rng.normal(size=(3,3))*.2;ds=rng.normal(size=(4,3))*.1
    kappa=.04;v0=58.;spring=.7
    def energy(h,s):
        v=np.linalg.det(h);r=s@h
        return .5*kappa*(v-v0)**2+.5*spring*np.sum(r*r)+pressure*v
    v=np.linalg.det(h);r=s@h
    sigma=kappa*(v-v0)*np.eye(3)+spring*(r.T@r)/v
    actual=configuration_work(h,sigma,dh,forces=-spring*r,fractional_direction=ds,pressure_eV_A3=pressure)
    step=1e-6
    expected=(energy(h+step*dh,s+step*ds)-energy(h-step*dh,s-step*ds))/(2*step)
    assert actual['total_eV_per_control']==pytest.approx(expected,rel=1e-7,abs=1e-7)


def test_rotated_oblique_row_cell_preserves_scalar_work_and_uses_no_factor_two():
    h=np.array([[3.,.3,.2],[.1,4.,.4],[.5,.2,5.]])
    sigma=np.array([[.1,.03,.02],[.03,-.02,.04],[.02,.04,.07]])
    dh=np.array([[.3,.1,.02],[.05,.2,.04],[.07,.03,.5]])
    theta=.63;r=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
    a=configuration_work(h,sigma,dh)
    b=configuration_work(h@r.T,r@sigma@r.T,dh@r.T)
    assert a['total_eV_per_control']==pytest.approx(b['total_eV_per_control'],abs=1e-12)
    shear=np.zeros((3,3));shear[0,1]=.2
    h=np.eye(3)
    assert configuration_work(h,sigma,shear)['cell_eV_per_control']==pytest.approx(.006)


def test_orthogonal_biaxial_scaling_and_rotation_work():
    strain=.01;h=np.diag([3*(1+strain),4*(1+strain),5.])
    sigma=np.diag([.1,.2,.3]);dh=np.zeros((3,3));dh[:2]=h[:2]/(1+strain)
    assert configuration_work(h,sigma,dh)['cell_eV_per_control']==pytest.approx(np.linalg.det(h)*.3/(1+strain))
    omega=np.array([[0.,-.3,.2],[.3,0.,-.4],[-.2,.4,0.]])
    assert configuration_work(h,sigma,h@omega.T)['total_eV_per_control']==pytest.approx(0.,abs=1e-12)


@pytest.mark.parametrize('bad',[np.zeros((3,3)),np.diag([-1.,1.,1.]),np.eye(2),np.full((3,3),np.nan)])
def test_invalid_cells_are_rejected(bad):
    with pytest.raises(ValueError):configuration_work(bad,np.eye(3),np.eye(3))


def test_wrong_stress_pressure_or_partial_atomic_arguments_are_rejected():
    for kwargs in ({'forces':np.ones((2,3))},{'pressure_eV_A3':True},{'pressure_eV_A3':np.nan},
                   {'pressure_eV_A3':[.1,.2]},{'pressure_eV_A3':'.1'},
                   {'forces':np.ones((2,3)),'fractional_direction':np.ones((1,3))}):
        with pytest.raises(ValueError):configuration_work(np.eye(3),np.eye(3),np.eye(3),**kwargs)
    with pytest.raises(ValueError):configuration_work(np.eye(3),np.arange(9).reshape(3,3),np.eye(3))
