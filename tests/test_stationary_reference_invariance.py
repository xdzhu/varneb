"""A complete harmonic release cannot gain accuracy from a basis rename."""
import numpy as np
import pytest

from vcneb.stationary_branch import stationary_quadratic_branch


@pytest.mark.parametrize('rotated',[False,True])
@pytest.mark.parametrize('rank',[1,3,4])
def test_complete_internal_release_matches_basis_independent_fixed_control_solve(rotated,rank):
    rng=np.random.default_rng(151)
    rotation=np.linalg.qr(rng.normal(size=(7,7)))[0] if rotated else np.eye(7)
    internal=rotation@np.diag([-2.,3.,4.,5.,6.,7.,8.])@rotation.T
    h=np.zeros((8,8));h[:7,:7]=internal;h[:7,7]=h[7,:7]=rng.normal(size=7)*.2;h[7,7]=-.3
    g=rng.normal(size=8)*.04
    local=np.eye(7);angle=.1*rank
    local[0,0]=local[6,6]=np.cos(angle);local[6,0]=np.sin(angle);local[0,6]=-np.sin(angle)
    basis=np.pad(rotation@local,((0,1),(0,0)))
    allowed=np.eye(8)[:,:7];control=np.eye(8)[:,7]
    branch=stationary_quadratic_branch(h,g,-10000.,basis[:,:rank],basis[:,rank:],control,
        admissible_internal_basis=allowed,expected_index=1,stability_floor=.01,internal_curvature_floor=.01)
    for t in (-.03,0.,.06):
        point=branch.evaluate(t)
        direct=np.r_[-np.linalg.solve(internal,g[:7]+h[:7,7]*t),t]
        np.testing.assert_allclose(point.full_displacement,direct,atol=1e-13,rtol=0)
        assert point.energy_change==pytest.approx(float(g@direct+.5*direct@h@direct),abs=1e-13)
        assert point.unrepresented_internal_gradient_norm<1e-13
        assert point.control_gradient==pytest.approx(float((g+h@direct)[7]),abs=1e-13)
    exact=float(h[7,7]-h[7,:7]@np.linalg.solve(internal,h[:7,7]))
    assert branch.control_curvature==pytest.approx(exact,abs=1e-13)
