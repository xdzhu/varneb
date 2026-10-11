"""Physical covector work identities; not a harmonic accuracy benchmark."""
import numpy as np
import pytest
from vcneb.path_work import sampled_path_work


@pytest.mark.parametrize('seed',range(8))
def test_quadratic_work_exact_for_every_polygon_segment_and_subspace_rotation(seed):
    rng=np.random.default_rng(seed)
    x=rng.normal(size=(11,7))
    h=rng.normal(size=(7,7));h=.5*(h+h.T)
    offset=rng.normal(size=7)
    g=x@h+offset
    energy=.5*np.einsum('ni,ni->n',x@h,x)+x@offset
    basis=np.linalg.qr(rng.normal(size=(7,7)))[0]
    parts=dict(retained=basis[:,:3],omitted=basis[:,3:6],cell=basis[:,6:])
    result=sampled_path_work(x,g,energy,parts)
    np.testing.assert_allclose(result['segment_full_work'],np.diff(energy),atol=1e-13,rtol=0)
    np.testing.assert_allclose(result['cumulative_full_work'],energy-energy[0],atol=1e-13,rtol=0)
    q=np.linalg.qr(rng.normal(size=(3,3)))[0]
    rotated=dict(retained=parts['retained']@q,omitted=parts['omitted']@q,cell=parts['cell'])
    other=sampled_path_work(x,g,energy,rotated)
    for name in parts:
        np.testing.assert_allclose(result['segment_partition_work'][name],other['segment_partition_work'][name],atol=1e-13,rtol=0)
    # A simultaneous orthogonal frame change is a dual coordinate transform.
    frame=np.linalg.qr(rng.normal(size=(7,7)))[0]
    transformed=sampled_path_work(x@frame,g@frame,energy,{n:frame.T@b for n,b in parts.items()})
    for name in parts:
        np.testing.assert_allclose(result['segment_partition_work'][name],transformed['segment_partition_work'][name],atol=1e-13,rtol=0)
    assert result['physical_error_bound_measured'] is False


def test_signed_work_not_clipped_and_empty_subspace_allowed():
    result=sampled_path_work([[0.,0.],[1.,2.]],[[1.,-2.],[1.,-2.]],[4.,1.],
        dict(up=np.eye(2)[:,:1],down=np.eye(2)[:,1:],empty=np.zeros((2,0))))
    assert result['segment_partition_work']==dict(up=[1.],down=[-4.],empty=[0.])
    assert result['segment_full_work']==[-3.]
    assert result['endpoint_energy_minus_work']==0.


def test_nonquadratic_coarse_error_is_explicit_and_refinement_changes_it():
    errors=[]
    for n in (3,5,9):
        x=np.linspace(0.,1.,n)[:,None]
        result=sampled_path_work(x,4*x**3,x[:,0]**4,dict(all=np.eye(1)))
        errors.append(abs(result['endpoint_energy_minus_work']))
        assert not result['physical_error_bound_measured']
    np.testing.assert_allclose(errors,[.25,.0625,.015625],atol=1e-13,rtol=0)


@pytest.mark.parametrize('bad', ['one','wrong_gradient','wrong_energy','nonfinite_x','nonfinite_g','nonfinite_E',
                                'no_partition','missing','overlap','nonorthogonal','nonfinite_basis','complex_x','complex_basis','bad_name'])
def test_malformed_or_incomplete_data_is_not_silently_repaired(bad):
    x=np.array([[0.,0.],[1.,2.]])
    g=x.copy();energy=np.array([0.,2.5]);parts=dict(all=np.eye(2))
    if bad=='one':x=x[:1];g=g[:1];energy=energy[:1]
    elif bad=='wrong_gradient':g=g[:1]
    elif bad=='wrong_energy':energy=energy[:1]
    elif bad=='nonfinite_x':x[0,0]=np.nan
    elif bad=='nonfinite_g':g[0,0]=np.inf
    elif bad=='nonfinite_E':energy[0]=np.nan
    elif bad=='no_partition':parts={}
    elif bad=='missing':parts=dict(all=np.eye(2)[:,:1])
    elif bad=='overlap':parts=dict(a=np.eye(2)[:,:1],b=np.eye(2)[:,:1])
    elif bad=='nonorthogonal':parts=dict(all=2*np.eye(2))
    elif bad=='nonfinite_basis':parts['all'][0,0]=np.nan
    elif bad=='complex_x':x=x.astype(complex)
    elif bad=='complex_basis':parts['all']=parts['all'].astype(complex)
    elif bad=='bad_name':parts={' ':np.eye(2)}
    with pytest.raises(ValueError):sampled_path_work(x,g,energy,parts)


def test_inputs_are_owned_and_not_mutated():
    x=np.array([[0.,0.],[1.,2.]]);g=x.copy();energy=np.array([0.,2.5]);basis=np.eye(2)
    original=[v.copy() for v in (x,g,energy,basis)]
    sampled_path_work(x,g,energy,dict(all=basis))
    for observed,old in zip((x,g,energy,basis),original):np.testing.assert_array_equal(observed,old)
