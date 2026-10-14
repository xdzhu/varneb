"""Independent tangent/work checks; not source-mode or material validation."""
import numpy as np
import pytest
from ase import Atoms

from vcneb import (ReferenceCellCoordinates, ClampedModeCoordinates,
    BiaxialClampedModeCoordinates, clamped_plane_vcneb_boundary,
    transport_atomic_reference_span)


def pair(kind, cell_basis):
    source=Atoms('CuAgCuAg',cell=[[4.1,.2,.1],[.1,3.9,.2],[.2,.1,4.3]],
        scaled_positions=[[0,0,0],[.5,0,.5],[0,.5,.5],[.5,.5,0]],pbc=True)
    rotation=np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]])
    target=source.copy(); target.set_cell(cell_basis@source.cell.array@rotation.T,scale_atoms=False)
    target.set_scaled_positions(source.get_scaled_positions(wrap=False)@np.linalg.inv(cell_basis))
    target.cell[2]+=[.03,-.02,.01]
    if kind=='free': chart=ReferenceCellCoordinates(target)
    else:
        boundary=clamped_plane_vcneb_boundary(4,target.cell.array,allow_tilt=True)
        chart=(ClampedModeCoordinates(target,boundary=boundary,cell_scale_A=4.) if kind=='clamped'
            else BiaxialClampedModeCoordinates(target,boundary=boundary,cell_scale_A=4.,
                                              anchor_strain=0.,strain_scale_A=4.))
    return source,chart


@pytest.mark.parametrize('kind',['free','clamped','controlled'])
@pytest.mark.parametrize('p',[np.eye(3),np.array([[1,1,0],[0,1,0],[0,0,1]]),
                              np.array([[0,0,1],[1,0,0],[0,1,0]])])
def test_ordered_fractional_tangent_mass_scale_and_full_virtual_work(kind,p):
    source,chart=pair(kind,p);rng=np.random.default_rng(46)
    values=np.linalg.qr(rng.normal(size=(12,4)))[0]; masses=np.array([63.5,107.9,63.5,107.9])
    before=values.copy();positions=source.positions.copy();cell=source.cell.array.copy()
    out=transport_atomic_reference_span(source,values,masses,chart,source_to_target_cell_basis=p)
    expected=np.empty((12,4))
    for j in range(4):
        fractional=(values[:,j].reshape(4,3)/np.sqrt(masses)[:,None])@np.linalg.inv(source.cell.array)
        target_fractional=np.linalg.solve(p.T,fractional.T).T
        expected[:,j]=(target_fractional@chart.reference_atoms.cell.array).ravel()
    np.testing.assert_allclose(out.raw_mapped_columns[:12],expected,atol=1e-14)
    np.testing.assert_allclose(out.local_basis.T@out.local_basis,np.eye(4),atol=1e-13)
    np.testing.assert_allclose(out.local_basis@out.source_to_local+out.excluded_translation_columns,
                               out.raw_mapped_columns,atol=1e-14)
    assert np.array_equal(out.local_basis[12:],np.zeros_like(out.local_basis[12:]))
    assert out.rank==4
    g=rng.normal(size=chart.coordinate_count);a=rng.normal(size=4)
    direct=float(g[:12]@expected@a)
    internal=float(g@out.local_basis@out.source_to_local@a)
    translation=float(g@out.excluded_translation_columns@a)
    assert direct==pytest.approx(internal+translation,abs=1e-13)
    assert abs(translation)>1e-5  # do not erase the deliberately raw net-force work
    atom_shifts=(out.raw_mapped_columns@a)[:12].reshape(4,3)
    gauge=(out.excluded_translation_columns@a)[:12].reshape(4,3)
    np.testing.assert_allclose(gauge,np.tile(atom_shifts.mean(axis=0),(4,1)),atol=1e-14)
    np.testing.assert_array_equal(values,before);np.testing.assert_array_equal(source.positions,positions)
    np.testing.assert_array_equal(source.cell.array,cell)
    for value in (out.local_basis,out.raw_mapped_columns,out.excluded_translation_columns,out.source_to_local):
        assert not value.flags.writeable


def test_basis_mixing_preserves_span_not_named_modes():
    source,chart=pair('controlled',np.eye(3));rng=np.random.default_rng(72)
    v=rng.normal(size=(12,3));r=np.array([[2.,.1,0],[0,1,.2],[.1,0,3.]])
    a=transport_atomic_reference_span(source,v,source.get_masses(),chart,source_to_target_cell_basis=np.eye(3))
    b=transport_atomic_reference_span(source,v@r,source.get_masses(),chart,source_to_target_cell_basis=np.eye(3))
    np.testing.assert_allclose(a.local_basis@a.local_basis.T,b.local_basis@b.local_basis.T,atol=1e-13)


@pytest.mark.parametrize('fault',['complex','mass','order','integer','orientation','dependent','translation','nan'])
def test_invalid_source_and_rank_loss_are_not_silently_repaired(fault):
    source,chart=pair('controlled',np.eye(3));masses=source.get_masses();p=np.eye(3)
    v=np.random.default_rng(19).normal(size=(12,3))
    if fault=='complex': v=v.astype(complex)+1j*1e-12
    if fault=='mass': masses[0]=0
    if fault=='order': source.numbers[0]=8
    if fault=='integer': p[0,1]=.5
    if fault=='orientation': p[0,0]=-1
    if fault=='dependent': v[:,1]=v[:,0]
    if fault=='translation': v=(np.sqrt(masses)[:,None]*np.tile([1.,0,0],(4,1))).reshape(12,1)
    if fault=='nan': v[0,0]=np.nan
    with pytest.raises(ValueError):
        transport_atomic_reference_span(source,v,masses,chart,source_to_target_cell_basis=p)
