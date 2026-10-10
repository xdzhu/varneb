"""Measured-column controls: no missing-block/stability or material claims."""
from dataclasses import replace

import numpy as np
import pytest

from vcneb import (MeasuredStationaryCentre,nested_stationary_controls,
                   nested_stationary_gap_responses,StationaryResponseContract)


def centre(index, *, partial=False, rotate=False, unstable_cell=False, **changes):
    h=np.array([[2. if index==0 or unstable_cell else -2., .2,.3,.4,0],
                [.2,4.,.1,.8,0],[.3,.1,-.6 if unstable_cell else 3.,1.1,0],
                [.4,.8,1.1,6.,0],[0,0,0,0,-7.]])
    g=np.array([.1,-.2,.05,.06,2.])
    u=np.linalg.qr(np.random.default_rng(192+index).normal(size=(5,5)))[0] if rotate else np.eye(5)
    e=u.T
    m=e[:,[0,3]] if partial else e[:,:4]
    values=dict(measured_basis=m,hessian_columns=u.T@h@u@m,gradient=u.T@g,
        reference_potential_eV_cell=-10000.+.8*index,retained_internal_basis=e[:,:1],
        admissible_internal_basis=e[:,:3],atomic_internal_basis=e[:,:2],control_direction=e[:,3],
        expected_index=index,stability_floor_eV_A2=.01,internal_curvature_floor_eV_A2=.01,
        symmetry_tolerance_eV_A2=1e-8,source_audit_sha256=('a' if index==0 else 'b')*64,
        measurement_ids=('shared-reference-basis-evaluation',f'centre-{index}-plus',f'centre-{index}-minus'))
    values.update(changes)
    return MeasuredStationaryCentre(**values)


def result_map(c):
    return {r.control:r for r in nested_stationary_controls(c)}


def contract():
    return StationaryResponseContract({'INPUT':'a'*64},'synthetic_clamped','c'*64,
        'biaxial_strain','1',0.,'synthetic_zero',4,0.,'E')


def domain(scale):
    return dict(scale_A_per_parameter_unit=scale,parameter_interval=(-.03,.03),maximum_full_displacement_A=2.)


def paired(a,b,x,**changes):
    args=dict(contract=contract(),initial_domain=domain(3.),bottleneck_domain=domain(5.),
              reference_gap_eV_cell=.8,reference_gap_identity_tolerance_eV_cell=1e-8)
    args.update(changes)
    return nested_stationary_gap_responses(a,b,x,**args)


def direct(c,internal,t):
    h=c.action_extension()
    y=internal@(-np.linalg.solve(internal.T@h@internal,
               internal.T@(c.gradient+h@c.control_direction*t)))+c.control_direction*t
    return float(c.gradient@y+.5*y@h@y),y


@pytest.mark.parametrize('index',[0,1])
@pytest.mark.parametrize('rotate',[False,True])
def test_nested_controls_match_independent_stationary_solves(index,rotate):
    c=centre(index,rotate=rotate)
    results=result_map(c)
    expected_dims=((1,0),(1,1),(1,2),(1,2))
    for r,(q,rn) in zip(results.values(),expected_dims):
        assert r.branch is not None and r.reason is None and r.promoted_dimension==0
        assert (r.retained_dimension,r.released_dimension)==(q,rn)
        assert r.unmeasured_internal_dimension==0
        internal=np.column_stack((r.branch.conditional.retained_basis[:,:-1],r.branch.conditional.eliminated_basis))
        for t in (-.04,0.,.06):
            e,y=direct(c,internal,t)
            point=r.branch.evaluate(t)
            assert point.energy_change==pytest.approx(e,abs=1e-14)
            np.testing.assert_allclose(point.full_displacement,y,atol=1e-14)
            assert point.declared_internal_gradient_norm<1e-13
            assert point.clamped_gradient_norm==pytest.approx(2.)
    assert results['B2_frozen'].branch.evaluate(.06).unrepresented_internal_gradient_norm>0
    assert results['B3_atomic_release'].branch.evaluate(.06).unrepresented_internal_gradient_norm>0
    assert results['B4_joint_release'].branch.evaluate(.06).unrepresented_internal_gradient_norm<1e-13
    assert results['B4_joint_release'].branch.control_curvature==pytest.approx(
        results['B5_training_instability_promotion'].branch.control_curvature,abs=1e-14)


def test_partial_columns_support_frozen_control_but_never_invent_stable_complement():
    c=centre(1,partial=True)
    results=result_map(c)
    frozen=results['B2_frozen']
    assert frozen.branch is not None and frozen.unmeasured_internal_dimension==2
    # Off-measured rows of the physical-gradient response are measured, not dropped.
    p=frozen.branch.evaluate(.02)
    physical=c.gradient+c.hessian_columns@np.linalg.lstsq(c.measured_basis,p.full_displacement,rcond=None)[0]
    assert p.unrepresented_internal_gradient_norm==pytest.approx(np.linalg.norm(physical[1:3]))
    for name in ('B3_atomic_release','B4_joint_release','B5_training_instability_promotion'):
        assert results[name].branch is None and results[name].status=='unavailable_unmeasured_curvature'


def test_mixed_measured_direction_is_not_full_control_or_full_internal_coverage():
    base=centre(0)
    m=np.array([[1],[0],[0],[1],[0]])/np.sqrt(2)
    c=replace(base,measured_basis=m,hessian_columns=base.action_extension()@m)
    results=result_map(c)
    assert all(r.branch is None and r.unmeasured_internal_dimension==3 for r in results.values())


def test_B5_promotes_training_unstable_cell_without_pseudoinverse_or_extra_data():
    c=centre(1,unstable_cell=True)
    results=result_map(c)
    assert all(results[n].branch is None for n in ('B2_frozen','B3_atomic_release','B4_joint_release'))
    b5=results['B5_training_instability_promotion']
    assert b5.branch is not None and b5.promoted_dimension==1
    assert b5.retained_dimension==2 and b5.released_dimension==1
    assert b5.measurement_ids==c.measurement_ids
    e,y=direct(c,c.admissible_internal_basis,.025)
    p=b5.branch.evaluate(.025)
    assert p.energy_change==pytest.approx(e,abs=1e-14)
    np.testing.assert_allclose(p.full_displacement,y,atol=1e-14)


def test_zero_or_additional_negative_full_direction_is_not_hidden():
    c=centre(1)
    h=c.action_extension(); h[2,:]=h[:,2]=0
    c=replace(c,hessian_columns=h@c.measured_basis)
    assert result_map(c)['B5_training_instability_promotion'].branch is None
    h[2,2]=-3.
    c=replace(c,hessian_columns=h@c.measured_basis)
    assert result_map(c)['B5_training_instability_promotion'].branch is None


def test_mixed_retained_atomic_cell_basis_cannot_masquerade_as_atomic_only():
    c=centre(0)
    q=np.array([[1],[0],[1],[0],[0]])/np.sqrt(2)
    c=replace(c,retained_internal_basis=q)
    result=result_map(c)['B3_atomic_release']
    assert result.branch is None and result.status=='unavailable_incompatible_atomic_ablation'


@pytest.mark.parametrize('rotate',[False,True])
def test_paired_gap_keeps_IS_response_offsets_and_distinct_control_scales(rotate):
    a,b=centre(0,rotate=rotate),centre(1,rotate=rotate)
    results=paired(a,b,.005)
    for r in results:
        assert r.response is not None
        p=r.response
        assert abs(p.initial_response_eV_cell)>1e-4
        assert abs(p.stationary_anchor_correction_eV_cell)>1e-4
        assert len(r.available_measurement_ids)==5  # shared reference preparation not counted twice
        dx=1e-5
        lo={r.control:r for r in paired(a,b,.005-dx)}[r.control].response
        hi={r.control:r for r in paired(a,b,.005+dx)}[r.control].response
        assert (hi.signed_stationary_gap_eV_fu-lo.signed_stationary_gap_eV_fu)/(2*dx)==pytest.approx(
            p.gap_derivative_eV_fu_per_parameter_unit,abs=1e-10)
        assert (hi.gap_derivative_eV_fu_per_parameter_unit-lo.gap_derivative_eV_fu_per_parameter_unit)/(2*dx)==pytest.approx(
            p.gap_curvature_eV_fu_per_parameter_unit2,abs=1e-9)
    b4,b5=results[-2:]
    assert b4.response.gap_response_eV_fu==pytest.approx(b5.response.gap_response_eV_fu,abs=1e-14)


def test_paired_controls_report_coverage_and_domain_failure_not_zero_prediction():
    partial=paired(centre(0,partial=True),centre(1,partial=True),.005)
    assert partial[0].response is not None and all(r.response is None for r in partial[1:])
    outside=paired(centre(0),centre(1),.04)
    assert all(r.response is None and r.status=='unavailable_local_domain' for r in outside)


def test_wrong_raw_gap_and_unregistered_domain_fields_fail_closed():
    for change in ({'reference_gap_eV_cell':.9},{'initial_domain':{'scale_A_per_parameter_unit':3.}},
                   {'reference_gap_identity_tolerance_eV_cell':-.1}):
        with pytest.raises(ValueError):
            paired(centre(0),centre(1),.005,**change)


@pytest.mark.parametrize('domain_change',[{'scale_A_per_parameter_unit':0},
    {'maximum_full_displacement_A':np.nan},{'parameter_interval':(.01,.02)},
    {'parameter_interval':(True,.02)}])
def test_invalid_domain_is_rejected_even_when_all_controls_are_unavailable(domain_change):
    d=domain(3.);d.update(domain_change)
    with pytest.raises(ValueError):
        paired(centre(0,partial=True),centre(1,partial=True),.005,initial_domain=d)


@pytest.mark.parametrize('fault',['complex','reciprocity','forbidden_measurement','forbidden_retained',
                                 'duplicate_ids','index','scale_basis','nan','control'])
def test_invalid_measured_data_rejected_without_repair(fault):
    c=centre(0)
    changes={}
    if fault=='complex': changes['hessian_columns']=c.hessian_columns.astype(complex)+1j
    elif fault=='reciprocity':
        d=c.hessian_columns.copy(); d[0,1]+=1e-3; changes['hessian_columns']=d
    elif fault=='forbidden_measurement':
        m=c.measured_basis.copy(); m[:,0]=np.eye(5)[:,4]; changes['measured_basis']=m
    elif fault=='forbidden_retained': changes['retained_internal_basis']=np.eye(5)[:,4:5]
    elif fault=='duplicate_ids': changes['measurement_ids']=('same','same')
    elif fault=='index': changes['expected_index']=True
    elif fault=='scale_basis': changes['measured_basis']=c.measured_basis*2
    elif fault=='nan': changes['gradient']=np.array([np.nan,0,0,0,0])
    elif fault=='control': changes['control_direction']=np.eye(5)[:,0]
    with pytest.raises(ValueError): replace(c,**changes)


def test_data_owned_not_mutated_by_control_selection():
    c=centre(1,unstable_cell=True)
    before=c.hessian_columns.copy()
    nested_stationary_controls(c)
    np.testing.assert_array_equal(c.hessian_columns,before)
    assert not c.hessian_columns.flags.writeable


def test_rotating_measured_column_basis_preserves_actions_and_controls():
    c=centre(1)
    rotation=np.linalg.qr(np.random.default_rng(98).normal(size=(4,4)))[0]
    rotated=replace(c,measured_basis=c.measured_basis@rotation,
                    hessian_columns=c.hessian_columns@rotation)
    for a,b in zip(nested_stationary_controls(c),nested_stationary_controls(rotated)):
        assert a.status==b.status and a.unmeasured_internal_dimension==b.unmeasured_internal_dimension
        assert a.branch.evaluate(.02).energy_change==pytest.approx(b.branch.evaluate(.02).energy_change,abs=1e-14)
        np.testing.assert_allclose(a.branch.evaluate(.02).full_displacement,
                                   b.branch.evaluate(.02).full_displacement,atol=1e-14)


def test_negative_external_control_curvature_is_not_an_extra_internal_instability():
    c=centre(1)
    h=c.action_extension(); h[3,3]=-8.
    c=replace(c,hessian_columns=h@c.measured_basis)
    assert all(r.branch is not None for r in nested_stationary_controls(c))
