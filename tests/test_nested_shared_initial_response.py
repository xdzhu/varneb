"""One visible well feature across B2--B5, without an initial Hessian.

The oracle solves the original physical blocks independently of the returned
branches. Synthetic checks are not HfO2 errors or prospective material gains.
"""
from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

from vcneb import (
    MeasuredStationaryCentre, SharedInitialEnergyResponse,
    StationaryResponseContract, nested_stationary_gaps_with_shared_initial_energy,
)


def feature(parameter=.015, *, shift=0., increment=.003):
    contract=StationaryResponseContract({'INPUT':'a'*64},'synthetic_clamped','b'*64,
        'biaxial_strain','1',0.,'synthetic_common_zero',4,0.,'E')
    anchor=-10000.04+shift
    if parameter==0.:
        increment=0.
    return SharedInitialEnergyResponse(contract,parameter,anchor,anchor+increment,
        increment,1e-7,'c'*64,'d'*64,7)


def centre(*, rotate=False, unstable_cell=False, partial=False, shift=0.):
    h=np.array([[2. if unstable_cell else -2.,.2,.3,.4,0.],
                [.2,4.,.1,.8,0.],[.3,.1,-.6 if unstable_cell else 3.,1.1,0.],
                [.4,.8,1.1,6.,0.],[0.,0.,0.,0.,-7.]])
    g=np.array([.1,-.2,.05,.06,2.])
    rotation=np.linalg.qr(np.random.default_rng(64).normal(size=(5,5)))[0] if rotate else np.eye(5)
    e=rotation.T
    h,g=rotation.T@h@rotation,rotation.T@g
    m=e[:,[0,3]] if partial else e[:,:4]
    c=MeasuredStationaryCentre(m,h@m,g,-9999.2+shift,e[:,:1],e[:,:3],e[:,:2],
        e[:,3],1,.01,.01,1e-8,'e'*64,('bottleneck-centre','bottleneck-plus','bottleneck-minus'))
    return c,h


def domain(**changes):
    d=dict(scale_A_per_parameter_unit=5.,parameter_interval=(-.03,.03),maximum_full_displacement_A=2.)
    d.update(changes)
    return d


def evaluate(f,c,**changes):
    args=dict(bottleneck_domain=domain(),reference_gap_eV_cell=.84,
        reference_gap_identity_tolerance_eV_cell=1e-7)
    args.update(changes)
    return nested_stationary_gaps_with_shared_initial_energy(f,c,**args)


def direct(h,c,internal,parameter):
    t=5.*parameter
    y=internal@(-np.linalg.solve(internal.T@h@internal,
        internal.T@(c.gradient+h@c.control_direction*t)))+c.control_direction*t
    return float(c.gradient@y+.5*y@h@y),y


@pytest.mark.parametrize('rotate',[False,True])
@pytest.mark.parametrize('parameter',[-.012,0.,.019])
def test_four_arms_match_original_block_oracle_and_share_one_well(rotate,parameter):
    c,h=centre(rotate=rotate)
    f=feature(parameter)
    results=evaluate(f,c)
    internals=(c.retained_internal_basis,c.atomic_internal_basis,
               c.admissible_internal_basis,c.admissible_internal_basis)
    assert [r.control for r in results]==[
        'B2_frozen','B3_atomic_release','B4_joint_release','B5_training_instability_promotion']
    for arm,internal in zip(results,internals):
        response=arm.response
        assert response is not None and arm.reason is None
        assert arm.status=='available_shared_initial_restricted_quadratic_gap_only'
        es,ys=direct(h,c,internal,parameter)
        es0,ys0=direct(h,c,internal,0.)
        assert response.signed_stationary_gap_eV_fu==pytest.approx(
            (.84+es-f.recorded_energy_change_eV_cell)/4,abs=1e-14)
        assert response.gap_response_eV_fu==pytest.approx(
            (es-es0-f.recorded_energy_change_eV_cell)/4,abs=1e-14)
        assert response.stationary_anchor_gap_eV_fu==pytest.approx((.84+es0)/4,abs=1e-14)
        assert response.stationary_anchor_correction_eV_cell==pytest.approx(es0,abs=1e-14)
        np.testing.assert_allclose(response.bottleneck_anchor_point.full_displacement,ys0,atol=1e-14)
        np.testing.assert_allclose(response.bottleneck_point.full_displacement,ys,atol=1e-14)
        assert arm.shared_initial_feature is f and response.shared_initial_feature is f
        assert arm.available_measurement_ids==c.measurement_ids
        assert arm.bottleneck.source_audit_sha256==c.source_audit_sha256
        assert response.shared_initial_feature.target_endpoint_DFT_calls==7
        assert not hasattr(response,'initial_point')
        assert not hasattr(response,'gap_derivative_eV_fu_per_parameter_unit')
        assert not hasattr(response,'gap_curvature_eV_fu_per_parameter_unit2')
    assert results[2].response.signed_stationary_gap_eV_fu==pytest.approx(
        results[3].response.signed_stationary_gap_eV_fu,abs=1e-14)


def test_omitted_movables_and_clamped_reactions_survive_at_both_conditions():
    results=evaluate(feature(),centre()[0])
    for arm in results:
        for point in (arm.response.bottleneck_anchor_point,arm.response.bottleneck_point):
            assert point.declared_internal_gradient_norm<1e-13
            assert point.clamped_gradient_norm==pytest.approx(2.)
            if arm.control in ('B2_frozen','B3_atomic_release'):
                assert point.unrepresented_internal_gradient_norm>0.
            else:
                assert point.unrepresented_internal_gradient_norm<1e-13


def test_partial_coverage_keeps_all_arms_and_same_feature_without_fabricating_curvature():
    c,_=centre(partial=True)
    f=feature()
    results=evaluate(f,c)
    assert results[0].response is not None
    for arm in results[1:]:
        assert arm.response is None and arm.reason==arm.bottleneck.reason
        assert 'outside the measured space' in arm.reason
        assert arm.bottleneck.status=='unavailable_unmeasured_curvature'
        assert arm.status=='unavailable_control_assumption_or_measurement'
        assert arm.shared_initial_feature is f and arm.available_measurement_ids==c.measurement_ids


@pytest.mark.parametrize('rotate',[False,True])
def test_B5_promotes_existing_training_columns_and_matches_full_internal_oracle(rotate):
    c,h=centre(unstable_cell=True,rotate=rotate)
    f=feature()
    results=evaluate(f,c)
    assert all(r.response is None for r in results[:3])
    promoted=results[3]
    assert promoted.bottleneck.promoted_dimension==1
    es,ys=direct(h,c,c.admissible_internal_basis,f.parameter)
    assert promoted.response.signed_stationary_gap_eV_fu==pytest.approx(
        (.84+es-f.recorded_energy_change_eV_cell)/4,abs=1e-14)
    np.testing.assert_allclose(promoted.response.bottleneck_point.full_displacement,ys,atol=1e-14)
    assert all(r.shared_initial_feature is f for r in results)
    assert promoted.available_measurement_ids==c.measurement_ids


def test_mixed_atomic_cell_retained_space_does_not_redefine_atomic_control():
    c,_=centre()
    q=np.array([[1.],[0.],[1.],[0.],[0.]])/np.sqrt(2.)
    c=replace(c,retained_internal_basis=q)
    result=evaluate(feature(),c)[1]
    assert result.response is None
    assert result.bottleneck.status=='unavailable_incompatible_atomic_ablation'


def test_negative_gaps_and_nonlinear_initial_increment_are_not_clipped_or_recalibrated():
    c,_=centre()
    f=feature(increment=1.1)
    results=evaluate(f,c)
    assert all(r.response.signed_stationary_gap_eV_fu<0. for r in results)
    assert all(r.response.initial_response_eV_cell==1.1 for r in results)
    assert all(not hasattr(r.response,'certified_barrier_eV_fu') for r in results)


def test_common_raw_energy_shift_does_not_change_recorded_small_response():
    original=evaluate(feature(),centre()[0])
    shifted=evaluate(feature(shift=1e8),centre(shift=1e8)[0])
    for before,after in zip(original,shifted):
        assert after.response.gap_response_eV_fu==before.response.gap_response_eV_fu
        assert after.response.signed_stationary_gap_eV_fu==before.response.signed_stationary_gap_eV_fu


@pytest.mark.parametrize('changes',[
    dict(bottleneck_domain=domain(scale_A_per_parameter_unit=0.)),
    dict(bottleneck_domain=domain(maximum_full_displacement_A=np.nan)),
    dict(bottleneck_domain=domain(parameter_interval=(.01,.02))),
    dict(bottleneck_domain=domain(parameter_interval=(True,.02))),
    dict(bottleneck_domain=domain(extra_field=1)),
    dict(bottleneck_domain={}),dict(reference_gap_eV_cell=.9),
    dict(reference_gap_eV_cell=True),dict(reference_gap_eV_cell=np.inf),
    dict(reference_gap_identity_tolerance_eV_cell=-1.),
    dict(reference_gap_identity_tolerance_eV_cell=True),
])
def test_invalid_global_declarations_fail_even_when_no_arm_can_evaluate(changes):
    c,h=centre()
    m=(np.eye(5)[:,0:1]+np.eye(5)[:,3:4])/np.sqrt(2.)
    c=replace(c,measured_basis=m,hessian_columns=h@m)
    assert all(r.response is None for r in evaluate(feature(),c))
    with pytest.raises(ValueError):
        evaluate(feature(),c,**changes)


@pytest.mark.parametrize('parameter,radius',[ (.04,2.),(.015,.001)])
def test_local_domain_failure_is_abstention_not_zero_or_automatic_recentering(parameter,radius):
    results=evaluate(feature(parameter),centre()[0],
        bottleneck_domain=domain(maximum_full_displacement_A=radius))
    assert all(r.response is None and r.status=='unavailable_local_domain' for r in results)
    assert all(r.reason and r.shared_initial_feature.parameter==parameter for r in results)


def test_wrong_types_index_and_owned_data_are_guarded():
    c,h=centre()
    f=feature()
    with pytest.raises(TypeError): evaluate(None,c)
    with pytest.raises(TypeError): evaluate(f,None)
    with pytest.raises(ValueError,match='index1'): evaluate(f,replace(c,expected_index=0))
    before=c.hessian_columns.copy()
    results=evaluate(f,c)
    np.testing.assert_array_equal(c.hessian_columns,before)
    assert not c.hessian_columns.flags.writeable
    with pytest.raises(FrozenInstanceError): results[0].status='changed'
    with pytest.raises(FrozenInstanceError): f.parameter=.02
