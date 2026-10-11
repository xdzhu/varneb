"""Shared visible well response without an extra initial Hessian."""

from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

from vcneb import (
    ControlledStationaryModel, SharedInitialEnergyResponse,
    StationaryResponseContract, stationary_gap_with_shared_initial_energy,
    stationary_quadratic_branch,
)


def context(**changes):
    values = dict(physical_contract={'INPUT': 'a'*64, 'KPT': 'b'*64},
                  mechanical_family='synthetic_shared_substrate', boundary_definition_sha256='c'*64,
                  parameter='biaxial_strain', parameter_unit='1', anchor_parameter=0.,
                  energy_zero_id='synthetic_common_zero', formula_units=4,
                  pressure_eV_A3=0., potential='E')
    values.update(changes)
    return StationaryResponseContract(**values)


def model(index, *, release=True, rotate=False, scale=5., shift=0.):
    h = np.array([[2. if index == 0 else -2., .2, .4, 0.],
                  [.2, 4., .3, 0.], [.4, .3, 3., 0.], [0., 0., 0., -7.]])
    g = np.array([.15, -.12, .09, 2.])
    u = np.linalg.qr(np.random.default_rng(51+index).normal(size=(4, 4)))[0] if rotate else np.eye(4)
    basis = u.T
    branch = stationary_quadratic_branch(u.T@h@u, u.T@g, -10000.+.8*index+shift,
        basis[:, :1], basis[:, 1:2] if release else np.empty((4, 0)), basis[:, 2],
        admissible_internal_basis=basis[:, :2], expected_index=index,
        stability_floor=.01, internal_curvature_floor=.01)
    return ControlledStationaryModel(branch, context(), scale, (-.03, .03), 2., 'd'*64)


def direct_energy(m, parameter):
    b = m.branch
    internal = np.column_stack((b.conditional.retained_basis[:, :-1], b.conditional.eliminated_basis))
    control = b.conditional.retained_basis[:, -1]
    t = m.scale_A_per_parameter_unit*(parameter-m.contract.anchor_parameter)
    y = internal@(-np.linalg.solve(internal.T@b.hessian@internal,
                                  internal.T@(b.gradient+b.hessian@control*t))) + control*t
    return float(b.gradient@y + .5*y@b.hessian@y)


def feature(parameter=.015, *, release=True, rotate=False, shift=0., **changes):
    initial = model(0, release=release, rotate=rotate, scale=3., shift=shift)
    anchor = -10000.+shift+direct_energy(initial, 0.)
    change = direct_energy(initial, parameter)-direct_energy(initial, 0.)
    values = dict(contract=context(), parameter=parameter, anchor_energy_eV_cell=anchor,
                  target_energy_eV_cell=anchor+change, recorded_energy_change_eV_cell=change,
                  energy_identity_tolerance_eV_cell=1e-7, anchor_audit_sha256='e'*64,
                  target_audit_sha256='f'*64, target_endpoint_DFT_calls=7)
    values.update(changes)
    return SharedInitialEnergyResponse(**values)


def evaluate(f, saddle, *, gap=None, tolerance=1e-7):
    if gap is None:
        # Reference identity only; expected response is checked independently below.
        gap = saddle.branch.conditional.reference_energy-f.anchor_energy_eV_cell
    return stationary_gap_with_shared_initial_energy(f, saddle,
        reference_gap_eV_cell=gap, reference_gap_identity_tolerance_eV_cell=tolerance)


@pytest.mark.parametrize('release', [False, True])
@pytest.mark.parametrize('rotate', [False, True])
@pytest.mark.parametrize('parameter', [-.012, 0., .019])
def test_independent_initial_and_bottleneck_solves(release, rotate, parameter):
    initial = model(0, release=release, rotate=rotate, scale=3.)
    saddle = model(1, release=release, rotate=rotate, scale=5.)
    f = feature(parameter, release=release, rotate=rotate)
    result = evaluate(f, saddle, gap=.8-direct_energy(initial, 0.))
    ei, es = direct_energy(initial, parameter), direct_energy(saddle, parameter)
    ei0, es0 = direct_energy(initial, 0.), direct_energy(saddle, 0.)
    assert result.signed_stationary_gap_eV_fu == pytest.approx((.8+es-ei)/4, abs=1e-14)
    assert result.stationary_anchor_gap_eV_fu == pytest.approx((.8+es0-ei0)/4, abs=1e-14)
    assert result.gap_response_eV_fu == pytest.approx((es-es0-ei+ei0)/4, abs=1e-14)
    assert result.stationary_anchor_correction_eV_cell == pytest.approx(es0, abs=1e-14)
    assert result.shared_initial_feature is f


def test_common_energy_shift_preserves_small_recorded_increments():
    initial0 = model(0, scale=3.)
    gap = .8-direct_energy(initial0, 0.)
    original = evaluate(feature(), model(1), gap=gap)
    shifted = evaluate(feature(shift=1e8), model(1, shift=1e8), gap=gap)
    assert shifted.gap_response_eV_fu == original.gap_response_eV_fu
    assert shifted.signed_stationary_gap_eV_fu == original.signed_stationary_gap_eV_fu


def test_nonlinear_measured_well_has_no_invented_derivative_or_hessian():
    f = feature()
    nonlinear_change = f.recorded_energy_change_eV_cell+.004
    f = replace(f, target_energy_eV_cell=f.anchor_energy_eV_cell+nonlinear_change,
                recorded_energy_change_eV_cell=nonlinear_change)
    result = evaluate(f, model(1))
    assert result.initial_response_eV_cell == nonlinear_change
    assert not hasattr(result, 'initial_point')
    assert not hasattr(result, 'gap_derivative_eV_fu_per_parameter_unit')
    assert not hasattr(result, 'gap_curvature_eV_fu_per_parameter_unit2')
    assert not hasattr(result, 'certified_barrier_eV_fu')


def test_frozen_omissions_and_clamped_reactions_are_retained_at_both_conditions():
    result = evaluate(feature(release=False), model(1, release=False))
    for point in (result.bottleneck_anchor_point, result.bottleneck_point):
        assert point.declared_internal_gradient_norm < 1e-14
        assert point.unrepresented_internal_gradient_norm > 0.
        assert point.clamped_gradient_norm == pytest.approx(2.)


def test_negative_gap_and_shared_cost_are_not_clipped_or_multiplied():
    f = feature()
    saddle = model(1)
    saddle = replace(saddle, branch=replace(saddle.branch,
        conditional=replace(saddle.branch.conditional, reference_energy=f.anchor_energy_eV_cell-1.)))
    results = [evaluate(f, saddle, gap=-1.) for _ in range(4)]
    assert all(r.signed_stationary_gap_eV_fu < 0 for r in results)
    assert all(r.shared_initial_feature is f and r.shared_initial_feature.target_endpoint_DFT_calls == 7 for r in results)
    with pytest.raises(FrozenInstanceError): f.parameter = .02


@pytest.mark.parametrize('changes', [dict(mechanical_family='other'), dict(boundary_definition_sha256='1'*64),
    dict(parameter_unit='percent'), dict(anchor_parameter=.001), dict(energy_zero_id='other'),
    dict(formula_units=2), dict(physical_contract={'INPUT': '1'*64})])
def test_incompatible_physical_boundary_reference_and_normalization(changes):
    with pytest.raises(ValueError, match='same physical'):
        evaluate(replace(feature(), contract=context(**changes)), model(1))


@pytest.mark.parametrize('changes', [dict(parameter=True), dict(parameter=np.nan),
    dict(anchor_energy_eV_cell=np.inf), dict(target_energy_eV_cell=True),
    dict(recorded_energy_change_eV_cell=True), dict(energy_identity_tolerance_eV_cell=-1.),
    dict(energy_identity_tolerance_eV_cell=np.nan), dict(anchor_audit_sha256='bad'),
    dict(target_audit_sha256='A'*64), dict(target_endpoint_DFT_calls=0),
    dict(target_endpoint_DFT_calls=True), dict(target_endpoint_DFT_calls=1.5),
    dict(target_energy_eV_cell=-9990.), dict(recorded_energy_change_eV_cell=.9)])
def test_invalid_feature_metadata_and_energy_identity(changes):
    with pytest.raises(ValueError): feature(**changes)


def test_no_pressure_or_same_condition_response_silent_repair():
    with pytest.raises(ValueError, match='P=0 only'):
        feature(contract=context(pressure_eV_A3=.01, potential='H=E+PV'))
    with pytest.raises(ValueError, match='P=0 only'):
        feature(contract=context(potential='H=E+PV'))
    with pytest.raises(ValueError, match='same-condition'):
        feature(parameter=0., target_energy_eV_cell=-10000.001, recorded_energy_change_eV_cell=.002,
                energy_identity_tolerance_eV_cell=1.)


def test_bottleneck_index_reference_and_domain_guards():
    f = feature()
    with pytest.raises(ValueError, match='index1'): evaluate(f, model(0))
    with pytest.raises(ValueError, match='contradicts'): evaluate(f, model(1), gap=.9)
    with pytest.raises(ValueError, match='nonnegative'): evaluate(f, model(1), tolerance=-1.)
    with pytest.raises(ValueError): evaluate(f, model(1), gap=True)
    with pytest.raises(ValueError, match='parameter domain'):
        evaluate(feature(.04), model(1))
    with pytest.raises(ValueError, match='declared radius'):
        evaluate(f, replace(model(1), maximum_full_displacement_A=.001))
    with pytest.raises(TypeError): evaluate(f, None, gap=.8)
