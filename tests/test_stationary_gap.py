"""Paired local response: shared contracts, offsets, scales and explicit limits."""

from dataclasses import replace

import numpy as np
import pytest

from vcneb import (
    ControlledStationaryModel, StationaryResponseContract,
    restricted_stationary_gap_response, stationary_quadratic_branch,
)


def contract(**changes):
    values = dict(physical_contract={'INPUT': 'a'*64, 'KPT': 'b'*64},
                  mechanical_family='test_substrate_tilt_open', boundary_definition_sha256='c'*64,
                  parameter='biaxial_strain', parameter_unit='1', anchor_parameter=0.,
                  energy_zero_id='synthetic_common_zero', formula_units=4,
                  pressure_eV_A3=0., potential='E')
    values.update(changes)
    return StationaryResponseContract(**values)


def model(index, energy, scale, *, release=True, rotate=False, coordinate_rescale=1., **changes):
    h = np.array([[2. if index == 0 else -2., .2, .4, 0.],
                  [.2, 4., .3, 0.], [.4, .3, 3., 0.], [0., 0., 0., -7.]])
    g = np.array([.15, -.12, .09, 2.])
    transform = np.eye(4); transform[2, 2] /= coordinate_rescale
    h, g = transform.T @ h @ transform, transform.T @ g
    u = np.linalg.qr(np.random.default_rng(51+index).normal(size=(4, 4)))[0] if rotate else np.eye(4)
    basis = u.T
    branch = stationary_quadratic_branch(u.T @ h @ u, u.T @ g, energy,
        basis[:, :1], basis[:, 1:2] if release else np.empty((4, 0)), basis[:, 2],
        admissible_internal_basis=basis[:, :2], expected_index=index,
        stability_floor=.01, internal_curvature_floor=.01)
    values = dict(branch=branch, contract=contract(), scale_A_per_parameter_unit=scale*coordinate_rescale,
                  parameter_interval=(-.03, .03), maximum_full_displacement_A=2., source_audit_sha256='d'*64)
    values.update(changes)
    return ControlledStationaryModel(**values)


def pair(**changes):
    return model(0, -10000., 3., **changes), model(1, -9999.2, 5., **changes)


def evaluate(initial, saddle, x, gap=.8, tolerance=1e-8):
    return restricted_stationary_gap_response(initial, saddle, x, reference_gap_eV_cell=gap,
                                              reference_gap_identity_tolerance_eV_cell=tolerance)


def direct(model, x):
    branch = model.branch
    internal = np.column_stack((branch.conditional.retained_basis[:, :-1], branch.conditional.eliminated_basis))
    c = branch.conditional.retained_basis[:, -1]
    t = model.scale_A_per_parameter_unit*(x-model.contract.anchor_parameter)
    solution = -np.linalg.solve(internal.T @ branch.hessian @ internal,
                               internal.T @ (branch.gradient + branch.hessian @ c*t))
    y = internal @ solution + c*t
    energy = float(branch.gradient @ y + .5*y @ branch.hessian @ y)
    gradient = float(c @ (branch.gradient + branch.hessian @ y))
    return energy, model.scale_A_per_parameter_unit*gradient, y


@pytest.mark.parametrize('release', [False, True])
@pytest.mark.parametrize('rotate', [False, True])
def test_separate_stationary_solves_and_finite_gap_derivatives(release, rotate):
    initial, saddle = pair(release=release, rotate=rotate)
    for x in (-.012, 0., .019):
        result = evaluate(initial, saddle, x)
        ei, gi, yi = direct(initial, x); es, gs, ys = direct(saddle, x)
        ei0, _, _ = direct(initial, 0.); es0, _, _ = direct(saddle, 0.)
        assert result.stationary_anchor_correction_eV_cell == pytest.approx(es0-ei0, abs=1e-14)
        assert result.stationary_anchor_gap_eV_fu == pytest.approx((.8+es0-ei0)/4, abs=1e-14)
        assert result.initial_response_eV_cell == pytest.approx(ei-ei0, abs=1e-14)
        assert result.bottleneck_response_eV_cell == pytest.approx(es-es0, abs=1e-14)
        assert result.signed_stationary_gap_eV_fu == pytest.approx((.8+es-ei)/4, abs=1e-14)
        assert result.gap_response_eV_fu == pytest.approx((es-es0-ei+ei0)/4, abs=1e-14)
        assert result.gap_derivative_eV_fu_per_parameter_unit == pytest.approx((gs-gi)/4, abs=1e-14)
        np.testing.assert_allclose(result.initial_point.full_displacement, yi, atol=1e-14)
        np.testing.assert_allclose(result.bottleneck_point.full_displacement, ys, atol=1e-14)
        dx = 1e-5
        a, b = evaluate(initial, saddle, x-dx), evaluate(initial, saddle, x+dx)
        assert (b.signed_stationary_gap_eV_fu-a.signed_stationary_gap_eV_fu)/(2*dx) == pytest.approx(
            result.gap_derivative_eV_fu_per_parameter_unit, abs=1e-10)
        assert (b.gap_derivative_eV_fu_per_parameter_unit-a.gap_derivative_eV_fu_per_parameter_unit)/(2*dx) == pytest.approx(
            result.gap_curvature_eV_fu_per_parameter_unit2, abs=1e-9)


def test_initial_response_cannot_be_omitted_or_reference_offset_hidden():
    initial, saddle = pair(); result = evaluate(initial, saddle, .02)
    assert abs(result.initial_response_eV_cell) > 1e-4
    assert abs(result.stationary_anchor_correction_eV_cell) > 1e-4
    wrong_saddle_only = result.bottleneck_response_eV_cell/4
    assert abs(wrong_saddle_only-result.gap_response_eV_fu) > 1e-5
    assert result.stationary_anchor_gap_eV_fu != pytest.approx(.8/4, abs=1e-4)


def test_common_energy_shift_and_independent_control_scales():
    original = pair(); shifted = pair()
    shifted = tuple(replace(m, branch=replace(m.branch, conditional=replace(m.branch.conditional,
                    reference_energy=m.branch.conditional.reference_energy+1e8))) for m in shifted)
    rescaled = (model(0, -10000., 3., coordinate_rescale=3.), model(1, -9999.2, 5., coordinate_rescale=4.))
    for x in (-.01, .015):
        expected = evaluate(*original, x)
        for models, tolerance in ((shifted, 1e-7), (rescaled, 1e-8)):
            actual = evaluate(*models, x, tolerance=tolerance)
            assert actual.signed_stationary_gap_eV_fu == pytest.approx(expected.signed_stationary_gap_eV_fu, abs=1e-14)
            assert actual.gap_response_eV_fu == pytest.approx(expected.gap_response_eV_fu, abs=1e-14)
            assert actual.gap_derivative_eV_fu_per_parameter_unit == pytest.approx(expected.gap_derivative_eV_fu_per_parameter_unit, abs=1e-14)
            assert actual.gap_curvature_eV_fu_per_parameter_unit2 == pytest.approx(expected.gap_curvature_eV_fu_per_parameter_unit2, abs=1e-13)


def test_partial_space_residuals_and_clamped_reactions_are_preserved():
    models = pair(release=False); result = evaluate(*models, .02)
    for model, point in zip(models, (result.initial_point, result.bottleneck_point)):
        _, _, y = direct(model, .02)
        # Explicit independent omitted-coordinate gradient, not a guessed bound.
        omitted = model.branch.gradient[1] + model.branch.hessian[1] @ y
        assert point.declared_internal_gradient_norm < 1e-14
        assert point.unrepresented_internal_gradient_norm == pytest.approx(abs(omitted), abs=1e-14)
        assert abs(omitted) > 1e-14
        assert point.clamped_gradient_norm == pytest.approx(2.)
    assert not hasattr(result, 'certified_barrier_eV_fu')


def test_signed_negative_gap_is_not_repaired_into_activation_barrier():
    result = evaluate(model(0, -10000., 3.), model(1, -10001., 5.), .01, gap=-1.)
    assert result.signed_stationary_gap_eV_fu < 0


def test_contract_and_domain_ownership():
    hashes = {'INPUT': 'a'*64}; context = contract(physical_contract=hashes)
    hashes['INPUT'] = 'b'*64
    assert context.physical_contract['INPUT'] == 'a'*64
    with pytest.raises(TypeError): context.physical_contract['INPUT'] = 'c'*64
    bounds = [-.03, .03]; owned = model(0, -1., 3., parameter_interval=bounds)
    bounds[0] = .04
    assert owned.parameter_interval == (-.03, .03)


@pytest.mark.parametrize('changes', [dict(mechanical_family='other'), dict(boundary_definition_sha256='e'*64),
    dict(parameter='other'), dict(parameter_unit='percent'), dict(anchor_parameter=.001),
    dict(energy_zero_id='other'), dict(formula_units=2), dict(potential='H=E+PV'),
    dict(physical_contract={'INPUT':'e'*64}), dict(potential='H=E+PV', pressure_eV_A3=.01)])
def test_incompatible_contracts_are_rejected(changes):
    initial, saddle = pair()
    with pytest.raises(ValueError, match='same physical'):
        evaluate(initial, replace(saddle, contract=contract(**changes)), .01)


@pytest.mark.parametrize('changes', [dict(physical_contract={}), dict(boundary_definition_sha256='bad'),
    dict(parameter=''), dict(formula_units=True), dict(formula_units=0), dict(anchor_parameter=np.nan),
    dict(pressure_eV_A3=.01), dict(pressure_eV_A3=True), dict(potential='F')])
def test_invalid_contracts(changes):
    with pytest.raises(ValueError): contract(**changes)


def test_finite_pressure_requires_declared_same_enthalpy_potential():
    context = contract(pressure_eV_A3=.01, potential='H=E+PV')
    result = evaluate(*pair(contract=context), .01)
    assert np.isfinite(result.signed_stationary_gap_eV_fu)  # analytic H inputs, no new material pressure


@pytest.mark.parametrize('changes', [dict(scale_A_per_parameter_unit=True), dict(scale_A_per_parameter_unit=0.),
    dict(maximum_full_displacement_A=-1.), dict(source_audit_sha256='bad'), dict(parameter_interval=(.01,.02)),
    dict(parameter_interval=(-.01,)), dict(parameter_interval=(-.01,np.inf)), dict(parameter_interval=(.03,-.03))])
def test_invalid_local_domain(changes):
    with pytest.raises(ValueError): model(0, -1., 3., **changes)


@pytest.mark.parametrize('x', [True, np.nan, .04, -.04, .1j])
def test_invalid_target_parameter(x):
    with pytest.raises(ValueError): evaluate(*pair(), x)


def test_anchor_and_target_displacement_domains_and_indices_are_checked():
    initial, saddle = pair()
    with pytest.raises(ValueError, match='declared radius'):
        evaluate(replace(initial, maximum_full_displacement_A=.001), saddle, 0.)
    with pytest.raises(ValueError, match='declared radius'):
        evaluate(initial, replace(saddle, maximum_full_displacement_A=.1), .03)
    with pytest.raises(ValueError, match='index0'):
        evaluate(saddle, initial, .01, gap=-.8)
    with pytest.raises(ValueError, match='contradicts'):
        evaluate(initial, saddle, .01, gap=.9)
    with pytest.raises(ValueError): evaluate(initial, saddle, .01, tolerance=-.01)
    with pytest.raises(ValueError): evaluate(initial, saddle, .01, gap=True)


def test_nonfinite_chain_rule_is_rejected_without_numeric_repair():
    initial, saddle = pair()
    with pytest.raises(ValueError, match='nonfinite paired'):
        evaluate(initial, replace(saddle, scale_A_per_parameter_unit=1e200), 0.)


def test_actual_biaxial_centres_share_the_prescribed_substrate_not_local_cell_scale():
    from ase.build import bulk
    from vcneb import BiaxialClampedCurvatureCoordinates, clamped_plane_vcneb_boundary
    reference = bulk('Cu', 'fcc', a=3.6, cubic=True)
    common = contract(); models = []; charts = []
    for index, cellscale, strainscale in ((0, 2., 3.), (1, 4., 5.)):
        centre = reference.copy(); cell = centre.cell.array.copy(); cell[2,2] *= (1+.03*index)
        centre.set_cell(cell, scale_atoms=True)
        boundary = clamped_plane_vcneb_boundary(len(centre), reference.cell.array, allow_tilt=True)
        chart = BiaxialClampedCurvatureCoordinates(centre, boundary=boundary, cell_scale_A=cellscale,
                                                anchor_strain=0., strain_scale_A=strainscale)
        allowed = chart.internal_relaxation_basis(); c = chart.controlled_direction()
        h = allowed @ np.diag([(-2. if index else 2.)]+[4.]*(allowed.shape[1]-1)) @ allowed.T
        h += np.outer(allowed[:,0],c)+np.outer(c,allowed[:,0])+3.*np.outer(c,c)
        branch = stationary_quadratic_branch(h, .1*c, -10000.+.8*index, allowed[:,:1], allowed[:,1:], c,
            admissible_internal_basis=allowed, expected_index=index, stability_floor=.01, internal_curvature_floor=.01)
        models.append(ControlledStationaryModel(branch, common, strainscale, (-.03,.03), 2., 'd'*64))
        charts.append(chart)
    response = evaluate(*models, .015)
    for chart, point in zip(charts, (response.initial_point, response.bottleneck_point)):
        atoms = chart.displaced(point.full_displacement)
        np.testing.assert_allclose(atoms.cell.array[:2], reference.cell.array[:2]*1.015, atol=1e-14)


def test_analytic_benchmark_and_existing_output_guard(tmp_path, monkeypatch):
    from scripts import check_stationary_gap as checker
    result = checker.run_benchmark()
    assert result['groups'] == 8 and result['points'] == 24
    assert result['DFT_calls'] == result['calculator_calls'] == 0
    assert not result['B2_to_B5_material_forecast_frozen'] and not result['HfO2_advantage_proven']
    output = tmp_path/'keep.json'; output.write_text('original bytes')
    monkeypatch.setattr('sys.argv', ['checker','--output',str(output)])
    with pytest.raises(FileExistsError): checker.main()
    assert output.read_text() == 'original bytes'
