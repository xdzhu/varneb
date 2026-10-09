"""Fixed control, saddle sign, nonstationary offsets and explicit restrictions."""

import numpy as np
import pytest

from vcneb import stationary_quadratic_branch


def example():
    h = np.array([
        [-2., .1, .3, 0., .5, .2], [.1, 3., .2, .1, -.2, 0.],
        [.3, .2, 4., .3, .4, 0.], [0., .1, .3, 5., .2, 0.],
        [.5, -.2, .4, .2, -4., 0.], [.2, 0., 0., 0., 0., -6.],
    ])
    g = np.array([.12, -.11, .08, .22, .04, 3.])
    basis = np.eye(6)
    return h, g, basis[:, :2], basis[:, 2:4], basis[:, 4]


def build(h, g, q, r, c, admissible_internal_basis=None, **kwargs):
    allowed = np.column_stack((q, r)) if admissible_internal_basis is None else admissible_internal_basis
    return stationary_quadratic_branch(h, g, -10000., q, r, c,
                                      admissible_internal_basis=allowed,
                                      expected_index=1, stability_floor=.01,
                                      internal_curvature_floor=.01, **kwargs)


def test_matches_independent_fixed_control_stationary_linear_solve():
    h, g, q, r, c = example()
    branch = build(h, g, q, r, c, admissible_internal_basis=np.eye(6)[:, [0, 1, 2, 3, 5]])
    internal = np.column_stack((q, r))
    for t in (-.15, 0., .2):
        point = branch.evaluate(t)
        coefficients = -np.linalg.solve(internal.T @ h @ internal, internal.T @ (g + h @ c * t))
        expected = internal @ coefficients + c * t
        np.testing.assert_allclose(point.full_displacement, expected, atol=1e-14)
        np.testing.assert_allclose(point.retained_internal_coordinates, coefficients[:2], atol=1e-14)
        np.testing.assert_allclose(point.eliminated_coordinates, coefficients[2:], atol=1e-14)
        assert point.energy_change == pytest.approx(g @ expected + .5 * expected @ h @ expected, abs=1e-14)
        assert point.control_gradient == pytest.approx(c @ (g + h @ expected), abs=1e-14)
        assert point.declared_internal_gradient_norm < 1e-14
        assert point.full_displacement[4] == t  # external coordinate stays prescribed
        assert point.full_displacement[5] == 0  # unmeasured direction stays frozen
        assert point.unrepresented_internal_gradient_norm > 2.9  # not full physical stationarity
        assert point.clamped_gradient_norm < 1e-14
        dt = 1e-5
        derivative = (branch.evaluate(t+dt).energy_change-branch.evaluate(t-dt).energy_change)/(2*dt)
        curvature = (branch.evaluate(t+dt).control_gradient-branch.evaluate(t-dt).control_gradient)/(2*dt)
        assert derivative == pytest.approx(point.control_gradient, abs=1e-10)
        assert curvature == pytest.approx(branch.control_curvature, abs=1e-10)


def test_saddle_control_response_can_harden_unlike_a_minimum():
    q, r, c = np.eye(3)[:, :1], np.eye(3)[:, 1:2], np.eye(3)[:, 2]
    saddle = np.array([[-2., 0., 1.], [0., 4., 0.], [1., 0., -3.]])
    minimum = saddle.copy(); minimum[0, 0] = 2.
    ts = build(saddle, np.zeros(3), q, r, c)
    well = stationary_quadratic_branch(minimum, np.zeros(3), 0., q, r, c,
                                      admissible_internal_basis=np.eye(3)[:, :2],
                                      expected_index=0, stability_floor=.01, internal_curvature_floor=.01)
    assert np.count_nonzero(np.linalg.eigvalsh(saddle)<0) == 2
    assert ts.expected_index == 1  # controlled negative curvature is not an internal instability
    assert ts.control_curvature == pytest.approx(-2.5)
    assert well.control_curvature == pytest.approx(-3.5)
    assert ts.control_curvature > saddle[-1, -1]
    assert well.control_curvature < minimum[-1, -1]


def test_rotation_covariance_including_retained_internal_frame():
    h, g, q, r, c = example(); original = build(h, g, q, r, c)
    rng = np.random.default_rng(512)
    u, _ = np.linalg.qr(rng.normal(size=h.shape))
    v, _ = np.linalg.qr(rng.normal(size=(2, 2)))
    rotated = build(u.T @ h @ u, u.T @ g, u.T @ q @ v, u.T @ r, u.T @ c)
    for t in (-.07, .1):
        a, b = original.evaluate(t), rotated.evaluate(t)
        np.testing.assert_allclose(b.full_displacement, u.T @ a.full_displacement, atol=1e-14)
        np.testing.assert_allclose(b.retained_internal_coordinates, v.T @ a.retained_internal_coordinates, atol=1e-14)
        assert b.energy_change == pytest.approx(a.energy_change, abs=1e-14)
        assert b.control_gradient == pytest.approx(a.control_gradient, abs=1e-14)
        assert b.unrepresented_internal_gradient_norm == pytest.approx(a.unrepresented_internal_gradient_norm)
        assert b.clamped_gradient_norm == pytest.approx(a.clamped_gradient_norm)
    assert rotated.control_curvature == pytest.approx(original.control_curvature)


def test_frozen_release_and_reference_energy_do_not_hide_anchor_offsets():
    h, g, q, _, c = example()
    branch = build(h, g, q, np.empty((6, 0)), c)
    other = stationary_quadratic_branch(h, g, 1e9, q, np.empty((6, 0)), c,
                                       admissible_internal_basis=q,
                                       expected_index=1, stability_floor=.01, internal_curvature_floor=.01)
    for t in (0., .01):
        point = branch.evaluate(t)
        assert point.eliminated_coordinates.size == 0
        np.testing.assert_array_equal(point.full_displacement[2:4], [0., 0.])
        assert point.energy_change == other.evaluate(t).energy_change
    assert branch.evaluate(0.).energy_change != 0.
    assert np.linalg.norm(branch.internal_offset) > 0.


def test_owned_immutable_arrays_and_input_changes():
    args = example(); branch = build(*args); expected = branch.evaluate(.02)
    for array in args: array[:] = 0.
    np.testing.assert_array_equal(branch.evaluate(.02).full_displacement, expected.full_displacement)
    for array in (branch.hessian, branch.gradient, branch.internal_offset, branch.internal_response,
                  branch.retained_internal_eigenvalues, expected.full_displacement,
                  branch.admissible_internal_basis,
                  expected.retained_internal_coordinates, expected.eliminated_coordinates):
        with pytest.raises(ValueError): array.flat[0] = 9.


@pytest.mark.parametrize('problem', ['unstable_release','unresolved_release','minimum_not_saddle',
                                   'two_negative','unresolved_internal','control_in_Q','control_in_R',
                                   'nonunit_control','complex_H','bad_control_shape','no_internal'])
def test_invalid_stationary_or_elimination_space_is_rejected(problem):
    h, g, q, r, c = example()
    if problem == 'unstable_release': h[2, 2] = -4.
    elif problem == 'unresolved_release': h[2:4, 2:4] = np.diag([.01, 5.])
    elif problem == 'minimum_not_saddle': h[0, 0] = 2.
    elif problem == 'two_negative': h[1, 1] = -3.
    elif problem == 'unresolved_internal':
        h[:2, :] = 0.; h[:, :2] = 0.; h[:2, :2] = np.diag([-.01, 3.])
    elif problem == 'control_in_Q': c = q[:, 0]
    elif problem == 'control_in_R': c = r[:, 0]
    elif problem == 'nonunit_control': c *= 2.
    elif problem == 'complex_H': h = h.astype(complex); h[0, 0] += 1j
    elif problem == 'bad_control_shape': c = c[:, None]
    elif problem == 'no_internal': q = np.empty((6, 0))
    with pytest.raises(ValueError): build(h, g, q, r, c)


@pytest.mark.parametrize('index', [True, 1., 2, -1])
def test_invalid_requested_index(index):
    with pytest.raises(ValueError):
        stationary_quadratic_branch(*example()[:2], 0., *example()[2:],
                                    admissible_internal_basis=np.eye(6)[:, [0,1,2,3,5]],
                                    expected_index=index, stability_floor=.01, internal_curvature_floor=.01)


@pytest.mark.parametrize('value', [True, np.nan, np.inf, [.1], .1j])
def test_invalid_control_value(value):
    with pytest.raises(ValueError): build(*example()).evaluate(value)


@pytest.mark.parametrize('name,value', [('stability_floor',True),('internal_curvature_floor',-.1),
                                      ('internal_curvature_floor',np.nan),('symmetry_tolerance',0.)])
def test_invalid_numeric_policy(name, value):
    h, g, q, r, c = example()
    policy=dict(expected_index=1,stability_floor=.01,internal_curvature_floor=.01,symmetry_tolerance=1e-8)
    policy[name]=value
    with pytest.raises(ValueError):
        stationary_quadratic_branch(h,g,0.,q,r,c,admissible_internal_basis=np.column_stack((q,r)),**policy)


@pytest.mark.parametrize('problem', ['missing_Q','missing_R','control_admissible','nonorthogonal'])
def test_actual_mechanical_space_is_required(problem):
    h,g,q,r,c=example()
    allowed=np.column_stack((q,r))
    if problem=='missing_Q': allowed=allowed[:,1:]
    elif problem=='missing_R': allowed=allowed[:,:3]
    elif problem=='control_admissible': allowed=np.column_stack((allowed,c))
    else: allowed[:,0]*=2.
    with pytest.raises(ValueError): build(h,g,q,r,c,admissible_internal_basis=allowed)


def test_clamped_reactions_are_not_unrepresented_internal_forces():
    h,g,q,r,c=example()
    point=build(h,g,q,r,c).evaluate(.03)
    assert point.declared_internal_gradient_norm < 1e-14
    assert point.unrepresented_internal_gradient_norm < 1e-14
    assert point.clamped_gradient_norm > 2.9


def test_actual_biaxial_chart_keeps_the_prescribed_substrate():
    from ase.build import bulk
    from vcneb import BiaxialClampedCurvatureCoordinates, clamped_plane_vcneb_boundary
    atoms=bulk('Cu','fcc',a=3.6,cubic=True)
    boundary=clamped_plane_vcneb_boundary(len(atoms),atoms.cell.array,allow_tilt=True)
    chart=BiaxialClampedCurvatureCoordinates(atoms,boundary=boundary,cell_scale_A=2.,
                                           anchor_strain=.01,strain_scale_A=3.)
    allowed=chart.internal_relaxation_basis(); c=chart.controlled_direction()
    h=allowed @ np.diag([-2.]+[4.]*(allowed.shape[1]-1)) @ allowed.T
    h+=np.outer(allowed[:,0],c)+np.outer(c,allowed[:,0])-3.*np.outer(c,c)
    branch=stationary_quadratic_branch(h,.1*c,0.,allowed[:,:1],allowed[:,1:],c,
                                      admissible_internal_basis=allowed,expected_index=1,
                                      stability_floor=.01,internal_curvature_floor=.01)
    point=branch.evaluate(.03); result=chart.displaced(point.full_displacement)
    np.testing.assert_allclose(result.cell.array[:2],atoms.cell.array[:2]*1.02/1.01,atol=1e-14)
    assert point.full_displacement[chart.controlled_index] == .03
    assert point.declared_internal_gradient_norm < 1e-13


def test_external_scale_chain_rule_does_not_change_stationary_energy():
    h,g,q,r,c=example(); original=build(h,g,q,r,c)
    scale=np.eye(6); scale[4,4]=1/3.
    rescaled=build(scale.T @ h @ scale,scale.T @ g,q,r,c)
    for t in (-.08,.12):
        a,b=original.evaluate(t),rescaled.evaluate(3*t)
        assert b.energy_change == pytest.approx(a.energy_change,abs=1e-14)
        np.testing.assert_allclose(scale @ b.full_displacement,a.full_displacement,atol=1e-14)
        assert 3*b.control_gradient == pytest.approx(a.control_gradient,abs=1e-14)
    assert 9*rescaled.control_curvature == pytest.approx(original.control_curvature)


def test_bounded_benchmark_is_analytic_not_material_evidence():
    from scripts.check_stationary_branch import run_benchmark
    result=run_benchmark()
    assert result['groups']==8 and result['points']==24
    assert result['DFT_calls']==0 and result['calculator_calls']==0
    assert not result['HfO2_curvature_or_independent_predictions_measured']
    assert result['maximum_errors']['internal_gradient'] < 1e-12
    assert all(len(digest)==64 for digest in result['source_sha256'].values())


def test_benchmark_refuses_existing_output_before_computation(tmp_path,monkeypatch):
    import scripts.check_stationary_branch as checker
    out=tmp_path/'keep.json';out.write_text('original bytes')
    monkeypatch.setattr('sys.argv',['checker','--output',str(out)])
    with pytest.raises(FileExistsError):checker.main()
    assert out.read_text()=='original bytes'
