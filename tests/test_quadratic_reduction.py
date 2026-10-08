"""Nonstationary offsets, gradients and restricted release are explicit."""

import numpy as np
import pytest

from vcneb.quadratic_reduction import condition_quadratic_energy


def test_nonstationary_offset_energy_and_envelope_gradient():
    h = np.array([[5., 2.], [2., 4.]])
    g = np.array([.3, -.8])
    model = condition_quadratic_energy(h, g, -10000., np.eye(2)[:, :1],
                                      np.eye(2)[:, 1:], stability_floor=.01)
    assert model.eliminated_offset[0] == pytest.approx(.2)
    assert model.reference_energy_change == pytest.approx(-.08)
    assert model.retained_gradient_at_zero[0] == pytest.approx(.7)
    assert model.curvature.relaxed[0, 0] == pytest.approx(4.)
    for value in (-.3, 0., .2):
        q = np.array([value])
        y = model.full_displacement(q)
        assert (g + h @ y)[1] == pytest.approx(0., abs=1e-14)
        assert model.energy_change(q) == pytest.approx(g @ y + .5 * y @ h @ y)
        assert model.retained_gradient(q)[0] == pytest.approx((g + h @ y)[0])
        dq = np.array([1e-6])
        derivative = (model.energy_change(q + dq) - model.energy_change(q - dq)) / 2e-6
        assert derivative == pytest.approx(model.retained_gradient(q)[0], abs=1e-9)


def test_covariance_multidimensional_release_and_omitted_coordinates():
    rng = np.random.default_rng(82)
    a = rng.normal(size=(7, 7))
    h = a.T @ a + np.eye(7)
    g = rng.normal(size=7)
    q, r = np.eye(7)[:, :2], np.eye(7)[:, 2:5]
    expected = condition_quadratic_energy(h, g, 2., q, r, stability_floor=.01)
    rotation, _ = np.linalg.qr(rng.normal(size=(7, 7)))
    actual = condition_quadratic_energy(rotation.T @ h @ rotation, rotation.T @ g,
                                      2., rotation.T @ q, rotation.T @ r, stability_floor=.01)
    values = np.array([.2, -.1])
    np.testing.assert_allclose(actual.eliminated_coordinates(values), expected.eliminated_coordinates(values), atol=1e-13)
    np.testing.assert_allclose(actual.full_displacement(values), rotation.T @ expected.full_displacement(values), atol=1e-13)
    np.testing.assert_allclose(actual.retained_gradient(values), expected.retained_gradient(values), atol=1e-13)
    assert actual.energy_change(values) == pytest.approx(expected.energy_change(values))
    np.testing.assert_array_equal(expected.full_displacement(values)[5:], np.zeros(2))
    np.testing.assert_allclose(r.T @ (g + h @ expected.full_displacement(values)), 0., atol=1e-13)


def test_no_elimination_and_reference_energy_shift():
    h, g, q = np.diag([2., -1., -3.]), np.array([.4, -.7, 10.]), np.eye(3)[:, :2]
    model = condition_quadratic_energy(h, g, -5., q, np.empty((3, 0)), stability_floor=.01)
    shifted = condition_quadratic_energy(h, g, 500., q, np.empty((3, 0)), stability_floor=.01)
    value = np.array([.3, -.2])
    assert model.reference_energy_change == 0.
    assert model.eliminated_coordinates(value).shape == (0,)
    assert model.full_displacement(value)[2] == 0.
    assert shifted.reference_energy == 500.
    assert shifted.energy_change(value) == model.energy_change(value)
    assert model.energy_change(value) == pytest.approx(g[:2] @ value + .5 * value @ h[:2, :2] @ value)


def test_stationary_reference_recovers_curvature_only_model():
    model = condition_quadratic_energy(np.array([[-1., 2.], [2., 4.]]), np.zeros(2),
                                      10., np.eye(2)[:, :1], np.eye(2)[:, 1:], stability_floor=.01)
    np.testing.assert_array_equal(model.eliminated_offset, [0.])
    np.testing.assert_array_equal(model.retained_gradient_at_zero, [0.])
    assert model.reference_energy_change == 0.
    assert model.energy_change(np.array([.2])) == pytest.approx(-.04)


def test_model_arrays_are_copies_and_readonly():
    h, g, q, r = np.eye(2), np.array([1., 2.]), np.eye(2)[:, :1], np.eye(2)[:, 1:]
    model = condition_quadratic_energy(h, g, 0., q, r, stability_floor=.01)
    g[:] = 9.; q[:] = 0.; r[:] = 0.; h[:] = 5.
    assert model.energy_change(np.zeros(1)) == -2.
    for array in (model.retained_basis, model.eliminated_basis, model.eliminated_offset,
                  model.retained_gradient_at_zero, model.curvature.relaxed,
                  model.curvature.orthogonal_response):
        with pytest.raises(ValueError):
            array.flat[0] = 3.


@pytest.mark.parametrize("problem", ["gradient_size", "gradient_nan", "energy_nan", "unstable", "unresolved", "asymmetric", "overlap"])
def test_invalid_reference_or_release_is_rejected(problem):
    h, g, q, r, e = np.eye(2), np.ones(2), np.eye(2)[:, :1], np.eye(2)[:, 1:], 0.
    if problem == "gradient_size":
        g = np.ones(3)
    elif problem == "gradient_nan":
        g[0] = np.nan
    elif problem == "energy_nan":
        e = np.nan
    elif problem == "unstable":
        h[1, 1] = -1.
    elif problem == "unresolved":
        h[1, 1] = .01
    elif problem == "asymmetric":
        h[0, 1] = 1.
    else:
        r = q
    with pytest.raises(ValueError):
        condition_quadratic_energy(h, g, e, q, r, stability_floor=.01)


@pytest.mark.parametrize("value", [0., [1., 2.], [float("nan")]])
def test_bad_evaluation_coordinates_are_rejected(value):
    model = condition_quadratic_energy(np.eye(2), np.ones(2), 0., np.eye(2)[:, :1],
                                      np.eye(2)[:, 1:], stability_floor=.01)
    for method in (model.energy_change, model.retained_gradient, model.eliminated_coordinates, model.full_displacement):
        with pytest.raises(ValueError):
            method(value)
