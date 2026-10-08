"""Analytic/coordinate tests for the conditional-curvature foundation."""

import numpy as np
import pytest

from vcneb.relaxed_curvature import relax_orthogonal_curvature


def test_coupled_harmonic_minimum_and_saddle():
    # E=.5*a*q^2 + b*q*r + .5*c*r^2; r*=-b*q/c.
    for a in (5.0, -1.0):
        result = relax_orthogonal_curvature(
            np.array([[a, 2.0], [2.0, 4.0]]), np.eye(2)[:, :1],
            np.eye(2)[:, 1:], stability_floor=0.01,
        )
        assert result.relaxed[0, 0] == pytest.approx(a - 1.0)
        assert result.orthogonal_response[0, 0] == pytest.approx(-0.5)
        assert result.softening[0, 0] == pytest.approx(1.0)
        for q in (-0.3, 0.2):
            r = result.orthogonal_response[0, 0] * q
            full = np.array([q, r])
            assert full @ np.array([[a, 2.0], [2.0, 4.0]]) @ full == pytest.approx(
                result.relaxed[0, 0] * q * q
            )


def test_softening_is_psd_and_covariant_under_coordinate_rotation():
    rng = np.random.default_rng(23)
    factor = rng.normal(size=(7, 7))
    hessian = factor.T @ factor + np.eye(7)
    q, r = np.eye(7)[:, :2], np.eye(7)[:, 2:6]
    expected = relax_orthogonal_curvature(hessian, q, r, stability_floor=0.01)
    rotation, _ = np.linalg.qr(rng.normal(size=(7, 7)))
    actual = relax_orthogonal_curvature(
        rotation.T @ hessian @ rotation, rotation.T @ q, rotation.T @ r,
        stability_floor=0.01,
    )
    np.testing.assert_allclose(actual.relaxed, expected.relaxed, atol=1e-12)
    np.testing.assert_allclose(actual.orthogonal_response, expected.orthogonal_response, atol=1e-12)
    assert np.linalg.eigvalsh(expected.softening).min() >= -1e-12


def test_omitted_directions_remain_clamped():
    hessian = np.array([[4.0, 1.0, 2.0], [1.0, 2.0, 0.0], [2.0, 0.0, 3.0]])
    partial = relax_orthogonal_curvature(
        hessian, np.eye(3)[:, :1], np.eye(3)[:, 1:2], stability_floor=0.01,
    )
    assert partial.relaxed[0, 0] == pytest.approx(3.5)
    frozen = relax_orthogonal_curvature(
        hessian, np.eye(3)[:, :1], np.empty((3, 0)), stability_floor=0.01,
    )
    np.testing.assert_array_equal(frozen.relaxed, [[4.0]])
    assert frozen.minimum_eliminated_curvature is None


@pytest.mark.parametrize("value", [-1.0, 0.0, 0.0001, 0.01])
def test_unstable_or_unresolved_elimination_is_not_pseudo_inverted(value):
    with pytest.raises(ValueError, match="unstable or unresolved"):
        relax_orthogonal_curvature(
            np.diag([2.0, value]), np.eye(2)[:, :1], np.eye(2)[:, 1:],
            stability_floor=0.01,
        )


@pytest.mark.parametrize("problem", ["nonfinite", "asymmetric", "overlap", "scale", "dimension", "floor"])
def test_invalid_contract_is_rejected(problem):
    h, q, r, floor = np.eye(3), np.eye(3)[:, :1], np.eye(3)[:, 1:], 0.01
    if problem == "nonfinite":
        h[0, 0] = np.nan
    elif problem == "asymmetric":
        h[0, 1] = 1
    elif problem == "overlap":
        r = q
    elif problem == "scale":
        q = q * 2
    elif problem == "dimension":
        r = np.ones((4, 1))
    else:
        floor = -0.1
    with pytest.raises(ValueError):
        relax_orthogonal_curvature(h, q, r, stability_floor=floor)
