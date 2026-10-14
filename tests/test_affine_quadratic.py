"""Standard constrained quadratic identities; all fixtures are analytic, not DFT."""
import numpy as np
import pytest

from vcneb.affine_quadratic import constrain_quadratic
from vcneb.quadratic_reduction import condition_quadratic_energy


@pytest.mark.parametrize('seed', range(6))
def test_full_reference_kkt_equals_adaptive_schur_with_affine_gradient(seed):
    rng = np.random.default_rng(seed)
    rot = np.linalg.qr(rng.normal(size=(7, 7)))[0]
    h = rot @ np.diag([-2., -1., .5, 1., 2., 3., 4.]) @ rot.T
    g = rng.normal(size=7)
    q, r = rot[:, :2], rot[:, 2:]
    schur = condition_quadratic_energy(h, g, -100., q, r, stability_floor=0.)
    full = constrain_quadratic(h, g, np.eye(7), q, stability_floor=0.)
    another_basis = constrain_quadratic(h, g, rot, q, stability_floor=0.)
    for coords in [np.zeros(2), np.array([-.1, .2]), np.array([.3, -.2])]:
        np.testing.assert_allclose(full.coordinates(coords), schur.full_displacement(coords), atol=2e-12, rtol=0)
        assert full.energy_change(coords) == pytest.approx(schur.energy_change(coords), abs=2e-12, rel=0)
        np.testing.assert_allclose(full.retained_gradient(coords), schur.retained_gradient(coords), atol=2e-12, rtol=0)
        np.testing.assert_allclose(full.coordinates(coords), another_basis.coordinates(coords), atol=2e-12, rtol=0)
        assert max(full.residuals(coords).values()) < 2e-12


def test_mixed_observables_keep_their_declared_scale_and_gradient_sign():
    h = np.diag([2., 3., 4.])
    g = np.array([.2, -.1, .3])
    a = np.array([[2., 1.], [0., 1.], [1., 0.]])
    model = constrain_quadratic(h, g, np.eye(3), a, stability_floor=0.)
    q = np.array([.1, -.2])
    np.testing.assert_allclose(a.T @ model.coordinates(q), q, atol=2e-12, rtol=0)
    numerical = []
    for i in range(2):
        delta = np.eye(2)[i] * 1e-6
        numerical.append((model.energy_change(q+delta) - model.energy_change(q-delta)) / 2e-6)
    np.testing.assert_allclose(model.retained_gradient(q), numerical, atol=2e-11, rtol=0)


def test_clamped_outside_reference_is_not_silently_released():
    h = np.eye(3)
    g = np.array([0., 1., 2.])
    a = np.eye(3)[:, :1]
    restricted = constrain_quadratic(h, g, np.eye(3)[:, :2], a, stability_floor=0.)
    full = constrain_quadratic(h, g, np.eye(3), a, stability_floor=0.)
    assert restricted.coordinates([.1])[2] == 0.
    assert full.energy_change([.1]) < restricted.energy_change([.1])


def test_empty_free_space_is_unique_constraint_not_physical_stability():
    model = constrain_quadratic([[-2.]], [.3], [[1.]], [[2.]], stability_floor=0.)
    assert model.minimum_free_curvature is None
    np.testing.assert_array_equal(model.coordinates([.4]), [.2])


@pytest.mark.parametrize('curvature', [-1., 0., .1])
def test_nonpositive_or_unresolved_free_block_rejected_without_clipping(curvature):
    with pytest.raises(ValueError, match='unstable or unresolved'):
        constrain_quadratic(np.diag([-2., curvature]), [0., .1], np.eye(2), np.eye(2)[:, :1],
                            stability_floor=.1)


@pytest.mark.parametrize('fault', ['complex', 'rank', 'metric', 'symmetry', 'nan'])
def test_invalid_reference_constraint_or_input_is_not_repaired(fault):
    h, g, u, a = np.eye(3), np.ones(3), np.eye(3)[:, :2], np.eye(3)[:, :1]
    if fault == 'complex': h = h.astype(complex)
    elif fault == 'rank': a = np.eye(3)[:, 2:]
    elif fault == 'metric': u *= 2
    elif fault == 'symmetry': h[0, 1] = .1
    else: g[0] = np.nan
    with pytest.raises(ValueError): constrain_quadratic(h, g, u, a, stability_floor=0.)


def test_model_arrays_do_not_alias_caller_or_allow_mutation():
    h, g = np.eye(2), np.array([.1, .2])
    model = constrain_quadratic(h, g, np.eye(2), np.eye(2)[:, :1], stability_floor=0.)
    old = model.energy_change([.3])
    h *= 10
    g[:] = 5
    assert model.energy_change([.3]) == old
    with pytest.raises(ValueError): model.offset[0] = 5.


def test_wrong_or_nonfinite_retained_coordinates_refused():
    model = constrain_quadratic(np.eye(2), [0., 0.], np.eye(2), np.eye(2)[:, :1], stability_floor=0.)
    for q in [[], [0., 1.], [float('nan')], [1.j]]:
        with pytest.raises(ValueError): model.energy_change(q)
