"""Path-SVD directions are translation-invariant empirical axes, not phonons."""

import numpy as np
import pytest

from vcneb.path_modes import fit_path_atomic_svd, leave_one_image_out_residuals


def _two_mode_path():
    masses = np.array([1.0, 4.0, 9.0])
    one = np.array([[1.0, 0.0, 0.0], [-0.25, 0.0, 0.0], [0.0, 0.0, 0.0]])
    two = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, -1.0 / 9.0, 0.0]])
    amplitudes = np.array([[0.0, 0.0], [1.0, 0.5], [2.0, -1.0], [1.5, 2.0], [3.0, -0.5]])
    path = amplitudes[:, 0, None, None] * one + amplitudes[:, 1, None, None] * two
    return path, masses


def test_empirical_two_mode_fit_and_holdout_are_exact():
    path, masses = _two_mode_path()
    fit = fit_path_atomic_svd(path, masses, n_components=2)
    assert fit.captured_squared_norm_fraction == pytest.approx(1.0, abs=1e-12)
    assert np.max(fit.residuals_sqrt_amu_A) < 1e-12
    assert np.max(leave_one_image_out_residuals(path, masses, n_components=2)) < 1e-12
    np.testing.assert_allclose(fit.basis_mass_weighted.T @ fit.basis_mass_weighted, np.eye(2), atol=1e-12)


def test_independent_rigid_translations_cannot_change_internal_basis():
    path, masses = _two_mode_path()
    shifted = path + np.arange(len(path))[:, None, None] * np.array([0.2, -0.3, 0.4])
    original = fit_path_atomic_svd(path, masses, n_components=2)
    moved = fit_path_atomic_svd(shifted, masses, n_components=2)
    np.testing.assert_allclose(moved.basis_mass_weighted, original.basis_mass_weighted, atol=1e-12)
    np.testing.assert_allclose(moved.coefficients_sqrt_amu_A, original.coefficients_sqrt_amu_A, atol=1e-12)
    assert np.max(moved.residuals_sqrt_amu_A) < 1e-12


def test_rank_deficiency_and_invalid_inputs_are_rejected():
    path, masses = _two_mode_path()
    with pytest.raises(ValueError, match="rank deficient"):
        fit_path_atomic_svd(path, masses, n_components=3)
    with pytest.raises(ValueError, match="positive finite"):
        fit_path_atomic_svd(path, [1.0, 0.0, 9.0], n_components=1)
    with pytest.raises(ValueError, match="n_components"):
        fit_path_atomic_svd(path, masses, n_components=20)
