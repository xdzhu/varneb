"""Calculator-free tests for Gamma-mode decomposition."""

from __future__ import annotations

import numpy as np
import pytest

from vcneb.phonons import (
    diagonalize_gamma_modes,
    mass_weighted_dynamical_matrix,
    project_displacements_onto_gamma_modes,
    tangent_mode_overlaps,
    translation_basis,
)


def test_mass_weighted_modes_and_signed_frequencies() -> None:
    force_constants = np.diag([4.0, 9.0, -16.0])
    modes = diagonalize_gamma_modes(force_constants, [4.0], project_translations=False)
    assert np.allclose(modes.eigenvalues_eV_per_A2_amu, [-4.0, 1.0, 2.25])
    assert np.allclose(modes.frequencies_cm1 / 521.4708986, [-2.0, 1.0, 1.5])


def test_translation_projection_removes_acoustic_component() -> None:
    masses = np.array([2.0, 2.0])
    force_constants = np.zeros((6, 6))
    force_constants[0, 0] = force_constants[3, 3] = 1.0
    force_constants[0, 3] = force_constants[3, 0] = -1.0
    modes = diagonalize_gamma_modes(force_constants, masses)
    translations = translation_basis(masses)
    assert np.allclose(translations.T @ translations, np.eye(3))
    assert np.count_nonzero(np.abs(modes.eigenvalues_eV_per_A2_amu) < 1e-12) == 5
    assert np.isclose(modes.eigenvalues_eV_per_A2_amu[-1], 1.0)


def test_known_mode_mixture_is_recovered() -> None:
    modes = diagonalize_gamma_modes(np.diag([1.0, 4.0, 9.0]), [1.0], project_translations=False)
    displacement = np.array([[0.25, -0.5, 0.75]])
    coordinates = project_displacements_onto_gamma_modes(displacement, modes)
    assert np.allclose(np.abs(coordinates), [0.25, 0.5, 0.75])
    path = np.array([coordinates * 0.0, coordinates * 0.5, coordinates])
    overlaps = tangent_mode_overlaps(path)
    assert overlaps.shape == (2, 3)
    assert np.allclose(overlaps[0], overlaps[1])


def test_force_constant_and_path_shapes_are_checked() -> None:
    with pytest.raises(ValueError):
        mass_weighted_dynamical_matrix(np.eye(4), [1.0])
    modes = diagonalize_gamma_modes(np.eye(3), [1.0], project_translations=False)
    with pytest.raises(ValueError):
        project_displacements_onto_gamma_modes(np.zeros((2, 2)), modes)
