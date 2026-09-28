"""Geometry of the third cubic soft-mode branch in a fixed-Q soft/soft plane."""

import numpy as np
import pytest

from examples.preflight_bto_transverse_soft_conditional import remaining_soft_y_direction
from vcneb.mode_surface import ModePlane, _orthogonal_directions, relax_orthogonal_at_q
from vcneb.phonons import GammaModes


@pytest.mark.parametrize("q_transverse", [0.0, 0.3])
def test_remaining_soft_y_is_open_and_mass_metric_unit(q_transverse: float) -> None:
    order = [5, 3, 4] + [index for index in range(15) if index not in (5, 3, 4)]
    modes = GammaModes(
        masses_amu=np.ones(5),
        eigenvalues_eV_per_A2_amu=np.r_[-1.0, -1.0, -1.0, np.ones(12)],
        frequencies_cm1=np.r_[-1.0, -1.0, -1.0, np.ones(12)],
        eigenvectors=np.eye(15)[:, order],
        translations_projected=True,
    )
    axes = np.zeros((15, 2))
    axes[0, 0], axes[1, 1] = 1.0, 1.0
    plane = ModePlane.from_gamma_modes_with_strain(
        modes, axes, strain_metric_weights_amu_A2=np.ones(6),
        axis_labels=("z soft", "x soft"), reference_id="synthetic",
    )
    direction = remaining_soft_y_direction(modes, plane)
    assert direction.shape == (21,)
    assert direction[4] == pytest.approx(1.0)
    assert np.count_nonzero(direction[15:]) == 0
    assert np.dot(direction * plane.metric_weights, direction) == pytest.approx(1.0)
    for sign in (-1.0, 1.0):
        target_q = [0.6, q_transverse]
        assert np.allclose(plane.project(plane.frozen_coordinates(target_q) + sign * 0.4 * direction),
                           target_q)


def test_remaining_soft_y_rejects_atomic_only_plane() -> None:
    modes = GammaModes(
        masses_amu=np.ones(5), eigenvalues_eV_per_A2_amu=np.ones(15),
        frequencies_cm1=np.ones(15), eigenvectors=np.eye(15),
        translations_projected=False,
    )
    axes = np.eye(15, 2)
    plane = ModePlane.from_gamma_modes(modes, axes, axis_labels=("a", "b"), reference_id="synthetic")
    with pytest.raises(ValueError, match="six open strain"):
        remaining_soft_y_direction(modes, plane)


def test_restricted_soft_sheet_fixes_qy_without_changing_qz_qx() -> None:
    order = [5, 3, 4] + [index for index in range(15) if index not in (5, 3, 4)]
    modes = GammaModes(
        masses_amu=np.ones(5),
        eigenvalues_eV_per_A2_amu=np.r_[-1.0, -1.0, -1.0, np.ones(12)],
        frequencies_cm1=np.r_[-1.0, -1.0, -1.0, np.ones(12)],
        eigenvectors=np.eye(15)[:, order],
        translations_projected=True,
    )
    axes = np.zeros((15, 2))
    axes[0, 0], axes[1, 1] = 1.0, 1.0
    plane = ModePlane.from_gamma_modes_with_strain(
        modes, axes, strain_metric_weights_amu_A2=np.ones(6),
        axis_labels=("z soft", "x soft"), reference_id="synthetic",
    )
    third = remaining_soft_y_direction(modes, plane)
    translations = np.eye(21)[:, :3]
    frozen = np.column_stack([translations, third])
    assert _orthogonal_directions(plane, translations).shape == (21, 16)
    assert _orthogonal_directions(plane, frozen).shape == (21, 15)
    q = np.array([0.6, 0.3])
    base = plane.frozen_coordinates(q)

    def energy_and_gradient(coordinates: np.ndarray) -> tuple[float, np.ndarray]:
        delta = coordinates - base
        y = float(third @ delta)
        energy = float(0.5 * np.dot(delta, delta) - y**2 + 0.25 * y**4)
        gradient = delta + (-2.0 * y + y**3) * third
        return energy, gradient

    unrestricted = relax_orthogonal_at_q(
        plane, q, energy_and_gradient, starts=[base + 0.5 * third],
        frozen_directions=translations, gradient_tolerance=1e-6,
    )
    restricted = relax_orthogonal_at_q(
        plane, q, energy_and_gradient, frozen_directions=frozen,
        gradient_tolerance=1e-6,
    )
    assert abs(float(third @ unrestricted.coordinates)) > 0.5
    assert abs(float(third @ restricted.coordinates)) < 1e-10
    assert np.allclose(plane.project(restricted.coordinates), q)
    assert restricted.energy > unrestricted.energy
