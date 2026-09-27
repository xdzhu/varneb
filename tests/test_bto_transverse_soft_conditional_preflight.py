"""Geometry of the third cubic soft-mode branch in a fixed-Q soft/soft plane."""

import numpy as np
import pytest

from examples.preflight_bto_transverse_soft_conditional import remaining_soft_y_direction
from vcneb.mode_surface import ModePlane
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
