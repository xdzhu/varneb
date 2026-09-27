import numpy as np
import pytest

from scripts.analyze_gan_2d_joint_mode_feasibility import plane_coverage


def test_exact_two_mode_plane_has_zero_residual():
    path = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [-1.0, 1.0, 0.0]])
    modes = np.eye(3)[:, :2]
    result = plane_coverage(path, modes)
    assert result["maximum_off_plane_norm_A"] == pytest.approx(0.0)
    assert result["captured_squared_path_fraction"] == pytest.approx(1.0)


def test_off_plane_residual_and_fraction_are_explicit():
    path = np.array([[1.0, 0.0, 0.2], [0.0, 1.0, -0.1]])
    result = plane_coverage(path, np.eye(3)[:, :2])
    assert result["off_plane_norm_A"] == pytest.approx([0.2, 0.1])
    assert result["captured_squared_path_fraction"] == pytest.approx(2.0 / 2.05)


def test_nonorthonormal_axes_rejected():
    with pytest.raises(ValueError, match="nonorthonormal"):
        plane_coverage(np.ones((2, 3)), np.ones((3, 2)))
