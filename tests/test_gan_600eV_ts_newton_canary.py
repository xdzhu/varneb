import numpy as np
import pytest

from scripts.prepare_gan_600eV_ts_newton_canary import basal_projection_A


def test_orthogonal_cell_has_zero_basal_projection():
    assert basal_projection_A(np.diag([3.0, 4.0, 5.0])) == pytest.approx(0.0)


def test_tiny_third_vector_shear_is_measured_in_angstrom():
    cell = np.diag([3.0, 4.0, 5.0])
    cell[2, 0] = 2e-4
    cell[2, 1] = -4e-5
    assert basal_projection_A(cell) == pytest.approx(2e-4)


def test_invalid_cell_rejected():
    with pytest.raises(ValueError, match="invalid cell"):
        basal_projection_A(np.ones((2, 3)))
