import pytest

from scripts.audit_gan_600eV_ts_2d_pilot import centered_curvature
from scripts.prepare_gan_600eV_ts_2d_pilot import sampling_coordinates


def test_2d_pilot_has_eight_grid_points_and_four_holdouts():
    points = sampling_coordinates()
    assert len(points) == 12
    assert len({name for name, _, _ in points}) == 12
    assert sum(name.startswith("grid_") for name, _, _ in points) == 8
    assert sum(name.startswith("hold_") for name, _, _ in points) == 4


def test_centered_curvature_preserves_sign():
    assert centered_curvature(-0.0002, 0, -0.0002, 0.01) == pytest.approx(-4.0)
    assert centered_curvature(0.0001, 0, 0.0001, 0.01) == pytest.approx(2.0)


def test_centered_curvature_rejects_zero_step():
    with pytest.raises(ValueError, match="nonpositive step"):
        centered_curvature(1, 0, 1, 0)
