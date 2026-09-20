"""Calculator-free tests for path geometry and modal contribution helpers."""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms

from vcneb.analysis import (
    dominant_mode_indices,
    mode_contribution_fractions,
    path_reaction_coordinate,
)


def test_path_reaction_coordinate_unwraps_periodic_crossing_and_cell_change() -> None:
    first = Atoms("H", scaled_positions=[[0.95, 0.0, 0.0]], cell=np.eye(3) * 4, pbc=True)
    second = Atoms("H", scaled_positions=[[1.05, 0.0, 0.0]], cell=np.eye(3) * 4.2, pbc=True)
    coordinate, segments = path_reaction_coordinate([first, second])
    assert coordinate.shape == (2,)
    assert segments.shape == (1,)
    assert coordinate[-1] > 0.0
    assert np.isclose(coordinate[-1], segments[0])
    assert coordinate[-1] < 1.0


def test_mode_contribution_fractions_and_dominant_modes() -> None:
    values = np.array([[0.0, 0.0, 0.0], [2.0, -1.0, 0.0]])
    fractions = mode_contribution_fractions(values)
    assert np.allclose(fractions.sum(axis=1), [0.0, 1.0])
    assert dominant_mode_indices(values, top=2) == [[0, 1], [0, 1]]
    with pytest.raises(ValueError):
        mode_contribution_fractions(np.zeros(3))
