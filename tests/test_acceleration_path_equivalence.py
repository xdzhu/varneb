"""Tests for common-arc comparison of first-threshold benchmark chains."""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from scripts.audit_acceleration_path_equivalence import audit_pair


def _chain(interior_x: float, *, final_x: float = 0.2) -> list[Atoms]:
    images = []
    for x in (0.0, interior_x, final_x):
        image = Atoms("H", positions=[[x, 0.0, 0.0]], cell=np.eye(3) * 3.0, pbc=True)
        image.calc = SinglePointCalculator(image, energy=-1.0 + x)
        images.append(image)
    return images


def test_common_arc_removes_image_redistribution() -> None:
    result = audit_pair(_chain(0.1), _chain(0.05), n_formula=1, pressure_gpa=0.0)
    assert result["max_same_index_joint_geometry_difference_A"] == pytest.approx(0.05)
    assert result["max_common_arc_joint_geometry_difference_A"] < 1e-12
    assert result["max_abs_common_arc_energy_difference_meV_per_formula"] < 1e-9


def test_common_arc_rejects_different_endpoints() -> None:
    with pytest.raises(ValueError, match="ordered endpoints differ"):
        audit_pair(_chain(0.1), _chain(0.05, final_x=0.21),
                   n_formula=1, pressure_gpa=0.0)
