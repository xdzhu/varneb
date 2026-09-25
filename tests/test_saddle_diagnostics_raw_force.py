"""Do not confuse an ordinary NEB tangent residual with TS stationarity."""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from vcneb import VCNEB


class ShiftedDoubleWell(Calculator):
    implemented_properties = ["energy", "forces", "stress"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x = float(atoms.positions[0, 0]) - 3.0
        self.results = {
            "energy": (x * x - 1.0) ** 2,
            "forces": np.array([[-4.0 * x * (x * x - 1.0), 0.0, 0.0]]),
            "stress": np.zeros((3, 3)),
        }


def test_ordinary_neb_tangent_residual_does_not_certify_a_saddle() -> None:
    images = [
        Atoms("Ar", positions=[[3.0 + x, 3.0, 3.0]],
              cell=np.diag([6.0] * 3), pbc=True)
        for x in (-1.0, -0.6, 0.2, 0.6, 1.0)
    ]
    for image in images:
        image.calc = ShiftedDoubleWell()
    chain = VCNEB(images, k=0.0, climb=False, cell_mask=np.zeros((3, 3)))

    diagnostic = chain.saddle_diagnostics()
    assert diagnostic["image_index"] == 2
    assert diagnostic["has_interior_barrier"]
    assert diagnostic["tangent_curvature_eV_per_A2"] < 0.0
    assert diagnostic["tangential_force_eV_per_A"] == pytest.approx(0.0, abs=1e-12)
    assert abs(diagnostic["true_tangential_force_eV_per_A"]) > 0.5
    assert diagnostic["true_generalized_force_max_vector_eV_per_A"] > 0.5
