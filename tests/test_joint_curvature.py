"""The variable-cell saddle coordinates must be energy/force consistent."""

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT

from vcneb.joint_curvature import (
    JointCurvatureCoordinates, central_difference_hessian,
    symmetric_strain_basis,
)


@pytest.mark.parametrize("deformed", [False, True])
def test_joint_gradient_matches_enthalpy_finite_difference(deformed: bool) -> None:
    reference = bulk("Cu", "fcc", a=3.6, cubic=True)
    coordinates = JointCurvatureCoordinates(reference, cell_scale_A=3.6)
    pressure = 0.03
    center = np.zeros(coordinates.size)
    if deformed:
        center[0] = 0.013
        center[3 * len(reference) + 0] = 0.025
        center[3 * len(reference) + 4] = -0.011
    at_center = coordinates.displaced(center)
    at_center.calc = EMT()
    gradient = coordinates.enthalpy_gradient(
        at_center, at_center.get_forces(), at_center.get_stress(voigt=False), pressure,
    )
    numerical = []
    for index in range(coordinates.size):
        pair = []
        for sign in (+1, -1):
            displacement = center.copy()
            displacement[index] += sign * 1e-5
            atoms = coordinates.displaced(displacement)
            assert len(atoms) == len(reference)
            atoms.calc = EMT()
            pair.append(atoms.get_potential_energy() + pressure * atoms.get_volume())
        numerical.append((pair[0] - pair[1]) / (2e-5))
    np.testing.assert_allclose(gradient, numerical, rtol=0, atol=3e-5)
    basis = coordinates.translation_free_basis()
    assert basis.shape == (coordinates.size, coordinates.size - 3)
    np.testing.assert_allclose(basis.T @ basis, np.eye(coordinates.size - 3), atol=1e-12)


def test_strain_basis_and_hessian_reciprocity() -> None:
    basis = symmetric_strain_basis()
    np.testing.assert_allclose(np.einsum("aij,bij->ab", basis, basis), np.eye(6), atol=1e-12)
    true_hessian = np.diag([-2.0, 1.0, 3.0, 4.0])
    true_hessian[0, 1] = true_hessian[1, 0] = 0.25
    h = 0.02
    plus = np.array([true_hessian @ (h * np.eye(4)[i]) for i in range(4)])
    minus = -plus
    recovered, defect = central_difference_hessian(plus, minus, h)
    np.testing.assert_allclose(recovered, true_hessian, atol=1e-12)
    assert defect == pytest.approx(0, abs=1e-12)
    plus[1, 0] += 0.001
    _, defect = central_difference_hessian(plus, minus, h)
    assert defect > 0


def test_joint_coordinate_rejects_invalid_geometry() -> None:
    reference = bulk("Cu", "fcc", a=3.6, cubic=True)
    coordinates = JointCurvatureCoordinates(reference, cell_scale_A=3.6)
    with pytest.raises(ValueError, match="invalid shape"):
        coordinates.displaced(np.zeros(coordinates.size - 1))
    with pytest.raises(ValueError, match="non-positive"):
        delta = np.zeros(coordinates.size)
        delta[3 * len(reference)] = -4.0
        coordinates.displaced(delta)
