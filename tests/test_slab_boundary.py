"""Two-dimensional symmetric cell freedom before any hBN production run."""

import numpy as np
import pytest
from ase import Atoms

from vcneb import VCNEB, symmetric_inplane_vcneb_boundary
from vcneb.core import cell_force_from_arrays, cell_from_deformation, deformation_from_cell


def _cell():
    return np.diag([6.0, 6.0, 20.0])


def test_symmetric_inplane_boundary_has_no_rotation_or_vacuum_direction():
    reference = _cell()
    mask, basis = symmetric_inplane_vcneb_boundary(1, reference)
    assert mask.shape == (3, 3)
    assert np.array_equal(mask, np.array([[1, 1, 0], [1, 1, 0], [0, 0, 0]]))
    assert basis.shape == (12, 6)
    assert np.allclose(basis.T @ basis, np.eye(6))
    images = [
        Atoms("Ar", positions=[[x, 3.0, 10.0]], cell=reference, pbc=True)
        for x in (2.0, 3.0, 4.0)
    ]
    neb = VCNEB(
        images, cell_scale=1.0, cell_mask=mask, mode_basis=basis,
        constraint_mode="subspace",
    )
    proposed = neb.get_x()
    proposed[4] = 0.6  # F_xy
    proposed[6] = -0.4  # F_yx: antisymmetric component should be removed
    proposed[11] = 0.5  # F_zz: vacuum strain should be removed
    candidate = neb.candidate_images_from_x(proposed)[1]
    deformation = deformation_from_cell(candidate.cell.array, reference)
    assert deformation[0, 1] == pytest.approx(0.1)
    assert deformation[1, 0] == pytest.approx(0.1)
    assert np.allclose(deformation[2], [0.0, 0.0, 1.0])
    assert np.allclose(deformation[:, 2], [0.0, 0.0, 1.0])
    assert candidate.cell[2, 2] == pytest.approx(20.0)


def test_symmetric_shear_stress_matches_energy_finite_difference():
    reference = _cell()
    mask, basis = symmetric_inplane_vcneb_boundary(1, reference)
    shear_direction = basis[3:, -1].reshape(3, 3)
    cell_scale = 6.0
    stiffness, preferred_shear = 8.0, 0.12

    def model(t):
        deformation = np.eye(3) + (t / cell_scale) * shear_direction
        atoms = Atoms(
            "Ar", positions=[[3.0, 3.0, 10.0]],
            cell=cell_from_deformation(deformation, reference), pbc=True,
        )
        shear = 0.5 * (deformation[0, 1] + deformation[1, 0])
        energy = 0.5 * stiffness * (shear - preferred_shear) ** 2
        first_piola = np.zeros((3, 3))
        first_piola[0, 1] = first_piola[1, 0] = 0.5 * stiffness * (shear - preferred_shear)
        stress = first_piola @ deformation.T / atoms.get_volume()
        assert np.allclose(stress, stress.T, atol=1e-12)
        force = cell_force_from_arrays(stress, atoms, reference, mask=mask)
        return energy, -float(np.sum(force * shear_direction)) / cell_scale

    for t in (-0.2, 0.0, 0.3):
        h = 1e-4
        energy_plus, _ = model(t + h)
        energy_minus, _ = model(t - h)
        _, stress_derivative = model(t)
        assert stress_derivative == pytest.approx(
            (energy_plus - energy_minus) / (2 * h), abs=1e-10,
        )


def test_boundary_rejects_tilted_vacuum_or_bad_atom_count():
    tilted = _cell()
    tilted[2, 0] = 0.5
    with pytest.raises(ValueError, match="z-directed vacuum"):
        symmetric_inplane_vcneb_boundary(4, tilted)
    with pytest.raises(ValueError, match="positive integer"):
        symmetric_inplane_vcneb_boundary(0, _cell())
