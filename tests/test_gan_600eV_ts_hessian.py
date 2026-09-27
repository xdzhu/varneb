import numpy as np
import pytest
from ase import Atoms

from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight


def _four_atom_cell(shear_x=0.0):
    cell = np.diag([5.0, 5.0, 5.0])
    cell[2, 0] = shear_x
    return Atoms(
        "Ga2N2",
        positions=[[0, 0, 0], [2.5, 2.5, 0], [2.5, 0, 2.5], [0, 2.5, 2.5]],
        cell=cell,
        pbc=True,
    )


def test_geometry_preflight_records_positive_volume_and_distance():
    result = geometry_preflight(_four_atom_cell())
    assert result["volume_A3"] == pytest.approx(125.0)
    assert result["minimum_distance_A"] > 3.0
    assert result["empirical_near_symmetry_warning"] is False


def test_empirical_warning_is_explicit_not_a_vasp_prediction():
    result = geometry_preflight(_four_atom_cell(2e-4))
    assert result["basal_projection_A"] == pytest.approx(2e-4)
    assert result["empirical_near_symmetry_warning"] is True


def test_overlapping_atoms_rejected():
    atoms = _four_atom_cell()
    atoms.positions[1] = atoms.positions[0] + [0.5, 0.0, 0.0]
    with pytest.raises(ValueError, match="unsafe"):
        geometry_preflight(atoms)
