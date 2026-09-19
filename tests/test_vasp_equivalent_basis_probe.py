from pathlib import Path
import runpy

import numpy as np
import pytest
from ase import Atoms

MODULE = runpy.run_path(str(Path(__file__).resolve().parents[1] / "scripts" / "prepare_vasp_equivalent_basis_probe.py"))


@pytest.mark.parametrize("transform", MODULE["TRANSFORMS"].values())
def test_probe_preserves_periodic_geometry_and_full_gamma_mesh(transform):
    atoms = Atoms("Ga2N2", cell=[[1.523621257, -2.639367291, 0], [1.523621261, 2.639367294, 0], [0, 0, 4.608269118]],
                  scaled_positions=[[2/3,1/3,0],[1/3,2/3,.5],[2/3,1/3,.488],[1/3,2/3,.988]],pbc=True)
    changed = MODULE["equivalent_atoms"](atoms, transform, [8,8,6])
    assert np.array_equal(changed.positions, atoms.positions)
    assert changed.get_volume() == pytest.approx(atoms.get_volume())
    assert MODULE["mesh_equivalence"](transform,[8,8,6]) == 384


def test_probe_rejects_changed_mesh_or_nonunimodular_cell():
    atoms = Atoms("GaN", cell=[4,4,4], pbc=True)
    with pytest.raises(ValueError,match="equal basal"):
        MODULE["equivalent_atoms"](atoms,MODULE["TRANSFORMS"]["basal_sum_b"],[8,6,6])
    with pytest.raises(ValueError,match="unimodular"):
        MODULE["equivalent_atoms"](atoms,np.diag([2,1,1]),[8,8,6])
