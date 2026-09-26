"""The GaN Gamma workflow must never create an expanded supercell."""

import numpy as np
import pytest
from ase import Atoms

pytest.importorskip("phonopy")
from scripts.prepare_gan_gamma_phonopy_1x1x1 import make_phonopy


def test_gamma_displacements_keep_four_atom_cell_and_central_pairs():
    atoms = Atoms(
        "Ga2N2",
        cell=[[3.1, 0, 0], [-1.55, 2.68467875, 0], [0, 0, 5.0]],
        scaled_positions=[[0, 0, 0], [2 / 3, 1 / 3, 1 / 2],
                          [0, 0, 3 / 8], [2 / 3, 1 / 3, 7 / 8]],
        pbc=True,
    )
    phonon = make_phonopy(atoms, 0.01)
    assert np.array_equal(phonon.supercell_matrix, np.eye(3, dtype=int))
    assert all(len(cell) == 4 for cell in phonon.supercells_with_displacements)
    displacements = phonon.dataset["first_atoms"]
    assert len(displacements) >= 4 and len(displacements) % 2 == 0
    for first, second in zip(displacements[::2], displacements[1::2]):
        assert first["number"] == second["number"]
        np.testing.assert_allclose(first["displacement"], -np.asarray(second["displacement"]), atol=1e-12)
