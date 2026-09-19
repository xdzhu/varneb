"""Calculator-free route/mapping tests; no licensed PAW data required."""
from pathlib import Path
import runpy

import numpy as np
import pytest
from ase import Atoms
from ase.io import write

ROOT = Path(__file__).resolve().parents[1]
SUITE = runpy.run_path(str(ROOT / "scripts" / "prepare_gan_qian_suite.py"))
DRIVER = runpy.run_path(str(ROOT / "examples" / "run_vcneb_vasp.py"))


def test_conventional_b3_b1_have_same_count_and_expected_symmetry():
    for offset, number in ((.25, 216), (.5, 225)):
        atoms = SUITE["conventional_gan"](4.2, offset)
        assert atoms.get_chemical_symbols() == ["Ga"] * 4 + ["N"] * 4
        assert all(record["number"] == number for record in SUITE["symmetry_audit"](atoms))


def test_b3_diagonal_path_preserves_all_three_components_and_winding():
    initial = SUITE["conventional_gan"](4.3, .25)
    final = SUITE["conventional_gan"](4.05, .5)
    chain = SUITE["body_diagonal_chain"](initial, final)
    displacement = chain[-1].get_scaled_positions(wrap=False) - chain[0].get_scaled_positions(wrap=False)
    assert len(chain) == 29
    assert np.allclose(displacement[:4], 0)
    assert np.allclose(displacement[4:], .25)
    assert np.max(chain[-1].get_scaled_positions(wrap=False)) == pytest.approx(1)


def test_hexagonal_anchor_has_h_mgo_symmetry_and_half_period():
    # Audited 45.7 GPa reference geometry; never run a DFT calculation here.
    initial = Atoms("Ga2N2", pbc=True,
                    cell=[[1.514987907424, -2.624036028511, 0],
                          [1.514987907424, 2.624036028511, 0], [0, 0, 4.93588733854]],
                    scaled_positions=[[2/3, 1/3, .499352376464], [1/3, 2/3, .999352376464],
                                      [2/3, 1/3, .875647623561], [1/3, 2/3, .375647623561]])
    final = Atoms("Ga2N2", pbc=True,
                  cell=[[2.029111234228, -2.029111234228, 0],
                        [2.029111234228, 2.029111234228, 0], [0, 0, 4.058064931312]],
                  scaled_positions=[[2/3, 1/3, .499352376464], [1/6, 5/6, .999352376464],
                                    [2/3, 1/3, .999352376464], [1/6, 5/6, .499352376464]])
    chain = SUITE["hexagonal_chain"](initial, final)
    anchor = chain[11]
    assert SUITE["symmetry_audit"](anchor)[1]["number"] == 194
    q = anchor.get_scaled_positions(wrap=False)
    assert np.abs(q[2:, 2] - q[:2, 2]) == pytest.approx([.5, .5])
    assert anchor.cell.angles()[2] == pytest.approx(initial.cell.angles()[2])
    assert np.allclose(chain[-1].cell, final.cell)


def test_explicit_chain_preserves_unwrapped_motion_and_rejects_extra_frames(tmp_path):
    initial = Atoms("Ba", cell=[4, 4, 4], pbc=True)
    final = initial.copy()
    final.set_scaled_positions([[1, 0, 0]])
    middle = initial.copy()
    middle.set_scaled_positions([[.6, 0, 0]])
    path = tmp_path / "seed.traj"
    write(path, [initial, middle, final])
    chain = DRIVER["read_explicit_initial_chain"](path, initial, initial, 3)
    assert chain[-1].get_scaled_positions(wrap=False)[0, 0] == pytest.approx(1)
    with pytest.raises(ValueError, match="exactly one"):
        DRIVER["read_explicit_initial_chain"](path, initial, initial, 2)


def test_explicit_chain_rejects_different_endpoint_mapping(tmp_path):
    initial = Atoms("Ba", cell=[4, 4, 4], pbc=True)
    wrong = initial.copy()
    wrong.set_scaled_positions([[.2, 0, 0]])
    path = tmp_path / "wrong.traj"
    write(path, [initial, initial, wrong])
    with pytest.raises(ValueError, match="endpoint atom mapping"):
        DRIVER["read_explicit_initial_chain"](path, initial, initial, 3)


def test_explicit_chain_rejects_final_species_identity_change(tmp_path):
    initial = Atoms("Ba", cell=[4, 4, 4], pbc=True)
    final = Atoms("Sr", cell=[4, 4, 4], pbc=True)
    path = tmp_path / "species.traj"
    write(path, [initial, initial, initial])
    with pytest.raises(ValueError, match="endpoint order/species"):
        DRIVER["read_explicit_initial_chain"](path, initial, final, 3)
