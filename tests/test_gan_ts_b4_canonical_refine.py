import numpy as np
import pytest
from ase import Atoms

from scripts.prepare_gan_ts_b4_canonical_refine import canonicalize_near_b4


def b4_structure():
    a, c = 3.03, 4.94
    cell = [[a / 2, -np.sqrt(3) * a / 2, 0],
            [a / 2, np.sqrt(3) * a / 2, 0], [0, 0, c]]
    scaled = [[2 / 3, 1 / 3, 0.5], [1 / 3, 2 / 3, 0.0],
              [2 / 3, 1 / 3, 0.875], [1 / 3, 2 / 3, 0.375]]
    return Atoms('Ga2N2', cell=cell, scaled_positions=scaled, pbc=True)


def test_b4_endpoint_canonicalization_is_small_and_phase_preserving():
    reference = b4_structure()
    candidate = reference.copy()
    cell = candidate.cell.array.copy()
    cell[0, 2] += 0.0003
    cell[1, 0] += 0.0002
    candidate.set_cell(cell, scale_atoms=True)
    scaled = candidate.get_scaled_positions()
    scaled += np.array([0.01, -0.01, 0.02])
    scaled[1, 2] += 0.0001
    scaled[2, 0] += 0.0001
    candidate.set_scaled_positions(scaled)
    canonical, metrics = canonicalize_near_b4(candidate, reference)
    assert canonical.get_chemical_symbols() == candidate.get_chemical_symbols()
    assert metrics['canonical_spacegroup_at_1e_minus_4A'] == 'P6_3mc'
    assert metrics['max_atomic_position_change_A'] < 0.005
    assert metrics['max_cell_vector_change_A'] < 0.005


def test_b4_endpoint_canonicalization_rejects_large_perturbation():
    reference = b4_structure()
    candidate = reference.copy()
    cell = candidate.cell.array.copy()
    cell[0, 2] += 0.03
    candidate.set_cell(cell, scale_atoms=True)
    with pytest.raises(ValueError):
        canonicalize_near_b4(candidate, reference)
