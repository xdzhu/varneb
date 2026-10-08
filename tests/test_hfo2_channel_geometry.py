"""Small synthetic ordered paths test geometry bookkeeping, not HfO2 physics."""

import numpy as np
import pytest
from ase import Atoms
from ase.io import write

from scripts.compare_hfo2_channel_geometries import chain_record


def write_path(directory, *, excursion=0.0):
    directory.mkdir()
    atoms = Atoms("Hf4O8", cell=np.diag([5.1, 5.2, 5.3]), pbc=True)
    initial = np.linspace(0.05, 0.9, 36).reshape(12, 3)
    for image in range(7):
        t = image / 6
        q = initial.copy()
        q[4, 0] += 0.08 * t + excursion * np.sin(np.pi * t)
        atoms.set_scaled_positions(q)
        write(directory / f"POSCAR_{image:02d}", atoms, format="vasp", direct=True, vasp5=True)


def test_same_endpoints_different_intermediate_geometry(tmp_path):
    ordinary, guided = tmp_path / "ordinary", tmp_path / "guided"
    write_path(ordinary)
    write_path(guided, excursion=0.02)
    left, right = chain_record(ordinary), chain_record(guided)
    assert [a["sha256"] == b["sha256"] for a, b in zip(left["endpoints"], right["endpoints"])] == [True, True]
    assert left["final_unwrap_integers"] == right["final_unwrap_integers"]
    assert left["maximum_segment_A"] > 0
    assert not np.allclose(left["translation_free_atomic_displacements_reference_A"],
                           right["translation_free_atomic_displacements_reference_A"])


def test_incomplete_chain_is_not_interpreted(tmp_path):
    chain = tmp_path / "chain"
    write_path(chain)
    (chain / "POSCAR_03").unlink()
    with pytest.raises(ValueError, match="complete"):
        chain_record(chain)
