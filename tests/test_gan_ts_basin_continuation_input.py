"""The GaN basin continuation preserves the VASP physical contract."""

import numpy as np
import pytest
from ase import Atoms

from scripts.prepare_gan_ts_basin_continuation import (
    continuation_incar,
    geometry_metrics,
)


def test_continuation_incar_changes_only_bounded_nsw():
    original = (
        b" ENCUT = 1000.000000\n EDIFF = 1.00e-07\n ISYM = -1\n"
        b" SYMPREC = 1.00e-04\n IBRION = 2\n ISIF = 3\n NSW = 10\n"
        b" PSTRESS = 457.0\n EDIFFG = -0.02\n POTIM = 0.25\n"
    )
    expected = original.replace(b" NSW = 10\n", b" NSW = 40\n")
    assert continuation_incar(original) == expected
    with pytest.raises(ValueError, match="LF-only"):
        continuation_incar(original.replace(b"\n", b"\r\n"))
    with pytest.raises(ValueError, match="setting changed"):
        continuation_incar(original.replace(b"ENCUT = 1000.000000", b"ENCUT = 600"))


def test_geometry_preflight_rejects_collisions_and_bad_cell():
    atoms = Atoms("Ga2N2", positions=[[0, 0, 0], [2, 2, 0], [1, 1, 1.5], [3, 3, 1.5]],
                  cell=np.diag([4.0, 4.0, 3.0]), pbc=True)
    metrics = geometry_metrics(atoms)
    assert metrics["volume_A3"] == pytest.approx(48.0)
    atoms.positions[2] = [0.1, 0.1, 0.1]
    with pytest.raises(ValueError, match="unsafe"):
        geometry_metrics(atoms)
