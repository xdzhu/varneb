"""The conditional-mode audit must read forces/stress from the ABACUS log."""

import numpy as np
import pytest
from ase import Atoms

from scripts.audit_bto_q1q2_conditional_pilot import _gradient_from_raw, _raw_force_stress
from vcneb.reference_cell import ReferenceCellCoordinates


LOG = """
 TOTAL-FORCE (eV/Angstrom)
---------------------------
 Ba1  0.0000000000  0.0062229451  0.0181086196
 Ti1  0.0000000000  0.2996816268  0.2343792499
 O1   0.0000000000 -0.0998218465 -0.1300343658
 O2   0.0000000000 -0.1066657311 -0.1342345331
 O3   0.0000000000 -0.0994169943  0.0117810294
---------------------------
 TOTAL-STRESS (KBAR)
---------------------------
 5.9689169599 0.0000000000 0.0000000000
 0.0000000000 18.4116300334 2.2974088112
 0.0000000000 2.2974088112 34.6977351888
---------------------------
 !FINAL_ETOT_IS -3737.8763209236112743 eV
"""


def test_parse_raw_abacus_force_stress() -> None:
    forces, stress = _raw_force_stress(LOG, symbols=["Ba", "Ti", "O", "O", "O"])
    assert forces.shape == (5, 3)
    assert stress.shape == (3, 3)
    assert np.isclose(forces[1, 1], 0.2996816268)
    assert np.isclose(stress[2, 2], 34.6977351888)


@pytest.mark.parametrize("bad_log", [LOG.replace(" O3   ", " Ba2  "),
                                           LOG.replace(" 0.0000000000 2.2974088112 34.6977351888", ""),
                                           LOG + LOG])
def test_reject_raw_abacus_force_stress_mismatch(bad_log: str) -> None:
    with pytest.raises(ValueError):
        _raw_force_stress(bad_log, symbols=["Ba", "Ti", "O", "O", "O"])


def test_raw_gradient_includes_deformation_and_shear_factor() -> None:
    chart = ReferenceCellCoordinates(
        Atoms("H", positions=[[0.0, 0.0, 0.0]], cell=np.diag([2.0, 3.0, 4.0]), pbc=True),
    )
    coordinates = np.zeros(9)
    coordinates[-6] = 0.1  # xx strain
    forces = np.array([[1.0, 2.0, 3.0]])
    stress = np.zeros((3, 3))
    stress[1, 2] = stress[2, 1] = 0.05
    gradient = _gradient_from_raw(chart, coordinates, forces, stress)
    assert np.allclose(gradient[:3], [-1.1, -2.0, -3.0])
    assert np.allclose(gradient[3:6], 0.0)
    assert gradient[6] == pytest.approx(2 * (2.0 * 3.0 * 4.0 * 1.1) * 0.05)
    assert np.allclose(gradient[7:], 0.0)
