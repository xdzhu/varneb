import pytest

from scripts.prepare_gan_ts_basin_pilot import relaxation_incar


STATIC = """INCAR created by Atomic Simulation Environment
 ENCUT = 1000.000000
 SIGMA = 0.050000
 EDIFF = 1.00e-07
 SYMPREC = 1.00e-04
 GGA = PE
 PREC = Accurate
 IBRION = -1
 ISIF = 2
 ISMEAR = 0
 ISYM = -1
 NSW = 0
 LCHARG = .FALSE.
 LWAVE = .FALSE.
"""


def test_pressure_relax_changes_only_declared_ionic_tags():
    result = relaxation_incar(STATIC)
    assert " ENCUT = 1000.000000\n" in result
    assert " SYMPREC = 1.00e-04\n" in result
    assert " IBRION = 2\n" in result
    assert " ISIF = 3\n" in result
    assert " NSW = 10\n" in result
    assert " PSTRESS = 457.0\n" in result
    assert " EDIFFG = -0.02\n" in result
    assert " POTIM = 0.25\n" in result
    assert "ISYM = -1" in result


def test_pressure_relax_rejects_changed_calculator_contract():
    with pytest.raises(ValueError, match="contract changed"):
        relaxation_incar(STATIC.replace("1000.000000", "600.000000"))
    with pytest.raises(ValueError, match="unexpectedly"):
        relaxation_incar(STATIC + " PSTRESS = 457.0\n")
