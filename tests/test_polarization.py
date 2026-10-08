import numpy as np
import pytest

from vcneb.polarization import quantum_lattice, modular_difference, parse_abacus_berry, sampled_band_gap


def test_quantum_lattice_units_and_triclinic_rows():
    np.testing.assert_allclose(quantum_lattice(np.eye(3) * 5), np.eye(3) * (16.02176634 / 25))
    cell = np.array([[4., 1., 0.], [0., 5., 0.], [.3, .2, 6.]])
    np.testing.assert_allclose(quantum_lattice(cell), cell * 16.02176634 / 120)
    for invalid in (np.zeros((3, 3)), np.diag([-1., 1., 1.]), np.ones((2, 3)), np.eye(3)*np.nan):
        with pytest.raises(ValueError):
            quantum_lattice(invalid)


def test_native_scalar_modulus_not_automatically_halved():
    body = """The calculated polarization direction is in R3 direction
    P = 0.5 (mod 1.2) (0.0, 0.0, 0.5) C/m^2
    """
    result = parse_abacus_berry(body)
    assert result["reported_modulus_C_m2"] == 1.2
    assert result["value_C_m2"] == .5 and not result["branch_selected"]
    assert modular_difference(.7, -.5, 1.2) == pytest.approx(0)
    assert modular_difference(.7, .5, 1.2) == pytest.approx(.2)
    for malformed in (body + body, body.replace("C/m^2", "e/bohr^2"), body.replace("mod 1.2", "mod -1.2"),
                      body.replace("0.0, 0.0, 0.5", "0.0, 0.0, 0.6"), body.replace("R3", "R4")):
        with pytest.raises(ValueError):
            parse_abacus_berry(malformed)
    for period in (0, -1, np.nan):
        with pytest.raises(ValueError):
            modular_difference(0, 1, period)


def test_sampled_indirect_gap_handles_k_weighted_occupations():
    body = """BAND Energy(ev) Occupation Kpoint = 1 (0 0 0)
1 -2 .25
2 -1 .25
3 3 0
BAND Energy(ev) Occupation Kpoint = 2 (0.5 0 0)
1 -2 .25
2 -.5 .25
3 2 0
"""
    result = sampled_band_gap(body, 2)
    assert result["sampled_indirect_gap_eV"] == 2.5 and result["n_kpoints"] == 2
    for invalid in (body.replace("3 2 0", "4 2 0"), body.replace("3 2 0", "3 nan 0"),
                    body.replace("3 2 0", "3 -.6 0"), ""):
        with pytest.raises(ValueError):
            sampled_band_gap(invalid, 2)
    with pytest.raises(ValueError):
        sampled_band_gap(body, 3)
