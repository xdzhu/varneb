"""Unit and sign guards for the two-basin GaN enthalpy summary."""

import pytest

from scripts.summarize_gan_ts_bilateral_evidence import (
    enthalpy_barriers_per_formula_unit,
)


def test_barriers_are_per_formula_unit_and_obey_cycle_identity():
    forward, reverse, relative = enthalpy_barriers_per_formula_unit(
        center_eV_per_cell=-11.5,
        initial_eV_per_cell=-12.2,
        final_eV_per_cell=-12.3,
        formula_units_per_cell=2,
    )
    assert forward == pytest.approx(0.35)
    assert reverse == pytest.approx(0.4)
    assert relative == pytest.approx(-0.05)
    assert forward - reverse == pytest.approx(relative)


@pytest.mark.parametrize(
    'center,initial,final,count',
    [(-12.3, -12.2, -12.4, 2), (-11.5, -12.2, -12.3, 0),
     (float('nan'), -12.2, -12.3, 2)],
)
def test_invalid_or_nonstationary_inputs_rejected(center, initial, final, count):
    with pytest.raises(ValueError):
        enthalpy_barriers_per_formula_unit(center, initial, final, count)
