"""Pressure residuals are not raw VASP compressive stresses."""

import numpy as np
import pytest
from ase.units import GPa

from scripts.audit_gan_ts_basin_followups import residual_stress_kbar


def test_pressure_target_is_minus_45p7_GPa_in_ase_stress_sign():
    target = -45.7 * GPa * np.eye(3)
    assert residual_stress_kbar(target, 45.7) == pytest.approx(0.0)
    shifted = target.copy()
    shifted[0, 0] += 0.1 * GPa
    shifted[0, 1] = shifted[1, 0] = 0.03 * GPa
    assert residual_stress_kbar(shifted, 45.7) == pytest.approx(1.0)
    with pytest.raises(ValueError, match="invalid"):
        residual_stress_kbar(np.zeros(6), 45.7)
