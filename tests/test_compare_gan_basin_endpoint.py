"""GaN endpoint comparison removes only a common periodic translation."""

import numpy as np
import pytest
from ase import Atoms

from scripts.compare_gan_basin_endpoint import relative_displacement_metrics


def test_common_translation_is_not_a_structural_mismatch():
    reference = Atoms("Ga2N2", scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5],
                                                  [0, 0, 0.5], [0.5, 0.5, 0]],
                      cell=np.diag([3.0, 3.0, 4.0]), pbc=True)
    candidate = reference.copy()
    candidate.set_scaled_positions(reference.get_scaled_positions() + [0.27, -0.18, 0.11])
    metrics = relative_displacement_metrics(candidate, reference)
    assert metrics["relative_displacement_rms_A"] == pytest.approx(0, abs=1e-12)
    candidate.positions[2, 2] += 0.02
    metrics = relative_displacement_metrics(candidate, reference)
    assert metrics["relative_displacement_max_A"] == pytest.approx(0.02)
    candidate[0].symbol = "N"
    with pytest.raises(ValueError, match="mapped"):
        relative_displacement_metrics(candidate, reference)
