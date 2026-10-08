from pathlib import Path

import numpy as np
import pytest

from scripts.analyze_hfo2_atomic_reduction import analyze


def test_real_matrix_reduction_io_preserves_coordinate_definition(tmp_path):
    # Software/I/O regression on archived DFT matrices, NOT independent DFT proof.
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"
    output = tmp_path / "analysis.json"
    result = analyze(root / "reference_variants", root / "gamma_analysis", output)
    q = np.array(result["retained_basis_cartesian"])
    np.testing.assert_allclose(q.T @ q, np.eye(3), atol=1e-12)
    assert result["eliminated_dimensions"] == 30 and result["clamped_cell"]
    for family in result["families"].values():
        response = np.array(family["full_atomic_response"])
        np.testing.assert_allclose(q.T @ response, 0, atol=1e-12)
        np.testing.assert_allclose(response.reshape(12, 3, 3).sum(axis=0), 0, atol=1e-12)
        assert family["minimum_eliminated_curvature_eV_A2"] > result["observed_two_step_operator_spread_eV_A2"]
        assert .58 < family["x_softening_fraction"] < .60
    with pytest.raises(FileExistsError):
        analyze(root / "reference_variants", root / "gamma_analysis", output)
