"""Phase-fixed Phonopy eigenvectors preserve a real orthonormal Gamma basis."""

import numpy as np
import pytest

pytest.importorskip("phonopy")
from scripts.audit_gan_gamma_phonopy_1x1x1 import real_gamma_eigenvectors


def test_phase_fix_recovers_real_columns():
    matrix = np.array([[0, 1j], [-1, 0]], dtype=complex)
    real = real_gamma_eigenvectors(matrix)
    np.testing.assert_allclose(real.T @ real, np.eye(2), atol=1e-12)
    assert np.isrealobj(real)


def test_phase_fix_rejects_intrinsically_complex_modes():
    matrix = np.array([[1, 1j], [1j, 1]], dtype=complex) / np.sqrt(2)
    with pytest.raises(ValueError, match="complex"):
        real_gamma_eigenvectors(matrix)
