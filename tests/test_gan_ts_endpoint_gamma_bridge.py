"""The joint TS direction has a mass-weighted, gauge-clean Gamma projection."""

from pathlib import Path

import numpy as np
import pytest

from scripts.audit_gan_ts_endpoint_gamma_bridge import (
    analyze, block_spectrum, project_atomic_direction,
)


def _two_atom_basis(masses):
    first, second = np.sqrt(masses)
    total = np.sqrt(np.sum(masses))
    basis = np.zeros((6, 6))
    for direction in range(3):
        basis[direction, direction] = first / total
        basis[3 + direction, direction] = second / total
        basis[direction, 3 + direction] = second / total
        basis[3 + direction, 3 + direction] = -first / total
    return basis


def test_mass_center_removal_cleans_spurious_acoustic_overlap():
    masses = np.array([4.0, 1.0])
    basis = _two_atom_basis(masses)
    direction = np.array([[1.0, 0, 0], [-1.0, 0, 0]])
    projected = project_atomic_direction(direction, masses, basis, [[3], [4], [5]])
    assert projected["raw_acoustic_squared_fraction_before_com_removal"] > 0.1
    assert projected["acoustic_squared_fraction_after_com_removal"] < 1e-25
    np.testing.assert_allclose(
        projected["group_fractions_of_mass_weighted_atomic_optical_direction"],
        [1.0, 0.0, 0.0], atol=1e-14,
    )
    shifted = project_atomic_direction(
        direction + [0.7, -0.2, 0.4], masses, basis, [[3], [4], [5]],
    )
    np.testing.assert_allclose(
        shifted["group_fractions_of_mass_weighted_atomic_optical_direction"],
        projected["group_fractions_of_mass_weighted_atomic_optical_direction"],
        atol=1e-14,
    )


def test_degenerate_subspace_fraction_is_basis_rotation_invariant():
    masses = np.array([4.0, 1.0])
    basis = _two_atom_basis(masses)
    direction = np.array([[1.0, 2.0, 3.0], [-1.0, -2.0, -3.0]])
    original = project_atomic_direction(direction, masses, basis, [[3, 4], [5]])
    rotated = basis.copy()
    angle = 0.77
    rotation = np.array([[np.cos(angle), -np.sin(angle)],
                         [np.sin(angle), np.cos(angle)]])
    rotated[:, 3:5] = basis[:, 3:5] @ rotation
    transformed = project_atomic_direction(direction, masses, rotated, [[3, 4], [5]])
    np.testing.assert_allclose(
        original["group_fractions_of_mass_weighted_atomic_optical_direction"],
        transformed["group_fractions_of_mass_weighted_atomic_optical_direction"],
        atol=1e-14,
    )


def test_missing_atomic_direction_is_not_misreported_as_phonon():
    masses = np.array([4.0, 1.0])
    with pytest.raises(ValueError, match="only a common translation"):
        project_atomic_direction(
            np.ones((2, 3)), masses, _two_atom_basis(masses), [[3], [4], [5]],
        )


def test_positive_frozen_blocks_can_have_joint_coupling_instability():
    spectrum = block_spectrum(np.array([[2.0]]), np.array([[3.0]]), np.array([[3.0]]))
    assert spectrum["lowest_frozen_atomic_eV_per_A2"] == 2.0
    assert spectrum["lowest_frozen_strain_eV_per_A2"] == 3.0
    assert spectrum["lowest_joint_two_eV_per_A2"][0] < 0
    assert spectrum["lowest_strain_relaxed_atomic_schur_eV_per_A2"] == -1.0


def test_archived_gan_joint_gamma_bridge_recomputes_from_eigenvectors():
    root = Path(__file__).resolve().parents[1] / "outputs"
    result = analyze(
        root / "gan_b4_b1_gamma_1x1x1_20260926",
        root / "gan_b4_b1_joint_curvature_20260926",
    )
    assert result["joint_unstable_direction_cross_step_absolute_overlap"] > 0.99999
    assert result["joint_unstable_direction_strain_coordinate_squared_fraction"] > 0.6
    for blocks in result["joint_hessian_atomic_strain_block_diagnostics"].values():
        assert blocks["lowest_frozen_atomic_eV_per_A2"] > 2.4
        assert blocks["lowest_frozen_strain_eV_per_A2"] > 3.9
        assert blocks["lowest_joint_two_eV_per_A2"][0] < -4.0
        assert blocks["lowest_strain_relaxed_atomic_schur_eV_per_A2"] < -16.0
    for phase in ("B4", "B1"):
        for phonon_step in ("d0.01", "d0.02"):
            for hessian_step in ("0p01", "0p02"):
                projected = result["phases"][phase][phonon_step][hessian_step]
                assert projected["raw_acoustic_squared_fraction_before_com_removal"] > 0.03
                assert projected["acoustic_squared_fraction_after_com_removal"] < 1e-6
                assert projected["path_selected_three_group_atomic_optical_fraction"] > 0.99999
