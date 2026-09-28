"""Reciprocal sampling density is a cell-and-mesh property, not an FFT grid."""

from pathlib import Path

import numpy as np
import pytest

from scripts.audit_gan_k_sampling_geometry import audit, sampling_metrics


def test_equal_density_can_use_different_integer_meshes():
    cells = np.array([np.eye(3) * 4.0, np.eye(3) * 2.0])
    meshes = np.array([[4, 4, 4], [8, 8, 8]])
    report = sampling_metrics(cells, meshes)
    # A smaller direct cell needs a denser mesh for equal reciprocal sampling.
    assert report["bz_volume_per_kpoint_A_minus3"][0] == pytest.approx(
        report["bz_volume_per_kpoint_A_minus3"][1]
    )
    unequal = sampling_metrics(cells, np.array([[4, 4, 4], [4, 4, 4]]))
    assert unequal["bz_volume_per_kpoint_A_minus3"][0] != pytest.approx(
        unequal["bz_volume_per_kpoint_A_minus3"][1]
    )


def test_invalid_mesh_and_cell_are_rejected():
    with pytest.raises(ValueError, match="positive integer"):
        sampling_metrics(np.array([np.eye(3)]), np.array([[4, 0, 4]]))
    with pytest.raises(ValueError, match="positive"):
        sampling_metrics(np.array([-np.eye(3)]), np.array([[4, 4, 4]]))


def test_archived_vasp_chain_matches_hash_bound_sampling_report():
    root = Path(__file__).resolve().parents[1] / "paper" / "VARNEB_CPC" / "evidence"
    report = audit(
        root / "gan_45p7_final_chains_20260927" / "gan_vasp_45p7_final_chain.traj",
        root / "gan_45p7_final_chain_audit_20260927.json",
        (8, 8, 6),
    )
    assert report["n_images_total"] == 29
    assert report["max_adjacent_fractional_change"] == pytest.approx(0.01275501314528609)
    assert report["max_adjacent_image_pair_zero_based"] == [22, 23]
