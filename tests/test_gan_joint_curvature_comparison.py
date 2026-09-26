"""Two-step Hessian comparisons must retain the same underlying GaN source."""

import json

import numpy as np
import pytest

from scripts.compare_gan_joint_curvature_steps import compare


def _make_step(directory, step: float, lowest: float, *, source: str = "a") -> None:
    directory.mkdir(exist_ok=True)
    eigenvalues = np.arange(1.0, 16.0)
    eigenvalues[0] = lowest
    eigenvectors = np.eye(18)[:, :15]
    np.savez_compressed(
        directory / "joint_hessian.npz", hessian=np.eye(18),
        eigenvalues=eigenvalues, eigenvectors=eigenvectors,
    )
    report = {
        "status": "GaN_image15_local_joint_curvature_at_nonstationary_NEB_candidate",
        "n_static_displacements": 36,
        "step_A": step,
        "raw_hessian_reciprocity_relative_defect": 0.001,
        "translation_null_relative_defect": 0.001,
        "energy_gradient_max_abs_difference_eV_per_A": 0.01,
        "energy_hessian_diagonal_max_abs_difference_eV_per_A2": 0.1,
        "tangent_squared_overlap_with_eigenvectors": [0.95] + [0.0] * 14,
        "center_gradient_translation_free_eV_per_A": 0.07,
        "source_sha256": {key: source for key in (
            "trajectory", "summary", "joint_curvature", "center_OUTCAR",
        )},
    }
    (directory / "audit.json").write_text(json.dumps(report), encoding="utf-8")


def test_compares_negative_direction_without_claiming_ts(tmp_path) -> None:
    small, large = tmp_path / "small", tmp_path / "large"
    _make_step(small, 0.01, -4.0)
    _make_step(large, 0.02, -4.1)
    report = compare(small, large)
    assert report["negative_count_below_minus_0p1_eV_per_A2"] == [1, 1]
    assert report["lowest_mode_squared_overlap_between_steps"] == 1.0
    assert "not_TS_certificate" in report["status"]
    _make_step(large, 0.02, -4.1, source="changed")
    with pytest.raises(ValueError, match="different trajectory"):
        compare(small, large)
