"""Regression checks for the audited GaN split-path figure semantics."""

from pathlib import Path

import numpy as np
from ase.io.trajectory import Trajectory

from scripts.plot_gan_45p7_two_segment_on_local_cut import (
    clipped_ray,
    evidence_rows,
    project,
)


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "paper/VARNEB_CPC/evidence/gan_45p7_split_20261004"
SURFACE = ROOT / "paper/VARNEB_CPC/figures/gan_600eV_local_joint_dft289_20260930_v2_source_data.csv"


def test_clipped_ray_stays_inside_local_cut():
    assert np.allclose(clipped_ray((0.072, 0.0094), (0.02, 0.0125)),
                       (0.02, 0.0094 * 0.02 / 0.072))
    assert np.allclose(clipped_ray((-0.102, 0.0096), (0.02, 0.0125)),
                       (-0.02, 0.0096 * 0.02 / 0.102))


def test_center_projects_to_zero():
    with Trajectory(str(EVIDENCE / "left.traj"), "r") as trajectory:
        center = trajectory[len(trajectory) - 1]
    with np.load(EVIDENCE / "joint_hessian.npz", allow_pickle=False) as archive:
        modes = archive["eigenvectors"]
    assert np.allclose(project(center, center, modes), (0, 0, 0), atol=1e-12)


def test_audited_path_has_only_center_in_local_chart():
    path, surface, audit = evidence_rows(
        EVIDENCE / "raw_audit.json", EVIDENCE,
        EVIDENCE / "joint_hessian.npz", SURFACE,
    )
    assert len(surface) == 289
    assert len(path) == 29
    assert [r["stitched_index"] for r in path if r["inside_measured_local_rectangle"]] == [15]
    assert path[15]["enthalpy_meV_per_GaN_relative_B4"] == audit["forward_barrier_meV_per_GaN"]
