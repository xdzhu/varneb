"""Guard the GaN wide-q failure and the fixed-cell atomic alternative."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.analyze_gan_600eV_atomic_tube_feasibility import atomic_normals
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


ROOT = Path(__file__).resolve().parents[1]
WIDE = ROOT / "benchmarks/numerical_integrity/gan_600eV_wide_tube_canary_20260928.json"
ATOMIC = ROOT / "benchmarks/numerical_integrity/gan_600eV_atomic_tube_canary_20260928.json"
TRAJECTORY = ROOT / "paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927/gan_vasp_45p7_final_chain.traj"
HESSIAN = ROOT / "benchmarks/numerical_integrity/gan_600eV_ts_hessian_0p02_20260928/joint_hessian.npz"


def test_wide_joint_chart_records_real_vasp_bravais_failures() -> None:
    report = json.loads(WIDE.read_text(encoding="utf-8"))
    assert report["status"] == "GaN_600eV_wide_frozen_tube_canary_partial_Bravais_failure_raw_audited"
    assert report["n_completed_static_cases"] == 6
    assert report["n_failed_Bravais_cases"] == 2
    assert {item["case"] for item in report["failures"]} == {
        "image_05_q_minus", "image_05_q_plus"
    }
    assert all(item["input_sha256"]["INCAR"] == PRODUCTION_INPUT_SHA256["INCAR"]
               for item in [*report["cases"], *report["failures"]])


def test_atomic_chart_succeeds_but_wide_displacement_is_not_a_final_surface() -> None:
    report = json.loads(ATOMIC.read_text(encoding="utf-8"))
    assert report["status"] == "GaN_600eV_atomic_transverse_canary_all_raw_audited"
    assert report["n_completed_static_cases"] == 4
    assert report["n_failed_cases"] == 0
    image_five = [row for row in report["cases"] if row["image_index"] == 5]
    assert min(row["delta_enthalpy_meV_per_GaN"] for row in image_five) > 200
    assert max(row["maximum_atomic_force_eV_per_A"] for row in image_five) > 3
    frames = read(TRAJECTORY, index=":")
    with np.load(HESSIAN, allow_pickle=False) as archive:
        normals, _ = atomic_normals(frames, archive["eigenvectors"][:, 1])
    assert min(float(normals[i] @ normals[i + 1]) for i in range(5, 22)) > 0.95
    assert min(float(normals[i] @ normals[i + 1]) for i in range(28)) < 0.2
