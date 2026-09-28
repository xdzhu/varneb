"""The dense GaN array must add only missing same-protocol grid points."""

from __future__ import annotations

import json
from pathlib import Path

from scripts.audit_gan_600eV_atomic_tube_dense import validation_errors
from scripts.prepare_gan_600eV_atomic_tube_dense import (
    CENTRAL_IMAGES,
    PROSPECTIVE_HOLDOUT_IMAGES,
    TRANSVERSE_Q_A,
    missing_points,
)


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / "benchmarks/numerical_integrity"


def test_dense_grid_reuses_all_prior_statics() -> None:
    first = json.loads(
        (BENCHMARKS / "gan_600eV_atomic_tube_central_20260928.json").read_text()
    )
    refinement = json.loads(
        (BENCHMARKS / "gan_600eV_atomic_tube_refinement_20260928.json").read_text()
    )
    new = missing_points(first, refinement)
    existing = {
        (int(case["image_index"]), float(case["q_atom_A"]))
        for case in first["cases"] + refinement["new_cases"]
    }
    target = {
        (index, q)
        for index in CENTRAL_IMAGES
        for q in TRANSVERSE_Q_A
        if q != 0.0
    }
    assert len(CENTRAL_IMAGES) == 18
    assert len(target) == 72
    assert len(existing) == 28
    assert len(new) == 44
    assert set(new).isdisjoint(existing)
    assert set(new) | existing == target
    assert set(PROSPECTIVE_HOLDOUT_IMAGES) == {10, 16}


def test_dense_grid_prospective_holdouts_detect_unmodeled_ridge() -> None:
    arc = [index / 28 for index in range(29)]
    smooth = {
        index: {q: 4 * q + 2 * q**2 + 0.5 * arc[index] * q
                for q in TRANSVERSE_Q_A}
        for index in CENTRAL_IMAGES
    }
    s_errors, q_errors = validation_errors(arc, smooth)
    assert max(item["absolute_error_meV_per_GaN"] for item in s_errors) < 1e-12
    assert max(item["absolute_error_meV_per_GaN"] for item in q_errors) < 1e-12
    smooth[10][0.05] += 2.0
    s_errors, q_errors = validation_errors(arc, smooth)
    assert max(item["absolute_error_meV_per_GaN"] for item in s_errors) > 1.99
    assert max(item["absolute_error_meV_per_GaN"] for item in q_errors) >= 0.5
