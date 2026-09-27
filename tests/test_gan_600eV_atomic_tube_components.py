"""Reproduce GaN tube even/odd failure analysis from audited 600-eV data."""

from __future__ import annotations

import json
from pathlib import Path

from scripts.analyze_gan_600eV_atomic_tube_components import analyze


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "benchmarks/numerical_integrity"


def test_force_information_does_not_validate_contour() -> None:
    calculated = analyze(
        BASE / "gan_600eV_atomic_tube_central_20260928.json",
        BASE / "gan_600eV_atomic_tube_refinement_20260928.json",
        BASE / "gan_600eV_atomic_tube_transport_comparison_20260928.json",
    )
    saved = json.loads((BASE / "gan_600eV_atomic_tube_components_20260928.json").read_text(
        encoding="utf-8"))
    assert calculated == saved
    assert calculated["maximum_even_LOO_absolute_error_meV_per_GaN"] > 1.0
    assert calculated["maximum_odd_LOO_absolute_error_meV_per_GaN"] > 1.0
    assert calculated["maximum_force_informed_signed_LOO_absolute_error_meV_per_GaN"] > 1.0
    assert calculated["contour_certified"] is False
    assert max(calculated["leave_one_anchor_out_components"], key=lambda row: abs(
        row["even_interpolation_error_meV_per_GaN"]))["image_index"] == 8
