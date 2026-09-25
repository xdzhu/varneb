"""Regression that VARNEB with every cell frozen agrees with ordinary ASE NEB."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def test_fixed_cell_vcneb_matches_ase_improved_tangent(tmp_path: Path) -> None:
    output = tmp_path / "fixed_cell_comparison.json"
    subprocess.run(
        [sys.executable, str(ROOT / "examples" / "run_fixed_cell_ase_comparison.py"),
         "--images", "7", "--fmax", "0.002", "--steps", "300", "--output", str(output)],
        cwd=ROOT, check=True, capture_output=True, text=True,
    )
    report = json.loads(output.read_text(encoding="utf-8"))
    assert report["n_images"] == 7
    assert report["max_cell_deviation_from_reference_A"] < 1e-12
    assert report["barrier_absolute_difference_eV"] < 1e-6
    assert report["max_corresponding_image_position_difference_A"] < 1e-4
    assert report["ase_cineb"]["optimizer_steps"] == report["vcneb_cell_mask_zero"]["optimizer_steps"]
