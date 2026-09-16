"""Regression checks for the BTO Gamma-mode manuscript figure source data."""

from __future__ import annotations

import json
import runpy
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "plot_bto_gamma_mode_figure.py"
OUTPUTS = ROOT / "outputs" / "batio3_t_to_c_pbe100_dzp10au"


def test_bto_gamma_mode_figure_data_preserves_soft_subspace_and_monotonic_controls() -> None:
    module = runpy.run_path(str(SCRIPT))
    modes = json.loads((OUTPUTS / "bto_tetragonal_to_cubic_n7_gamma_modes.json").read_text(encoding="utf-8"))
    summary = json.loads((OUTPUTS / "bto_tetragonal_to_cubic_n7_static_audit_summary.json").read_text(encoding="utf-8"))
    rows, curves = module["extract_source_data"](modes, summary)
    assert len(rows) == 7
    assert len(curves) == 3
    assert curves[0][0].startswith("soft")
    assert np.all(np.diff([row["relative_enthalpy_meV_per_formula_unit"] for row in rows]) >= 0.0)
    assert np.all(np.diff([row["relative_volume_percent"] for row in rows]) <= 0.0)
