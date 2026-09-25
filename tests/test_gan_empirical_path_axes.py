"""Regression contract for the published, non-phonon GaN path coordinates."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "benchmarks/numerical_integrity/gan_b4_b1_tetragonal_empirical_atomic_strain_20260926.csv"
AUDIT = ROOT / "outputs/neb_literature_benchmarks/gan_tetragonal_empirical_atomic_axes_audit_20260926.json"


def test_published_gan_atomic_axes_are_traceable_and_not_called_phonons() -> None:
    report = json.loads(AUDIT.read_text(encoding="utf-8"))
    with CSV.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert report["status"] == "audited_GaN_path_geometric_empirical_axes_not_phonons"
    assert report["atomic_axis_kind"] == "mass_weighted_path_SVD_empirical_not_phonon"
    assert report["endpoint_identity"]["status"] == "passed"
    assert report["phonon_force_constants_available_in_this_analysis"] is False
    assert report["full_variable_cell_TS_certified"] is False
    assert len(rows) == report["n_images_total"] == len(report["rows"]) == 29
    assert report["n_formula_units"] == 2
    assert set(report["source_sha256"]) == {
        "summary", "trajectory", "generic_audit", "route_analysis", "ts_audit", "auditor"
    }
    assert all(len(value) == 64 for value in report["source_sha256"].values())

    s = np.array([float(row["reaction_coordinate_normalized"]) for row in rows])
    h = np.array([float(row["relative_enthalpy_eV_per_GaN"]) for row in rows])
    residual = np.array([float(row["rank3_held_out_residual_sqrt_amu_A"]) for row in rows])
    assert s[0] == pytest.approx(0.0)
    assert s[-1] == pytest.approx(1.0)
    assert np.all(np.diff(s) > 0.0)
    assert int(np.argmax(h)) == report["peak_image_index"] == 15
    assert h.max() == pytest.approx(report["peak_enthalpy_eV_per_GaN"], abs=1e-12)
    assert np.max(residual) == pytest.approx(report["max_rank3_held_out_residual_sqrt_amu_A"])
    assert report["captured_squared_atomic_path_norm_fraction_by_rank"]["3"] > 0.99999999
    assert report["max_cell_deformation_antisymmetry"] < 1e-5

    for parsed, stored in zip(rows, report["rows"]):
        assert int(parsed["image_index"]) == stored["image_index"]
        for key in parsed:
            if key != "image_index":
                assert float(parsed[key]) == pytest.approx(stored[key], abs=1e-13)
