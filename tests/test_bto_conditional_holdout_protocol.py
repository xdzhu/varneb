"""Keep the preregistered BTO center prediction tied to four measured corners."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "paper" / "VARNEB_CPC" / "figures"
SOURCE = FIGURES / "bto_conditional_four_point_stage_2026-09-27_source_data.csv"
QA = FIGURES / "bto_conditional_four_point_stage_2026-09-27_qa.json"
EXPECTED_SOURCE_SHA256 = "c1d2c39f0e8fae56660f0b2dc7aca385b916967b0c5fd73063cb0fcb16abbd63"


def test_center_holdout_prediction_uses_only_prior_four_corners() -> None:
    checksum = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    assert checksum == EXPECTED_SOURCE_SHA256
    qa = json.loads(QA.read_text(encoding="utf-8"))
    assert qa["source_csv_sha256"] == checksum
    assert qa["n_measured_Q_points"] == 4
    assert qa["interpolation_used"] is False
    assert qa["ABACUS_ecutwfc_Ry"] == 100
    assert qa["ABACUS_orbitals"] == "10 au DZP"
    assert qa["phonon_supercell"] == [1, 1, 1]
    assert qa["electronic_kpoints"] == [4, 4, 4]

    with SOURCE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 4
    coordinates = {
        (float(row["Q_z_sqrt_amu_A"]), float(row["Q_x_sqrt_amu_A"]))
        for row in rows
    }
    assert coordinates == {(0.6, 0.0), (0.6, 0.3), (0.9, 0.0), (0.9, 0.3)}

    def mean(field: str) -> float:
        return sum(float(row[field]) for row in rows) / 4

    assert mean("selected_E_minus_C_eV_per_BTO") == pytest.approx(
        -0.106411735308, abs=5e-13
    )
    assert mean("relaxed_Qy_zero_E_minus_C_eV_per_BTO") == pytest.approx(
        -0.063858619978, abs=5e-13
    )
    abs_qy = sum(abs(float(row["selected_Q_y_sqrt_amu_A"])) for row in rows) / 4
    assert abs_qy == pytest.approx(1.016318129613, abs=5e-13)
