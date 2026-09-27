"""A four-point stage figure must never masquerade as a continuous BTO PES."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest


FIGURES = Path(__file__).resolve().parents[1] / "paper/VARNEB_CPC/figures"
PREFIX = FIGURES / "bto_conditional_four_point_stage_2026-09-27"


def test_source_table_and_claim_limit():
    csv_path = Path(str(PREFIX) + "_source_data.csv")
    qa = json.loads(Path(str(PREFIX) + "_qa.json").read_text(encoding="utf-8"))
    with csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert qa["kind"] == "BTO_four_measured_fixed_Q_branch_candidates_not_continuous_PES_or_barrier"
    assert qa["n_measured_Q_points"] == len(rows) == 4
    assert qa["interpolation_used"] is False
    assert qa["phonon_supercell"] == [1, 1, 1]
    assert qa["electronic_kpoints"] == [4, 4, 4]
    assert qa["ABACUS_ecutwfc_Ry"] == 100
    assert qa["source_csv_sha256"] == hashlib.sha256(csv_path.read_bytes()).hexdigest()
    assert {(float(row["Q_z_sqrt_amu_A"]), float(row["Q_x_sqrt_amu_A"]))
            for row in rows} == {(0.6, 0.0), (0.6, 0.3), (0.9, 0.0), (0.9, 0.3)}
    expected_lowerings = [57.46661343, 49.98140050, 34.45637062, 28.30807677]
    for row, expected in zip(rows, expected_lowerings):
        selected = float(row["selected_E_minus_C_eV_per_BTO"])
        constrained = float(row["relaxed_Qy_zero_E_minus_C_eV_per_BTO"])
        assert 1000 * (constrained - selected) == pytest.approx(expected, abs=1e-6)
        assert len(row["result_sha256"]) == len(row["replay_sha256"]) == 64
    for extension in (".png", ".pdf", ".svg"):
        assert PREFIX.with_suffix(extension).is_file()
