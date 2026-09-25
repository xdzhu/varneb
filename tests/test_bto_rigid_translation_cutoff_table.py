"""Keep the public BTO numerical-sensitivity source table self-consistent."""

import csv
from pathlib import Path

import pytest


def test_bto_rigid_z_cutoff_table_matches_raw_audited_spans():
    path = (Path(__file__).resolve().parents[1] / "benchmarks" /
            "numerical_integrity" / "bto_q060_rigid_z_ecutwfc100_120.csv")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [float(row["shift_z_A"]) for row in rows] == [-0.02, -0.01, 0.0, 0.01, 0.02]
    for cutoff, expected in ((100, 0.02146751330656116),
                             (120, 0.020593628505594097)):
        values = [float(row[f"energy_{cutoff}Ry_eV"]) for row in rows]
        assert 1000 * (max(values) - min(values)) == pytest.approx(expected, abs=1e-6)
