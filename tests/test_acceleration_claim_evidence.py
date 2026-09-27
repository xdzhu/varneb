"""Reconstruct the CPC acceleration table from archived optimizer traces.

These are evaluation *estimates* from the recorded update count, not an
independent audit of every DFT launch. See the manuscript evidence audit.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "benchmarks" / "convergence" / "results"
ROW = re.compile(
    r"^(?:FIRE|BlockFIRE):\s+(\d+)\s+\d\d:\d\d:\d\d\s+"
    r"[-+\d.eE]+\s+([-+\d.eE]+)\s*$"
)


def _rows(name: str) -> list[tuple[int, float]]:
    rows = []
    for line in (RESULTS / name).read_text(encoding="utf-8").splitlines():
        if match := ROW.match(line):
            rows.append((int(match[1]), float(match[2])))
    assert rows and [step for step, _ in rows] == list(range(len(rows)))
    return rows


@pytest.mark.parametrize(
    ("log_name", "n_total", "first_step", "reported_evaluations", "initial_force"),
    [
        ("bto_fire_baseline.log", 9, 28, 205, 0.438473),
        ("bto_blockfire_002.log", 9, 10, 79, 0.438473),
        ("bto_blockfire_001.log", 9, 5, 44, 0.438473),
        ("bto_fire_scaled_002.log", 9, 8, 65, 0.438473),
        ("hfo2_fire_baseline.log", 7, 42, 217, 0.965655),
    ],
)
def test_first_crossing_and_reconstructed_evaluations(
    log_name: str, n_total: int, first_step: int,
    reported_evaluations: int, initial_force: float,
) -> None:
    rows = _rows(log_name)
    assert rows[0][1] == initial_force
    assert all(force > 0.10 for _, force in rows[:first_step])
    assert rows[first_step][1] <= 0.10
    # Two endpoints are evaluated once at startup, while only interior
    # images are evaluated at step 0 and each subsequent accepted update.
    estimated = 2 + (n_total - 2) * (first_step + 1)
    assert estimated == reported_evaluations


def test_staged_hfo2_count_and_continuity() -> None:
    coarse = _rows("hfo2_fire_scaled.log")
    refine = _rows("hfo2_fire_staged_refine23.log")
    assert coarse[23][1] == refine[0][1] == 0.110833
    assert all(force > 0.10 for _, force in coarse[:24])
    assert refine[1][1] == 0.092367
    # Both independent run startups are included in this historical estimate.
    n_interior = 7 - 2
    estimated = (2 + n_interior * (23 + 1)) + (2 + n_interior * (1 + 1))
    assert estimated == 134


def test_saved_percentages_are_case_specific() -> None:
    assert round(100 * (205 - 44) / 205, 1) == 78.5
    assert round(100 * (217 - 134) / 217, 1) == 38.2
