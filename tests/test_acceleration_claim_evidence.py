"""Reconstruct the CPC acceleration table from archived optimizer traces.

The formula here is an update-count estimate. The separate archived Slurm
accounting test corroborates actual ABACUS process launches; neither test
audits every SCF history. See the manuscript evidence audit.
"""

from __future__ import annotations

import json
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


def _arguments(name: str) -> list[str]:
    summary = json.loads((RESULTS / f"{name}_summary.json").read_text(encoding="utf-8"))
    return summary["command_line"]


def _values(arguments: list[str], option: str, count: int = 1) -> tuple[tuple[str, ...], ...]:
    return tuple(
        tuple(arguments[index + 1:index + 1 + count])
        for index, value in enumerate(arguments)
        if value == option
    )


@pytest.mark.parametrize(
    "names",
    [
        ("bto_fire_baseline", "bto_blockfire_001", "bto_blockfire_002", "bto_fire_scaled_002"),
        ("hfo2_fire_baseline", "hfo2_fire_scaled", "hfo2_fire_staged_refine23"),
    ],
)
def test_archived_calculator_contract_is_matched(names: tuple[str, ...]) -> None:
    arguments = [_arguments(name) for name in names]
    for option in (
        "--pseudo-dir", "--basis-dir", "--pp", "--basis", "--ecutwfc",
        "--scf-thr", "--scf-nmax", "--mixing-beta", "--command", "--n-images", "--k",
    ):
        reference = _values(arguments[0], option)
        assert reference
        assert all(_values(item, option) == reference for item in arguments)
    k_mesh = _values(arguments[0], "--kpts", 3)
    assert k_mesh
    assert all(_values(item, "--kpts", 3) == k_mesh for item in arguments)


def test_common_analysis_threshold_is_not_misreported_as_common_run_stop() -> None:
    assert _values(_arguments("bto_fire_baseline"), "--fmax") == (("0.02",),)
    assert _values(_arguments("hfo2_fire_baseline"), "--fmax") == (("0.05",),)
    for name in ("bto_blockfire_001", "bto_blockfire_002", "bto_fire_scaled_002",
                 "hfo2_fire_scaled", "hfo2_fire_staged_refine23"):
        assert _values(_arguments(name), "--fmax") == (("0.10",),)
    manuscript = (ROOT / "paper/VARNEB_CPC/varneb_CPC.tex").read_text(encoding="utf-8")
    assert "first crossing" in manuscript
    assert "historical" in manuscript
    assert "stopping rule fixed" not in manuscript


def test_hfo2_acceleration_reuses_the_linear_seed_despite_cli_label() -> None:
    baseline = json.loads((RESULTS / "hfo2_fire_baseline_summary.json").read_text(encoding="utf-8"))
    coarse = _arguments("hfo2_fire_scaled")
    refine = _arguments("hfo2_fire_staged_refine23")
    assert _values(baseline["command_line"], "--cell-interpolation") == (("linear",),)
    assert _values(coarse, "--cell-interpolation") == (("log_strain",),)
    assert _values(refine, "--cell-interpolation") == (("log_strain",),)
    assert "--resume" in coarse and "--resume" in refine
    assert _values(coarse, "--resume-trajectory") == (
        (baseline["workdir"] + "/initial-vcneb.traj",),
    )
    assert _values(refine, "--resume-trajectory")[0][0].endswith(
        "/hfo2_fire_scaled/snapshots/chain_step_0023.traj"
    )
