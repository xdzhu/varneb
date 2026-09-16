"""Regression coverage for the public ordinary-NEB force threshold."""

from __future__ import annotations

import inspect
import runpy
import sys
from pathlib import Path

from vcneb.core import run_vcneb


ROOT = Path(__file__).resolve().parents[1]
PUBLIC_CLI_DEFAULTS = (
    ("examples/run_vcneb_abacus.py", ("run", "--initial", "initial.stru", "--final", "final.stru")),
    ("examples/run_vcneb_qe.py", ("run", "--initial", "initial.stru", "--final", "final.stru")),
    ("examples/run_vcneb_vasp.py", ("run",)),
    ("examples/run_fixed_cell_ase_comparison.py", ("run",)),
    ("examples/run_vcneb_convergence.py", ("run",)),
    ("examples/run_vcneb_robustness.py", ("run",)),
    ("examples/compare_mode_path_variants.py", ("run",)),
)
PRODUCTION_TEMPLATES = (
    "cluster/hf_batio3_vcneb.slurm",
    "cluster/hf_batio3_vcneb_parallel.slurm",
    "cluster/hf_hfo2_vcneb.slurm",
    "cluster/hf_hfo2_vcneb_parallel.slurm",
    "cluster/hf_hfo2_vcneb_distributed.slurm",
)


def _parse_args(script: str, argv: tuple[str, ...], monkeypatch) -> float:
    namespace = runpy.run_path(str(ROOT / script))
    monkeypatch.setattr(sys, "argv", list(argv))
    return namespace["parse_args"]().fmax


def test_public_neb_default_force_threshold_is_0p10(monkeypatch) -> None:
    assert inspect.signature(run_vcneb).parameters["fmax"].default == 0.10
    for script, argv in PUBLIC_CLI_DEFAULTS:
        assert _parse_args(script, argv, monkeypatch) == 0.10, script


def test_production_neb_templates_default_to_0p10() -> None:
    for template in PRODUCTION_TEMPLATES:
        text = (ROOT / template).read_text(encoding="utf-8")
        assert "fmax=${FMAX:-0.10}" in text, template
