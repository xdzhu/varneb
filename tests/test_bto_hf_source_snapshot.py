"""The paper's BTO DFT provenance must retain its exact executed source."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / "benchmarks" / "numerical_integrity" / "bto_hf_source_snapshot_2026-09-27"


@pytest.mark.parametrize("relative,expected", [
    ("examples/run_bto_q1q2_conditional_abacus.py",
     "de0c1c3aeac376b9e4f4ea1f9ff742e85e70f9ba08823824c46c358d547a17bf"),
    ("scripts/audit_bto_q1q2_conditional_pilot.py",
     "afe51053a1c23f119b256246cbe21e97529c6da3c9e037dec45895b360708f29"),
    ("vcneb/mode_surface.py",
     "91b70d48a488b5a26f1e563bbce1ba0e75581a01d9145158d57108b8dac4bbfc"),
    ("vcneb/__init__.py",
     "25576062e4c39ea92e552fa78c9e7dd554ce00d8d7667d8a7404eb58245c1868"),
    ("vcneb/core.py",
     "c69e4066e314fa5f81b30c8eb3f66e3dfe327d9d62abdd718369a7932f920169"),
])
def test_executed_snapshot_is_byte_identical(relative: str, expected: str) -> None:
    assert hashlib.sha256((SNAPSHOT / relative).read_bytes()).hexdigest() == expected
