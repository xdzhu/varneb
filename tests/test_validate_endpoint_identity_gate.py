"""No-DFT tests for the production identity-gate checker."""

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "validate_endpoint_identity_gate.py"


def test_identity_gate_checker_accepts_only_two_passing_endpoint_records(tmp_path: Path) -> None:
    gate = tmp_path / "gate.json"
    gate.write_text(json.dumps({"matches": True, "endpoints": {"initial": {"matches": True}, "final": {"matches": True}}}), encoding="utf-8")
    passed = subprocess.run([sys.executable, str(SCRIPT), "--gate", str(gate)], text=True, capture_output=True, check=True)
    assert "passed" in passed.stdout
    gate.write_text(json.dumps({"matches": False, "endpoints": {"initial": {"matches": True}, "final": {"matches": False}}}), encoding="utf-8")
    blocked = subprocess.run([sys.executable, str(SCRIPT), "--gate", str(gate)], text=True, capture_output=True, check=False)
    assert blocked.returncode != 0
    assert "BLOCKED" in blocked.stderr
