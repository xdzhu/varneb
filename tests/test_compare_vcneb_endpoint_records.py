"""No-DFT test for strict cross-calculator endpoint comparison."""

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "compare_vcneb_endpoint_records.py"


def _preflight(path: Path, initial: str, final: str) -> None:
    path.write_text(
        json.dumps({"endpoint_structures": {"initial": {"sha256": initial}, "final": {"sha256": final}}}),
        encoding="utf-8",
    )


def test_endpoint_record_gate_writes_auditable_match_and_blocks_mismatch(tmp_path: Path) -> None:
    reference, candidate = tmp_path / "abacus.json", tmp_path / "qe.json"
    _preflight(reference, "a", "b")
    _preflight(candidate, "a", "b")
    output = tmp_path / "gate.json"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--reference-preflight", str(reference), "--candidate-preflight", str(candidate), "--output", str(output)],
        text=True,
        capture_output=True,
        check=True,
    )
    assert "identities match" in result.stdout
    assert json.loads(output.read_text(encoding="utf-8"))["matches"]
    _preflight(candidate, "a", "different")
    blocked = subprocess.run(
        [sys.executable, str(SCRIPT), "--reference-preflight", str(reference), "--candidate-preflight", str(candidate), "--output", str(output)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert blocked.returncode != 0
    assert not json.loads(output.read_text(encoding="utf-8"))["matches"]
