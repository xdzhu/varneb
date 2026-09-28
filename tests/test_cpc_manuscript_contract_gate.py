"""Include the CPC claim-boundary check in the normal repository test suite."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_cpc_manuscript_contract_gate() -> None:
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, str(root / "tests" / "check_cpc_manuscript_contract.py")],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "cpc_manuscript_contract=ok" in result.stdout
