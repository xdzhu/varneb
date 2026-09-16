"""No-DFT test for creating a cross-calculator endpoint reference."""

import json
import subprocess
import sys
from pathlib import Path

from ase import Atoms
from ase.io import write


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "record_vcneb_endpoint_identity.py"


def test_endpoint_identity_recorder_writes_both_ordered_endpoint_records(tmp_path: Path) -> None:
    initial, final = tmp_path / "initial.vasp", tmp_path / "final.vasp"
    write(initial, Atoms("Ba", cell=[4, 4, 4], pbc=True), format="vasp")
    write(final, Atoms("Ba", scaled_positions=[[0.1, 0.0, 0.0]], cell=[4.1, 4, 4], pbc=True), format="vasp")
    output = tmp_path / "reference.json"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--initial", str(initial), "--final", str(final), "--output", str(output), "--label", "BTO ABACUS accepted endpoints"],
        text=True,
        capture_output=True,
        check=True,
    )
    assert "endpoint identity" in result.stdout
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["status"] == "accepted_reference_identity"
    assert payload["endpoint_structures"]["initial"]["sha256"] != payload["endpoint_structures"]["final"]["sha256"]
