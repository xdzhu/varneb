"""Validate endpoint-mode derivation on the tracked 12-atom HfO2 inputs."""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import Mode
from scripts.derive_endpoint_mode import derive_endpoint_mode


def main() -> None:
    initial = ROOT / "validation" / "hfo2_t_to_po" / "T_HfO2_12.vasp"
    final = ROOT / "validation" / "hfo2_t_to_po" / "PO_HfO2_12_mapped.vasp"
    payload = derive_endpoint_mode(initial, final)
    if payload["n_atoms"] != 12 or payload["symbols"] != ["Hf"] * 4 + ["O"] * 8:
        raise SystemExit("endpoint mode atom metadata is inconsistent")
    if len(payload["atomic"]) != 12 or len(payload["cell"]) != 3:
        raise SystemExit("endpoint mode dimensions are invalid")
    if payload["semantics"] != "endpoint displacement diagnostic; not a phonon eigenvector":
        raise SystemExit("endpoint mode semantics were not recorded")
    Mode(payload["atomic"], cell=payload["cell"])
    if payload["mapping"] is None or len(payload["mapping"]) != 12:
        raise SystemExit("endpoint mode mapping metadata is missing")
    print("endpoint_mode_regression=ok")


if __name__ == "__main__":
    main()
