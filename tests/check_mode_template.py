"""Validate the copy-and-edit JSON mode template without running a calculator."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.modes import Mode


def main() -> None:
    mode = Mode.from_file(ROOT / "examples" / "mode_template.json")
    assert mode.label == "replace-with-mode-label"
    assert mode.frequency is None
    assert mode.atomic.shape == (1, 3)
    assert mode.cell is not None and mode.cell.shape == (3, 3)
    assert np.allclose(mode.atomic, [[1.0, 0.0, 0.0]])
    assert np.allclose(mode.cell, 0.0)
    print("mode_template_regression=ok")


if __name__ == "__main__":
    main()
