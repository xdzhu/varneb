"""Regression check for trajectory-only structural metric export."""

from __future__ import annotations

import csv
from pathlib import Path
import sys
import tempfile

import numpy as np
from ase import Atoms
from ase.io.trajectory import Trajectory

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.export_vcneb_structural_metrics import Pair, export_metrics
from vcneb.modes import Mode


def main() -> None:
    cell = np.diag([5.0, 5.0, 5.0])
    images = [
        Atoms("Ar2", positions=[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], cell=cell, pbc=True),
        Atoms("Ar2", positions=[[0.0, 0.0, 0.0], [2.5, 0.2, 0.0]], cell=cell, pbc=True),
        Atoms("Ar2", positions=[[0.0, 0.0, 0.0], [3.0, 0.4, 0.0]], cell=cell, pbc=True),
    ]
    mode = Mode([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    with tempfile.TemporaryDirectory(prefix="vcneb-geometry-metrics-") as tmp:
        root = Path(tmp)
        trajectory = root / "chain.traj"
        with Trajectory(str(trajectory), "w") as handle:
            for image in images:
                handle.write(image)
        output = root / "metrics.csv"
        metadata = export_metrics(
            trajectory,
            n_images=3,
            output=output,
            pairs=[Pair("Ar0-Ar1", 0, 1)],
            mode=mode,
            mode_source="test-mode",
        )
        with output.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != 3 or rows[1]["Ar0-Ar1"] == rows[0]["Ar0-Ar1"]:
            raise SystemExit("key-pair distances were not exported per image")
        if "mode_coefficient_0" not in rows[1] or metadata["status"] != "ok":
            raise SystemExit("mode projection or metadata was not exported")
        if not output.with_suffix(".json").exists():
            raise SystemExit("structural metric provenance JSON is missing")
    print("geometry_metrics_regression=ok")


if __name__ == "__main__":
    main()
