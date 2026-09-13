"""Regression checks for calculator-free summary metric export."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.export_vcneb_metrics import FIELDS, rows_from_summary, write_csv_atomic


def _summary() -> dict:
    images = []
    for index in range(3):
        images.append(
            {
                "image_index": index,
                "interior": index == 1,
                "enthalpy_eV": float(index),
                "relative_enthalpy_eV": float(index),
                "volume_A3": 100.0 + index,
                "cell_lengths_A": [5.0 + index, 5.1, 5.2],
                "cell_angles_deg": [90.0, 90.0, 90.0],
            }
        )
    return {
        "status": "completed",
        "n_images": 3,
        "path_diagnostics": {
            "n_images": 3,
            "geometry": {"segment_lengths_A": [2.0, 3.0]},
            "images": images,
        },
    }


def main() -> None:
    rows = rows_from_summary(_summary())
    if len(rows) != 3 or rows[1]["reaction_coordinate"] != 0.4:
        raise SystemExit(f"reaction coordinates were exported incorrectly: {rows}")
    if rows[0]["a_A"] != 5.0 or rows[2]["c_A"] != 5.2:
        raise SystemExit("cell metrics were not exported")
    with tempfile.TemporaryDirectory() as tmp:
        output = Path(tmp) / "metrics.csv"
        write_csv_atomic(output, rows)
        lines = output.read_text(encoding="utf-8").splitlines()
        if lines[0].split(",") != list(FIELDS) or len(lines) != 4:
            raise SystemExit("metrics CSV header or row count is invalid")
    try:
        rows_from_summary({**_summary(), "status": "cancelled"})
    except ValueError:
        pass
    else:
        raise SystemExit("incomplete summary was exported as final metrics")
    print("vcneb_metrics_regression=ok")


if __name__ == "__main__":
    main()
