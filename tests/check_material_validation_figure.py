"""Regression checks for calculator-free material-validation figure inputs."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.plot_material_validation_figure import (  # noqa: E402
    HFO2_FORMULA_UNITS_PER_CELL,
    read_completed_summary,
    write_source_data,
)


def _summary(*, climbing: bool = False) -> dict:
    return {
        "status": "completed",
        "n_images": 3,
        "path_diagnostics": {
            "geometry": {"segment_lengths_A": [2.0, 3.0]},
            "images": [
                {"image_index": 0, "relative_enthalpy_eV": 0.0, "volume_A3": 100.0},
                {"image_index": 1, "relative_enthalpy_eV": 0.12, "volume_A3": 101.0, "is_climbing_image": climbing},
                {"image_index": 2, "relative_enthalpy_eV": -0.2, "volume_A3": 102.0},
            ],
        },
    }


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        source = root / "summary.json"
        source.write_text(json.dumps(_summary(climbing=True)), encoding="utf-8")
        bto = read_completed_summary(source, label="BTO n3", material="BaTiO3")
        hfo2 = read_completed_summary(source, label="HfO2 CI n3", material="HfO2")
        if bto.coordinate != (0.0, 0.4, 1.0) or bto.barrier_eV != 0.12:
            raise SystemExit(f"path normalization failed: {bto}")
        if hfo2.climbing_image_index != 1:
            raise SystemExit("climbing-image evidence was lost")
        path_data, barrier_data = write_source_data(
            root / "out",
            bto_paths=(bto,),
            hfo2_paths=(hfo2,),
            hfo2_literature_barrier_eV_per_fu=0.032,
            bto_literature_barrier_kcal_per_mol=2.1,
        )
        if not path_data.exists() or not barrier_data.exists():
            raise SystemExit("source-data files were not written")
        hfo2_text = barrier_data.read_text(encoding="utf-8")
        expected = 0.12 / HFO2_FORMULA_UNITS_PER_CELL
        if f"{expected:.12g}" not in hfo2_text or "Literature VC-NEB" not in hfo2_text:
            raise SystemExit("per-formula-unit or literature records are missing")
    print("material_validation_figure_regression=ok")


if __name__ == "__main__":
    main()
