"""Regression check for the optional metrics plotting helper."""

from __future__ import annotations

import csv
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main() -> None:
    try:
        from scripts.plot_vcneb_metrics import plot_metrics
    except ImportError as exc:
        raise SystemExit(f"plotting import failed: {exc}") from exc

    with tempfile.TemporaryDirectory(prefix="vcneb-plot-") as temporary:
        directory = Path(temporary)
        metrics = directory / "metrics.csv"
        fields = [
            "reaction_coordinate",
            "relative_enthalpy_eV",
            "volume_A3",
            "a_A",
            "b_A",
            "c_A",
            "max_true_generalized_force_eV_per_A",
            "neb_residual_generalized_force_eV_per_A",
            "is_climbing_image",
        ]
        with metrics.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for coordinate, energy in ((0.0, 0.0), (0.5, 0.2), (1.0, 0.0)):
                writer.writerow(
                    {
                        "reaction_coordinate": coordinate,
                        "relative_enthalpy_eV": energy,
                        "volume_A3": 100.0 + coordinate,
                        "a_A": 4.0 + coordinate,
                        "b_A": 4.1 + coordinate,
                        "c_A": 4.2 + coordinate,
                        "max_true_generalized_force_eV_per_A": 0.1,
                        "neb_residual_generalized_force_eV_per_A": 0.05,
                        "is_climbing_image": coordinate == 0.5,
                    }
                )
        output = directory / "plot.png"
        plot_metrics([("toy", metrics)], output)
        if not output.exists() or output.stat().st_size <= 0:
            raise SystemExit("metrics plotting did not create a non-empty output")
    print("vcneb_plot_regression=ok")


if __name__ == "__main__":
    main()
