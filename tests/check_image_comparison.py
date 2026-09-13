"""Regression checks for image-count convergence diagnostics."""

from __future__ import annotations

import sys
import json
import tempfile

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.compare_vcneb_images import compare


def _summary(n_images: int, *, barrier: float, peak_index: int, segment: float) -> dict:
    lengths = [segment] * (n_images - 1)
    return {
        "n_images": n_images,
        "status": "completed",
        "barrier_enthalpy_eV": barrier,
        "reaction_enthalpy_eV": -0.3,
        "final_max_generalized_force_eV_per_A": 0.04,
        "fmax_target_eV_per_A": 0.05,
        "calculator_parameters": {"ecutwfc": 100.0, "kpts": [2, 2, 2]},
        "path_diagnostics": {
            "n_images": n_images,
            "highest_image_index": peak_index,
            "has_interior_barrier": True,
            "geometry": {
                "valid": True,
                "segment_lengths_A": lengths,
                "images": [
                    {"minimum_interatomic_distance_A": 2.1}
                    for _ in range(n_images)
                ],
            },
        },
    }


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        close_dirs = []
        for name, summary in (
            ("n7", _summary(7, barrier=0.1567, peak_index=2, segment=1.0)),
            ("n9", _summary(9, barrier=0.1562, peak_index=3, segment=0.75)),
        ):
            directory = root / name
            directory.mkdir()
            (directory / "vcneb_summary.json").write_text(json.dumps(summary), encoding="utf-8")
            close_dirs.append(directory)
        close = compare(close_dirs, barrier_tol=0.02, reaction_coordinate_tol=0.5)
        if close["status"] != "ok" or close["issues"]:
            raise SystemExit(f"discrete peak shift was rejected unexpectedly: {close}")

        far_summary = _summary(7, barrier=0.1567, peak_index=1, segment=1.0)
        (close_dirs[0] / "vcneb_summary.json").write_text(json.dumps(far_summary), encoding="utf-8")
        far = compare(close_dirs, barrier_tol=0.02, reaction_coordinate_tol=0.5)
        if far["status"] != "failed" or not any("coordinate spread" in issue for issue in far["issues"]):
            raise SystemExit("large saddle-coordinate shift was not rejected")
    print("image_comparison_regression=ok")


if __name__ == "__main__":
    main()
