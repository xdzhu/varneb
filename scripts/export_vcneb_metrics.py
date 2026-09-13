"""Export per-image VC-NEB geometry and force metrics from a summary JSON."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import tempfile


FIELDS = (
    "image_index",
    "interior",
    "reaction_coordinate",
    "enthalpy_eV",
    "relative_enthalpy_eV",
    "volume_A3",
    "a_A",
    "b_A",
    "c_A",
    "alpha_deg",
    "beta_deg",
    "gamma_deg",
    "max_atom_force_eV_per_A",
    "max_stress_eV_per_A3",
    "max_cell_force_eV",
    "max_true_generalized_force_eV_per_A",
    "neb_residual_generalized_force_eV_per_A",
    "is_climbing_image",
)


def _reaction_coordinates(segment_lengths: list[float], n_images: int) -> list[float]:
    if len(segment_lengths) != max(0, n_images - 1):
        raise ValueError("summary geometry segment_lengths_A does not match n_images")
    cumulative = [0.0]
    for length in segment_lengths:
        cumulative.append(cumulative[-1] + float(length))
    total = cumulative[-1]
    if total <= 0.0:
        raise ValueError("summary geometry has zero total path length")
    return [value / total for value in cumulative]


def rows_from_summary(summary: dict) -> list[dict[str, object]]:
    if summary.get("status") != "completed":
        raise ValueError("only completed VC-NEB summaries can be exported")
    diagnostics = summary.get("path_diagnostics") or {}
    images = diagnostics.get("images") or []
    n_images = int(summary.get("n_images", diagnostics.get("n_images", len(images))))
    if len(images) != n_images:
        raise ValueError("summary path_diagnostics.images does not match n_images")
    geometry = diagnostics.get("geometry") or {}
    coordinates = _reaction_coordinates(
        [float(value) for value in geometry.get("segment_lengths_A", [])], n_images
    )
    rows: list[dict[str, object]] = []
    for record in images:
        index = int(record["image_index"])
        lengths = list(record.get("cell_lengths_A") or [None, None, None])
        angles = list(record.get("cell_angles_deg") or [None, None, None])
        rows.append(
            {
                "image_index": index,
                "interior": bool(record.get("interior", 0 < index < n_images - 1)),
                "reaction_coordinate": coordinates[index],
                "enthalpy_eV": record.get("enthalpy_eV"),
                "relative_enthalpy_eV": record.get("relative_enthalpy_eV"),
                "volume_A3": record.get("volume_A3"),
                "a_A": lengths[0],
                "b_A": lengths[1],
                "c_A": lengths[2],
                "alpha_deg": angles[0],
                "beta_deg": angles[1],
                "gamma_deg": angles[2],
                "max_atom_force_eV_per_A": record.get("max_atom_force_eV_per_A"),
                "max_stress_eV_per_A3": record.get("max_stress_eV_per_A3"),
                "max_cell_force_eV": record.get("max_cell_force_eV"),
                "max_true_generalized_force_eV_per_A": record.get(
                    "max_true_generalized_force_eV_per_A"
                ),
                "neb_residual_generalized_force_eV_per_A": record.get(
                    "neb_residual_generalized_force_eV_per_A"
                ),
                "is_climbing_image": bool(record.get("is_climbing_image", False)),
            }
        )
    return rows


def write_csv_atomic(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path, help="completed vcneb_summary.json")
    parser.add_argument("--output", type=Path, help="output CSV (default: next to summary)")
    args = parser.parse_args()
    try:
        summary = json.loads(args.summary.read_text(encoding="utf-8"))
        rows = rows_from_summary(summary)
    except (OSError, json.JSONDecodeError, ValueError, KeyError) as exc:
        parser.error(str(exc))
    output = args.output or args.summary.with_name("vcneb_metrics.csv")
    write_csv_atomic(output, rows)
    print(f"wrote={output.resolve()} rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
