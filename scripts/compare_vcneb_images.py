"""Compare completed VC-NEB summaries across image counts.

The comparison is deliberately calculator-free: it checks that the summaries
use the same numerical setup, reports the barrier/profile changes, and makes
the image-count decision explicit for both barrier-bearing and monotonic paths.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdirs", nargs="+", type=Path, help="completed VC-NEB work directories")
    parser.add_argument("--barrier-tol", type=float, default=0.02, help="allowed barrier spread in eV")
    parser.add_argument(
        "--allow-duplicate-image-counts",
        action="store_true",
        help="allow variants with the same n_images (for interpolation/optimizer comparisons)",
    )
    parser.add_argument(
        "--reaction-coordinate-tol",
        type=float,
        default=0.5,
        help="allowed highest-image coordinate shift as a fraction of the local segment (default: half a segment)",
    )
    parser.add_argument("--output", type=Path, help="optional JSON report path")
    return parser.parse_args()


def _summary(path: Path) -> dict:
    summary_path = path / "vcneb_summary.json"
    try:
        return json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read {summary_path}: {exc}") from exc


def _reaction_coordinates(records: list[dict]) -> list[float]:
    lengths = []
    geometry = records[0].get("path_diagnostics", {}).get("geometry", {})
    lengths = [float(value) for value in geometry.get("segment_lengths_A", [])]
    cumulative = [0.0]
    for value in lengths:
        cumulative.append(cumulative[-1] + value)
    total = cumulative[-1] if cumulative else 0.0
    return [value / total for value in cumulative] if total > 0.0 else [0.0] * len(cumulative)


def _record(workdir: Path, summary: dict) -> dict:
    diagnostics = summary.get("path_diagnostics") or {}
    geometry = diagnostics.get("geometry") or {}
    images = diagnostics.get("images") or []
    coordinates = _reaction_coordinates([{"path_diagnostics": diagnostics}])
    highest = diagnostics.get("highest_image_index")
    highest_coordinate = coordinates[highest] if isinstance(highest, int) and highest < len(coordinates) else None
    segment = None
    if isinstance(highest, int) and 0 < highest < len(coordinates) - 1:
        segment = min(coordinates[highest] - coordinates[highest - 1], coordinates[highest + 1] - coordinates[highest])
    parameters = summary.get("calculator_parameters") or {}
    return {
        "workdir": str(workdir),
        "n_images": int(summary.get("n_images", diagnostics.get("n_images", 0))),
        "cell_interpolation": summary.get("cell_interpolation"),
        "optimizer": summary.get("optimizer"),
        "status": summary.get("status"),
        "barrier_enthalpy_eV": summary.get("barrier_enthalpy_eV"),
        "reaction_enthalpy_eV": summary.get("reaction_enthalpy_eV"),
        "final_max_generalized_force_eV_per_A": summary.get("final_max_generalized_force_eV_per_A"),
        "fmax_target_eV_per_A": summary.get("fmax_target_eV_per_A"),
        "highest_image_index": highest,
        "highest_image_coordinate": highest_coordinate,
        "local_segment_coordinate": segment,
        "has_interior_barrier": diagnostics.get("has_interior_barrier"),
        "geometry_valid": geometry.get("valid"),
        "minimum_path_distance_A": min(
            (entry.get("minimum_interatomic_distance_A") for entry in geometry.get("images", [])),
            default=None,
        ),
        "calculator_parameters": parameters,
    }


def compare(
    workdirs: list[Path],
    *,
    barrier_tol: float,
    reaction_coordinate_tol: float,
    require_unique_image_counts: bool = True,
) -> dict:
    records = [_record(workdir, _summary(workdir)) for workdir in workdirs]
    records.sort(key=lambda item: item["n_images"])
    issues: list[str] = []
    if require_unique_image_counts and len({record["n_images"] for record in records}) != len(records):
        issues.append("duplicate image counts supplied")
    reference = records[0]["calculator_parameters"]
    for record in records[1:]:
        if record["calculator_parameters"] != reference:
            issues.append(f"calculator parameters differ for {record['workdir']}")
    barriers = [record["barrier_enthalpy_eV"] for record in records if record["barrier_enthalpy_eV"] is not None]
    barrier_spread = max(barriers) - min(barriers) if barriers else None
    coordinates = [record["highest_image_coordinate"] for record in records if record["highest_image_coordinate"] is not None]
    coordinate_spread = max(coordinates) - min(coordinates) if coordinates else None
    all_barrierless = all(record["has_interior_barrier"] is False for record in records)
    if barrier_spread is not None and barrier_spread > barrier_tol:
        issues.append(f"barrier spread {barrier_spread:.6g} eV exceeds {barrier_tol:.6g} eV")
    # For a monotonic path the "highest image" is an endpoint-near sampling
    # point, not a saddle estimate; its fractional index necessarily changes
    # as images are added and must not fail the saddle-location gate.
    if coordinate_spread is not None and not all_barrierless:
        local_segments = [record["local_segment_coordinate"] for record in records if record["local_segment_coordinate"]]
        scale = min(local_segments) if local_segments else 1.0
        if coordinate_spread > reaction_coordinate_tol * scale:
            issues.append(
                f"highest-image coordinate spread {coordinate_spread:.6g} exceeds "
                f"{reaction_coordinate_tol:.6g} times the smallest local segment {scale:.6g}"
            )
    if any(record["geometry_valid"] is False for record in records):
        issues.append("at least one path has invalid geometry")
    status = "barrierless-consistent" if all_barrierless and not issues else ("ok" if not issues else "failed")
    return {
        "status": status,
        "issues": issues,
        "barrier_tolerance_eV": barrier_tol,
        "reaction_coordinate_tolerance_fraction": reaction_coordinate_tol,
        "image_count_uniqueness_required": require_unique_image_counts,
        "barrier_spread_eV": barrier_spread,
        "highest_image_coordinate_spread": coordinate_spread,
        "records": records,
    }


def write_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def main() -> int:
    args = parse_args()
    try:
        report = compare(
            args.workdirs,
            barrier_tol=args.barrier_tol,
            reaction_coordinate_tol=args.reaction_coordinate_tol,
            require_unique_image_counts=not args.allow_duplicate_image_counts,
        )
    except ValueError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    output = args.output or (args.workdirs[0] / "vcneb_image_comparison.json")
    write_atomic(output, report)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if report["status"] in {"ok", "barrierless-consistent"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
