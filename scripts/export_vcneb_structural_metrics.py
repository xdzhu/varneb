"""Export key bond lengths and optional mode projections from a VC-NEB trajectory.

The command reads only the latest complete chain snapshot. A supplied mode is
projected with the same extended atomic/cell metric as VC-NEB;
``--derive-endpoint-mode`` is a structural endpoint-displacement diagnostic,
not a phonon eigenvector and must not be described as one.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import os
from pathlib import Path
import sys
import tempfile
from typing import Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.core import deformation_from_cell, path_geometry_diagnostics, read_chain_trajectory
from vcneb.modes import Mode, project_path_onto_modes


@dataclass(frozen=True)
class Pair:
    label: str
    first: int
    second: int


def parse_pair(spec: str, *, n_atoms: int) -> Pair:
    """Parse ``label:i,j`` (zero-based atom indices)."""

    if ":" in spec:
        label, indices = spec.split(":", 1)
    else:
        indices = spec
        label = f"pair_{spec.replace(',', '_')}"
    fields = [field.strip() for field in indices.split(",")]
    if len(fields) != 2 or not label.strip():
        raise ValueError(f"pair must have the form label:i,j, got {spec!r}")
    try:
        first, second = (int(field) for field in fields)
    except ValueError as exc:
        raise ValueError(f"pair indices must be integers, got {spec!r}") from exc
    if first == second or min(first, second) < 0 or max(first, second) >= n_atoms:
        raise ValueError(f"pair indices out of range for {n_atoms} atoms: {spec!r}")
    return Pair(label.strip(), first, second)


def _reaction_coordinates(segments: Sequence[float]) -> list[float]:
    cumulative = [0.0]
    for length in segments:
        cumulative.append(cumulative[-1] + float(length))
    if cumulative[-1] <= 0.0:
        raise ValueError("trajectory has zero total extended-coordinate path length")
    return [value / cumulative[-1] for value in cumulative]


def _atomic_write_csv(path: Path, rows: list[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
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


def export_metrics(
    trajectory: str | Path,
    *,
    n_images: int,
    output: str | Path,
    pairs: Sequence[Pair] = (),
    mode: Mode | None = None,
    mode_source: str | None = None,
) -> dict[str, object]:
    """Export metrics and return the machine-readable provenance payload."""

    if n_images < 2:
        raise ValueError("n_images must be at least 2")
    images = read_chain_trajectory(trajectory, n_images=n_images, step=-1)
    geometry = path_geometry_diagnostics(images)
    coordinates = _reaction_coordinates(geometry["segment_lengths_A"])
    projections = None
    if mode is not None:
        projections = project_path_onto_modes(
            images,
            images[0],
            mode,
            cell_scale=float(geometry["cell_scale_A"]),
        )
    fields = ["image_index", "reaction_coordinate"]
    fields.extend(pair.label for pair in pairs)
    if projections is not None:
        fields.extend(f"mode_coefficient_{index}" for index in range(projections.shape[1]))

    rows: list[dict[str, object]] = []
    for image_index, image in enumerate(images):
        row: dict[str, object] = {
            "image_index": image_index,
            "reaction_coordinate": coordinates[image_index],
        }
        for pair in pairs:
            row[pair.label] = float(image.get_distance(pair.first, pair.second, mic=True))
        if projections is not None:
            for mode_index, coefficient in enumerate(projections[image_index]):
                row[f"mode_coefficient_{mode_index}"] = float(coefficient)
        rows.append(row)

    destination = Path(output)
    _atomic_write_csv(destination, rows, fields)
    metadata: dict[str, object] = {
        "status": "ok",
        "trajectory": str(Path(trajectory).resolve()),
        "n_images": n_images,
        "pairs": [
            {"label": pair.label, "first": pair.first, "second": pair.second}
            for pair in pairs
        ],
        "mode_source": mode_source,
        "mode_semantics": (
            "user-supplied mode projected in the extended atomic/cell metric"
            if mode is not None and mode_source != "endpoint_displacement"
            else "endpoint displacement diagnostic; not a phonon eigenvector"
            if mode is not None
            else None
        ),
        "cell_scale_A": float(geometry["cell_scale_A"]),
        "output_csv": str(destination.resolve()),
    }
    destination.with_suffix(".json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return metadata


def _derive_endpoint_mode(images) -> Mode:
    initial, final = images[0], images[-1]
    atomic = np.asarray(final.positions - initial.positions, dtype=float)
    deform = deformation_from_cell(final.cell.array, initial.cell.array) - np.eye(3)
    return Mode(atomic=atomic, cell=deform, label="endpoint-displacement")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path, help="VC-NEB trajectory containing complete image groups")
    parser.add_argument("--n-images", type=int, required=True, help="total images, including fixed endpoints")
    parser.add_argument("--pair", action="append", default=[], help="key pair as label:i,j with zero-based indices")
    parser.add_argument("--mode", type=Path, help="JSON/NPZ/text mode file to project")
    parser.add_argument(
        "--derive-endpoint-mode",
        action="store_true",
        help="project the endpoint displacement (structural diagnostic, not a phonon mode)",
    )
    parser.add_argument("--output", type=Path, required=True, help="output CSV")
    args = parser.parse_args()
    if args.mode is not None and args.derive_endpoint_mode:
        parser.error("--mode and --derive-endpoint-mode are mutually exclusive")
    try:
        images = read_chain_trajectory(args.trajectory, n_images=args.n_images, step=-1)
        pairs = [parse_pair(spec, n_atoms=len(images[0])) for spec in args.pair]
        mode = None
        mode_source = None
        if args.mode is not None:
            mode = Mode.from_file(args.mode, n_atoms=len(images[0]))
            mode_source = str(args.mode.resolve())
        elif args.derive_endpoint_mode:
            mode = _derive_endpoint_mode(images)
            mode_source = "endpoint_displacement"
        metadata = export_metrics(
            args.trajectory,
            n_images=args.n_images,
            output=args.output,
            pairs=pairs,
            mode=mode,
            mode_source=mode_source,
        )
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.error(str(exc))
    print(json.dumps(metadata, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
