"""Compare ordered historical chains; geometry alone does not identify modes."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read

from vcneb.analysis import path_reaction_coordinate
from vcneb.provenance import endpoint_structure_record


def chain_record(directory: Path) -> dict:
    files = sorted(directory.glob("POSCAR_*"))
    if len(files) != 7:
        raise ValueError("this audit requires a complete seven-total-image snapshot")
    images = [read(path, format="vasp") for path in files]
    if any(image.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8 for image in images):
        raise ValueError("chain must preserve ordered Hf4O8")
    if any(image.get_volume() <= 0 or not np.isfinite(image.positions).all() for image in images):
        raise ValueError("invalid chain geometry")
    wrapped = np.array([image.get_scaled_positions(wrap=True) for image in images])
    unwrapped = [wrapped[0].copy()]
    for previous, current in zip(wrapped[:-1], wrapped[1:]):
        delta = current - previous
        delta -= np.rint(delta)
        unwrapped.append(unwrapped[-1] + delta)
    unwrapped = np.array(unwrapped)
    winding = np.rint(unwrapped[-1] - wrapped[-1]).astype(int)
    reference_cell = images[0].cell.array
    displacement = (unwrapped - unwrapped[0]) @ reference_cell
    translation = np.average(displacement, weights=images[0].get_masses(), axis=1)
    displacement -= translation[:, None, :]
    coordinate, segment = path_reaction_coordinate(images)
    cells = np.array([image.cell.array for image in images])
    return {
        "source_files": [{"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                         for path in files],
        "endpoints": [endpoint_structure_record(images[0]), endpoint_structure_record(images[-1])],
        "n_images_total": len(images),
        "periodic_unwrap_policy": "nearest fractional increment at every link; no atom remapping",
        "maximum_fractional_link_component_before_unwrap": float(np.max(np.abs(np.diff(wrapped, axis=0)))),
        "final_unwrap_integers": winding.tolist(),
        "reaction_coordinate_A": coordinate.tolist(),
        "segment_lengths_A": segment.tolist(),
        "minimum_segment_A": float(segment.min()),
        "maximum_segment_A": float(segment.max()),
        "translation_free_atomic_displacements_reference_A": displacement.tolist(),
        "removed_mass_weighted_translations_A": translation.tolist(),
        "cells_A": cells.tolist(),
        "volumes_A3": [float(image.get_volume()) for image in images],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ordinary", type=Path, required=True)
    parser.add_argument("--guided", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    left, right = chain_record(args.ordinary), chain_record(args.guided)
    endpoint_equal = [a["sha256"] == b["sha256"] for a, b in zip(left["endpoints"], right["endpoints"])]
    sampled_arc = np.linspace(0, 1, 101)
    curves = []
    for record in (left, right):
        s = np.array(record["reaction_coordinate_A"])
        if s[-1] <= 0 or np.any(np.diff(s) <= 0):
            raise ValueError("degenerate geometric path cannot be compared by arc length")
        atomic = np.array(record["translation_free_atomic_displacements_reference_A"]).reshape(7, -1)
        curves.append(np.array([np.interp(sampled_arc, s / s[-1], atomic[:, index])
                                for index in range(atomic.shape[1])]).T.reshape(101, 12, 3))
    rms = np.sqrt(np.mean(np.sum((curves[0] - curves[1]) ** 2, axis=2), axis=1))
    report = {
        "schema_version": 1,
        "scope": "ordered geometry and nearest-link unwrap audit; no phonon labels or topology certificate",
        "ordered_endpoint_hashes_equal": endpoint_equal,
        "final_nearest_link_unwrap_integers_equal": left["final_unwrap_integers"] == right["final_unwrap_integers"],
        "maximum_atomic_RMS_difference_at_equal_normalized_arc_A": float(rms.max()),
        "arc_interpolation_is_geometric_only_not_new_DFT": True,
        "interpretation": "different geometries or winding are observations, not proof of X2 sign reversal",
        "ordinary": left, "guided": right,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "ordered_endpoint_hashes_equal", "final_nearest_link_unwrap_integers_equal",
        "maximum_atomic_RMS_difference_at_equal_normalized_arc_A")}, indent=2))


if __name__ == "__main__":
    main()
