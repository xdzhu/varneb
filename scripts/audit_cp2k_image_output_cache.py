"""Read-only comparison of cumulative CP2K FORCE_EVAL output and VARNEB cache.

Run with the original per-image ``cp2k.out`` and its parent run's
``image_cache`` directory. Cache entries are ordered by modification time;
the order is diagnostic only and is not assumed to be a guaranteed call ID.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.units import Hartree


ENERGY_RE = re.compile(
    r"^\s*ENERGY\| Total FORCE_EVAL \( QS \) energy \[a\.u\.\]:\s*"
    r"([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s*$",
    re.MULTILINE,
)


def compare(output: Path, cache_dir: Path, image_index: int, tolerance_ev: float) -> dict:
    raw = output.read_bytes()
    text = raw.decode("utf-8", errors="replace")
    output_energies_ha = [float(value) for value in ENERGY_RE.findall(text)]
    if not output_energies_ha:
        raise ValueError(f"no CP2K FORCE_EVAL energy in {output}")
    output_energies = [value * Hartree for value in output_energies_ha]

    cache_paths = sorted(
        cache_dir.glob(f"image_{image_index:04d}_*.npz"),
        key=lambda path: (path.stat().st_mtime_ns, path.name),
    )
    if not cache_paths:
        raise ValueError(f"no image {image_index} cache entries in {cache_dir}")
    cached = []
    for path in cache_paths:
        with np.load(path, allow_pickle=False) as data:
            energy = float(np.asarray(data["energy"]).reshape(()))
        cached.append({
            "name": path.name,
            "energy_ev": energy,
            "mtime_ns": path.stat().st_mtime_ns,
        })

    out = np.asarray(output_energies)
    missing = []
    for row in cached:
        nearest = int(np.argmin(np.abs(out - row["energy_ev"])))
        difference = float(out[nearest] - row["energy_ev"])
        if abs(difference) > tolerance_ev:
            missing.append({**row, "nearest_output_record": nearest + 1,
                            "nearest_output_minus_cache_ev": difference})
    cache_values = np.asarray([row["energy_ev"] for row in cached])
    unmatched_output = []
    for index, energy in enumerate(output_energies, start=1):
        nearest = int(np.argmin(np.abs(cache_values - energy)))
        difference = float(cache_values[nearest] - energy)
        if abs(difference) > tolerance_ev:
            unmatched_output.append({
                "record": index,
                "energy_ev": energy,
                "nearest_cache": cached[nearest]["name"],
                "nearest_cache_minus_output_ev": difference,
            })
    aligned = min(len(cached), len(out))
    order_differences = [float(out[i] - cached[i]["energy_ev"])
                         for i in range(aligned)]
    paired_ha = np.asarray(output_energies_ha[:aligned])
    paired_cache = cache_values[:aligned]
    effective_factor = float(np.dot(paired_ha, paired_cache)
                             / np.dot(paired_ha, paired_ha))
    factor_residuals = paired_cache - effective_factor * paired_ha
    return {
        "output_sha256": hashlib.sha256(raw).hexdigest(),
        "image_index": image_index,
        "tolerance_ev": tolerance_ev,
        "n_output_records": len(output_energies),
        "n_unique_cache_entries": len(cached),
        "n_matching_cache_entries_with_ase_constant": len(cached) - len(missing),
        "cache_entries_without_ase_hartree_match": len(missing),
        "first_missing_cache_sample": missing[:1],
        "output_records_without_ase_hartree_match": len(unmatched_output),
        "first_unmatched_output_sample": unmatched_output[:1],
        "n_order_aligned_matches": sum(abs(value) <= tolerance_ev
                                       for value in order_differences),
        "ase_hartree_to_ev": Hartree,
        "effective_order_paired_hartree_to_ev": effective_factor,
        "effective_minus_ase_hartree_to_ev": effective_factor - Hartree,
        "max_absolute_energy_residual_after_factor_ev": float(
            np.max(np.abs(factor_residuals))
        ),
        "first_three_order_differences_ev": order_differences[:3],
        "last_output_energy_ev": output_energies[-1],
        "last_cache_by_mtime": cached[-1],
        "last_three_order_differences_ev": order_differences[-3:],
        "last_output_has_program_end": (
            text.rfind("PROGRAM ENDED AT") > text.rfind("ENERGY| Total FORCE_EVAL")
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("cache_dir", type=Path)
    parser.add_argument("image_index", type=int)
    parser.add_argument("--tolerance-ev", type=float, default=1e-8)
    args = parser.parse_args()
    print(json.dumps(compare(args.output, args.cache_dir, args.image_index,
                             args.tolerance_ev), indent=2))


if __name__ == "__main__":
    main()
