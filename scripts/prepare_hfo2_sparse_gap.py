"""Three fixed-contract statics inside the sparsest guided-chain segment."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from ase.io import read
import numpy as np

from scripts.audit_hfo2_static_replica import INPUT_FILES, audited_results, sha256
from scripts.prepare_hfo2_channel_work_probes import write_structure
from scripts.prepare_hfo2_gamma_and_seeds import geometry_roundtrip
from vcneb import interpolate_vcneb


def gap_geometries(left, right):
    return interpolate_vcneb(left, right, 5, align_cells=False, align_translation=False,
                             mic=True, mapping=None, minimum_distance=1.6)[1:-1]


def prepare(source_root, chain, baseline_source, output):
    if output.exists():
        raise FileExistsError("refusing existing gap namespace")
    endpoint_records, endpoints = [], []
    for index in (1, 2):
        source = source_root / f"{index:02d}"
        actual = read(source / "STRU", format="abacus")
        path = chain / f"POSCAR_{index:02d}"
        geometry = read(path, format="vasp")
        delta = geometry.get_scaled_positions() - actual.get_scaled_positions()
        delta -= np.rint(delta)
        if not np.allclose(geometry.cell.array, actual.cell.array, atol=1e-9) or np.max(np.abs(delta)) > 1e-9:
            raise ValueError("gap source static differs from final guided snapshot")
        raw = audited_results(source)
        endpoint_records.append({"index": index, "energy_eV_cell": float(raw["energy"]),
                                 "STRU_sha256": sha256(source / "STRU"), "POSCAR_sha256": sha256(path),
                                 "log_sha256": sha256(source / "OUT.ABACUS/running_scf.log")})
        endpoints.append(geometry)
    source = source_root / "02"
    for name in INPUT_FILES:
        if name != "STRU" and sha256(source / name) != sha256(baseline_source / name):
            raise ValueError("gap and T baseline electronic contracts differ")
    baseline = audited_results(baseline_source)
    output.mkdir(parents=True)
    points = []
    for index, geometry in enumerate(gap_geometries(*endpoints)):
        target = output / "points" / f"{index:02d}"
        target.mkdir(parents=True)
        for name in INPUT_FILES:
            if name != "STRU":
                shutil.copyfile(source / name, target / name)
        write_structure(target / "STRU", geometry)
        minimum = geometry_roundtrip(target / "STRU", geometry)
        points.append({"index": index, "fraction": (index + 1) / 4,
                       "minimum_distance_A": minimum,
                       "input_sha256": {n: sha256(target / n) for n in INPUT_FILES}})
    manifest = {"schema_version": 1, "purpose": "HfO2_guided_step300_gap01_02_three_statics",
                "source_root": str(source_root), "chain": str(chain),
                "baseline_source": str(baseline_source), "T_energy_eV_cell": float(baseline["energy"]),
                "baseline_log_sha256": sha256(baseline_source / "OUT.ABACUS/running_scf.log"),
                "source_endpoints": endpoint_records, "points": points,
                "interpolation": "linear cell/fractional coordinates, nearest link unwrap, unchanged atom order, no rotation/translation alignment",
                "limitations": "static sampled interpolation only; neither a new relaxed MEP nor certification of its global saddle"}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source-root", "chain", "baseline-source", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    manifest = prepare(args.source_root, args.chain, args.baseline_source, args.output)
    print(json.dumps({"status": "staged_no_DFT", "n_points": len(manifest["points"])}))


if __name__ == "__main__":
    main()
