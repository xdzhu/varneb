"""Stage only the 22 unsampled nodes of the audited frozen BTO 9x9 grid.

The old 59 SCF results are reused by coordinate and structure hash. No
electronic calculation or conditional relaxation occurs in this preparer.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.bto_q1q2_reference import sha256


DENSE_KIND = "inert_frozen_cubic_cell_transverse_soft_mode_points_not_T_to_C_barrier"
OLD_KIND = "BTO_transverse_soft_frozen_C_cell_59_real_DFT_points_not_conditional_PES_or_MEP"
NEW_KIND = "inert_frozen_cubic_cell_transverse_soft_missing22_not_T_to_C_barrier"


def stage(dense_dir: Path, analysis59_path: Path, output_dir: Path) -> dict:
    if output_dir.exists():
        raise FileExistsError(output_dir)
    dense_path = dense_dir / "manifest.json"
    dense = json.loads(dense_path.read_text(encoding="utf-8"))
    old = json.loads(analysis59_path.read_text(encoding="utf-8"))
    if (dense.get("kind") != DENSE_KIND or dense.get("n_points") != 81
            or len(dense.get("points", [])) != 81
            or old.get("kind") != OLD_KIND
            or old.get("status") != "eighteen_independent_adaptive_edge_holdouts_scored"
            or old.get("n_real_DFT_points") != 59
            or len(old.get("samples", [])) != 59
            or old.get("source_sha256", {}).get("dense_manifest") != sha256(dense_path)
            or old.get("axis_labels") != dense.get("axis_labels")
            or old.get("axis_mode_weights") != dense.get("axis_mode_weights")
            or old.get("boundary_condition") != dense.get("boundary_condition")):
        raise ValueError("81-point geometry and 59-point DFT audit are not aligned")
    by_index = {tuple(point["grid_index_q1_q2"]): point for point in dense["points"]}
    old_by_index = {(point["dense_i_q1"], point["dense_j_q2"]): point
                    for point in old["samples"]}
    if len(by_index) != 81 or len(old_by_index) != 59:
        raise ValueError("duplicate dense or previously measured coordinate")
    for index, measured in old_by_index.items():
        candidate = by_index.get(index)
        if (candidate is None
                or candidate["structure_sha256"] != measured["source_structure_sha256"]
                or candidate["q1"] != measured["q_parallel_sqrt_amu_A"]
                or candidate["q2"] != measured["q_transverse_sqrt_amu_A"]):
            raise ValueError(f"old DFT point does not match this plane: {index}")
    missing = [point for point in dense["points"]
               if tuple(point["grid_index_q1_q2"]) not in old_by_index]
    if len(missing) != 22:
        raise ValueError(f"expected exactly 22 new nodes, got {len(missing)}")
    for point in missing:
        if sha256(dense_dir / point["structure"]) != point["structure_sha256"]:
            raise ValueError(f"unsampled structure hash changed: {point['name']}")
    output_dir.mkdir(parents=True)
    for point in missing:
        target = output_dir / point["structure"]
        target.parent.mkdir(parents=True)
        shutil.copy2(dense_dir / point["structure"], target)
        if sha256(target) != point["structure_sha256"]:
            raise RuntimeError(f"staged structure changed: {point['name']}")
    result = {
        "kind": NEW_KIND, "n_points": 22, "n_atoms": dense["n_atoms"],
        "amplitude_unit": dense["amplitude_unit"],
        "axis_labels": dense["axis_labels"],
        "axis_definition": dense["axis_definition"],
        "axis_mode_weights": dense["axis_mode_weights"],
        "reference_id": dense["reference_id"],
        "boundary_condition": dense["boundary_condition"],
        "grid": None, "no_dft_launched": True,
        "input_sha256": dense["input_sha256"],
        "source_dense_manifest_sha256": sha256(dense_path),
        "source_59_analysis_sha256": sha256(analysis59_path),
        "selection_rule": "all 9x9 nodes not already among the 59 hash-matched DFT samples",
        "points": missing,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dense-dir", type=Path, required=True)
    parser.add_argument("--analysis59", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = stage(args.dense_dir, args.analysis59, args.output_dir)
    print(json.dumps({"status": "inputs_staged_no_DFT", "n_points": result["n_points"]}))


if __name__ == "__main__":
    main()
