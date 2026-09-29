"""Preflight only missing nodes of a 9x3 symmetry-restricted BTO sheet.

The 12 already measured grid nodes are reused. Fifteen missing Qy=0
conditional optimizations use the unchanged 100-Ry/10-au-DZP ABACUS source.
The old nine-node even-mode fit is frozen before any new DFT evaluation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.preflight_bto_transverse_soft_conditional import build_preflight
from scripts.plan_bto_qy0_even_mode_validation import even_basis


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(source_dir: Path, measured_csv: Path, frozen_plan: Path,
            output_dir: Path) -> dict:
    if output_dir.exists():
        raise FileExistsError(output_dir)
    model = json.loads(frozen_plan.read_text(encoding="utf-8"))
    if (model.get("status") != "two_model_predictions_frozen_before_new_DFT"
            or model.get("n_training_grid_nodes") != 9
            or model.get("fixed_third_soft_mode") != "Q_y=0"
            or model.get("energy_absolute_error_gate_meV_per_BTO") != 2.0):
        raise ValueError("the previously validated even-mode model changed")
    measured = {}
    with measured_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["class"] == "T_endpoint_retrospective":
                continue
            q = (round(float(row["qz"]), 8), round(float(row["qx"]), 8))
            if q in measured:
                raise ValueError(f"duplicate measured Q point {q}")
            measured[q] = row
    grid_qz = [round(0.15 * i, 8) for i in range(9)]
    grid_qx = [0.0, 0.15, 0.3]
    reused = []
    planned = []
    inputs = {
        "report": source_dir / "bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json",
        "reference": source_dir / "cubic_CONTCAR",
        "force_constants": source_dir / "bto_cubic_gamma_force_constants.npz",
        "phonopy_eigenpairs": source_dir / "bto_cubic_phonopy_gamma_eigenpairs.npz",
        "gamma_provenance": source_dir / "bto_cubic_gamma_phonon_provenance.json",
        "force_sets": source_dir / "bto_cubic_gamma_FORCE_SETS",
        "eigenpairs_provenance": source_dir / "bto_cubic_phonopy_gamma_eigenpairs_provenance.json",
        "grid_result_manifest": source_dir / "result_manifest_hf_27777454.json",
    }
    if any(not path.is_file() for path in inputs.values()):
        raise FileNotFoundError("BTO source INPUT/phonon/manifest package incomplete")
    coefficient = np.asarray(model["coefficients_eV_per_BTO"], dtype=float)
    output_dir.mkdir(parents=True)
    (output_dir / "preflights").mkdir()
    for i, qz in enumerate(grid_qz):
        for j, qx in enumerate(grid_qx):
            q = (qz, qx)
            if q in measured:
                reused.append({
                    "grid_index": [i, j], "q_sqrt_amu_A": list(q),
                    "measured_minus_C_meV_per_BTO": float(measured[q]["observed_meV_per_BTO"]),
                    "class": measured[q]["class"],
                })
                continue
            name = f"qz{i:02d}_qx{j:02d}"
            args = argparse.Namespace(
                **inputs, q_parallel=qz, q_transverse=qx,
                seed_amplitude=0.4, lock_third_soft_at_zero=True,
                minimum_distance_A=1.6,
            )
            preflight = build_preflight(args)
            if (preflight.get("kind")
                    != "bto_symmetry_restricted_soft_qy_zero_preflight_no_dft"
                    or preflight.get("phonon_supercell") != [1, 1, 1]
                    or preflight.get("electronic_kpoints") != [4, 4, 4]
                    or preflight.get("q_parallel_q_transverse_sqrt_amu_A") != [qz, qx]):
                raise ValueError(f"Qy=0 preflight failed at {name}")
            target = output_dir / "preflights" / f"{name}.json"
            target.write_text(json.dumps(preflight, indent=2) + "\n", encoding="utf-8")
            planned.append({
                "name": name, "grid_index": [i, j], "q_sqrt_amu_A": list(q),
                "frozen_model_prediction_minus_C_meV_per_BTO": float(
                    1000 * (even_basis(np.asarray(q))[0] @ coefficient)),
                "preflight_sha256": sha256(target),
                "minimum_seed_distance_A": min(
                    row["minimum_atomic_distance_A"] for row in preflight["branch_starts"]),
            })
    if len(reused) != 12 or len(planned) != 15:
        raise ValueError(f"expected 12 reused and 15 new nodes, found {len(reused)}, {len(planned)}")
    contract_ids = {json.loads((output_dir / "preflights" / f"{row['name']}.json")
                               .read_text(encoding="utf-8"))["calculator_id"]
                    for row in planned}
    if len(contract_ids) != 1:
        raise ValueError("new Qy=0 points have different ABACUS contracts")
    result = {
        "status": "15_missing_Qy0_inputs_preflighted_no_DFT",
        "claim_limit": "one restricted local branch, Qy=0; not a full conditional PES",
        "grid_qz_sqrt_amu_A": grid_qz, "grid_qx_sqrt_amu_A": grid_qx,
        "n_grid_nodes": 27, "n_reused_measured_nodes": 12, "n_new_conditional_points": 15,
        "new_point_energy_prediction_gate_meV_per_BTO": 2.0,
        "gradient_tolerance_eV_per_sqrt_amu_A": 0.003,
        "maximum_absolute_stress_kbar": 2.0,
        "calculator_id": next(iter(contract_ids)),
        "reused": reused, "new": planned,
        "source_sha256": {
            "measured_csv": sha256(measured_csv), "frozen_plan": sha256(frozen_plan),
            **{key: sha256(path) for key, path in inputs.items()},
            "preparer": sha256(Path(__file__)),
        },
    }
    (output_dir / "plan.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source-dir", "measured-csv", "frozen-plan", "output-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.source_dir, args.measured_csv, args.frozen_plan,
                     args.output_dir)
    print(json.dumps({key: result[key] for key in (
        "status", "n_reused_measured_nodes", "n_new_conditional_points", "calculator_id",
    )}))


if __name__ == "__main__":
    main()
