"""Audit independent holdout SCFs against the preregistered harmonic model."""

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.audit_hfo2_static_replica import audited_results, sha256


def audit_point(root, index):
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    if manifest["purpose"] != "HfO2_T_atomic_reduction_8_independent_holdouts":
        raise ValueError("unexpected holdout purpose")
    item = manifest["points"][index]
    directory = root / "calculations" / f"{index:02d}"
    if item["index"] != index or any(sha256(directory / n) != h for n, h in item["input_sha256"].items()):
        raise ValueError("holdout index/input mismatch")
    raw = audited_results(directory)
    vector = np.array(manifest["direction_vectors"][item["kind"]])
    q = np.array(manifest["retained_basis_cartesian"])
    force = raw["forces"].copy()
    force -= force.mean(axis=0)
    force = force.reshape(-1)
    orthogonal_force = force - q @ (q.T @ force)
    report = {**item, "status": "audited", "energy_eV_cell": float(raw["energy"]),
              "directional_gradient_eV_A": float(-raw["forces"].reshape(-1) @ vector),
              "orthogonal_nontranslational_force_norm_eV_A": float(np.linalg.norm(orthogonal_force)),
              "results": {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in raw.items()},
              "raw_log_sha256": sha256(directory / "OUT.ABACUS/running_scf.log")}
    (directory / "point_audit.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def score_pairs(manifest, results):
    pairs = []
    for kind, predicted in manifest["prediction_fixed_before_DFT"].items():
        if not np.isfinite(predicted) or predicted <= 0:
            raise ValueError("this positive-curvature pilot requires a finite positive prediction")
        for h in (.05, .10):
            pair = {p["sign"]: results[p["index"]] for p in manifest["points"]
                    if p["kind"] == kind and p["amplitude_A"] == h and p["index"] in results}
            if set(pair) != {-1, 1}:
                continue
            ke = (pair[1]["energy_eV_cell"] + pair[-1]["energy_eV_cell"] - 2 * manifest["T_energy_eV_cell"]) / h ** 2
            kf = (pair[1]["directional_gradient_eV_A"] - pair[-1]["directional_gradient_eV_A"]) / (2 * h)
            errors = {"energy_prediction": abs(ke - predicted) / predicted,
                      "force_prediction": abs(kf - predicted) / predicted,
                      "energy_force_disagreement": abs(ke - kf) / predicted}
            pairs.append({"kind": kind, "amplitude_A": h, "predicted_curvature_eV_A2": predicted,
                          "energy_curvature_eV_A2": ke, "force_curvature_eV_A2": kf,
                          "relative_errors": errors,
                          "passes_preregistered_10percent_gate": max(errors.values()) <= manifest["relative_curvature_acceptance"]})
    return pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    if args.index is not None:
        report = audit_point(args.root, args.index)
        print(json.dumps({k: report[k] for k in ("index", "status", "energy_eV_cell")}))
    else:
        manifest = json.loads((args.root / "manifest.json").read_text(encoding="utf-8"))
        results = {i: audit_point(args.root, i) for i in range(8)}
        pairs = score_pairs(manifest, results)
        report = {"status": "all8_independent_SCFs_audited", "manifest_sha256": sha256(args.root / "manifest.json"),
                  "pairs": pairs, "all_four_pairs_pass": len(pairs) == 4 and all(p["passes_preregistered_10percent_gate"] for p in pairs),
                  "point_audits": list(results.values()), "limitations": manifest["limitations"]}
        (args.root / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({k: v for k, v in report.items() if k != "point_audits"}))


if __name__ == "__main__":
    main()
