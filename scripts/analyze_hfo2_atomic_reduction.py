"""Apply the stable-subspace reducer to two DFT-derived T Gamma matrices.

This is fixed-cell atomic harmonic reduction only. Its predictions are to be
tested by independent SCFs, not certified by fitting the same force constants.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ase.io import read
import numpy as np
from scipy.linalg import null_space

from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.audit_hfo2_static_replica import sha256
from vcneb.relaxed_curvature import relax_orthogonal_curvature


def analyze(variants, gamma, output):
    if output.exists():
        raise FileExistsError("refusing existing atomic-reduction analysis")
    t = read(variants / "T.vasp", format="vasp")
    _, patterns, _ = rotated_t_triplet(t)
    q = patterns.reshape(3, 36).T
    translations = np.tile(np.eye(3), (12, 1)) / np.sqrt(12)
    orthogonal = null_space(np.column_stack([q, translations]).T)
    if orthogonal.shape != (36, 30):
        raise ValueError("expected30 relaxed atomic directions after excluding translations and3 patterns")
    matrices, source_hashes = [], {}
    for size in (.01, .02):
        path = gamma / f"T_d{size:.2f}.npz"
        with np.load(path, allow_pickle=False) as data:
            hessian = data["force_constants_symmetrized_eV_A2"].transpose(0, 2, 1, 3).reshape(36, 36)
        matrices.append(hessian)
        source_hashes[path.name] = sha256(path)
    # This observed spread is an operational gate, NOT a proven error bound.
    spread = float(np.linalg.norm(matrices[0] - matrices[1], ord=2))
    families = {}
    for size, hessian in zip((.01, .02), matrices):
        reduced = relax_orthogonal_curvature(hessian, q, orthogonal, stability_floor=spread)
        families[f"d{size:.2f}"] = {
            "frozen_curvature_eV_A2": reduced.frozen.tolist(),
            "relaxed_curvature_eV_A2": reduced.relaxed.tolist(),
            "softening_eV_A2": reduced.softening.tolist(),
            "full_atomic_response": (orthogonal @ reduced.orthogonal_response).tolist(),
            "minimum_eliminated_curvature_eV_A2": reduced.minimum_eliminated_curvature,
            "eliminated_condition_number": reduced.eliminated_condition_number,
            "x_softening_fraction": float(reduced.softening[0, 0] / reduced.frozen[0, 0])}
    report = {"schema_version": 1, "source_sha256": source_hashes,
              "T_reference_sha256": sha256(variants / "T.vasp"), "retained_basis_cartesian": q.tolist(),
              "eliminated_dimensions": 30, "clamped_cell": True,
              "observed_two_step_operator_spread_eV_A2": spread, "families": families,
              "holdout_design": "independent frozen/linearly responded x-pattern statics at +/-0.05 and +/-0.10 Angstrom",
              "limitations": ["a Schur-complement material pilot, not a new theorem",
                              "two-step spread is not a rigorous total DFT error bound",
                              "response is harmonic, not a fully relaxed conditional surface",
                              "no cell release, barrier prediction, or transition-state certificate"]}
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("variants", "gamma", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.variants, args.gamma, args.output)
    print(json.dumps({"two_step_operator_spread_eV_A2": result["observed_two_step_operator_spread_eV_A2"],
                      "families": {k: {a: v[a] for a in ("minimum_eliminated_curvature_eV_A2", "x_softening_fraction")}
                                   for k, v in result["families"].items()}}))


if __name__ == "__main__":
    main()
