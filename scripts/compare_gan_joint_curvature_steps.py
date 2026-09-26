"""Compare two audited GaN image-15 joint-Hessian displacement sizes.

This tests numerical stability of a *nonstationary* VCNEB maximum image; it
never upgrades that image to an independently verified saddle point.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(first_dir: Path, second_dir: Path) -> dict:
    reports = [json.loads((directory / "audit.json").read_text(encoding="utf-8"))
               for directory in (first_dir, second_dir)]
    archives = [np.load(directory / "joint_hessian.npz", allow_pickle=False)
                for directory in (first_dir, second_dir)]
    try:
        for report in reports:
            if (report["status"] != "GaN_image15_local_joint_curvature_at_nonstationary_NEB_candidate"
                    or report["n_static_displacements"] != 36):
                raise ValueError("a source step has not passed the full static audit")
        if not reports[0]["step_A"] < reports[1]["step_A"]:
            raise ValueError("steps must be distinct and passed in ascending size")
        for key in ("trajectory", "summary", "joint_curvature", "center_OUTCAR"):
            if reports[0]["source_sha256"][key] != reports[1]["source_sha256"][key]:
                raise ValueError(f"two step-size audits have different {key} sources")
        h0, h1 = (archive["hessian"] for archive in archives)
        ev0, ev1 = (archive["eigenvalues"] for archive in archives)
        modes0, modes1 = (archive["eigenvectors"] for archive in archives)
        if (h0.shape != h1.shape or h0.shape != (18, 18)
                or ev0.shape != ev1.shape or ev0.shape != (15,)
                or modes0.shape != modes1.shape or modes0.shape != (18, 15)
                or not all(np.isfinite(item).all() for item in (h0, h1, ev0, ev1, modes0, modes1))):
            raise ValueError("incompatible or non-finite joint-Hessian arrays")
        negative0 = int(np.count_nonzero(ev0 < -0.1))
        negative1 = int(np.count_nonzero(ev1 < -0.1))
        lowest_mode_overlap = float(abs(modes0[:, 0] @ modes1[:, 0]) ** 2)
        report = {
            "status": "GaN_nonstationary_joint_curvature_two_step_numerical_comparison_not_TS_certificate",
            "step_A": [reports[0]["step_A"], reports[1]["step_A"]],
            "negative_count_below_minus_0p1_eV_per_A2": [negative0, negative1],
            "lowest_eigenvalue_eV_per_A2": [float(ev0[0]), float(ev1[0])],
            "next_eigenvalue_eV_per_A2": [float(ev0[1]), float(ev1[1])],
            "lowest_mode_squared_overlap_between_steps": lowest_mode_overlap,
            "maximum_sorted_eigenvalue_difference_eV_per_A2": float(np.max(np.abs(ev0 - ev1))),
            "relative_hessian_frobenius_difference": float(
                np.linalg.norm(h0 - h1) / max(np.linalg.norm(h0), 1e-30)
            ),
            "raw_reciprocity_defect": [r["raw_hessian_reciprocity_relative_defect"] for r in reports],
            "translation_null_defect": [r["translation_null_relative_defect"] for r in reports],
            "energy_gradient_max_abs_difference_eV_per_A": [
                r["energy_gradient_max_abs_difference_eV_per_A"] for r in reports
            ],
            "energy_hessian_diagonal_max_abs_difference_eV_per_A2": [
                r["energy_hessian_diagonal_max_abs_difference_eV_per_A2"] for r in reports
            ],
            "tangent_squared_overlap_with_lowest_mode": [
                r["tangent_squared_overlap_with_eigenvectors"][0] for r in reports
            ],
            "center_gradient_translation_free_eV_per_A": [
                r["center_gradient_translation_free_eV_per_A"] for r in reports
            ],
            "source_sha256": {
                "audits": [sha256(directory / "audit.json")
                           for directory in (first_dir, second_dir)],
                "hessian_archives": [sha256(directory / "joint_hessian.npz")
                                     for directory in (first_dir, second_dir)],
                "comparator": sha256(Path(__file__)),
            },
            "limitation": (
                "Step-size-stable curvature at a nonstationary candidate does not prove a "
                "variable-cell transition state or the two basin connections."
            ),
        }
        return report
    finally:
        for archive in archives:
            archive.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small", type=Path, required=True)
    parser.add_argument("--large", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = compare(args.small, args.large)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "status", "lowest_eigenvalue_eV_per_A2",
        "lowest_mode_squared_overlap_between_steps",
        "relative_hessian_frobenius_difference",
    )}))


if __name__ == "__main__":
    main()
