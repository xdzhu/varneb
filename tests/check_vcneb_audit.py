"""Regression checks for calculator-free VC-NEB result auditing."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_vcneb_result import audit


def _summary(stress_eV_A3: float) -> dict:
    return {
        "status": "completed",
        "n_images": 3,
        "fmax_target_eV_per_A": 0.02,
        "final_max_generalized_force_eV_per_A": 0.01,
        "path_diagnostics": {
            "n_images": 3,
            "geometry": {"valid": True, "images": []},
            "images": [
                {
                    "image_index": i,
                    "volume_A3": 100.0,
                    "max_atom_force_eV_per_A": 0.01,
                    "max_stress_eV_per_A3": stress_eV_A3,
                }
                for i in range(3)
            ],
        },
    }


def main() -> None:
    good = audit(
        _summary(1.0e-5),
        minimum_distance=1.0,
        maximum_deformation=None,
        maximum_stress_kbar=0.1,
    )
    if good["status"] != "ok" or good["maximum_stress_kbar"] is None:
        raise SystemExit("stress-gated VCNEB audit rejected a valid summary")

    bad = audit(
        _summary(2.0e-4),
        minimum_distance=1.0,
        maximum_deformation=None,
        maximum_stress_kbar=0.1,
    )
    if bad["status"] != "failed" or not any("stress" in issue for issue in bad["issues"]):
        raise SystemExit("stress-gated VCNEB audit missed an over-threshold path")
    loose_summary = _summary(1.0e-5)
    loose_summary["final_max_generalized_force_eV_per_A"] = 0.08
    strict = audit(
        loose_summary,
        minimum_distance=1.0,
        maximum_deformation=None,
    )
    if strict["status"] != "failed" or strict["fmax_target_eV_per_A"] != 0.02:
        raise SystemExit("VCNEB audit did not enforce the recorded force target")
    loose = audit(
        loose_summary,
        minimum_distance=1.0,
        maximum_deformation=None,
        fmax_target=0.10,
    )
    if loose["status"] != "ok" or loose["recorded_fmax_target_eV_per_A"] != 0.02:
        raise SystemExit("VCNEB audit did not allow an explicit loose force target")
    incomplete_summary = _summary(1.0e-5)
    incomplete_summary["status"] = "cancelled"
    incomplete = audit(
        incomplete_summary,
        minimum_distance=1.0,
        maximum_deformation=None,
    )
    if incomplete["status"] != "failed" or not any("completed" in issue for issue in incomplete["issues"]):
        raise SystemExit("VCNEB audit accepted a non-completed summary")
    print("vcneb_audit_regression=ok")


if __name__ == "__main__":
    main()
