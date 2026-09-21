"""Validate a calculator-specific ASE endpoint static summary before VC-NEB."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


EV_A3_TO_KBAR = 1602.176634


def validate_summary(
    summary: dict,
    *,
    fmax_eV_per_A: float = 0.10,
    stress_kbar: float = 0.10,
) -> dict:
    """Return a machine-readable gate report or raise on invalid endpoints."""

    issues: list[str] = []
    if summary.get("status") != "static_completed":
        issues.append("summary status is not static_completed")
    endpoints = summary.get("static_endpoints")
    if not isinstance(endpoints, list) or len(endpoints) != 2:
        issues.append("summary must contain exactly two static endpoints")
        endpoints = []

    records = []
    for endpoint in endpoints:
        label = str(endpoint.get("label", endpoint.get("index", "unknown")))
        force = endpoint.get("max_generalized_force_eV_per_A")
        if force is None:
            force = endpoint.get("max_force_eV_per_A")
        stress = endpoint.get("stress_eV_per_A3")
        try:
            force_value = float(force)
        except (TypeError, ValueError):
            issues.append(f"{label}: missing finite endpoint force")
            continue
        if not force_value < fmax_eV_per_A:
            issues.append(f"{label}: force {force_value:.8g} > {fmax_eV_per_A:.8g} eV/A")
        try:
            stress_value = max(abs(float(value)) for row in stress for value in row)
        except (TypeError, ValueError):
            issues.append(f"{label}: missing finite endpoint stress")
            continue
        stress_value_kbar = stress_value * EV_A3_TO_KBAR
        if not stress_value_kbar < stress_kbar:
            issues.append(
                f"{label}: stress {stress_value_kbar:.8g} > {stress_kbar:.8g} kbar"
            )
        records.append(
            {
                "label": label,
                "max_generalized_force_eV_per_A": force_value,
                "max_abs_stress_kbar": stress_value_kbar,
                "n_atoms": endpoint.get("n_atoms"),
            }
        )

    return {
        "status": "passed" if not issues else "failed",
        "fmax_target_eV_per_A": fmax_eV_per_A,
        "stress_target_kbar": stress_kbar,
        "endpoints": records,
        "issues": issues,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--fmax", type=float, default=0.10)
    parser.add_argument("--stress-kbar", type=float, default=0.10)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    report = validate_summary(
        summary,
        fmax_eV_per_A=args.fmax,
        stress_kbar=args.stress_kbar,
    )
    if args.output is not None:
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if report["status"] != "passed":
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
