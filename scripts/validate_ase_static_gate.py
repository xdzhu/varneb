"""Validate a calculator-specific ASE endpoint static summary before VC-NEB."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


EV_A3_TO_KBAR = 1602.176634
GPA_PER_EV_A3 = 160.21766208


def validate_summary(
    summary: dict,
    *,
    fmax_eV_per_A: float = 0.10,
    stress_kbar: float = 1.0,
    pressure_gpa: float = 0.0,
) -> dict:
    """Return a machine-readable gate report or raise on invalid endpoints."""

    issues: list[str] = []
    if not math.isfinite(pressure_gpa):
        raise ValueError("target pressure must be finite")
    declared_pressure = summary.get("external_pressure_gpa")
    if pressure_gpa != 0.0 and declared_pressure is None:
        issues.append("nonzero target pressure requires external_pressure_gpa provenance")
    elif declared_pressure is not None:
        try:
            if not math.isclose(float(declared_pressure), pressure_gpa, abs_tol=1e-6):
                issues.append("summary external pressure differs from gate target")
        except (TypeError, ValueError):
            issues.append("summary external pressure is invalid")
    if summary.get("status") != "static_completed":
        issues.append("summary status is not static_completed")
    endpoints = summary.get("static_endpoints")
    if not isinstance(endpoints, list) or len(endpoints) != 2:
        issues.append("summary must contain exactly two static endpoints")
        endpoints = []

    identity = summary.get("endpoint_structures")
    if not isinstance(identity, dict):
        issues.append("summary must contain endpoint_structures identity records")
        identity = {}
    initial_identity = identity.get("initial")
    final_identity = identity.get("final")
    if not isinstance(initial_identity, dict) or not isinstance(final_identity, dict):
        issues.append("endpoint_structures must contain initial and final records")
    elif initial_identity.get("n_atoms") != final_identity.get("n_atoms"):
        issues.append("initial/final endpoint atom counts differ")
    elif initial_identity.get("composition") != final_identity.get("composition"):
        issues.append("initial/final endpoint compositions differ")

    records = []
    labels_seen: set[str] = set()
    for endpoint in endpoints:
        if not isinstance(endpoint, dict):
            issues.append("each static endpoint must be an object")
            continue
        label = str(endpoint.get("label", endpoint.get("index", "unknown")))
        labels_seen.add(label)
        force = endpoint.get("max_generalized_force_eV_per_A")
        if force is None:
            force = endpoint.get("max_force_eV_per_A")
        stress = endpoint.get("stress_eV_per_A3")
        try:
            force_value = float(force)
        except (TypeError, ValueError):
            issues.append(f"{label}: missing finite endpoint force")
            continue
        if not math.isfinite(force_value) or not force_value < fmax_eV_per_A:
            issues.append(f"{label}: force {force_value:.8g} > {fmax_eV_per_A:.8g} eV/A")
        try:
            rows = [[float(value) for value in row] for row in stress]
            if len(rows) != 3 or any(len(row) != 3 for row in rows):
                raise ValueError("stress must be 3x3")
            if any(not math.isfinite(value) for row in rows for value in row):
                raise ValueError("stress must be finite")
        except (TypeError, ValueError):
            issues.append(f"{label}: missing finite endpoint stress")
            continue
        stress_value_kbar = max(abs(value) for row in rows for value in row) * EV_A3_TO_KBAR
        target_stress = -pressure_gpa / GPA_PER_EV_A3
        residual_kbar = max(
            abs(rows[i][j] - (target_stress if i == j else 0.0))
            for i in range(3) for j in range(3)
        ) * EV_A3_TO_KBAR
        if not residual_kbar < stress_kbar:
            issues.append(
                f"{label}: stress residual {residual_kbar:.8g} > {stress_kbar:.8g} kbar"
            )
        records.append(
            {
                "label": label,
                "max_generalized_force_eV_per_A": force_value,
                "max_abs_stress_kbar": stress_value_kbar,
                "max_abs_stress_residual_kbar": residual_kbar,
                "n_atoms": endpoint.get("n_atoms"),
            }
        )
    if labels_seen != {"initial", "final"}:
        issues.append("static endpoints must be labelled initial and final exactly once")
    expected_atoms = initial_identity.get("n_atoms") if isinstance(initial_identity, dict) else None
    if expected_atoms is not None:
        for record in records:
            if record["n_atoms"] != expected_atoms:
                issues.append(f"{record['label']}: atom count differs from endpoint identity")

    return {
        "status": "passed" if not issues else "failed",
        "fmax_target_eV_per_A": fmax_eV_per_A,
        "stress_target_kbar": stress_kbar,
        "external_pressure_gpa": pressure_gpa,
        "endpoints": records,
        "issues": issues,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--fmax", type=float, default=0.10)
    parser.add_argument("--stress-kbar", type=float, default=1.0)
    parser.add_argument("--pressure-gpa", type=float, default=0.0)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    report = validate_summary(
        summary,
        fmax_eV_per_A=args.fmax,
        stress_kbar=args.stress_kbar,
        pressure_gpa=args.pressure_gpa,
    )
    if args.output is not None:
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if report["status"] != "passed":
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
