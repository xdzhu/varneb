"""Audit a three-or-more-point QE fixed-endpoint cutoff/SCF series.

The inputs are ``qe_static_summary.json`` files produced by
``examples/run_vcneb_qe.py --static-only``.  This program deliberately does
not choose tolerances: the energy, force and stress tolerances must be stated
on the command line before the series is inspected.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import sys
import tempfile
from typing import Any

import numpy as np


EV_A3_TO_GPA = 160.2176634


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", action="append", required=True, type=Path, help="QE static summary; repeat at least three times")
    parser.add_argument("--energy-tol-mev-per-atom", required=True, type=float)
    parser.add_argument("--force-tol-ev-per-a", required=True, type=float)
    parser.add_argument("--stress-tol-gpa", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def _read_summary(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} does not contain a JSON object")
    payload["_source"] = str(path.resolve())
    return payload


def _require_static_summary(summary: dict[str, Any]) -> None:
    source = summary["_source"]
    if summary.get("status") != "completed":
        raise ValueError(f"{source}: status must be 'completed'")
    if summary.get("execution_mode") != "fixed_initial_endpoint_static_scf":
        raise ValueError(f"{source}: not a fixed-initial-endpoint static SCF summary")
    if summary.get("evaluated_image_index") != 0:
        raise ValueError(f"{source}: static baseline must evaluate image index 0")
    if summary.get("n_images") != 7 or summary.get("n_interior_images") != 5:
        raise ValueError(f"{source}: expected the BTO 7-total-image layout")
    for key in ("potential_energy_eV", "forces_eV_per_A", "stress_eV_per_A3_voigt"):
        if key not in summary:
            raise ValueError(f"{source}: missing {key}")


def _settings_signature(summary: dict[str, Any]) -> dict[str, Any]:
    """Return the conditions that must not change across cutoff points."""

    parameters = deepcopy(summary.get("calculator_parameters") or {})
    input_data = parameters.get("input_data") or {}
    system = dict(input_data.get("system") or {})
    system.pop("ecutwfc", None)
    system.pop("ecutrho", None)
    input_data["system"] = system
    parameters["input_data"] = input_data
    pseudo = []
    for item in summary.get("pseudopotential_reports") or []:
        pseudo.append({key: item.get(key) for key in ("species", "filename", "md5", "expected_md5")})
    pseudo.sort(key=lambda item: str(item["species"]))
    endpoints = summary.get("endpoint_structures") or {}
    return {
        "git_revision": summary.get("git_revision"),
        "n_images": summary.get("n_images"),
        "n_interior_images": summary.get("n_interior_images"),
        "initial_endpoint_sha256": (endpoints.get("initial") or {}).get("sha256"),
        "final_endpoint_sha256": (endpoints.get("final") or {}).get("sha256"),
        "calculator_parameters_except_cutoffs": parameters,
        "pseudopotentials": pseudo,
    }


def _cutoffs(summary: dict[str, Any]) -> tuple[float, float]:
    system = ((summary.get("calculator_parameters") or {}).get("input_data") or {}).get("system") or {}
    try:
        wavefunction = float(system["ecutwfc"])
        density = float(system["ecutrho"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"{summary['_source']}: missing numeric QE ecutwfc/ecutrho") from exc
    if wavefunction <= 0.0 or density < wavefunction:
        raise ValueError(f"{summary['_source']}: invalid cutoffs ecutwfc={wavefunction}, ecutrho={density}")
    return wavefunction, density


def _compare(lower: dict[str, Any], higher: dict[str, Any], *, natoms: int) -> dict[str, Any]:
    energy_difference = abs(float(higher["potential_energy_eV"]) - float(lower["potential_energy_eV"]))
    forces_lower = np.asarray(lower["forces_eV_per_A"], dtype=float)
    forces_higher = np.asarray(higher["forces_eV_per_A"], dtype=float)
    if forces_lower.shape != forces_higher.shape or forces_lower.shape != (natoms, 3):
        raise ValueError(f"force-array shape mismatch between {lower['_source']} and {higher['_source']}")
    stress_lower = np.asarray(lower["stress_eV_per_A3_voigt"], dtype=float)
    stress_higher = np.asarray(higher["stress_eV_per_A3_voigt"], dtype=float)
    if stress_lower.shape != stress_higher.shape or stress_lower.size != 6:
        raise ValueError(f"stress-array shape mismatch between {lower['_source']} and {higher['_source']}")
    return {
        "lower_cutoffs_Ry": _cutoffs(lower),
        "higher_cutoffs_Ry": _cutoffs(higher),
        "energy_difference_meV_per_atom": 1000.0 * energy_difference / natoms,
        "max_atom_force_difference_eV_per_A": float(np.linalg.norm(forces_higher - forces_lower, axis=1).max()),
        "max_stress_component_difference_GPa": float(np.abs(stress_higher - stress_lower).max() * EV_A3_TO_GPA),
    }


def audit(
    summaries: list[dict[str, Any]],
    *,
    energy_tol_mev_per_atom: float,
    force_tol_ev_per_a: float,
    stress_tol_gpa: float,
) -> dict[str, Any]:
    """Check a fixed-geometry QE series and accept only its final cutoff pair."""

    tolerances = {
        "energy_meV_per_atom": float(energy_tol_mev_per_atom),
        "force_eV_per_A": float(force_tol_ev_per_a),
        "stress_GPa": float(stress_tol_gpa),
    }
    if any(value < 0.0 for value in tolerances.values()):
        raise ValueError("all convergence tolerances must be non-negative")
    if len(summaries) < 3:
        raise ValueError("at least three static summaries are required for a cutoff series")
    for summary in summaries:
        _require_static_summary(summary)
    signatures = [_settings_signature(summary) for summary in summaries]
    if any(signature != signatures[0] for signature in signatures[1:]):
        raise ValueError("static summaries do not share identical endpoints, code, UPFs and non-cutoff QE settings")
    endpoint = (summaries[0].get("endpoint_structures") or {}).get("initial") or {}
    natoms = endpoint.get("n_atoms")
    if not isinstance(natoms, int) or natoms <= 0:
        raise ValueError("static summaries lack a valid initial endpoint atom count")
    ordered = sorted(summaries, key=_cutoffs)
    cutoff_pairs = [_cutoffs(summary) for summary in ordered]
    if len(set(cutoff_pairs)) != len(cutoff_pairs):
        raise ValueError("each cutoff point must have a unique (ecutwfc, ecutrho) pair")
    comparisons = [_compare(low, high, natoms=natoms) for low, high in zip(ordered, ordered[1:])]
    final = comparisons[-1]
    accepted = (
        final["energy_difference_meV_per_atom"] <= tolerances["energy_meV_per_atom"]
        and final["max_atom_force_difference_eV_per_A"] <= tolerances["force_eV_per_A"]
        and final["max_stress_component_difference_GPa"] <= tolerances["stress_GPa"]
    )
    return {
        "status": "ok" if accepted else "failed",
        "acceptance_policy": "The highest and penultimate cutoff points must agree within declared tolerances; every earlier adjacent comparison is retained for inspection.",
        "declared_tolerances": tolerances,
        "fixed_conditions": signatures[0],
        "n_static_points": len(ordered),
        "cutoff_points_Ry": cutoff_pairs,
        "adjacent_comparisons": comparisons,
        "accepted_final_pair": final,
        "sources": [summary["_source"] for summary in ordered],
    }


def write_report(path: Path, report: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def main() -> int:
    args = parse_args()
    try:
        summaries = [_read_summary(path) for path in args.summary]
        report = audit(
            summaries,
            energy_tol_mev_per_atom=args.energy_tol_mev_per_atom,
            force_tol_ev_per_a=args.force_tol_ev_per_a,
            stress_tol_gpa=args.stress_tol_gpa,
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    write_report(args.output, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
