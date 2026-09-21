"""Combine two calculator-agnostic endpoint relax summaries for the static gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def _load(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"summary must be a JSON object: {path}")
    if payload.get("status") != "completed" or payload.get("converged") is not True:
        raise ValueError(f"endpoint relaxation is not converged: {path}")
    endpoint = payload.get("endpoint")
    if not isinstance(endpoint, dict):
        raise ValueError(f"endpoint identity is missing: {path}")
    stress = payload.get("stress_eV_per_A3")
    if not isinstance(stress, list):
        raise ValueError(f"endpoint stress is missing: {path}")
    return payload


def build_summary(initial_path: str | Path, final_path: str | Path) -> dict:
    initial = _load(Path(initial_path).resolve())
    final = _load(Path(final_path).resolve())
    initial_identity = initial["endpoint"]
    final_identity = final["endpoint"]
    if initial_identity.get("n_atoms") != final_identity.get("n_atoms"):
        raise ValueError("relaxed endpoints have different atom counts")
    if initial_identity.get("composition") != final_identity.get("composition"):
        raise ValueError("relaxed endpoints have different compositions")

    def record(payload: dict, label: str) -> dict:
        return {
            "index": 0 if label == "initial" else 1,
            "label": label,
            "energy_eV": payload["potential_energy_eV"],
            "max_force_eV_per_A": payload["max_atomic_force_eV_per_A"],
            "max_generalized_force_eV_per_A": payload["max_generalized_force_eV_per_A"],
            "stress_eV_per_A3": payload["stress_eV_per_A3"],
            "n_atoms": payload["endpoint"]["n_atoms"],
        }

    return {
        "status": "static_completed",
        "endpoint_evaluation_policy": "independent_bfgs_relaxations",
        "endpoint_structures": {"initial": initial_identity, "final": final_identity},
        "static_endpoints": [record(initial, "initial"), record(final, "final")],
        "source_relax_summaries": {
            "initial": str(Path(initial_path).resolve()),
            "final": str(Path(final_path).resolve()),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial", required=True, type=Path)
    parser.add_argument("--final", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    summary = build_summary(args.initial, args.final)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
